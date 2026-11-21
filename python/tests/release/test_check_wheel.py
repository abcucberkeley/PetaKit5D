"""The release checks of the wheels and the source distribution (tools/check_wheel.py), on small
wheels made here with ELF, Mach-O and PE files laid out by hand, so that every platform's checks
run on every platform; and its readers on this build's own modules."""

from __future__ import annotations

import importlib.util
import io
import platform
import re
import struct
import sys
import tarfile
import zipfile
import zlib
from pathlib import Path

import pytest

TOOL = Path(__file__).resolve().parents[2] / "tools" / "check_wheel.py"


def _tool():
    # (registered in sys.modules, as its dataclasses need)
    spec = importlib.util.spec_from_file_location("pk5d_check_wheel", TOOL)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


cw = _tool()

# --- binaries, laid out as their formats define them ---


def _pad(b: bytearray, n: int) -> None:
    b.extend(bytes(-len(b) % n))


def elf(needed=(), runpath=None, rpath=None, soname=None, symbols=(), extra=b"", debug=b"",
        machine=62) -> bytes:  # fmt: skip
    """A 64-bit little-endian shared object: .dynstr, .dynamic, .dynsym (the symbols, defined,
    and an undefined fftw symbol), .rodata (extra) and, with debug, a zlib-compressed
    .debug_str."""
    dynstr = bytearray(b"\0")

    def s(text: str) -> int:
        dynstr.extend(text.encode() + b"\0")
        return len(dynstr) - len(text) - 1

    dyn = [(1, s(n)) for n in needed]
    dyn += [(14, s(soname))] if soname else []
    dyn += [(15, s(rpath))] if rpath is not None else []
    dyn += [(29, s(runpath))] if runpath is not None else []
    dynamic = b"".join(struct.pack("<qQ", t, v) for t, v in [*dyn, (0, 0)])
    syms = [bytes(24)] + [struct.pack("<IBBHQQ", s(n), 0x12, 0, 4, 0, 0) for n in symbols]
    syms.append(struct.pack("<IBBHQQ", s("fftw_plan_elsewhere"), 0x12, 0, 0, 0, 0))
    sections = [
        (".dynstr", 3, 0, bytes(dynstr), 0),
        (".dynamic", 6, 0, dynamic, 1),
        (".dynsym", 11, 0, b"".join(syms), 1),
        (".rodata", 1, 0, extra, 0),
    ]
    if debug:
        chdr = struct.pack("<IIQQ", 1, 0, len(debug), 1)  # ELFCOMPRESS_ZLIB
        sections.append((".debug_str", 1, 0x800, chdr + zlib.compress(debug), 0))
    sections.append((".shstrtab", 3, 0, None, 0))
    out, heads, shstr = bytearray(64), [bytes(64)], bytearray(b"\0")
    for name, typ, flags, blob, link in sections:
        shstr.extend(name.encode() + b"\0")
        name_at = len(shstr) - len(name) - 1
        blob = bytes(shstr) if blob is None else blob
        _pad(out, 8)
        heads.append(
            struct.pack("<IIQQQQIIQQ", name_at, typ, flags, 0, len(out), len(blob), link, 0, 8, 0)
        )
        out.extend(blob)
    _pad(out, 8)
    shoff = len(out)
    out.extend(b"".join(heads))
    n = len(heads)
    header = struct.pack("<HHIQQQIHHHHHH", 3, machine, 1, 0, 0, shoff, 0, 64, 0, 0, 64, n, n - 1)
    out[:64] = b"\x7fELF" + bytes([2, 1, 1]) + bytes(9) + header
    return bytes(out)


def _lc(cmd: int, fields: bytes, name: str) -> bytes:
    """A load command whose string follows its fields, padded to 8 bytes."""
    raw = fields + name.encode() + b"\0"
    raw += bytes(-(len(raw) + 8) % 8)
    return struct.pack("<II", cmd, len(raw) + 8) + raw


def macho(dylibs=(), rpaths=(), ident=None, minos=(13, 0, 0), symbols=(), extra=b"",
          cpu=0x0100000C) -> bytes:  # fmt: skip
    """A 64-bit little-endian Mach-O dylib: LC_ID_DYLIB, LC_LOAD_DYLIB, LC_RPATH, LC_BUILD_VERSION
    and LC_SYMTAB (the symbols, defined, and an undefined fftw symbol)."""
    versions = struct.pack("<IIII", 24, 0, 0x10000, 0x10000)
    cmds = [_lc(0xD, versions, ident)] if ident else []
    cmds += [_lc(0xC, versions, d) for d in dylibs]
    cmds += [_lc(0x8000001C, struct.pack("<I", 12), r) for r in rpaths]
    if minos:
        v = minos[0] << 16 | minos[1] << 8 | minos[2]
        cmds.append(struct.pack("<IIIIII", 0x32, 24, 1, v, v, 0))
    strtab, nlist = bytearray(b" \0"), []
    for name, ntype in [*((n, 0x0F) for n in symbols), ("_fftw_plan_elsewhere", 0x01)]:
        nlist.append(struct.pack("<IBBHQ", len(strtab), ntype, ntype == 0x0F, 0, 0))
        strtab.extend(name.encode() + b"\0")
    sizeofcmds = sum(map(len, cmds)) + 24
    symoff = 32 + sizeofcmds
    stroff = symoff + 16 * len(nlist)
    cmds.append(struct.pack("<IIIIII", 0x2, 24, symoff, len(nlist), stroff, len(strtab)))
    head = struct.pack("<IiiIIIII", 0xFEEDFACF, cpu, 0, 6, len(cmds), sizeofcmds, 0, 0)
    return head + b"".join(cmds) + b"".join(nlist) + bytes(strtab) + extra


def fat(*slices: tuple[int, bytes]) -> bytes:
    """A universal binary of (CPU type, Mach-O) slices, each at a 4096-byte boundary."""
    out = bytearray(struct.pack(">II", 0xCAFEBABE, len(slices)) + bytes(20 * len(slices)))
    for i, (cpu, data) in enumerate(slices):
        _pad(out, 4096)
        struct.pack_into(">iiIII", out, 8 + 20 * i, cpu, 0, len(out), len(data), 12)
        out.extend(data)
    return bytes(out)


def pe(imports=(), delay=(), exports=(), extra=b"", machine=0x8664, signer=None,
       linker=(14, 44)) -> bytes:  # fmt: skip
    """A PE32+ DLL with one section holding its import, delay-import and export tables, and with
    signer, a certificate table (an Authenticode entry whose blob names the signer); linker: the
    linker's version in the optional header (MSVC's 14.x; GNU ld writes its own, 2.x)."""
    rva0, raw0 = 0x1000, 0x400
    sec = bytearray()

    def put(b: bytes) -> int:
        sec.extend(b)
        return rva0 + len(sec) - len(b)

    imp = [put(n.encode() + b"\0") for n in imports]
    dly = [put(n.encode() + b"\0") for n in delay]
    exp = [put(n.encode() + b"\0") for n in exports]
    _pad(sec, 8)
    imp_at = put(b"".join(struct.pack("<5I", 0, 0, 0, n, 0x2000) for n in imp) + bytes(20))
    dly_at = put(b"".join(struct.pack("<8I", 1, n, 0, 0, 0, 0, 0, 0) for n in dly) + bytes(32))
    names_at = put(b"".join(struct.pack("<I", n) for n in exp))
    exp_dir = struct.pack("<IIHHIIIIIII", 0, 0, 0, 0, 0, 1, len(exp), len(exp), 0, names_at, 0)
    exp_at = put(exp_dir)
    sec.extend(extra)
    dirs = [(0, 0)] * 16
    if exports:
        dirs[0] = (exp_at, 40)
    if imports:
        dirs[1] = (imp_at, 20 * (len(imp) + 1))
    if delay:
        dirs[13] = (dly_at, 32 * (len(dly) + 1))
    opt = struct.pack("<HBB", 0x20B, *linker) + bytes(104) + struct.pack("<I", 16)
    opt += b"".join(struct.pack("<II", *d) for d in dirs)
    dos = bytearray(0x40)
    dos[:2] = b"MZ"
    struct.pack_into("<I", dos, 0x3C, 0x40)
    coff = struct.pack("<HHIIIHH", machine, 1, 0, 0, 0, len(opt), 0x2022)
    section = struct.pack("<8sIIIIIIHHI", b".rdata", len(sec), rva0, len(sec), raw0, 0, 0, 0, 0, 0)
    head = bytes(dos) + b"PE\0\0" + coff + opt + section
    body = head + bytes(raw0 - len(head)) + bytes(sec)
    if signer is None:
        return body
    # WIN_CERTIFICATE: length, revision 2.0, PKCS#7 signed data (here only the signer's name)
    blob = b"\x30\x82" + b"O=" + signer.encode() + bytes(-len(signer) % 8)
    cert = struct.pack("<IHH", 8 + len(blob), 0x200, 2) + blob
    out = bytearray(body + bytes(-len(body) % 8))
    struct.pack_into("<II", out, 0x40 + 24 + 112 + 8 * 4, len(out), len(cert))  # entry 4
    return bytes(out) + cert


# --- wheels and source distributions ---

LICENSES = ("LICENSE", "THIRD_PARTY_NOTICES.md", "cpp/fft/third_party/pocketfft/LICENSE.md")
HASH = "0123456789abcdef0123456789abcdef"
GOMP = "libgomp-e985bcbb.so.1.0.0"
PLATFORMS = {
    "linux": "manylinux_2_24_x86_64.manylinux_2_28_x86_64",
    "macos": "macosx_13_0_arm64",
    "windows": "win_amd64",
}
OMP = "@loader_path/../.dylibs/libomp.dylib"
# a sample home folder, composed so that the repository's text scans don't take it for a real one
HOME = "/".join(("", "home", "builder", ""))


def _metadata(dist, version, license_expr, license_files) -> bytes:
    text = f"Metadata-Version: 2.4\nName: {dist}\nVersion: {version}\n"
    text += f"License-Expression: {license_expr}\n" if license_expr else ""
    return (text + "".join(f"License-File: {f}\n" for f in license_files)).encode()


def linux_module(**kw) -> bytes:
    return elf(**{"needed": [GOMP, "libm.so.6", "libc.so.6"], "runpath": "$ORIGIN", **kw})


# delvewheel's copy of the OpenMP runtime, with a mangled name, as the windows-2025 runners'
# toolset leaves it to delvewheel
OMP140 = f"libomp140.x86_64-{HASH[::-1]}.dll"
MS = "Microsoft Corporation"
MS_PATH = b"D:\\a\\_work\\1\\s\\src\\vctools\\crt\\github\\stl\\src\\x.cpp\0"


def windows_module(*imports, omp=OMP140, **kw) -> bytes:
    std = [omp, f"msvcp140-{HASH}.dll", "KERNEL32.dll", "VCRUNTIME140.dll", "VCRUNTIME140_1.dll",
           "api-ms-win-crt-runtime-l1-1-0.dll"]  # fmt: skip
    return pe(imports=[*(imports or ["python3.dll"]), *std], exports=["PyInit_core"], **kw)


def petakit5d_files(system: str) -> dict[str, bytes]:
    """The members of a clean petakit5d cp312-abi3 wheel of the system."""
    files = {"petakit5d/__init__.py": b'"""petakit5d"""\n'}
    if system == "linux":
        files["petakit5d/_ext/core.abi3.so"] = linux_module(symbols=["PyInit_core"])
        files[f"petakit5d/_ext/{GOMP}"] = elf(needed=["libc.so.6"], soname=GOMP)
    elif system == "macos":
        libcxx = "/usr/lib/libc++.1.dylib"
        files["petakit5d/_ext/core.abi3.so"] = macho(dylibs=[OMP, libcxx], symbols=["_PyInit_core"])
        files["petakit5d/.dylibs/libomp.dylib"] = macho(
            ident="/DLC/petakit5d/.dylibs/libomp.dylib", dylibs=["/usr/lib/libSystem.B.dylib"]
        )
    else:
        files["petakit5d/_ext/core.pyd"] = windows_module()
        # Microsoft's redistributables, signed, naming Microsoft's own build folders
        omp = pe(imports=["KERNEL32.dll", "PSAPI.DLL"], extra=MS_PATH, signer=MS)
        files[f"petakit5d.libs/{OMP140}"] = omp
        files[f"petakit5d.libs/msvcp140-{HASH}.dll"] = pe(extra=MS_PATH, signer=MS)
    return files


def write_wheel(folder: Path, name: str, files: dict[str, bytes], metadata=None, tags=None,
                licenses=LICENSES) -> Path:  # fmt: skip
    dist, version, pys, abis, plats = name[: -len(".whl")].split("-")
    di = f"{dist}-{version}.dist-info"
    if tags is None:
        tags = [f"{p}-{a}-{q}" for p in pys.split(".") for a in abis.split(".")
                for q in plats.split(".")]  # fmt: skip
    expr = "GPL-2.0-or-later" if dist == "petakit5d_fftw" else "BSD-3-Clause"
    members = dict(files)
    wheel = "Wheel-Version: 1.0\nRoot-Is-Purelib: false\n" + "".join(f"Tag: {t}\n" for t in tags)
    members[f"{di}/WHEEL"] = wheel.encode()
    members[f"{di}/METADATA"] = metadata or _metadata(dist, version, expr, licenses)
    members.update({f"{di}/licenses/{f}": b"license text\n" for f in licenses})
    members[f"{di}/RECORD"] = "".join(f"{m},,\n" for m in members).encode()
    path = folder / name
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        for m, data in members.items():
            z.writestr(m, data)
    return path


def wheel_of(folder: Path, system: str, change=None, name=None, **kw) -> Path:
    files = petakit5d_files(system)
    if change:
        change(files)
    name = name or f"petakit5d-0.1.0-cp312-abi3-{PLATFORMS[system]}.whl"
    return write_wheel(folder, name, files, **kw)


def run(capsys, *args) -> tuple[int, str]:
    rc = cw.main([str(a) for a in args])
    out = capsys.readouterr()
    return rc, out.out + out.err


def put(path, data):
    return lambda files: files.__setitem__(path, data)


def drop(path):
    return lambda files: files.pop(path)


def rename(old, new):
    return lambda files: files.__setitem__(new, files.pop(old))


# --- the checks ---


@pytest.mark.parametrize("system", ["linux", "macos", "windows"])
def test_a_clean_wheel_passes(tmp_path, capsys, system):
    rc, out = run(capsys, wheel_of(tmp_path, system))
    assert rc == 0, out
    assert "ok, 0 findings" in out
    assert "1 archives, 0 with findings" in out


def test_the_inventory_names_the_libraries_and_search_paths(tmp_path, capsys):
    rc, out = run(capsys, wheel_of(tmp_path, "linux"))
    assert rc == 0, out
    assert f"library petakit5d/_ext/{GOMP}: " in out and f"SONAME {GOMP}" in out
    assert f"1 modules (core): needs {GOMP} libm.so.6 libc.so.6 | search path $ORIGIN" in out
    rc, out = run(capsys, "-q", wheel_of(tmp_path, "linux"))
    assert "library" not in out and "ok, 0 findings" in out


MOD = "petakit5d/_ext/core.abi3.so"


@pytest.mark.parametrize(
    "change, message",
    [
        pytest.param(
            put(MOD, linux_module(runpath="$ORIGIN:/opt/gcc/lib64")),
            r"RUNPATH '/opt/gcc/lib64' \(only \$ORIGIN\)",
            id="runpath",
        ),
        pytest.param(
            put(MOD, elf(needed=[GOMP], rpath="$ORIGIN")), r"RPATH '\$ORIGIN'", id="rpath"
        ),
        pytest.param(
            put(MOD, linux_module(needed=[GOMP, "libstdc++.so.6"])),
            r"needs libstdc\+\+\.so\.6, which is neither in the wheel",
            id="libstdc++",
        ),
        pytest.param(
            put(MOD, linux_module(needed=["libgomp.so.1"])),
            r"needs libgomp\.so\.1, which is neither",
            id="toolchain-libgomp",
        ),
        pytest.param(
            put(MOD, linux_module(runpath="/opt/gcc/lib64")),
            r"needs libgomp-e985bcbb\.so\.1\.0\.0, which is neither",
            id="runpath-elsewhere",
        ),
        pytest.param(drop(f"petakit5d/_ext/{GOMP}"), r"0 libraries match .*libgomp", id="no-omp"),
        pytest.param(
            put(
                "petakit5d/_ext/libgomp-0123abcd.so.1.0.0", elf(soname="libgomp-0123abcd.so.1.0.0")
            ),
            r"2 libraries match .*libgomp.* \(expected one\)",
            id="two-runtimes",
        ),
        pytest.param(
            put("petakit5d.libs/libgomp-0badc0de.so.1.0.0", elf(soname="x")),
            r"petakit5d\.libs/libgomp-0badc0de\.so\.1\.0\.0: a library the wheel should not carry",
            id="grafted",
        ),
        pytest.param(
            put("petakit5d/tmp/build-env-ab12/x.h", b""),
            r"x\.h \(name\): /petakit5d/tmp/build-env-ab12/x\.h \[temporary build environment\]",
            id="member-name",
        ),
        pytest.param(
            put(MOD, linux_module(symbols=["fftw_execute"])),
            r"defines FFTW symbols \(fftw_execute\)",
            id="fftw-symbol",
        ),
        pytest.param(
            put(MOD, linux_module(extra=b"fftw-3.3.8-sse2\0")),
            r"has FFTW's version string \(fftw-3\.3\.8\)",
            id="fftw-string",
        ),
        pytest.param(
            put(MOD, linux_module(extra=f"{HOME}src/cpp/x.cpp\0".encode())),
            re.escape(f"{HOME}src/cpp/x.cpp [POSIX home folder]"),
            id="home",
        ),
        pytest.param(
            put(MOD, linux_module(debug=b"/tmp/build-env-abc123/include/x.h\0")),
            r"core\.abi3\.so\[\.debug_str\]: /tmp/build-env-abc123/include/x\.h",
            id="compressed",
        ),
        pytest.param(
            put(MOD, linux_module(extra=b"/tmp/tmpq7x_0abc/build/_deps/a.c\0")),
            r"\[Python temporary folder\]",
            id="tempdir",
        ),
        pytest.param(
            put("petakit5d/x.py", b"# /project/python/cpp/core/a.cpp\n"),
            r"\[cibuildwheel's project folder\]",
            id="container",
        ),
        pytest.param(
            put("petakit5d/__pycache__/x.cpython-312.pyc", b"\0"),
            r"compiled Python in the wheel",
            id="pyc",
        ),
        pytest.param(
            put(MOD, linux_module(machine=183)), r"built for arm64, not x86_64", id="machine"
        ),
        pytest.param(put(MOD, macho()), r"a macho file in a linux wheel", id="macho"),
        pytest.param(
            rename(MOD, "petakit5d/_ext/core.cpython-312-x86_64-linux-gnu.so"),
            r"the module's name does not fit the ABI tag abi3",
            id="module-name",
        ),
        pytest.param(drop(MOD), r"no extension modules", id="no-modules"),
        pytest.param(put(MOD, linux_module()[:100]), r"can't be read as elf", id="truncated"),
    ],
)
def test_linux_findings(tmp_path, capsys, change, message):
    rc, out = run(capsys, wheel_of(tmp_path, "linux", change))
    assert rc == 1, out
    assert re.search(message, out), out


@pytest.mark.parametrize(
    "change, message",
    [
        pytest.param(
            put(MOD, macho(dylibs=[OMP], rpaths=["/opt/homebrew/opt/libomp/lib"])),
            r"LC_RPATH '/opt/homebrew/opt/libomp/lib' \(only @loader_path\)",
            id="rpath",
        ),
        pytest.param(
            put(MOD, macho(dylibs=["/opt/homebrew/opt/libomp/lib/libomp.dylib"])),
            r"loads /opt/homebrew/opt/libomp/lib/libomp\.dylib, which is neither",
            id="absolute",
        ),
        pytest.param(
            put(MOD, macho(dylibs=["@rpath/libomp.dylib"], rpaths=["@loader_path"])),
            r"loads @rpath/libomp\.dylib, which is neither",
            id="rpath-misses",
        ),
        pytest.param(
            # the executable's folder is not the module's
            put(MOD, macho(dylibs=["@rpath/libomp.dylib"], rpaths=["@executable_path/../.dylibs"])),
            r"loads @rpath/libomp\.dylib, which is neither",
            id="executable-path",
        ),
        pytest.param(
            put("petakit5d/.dylibs/libomp.dylib", macho(ident="/Users/runner/work/o/libomp.dylib")),
            r"install name '/Users/runner/work/o/libomp\.dylib' is a folder outside",
            id="install-name",
        ),
        pytest.param(
            put(MOD, macho(dylibs=[OMP], minos=(15, 0, 0))),
            r"deployment target 15\.0 is newer than 13\.0",
            id="target",
        ),
        pytest.param(
            put(MOD, fat((0x01000007, macho(cpu=0x01000007)), (0x0100000C, macho()))),
            r"built for x86_64, arm64, not arm64",
            id="universal",
        ),
        pytest.param(
            put(MOD, macho(dylibs=[OMP], symbols=["_fftwf_plan_dft_1d"])),
            r"defines FFTW symbols \(_fftwf_plan_dft_1d\)",
            id="fftw-symbol",
        ),
        pytest.param(
            put(MOD, macho(dylibs=[OMP], extra=b"/private/var/folders/x1/abc/T/a.o\0")),
            r"\[macOS per-user temporary folder\]",
            id="tempdir",
        ),
        pytest.param(
            put(MOD, macho(dylibs=[OMP], extra=b"/opt/homebrew/Cellar/libomp/21/x\0")),
            r"\[Homebrew prefix\]",
            id="homebrew",
        ),
        pytest.param(
            drop("petakit5d/.dylibs/libomp.dylib"),
            r"loads @loader_path/\.\./\.dylibs/libomp\.dylib, which is neither",
            id="no-omp",
        ),
    ],
)
def test_macos_findings(tmp_path, capsys, change, message):
    rc, out = run(capsys, wheel_of(tmp_path, "macos", change))
    assert rc == 1, out
    assert re.search(message, out), out


def test_macos_resolves_through_a_relative_rpath(tmp_path, capsys):
    mod = macho(dylibs=["@rpath/libomp.dylib"], rpaths=["@loader_path/../.dylibs"])
    rc, out = run(capsys, wheel_of(tmp_path, "macos", put(MOD, mod)))
    assert rc == 0, out


WIN = "petakit5d/_ext/core.pyd"


@pytest.mark.parametrize(
    "change, message",
    [
        pytest.param(
            put(WIN, windows_module("python312.dll")),
            r"loads python312\.dll, not python3\.dll",
            id="python-dll",
        ),
        pytest.param(
            put(WIN, windows_module("python3.dll", "libgomp-1.dll")),
            r"loads libgomp-1\.dll, which is neither in the wheel",
            id="missing",
        ),
        pytest.param(
            put(WIN, windows_module(delay=["cudart64_12.dll"])),
            r"loads cudart64_12\.dll, which is neither",
            id="delayed",
        ),
        pytest.param(
            put(
                WIN, windows_module(extra=b"C:\\Users\\builder\\AppData\\Local\\Temp\\pip-x\\a.h\0")
            ),
            r"\[Windows user folder, temporary build environment\]",
            id="user",
        ),
        pytest.param(
            put(WIN, windows_module(extra=b"D:\\a\\PetaKit5D\\PetaKit5D\\python\\build\\c.pdb\0")),
            r"\[CI runner workspace\]",
            id="workspace",
        ),
        pytest.param(
            put(WIN, pe(imports=["python3.dll"], exports=["fftw_execute"])),
            r"defines FFTW symbols \(fftw_execute\)",
            id="fftw-export",
        ),
        pytest.param(
            put(f"petakit5d.libs/libwinpthread-1-{HASH}.dll", pe()),
            r"libwinpthread-1-[0-9a-f]+\.dll: a library the wheel should not carry",
            id="unexpected",
        ),
        pytest.param(
            put(WIN, windows_module(machine=0xAA64)), r"built for arm64, not x86_64", id="machine"
        ),
        pytest.param(
            rename(WIN, "petakit5d/_ext/core.cp312-win_amd64.pyd"),
            r"the module's name does not fit the ABI tag abi3",
            id="module-name",
        ),
        pytest.param(
            rename(f"petakit5d.libs/{OMP140}", f"petakit5d/x/{OMP140}"),
            r"loads libomp140\.x86_64-[0-9a-f]+\.dll, which is neither",
            id="dll-elsewhere",
        ),
        pytest.param(
            put("petakit5d/_ext/libomp140.x86_64.dll", pe(imports=["KERNEL32.dll"])),
            r"2 libraries match .*libomp140.* \(expected one\)",
            id="two-runtime-files",
        ),
        pytest.param(
            put("petakit5d/_ext/io.pyd", windows_module(omp="vcomp140.dll")),
            r"the modules load 2 OpenMP runtimes: libomp140\.x86_64-[0-9a-f]+\.dll, vcomp140\.dll",
            id="two-runtimes",
        ),
        # Microsoft's build folders count anywhere but in its signed redistributables
        pytest.param(
            put(f"petakit5d.libs/msvcp140-{HASH}.dll", pe(extra=MS_PATH)),
            r"msvcp140-[0-9a-f]+\.dll: D:\\a\\_work\\.* \[CI runner workspace\]",
            id="unsigned-redistributable",
        ),
        pytest.param(
            put(f"petakit5d.libs/msvcp140-{HASH}.dll", pe(extra=MS_PATH, signer="Example Corp")),
            r"msvcp140-[0-9a-f]+\.dll: D:\\a\\_work\\.* \[CI runner workspace\]",
            id="other-signer",
        ),
        pytest.param(
            put(WIN, windows_module(extra=MS_PATH)),
            r"core\.pyd: D:\\a\\_work\\.* \[CI runner workspace\]",
            id="module-with-azure-folder",
        ),
    ],
)
def test_windows_findings(tmp_path, capsys, change, message):
    rc, out = run(capsys, wheel_of(tmp_path, "windows", change))
    assert rc == 1, out
    assert re.search(message, out), out


def test_windows_runtime_next_to_the_modules(tmp_path, capsys):
    # a build that copies libomp140 next to the modules, and delvewheel that leaves it there
    def copied(files):
        files.pop(f"petakit5d.libs/{OMP140}")
        files["petakit5d/_ext/libomp140.x86_64.dll"] = pe(imports=["KERNEL32.dll"])
        files[WIN] = windows_module(omp="libomp140.x86_64.dll")

    rc, out = run(capsys, wheel_of(tmp_path, "windows", copied))
    assert rc == 0, out


# a wheel of modules built with MinGW-w64's GCC: its libgomp-1.dll and the DLLs that loads, next to
# the modules under their own names (no delvewheel copies, no Microsoft runtime)
MINGW_DLLS = {
    "petakit5d/_ext/libgomp-1.dll": [
        "libgcc_s_seh-1.dll",
        "KERNEL32.dll",
        "api-ms-win-crt-runtime-l1-1-0.dll",
        "libwinpthread-1.dll",
    ],
    "petakit5d/_ext/libwinpthread-1.dll": ["KERNEL32.dll", "api-ms-win-crt-runtime-l1-1-0.dll"],
    "petakit5d/_ext/libgcc_s_seh-1.dll": ["KERNEL32.dll", "libwinpthread-1.dll"],
}


GNU_LD = (2, 39)
# the toolchain's notices that a MinGW-w64 build copies in (cmake/Pk5dExtension.cmake)
NOTICES = [
    "petakit5d/_ext/licenses/gcc/COPYING3",
    "petakit5d/_ext/licenses/gcc/COPYING.RUNTIME",
    "petakit5d/_ext/licenses/mingw-w64/COPYING.MinGW-w64-runtime.txt",
    "petakit5d/_ext/licenses/winpthreads/COPYING",
]


def mingw_module(*imports, omp="libgomp-1.dll", **kw) -> bytes:
    std = ["KERNEL32.dll", "api-ms-win-crt-math-l1-1-0.dll", "libwinpthread-1.dll", omp]
    kw.setdefault("linker", GNU_LD)
    return pe(imports=[*(imports or ["python3.dll"]), *std], exports=["PyInit_core"], **kw)


def mingw(*changes):
    """the MinGW-w64 wheel, then the changes"""

    def change(files):
        for name in [n for n in files if n.startswith("petakit5d.libs/")]:
            files.pop(name)
        files[WIN] = mingw_module()
        for name, imports in MINGW_DLLS.items():
            files[name] = pe(imports=imports, linker=GNU_LD)
        for n in NOTICES:
            files[n] = b"notice text, as the toolchain ships it\n"
        for c in changes:
            c(files)

    return change


def test_a_mingw_wheel_passes(tmp_path, capsys):
    rc, out = run(capsys, wheel_of(tmp_path, "windows", mingw()))
    assert rc == 0, out
    assert "OpenMP runtime of the modules: libgomp-1.dll" in out
    assert "library petakit5d/_ext/libgcc_s_seh-1.dll" in out
    assert "linked by GNU ld (MinGW-w64): core.pyd libgcc_s_seh-1.dll" in out
    assert f"{NOTICES[2]}: the toolchain's notice, 39 bytes, sha256" in out


def test_the_mingw_notices_are_left_out_of_the_search_by_their_paths(tmp_path, capsys):
    # the toolchain's runtime notice cites a URL that a site's pattern of node names can match:
    # here a made-up site pattern and text
    url = b"<http://lists.example.org/legal/2004/12/rack07.html>\n"
    forbid = r"rack[0-9]{2}"
    path = wheel_of(tmp_path, "windows", mingw(put(NOTICES[2], url)))
    rc, out = run(capsys, "--forbid", forbid, path)
    assert rc == 0, out
    assert f"{NOTICES[2]}: a notice of MinGW-w64's toolchain, as it ships it (not searched)" in out
    # anywhere else it is searched: another file beside them, a module, a notice's name elsewhere
    for where in ("petakit5d/_ext/licenses/README.txt", "petakit5d/_ext/licenses/gcc/COPYING",
                  "petakit5d_extra/licenses/mingw-w64/COPYING.MinGW-w64-runtime.txt"):  # fmt: skip
        path = wheel_of(tmp_path, "windows", mingw(put(where, url)))
        rc, out = run(capsys, "--forbid", forbid, path)
        assert rc == 1 and f"{where}: " in out and "[--forbid 1]" in out, (where, out)
    path = wheel_of(tmp_path, "windows", mingw(put(WIN, mingw_module(extra=url))))
    rc, out = run(capsys, "--forbid", forbid, path)
    assert rc == 1 and "core.pyd: " in out and "[--forbid 1]" in out, out
    # the folders of a build machine too: in a notice's path they are not searched for either,
    # but in a binary at that path they are
    user = b"C:\\Users\\builder\\x.c\0"
    rc, out = run(capsys, wheel_of(tmp_path, "windows", mingw(put(NOTICES[0], pe(extra=user)))))
    assert rc == 1 and "[Windows user folder]" in out, out


@pytest.mark.parametrize(
    "change, message",
    [
        pytest.param(drop(NOTICES[3]), rf"{NOTICES[3]}: missing", id="missing"),
        pytest.param(put(NOTICES[1], b"text\0more"), rf"{NOTICES[1]}: not a text file", id="nul"),
    ],
)
def test_mingw_notice_findings(tmp_path, capsys, change, message):
    rc, out = run(capsys, wheel_of(tmp_path, "windows", mingw(change)))
    assert rc == 1, out
    assert re.search(message, out), out


def test_an_msvc_wheel_needs_no_mingw_notices(tmp_path, capsys):
    rc, out = run(capsys, wheel_of(tmp_path, "windows"))
    assert rc == 0, out
    assert "linked by GNU ld" not in out and "notice" not in out


@pytest.mark.parametrize(
    "change, message",
    [
        # delvewheel's renamed copy would be a runtime of its own, beside the I/O packages'
        pytest.param(
            mingw(
                rename("petakit5d/_ext/libgomp-1.dll", f"petakit5d.libs/libgomp-1-{HASH}.dll"),
                put(WIN, mingw_module(omp=f"libgomp-1-{HASH}.dll")),
            ),
            r"libgomp-1-[0-9a-f]+\.dll: a library the wheel should not carry",
            id="renamed",
        ),
        pytest.param(
            mingw(drop("petakit5d/_ext/libgcc_s_seh-1.dll")),
            r"libgomp-1\.dll: loads libgcc_s_seh-1\.dll, which is neither in the wheel",
            id="libgcc-missing",
        ),
        pytest.param(
            mingw(drop("petakit5d/_ext/libwinpthread-1.dll")),
            r"0 libraries match petakit5d/_ext/libwinpthread-1\\\.dll \(expected one\)",
            id="winpthreads-missing",
        ),
        pytest.param(
            mingw(put("petakit5d/_ext/libstdc++-6.dll", pe(imports=["KERNEL32.dll"]))),
            r"libstdc\+\+-6\.dll: a library the wheel should not carry",
            id="libstdc++",
        ),
        pytest.param(
            mingw(put(WIN, mingw_module("python3.dll", "libstdc++-6.dll"))),
            r"loads libstdc\+\+-6\.dll, which is neither in the wheel",
            id="libstdc++-dynamic",
        ),
        pytest.param(
            mingw(put("petakit5d/_ext/io.pyd", windows_module())),
            r"the modules load 2 OpenMP runtimes",
            id="two-toolchains",
        ),
        pytest.param(
            mingw(put("petakit5d/_ext/libgomp-1.dll", pe(extra=b"C:\\Users\\builder\\gcc\\x.c\0"))),
            r"libgomp-1\.dll: C:\\Users\\builder.* \[Windows user folder\]",
            id="build-folder",
        ),
    ],
)
def test_mingw_findings(tmp_path, capsys, change, message):
    rc, out = run(capsys, wheel_of(tmp_path, "windows", change))
    assert rc == 1, out
    assert re.search(message, out), out


def test_signed_redistributables_get_the_forbid_patterns_only(tmp_path, capsys):
    rc, out = run(capsys, wheel_of(tmp_path, "windows"))
    assert rc == 0, out
    assert f"petakit5d.libs/{OMP140}: Microsoft's signed redistributable" in out
    assert f"OpenMP runtime of the modules: {OMP140}" in out
    site = pe(extra=MS_PATH + b"built on node42\0", signer=MS)
    path = wheel_of(tmp_path, "windows", put(f"petakit5d.libs/msvcp140-{HASH}.dll", site))
    rc, out = run(capsys, "--forbid", r"node\d+", path)
    assert rc == 1 and "built on node42 [--forbid 1]" in out and "CI runner workspace" not in out


def test_readers_on_signatures():
    assert cw.parse_pe(pe(signer=MS)).signer == MS
    assert cw.parse_pe(pe(signer="Example Corp")).signer == "another publisher"
    assert cw.parse_pe(pe()).signer == ""


def test_a_cp311_wheel_has_its_versions_modules(tmp_path, capsys):
    name = "petakit5d-0.1.0-cp311-cp311-win_amd64.whl"
    cp311 = "petakit5d/_ext/core.cp311-win_amd64.pyd"
    ok = put(cp311, windows_module("python311.dll"))
    rc, out = run(capsys, wheel_of(tmp_path, "windows", lambda f: (f.pop(WIN), ok(f)), name=name))
    assert rc == 0, out
    stable = put(cp311, windows_module("python3.dll"))
    rc, out = run(
        capsys, wheel_of(tmp_path, "windows", lambda f: (f.pop(WIN), stable(f)), name=name)
    )
    assert rc == 1 and "loads python3.dll, not python311.dll" in out
    linux = rename(MOD, "petakit5d/_ext/core.cpython-311-x86_64-linux-gnu.so")
    name = "petakit5d-0.1.0-cp311-cp311-manylinux_2_28_x86_64.whl"
    rc, out = run(capsys, wheel_of(tmp_path, "linux", linux, name=name))
    assert rc == 0, out


def test_tags_and_license(tmp_path, capsys):
    rc, out = run(capsys, wheel_of(tmp_path, "linux", tags=["cp312-abi3-manylinux_2_28_x86_64"]))
    assert rc == 1 and "WHEEL's tags" in out
    meta = _metadata("petakit5d", "0.1.0", "MIT", LICENSES)
    rc, out = run(capsys, wheel_of(tmp_path, "linux", metadata=meta))
    assert rc == 1 and "License-Expression is 'MIT', not 'BSD-3-Clause'" in out
    rc, out = run(capsys, wheel_of(tmp_path, "linux", licenses=LICENSES[:2]))
    assert rc == 1
    assert "the license file cpp/fft/third_party/pocketfft/LICENSE.md is not declared" in out
    assert "is not in petakit5d-0.1.0.dist-info/licenses/" in out


def test_forbidden_patterns(tmp_path, capsys):
    path = wheel_of(tmp_path, "linux", put("petakit5d/site.py", b"# built on node42.example.org\n"))
    assert run(capsys, path)[0] == 0
    rc, out = run(capsys, "--forbid", r"node\d+\.example", path)
    assert rc == 1 and "[--forbid 1]" in out
    patterns = tmp_path / "site.txt"
    patterns.write_text("# the site\n\nexample\\.org\n", newline="\n")
    rc, out = run(capsys, "--forbid-file", patterns, path)
    assert rc == 1 and "petakit5d/site.py: # built on node42.example.org [--forbid 1]" in out
    assert run(capsys, "--forbid", "(", path)[0] == 2


ADDON = {
    "linux": ("manylinux_2_28_x86_64", "libpk5d_fftw3{}.so"),
    "macos": ("macosx_11_0_arm64", "libpk5d_fftw3{}.dylib"),
    "windows": ("win_amd64", "pk5d_fftw3{}.dll"),
}


def addon_files(system: str) -> dict[str, bytes]:
    def lib(symbol):
        if system == "linux":
            return elf(needed=["libm.so.6", "libc.so.6"], symbols=[symbol], extra=b"fftw-3.3.8\0")
        if system == "macos":
            return macho(
                dylibs=["/usr/lib/libSystem.B.dylib"], symbols=["_" + symbol], minos=(11, 0, 0)
            )
        return pe(imports=["KERNEL32.dll", "VCRUNTIME140.dll"], exports=[symbol])

    name = "petakit5d_fftw/" + ADDON[system][1]
    return {
        "petakit5d_fftw/__init__.py": b"",
        "petakit5d_fftw/FFTW-COPYING": b"GPL\n",
        name.format("f"): lib("fftwf_execute"),
        name.format(""): lib("fftw_execute"),
    }


@pytest.mark.parametrize("system", ["linux", "macos", "windows"])
def test_the_addon(tmp_path, capsys, system):
    name = f"petakit5d_fftw-3.3.8.0-py3-none-{ADDON[system][0]}.whl"
    rc, out = run(capsys, write_wheel(tmp_path, name, addon_files(system), licenses=("COPYING",)))
    assert rc == 0, out
    files = addon_files(system)
    files["petakit5d_fftw/extra.so"] = elf(needed=["libc.so.6"])
    files["petakit5d_fftw/m.abi3.so"] = elf()
    rc, out = run(capsys, write_wheel(tmp_path, name, files, licenses=()))
    assert rc == 1
    assert "petakit5d_fftw/extra.so: a library the wheel should not carry" in out
    assert "extension modules in a library wheel: petakit5d_fftw/m.abi3.so" in out
    assert "the license file COPYING is not declared" in out


def test_the_addon_built_with_mingw(tmp_path, capsys):
    # libraries that GNU ld linked (libgcc and MinGW-w64's runtime linked in) carry the
    # toolchain's notices of those, under petakit5d_fftw/licenses/
    name = "petakit5d_fftw-3.3.8.0-py3-none-win_amd64.whl"
    files = addon_files("windows")
    for v in ("f", ""):
        files[f"petakit5d_fftw/pk5d_fftw3{v}.dll"] = pe(
            imports=["KERNEL32.dll", "api-ms-win-crt-math-l1-1-0.dll"],
            exports=[f"fftw{v}_execute"],
            linker=GNU_LD,
        )
    rc, out = run(capsys, write_wheel(tmp_path, name, files, licenses=("COPYING",)))
    assert rc == 1
    for n in ("gcc/COPYING3", "gcc/COPYING.RUNTIME", "mingw-w64/COPYING.MinGW-w64-runtime.txt"):
        assert f"petakit5d_fftw/licenses/{n}: missing" in out
        files[f"petakit5d_fftw/licenses/{n}"] = b"rack07\n"
    assert "winpthreads" not in out
    path = write_wheel(tmp_path, name, files, licenses=("COPYING",))
    rc, out = run(capsys, "--forbid", r"rack[0-9]{2}", path)
    assert rc == 0, out


def write_sdist(folder: Path, files: dict[str, bytes], owner: str = "") -> Path:
    path = folder / "petakit5d-0.1.0.tar.gz"
    members = {
        "PKG-INFO": _metadata("petakit5d", "0.1.0", "BSD-3-Clause", LICENSES),
        "pyproject.toml": b"[project]\n",
        "src/petakit5d/__init__.py": b"",
        **dict.fromkeys(LICENSES, b"license text\n"),
        **files,
    }
    with tarfile.open(path, "w:gz") as t:
        for m, data in members.items():
            info = tarfile.TarInfo("petakit5d-0.1.0/" + m)
            info.size, info.uname, info.gname = len(data), owner, owner
            t.addfile(info, io.BytesIO(data))
    return path


def test_the_sdist(tmp_path, capsys):
    # a source may name a build machine's folders (a search hint): only --forbid finds them there
    hint = {"cmake/omp.cmake": b"find_path(OMP omp.h HINTS /opt/homebrew/opt/libomp)\n"}
    rc, out = run(capsys, write_sdist(tmp_path, hint))
    assert rc == 0, out
    rc, out = run(capsys, "--forbid", "/opt/homebrew", write_sdist(tmp_path, hint))
    assert rc == 1 and "cmake/omp.cmake: find_path(OMP" in out and "[--forbid 1]" in out
    bad = {
        "cpp/fft/exact/CMakeLists.txt": b"",
        "build/cp312/x.txt": b"",
        "src/petakit5d/__pycache__/x.cpython-312.pyc": b"",
        "src/petakit5d/_ext/core.abi3.so": elf(),
        "tests/data/blob.bin": macho(),
        "tmp/build-env-ab12/x.h": b"",
    }
    rc, out = run(capsys, write_sdist(tmp_path, bad, owner="builder"))
    assert rc == 1
    for text in (
        "cpp/fft/exact/CMakeLists.txt: the exact FFT add-on",
        "build/cp312/x.txt: a build folder",
        "x.cpython-312.pyc: compiled Python",
        "core.abi3.so: a binary",
        "tests/data/blob.bin: a binary (ELF, Mach-O or PE)",
        "x.h (name): /petakit5d-0.1.0/tmp/build-env-ab12/x.h [temporary build environment]",
        "owned by the name 'builder'",
    ):
        assert text in out, (text, out)
    meta = _metadata("petakit5d", "0.1.0", "BSD-3-Clause", LICENSES) + b"Home-page: /Users/ci/x/\n"
    rc, out = run(capsys, write_sdist(tmp_path, {"PKG-INFO": meta}))
    assert rc == 1 and "PKG-INFO: Home-page: /Users/ci/x/ [macOS user folder]" in out


def test_folders_and_unreadable_archives(tmp_path, capsys):
    wheel_of(tmp_path, "linux")
    write_sdist(tmp_path, {})
    rc, out = run(capsys, tmp_path)
    assert rc == 0 and "2 archives, 0 with findings" in out
    other = tmp_path / "other"
    other.mkdir()
    assert run(capsys, other)[0] == 2
    assert run(capsys, tmp_path / "missing.whl")[0] == 2
    broken = other / "petakit5d-0.1.0-cp312-abi3-win_amd64.whl"
    broken.write_bytes(b"not a zip")
    assert run(capsys, broken)[0] == 2
    rc, out = run(capsys, write_wheel(other, "numpy-2.0-cp312-cp312-win_amd64.whl", {}))
    assert rc == 2 and "no profile for the distribution numpy" in out


def test_a_string_matching_several_patterns_is_one_finding():
    pats = [(label, re.compile(rx)) for label, rx in cw.PATH_PATTERNS]
    data = (
        b"\0C:\\Users\\builder\\AppData\\Local\\Temp\\build-env-12\\x.h\0" + HOME.encode() + b"y\0"
    )
    shown, total = cw.scan_bytes(data, pats)
    assert total == 2
    assert shown[0].endswith("[Windows user folder, temporary build environment]")
    # placeholders are no one's folders
    placeholders = HOME.replace("builder", "<user>") + "x /Users/username/y C:\\Users\\name\\z"
    assert cw.scan_bytes(placeholders.encode(), pats) == ([], 0)


def test_readers_on_this_builds_modules():
    # this platform's installed (or built) modules, read as the wheels' are
    from petakit5d import _ext

    if not _ext.available("core"):
        pytest.skip("the extension modules aren't built")
    kind = {"linux": "elf", "darwin": "macho", "win32": "pe"}.get(sys.platform)
    if kind is None:
        pytest.skip(f"no reader for {sys.platform}")
    b = cw.parse_binary(Path(_ext.load("core").__file__).read_bytes())
    assert b is not None and b.kind == kind and not b.problems
    machine = platform.machine().lower()
    assert b.machines == [{"amd64": "x86_64", "aarch64": "arm64"}.get(machine, machine)]
    assert any(s.lstrip("_") == "PyInit_core" for s in b.defined)
    if kind == "pe":
        assert any(re.fullmatch(r"python3\d*\.dll", n.lower()) for n in b.needed), b.needed
    else:
        assert any("libc" in n or "libSystem" in n for n in b.needed), b.needed
