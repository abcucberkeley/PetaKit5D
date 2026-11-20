"""The release checks of petakit5d's wheels and source distribution (docs/release.md).

    python tools/check_wheel.py [--forbid REGEX]... [--forbid-file FILE] [-q] ARCHIVE|FOLDER...

A folder means its *.whl and *.tar.gz files. The wheels' binaries (ELF, Mach-O, PE) are read
with the standard library alone, so the checks are the same on every platform and a wheel of any
platform can be checked on any machine. The distribution is taken from the file name: petakit5d,
or petakit5d_fftw (the exact FFT add-on, GPL). For each wheel:

- the tags: the WHEEL file has the file name's tags; every binary is of the platform's format
  and architecture; the extension modules are built for the tags (abi3: `<group>.abi3.so`, or a
  `.pyd` that loads python3.dll; cp3XX: that version's suffix and DLL); a Mach-O file's
  deployment target is not newer than the wheel's;
- the bundled libraries: the libraries besides the modules are the expected ones (petakit5d's
  OpenMP runtime: on Linux the libgomp it shares with cpp-zarr, next to the modules; on macOS
  libomp in .dylibs; on Windows libomp140 next to the modules, and the C++ runtime delvewheel
  adds; the add-on's two FFTW libraries), and every library a binary loads is in the wheel or
  the system's;
- the search paths: RUNPATH (Linux) only $ORIGIN; LC_RPATH (macOS) only @loader_path; no install
  name outside the wheel and the system's folders;
- no folder of the build machine named in any member (cpp-tiff's check_wheel_paths.py approach:
  every member's bytes, with ELF's compressed sections decompressed, and its name, against home
  folders, temporary build environments, cibuildwheel's and the CI runners' folders and the
  --forbid patterns), and no compiled Python file;
- no FFTW (GPL) in petakit5d's binaries: no fftw symbol defined, no FFTW version string;
- the license: METADATA's License-Expression, and every License-File (the project's own among
  them) in .dist-info/licenses.

For the source distribution: no cpp/fft/exact/ (the add-on), no build folders, compiled files or
binaries, no tar owner names, the license, and the folder check of the member names and PKG-INFO
(the other files are the repository's sources, where a build folder's name can be legitimate:
their contents get the --forbid patterns only).

Exit status: 0 if every archive is clean, 1 on any finding, 2 on bad arguments or an archive that
can't be read.
"""

from __future__ import annotations

import argparse
import email.parser
import os
import posixpath
import re
import struct
import sys
import tarfile
import zipfile
import zlib
from dataclasses import dataclass, field
from pathlib import PurePosixPath

# --- what each distribution holds ---


@dataclass(frozen=True)
class Profile:
    license: str
    license_files: tuple[str, ...]
    # per system, the libraries besides the modules: (pattern of the path, required)
    libraries: dict[str, tuple[tuple[str, bool], ...]]
    modules: bool
    fftw: bool


PROFILES = {
    "petakit5d": Profile(
        license="BSD-3-Clause",
        license_files=(
            "LICENSE",
            "THIRD_PARTY_NOTICES.md",
            "cpp/fft/third_party/pocketfft/LICENSE.md",
        ),
        libraries={
            # the libgomp of the cpp-zarr and cpp-tiff wheels, as auditwheel named it there
            # (docs/decisions/0008-openmp-runtime.md)
            "linux": ((r"petakit5d/_ext/libgomp-[0-9a-f]{8}\.so\.1\.0\.0", True),),
            "macos": ((r"petakit5d/\.dylibs/libomp\.dylib", True),),
            "windows": (
                (r"petakit5d/_ext/libomp140\.x86_64\.dll", True),
                # delvewheel's copies of MSVC's C++ runtime, with mangled names
                (r"petakit5d\.libs/msvcp140(?:_[0-9a-z_]+)?-[0-9a-f]{32}\.dll", False),
            ),
        },
        modules=True,
        fftw=False,
    ),
    "petakit5d_fftw": Profile(
        license="GPL-2.0-or-later",
        license_files=("COPYING",),
        libraries={
            "linux": (
                (r"petakit5d_fftw/libpk5d_fftw3f\.so", True),
                (r"petakit5d_fftw/libpk5d_fftw3\.so", True),
            ),
            "macos": (
                (r"petakit5d_fftw/libpk5d_fftw3f\.dylib", True),
                (r"petakit5d_fftw/libpk5d_fftw3\.dylib", True),
            ),
            "windows": (
                (r"petakit5d_fftw/pk5d_fftw3f\.dll", True),
                (r"petakit5d_fftw/pk5d_fftw3\.dll", True),
            ),
        },
        modules=False,
        fftw=True,
    ),
}

# libraries of the system, which a wheel does not carry
ELF_SYSTEM = {
    "libc.so.6",
    "libm.so.6",
    "libpthread.so.0",
    "libdl.so.2",
    "librt.so.1",
    "ld-linux-x86-64.so.2",
    "ld-linux-aarch64.so.1",
}
MACHO_SYSTEM = ("/usr/lib/", "/System/Library/")
PE_SYSTEM = {
    "kernel32.dll",
    "user32.dll",
    "advapi32.dll",
    "ntdll.dll",
    "ws2_32.dll",
    "bcrypt.dll",
    "ole32.dll",
    "oleaut32.dll",
    "shell32.dll",
    "shlwapi.dll",
    "psapi.dll",
    "dbghelp.dll",
    "version.dll",
    "winmm.dll",
    "ucrtbase.dll",
    # MSVC's C runtime, which every CPython for Windows ships
    "vcruntime140.dll",
    "vcruntime140_1.dll",
}
PE_SYSTEM_PREFIXES = ("api-ms-win-", "ext-ms-win-")

# FFTW's symbols (also with Mach-O's leading underscore) and its version string ("fftw-3.3.8")
FFTW_SYMBOL = re.compile(r"_?[ds]?fftw[flq]?_\w")
FFTW_STRING = re.compile(rb"fftw-3\.\d+\.\d+")

# --- folders of a build machine (the approach of cpp-tiff's check_wheel_paths.py) ---

# A name in a path: no separator, NUL, control byte or Windows-reserved character (so that a
# placeholder such as <user> is no name); POSIX names here end at a space
_WIN_NAME = rb'[^\\/\x00-\x1f"<>|:*?]+'
_POSIX_NAME = rb"[A-Za-z0-9_.][A-Za-z0-9_.-]*"
# a placeholder of a docstring (/Users/username/) is no user's folder
_NOT_PLACEHOLDER = rb"(?!(?i:user|username|yourname|your_?username|name)[\\/])"
_WIN_USER = _NOT_PLACEHOLDER + _WIN_NAME
_POSIX_USER = _NOT_PLACEHOLDER + _POSIX_NAME
# the start of an absolute POSIX path, not the middle of one; the folder of home folders
_ROOT = rb"(?<![A-Za-z0-9_.-])/"
_HOMES = rb"home/"
PATH_PATTERNS = [
    ("Windows user folder", rb"(?i)[A-Z]:[\\/]+Users[\\/]+" + _WIN_USER + rb"[\\/]"),
    ("POSIX home folder", _ROOT + _HOMES + _POSIX_USER + rb"/"),
    (
        "POSIX home folder under a site root",
        _ROOT + _POSIX_NAME + rb"/" + _HOMES + _POSIX_USER + rb"/",
    ),
    ("macOS user folder", rb"(?<![A-Za-z0-9_.:-])/Users/" + _POSIX_USER + rb"/"),
    ("temporary build environment", rb"/tmp/+(?:build-env-|pip-)"),
    ("temporary build environment", rb"(?i)[\\/]Temp[\\/]+(?:build-env-|pip-)"),
    ("Python temporary folder", _ROOT + rb"tmp/+tmp[a-z0-9_]{8}/"),
    ("Python temporary folder", rb"(?i)[\\/]Temp[\\/]+tmp[a-z0-9_]{8}[\\/]"),
    ("macOS per-user temporary folder", rb"/private/var/folders/"),
    ("macOS per-user temporary folder", rb"/var/folders/[A-Za-z0-9_+-]{2}/[A-Za-z0-9_+-]{2,}/"),
    ("cibuildwheel folder", rb"cibw-run-|/tmp/cibuildwheel/"),
    ("cibuildwheel cache", rb"(?i)pypa[\\/]+cibuildwheel[\\/]+Cache"),
    # the project inside cibuildwheel's Linux container, with python/ or as the root
    ("cibuildwheel's project folder", _ROOT + rb"project/(?:python/)?(?:src|cpp|cmake|build)/"),
    # GitHub's runners: the workspace (D:\a\<repo>\<repo> on Windows) and the tool cache
    ("CI runner workspace", rb"(?i)(?<![A-Za-z0-9])[A-Z]:[\\/]+a[\\/]+[A-Za-z0-9_.-]+[\\/]"),
    ("CI runner tool cache", rb"(?i)hostedtoolcache"),
    # the build machine's Homebrew, which delocate's copies must not name
    ("Homebrew prefix", rb"/opt/homebrew/|/usr/local/(?:Cellar|opt)/|@@HOMEBREW_"),
    ("Spack build stage", rb"spack-stage-|[\\/]spack-(?:temporary|stage)[\\/]"),
]
MAX_PER_MEMBER = 5
SHOW = 120


def _printable(c: int) -> bool:
    return 0x20 <= c < 0x7F


def scan_bytes(data: bytes, patterns: list[tuple[str, re.Pattern]]) -> tuple[list[str], int]:
    """The first MAX_PER_MEMBER printable strings with matches (with their pattern names), and
    the number of such strings (a string that several patterns match is one finding)."""
    found: dict[int, tuple[int, int, list[str]]] = {}
    for label, rx in patterns:
        for m in rx.finditer(data):
            lo = m.start()
            while lo > 0 and m.start() - lo < 4096 and _printable(data[lo - 1]):
                lo -= 1
            start, end, labels = found.setdefault(lo, (m.start(), m.end(), []))
            if label not in labels:
                labels.append(label)
            if m.start() < start:
                found[lo] = (m.start(), m.end(), labels)
    shown = []
    for lo in sorted(found)[:MAX_PER_MEMBER]:
        start, end, labels = found[lo]
        lo = max(lo, start - 40)
        hi = end
        while hi < len(data) and hi - lo < SHOW and _printable(data[hi]):
            hi += 1
        text = "".join(chr(c) if _printable(c) else f"\\x{c:02x}" for c in data[lo:hi])
        shown.append(f"{text} [{', '.join(labels)}]")
    return shown, len(found)


# --- binaries ---


@dataclass
class Binary:
    kind: str  # elf, macho or pe
    machines: list[str] = field(default_factory=list)
    needed: list[str] = field(default_factory=list)
    runpath: list[str] = field(default_factory=list)  # DT_RUNPATH; LC_RPATH
    rpath: list[str] = field(default_factory=list)  # DT_RPATH (ELF)
    soname: str = ""  # DT_SONAME; LC_ID_DYLIB
    defined: set[str] = field(default_factory=set)
    minos: list[tuple[int, int, int]] = field(default_factory=list)
    sections: list[tuple[str, bytes | None]] = field(default_factory=list)  # compressed (ELF)
    problems: list[str] = field(default_factory=list)


def _cstr(data: bytes, off: int, end: int | None = None) -> str:
    stop = data.find(b"\0", off, len(data) if end is None else end)
    return data[off : (len(data) if end is None else end) if stop < 0 else stop].decode(
        "utf-8", "surrogateescape"
    )


ELF_MACHINES = {3: "x86", 62: "x86_64", 183: "arm64"}
SHT_SYMTAB, SHT_DYNAMIC, SHT_NOBITS, SHT_DYNSYM = 2, 6, 8, 11
DT_NEEDED, DT_SONAME, DT_RPATH, DT_RUNPATH = 1, 14, 15, 29
SHF_COMPRESSED = 0x800


def parse_elf(data: bytes) -> Binary:
    b = Binary("elf")
    is64 = data[4] == 2
    end = "<" if data[5] == 1 else ">"
    machine = struct.unpack_from(end + "H", data, 18)[0]
    b.machines.append(ELF_MACHINES.get(machine, f"ELF machine {machine}"))
    if is64:
        shoff = struct.unpack_from(end + "Q", data, 0x28)[0]
        shentsize, shnum, shstrndx = struct.unpack_from(end + "HHH", data, 0x3A)
        shfmt, symfmt, dynfmt, chfmt = end + "IIQQQQIIQQ", end + "IBBHQQ", end + "qQ", end + "IIQQ"
    else:
        shoff = struct.unpack_from(end + "I", data, 0x20)[0]
        shentsize, shnum, shstrndx = struct.unpack_from(end + "HHH", data, 0x2E)
        shfmt, symfmt, dynfmt, chfmt = end + "IIIIIIIIII", end + "IIIBBH", end + "iI", end + "III"
    secs = [struct.unpack_from(shfmt, data, shoff + i * shentsize) for i in range(shnum)]

    def body(i: int) -> bytes:  # (name, type, flags, addr, offset, size, link, ...)
        s = secs[i]
        return b"" if s[1] == SHT_NOBITS else data[s[4] : s[4] + s[5]]

    names = body(shstrndx) if shstrndx < len(secs) else b""
    for i, s in enumerate(secs):
        typ, flags, link = s[1], s[2], s[6]
        if typ == SHT_DYNAMIC and link < len(secs):
            strs, dyn = body(link), body(i)
            size = struct.calcsize(dynfmt)
            for j in range(0, len(dyn) - size + 1, size):
                tag, val = struct.unpack_from(dynfmt, dyn, j)
                if tag == 0:
                    break
                if tag == DT_NEEDED:
                    b.needed.append(_cstr(strs, val))
                elif tag == DT_SONAME:
                    b.soname = _cstr(strs, val)
                elif tag == DT_RUNPATH:
                    b.runpath += _cstr(strs, val).split(":")
                elif tag == DT_RPATH:
                    b.rpath += _cstr(strs, val).split(":")
        elif typ in (SHT_SYMTAB, SHT_DYNSYM) and link < len(secs):
            strs, syms = body(link), body(i)
            size = struct.calcsize(symfmt)
            for j in range(0, len(syms) - size + 1, size):
                f = struct.unpack_from(symfmt, syms, j)
                name, shndx = f[0], (f[3] if is64 else f[5])
                if name and shndx != 0:  # SHN_UNDEF: needed from another library
                    b.defined.add(_cstr(strs, name))
        name = _cstr(names, s[0])
        raw = body(i)
        try:
            if flags & SHF_COMPRESSED:  # an Elf_Chdr, then the data
                if struct.unpack_from(chfmt, raw, 0)[0] == 1:  # ELFCOMPRESS_ZLIB
                    b.sections.append((name, zlib.decompress(raw[struct.calcsize(chfmt) :])))
                else:  # zstd, which the standard library reads only from 3.14
                    b.sections.append((name, None))
            elif name.startswith(".zdebug") and raw[:4] == b"ZLIB":
                b.sections.append((name, zlib.decompress(raw[12:])))
        except (struct.error, zlib.error):
            b.sections.append((name, None))
    return b


MACHO_CPUS = {7: "x86", 0x01000007: "x86_64", 0x0100000C: "arm64", 0x0200000C: "arm64_32"}
LC_SYMTAB, LC_LOAD_DYLIB, LC_ID_DYLIB, LC_RPATH = 0x2, 0xC, 0xD, 0x8000001C
LC_LOAD_DYLIBS = {LC_LOAD_DYLIB, 0x80000018, 0x8000001F, 0x20, 0x80000023}  # weak, reexport, ...
LC_BUILD_VERSION, LC_VERSION_MIN_MACOSX = 0x32, 0x24


def _macho_version(v: int) -> tuple[int, int, int]:
    return v >> 16, (v >> 8) & 0xFF, v & 0xFF


def _parse_macho_slice(data: bytes, off: int, b: Binary) -> None:
    magic = struct.unpack_from("<I", data, off)[0]
    end = "<" if magic in (0xFEEDFACF, 0xFEEDFACE) else ">"
    is64 = magic in (0xFEEDFACF, 0xCFFAEDFE)
    cpu, _, _, ncmds, _, _ = struct.unpack_from(end + "iiIIII", data, off + 4)
    b.machines.append(MACHO_CPUS.get(cpu & 0xFFFFFFFF, f"Mach-O CPU {cpu:#x}"))
    p = off + (32 if is64 else 28)
    for _ in range(ncmds):
        cmd, size = struct.unpack_from(end + "II", data, p)
        if size < 8:
            raise struct.error("a load command of less than 8 bytes")
        if cmd in LC_LOAD_DYLIBS or cmd in (LC_ID_DYLIB, LC_RPATH):
            name = _cstr(data, p + struct.unpack_from(end + "I", data, p + 8)[0], p + size)
            if cmd == LC_ID_DYLIB:
                b.soname = name
            elif cmd == LC_RPATH:
                b.runpath.append(name)
            else:
                b.needed.append(name)
        elif cmd == LC_BUILD_VERSION:
            b.minos.append(_macho_version(struct.unpack_from(end + "II", data, p + 8)[1]))
        elif cmd == LC_VERSION_MIN_MACOSX:
            b.minos.append(_macho_version(struct.unpack_from(end + "I", data, p + 8)[0]))
        elif cmd == LC_SYMTAB:
            symoff, nsyms, stroff, _ = struct.unpack_from(end + "IIII", data, p + 8)
            entry = 16 if is64 else 12
            for k in range(nsyms):
                strx, ntype = struct.unpack_from(end + "IB", data, off + symoff + k * entry)
                if not ntype & 0xE0 and ntype & 0x0E:  # not a debugging entry, not undefined
                    b.defined.add(_cstr(data, off + stroff + strx))
        p += size


def parse_macho(data: bytes) -> Binary:
    b = Binary("macho")
    magic = struct.unpack_from(">I", data, 0)[0]
    if magic in (0xCAFEBABE, 0xCAFEBABF):  # a universal binary: its slices
        n = struct.unpack_from(">I", data, 4)[0]
        fmt = ">iiIII" if magic == 0xCAFEBABE else ">iiQQII"
        for i in range(n):
            arch = struct.unpack_from(fmt, data, 8 + i * struct.calcsize(fmt))
            _parse_macho_slice(data, arch[2], b)
    else:
        _parse_macho_slice(data, 0, b)
    return b


PE_MACHINES = {0x14C: "x86", 0x8664: "x86_64", 0xAA64: "arm64"}


def parse_pe(data: bytes) -> Binary:
    b = Binary("pe")
    pe = struct.unpack_from("<I", data, 0x3C)[0]
    machine, nsec, _, _, _, optsize, _ = struct.unpack_from("<HHIIIHH", data, pe + 4)
    b.machines.append(PE_MACHINES.get(machine, f"PE machine {machine:#x}"))
    opt = pe + 24
    magic = struct.unpack_from("<H", data, opt)[0]
    ndirs_at = opt + (108 if magic == 0x20B else 92)  # PE32+, PE32
    ndirs = min(struct.unpack_from("<I", data, ndirs_at)[0], 16)
    dirs = [struct.unpack_from("<II", data, ndirs_at + 4 + 8 * i) for i in range(ndirs)]
    secs = []
    for i in range(nsec):
        _, vsize, vaddr, rawsize, rawptr = struct.unpack_from(
            "<8sIIII", data, opt + optsize + 40 * i
        )
        secs.append((vaddr, max(vsize, rawsize), rawptr))

    def at(rva: int) -> int:
        for vaddr, size, rawptr in secs:
            if vaddr <= rva < vaddr + size:
                return rawptr + rva - vaddr
        raise struct.error(f"RVA {rva:#x} outside the sections")

    def directory(i: int) -> int | None:
        return at(dirs[i][0]) if i < len(dirs) and dirs[i][0] and dirs[i][1] else None

    if (p := directory(1)) is not None:  # imports: 20-byte descriptors, the name at 12
        while any(struct.unpack_from("<5I", data, p)):
            b.needed.append(_cstr(data, at(struct.unpack_from("<I", data, p + 12)[0])))
            p += 20
    if (p := directory(13)) is not None:  # delayed imports: 32-byte descriptors, the name at 4
        while name := struct.unpack_from("<I", data, p + 4)[0]:
            b.needed.append(_cstr(data, at(name)))
            p += 32
    if (p := directory(0)) is not None:  # exports
        nnames, names = (
            struct.unpack_from("<I", data, p + 24)[0],
            struct.unpack_from("<I", data, p + 32)[0],
        )
        for k in range(nnames):
            b.defined.add(_cstr(data, at(struct.unpack_from("<I", data, at(names) + 4 * k)[0])))
    return b


def _u32be(data: bytes, off: int) -> int:
    return struct.unpack_from(">I", data, off)[0] if len(data) >= off + 4 else 0


def parse_binary(data: bytes) -> Binary | None:
    """The libraries, search paths, symbols and targets of an ELF, Mach-O or PE file, or None
    for any other file. A file that looks like one but can't be read has problems."""
    head = data[:4]
    if head == b"\x7fELF":
        parse = parse_elf
    elif head in (
        b"\xcf\xfa\xed\xfe",
        b"\xce\xfa\xed\xfe",
        b"\xfe\xed\xfa\xcf",
        b"\xfe\xed\xfa\xce",
    ):
        parse = parse_macho
    elif head in (b"\xca\xfe\xba\xbe", b"\xca\xfe\xba\xbf") and 0 < _u32be(data, 4) < 32:
        parse = parse_macho  # (a Java class file has the same magic, with a large count)
    elif data[:2] == b"MZ" and len(data) > 0x40:
        pe = struct.unpack_from("<I", data, 0x3C)[0]
        if data[pe : pe + 4] != b"PE\0\0":
            return None
        parse = parse_pe
    else:
        return None
    try:
        return parse(data)
    except (struct.error, IndexError, ValueError) as e:
        b = Binary(parse.__name__[6:])
        b.problems.append(f"can't be read as {b.kind}: {e}")
        return b


# --- wheels ---


def _tags(name: str) -> tuple[str, set[str], set[str], set[str]]:
    """The distribution and the Python, ABI and platform tags of a wheel's file name."""
    parts = name[: -len(".whl")].split("-")
    if len(parts) not in (5, 6):
        raise ValueError(f"not a wheel file name: {name}")
    return parts[0], set(parts[-3].split(".")), set(parts[-2].split(".")), set(parts[-1].split("."))


def _system(plats: set[str]) -> str:
    if all(p.startswith(("manylinux", "musllinux", "linux")) for p in plats):
        return "linux"
    if all(p.startswith("macosx_") for p in plats):
        return "macos"
    if all(p.startswith("win") for p in plats):
        return "windows"
    return "?"


def _arch(plats: set[str]) -> str:
    """The architecture of the platform tags (as the binaries' machines are named), or ?."""
    found = set()
    for p in plats:
        a = "x86" if p == "win32" else next((x for x in ARCHS if p.endswith("_" + x)), "?")
        found.add(ARCHS.get(a, a))
    return found.pop() if len(found) == 1 else "?"


ARCHS = {"x86_64": "x86_64", "amd64": "x86_64", "aarch64": "arm64", "arm64": "arm64", "i686": "x86"}
MODULE_SUFFIX = re.compile(r"\.(?:abi3|cpython-\d+[a-z]*-[\w-]+|cp\d+[a-z]*-win\w*)\.(?:so|pyd)$")


def _is_module(path: str) -> bool:
    return bool(MODULE_SUFFIX.search(path)) or path.endswith(".pyd")


def _python_dll(abi: set[str], pys: set[str]) -> set[str]:
    """The Python DLLs a Windows module may load: python3.dll for abi3, python3XX.dll for cpXX."""
    if "abi3" in abi:
        return {"python3.dll"}
    return {f"python{a[2:]}.dll" for a in abi if a.startswith("cp")} | {
        f"python{p[2:]}.dll" for p in pys if p.startswith("cp") and "none" not in abi
    }


class Report:
    def __init__(self, path: str, quiet: bool) -> None:
        self.path, self.quiet, self.findings = path, quiet, 0

    def info(self, text: str) -> None:
        if not self.quiet:
            print(f"  {text}")

    def fail(self, text: str) -> None:
        self.findings += 1
        print(f"{self.path}: {text}")


def _licenses(r: Report, meta: str, profile: Profile, files: set[str], prefix: str) -> None:
    msg = email.parser.HeaderParser().parsestr(meta)
    expr = msg.get("License-Expression")
    if expr != profile.license:
        r.fail(f"License-Expression is {expr!r}, not {profile.license!r}")
    declared = msg.get_all("License-File") or []
    for f in sorted(set(declared) | set(profile.license_files)):
        if f not in declared:
            r.fail(f"the license file {f} is not declared (License-File)")
        if prefix + f not in files:
            r.fail(f"the license file {f} is not in {prefix}")
    r.info(f"license {expr}, files: {' '.join(declared) or '-'}")


def _scan(r: Report, name: str, data: bytes, patterns, sections=()) -> None:
    for part, body in [(name, data), *((f"{name}[{s}]", d) for s, d in sections)]:
        if body is None:
            r.fail(f"{part}: a compressed section that can't be read here (not searched)")
            continue
        shown, total = scan_bytes(body, patterns)
        for s in shown:
            r.fail(f"{part}: {s}")
        if total > len(shown):
            r.fail(f"{part}: ... and {total - len(shown)} more")


def _resolve_elf(path: str, needed: str, b: Binary, members: set[str]) -> bool:
    folder = str(PurePosixPath(path).parent)
    for entry in b.runpath:
        entry = entry.replace("${ORIGIN}", "$ORIGIN")
        if entry.startswith("$ORIGIN"):
            p = posixpath.normpath(folder + entry[len("$ORIGIN") :] + "/" + needed)
            if p in members:
                return True
    return False


def _resolve_macho(path: str, name: str, b: Binary, members: set[str]) -> bool:
    folder = str(PurePosixPath(path).parent)
    if name.startswith("@loader_path/"):
        candidates = [name]
    elif name.startswith("@rpath/"):
        candidates = [e + name[len("@rpath") :] for e in b.runpath if e.startswith("@loader_path")]
    else:
        return False
    for c in candidates:
        p = posixpath.normpath(folder + c[len("@loader_path") :])
        if p in members:
            return True
    return False


def _pe_search_path(path: str, binaries: dict[str, Binary], dist: str) -> set[str]:
    """The DLLs a Windows binary finds in the wheel: those next to it and in <dist>.libs/, the
    folder that delvewheel's patch adds to the DLL search path."""
    here = PurePosixPath(path).parent
    return {
        PurePosixPath(x).name.lower()
        for x in binaries
        if PurePosixPath(x).parent == here or x.startswith(f"{dist}.libs/")
    }


def check_wheel(path: str, patterns, profile_name: str | None, quiet: bool) -> int:
    name = os.path.basename(path)
    dist, pys, abis, plats = _tags(name)
    profile = PROFILES.get(profile_name or dist)
    if profile is None:
        raise ValueError(f"no profile for the distribution {dist} (one of {', '.join(PROFILES)})")
    r = Report(path, quiet)
    system, arch = _system(plats), _arch(plats)
    print(f"{path}: {dist}, {'.'.join(sorted(pys))}-{'.'.join(sorted(abis))}, {system} {arch}")
    with zipfile.ZipFile(path) as z:
        infos = [i for i in z.infolist() if not i.is_dir()]
        members = {i.filename for i in infos}
        data = {i.filename: z.read(i) for i in infos}
    distinfo = [m for m in members if re.fullmatch(r"[^/]+\.dist-info/WHEEL", m)]
    if len(distinfo) != 1:
        r.fail(f"{len(distinfo)} .dist-info/WHEEL files")
        return r.findings
    di = distinfo[0].rsplit("/", 1)[0] + "/"

    # the tags
    wheel = email.parser.HeaderParser().parsestr(data[di + "WHEEL"].decode("utf-8", "replace"))
    listed = set(wheel.get_all("Tag") or [])
    named = {f"{p}-{a}-{q}" for p in pys for a in abis for q in plats}
    if listed != named:
        r.fail(f"WHEEL's tags {sorted(listed)} are not the file name's {sorted(named)}")
    if wheel.get("Root-Is-Purelib", "").lower() != "false":
        r.fail("WHEEL's Root-Is-Purelib is not false")
    if system == "?":
        r.fail(f"no single system in the platform tags {sorted(plats)}")
    want_kind = {"linux": "elf", "macos": "macho", "windows": "pe"}.get(system)
    target = None
    if system == "macos":
        target = min(tuple(int(x) for x in p.split("_")[1:3]) for p in plats)

    # members: compiled Python, the folder check
    for m in sorted(members):
        if m.endswith((".pyc", ".pyo")) or "/__pycache__/" in f"/{m}":
            r.fail(f"{m}: compiled Python in the wheel")
        _scan(r, f"{m} (name)", b"/" + m.encode("utf-8", "surrogateescape"), patterns)

    # binaries
    binaries = {}
    for m in sorted(members):
        b = parse_binary(data[m])
        if b is not None:
            binaries[m] = b
    modules = sorted(m for m in binaries if _is_module(m))
    libs = sorted(m for m in binaries if not _is_module(m))
    pydlls = _python_dll(abis, pys)
    for m, b in binaries.items():
        for p in b.problems:
            r.fail(f"{m}: {p}")
        if want_kind and b.kind != want_kind:
            r.fail(f"{m}: a {b.kind} file in a {system} wheel")
        if arch != "?" and any(a != arch for a in b.machines):
            r.fail(f"{m}: built for {', '.join(b.machines)}, not {arch}")
        if target is not None:
            wheel_target = f"{target[0]}.{target[1]}"
            for v in b.minos:
                if v[:2] > target:
                    r.fail(f"{m}: deployment target {v[0]}.{v[1]} is newer than {wheel_target}")
        if b.kind == "elf":
            for e in b.rpath:
                r.fail(f"{m}: RPATH {e!r} (RUNPATH $ORIGIN only)")
            for e in b.runpath:
                if e not in ("$ORIGIN", "${ORIGIN}"):
                    r.fail(f"{m}: RUNPATH {e!r} (only $ORIGIN)")
            for n in b.needed:
                if n not in ELF_SYSTEM and not _resolve_elf(m, n, b, members):
                    r.fail(
                        f"{m}: needs {n}, which is neither in the wheel (RUNPATH) nor the system's"
                    )
        elif b.kind == "macho":
            for e in b.runpath:
                if not e.startswith("@loader_path"):
                    r.fail(f"{m}: LC_RPATH {e!r} (only @loader_path)")
            # (delocate names its copies /DLC/<package>/...: a unique name, no folder)
            if b.soname.startswith("/") and not b.soname.startswith((*MACHO_SYSTEM, "/DLC/")):
                r.fail(f"{m}: install name {b.soname!r} is a folder outside the wheel")
            for n in b.needed:
                if not n.startswith(MACHO_SYSTEM) and not _resolve_macho(m, n, b, members):
                    r.fail(f"{m}: loads {n}, which is neither in the wheel nor the system's")
        elif b.kind == "pe":
            near = _pe_search_path(m, binaries, dist)
            for n in b.needed:
                low = n.lower()
                if re.fullmatch(r"python3\d*t?\.dll", low):
                    if low not in pydlls:
                        want = " or ".join(sorted(pydlls)) or "no Python DLL"
                        r.fail(f"{m}: loads {n}, not {want}")
                elif not (low in PE_SYSTEM or low.startswith(PE_SYSTEM_PREFIXES) or low in near):
                    r.fail(f"{m}: loads {n}, which is neither in the wheel nor the system's")
        if not profile.fftw:
            fftw = sorted(s for s in b.defined if FFTW_SYMBOL.match(s))
            if fftw:
                more = " ..." if len(fftw) > 5 else ""
                r.fail(f"{m}: defines FFTW symbols ({', '.join(fftw[:5])}{more})")
            hits = {x.decode() for x in FFTW_STRING.findall(data[m])}
            if hits:
                r.fail(f"{m}: has FFTW's version string ({', '.join(sorted(hits))})")

    # the modules' names
    if profile.modules and not modules:
        r.fail("no extension modules")
    if not profile.modules and modules:
        r.fail(f"extension modules in a library wheel: {', '.join(modules)}")
    for m in modules:
        base = PurePosixPath(m).name
        if "abi3" in abis:
            ok = base.endswith(".abi3.so") if system != "windows" else base.count(".") == 1
        else:
            vers = {a[2:] for a in abis if a.startswith("cp")}
            ok = any(f".cpython-{v}-" in base or f".cp{v}-" in base for v in vers)
        if not ok:
            r.fail(f"{m}: the module's name does not fit the ABI tag {'.'.join(sorted(abis))}")

    # the bundled libraries
    expected = profile.libraries.get(system, ())
    for m in libs:
        if not any(re.fullmatch(p, m) for p, _ in expected):
            r.fail(f"{m}: a library the wheel should not carry")
    for p, required in expected:
        n = sum(1 for m in libs if re.fullmatch(p, m))
        if required and n != 1:
            r.fail(f"{n} libraries match {p} (expected one)")
    for m in libs:
        b = binaries[m]
        name = "install name" if b.kind == "macho" else "SONAME"
        needs = " ".join(b.needed) or "-"
        r.info(f"library {m}: {len(data[m])} bytes, {name} {b.soname or '-'}, needs {needs}")
    sets: dict[str, list[str]] = {}
    for m in modules:
        b = binaries[m]
        paths = " ".join(b.runpath + [f"RPATH:{e}" for e in b.rpath]) or "-"
        sets.setdefault(f"needs {' '.join(b.needed) or '-'} | search path {paths}", []).append(
            PurePosixPath(m).name.split(".")[0]
        )
    for k, v in sets.items():
        r.info(f"{len(v)} modules ({' '.join(v)}): {k}")
    newest = max((v for b in binaries.values() for v in b.minos), default=None)
    if target is not None and newest:
        r.info(f"deployment target {target[0]}.{target[1]}, the binaries' {newest[0]}.{newest[1]}")

    # the folder check of the binaries' compressed sections and every member's bytes
    for m in sorted(members):
        _scan(r, m, data[m], patterns, binaries[m].sections if m in binaries else ())

    # the license
    meta = data.get(di + "METADATA")
    if meta is None:
        r.fail("no METADATA")
    else:
        _licenses(r, meta.decode("utf-8", "replace"), profile, members, di + "licenses/")
    return r.findings


# --- the source distribution ---

SDIST_FORBIDDEN = [
    ("the exact FFT add-on (GPL), a package of its own", re.compile(r"cpp/fft/exact(?:/|$)")),
    ("a build folder", re.compile(r"(?:build|dist|_skbuild|\.venv|[^/]+\.egg-info)(?:/|$)")),
    ("compiled Python", re.compile(r"(?:.*/)?(?:__pycache__/.*|[^/]+\.py[co])$")),
    (
        "a binary",
        re.compile(r"(?:.*/)?[^/]+\.(?:so(?:\.\d+)*|pyd|dylib|dll|exe|o|obj|a|lib|whl)$", re.I),
    ),
]


def check_sdist(path: str, patterns, profile_name: str | None, quiet: bool, forbid=()) -> int:
    name = os.path.basename(path)
    m = re.fullmatch(r"([A-Za-z0-9_.]+)-([^-]+)\.tar\.gz", name)
    if not m:
        raise ValueError(f"not a source distribution's file name: {name}")
    dist = m.group(1).replace(".", "_").replace("-", "_").lower()
    profile = PROFILES.get(profile_name or dist)
    if profile is None:
        raise ValueError(f"no profile for the distribution {dist} (one of {', '.join(PROFILES)})")
    r = Report(path, quiet)
    print(f"{path}: {dist} source distribution")
    base = name[: -len(".tar.gz")] + "/"
    files: dict[str, bytes] = {}
    owners = set()
    with tarfile.open(path, "r:gz") as t:
        for info in t:
            if not info.name.startswith(base) and info.name != base.rstrip("/"):
                r.fail(f"{info.name}: outside the folder {base}")
            _scan(
                r,
                f"{info.name} (name)",
                b"/" + info.name.encode("utf-8", "surrogateescape"),
                patterns,
            )
            if info.linkname:
                _scan(
                    r,
                    f"{info.name} (link)",
                    info.linkname.encode("utf-8", "surrogateescape"),
                    patterns,
                )
            for owner in (info.uname, info.gname):
                if owner and owner != "root":
                    owners.add(owner)
            if info.isfile():
                files[info.name] = t.extractfile(info).read()
    for o in sorted(owners):
        r.fail(f"members are owned by the name {o!r} (no owner names)")
    for f in sorted(files):
        rel = f[len(base) :]
        for why, rx in SDIST_FORBIDDEN:
            if rx.match(rel):
                r.fail(f"{f}: {why}")
        if parse_binary(files[f]) is not None:
            r.fail(f"{f}: a binary (ELF, Mach-O or PE)")
        # the sources are the repository's files, in which build folders' names can be legitimate
        # (search hints, this check's own patterns and cases): their contents get the --forbid
        # patterns only, the generated PKG-INFO every pattern
        _scan(r, f, files[f], patterns if rel == "PKG-INFO" else forbid)
    for need in ("PKG-INFO", "pyproject.toml"):
        if base + need not in files:
            r.fail(f"no {need}")
    if base + "PKG-INFO" in files:
        _licenses(r, files[base + "PKG-INFO"].decode("utf-8", "replace"), profile, set(files), base)
    r.info(f"{len(files)} files")
    return r.findings


# --- the command ---


def _expand(paths: list[str]) -> list[str]:
    out = []
    for p in paths:
        if os.path.isdir(p):
            found = sorted(
                os.path.join(p, f)
                for f in os.listdir(p)
                if f.endswith((".whl", ".tar.gz")) and os.path.isfile(os.path.join(p, f))
            )
            if not found:
                raise FileNotFoundError(f"no wheels or source distributions in {p}")
            out += found
        elif os.path.isfile(p):
            out.append(p)
        else:
            raise FileNotFoundError(f"no archive or folder {p}")
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("paths", nargs="+", metavar="ARCHIVE|FOLDER")
    ap.add_argument(
        "--forbid",
        action="append",
        default=[],
        metavar="REGEX",
        help="another pattern to fail on, such as a build machine's or a site's names",
    )
    ap.add_argument(
        "--forbid-file", metavar="FILE", help="patterns to fail on, one per line (# comments)"
    )
    ap.add_argument("--profile", choices=sorted(PROFILES), help="instead of the file name's")
    ap.add_argument("-q", "--quiet", action="store_true", help="print the findings only")
    a = ap.parse_args(argv)
    extra = list(a.forbid)
    try:
        if a.forbid_file:
            with open(a.forbid_file, encoding="utf-8") as f:
                extra += [ln.strip() for ln in f if ln.strip() and not ln.lstrip().startswith("#")]
        forbid = [(f"--forbid {i + 1}", re.compile(rx.encode())) for i, rx in enumerate(extra)]
        patterns = [(label, re.compile(rx)) for label, rx in PATH_PATTERNS] + forbid
        archives = _expand(a.paths)
    except (OSError, re.error) as e:
        print(f"check_wheel.py: {e}", file=sys.stderr)
        return 2
    bad = 0
    for p in archives:
        try:
            if p.endswith(".whl"):
                n = check_wheel(p, patterns, a.profile, a.quiet)
            else:
                n = check_sdist(p, patterns, a.profile, a.quiet, forbid)
        except (OSError, ValueError, KeyError, zipfile.BadZipFile, tarfile.TarError, EOFError) as e:
            print(f"{p}: can't be read: {e}", file=sys.stderr)
            return 2
        print(f"{p}: {'FAIL' if n else 'ok'}, {n} findings")
        bad += n > 0
    print(f"{len(archives)} archives, {bad} with findings")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
