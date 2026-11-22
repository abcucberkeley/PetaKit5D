"""The toolchain of a source build on Windows (docs/platforms.md, "Windows: the toolchain of a
source build"): pyproject.toml's overrides, read by scikit-build-core as on Windows, and
cmake/Pk5dToolchain.cmake's choice and check, run in CMake's script mode with stand-in compilers
(POSIX shell scripts that answer as MinGW-w64's gcc and g++ do), so that both run off Windows. The
real toolchains are the CI's: the Python workflow's Windows jobs build with MinGW-builds' GCC 12.2.0
where they put it and with MSVC where the runners' own MinGW-w64 is left out, and check which one
the build took and why."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import tomllib
from pathlib import Path

import pytest

PY = Path(__file__).resolve().parents[2]
MODULE = PY / "cmake" / "Pk5dToolchain.cmake"
GET_MINGW = PY / "tools" / "get_mingw.sh"
# petakit5d and the exact add-on, whose source tree stands alone with a copy of the module
ADDON = PY / "cpp" / "fft" / "exact"
PROJECTS = {"petakit5d": PY, "petakit5d-fftw": ADDON}

# the banners of `gcc --version`: tools/get_mingw.sh's toolchain, the GitHub windows-2025 image's
# own (MinGW-builds' 15.2.0, also UCRT, in C:\mingw64\bin on PATH) and MSYS2's
PINNED = "gcc.exe (x86_64-posix-seh-rev2, Built by MinGW-W64 project) 12.2.0"
RUNNER = "gcc.exe (x86_64-posix-seh-rev0, Built by MinGW-Builds project) 15.2.0"
MSYS2 = "gcc.exe (Rev2, Built by MSYS2 project) 14.2.0"
UCRT, MSVCRT = "0xE00", "0x700"
FAIL = "cc1.exe: libwinpthread-1.dll was not found"
# the folder of tools/get_mingw.sh's toolchain under LOCALAPPDATA
GET_MINGW_DIR = "petakit5d/mingw-12.2.0"

# --- pyproject.toml's overrides, as scikit-build-core reads them ---

ENV_NAMES = ("CC", "CXX", "CMAKE_GENERATOR", "CMAKE_ARGS", "CMAKE_TOOLCHAIN_FILE")


@pytest.fixture(params=sorted(PROJECTS))
def project(request) -> Path:
    return PROJECTS[request.param]


def _settings(monkeypatch, project: Path, platform: str, env: dict, retry: bool = False):
    """scikit-build-core's settings of the project's pyproject.toml on `platform` with the
    environment `env` alone, at the first attempt or (retry) the second: (cmake.args,
    cmake.define, whether a failed first attempt is tried again, messages.after-failure)"""
    skb = pytest.importorskip("scikit_build_core.settings.skbuild_read_settings")
    pyproject = tomllib.loads((project / "pyproject.toml").read_text(encoding="utf-8"))
    with monkeypatch.context() as m:
        m.setattr(sys, "platform", platform)
        r = skb.SettingsReader(pyproject, {}, state="wheel", env=env, retry=retry)
    s = r.settings
    return (
        list(s.cmake.args),
        dict(s.cmake.define),
        "failed" in r.overrides,
        s.messages.after_failure,
    )


def test_linux_and_macos_builds_are_unchanged(monkeypatch, project):
    for platform in ("linux", "darwin"):
        for retry in (False, True):
            args, define, again, _ = _settings(monkeypatch, project, platform, {}, retry)
            assert (args, define, again) == ([], {}, False)


def test_windows_tries_mingw_first_then_msvc(monkeypatch, project):
    args, define, _, _ = _settings(monkeypatch, project, "win32", {})
    assert args == ["-GNinja"]
    assert define == {"PK5D_WINDOWS_TOOLCHAIN": "default"}
    # the second attempt: scikit-build-core's own default (Visual Studio's generator, MSVC)
    args, define, again, msg = _settings(monkeypatch, project, "win32", {}, retry=True)
    assert (args, define, again) == ([], {"PK5D_WINDOWS_TOOLCHAIN": "fallback"}, True)
    assert "MinGW-builds' GCC 12.2.0" in msg and "MSVC" in msg
    # an empty variable is an unset one, for CMake too
    env = {"CC": "", "CMAKE_ARGS": ""}
    args, define, _, _ = _settings(monkeypatch, project, "win32", env)
    assert (args, define) == (["-GNinja"], {"PK5D_WINDOWS_TOOLCHAIN": "default"})


@pytest.mark.parametrize("name", ENV_NAMES)
def test_a_toolchain_the_environment_names_wins(monkeypatch, project, name):
    value = {"CC": "gcc", "CXX": "cl", "CMAKE_GENERATOR": "Ninja"}.get(name, "x")
    for retry in (False, True):
        for switch in ({}, {"PETAKIT5D_WINDOWS_TOOLCHAIN": "mingw"}):
            env = {name: value, **switch}
            assert _settings(monkeypatch, project, "win32", env, retry)[:3] == ([], {}, False)


def test_the_switch(monkeypatch, project):
    on = {"PETAKIT5D_WINDOWS_TOOLCHAIN": "mingw"}
    for retry in (False, True):
        # mingw: that toolchain or an error, no second attempt
        args, define, again, _ = _settings(monkeypatch, project, "win32", on, retry)
        assert (args, define, again) == (["-GNinja"], {"PK5D_WINDOWS_TOOLCHAIN": "mingw"}, False)
        # msvc, and any other value (which the CMake side refuses): scikit-build-core's default
        for value in ("msvc", "MinGW", "auto", "xmingw"):
            env = {"PETAKIT5D_WINDOWS_TOOLCHAIN": value}
            assert _settings(monkeypatch, project, "win32", env, retry)[:3] == ([], {}, False)


def test_the_addon_takes_the_same_overrides_and_module():
    # its own copy of the module (its source tree stands alone), the same bytes, and the same
    # overrides but for the failure's message
    assert (ADDON / "cmake" / "Pk5dToolchain.cmake").read_bytes() == MODULE.read_bytes()

    def overrides(project):
        t = tomllib.loads((project / "pyproject.toml").read_text(encoding="utf-8"))
        rows = t["tool"]["scikit-build"]["overrides"]
        for r in rows:
            r.get("messages", {}).pop("after-failure", None)
        return rows

    assert overrides(ADDON) == overrides(PY)
    text = (ADDON / "CMakeLists.txt").read_text(encoding="utf-8")
    call = text.index("pk5d_windows_toolchain()")
    assert text.index("include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Pk5dToolchain.cmake)") < call
    assert call < text.index("project(petakit5d_fftw LANGUAGES C)")


# --- cmake/Pk5dToolchain.cmake in script mode, with stand-in compilers ---


def _cmake() -> str:
    exe = shutil.which("cmake")
    if exe is None:
        pytest.skip("no cmake on PATH")
    if os.name == "nt":
        pytest.skip("the stand-in compilers are POSIX shell scripts (Windows: the CI's real ones)")
    out = subprocess.run([exe, "--version"], capture_output=True, text=True, check=True).stdout
    if tuple(int(x) for x in re.search(r"(\d+)\.(\d+)", out).groups()) < (3, 26):
        pytest.skip(f"cmake older than the build's 3.26: {out.splitlines()[0]}")
    return exe


def _compiler(folder: Path, banner: str | None = PINNED, crt: str | None = UCRT, gxx=True) -> Path:
    """gcc (and g++) in folder: `--version` prints banner (None: the pinned one, then fails),
    `-x c -dM -E <file>` the C library's macros (__MSVCRT_VERSION__ crt; None: none; FAIL: fails,
    as when cc1 can't load its DLLs)"""
    folder.mkdir(parents=True, exist_ok=True)
    version = (
        f'echo "{banner}"; echo "Copyright (C) the FSF"' if banner else f'echo "{PINNED}"; exit 1'
    )
    macros = f'echo "#define _UCRT"; echo "#define __MSVCRT_VERSION__ {crt}"' if crt else "true"
    if crt == FAIL:
        macros = f'echo "{FAIL}" >&2; exit 1'
    text = f'#!/bin/sh\ncase "$1" in\n  --version) {version};;\n  -x) {macros};;\n'
    text += "  *) exit 1;;\nesac\n"
    for name in ("gcc", "g++") if gxx else ("gcc",):
        (folder / name).write_text(text, newline="\n")
        (folder / name).chmod(0o755)
    return folder / "gcc"


def _choose(tmp_path, *, path=(), appdata=None, mode="default", generator="Ninja", switch=None,
            pre="", toolchain_file=None):  # fmt: skip
    """pk5d_windows_toolchain() with PATH `path`, LOCALAPPDATA `appdata`, PK5D_WINDOWS_TOOLCHAIN
    `mode`, the generator and PETAKIT5D_WINDOWS_TOOLCHAIN `switch`, after the CMake code `pre`:
    (exit status, {CC, CXX, WHY}, output)"""
    exe = _cmake()
    run = tmp_path / "run"
    run.mkdir(exist_ok=True)
    script = tmp_path / "choose.cmake"
    script.write_text(
        "cmake_minimum_required(VERSION 3.26)\n"
        f'include("{MODULE.as_posix()}")\n'
        f'set(CMAKE_GENERATOR "{generator}")\n'
        f'set(PK5D_WINDOWS_TOOLCHAIN "{mode}")\n'
        f"{pre}\n"
        "pk5d_windows_toolchain()\n"
        'message(STATUS "CC=${CMAKE_C_COMPILER}")\n'
        'message(STATUS "CXX=${CMAKE_CXX_COMPILER}")\n'
        'message(STATUS "WHY=${PK5D_TOOLCHAIN_WHY}")\n',
        newline="\n",
    )
    env = {"PATH": os.pathsep.join(str(p) for p in path) or str(tmp_path / "nothing")}
    if appdata is not None:
        env["LOCALAPPDATA"] = str(appdata)
    if switch is not None:
        env["PETAKIT5D_WINDOWS_TOOLCHAIN"] = switch
    if toolchain_file is not None:
        env["CMAKE_TOOLCHAIN_FILE"] = toolchain_file
    p = subprocess.run([exe, "-P", str(script)], cwd=run, env=env, capture_output=True, text=True)
    got = dict(re.findall(r"^-- (CC|CXX|WHY)=(.*)$", p.stdout, re.M))
    return p.returncode, got, p.stdout + p.stderr


def _flat(text: str) -> str:
    # (CMake wraps an error's message over lines)
    return " ".join(text.split())


def test_takes_the_pinned_toolchain_on_path(tmp_path):
    gcc = _compiler(tmp_path / "tc" / "bin")
    rc, got, out = _choose(tmp_path, path=[gcc.parent])
    assert rc == 0, out
    assert got == {
        "CC": gcc.as_posix(),
        "CXX": (gcc.parent / "g++").as_posix(),
        "WHY": "MinGW-builds' GCC 12.2.0 with the UCRT on PATH: the default of source builds",
    }
    rc, got, out = _choose(tmp_path, path=[gcc.parent], mode="mingw", switch="mingw")
    assert rc == 0, out
    assert got["WHY"].endswith("on PATH: PETAKIT5D_WINDOWS_TOOLCHAIN=mingw")


def test_takes_get_mingw_folder_when_the_first_gcc_on_path_is_another(tmp_path):
    other = _compiler(tmp_path / "runner" / "bin", RUNNER)
    appdata = tmp_path / "appdata"
    gcc = _compiler(appdata / GET_MINGW_DIR / "mingw64" / "bin")
    rc, got, out = _choose(tmp_path, path=[other.parent], appdata=appdata)
    assert rc == 0, out
    assert got["CC"] == gcc.as_posix()
    assert got["WHY"] == (
        "MinGW-builds' GCC 12.2.0 with the UCRT in tools/get_mingw.sh's folder: "
        "the default of source builds"
    )
    # the first gcc on PATH wins when it is the pinned one; one later on PATH is not looked at
    first = _compiler(tmp_path / "first" / "bin")
    assert _choose(tmp_path, path=[first.parent], appdata=appdata)[1]["CC"] == first.as_posix()
    rc, got, out = _choose(tmp_path, path=[other.parent, first.parent])
    assert rc != 0 and "CC" not in got


@pytest.mark.parametrize(
    "banner, crt, gxx, why",
    [
        (RUNNER, UCRT, True, f"is {RUNNER}"),
        (MSYS2, UCRT, True, f"is {MSYS2}"),
        # MinGW-builds' other builds: another revision, win32 threads, another 12
        (PINNED.replace("rev2", "rev1"), UCRT, True, "is " + PINNED.replace("rev2", "rev1")),
        (PINNED.replace("posix", "win32"), UCRT, True, "is " + PINNED.replace("posix", "win32")),
        (
            PINNED.replace("12.2.0", "12.1.0"),
            UCRT,
            True,
            "is " + PINNED.replace("12.2.0", "12.1.0"),
        ),
        (PINNED, MSVCRT, True, "is MinGW-builds' 12.2.0 for msvcrt.dll, not the UCRT"),
        (PINNED, None, True, "(12.2.0) names no C library (no __MSVCRT_VERSION__)"),
        (PINNED, FAIL, True, f"(12.2.0) fails to preprocess a C file (1: {FAIL})"),
        (None, UCRT, True, "does not run"),
        (PINNED, UCRT, False, "has no g++ next to it"),
    ],
)
def test_leaves_out_another_toolchain(tmp_path, banner, crt, gxx, why):
    gcc = _compiler(tmp_path / "tc" / "bin", banner, crt, gxx)
    appdata = tmp_path / "appdata"
    rc, got, out = _choose(tmp_path, path=[gcc.parent], appdata=appdata)
    flat = _flat(out)
    assert rc != 0 and "CC" not in got and "WHY" not in got, out
    assert f"gcc on PATH: {gcc.as_posix()} {why}" in flat
    assert f"none in tools/get_mingw.sh's folder {(appdata / GET_MINGW_DIR).as_posix()}" in flat
    assert "scikit-build-core builds again with MSVC" in flat
    # PETAKIT5D_WINDOWS_TOOLCHAIN=mingw: an error that says how to get it, no second attempt
    rc, got, out = _choose(tmp_path, path=[gcc.parent], mode="mingw", switch="mingw")
    flat = _flat(out)
    assert rc != 0 and "CC" not in got
    assert "PETAKIT5D_WINDOWS_TOOLCHAIN=mingw, but no MinGW-builds GCC 12.2.0" in flat
    assert "bash tools/get_mingw.sh fetches it" in flat
    assert "builds again" not in flat


def test_no_gcc_at_all(tmp_path):
    rc, got, out = _choose(tmp_path)
    assert rc != 0 and "CC" not in got
    assert "the default toolchain of source builds on Windows: no gcc on PATH" in _flat(out)


def test_a_decimal_runtime_version_counts(tmp_path):
    gcc = _compiler(tmp_path / "tc" / "bin", crt="3584")
    assert _choose(tmp_path, path=[gcc.parent])[1]["CC"] == gcc.as_posix()
    gcc = _compiler(tmp_path / "old" / "bin", crt="3583")
    assert _choose(tmp_path, path=[gcc.parent])[0] != 0


def test_named_compilers_and_generators_win(tmp_path):
    gcc = _compiler(tmp_path / "tc" / "bin")
    cases = [
        ({"pre": 'set(CMAKE_CXX_COMPILER "/opt/cl")'}, "/opt/cl"),
        ({"pre": 'set(CMAKE_C_COMPILER "/opt/cc")'}, ""),
        ({"pre": 'set(CMAKE_TOOLCHAIN_FILE "/opt/tc.cmake")'}, ""),
        ({"toolchain_file": "/opt/tc.cmake"}, ""),
        ({"generator": "Visual Studio 17 2022"}, ""),
        ({"generator": "Ninja Multi-Config"}, ""),
    ]
    for kw, cxx in cases:
        rc, got, out = _choose(tmp_path, path=[gcc.parent], **kw)
        assert rc == 0, out
        assert got["WHY"] == "as CMake's command line names it", kw
        assert got["CXX"] == cxx, kw


def test_without_the_overrides(tmp_path):
    gcc = _compiler(tmp_path / "tc" / "bin")
    rc, got, out = _choose(tmp_path, path=[gcc.parent], mode="")
    assert rc == 0, out
    assert got["CC"] == "" and got["WHY"].startswith("as CC, CXX, CMAKE_GENERATOR, CMAKE_ARGS")
    rc, got, out = _choose(tmp_path, path=[gcc.parent], mode="", switch="msvc")
    assert (rc, got["CC"], got["WHY"]) == (0, "", "PETAKIT5D_WINDOWS_TOOLCHAIN=msvc")
    rc, got, out = _choose(tmp_path, path=[gcc.parent], mode="fallback")
    assert rc == 0 and got["CC"] == ""
    assert got["WHY"] == "the fallback: the attempt with MinGW-w64 failed (its log is above)"


def test_a_configure_again_keeps_the_first_choice(tmp_path):
    pre = (
        'set(PK5D_TOOLCHAIN_WHY "chosen before" CACHE INTERNAL "")\n'
        'set(CMAKE_CXX_COMPILER "/opt/g++")'
    )
    rc, got, out = _choose(tmp_path, pre=pre)
    assert (rc, got["CXX"], got["WHY"]) == (0, "/opt/g++", "chosen before"), out


@pytest.mark.parametrize(
    "kw, msg",
    [
        (
            {"mode": "", "switch": "MinGW"},
            "PETAKIT5D_WINDOWS_TOOLCHAIN is 'MinGW': mingw, msvc or unset",
        ),
        ({"switch": "auto"}, "PETAKIT5D_WINDOWS_TOOLCHAIN is 'auto': mingw, msvc or unset"),
        ({"mode": "gcc"}, "PK5D_WINDOWS_TOOLCHAIN is 'gcc': default, fallback, mingw or unset"),
    ],
)
def test_refuses_unknown_values(tmp_path, kw, msg):
    rc, got, out = _choose(tmp_path, **kw)
    assert rc != 0 and msg in _flat(out)


def test_the_check_alone(tmp_path):
    exe = _cmake()
    good = _compiler(tmp_path / "a")
    bad = _compiler(tmp_path / "b", crt=MSVCRT)
    script = tmp_path / "check.cmake"
    script.write_text(
        "cmake_minimum_required(VERSION 3.26)\n"
        f'include("{MODULE.as_posix()}")\n'
        + "".join(
            f'pk5d_mingw_toolchain_check("{g.as_posix()}" ok why)\n'
            'message(STATUS "OK=${ok} ${why}")\n'
            for g in (good, bad)
        ),
        newline="\n",
    )
    out = subprocess.run(
        [exe, "-P", str(script)], cwd=tmp_path, capture_output=True, text=True, check=True
    ).stdout
    assert re.findall(r"^-- OK=(.*)$", out, re.M) == [
        "TRUE MinGW-builds' GCC 12.2.0 with the UCRT",
        f"FALSE {bad.as_posix()} is MinGW-builds' 12.2.0 for msvcrt.dll, not the UCRT",
    ]


# --- tools/get_mingw.sh and the check name the same toolchain ---


def test_get_mingw_fetches_the_toolchain_the_check_takes():
    sh = GET_MINGW.read_text(encoding="utf-8")
    cm = MODULE.read_text(encoding="utf-8")
    url = re.search(r"^url=(\S+)$", sh, re.M).group(1)
    version, rt, rev = re.search(r"/download/([\d.]+)-(rt_v\d+)-(rev\d+)/", url).groups()
    assert url.endswith(f"/x86_64-{version}-release-posix-seh-ucrt-{rt}-{rev}.7z")
    # the check's banner (CMake's regex, its escapes undone) is that build's, and no other's
    rx = re.search(r'_banner MATCHES "([^"]+)"', cm).group(1).replace("\\\\", "\\")
    assert re.search(rx, f"gcc.exe (x86_64-posix-seh-{rev}, Built by MinGW-W64 project) {version}")
    assert re.search(rx, PINNED)
    assert not re.search(rx, RUNNER) and not re.search(rx, MSYS2)
    # the default folder, where the build looks
    assert f"dest=$base/{GET_MINGW_DIR}" in sh
    # the compiler's programs find their libwinpthread-1.dll off PATH: a copy next to them
    libexec = f"$dest/mingw64/libexec/gcc/x86_64-w64-mingw32/{version}/"
    assert f'cp "$dest/mingw64/bin/libwinpthread-1.dll" "{libexec}"' in sh
    assert f'"$ENV{{LOCALAPPDATA}}/{GET_MINGW_DIR}"' in cm
