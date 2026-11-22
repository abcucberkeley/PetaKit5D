#!/bin/bash
# The toolchain of the Windows build with MinGW-w64 (docs/decisions/0011-windows-toolchain.md): the
# MinGW-builds project's GCC 12.2.0 (UCRT, posix threads, SEH; binutils 2.39), the toolchain of the
# release wheels and of the cpp-zarr and cpp-tiff Windows wheels, from its release by URL and SHA-256
# (the release publishes no digest), unpacked into <folder>/mingw64. Needs curl, sha256sum and 7z
# (GitHub's Windows runners, Git for Windows' bash with 7-Zip).
#
#   bash tools/get_mingw.sh [<folder>]
#
# Without <folder>: %LOCALAPPDATA%\petakit5d\mingw-12.2.0, where a source build looks for the
# toolchain when the first gcc on PATH is not it (cmake/Pk5dToolchain.cmake; docs/platforms.md).
# Prints the compiler's version and, last, its bin folder, to put first on PATH where wanted
# (on GitHub: >> "$GITHUB_PATH").
set -euo pipefail
dest=${1:-}
if [ -z "$dest" ]; then
  if [ -z "${LOCALAPPDATA:-}" ]; then
    echo "get_mingw.sh: no <folder>, and no LOCALAPPDATA for the default one" >&2
    exit 2
  fi
  base=$LOCALAPPDATA
  if command -v cygpath > /dev/null; then
    base=$(cygpath -m "$base")
  fi
  dest=$base/petakit5d/mingw-12.2.0
fi
url=https://github.com/niXman/mingw-builds-binaries/releases/download/12.2.0-rt_v10-rev2/x86_64-12.2.0-release-posix-seh-ucrt-rt_v10-rev2.7z
sha=dc1f299e8450c96f59f4b8f7c1290ebb78ee00ddb71f9ed4ca0efaf5462f8c4e
mkdir -p "$dest"
curl -sSL --fail --retry 3 -o "$dest/mingw.7z" "$url"
got=$(sha256sum < "$dest/mingw.7z" | cut -d' ' -f1)
if [ "$got" != "$sha" ]; then
  echo "get_mingw.sh: the archive's SHA-256 is $got, not $sha" >&2
  exit 1
fi
7z x -y -o"$dest" "$dest/mingw.7z" > /dev/null
rm -f "$dest/mingw.7z"
bin="$dest/mingw64/bin"
if command -v cygpath > /dev/null; then
  bin=$(cygpath -m "$bin")
fi
"$bin/g++" --version | head -n 1
echo "$bin"
