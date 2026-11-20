#!/bin/bash
# LLVM's OpenMP runtime (libomp) for the macOS wheels, built from LLVM's release sources for the
# wheels' deployment target: cibuildwheel's before-all ([tool.cibuildwheel.macos] in
# pyproject.toml), with MACOSX_DEPLOYMENT_TARGET and OpenMP_ROOT, where it is installed (the
# build's Pk5dOpenMP.cmake finds it there; delocate copies it into the wheel). Homebrew's libomp
# is built for the machine's macOS (15.0 on macos-15 runners), which would be the wheels' floor.
#
#   OpenMP_ROOT=<prefix> MACOSX_DEPLOYMENT_TARGET=13.0 bash tools/build_libomp.sh
#
# The sources are fetched by URL and SHA-256; the build folders are left out of the library's
# file names (its assertions name their source files). Nothing is done when <prefix> holds it.
set -euo pipefail
VERSION=21.1.8
SHA256_OPENMP=856b023748b41ac7b2c83fd8e9f765ff48a4df2fe6777d2811ef7c7ed8f2f977
SHA256_CMAKE=85735f20fd8c81ecb0a09abb0c267018475420e93b65050cc5b7634eab744de9
URL=https://github.com/llvm/llvm-project/releases/download/llvmorg-$VERSION
PREFIX=${OpenMP_ROOT:?set OpenMP_ROOT to the folder to install libomp in}
TARGET=${MACOSX_DEPLOYMENT_TARGET:?set MACOSX_DEPLOYMENT_TARGET}
if [ -f "$PREFIX/lib/libomp.dylib" ] && [ "$(cat "$PREFIX/TARGET" 2>/dev/null)" = "$VERSION $TARGET" ]; then
    echo "libomp $VERSION for macOS $TARGET is in $PREFIX"
    exit 0
fi
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
for part in openmp:$SHA256_OPENMP cmake:$SHA256_CMAKE; do
    name=${part%%:*}
    curl -fsSL --retry 3 -o "$T/$name.tar.xz" "$URL/$name-$VERSION.src.tar.xz"
    echo "${part#*:}  $T/$name.tar.xz" | shasum -a 256 -c -
    tar -xf "$T/$name.tar.xz" -C "$T"
done
# the standalone build takes LLVM's shared CMake modules from ../cmake
mv "$T/cmake-$VERSION.src" "$T/cmake"
map="-ffile-prefix-map=$T=llvm-$VERSION"
cmake -S "$T/openmp-$VERSION.src" -B "$T/build" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PREFIX" -DCMAKE_OSX_DEPLOYMENT_TARGET="$TARGET" \
    -DCMAKE_OSX_ARCHITECTURES="$(uname -m)" -DCMAKE_C_FLAGS="$map" -DCMAKE_CXX_FLAGS="$map" \
    -DLIBOMP_INSTALL_ALIASES=OFF -DOPENMP_ENABLE_LIBOMPTARGET=OFF
cmake --build "$T/build" --parallel "${CMAKE_BUILD_PARALLEL_LEVEL:-3}"
cmake --install "$T/build"
cp "$T/openmp-$VERSION.src/LICENSE.TXT" "$PREFIX/LICENSE.TXT"
echo "$VERSION $TARGET" > "$PREFIX/TARGET"
echo "libomp $VERSION for macOS $TARGET installed in $PREFIX"
