# Releasing petakit5d

The plan for the first release (Phase 1) and the checklist for every release: the version, the
wheels, the licenses of everything the wheels bundle or require, and the move of the package into
a repository of its own. What waits for a decision of the user is marked **(user)**.

## Version and metadata

- **Scheme:** PEP 440. The release candidates of Phase 1 are `0.1.0rc1`, `0.1.0rc2`, ...; the
  release is `0.1.0`. The version lives in `src/petakit5d/_version.py` alone (`pyproject.toml`
  reads it), so a release changes that one line.
- **The exact FFT add-on** `petakit5d-fftw` has a version of its own, FFTW's version and a
  packaging number (`3.3.8.0`), in `cpp/fft/exact/pyproject.toml`. It changes only when FFTW or its
  build changes, not with every petakit5d release.
- **Metadata** (`pyproject.toml`): license `BSD-3-Clause` with `LICENSE`,
  `THIRD_PARTY_NOTICES.md` and pocketfft's license in the wheel's `licenses/`; Python 3.11 and
  newer; the classifiers of the tested platforms and Pythons. `[project.urls]` (home page,
  documentation, issues) is added with the new repository **(user: its name and owner)**.
- **Dependencies:** NumPy, SciPy, scikit-image, h5py, cpp-zarr and cpp-tiff. The floors of
  cpp-zarr and cpp-tiff move to their next release (1.8.0) once it is published, in both
  `[project] dependencies` and `[build-system] requires` (the Linux build links cpp-zarr's
  OpenMP runtime, below).

## Wheels

### The matrix

| Platform | Wheel tag | Build | Repair | OpenMP runtime |
|---|---|---|---|---|
| Linux x86-64 | `manylinux_2_28_x86_64` (auditwheel also tags it `manylinux_2_24`) | GCC 14 of the image `quay.io/pypa/manylinux_2_28_x86_64:2026.06.04-1`, the image of the cpp-zarr and cpp-tiff wheels | auditwheel | the cpp-zarr wheel's libgomp, carried in the wheel (below) |
| macOS arm64 | `macosx_13_0_arm64` | Apple Clang | delocate | LLVM's libomp, built for macOS 13.0 from LLVM's sources and copied in by delocate |
| Windows x86-64 | `win_amd64` | MSVC (on the windows-2025 runners Visual Studio 2026's, toolset 14.51), `/openmp:llvm` | delvewheel `--no-diagnostic` | `libomp140.x86_64.dll`, which delvewheel copies in with the C++ runtime |

For each platform two wheels: **CPython 3.11** (`cp311-cp311`) and **CPython 3.12 and newer**
(`cp312-abi3`), with nanobind's stable ABI. The modules use no Python API outside the limited API
(a capsule is the only direct call), so one abi3 wheel serves 3.12, 3.13, 3.14 and later
versions. Free-threaded Pythons (`cp314t`) get no wheel: they have no stable ABI.
`[tool.cibuildwheel]` in `pyproject.toml` holds the matrix; `cibuildwheel` builds it on each
platform (the workflow below), or, on Linux machines with Apptainer only, the same recipe that
builds the cpp-zarr and cpp-tiff wheels (a script that runs cibuildwheel's steps in the same
image).

Each wheel is tested where it is built, in a fresh environment, by cibuildwheel's test command:
the smoke tests and the OpenMP runtime tests (`tests/test_smoke.py`,
`tests/test_openmp_runtime.py`); the abi3 wheel also on CPython 3.14 (the `cp314-*` entry of the
matrix, for which cibuildwheel tests the abi3 wheel instead of building one). The full suite runs
against an installed wheel in CI (below).

### The OpenMP runtime

- **Linux:** the modules link the libgomp that the cpp-zarr wheel bundles
  (`libgomp-e985bcbb.so.1.0.0`, the image's libgomp renamed by auditwheel) and carry a copy in
  `petakit5d/_ext/` with the RUNPATH `$ORIGIN`, so that a process has one OpenMP runtime with
  cpp-zarr and cpp-tiff ([decisions/0008](decisions/0008-openmp-runtime.md)). auditwheel leaves
  the copy in place and bundles nothing else: the modules need only the C library, libm,
  libpthread, libdl and the loader (libstdc++ and libgcc are linked statically and hidden). The
  name depends on the image, so the image is pinned in all three projects and changes in all three
  together.
- **macOS:** LLVM's libomp (21.1.8), built by `tools/build_libomp.sh` (cibuildwheel's
  before-all) from LLVM's release sources, fetched by URL and SHA-256, for the wheel's deployment
  target, macOS 13.0, the target of the cpp-zarr and cpp-tiff wheels; delocate copies it into
  `petakit5d/.dylibs/`. Homebrew's libomp would not do: a bottle is built for the machine's macOS
  (15.0 on the macOS 15 runners, and Homebrew has none for macOS 13 or 14), so the wheel's target
  would follow it, and a `macosx_15_0` wheel installs on no older macOS, though every dependency
  supports macOS 13.
- **Windows:** MSVC's LLVM OpenMP runtime `libomp140.x86_64.dll`. The build copies it next to
  the modules when the toolset's redistributable folders have it; Visual Studio 2026's (the
  windows-2025 runners') don't, so delvewheel copies it into `petakit5d.libs/` under a mangled
  name, with the C++ runtime (`msvcp140.dll`), from the first folder on `PATH` that has them:
  on the runners System32's libomp140, and for the C++ runtime, unless the workflow puts the
  toolset's own redistributable one first, any program's copy that comes earlier (a JDK's 14.40
  in the first trial, older than the 14.51 toolset, which Microsoft does not support and
  delvewheel warns of). All modules load the one copy of each.

### The GPU extra

`pip install "petakit5d[gpu]"` adds CuPy for CUDA 12, CUDA 12.9's runtime and NVRTC and cuFFT
11.0.8.103, MATLAB R2024b's, with which the GPU FFTs give MATLAB's bits
([decisions/0003](decisions/0003-gpu.md)). The CUDA kernels are compiled at run time from the
sources in `petakit5d/_ext/cuda/`, so the wheels carry no GPU binaries, and without the extra
every CPU path works. The cuFFT pin can conflict with other packages that pin their own cuFFT
(PyTorch's CUDA wheels): an environment with both needs one of them installed with `--no-deps`.

### The exact FFT add-on (GPL)

`cpp/fft/exact` is the separate package `petakit5d-fftw` (FFTW 3.3.8, GPL-2.0-or-later). Its
wheel holds two C libraries and no Python extension, so one `py3-none-<platform>` wheel per
platform serves every Python. It is built, published and installed apart from petakit5d. Nothing
of it is in petakit5d's wheels or source distribution: the checks below find no FFTW symbol in the
modules.

### Checks of every wheel before it is published

`tools/check_wheel.py` makes them, the same on every platform: it reads the binaries of all
three formats (ELF, Mach-O, PE) with the standard library alone, so any wheel can be checked on
any machine.

```
python tools/check_wheel.py [--forbid REGEX]... [--forbid-file FILE] [-q] <wheels, sdists or folders>
```

It prints each wheel's inventory (its bundled libraries, and the libraries and search paths of
its modules), then every finding, and exits with 0 if all are clean, 1 on a finding and 2 for an
archive it can't read. The distribution comes from the file name: petakit5d or petakit5d_fftw.

- **The tags:** the WHEEL file has the file name's tags; every binary is of the platform's format
  and architecture; the modules are built for the tags (abi3: `<group>.abi3.so`, or a `.pyd` that
  loads `python3.dll`; cp311: that version's suffix and `python311.dll`); no Mach-O file has a
  newer deployment target than the wheel's tag.
- **The bundled libraries:** the libraries besides the modules are the expected ones: Linux, the
  libgomp copy next to the modules alone; macOS, `.dylibs/libomp.dylib`; Windows, one
  `libomp140.x86_64.dll` (next to the modules, or delvewheel's copy in `petakit5d.libs/`) and
  delvewheel's copies of `msvcp140*.dll`; the add-on, its two FFTW libraries. Every library a
  binary loads is in the wheel (found through its search path, or on Windows next to it or in
  `petakit5d.libs/`) or the system's: the C library, macOS's `/usr/lib` and `/System`, Windows'
  system DLLs, API sets and the Python DLL and C runtime every CPython ships. All modules load
  the same OpenMP runtime, so that a process has one.
- **The search paths:** every RUNPATH (Linux) is `$ORIGIN`, with no old-style RPATH; every
  LC_RPATH (macOS) starts with `@loader_path`; no install name is a folder outside the wheel.
- **No folder of the build machine:** every member's bytes (with ELF's compressed sections
  decompressed) and its name are searched, the approach of cpp-tiff's `check_wheel_paths.py`,
  for home and user folders, temporary build environments and Python's temporary folders,
  macOS's per-user temporary folders, cibuildwheel's folders and container, the CI runners'
  workspace and tool cache, Homebrew's prefix and Spack's stages, and `--forbid` patterns: a
  site's own names, which are kept outside the repository. Microsoft's redistributable DLLs
  (delvewheel's copies of `libomp140` and `msvcp140*`, by name) name the folders of Microsoft's
  own builds (`D:\a\_work\1\s\...`): when Microsoft's Authenticode certificate is on them, they
  get the `--forbid` patterns only. No compiled Python file is carried.
- **No FFTW** (GPL) in petakit5d's binaries: no symbol `fftw*_` defined, no FFTW version string.
- **The license:** METADATA's `License-Expression` (`BSD-3-Clause`; the add-on's
  `GPL-2.0-or-later`), and every `License-File`, the project's three among them (the add-on's
  `COPYING`), in `.dist-info/licenses/`.

For the source distribution: no `cpp/fft/exact/`, no build folders, compiled files or binaries,
no owner names in the tar headers, the license as above, and the folder check of the member
names and of PKG-INFO. Its other files are the repository's sources, where a build folder can be
named for a reason (the libomp search hints of `cmake/Pk5dOpenMP.cmake`, the check's own patterns
and test cases), so their contents get the `--forbid` patterns only.

### The wheel workflow

`.github/workflows/wheels.yml` builds all of it and keeps it as the run's artifacts; it
publishes nothing, uploads nothing anywhere else and uses no secret. It runs for a tag `v*`, for
a push to the branch `wheels-trial` (trial runs) and by hand (`workflow_dispatch`, for the new
repository's default branch later). The `test` branch's pushes don't run it.

| Job | Runners | Steps |
|---|---|---|
| `wheels` | ubuntu-24.04, macos-15, windows-2025 | cibuildwheel 4.3.5, run from `python/` with `[tool.cibuildwheel]`: per Python the build, the repair, the audit of the abi3 wheel (`abi3audit --strict`) and the test command in a fresh environment; then the checks |
| `sdist` | ubuntu-24.04 | `python -m build --sdist`; a wheel built from the sdist alone (`pip wheel`, isolated), installed in a fresh environment and tested by the sdist's own smoke and OpenMP tests; the checks of both |
| `addon` | the three | cibuildwheel from `python/cpp/fft/exact`: one `py3-none` wheel per platform, built once (CPython 3.12; on Linux in petakit5d's image), its test (both libraries load and are `fftw-3.3.8`) and the checks |

The artifacts are `petakit5d-<system>-<architecture>`, `sdist` and
`petakit5d-fftw-<system>-<architecture>`. To check wheels anywhere else, run
`python tools/check_wheel.py <folder>` from `python/`; with the site's own names
(`--forbid-file`, a file outside the repository) on the lab's machines.

### First results (2026-11-19)

- **Linux x86-64, on the cluster** (Apptainer, the pinned image, the recipe of the cpp-zarr and
  cpp-tiff wheels, the runner-like Ubuntu 24.04 for the sdist): the `cp311` and `cp312-abi3`
  wheels passed their test commands (27 passed each; the abi3 wheel also on CPython 3.14.5,
  27 passed) and `abi3audit --strict`, and the checks found nothing, with the site's names too:
  one bundled library, `libgomp-e985bcbb.so.1.0.0`; every module's RUNPATH `$ORIGIN`; no build
  folder; no FFTW symbol. The sdist has 836 files and no `cpp/fft/exact/`, and the wheel built
  from it passed the checks. The add-on's wheel passed its test and the checks after its build
  left FFTW's source folder out of the file names its assertions carry (they named the build's
  temporary folder).
- **macOS arm64, natively on a macOS 13.7 machine** (cibuildwheel 4.3.5 with python.org's
  CPythons): libomp 21.1.8 built for macOS 13.0; the `cp311` and `cp312-abi3` wheels
  (`macosx_13_0_arm64`) passed their test commands (15 passed, 12 skipped each; the abi3 wheel
  also on CPython 3.14), `abi3audit --strict` and the checks: one bundled library,
  `.dylibs/libomp.dylib`, which every module loads through `@loader_path`; no LC_RPATH left;
  every binary's deployment target 13.0; no build folder, after the core module's build
  information named libomp by its file name instead of its path. The add-on's
  `macosx_11_0_arm64` wheel passed its test and the checks.
- **GitHub, the first trial run** (all three systems, the workflow as above): the sdist, the
  Linux and macOS wheels and the add-on's Linux and macOS wheels passed every step. The Windows
  `cp311` and `cp312-abi3` wheels were built, repaired, audited and passed their test commands
  (the abi3 wheel also on CPython 3.14); delvewheel had copied in `libomp140.x86_64.dll` (Visual
  Studio 2026 has none in its redistributable folders) and the system's `msvcp140.dll` 14.40,
  older than the toolset (14.51). The checks then expected the build's placement of libomp140
  and took the folders of Microsoft's own builds in those two DLLs for the runner's; both are
  handled as described above. The add-on's Windows job had run its steps in PowerShell.
- **GitHub, the second trial run:** every job passed but the add-on's Windows wheel. The Windows
  `cp311` and `cp312-abi3` wheels passed their test commands (15 passed, 12 skipped each; the
  abi3 wheel also on CPython 3.14), `abi3audit --strict` and the checks: libomp140 (System32's)
  and `msvcp140.dll` 14.51.36247 (the toolset's) bundled by delvewheel, Microsoft-signed, the one
  OpenMP runtime of all modules. The add-on's Windows wheel passed its test (`fftw-3.3.8-sse2-avx`
  from both libraries; FFTW's file names trimmed by MSVC too), but delvewheel had copied in
  `vcruntime140.dll`, which a `py3` wheel's delvewheel doesn't know CPython ships: it is excluded
  now.
- **GitHub, the third trial run:** every job passed, on all three systems: the wheels' test
  commands (Linux 27 passed for each Python; macOS and Windows 15 passed, 12 skipped), the abi3
  audits, the sdist with the wheel built from it, the add-on's three wheels, and the checks of
  every wheel and the sdist.
- Before: the abi3 wheel, installed into a fresh environment on a machine of another CPU
  generation than the build's, passed the full CPU suite with and without the golden store and
  with the add-on's wheel installed, and installed with its dependencies from PyPI on CPython
  3.14 it passed the FFT, core, deskew and decon tests.

## Licenses

petakit5d is BSD-3-Clause. **Nothing GPL-licensed is in its wheels**; the GCC runtime parts that
are linked or carried are under the GPL with the GCC Runtime Library Exception, which allows
shipping them with code under any license.

### In the wheels

| Component | Where | License |
|---|---|---|
| petakit5d (Python, C++, CUDA sources) | everything | BSD-3-Clause |
| pocketfft | compiled into the `fft` module | BSD-3-Clause |
| nanobind and its `tsl::robin_map` | compiled into every module | BSD-3-Clause; MIT |
| fdlibm and FreeBSD msun algorithms | `matlabcompat.elementary`, the `decon` module | the Sun notice ("freely granted, provided that this notice is preserved"), kept in the sources |
| normxcorr3 (method of J. P. Lewis) | `kernels.xcorr`, the `xcorr` module | BSD (`THIRD_PARTY_NOTICES.md`) |
| libstdc++, libgcc (Linux) | linked statically | GPL-3.0 with the GCC Runtime Library Exception |
| libgomp (Linux) | `petakit5d/_ext/libgomp-e985bcbb.so.1.0.0` | GPL-3.0 with the GCC Runtime Library Exception |
| libomp (macOS) | `petakit5d/.dylibs/` | Apache-2.0 with LLVM exceptions |
| libomp140, MSVC C++ runtime (Windows) | `petakit5d.libs/` (or `petakit5d/_ext/`) | Microsoft's redistributable terms (the Visual Studio that builds the wheel: 2026 on the windows-2025 runners) |

### Required packages

| Package | License | Bundles (their wheels) |
|---|---|---|
| NumPy | BSD-3-Clause AND 0BSD AND MIT AND Zlib AND CC0-1.0 | OpenBLAS (BSD-3-Clause) |
| SciPy | BSD-3-Clause | OpenBLAS (BSD-3-Clause), libgfortran (GPL-3.0 with the GCC Runtime Library Exception), libquadmath (LGPL-2.1) |
| scikit-image | BSD-3-Clause | - |
| networkx, tifffile, lazy-loader (scikit-image's) | BSD-3-Clause | - |
| imageio (scikit-image's) | BSD-2-Clause | - |
| pillow (scikit-image's) | MIT-CMU | libjpeg-turbo, libpng, libtiff, libwebp, libavif, freetype, harfbuzz, brotli, lcms2, zlib, libXau, ... (permissive) |
| packaging (scikit-image's) | Apache-2.0 OR BSD-2-Clause | - |
| h5py | BSD-3-Clause | HDF5 (BSD-3-Clause-like), libaec and libsz (BSD-2-Clause) |
| cpp-zarr, cpp-tiff | BSD-3-Clause (their repositories' `LICENSE`; the 1.7.3 wheels declare no license in their metadata) | their library, with its compression libraries linked in; libgomp (GPL-3.0 with the GCC Runtime Library Exception) |

### Optional extras

| Extra | Packages | Licenses |
|---|---|---|
| `gpu` | cupy-cuda12x, cuda-pathfinder; nvidia-cuda-runtime-cu12, nvidia-cuda-nvrtc-cu12, nvidia-cufft-cu12 | MIT; Apache-2.0; NVIDIA's proprietary license (the CUDA EULA: not open source, installed by the user from PyPI, never bundled) |
| `figures` | matplotlib | PSF-based (matplotlib's own) |
| `n5` | numcodecs | MIT |
| `formats` | nd2, czifile, imagecodecs | BSD-3-Clause |
| `basic` | basicpy (with `--no-deps`), torch, torch-dct, hyperactive, pooch, pydantic, dask, tqdm | MIT; BSD-3-Clause; MIT; MIT; BSD-3-Clause; MIT; BSD-3-Clause; MPL-2.0 and MIT |
| `test`, `dev` | pytest, pytest-timeout, ruff | MIT |
| (add-on) | petakit5d-fftw | **GPL-2.0-or-later**: a package of its own, never a dependency or an extra of petakit5d |

The table is checked at each release from the installed packages' metadata (`License-Expression`,
the license classifiers) and their bundled libraries (`*.libs/`, `.dylibs/`).

## The new repository

The package moves into a repository of its own, with its history. PetaKit5D (MATLAB) stays as it
is: it remains the reference that the oracle cases were made from.

### What travels, and how

| What | How |
|---|---|
| `python/` with its history | `git filter-repo` (below): every commit that touched it, with its authors, dates and messages; `python/` becomes the root |
| the CI workflow (`.github/workflows/python.yml`) | in the same filter, with its history; then one commit for the new layout: no `working-directory: python`, no `paths:` filters, the Slurm config's path without `python/` |
| the wheel workflow (`.github/workflows/wheels.yml`) | in the same filter, with its history; in the layout commit: no `working-directory: python` (the add-on's `python/cpp/fft/exact` becomes `cpp/fft/exact`; its `../../../tools/check_wheel.py` stays), and the artifacts' paths without `python/` |
| the oracle manifests (`tests/**/*.oracle.json`) | in `python/`, so with the history: they carry the MATLAB references that CI compares with |
| the oracle harness (`tests/oracle_matlab`: the MATLAB side and `ORACLE_COMMIT`) | in `python/`. Making new references needs a MATLAB PetaKit5D tree at the commit `ORACLE_COMMIT` names (`PETAKIT5D_ORACLE_DIR`): a checkout of the MATLAB repository, which stays public. The docs name that repository and the commit |
| the golden store (the full MATLAB outputs, `PETAKIT5D_GOLDEN_DIR`) | not in git; CI does not need it. kept local; the user's decision of 2026-11-19 (published as a versioned archive keyed by the oracle commit only if another site needs it) |
| the demo data (`PETAKIT5D_DEMO_DATA`) | not in git, as now |
| docs, examples, benchmarks, the lock files, the MATLAB test map | in `python/` |
| the exact FFT add-on (`cpp/fft/exact`) | in `python/`, as a package with its own license. it stays a folder of the new repository (the user, 2026-11-19) |
| paths that say `python/` (`pip install ./python`, `python/docs/...`) | one commit in the new repository |

### What stays in PetaKit5D

The MATLAB code, its tests (`tests/`), demos, `mcc/`, `setup.m`, `third_parties/` and the MATLAB
workflow (`tests.yml`). The compiled MATLAB package is retired after the port's first release;
`pk5d-master` and `pk5d mcc-shim` replace it for PetaKit5D-GUI and PyPetaKit5D. PetaKit5D's README gets no pointer to the new repository for now (the user, 2026-11-19), and
`python/`, `python.yml` and `wheels.yml` stay in PetaKit5D until the user decides otherwise.

### The steps (rehearsed on a throw-away clone)

```
git clone --no-local --single-branch --branch <branch> <PetaKit5D> petakit5d-split
cd petakit5d-split
git filter-repo --path python/ --path .github/workflows/python.yml \
    --path .github/workflows/wheels.yml --path-rename python/:
# check: the commit count, the authors, the trees (below)
# then, in the new repository: the workflows' layout commit, the paths commit, the remote
git remote add origin <the new repository>
```

The checks after the filter:

- the number of commits equals
  `git rev-list --count <branch> -- python .github/workflows/python.yml .github/workflows/wheels.yml`
  in PetaKit5D;
- with `--subdirectory-filter python` instead, every commit's tree is the `python/` tree of its
  original, one for one (`git rev-parse <commit>:python` against the new commit's tree);
- authors, committers and dates are unchanged; messages too, except the abbreviated hashes of
  kept commits, which filter-repo rewrites to the new hashes (hashes of MATLAB commits stay);
- the confidentiality scan over the whole history (every blob and message), not only the last
  tree.

## The release checklist

1. **The libraries (user):** cpp-tiff 1.8.0 and cpp-zarr 1.8.0 are published, from the pinned
   image. The user publishes them once the port's work is done. The candidates (cpp-tiff's build 6
   and cpp-zarr's release build) carry their license metadata; their Linux and macOS code is the
   code the port's suites passed with. Raise the floors to 1.8.0 (dependencies and build requirements), check that their Linux
   wheels still bundle `libgomp-e985bcbb.so.1.0.0`, and rerun the OpenMP runtime tests against
   them. On macOS this also ends the known crash of the 1.7.3 wheels.
2. **The version:** `src/petakit5d/_version.py` to `0.1.0rcN` (then `0.1.0`); the add-on's only
   if FFTW or its build changed.
3. **The suites, green:** the full CPU suite with and without the golden store, the exact add-on's
   suite, the GPU tests on a Volta, an A100-class and a Blackwell GPU, and the CI on Linux, macOS
   and Windows, all on the release commit.
4. **The wheels:** the wheel workflow's run for the release tag (or a `wheels-trial` run of the
   release commit) is green: every wheel of the matrix and of the add-on passes its test command
   and the checks; one Linux wheel installed into a fresh environment on another machine passes
   the full CPU suite; the wheels that are published are that run's artifacts.
5. **The source distribution:** the same run's `sdist` job: no `cpp/fft/exact/`, no build
   folders, and a wheel built from it alone passes the smoke tests.
6. **The licenses:** the tables above against the installed packages' metadata; new bundled
   libraries go into `THIRD_PARTY_NOTICES.md`.
7. **The new repository (user):** its name and owner (to be decided later); the split as above; the CI in its new
   layout green there; `[project.urls]` filled in.
8. **Publishing:** PyPI (petakit5d, and petakit5d-fftw apart from it), from CI with trusted
   publishing, first to TestPyPI, after the libraries (agreed with the user, 2026-11-19); the
   history's authors and their e-mail addresses stay as they are, and so does `[project] authors`; the tag `v0.1.0rcN` on the release commit.
9. **The switch (user):** PetaKit5D-GUI and PyPetaKit5D move to `pk5d mcc-shim` / `pk5d-master`
   ([migration.md](migration.md)); the compiled MATLAB package is retired after the first release.
10. **The old repository (user):** no pointer in PetaKit5D's README for now; whether PetaKit5D
    still takes the 1.8.0 vendoring of cpp-tiff and cpp-zarr and its Windows demos.
