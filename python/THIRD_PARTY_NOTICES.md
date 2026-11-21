# Third-party notices

PetaKit5D's Python package is BSD 3-Clause (`LICENSE`). Parts of it derive from, or bundle, the
works below; their notices follow, as their licenses require.

## pocketfft

The CPU FFTs' default backend, bundled in `cpp/fft/third_party/pocketfft`. Its license (BSD
3-Clause) is `cpp/fft/third_party/pocketfft/LICENSE.md`, which is included with the package.

## normxcorr3 (Daniel Eaton)

`petakit5d.kernels.xcorr` (`python/src/petakit5d/kernels/xcorr.py`, `python/cpp/xcorr/`) ports
PetaKit5D's `normxcorr3_fast`. That function is the lab's adaptation of Daniel Eaton's
`normxcorr3` ("normxcorr3 (Fast 3D NCC)", MATLAB File Exchange submission 73946, version 1.0.1,
2020), which computes the normalized cross-correlation by the method of J. P. Lewis, "Fast
Normalized Cross-Correlation", Vision Interface 1995. Its BSD license, as distributed with
PetaKit5D's copy of `normxcorr3` (`third_parties/normxcorr3/license.txt` in the MATLAB
repository):

```
Copyright (c) 2005, Daniel Eaton
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution
* Neither the name of  nor the names of its
  contributors may be used to endorse or promote products derived from this
  software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

## nanobind and tsl::robin_map

The extension modules are built with nanobind, whose support code is compiled into every module,
together with the hash table tsl::robin_map that nanobind carries (`ext/robin_map`). nanobind's
license (BSD 3-Clause):

```
Copyright (c) 2022 Wenzel Jakob <wenzel.jakob@epfl.ch>, All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software
   without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

tsl::robin_map's license (MIT):

```
MIT License

Copyright (c) 2017 Thibaut Goetghebuer-Planchon <tessil@gmx.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## fdlibm and FreeBSD's msun

`petakit5d.matlabcompat.elementary`, `petakit5d.matlabcompat.special` (`erf` and `erfc`, fdlibm's
`s_erf.c`), the decon module (`cpp/decon/omw.cpp`), the detection module's exponential
(`cpp/detection/module.cpp`) and the resize module's power (`cpp/resize/module.cpp`) follow the
algorithms and constants of fdlibm 5.3 and FreeBSD's msun, as their comments name them, to give
MATLAB's elementary functions bit for bit. Their notice, kept with that code:

```
Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.

Developed at SunSoft, a Sun Microsystems, Inc. business.
Permission to use, copy, modify, and distribute this
software is freely granted, provided that this notice
is preserved.
```

fdlibm 5.3's `k_tan.c`, `s_tan.c`, `e_exp.c` and `e_pow.c` carry the same notice with
"Copyright (C) 1993, 2004".

## MINPACK

The detection module's Levenberg-Marquardt fit (`cpp/detection/mexfit.hpp`: enorm, qrfac, qrsolv,
lmpar and lmder's iteration) follows MINPACK-1 (Moré, Garbow and Hillstrom, 1980). Its notice,
kept with that code:

```
Minpack Copyright Notice (1999) University of Chicago. All rights reserved
Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of
conditions and the following disclaimer in the documentation and/or other materials provided
with the distribution.
3. The end-user documentation included with the redistribution, if any, must include the
following acknowledgment: "This product includes software developed by the University of
Chicago, as Operator of Argonne National Laboratory." Alternately, this acknowledgment may
appear in the software itself, if and wherever such third-party acknowledgments normally
appear.
4. WARRANTY DISCLAIMER. THE SOFTWARE IS SUPPLIED "AS IS" WITHOUT WARRANTY OF ANY KIND. THE
COPYRIGHT HOLDER, THE UNITED STATES, THE UNITED STATES DEPARTMENT OF ENERGY, AND THEIR
EMPLOYEES: (1) DISCLAIM ANY WARRANTIES, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO ANY
IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE OR
NON-INFRINGEMENT, (2) DO NOT ASSUME ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY,
COMPLETENESS, OR USEFULNESS OF THE SOFTWARE, (3) DO NOT REPRESENT THAT USE OF THE SOFTWARE
WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS, (4) DO NOT WARRANT THAT THE SOFTWARE WILL FUNCTION
UNINTERRUPTED, THAT IT IS ERROR-FREE OR THAT ANY ERRORS WILL BE CORRECTED.
5. LIMITATION OF LIABILITY. IN NO EVENT WILL THE COPYRIGHT HOLDER, THE UNITED STATES, THE
UNITED STATES DEPARTMENT OF ENERGY, OR THEIR EMPLOYEES: BE LIABLE FOR ANY INDIRECT,
INCIDENTAL, CONSEQUENTIAL, SPECIAL OR PUNITIVE DAMAGES OF ANY KIND OR NATURE, INCLUDING BUT
NOT LIMITED TO LOSS OF PROFITS OR LOSS OF DATA, FOR ANY REASON WHATSOEVER, WHETHER SUCH
LIABILITY IS ASSERTED ON THE BASIS OF CONTRACT, TORT (INCLUDING NEGLIGENCE OR STRICT
LIABILITY), OR OTHERWISE, EVEN IF ANY OF SAID PARTIES HAS BEEN WARNED OF THE POSSIBILITY OF
SUCH LOSS OR DAMAGES.
```

This product includes software developed by the University of Chicago, as Operator of Argonne
National Laboratory.

## The GCC runtime libraries

The Linux modules link GCC's libstdc++ and libgcc statically, and the Linux wheels carry GCC's
OpenMP runtime libgomp (`petakit5d/_ext/libgomp-e985bcbb.so.1.0.0`, the copy that the cpp-zarr
wheel bundles: [decisions/0008](docs/decisions/0008-openmp-runtime.md)). The Windows modules built
with MinGW-w64's GCC link libstdc++ and libgcc statically too, and their wheels carry GCC's
`libgomp-1.dll` and `libgcc_s_seh-1.dll` (`petakit5d/_ext/`, the files of the MinGW-builds GCC
12.2.0 toolchain that the cpp-zarr and cpp-tiff wheels carry:
[decisions/0011](docs/decisions/0011-windows-toolchain.md)). These libraries are under the GNU
General Public License, version 3, with the GCC Runtime Library Exception, version 3.1, which
allows programs compiled with GCC and linked with them to be distributed under licenses of their
own (https://www.gnu.org/licenses/gcc-exception-3.1.html). The Windows wheels built with
MinGW-w64 carry both texts as the toolchain ships them, copied from it by the build:
`petakit5d/_ext/licenses/gcc/COPYING3` and `petakit5d/_ext/licenses/gcc/COPYING.RUNTIME`. GCC's
sources are at https://gcc.gnu.org.

## MinGW-w64's runtime and winpthreads (Windows, built with MinGW-w64)

The Windows modules built with MinGW-w64's GCC contain MinGW-w64's runtime: its start-up code and,
linked statically, functions of its C library (libmingwex: among them the math functions that
`cpp/common/mingw_libm.c` wraps). Their wheels carry MinGW-w64's POSIX threads library winpthreads
(`petakit5d/_ext/libwinpthread-1.dll`), which GCC's libgomp and libstdc++ use there. Their
licenses: the mingw-w64 project's BSD-style license for the runtime, with the notices of its parts
under other permissive licenses, which binary distributions must include; MIT for winpthreads,
with the BSD notice of the parts it derives from. The wheels carry the notices as the toolchain
ships them, copied from it by the build (the texts cpp-zarr's repository ships for its wheels):
`petakit5d/_ext/licenses/mingw-w64/COPYING.MinGW-w64-runtime.txt` and
`petakit5d/_ext/licenses/winpthreads/COPYING`. Built with MinGW-builds' GCC 12.2.0 (the toolchain
of `tools/get_mingw.sh`), the build checks the four files' SHA-256 against that toolchain's
(`cmake/Pk5dExtension.cmake`), and the release checks find them in the wheels
(`tools/check_wheel.py`). The exact add-on built with MinGW-w64 (`petakit5d_fftw`, libgcc and the
runtime linked in) carries the GCC and runtime notices in `petakit5d_fftw/licenses/`. MinGW-w64's
sources are at https://www.mingw-w64.org.
