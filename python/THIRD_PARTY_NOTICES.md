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
Normalized Cross-Correlation", Vision Interface 1995. File Exchange submissions without a linked
repository are distributed under its BSD license:

```
Copyright (c) 2020, Daniel Eaton
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```
