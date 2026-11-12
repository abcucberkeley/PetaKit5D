# pocketfft (C++ header)

`pocketfft_hdronly.h` and `LICENSE.md` unchanged from pocketfft's `cpp` branch,
<https://github.com/mreineck/pocketfft>, commit `4d865cba24359b05b37da4e75aee24498c3b78e7`
(2026-10-18; tarball SHA-256 `db236e44b1020dd0b415e5040ce10d8be604fefb7ce83ec0b32bc0f78975aa2f`).
License: BSD 3-Clause (`LICENSE.md`, and the header's own notice). NumPy bundles the same library.

The port's default CPU FFT backend ("permissive", `cpp/fft/pocketfft_api.cpp`, `cpp/fft/halfspec.cpp`).
To update: replace both files from a newer commit, note it here, and run `tests/fft` and
`tests/matlabcompat/test_fft.py`.
