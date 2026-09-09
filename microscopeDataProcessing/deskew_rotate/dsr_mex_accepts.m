function ok = dsr_mex_accepts(T)
% True if the 0-based backward map T (4x4, mex axis order: T(i,:) gives input axis i from
% [out1 out2 out3 1]) lies in the affine class the deskew/rotate warp kernels implement.
%
% The kernels (volume_deskew_rotate_warp_mex.cpp, skewed_space_interp_volume_deskew_rotate_warp_mex.cpp)
% do not apply a general affine. They copy input dim 1 to output dim 1 with unit stride (no
% matrix entry read; the only offset is the bbox start), compute input dim 2 from output
% dim 3 alone, and input dim 3 from output dims 2 and 3. Any other coefficient is silently
% dropped, so a transform outside this class must go through imwarp.
%
% Sample-scan deskew/rotate has this form (the shear and the rotation are both in the
% x-z plane and the in-sheet coordinate maps to height only), with or without an x or z
% resampling factor. A y factor or an objective-scan rotation does not.
%
% T(1,4) must be 0: the kernel takes its dim-1 offset from the bbox start it is given,
% not from the matrix, so a genuine dim-1 translation cannot be expressed to it.
%
% tol: structurally-zero entries come out ~1e-16 from the matrix products; live entries
% are O(0.1-1) (sin/cos of the skew angle, resampling factors).
%
% Author: Martin Alvarez-Kuglen (09/08/2026)

validateattributes(T, {'numeric'}, {'size', [4, 4], 'finite'}, mfilename, 'T');
tol = 1e-9;
ok = all(abs(T(1, :) - [1 0 0 0]) < tol) ...
    && abs(T(2, 1)) < tol && abs(T(2, 2)) < tol ...
    && abs(T(3, 1)) < tol ...
    && all(abs(T(4, :) - [0 0 0 1]) < tol);
end
