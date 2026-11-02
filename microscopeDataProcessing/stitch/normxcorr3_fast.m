function C = normxcorr3_fast(T, A, shape, bbox, fastFFT)
% C = normxcorr3(TEMPLATE, IMAGE, SHAPE, BBOX, FASTFFT)
%
%       TEMPLATE - type double, ndims==3, size <= size of image
%       IMAGE    - type double, ndims==3
%       SHAPE    - one of: 'valid', 'same', 'full'. same as conv2 shape parameter
%                  'full' by default
%       BBOX     - [s1, s2, s3, t1, t2, t3] (optional): only this part of the
%                  'full' output is computed and returned (SHAPE must be 'full')
%       FASTFFT  - true or false (default: false). With BBOX, let the mex use
%                  smaller FFTs (2, 3, 5, 7-smooth sizes without wrap-around
%                  in BBOX), which changes C in the last bits
%
%       C        - values in [-1,1]. size depends on SHAPE
%
% adapted from daniel eaton's version with support for 1d and 2d, and
% optimize performance with mex based image integral
%
% With BBOX, single TEMPLATE and IMAGE on Linux, normxcorr3_window_mex computes
% the correlation and the integral images only in BBOX, with the same FFTs as
% MATLAB's, so C is the same bit for bit (unless FASTFFT).
%
% Author: Xiongtao Ruan (04/27/2024)


persistent grain

if nargin<3
	shape = 'full';
end
if nargin < 4
    bbox = [];
end
if nargin < 5
    fastFFT = false;
end

if ndims(A)>3 || ndims(T)>3
	error('A and T must be no more than 3 dimensional matrices');
end

szT = size(T, 1 : 3);
szA = size(A, 1 : 3);

if any(szT>szA)
	error('template must be smaller than image');
end

pSzT = prod(szT);

szOut = szT + szA - 1;

% the mex takes real single data whose output has at least 2 elements along y
% and x, and a bbox inside the output
useMex = false;
if ~isempty(bbox)
    if ~strcmpi(shape, 'full')
        error('BBOX requires SHAPE ''full''');
    end
    bbox = bbox(:)';
    useMex = isunix && ~ismac && isa(T, 'single') && isa(A, 'single') && isreal(T) && isreal(A) ...
        && all(szOut(1 : 2) > 1) && numel(bbox) == 6 && all(bbox == round(bbox)) ...
        && all(bbox(1 : 3) >= 1) && all(bbox(1 : 3) <= bbox(4 : 6)) && all(bbox(4 : 6) <= szOut);
    if useMex && exist('normxcorr3_window_mex', 'file') ~= 3
        mex_fallback('normxcorr3_window_mex');
        useMex = false;
    end
end
if useMex
    try
        if isempty(grain)
            grain = matlab_fft_grain();
        end
        [corrTA, intImgA, intImgA2] = normxcorr3_window_mex('ncc', T, A, bbox, maxNumCompThreads, grain, fastFFT);
    catch ME
        mex_fallback('normxcorr3_window_mex', ME);
        useMex = false;
    end
end

% compute the numerator of the NCC
if ~useMex
    corrTA = real(ifftn(fftn(A,szOut).*fftn(T(end:-1:1, end:-1:1, end:-1:1), szOut)));
end

sumT = sum(T(:));
denomT = std(T(:));
% clear T;

% make the running-sum/integral-images of A and A^2, which are
% used to speed up the computation of the NCC denominator
if ~useMex
    intImgA = integral_image_3d(A,szT);
end

num = (corrTA - intImgA*sumT/pSzT ) / (pSzT-1);
clear corrTA; 

% compute the denominator of the NCC
if ~useMex
    intImgA2 = integral_image_3d(A.*A,szT);
end

denom = denomT * sqrt(max(intImgA2 - (intImgA.^2)/pSzT, 0) / (pSzT-1) );
clear intImgA intImgA2;

% compute the NCC
C = num ./ (denom + eps) .* (denom ~= 0);

if ~isempty(bbox)
    if ~useMex
        C = crop3d(C, bbox);
    end
    return;
end

switch( lower(shape) )
	case 'full'
	case 'same'
		szTp = fix((szT-1)/2);
		C = C( szTp(1)+1:szTp(1)+szA(1), szTp(2)+1:szTp(2)+szA(2), szTp(3)+1:szTp(3)+szA(3) );
	case 'valid'
		C = C(szT(1):end-szT(1)+1,szT(1):end-szT(2)+1,szT(3):end-szT(3)+1);
	otherwise
		error('unknown SHAPE %s, assuming FULL by default', shape);
end

end


function [grain] = matlab_fft_grain()
% MATLAB's FFTs transform the lines of each dimension in groups of four and
% some lines alone, which for some lengths changes the last bits. When MATLAB
% runs its FFTs in parallel (it decides that at its first FFT, from the number
% of threads then), it splits the lines into pieces of at most 32 lines and
% runs the last (length % 4) lines of each piece alone, otherwise only the last
% lines of a dimension. Find out which split this MATLAB session uses from the
% FFT of test data, so that the mex runs the same lines alone (grain 32: the
% split, 0: none). As in decon_lucy_rfft.

grains = [32, 0];
for n = [64, 128, 512]
    m = 2000;
    X = single(mod(reshape(0 : n * m - 1, n, m) * 0.6180339887, 1));
    Y = fft(X, [], 1);
    Y = Y(1 : n / 2 + 1, :);
    match = false(size(grains));
    for i = 1 : numel(grains)
        match(i) = isequal(normxcorr3_window_mex('fft1', X, maxNumCompThreads, grains(i)), Y);
    end
    if sum(match) == 1
        grain = grains(match);
        return;
    end
end
% the test data can't tell (or neither split matches): the parallel split
grain = grains(1);

end
