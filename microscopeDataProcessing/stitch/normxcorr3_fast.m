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
%       C        - values in [-1,1] (up to rounding), in the class of IMAGE.
%                  size depends on SHAPE
%
% adapted from daniel eaton's version with support for 1d and 2d, and
% optimize performance with mex based image integral
%
% The correlation is computed by FFTs in the class of the inputs, of the
% template and the image each less its mean; the window sums, the
% normalization and the NCC in double. Windows of the image whose variance is
% at most eps(class) times the image's (flat to the precision of the class)
% and a flat template give 0.
%
% With BBOX, single TEMPLATE and IMAGE on Linux, normxcorr3_window_mex computes
% the correlation and the window sums only in BBOX, with the same FFTs and sums
% as the MATLAB code, so C is the same bit for bit (unless FASTFFT).
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
    if ~(numel(bbox) == 6 && all(bbox == round(bbox)) && all(bbox(1 : 3) >= 1) ...
            && all(bbox(1 : 3) <= bbox(4 : 6)) && all(bbox(4 : 6) <= szOut))
        error('BBOX must hold integers [s1, s2, s3, t1, t2, t3] with 1 <= s <= t <= the size of the full output');
    end
    useMex = isunix && ~ismac && isa(T, 'single') && isa(A, 'single') && isreal(T) && isreal(A) ...
        && all(szOut(1 : 2) > 1);
    if useMex && exist('normxcorr3_window_mex', 'file') ~= 3
        mex_fallback('normxcorr3_window_mex');
        useMex = false;
    end
    win = bbox;
else
    win = [1, 1, 1, szOut];
end

% the template and the image, each less its mean (in their class): with the
% means, the correlation and the window sums keep a large constant part, whose
% rounding (and that of its removal) in single precision is larger than the
% NCC of faint or flat regions. The image stays zero-padded (its mean is
% accounted for where the template overlaps it).
cT = cast(mean(T(:), 'double'), 'like', T);
cA = cast(mean(A(:), 'double'), 'like', A);
T0 = T - cT;
A0 = A - cA;
% at the output position k (along each dimension) the template's elements
% loT .. hiT overlap the image's loA .. hiA
[loT, hiT, loA, hiA] = deal(cell(1, 3));
for d = 1 : 3
    k = win(d) : win(d + 3);
    loT{d} = max(1, szT(d) - k + 1);
    hiT{d} = min(szT(d), szT(d) - k + szA(d));
    loA{d} = max(1, k - szT(d) + 1);
    hiA{d} = min(k, szA(d));
end
% the correlation (in the class of the inputs) and the window sums of A0,
% A0 .* A0 and T0 (where it overlaps the image) in double
if useMex
    try
        if isempty(grain)
            grain = matlab_fft_grain();
        end
        [corrTA, sA, sA2, sT] = normxcorr3_window_mex('ncc', T0, A0, bbox, maxNumCompThreads, grain, fastFFT);
    catch ME
        mex_fallback('normxcorr3_window_mex', ME);
        useMex = false;
    end
end
if ~useMex
    corrTA = real(ifftn(fftn(A0, szOut) .* fftn(T0(end : -1 : 1, end : -1 : 1, end : -1 : 1), szOut)));
    if ~isempty(bbox)
        corrTA = crop3d(corrTA, bbox);
    end
    X = double(A0);
    sA = box_sums(X, loA, hiA);
    sA2 = box_sums(X .* X, loA, hiA);
    clear X;
    sT = box_sums(double(T0), loT, hiT);
end
clear A0;

% the rest in double, with the sums of the windows' elements less cA (the
% padding -cA)
cnt = cellfun(@(l, h) h - l + 1, loA, hiA, 'unif', 0);
nPad = pSzT - cnt{1}(:) .* cnt{2}(:)' .* reshape(cnt{3}, 1, 1, []);
cAd = double(cA);
sB = sA - nPad * cAd;
sB2 = sA2 + nPad * cAd ^ 2;
clear sA sA2 nPad;
sumT0 = sum(T0(:), 'double');
denomT = std(double(T0(:)));
clear T0;
num = (double(corrTA) - cAd * (sumT0 - sT) - sB * (sumT0 / pSzT)) / (pSzT - 1);
clear corrTA sT;
varA = max(sB2 - sB .^ 2 / pSzT, 0) / (pSzT - 1);
clear sB sB2;

% the NCC, 0 for windows that are flat to the precision of the image's class
% (variance at most eps times the image's) and for a flat template
valid = varA > eps(class(A)) * var(A(:), 1) & denomT > 0;
C = zeros(size(num), class(A));
C(valid) = num(valid) ./ (denomT * sqrt(varA(valid)));

if ~isempty(bbox)
    return;
end

switch( lower(shape) )
	case 'full'
	case 'same'
		szTp = fix((szT-1)/2);
		C = C( szTp(1)+1:szTp(1)+szA(1), szTp(2)+1:szTp(2)+szA(2), szTp(3)+1:szTp(3)+szA(3) );
	case 'valid'
		C = C(szT(1):end-szT(1)+1,szT(2):end-szT(2)+1,szT(3):end-szT(3)+1);
	otherwise
		error('unknown SHAPE %s, assuming FULL by default', shape);
end

end


function [S] = box_sums(X, lo, hi)
% sums of X over the boxes lo{d}(k) .. hi{d}(k) along each dimension d (a sum
% over every combination of the k of the dimensions), by cumulative sums
S = X;
for d = 1 : 3
    P = cumsum(S, d);
    idx = repmat({':'}, 1, 3);
    idx{d} = hi{d};
    W = P(idx{:});
    m = lo{d} > 1;
    if any(m)
        idx{d} = lo{d}(m) - 1;
        L = P(idx{:});
        idx{d} = find(m);
        W(idx{:}) = W(idx{:}) - L;
    end
    S = W;
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
