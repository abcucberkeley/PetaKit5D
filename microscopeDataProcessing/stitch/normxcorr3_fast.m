function C = normxcorr3_fast(T, A, shape)
% C = normxcorr3(TEMPLATE, IMAGE, SHAPE)
%
%       TEMPLATE - type double, ndims==3, size <= size of image
%       IMAGE    - type double, ndims==3
%       SHAPE    - one of: 'valid', 'same', 'full'. same as conv2 shape parameter
%                  'full' by default
%
%       C        - values in [-1,1]. size depends on SHAPE
%
% adapted from daniel eaton's version with support for 1d and 2d, and
% optimize performance with mex based image integral
%
% Author: Xiongtao Ruan (04/27/2024)


if nargin<3
	shape = 'full';
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

% compute the numerator of the NCC
corrTA = real(ifftn(fftn(A,szOut).*fftn(T(end:-1:1, end:-1:1, end:-1:1), szOut)));

sumT = sum(T(:));
denomT = std(T(:));
% clear T;

% make the running-sum/integral-images of A and A^2, which are
% used to speed up the computation of the NCC denominator
intImgA = integral_image_3d(A,szT);

num = (corrTA - intImgA*sumT/pSzT ) / (pSzT-1);
clear corrTA; 

% compute the denominator of the NCC
intImgA2 = integral_image_3d(A.*A,szT);

denom = denomT * sqrt(max(intImgA2 - (intImgA.^2)/pSzT, 0) / (pSzT-1) );
clear intImgA intImgA2;

% compute the NCC
C = num ./ (denom + eps) .* (denom ~= 0);

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


