function [IntA] = integral_image_3d(A, szT)
% integral of 3d image
%
% Author: Xiongtao Ruan (04/27/2024)


szA = size(A, 1 : 3);

if any(szA < szT)
    error('The size of A must be not smaller than size T!')
end

try
    IntA = integral_image_3d_mex(A, szT);
catch ME
    disp(ME)
    IntA = integralImage(A, szT);
end

end


function c = integralImage(A,szT)
% this is adapted from Matlab's normxcorr2

szA = size(A, 1 : 3);

B = zeros( szA+2*szT-1, class(A));
% B( szT(1)+1:szT(1)+szA(1), szT(2)+1:szT(2)+szA(2), szT(3)+1:szT(3)+szA(3) ) = A;
indexing3d_mex(B, [szT+1, szT+szA], A);

s = cumsum(B,1);
clear B;
% c = s(1+szT(1):end,:,:)-s(1:end-szT(1),:,:);
c = crop3d_mex(s, [1 + szT(1), 1, 1, size(s, 1 : 3)]) - crop3d_mex(s, [1, 1, 1, size(s, 1) - szT(1), size(s, 2 : 3)]);
s = cumsum(c,2);
% c = s(:,1+szT(2):end,:)-s(:,1:end-szT(2),:);
c = crop3d_mex(s, [1, 1 + szT(2), 1, size(s, 1 : 3)]) - crop3d_mex(s, [1, 1, 1, size(s, 1), size(s, 2) - szT(2), size(s, 3)]);
s = cumsum(c,3);
% integralImageA = s(:,:,1+szT(3):end)-s(:,:,1:end-szT(3));
c = crop3d_mex(s, [1, 1, 1 + szT(3), size(s, 1 : 3)]) - crop3d_mex(s, [1, 1, 1, size(s, 1), size(s, 2), size(s, 3) - szT(3)]);

end