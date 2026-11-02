function [max_off, max_corr, C] = normxcorr3_max_shift(T, A, maxShifts, fastFFT)
% compute max xcorr based on fast normalized xcorr
%
% fastFFT: true or false (default: false). Smaller FFTs for the xcorr within
% the max shifts (see normxcorr3_fast), which changes it in the last bits.
%
% Author: Xiongtao Ruan (10/15/2020)
% xruan (11/02/2020): add constraints for the range of max shifts
% xruan (12/07/2020): update constraints for the range of max shift, set as
% lower bound and upper bound
% xruan (02/04/2022): first crop C based on the maxshift constraint, and
% directly find the max xcorr. For big image, it may stuck in the while loop. 

if nargin < 4
    fastFFT = false;
end

if size(maxShifts, 1) == 1
    maxShifts = [-maxShifts; maxShifts];
end

sz_t = size(T, 1 : 3);
% size of the full xcorr
sz_c = sz_t + size(A, 1 : 3) - 1;

s = max(1, ceil(maxShifts(1, :) + sz_t));
t = min(floor(maxShifts(2, :) + sz_t), sz_c);

% tic
% C = normxcorr3(T, A);
% C = normxcorr3_updated(T, A);
% only the part within the max shifts
C = normxcorr3_fast(T, A, 'full', [s, t], fastFFT);
% toc

[max_corr, ind] = max(C(:));
[y, x, z] = ind2sub(size(C), ind);
max_inds = [y, x, z];

max_off = max_inds + (s - 1) - sz_t;


end
