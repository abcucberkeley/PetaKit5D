function [max_off, max_corr, C] = normxcorr3_max_shift(T, A, maxShifts)
% compute max xcorr based on fast normalized xcorr
%
% Author: Xiongtao Ruan (10/15/2020)
% xruan (11/02/2020): add constraints for the range of max shifts
% xruan (12/07/2020): update constraints for the range of max shift, set as
% lower bound and upper bound
% xruan (02/04/2022): first crop C based on the maxshift constraint, and
% directly find the max xcorr. For big image, it may stuck in the while loop. 

if size(maxShifts, 1) == 1
    maxShifts = [-maxShifts; maxShifts];
end

% tic
% C = normxcorr3(T, A);
% C = normxcorr3_updated(T, A);
C = normxcorr3_fast(T, A);
% toc

sz_t = size(T, 1 : 3);

s = max(1, ceil(maxShifts(1, :) + sz_t));
t = min(floor(maxShifts(2, :) + sz_t), size(C));

C = crop3d(C, [s, t]);

[max_corr, ind] = max(C(:));
[y, x, z] = ind2sub(size(C), ind);
max_inds = [y, x, z];

max_off = max_inds + (s - 1) - sz_t;


end
