function [max_off, max_corr, C] = normxcorr2_max_shift(T, A, maxShifts)
% compute max xcorr based on fast normalized xcorr
% 
% from normcorr3_max_shift.m
%
% Author: Xiongtao Ruan (05/09/2022)
% xruan (11/02/2020): add constraints for the range of max shifts
% xruan (12/07/2020): update constraints for the range of max shift, set as
% lower bound and upper bound
% xruan (02/04/2022): first crop C based on the maxshift constraint, and
% directly find the max xcorr. For big image, it may stuck in the while loop. 

if size(maxShifts, 1) == 1
    maxShifts = [-maxShifts; maxShifts];
end

tic
C = normxcorr2(T, A);
toc

sz_t = size(T);

s = max(1, ceil(maxShifts(1, 1 : 2) + sz_t));
t = min(floor(maxShifts(2, 1 : 2) + sz_t), size(C));

C = C(s(1) : t(1), s(2) : t(2));

[max_corr, ind] = max(C(:));
[y, x] = ind2sub(size(C), ind);
max_inds = [y, x];

max_off = [max_inds + (s - 1) - sz_t, 0];

% k = 100;
% mind = [];
% while isempty(mind)
%     [~, inds] = maxk(C(:), k);
%     [y, x, z] = ind2sub(size(C), inds);
%     offsets = [y, x, z] - size(T);
%     % mind = find(all(abs(offsets) <= maxShifts, 2), 1, 'first');
%     mind = find(all(offsets >= maxShifts(1, :) & offsets <= maxShifts(2, :), 2), 1, 'first');
%     k = k * 10;
% end
% 
% max_off = offsets(mind, :);
% max_corr = C(y(mind), x(mind), z(mind));

end
