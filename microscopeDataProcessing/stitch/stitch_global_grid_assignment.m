function [d_shift] = stitch_global_grid_assignment(nF, max_xcorr_mat, absolute_shift_mat, overlap_matrix, max_allow_shift, xcorr_thresh, tileIdx, axisWeight)
% solve weighted constrained linear least square problem for the assignment
% the weight is the function of max shift (currently just y=x). 
% max allowed shift is based on the maxShift parameter and the number of
% overlap between tiles.
% 02/03/2022: add threshold for objective to remove the samples with small xcorr

neq = sum(overlap_matrix(:));

R = zeros(neq, nF);

[n_i, n_j] = find(overlap_matrix);

inds_i = sub2ind(size(R), 1 : neq, n_i');
inds_j = sub2ind(size(R), 1 : neq, n_j');

R(inds_i) = -1; 
R(inds_j) = 1;

max_xcorr_mat_filt = max_xcorr_mat;
filt_inds = max_xcorr_mat_filt(:, 3) < xcorr_thresh;
max_xcorr_mat_filt(filt_inds, :) = [];
% w = max_xcorr_mat_filt(:, 3);
nP = size(max_xcorr_mat_filt, 1);
w = ones(nP, 1) * 0.00;
for i = 1 : nP
    s = max_xcorr_mat_filt(i, 1);
    t = max_xcorr_mat_filt(i, 2);
    
    % check if two tiles are neighboring tiles
    if sum(abs(tileIdx(s, 1 : 3) - tileIdx(t, 1 : 3))) == 1
        aind = find(abs(tileIdx(s, 1 : 3) - tileIdx(t, 1 : 3)));
        switch aind
            case 1
                % x
                w(i) = axisWeight(2);
            case 2
                % y
                w(i) = axisWeight(1);
            case 3
                % z
                w(i) = axisWeight(3);
        end
    end    
end

R_w = zeros(nP, nF);

np_i = max_xcorr_mat_filt(:, 1);
np_j = max_xcorr_mat_filt(:, 2);
inds_i = sub2ind(size(R_w), 1 : nP, np_i');
inds_j = sub2ind(size(R_w), 1 : nP, np_j');
R_w(inds_i) = -1; 
R_w(inds_j) = 1;

W = diag(w);
R_w = W.^0.5 * R_w;
% R_w = R;

% c_inds = (np_i - 1) * nF - np_i .* (np_i + 1) / 2 + np_j;
d_w = absolute_shift_mat(~filt_inds, 3 : 5);
d_w = W.^0.5 * d_w;

d_shift = zeros(nF, 3);
for i = 1 : 3
    C = R_w;
    d = d_w(:, i);
    A = [R; -R];

    l = max_allow_shift(:, i);
    u = max_allow_shift(:, 3 + i);
    b = [u; -l];

    [x,resnorm,residual,exitflag,output,lambda] = lsqlin(C,d,A,b);
    d_shift(:, i) = x;
end

% round to integers and normalize for the first tile.
d_shift = round(d_shift);
d_shift = d_shift - d_shift(1, :);

end
