function [xoff, zoff] = XZoffsetEstimation3D(frame_1, frame_2, zoff_mat)
% function to estimate offset between frame_1 and frame_2
% require to be a local maximum
% 
% Author: Xiongtao Ruan 11/24/2019
% 

if nargin < 3
    xoff_mat = -7 : 7;
    zoff_mat = -5 : 5;
end

[ny, nx, nz] = size(frame_1);
num_xoff = numel(xoff_mat);
num_zoff = numel(zoff_mat);

xcorr_mat = zeros(num_xoff, num_zoff);

for i = 1 : num_xoff
    xoff = xoff_mat(i);
    if xoff <= 0
        xa_1 = 1 : nx + xoff;
        xa_2 = 1 - xoff : nx;
    else
        xa_1 = 1 + xoff : nx;
        xa_2 = 1 : nx - xoff;
    end
    
    for j = 1 : num_zoff
        zoff = zoff_mat(j);        
        if zoff <= 0
            za_1 = 1 : nz + zoff;
            za_2 = 1 - zoff : nz;
        else
            za_1 = 1 + zoff : nz;
            za_2 = 1 : nz - zoff;
        end
        frame_1_xz = frame_1(:, xa_1, za_1);
        frame_2_xz = frame_2(:, xa_2, za_2);

        nonnan_inds = ~isnan(frame_1_xz + frame_2_xz);

        inds = nonnan_inds;

        xcorr_xz = corr(frame_1_xz(inds) , frame_2_xz(inds));
        xcorr_mat(i, j) = xcorr_xz;
    end
end

% [pks, locs] = findpeaks(xcorr_mat, 'minPeakProminence', 0.01);
BW = imregionalmax(xcorr_mat);
% remove boundary values
BW(1 : end, [1, end]) = false;
BW([1, end], 1 : end) = false;

[lmx, lmz] = find(BW);

% if there are multiple ones, use the one mostly close the center
% (conservative correction)
if numel(lmz) > 1
    [~, min_ind] = min((lmz - num_zoff / 2) .^ 2 + (lmx - num_xoff / 2) .^ 2);
    lmx = lmx(min_ind);      
    lmz = lmz(min_ind);
end

if ~isempty(lmz)
    xoff = xoff_mat(lmx);    
    zoff = zoff_mat(lmz);
else
    xoff = 0;
    zoff = 0;
end


end

