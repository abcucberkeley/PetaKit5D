function [] = XR_correctZoffsetData3D(data, varargin)
% adjust offsets for slave channels, based on the primary channel, the idea
% is to find the offset by the calcualte the cross correlation between the 
% frames of the same. The final offset is the average of those for each
% frame. 
% 
% Author: Xiongtao Ruan 11/18/2019


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('DS', false, @islogical); %
ip.parse(data, varargin{:});

nd = numel(data);
for i = 1 : nd
    nCh = numel(data(i).channels);
    nF = numel(data(i).framePaths{1});
    
    fprintf('Correct z-offset for ''%s''\n', getShortPath(data(i)));

    % estimate offset in z-axis. 
    if nF < 3
        Fchosen = 1 : nF;
    else
        Fchosen = [1, round((nF - 1) / 2), nF];
    end
    zoff_mat = zeros(numel(Fchosen), nCh - 1);
    zoff_done_mat = false(nCh, 1);
    for j = 1 : nCh 
        DSPath = fileparts(data(i).framePathsDS{j}{1});
        
        zoff_flag_filename = sprintf('%s/z_offset.mat', DSPath);
        if exist(zoff_flag_filename, 'file')
            zoff_done_mat(j) = true;
        end
    end
    
    if all(zoff_done_mat) && ~ip.Results.Overwrite
        continue;
    end
    
    for k = 1 : numel(Fchosen)
        FInd = Fchosen(k);
        prim_Path = data(i).framePathsDS{1}{FInd};
        frame_p = double(readtiff(prim_Path));
        frame_p(frame_p == 0) = nan;
        
        for j = 1 : nCh - 1
            ch_j = j + 1;
            Path_sj = data(i).framePathsDS{ch_j}{FInd};
            frame_sj = double(readtiff(Path_sj));
            frame_sj(frame_sj == 0) = nan;
            tic
            zoff = offsetEstimation3D(frame_p, frame_sj);
            toc
            zoff_mat(k, j) = zoff;
        end
    end
    
    zoff_move_mat = [0, fix(median(zoff_mat))];
    % [bottom, top] pad number
    pad_mat = abs(zoff_move_mat' - [min(zoff_move_mat), max(zoff_move_mat)]);
    nPad = pad_mat(1, 1) + pad_mat(1, 2);
    channels = data(1).channels;
    
    % perform the correction
    for j = 1 : nCh     
        DSPath = fileparts(data(i).framePathsDS{j}{1});
        zoff_flag_filename = sprintf('%s/z_offset.mat', DSPath);        
        if exist(zoff_flag_filename, 'file') && ~ip.Results.Overwrite
            continue;
        end
        
        cur_pad_mat = pad_mat(j, :);
        if all(cur_pad_mat == 0)
            save('-v7.3', zoff_flag_filename, 'zoff_mat', 'pad_mat', 'channels');
            continue;
        end
        
        for k = 1 : nF
            Path = data(i).framePathsDS{j}{k};
            im = readtiff(Path);
            im2 = zeros(size(im)+[0 0 nPad]);
            im2(:, :, cur_pad_mat(1) + 1 : end - cur_pad_mat(2)) = im;

            writetiff(single(im2), Path);           
        end
        
        save('-v7.3', zoff_flag_filename, 'zoff_mat', 'pad_mat', 'channels');
    end
    fprintf('Done!');    
end


end



function [zoff] = offsetEstimation3D(frame_1, frame_2, zoff_mat)
% function to estimate offset between frame_1 and frame_2
% require to be a local maximum

if nargin < 3
    zoff_mat = -6 : 6;
end

[ny, nx, nz] = size(frame_1);
num_off = numel(zoff_mat);

xcorr_mat = zeros(num_off, 1);

for i = 1 : num_off
    zoff = zoff_mat(i);
    
    if zoff <= 0
        frame_1_z = frame_1(:, :, 1 : nz + zoff);
        frame_2_z = frame_2(:, :, 1 - zoff : nz);
    else
        frame_1_z = frame_1(:, :, 1 + zoff : nz);
        frame_2_z = frame_2(:, :, 1 : nz - zoff);
    end
    nonnan_inds = ~isnan(frame_1_z + frame_2_z);
    
    inds = nonnan_inds;
    
    xcorr_z = corr(frame_1_z(inds) , frame_2_z(inds));
    xcorr_mat(i) = xcorr_z;
end

[pks, locs] = findpeaks(xcorr_mat, 'minPeakProminence', 0.01);
locs(locs == 1 | locs == num_off) = [];
if numel(locs) > 1
    [~, max_ind] = max(pks);
    locs = locs(max_ind);
end

if ~isempty(locs)
    zoff = zoff_mat(locs);
else
    zoff = 0;
end


end

