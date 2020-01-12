function [data] = XR_correctXZoffsetData3D(data, varargin)
% adjust offsets for slave channels, based on the primary channel, the idea
% is to find the offset by the calcualte the cross correlation between the 
% frames of the same. The final offset is the average of those for each
% frame. 
% 
% Inputs (optional):
%      data : structure returned by loadConditionData()
%
% Options ('specifier', value):
%         'Overwrite' : true|{false}, enables overwriting of previously de-skewed/rotated data
%                'DS' : use deskewed data as input images. Default: true
%   'sCMOSCameraFlip' : sCMOS camera flip. Default: false 
% 
% 
% Author: Xiongtao Ruan 11/24/2019
% xruan 11/25/2019 sorted channels so that different running will not
% conflict each other if the channels are permutated for different data. 
% xruan 11/27/2019 first write all corrected files to a new folder, than
% delete the old folder and move the new folder to DS. 

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('DS', true, @islogical); %
ip.addParameter('sCMOSCameraFlip', false, @islogical); %
ip.parse(data, varargin{:});

Overwrite = ip.Results.Overwrite;

nd = numel(data);
for i = 1 : nd
    data_i = data(i);
    channels = data_i.channels;
    % make sure the channels are sorted
    if ~issorted(channels)
        [~, sinds] = sort(channels);
        data_i.channels = data_i.channels(sinds);
        data_i.framePaths = data_i.framePaths(sinds);
        data_i.framePathsDS = data_i.framePathsDS(sinds);
        data_i.framePathsDSR = data_i.framePathsDSR(sinds);
    end
    
    nCh = numel(data_i.channels);
    nF = numel(data_i.framePaths{1});
    
    fprintf('Correct xz-offset for ''%s''\n', getShortPath(data_i));

    % estimate offset in z-axis. 
    if nF < 3
        Fchosen = 1 : nF;
    else
        Fchosen = [1, round((nF - 1) / 2), nF];
    end
    xzoff_mat = zeros(numel(Fchosen), nCh - 1, 2);
    xzoff_done_mat = false(nCh, 1);
    % only check Ch1 in case of conflict. 
    for j = 1
        DSPath = fileparts(data_i.framePathsDS{j}{1});
        
        xzoff_flag_filename = sprintf('%s/xz_offset.mat', DSPath);
        if exist(xzoff_flag_filename, 'file')            
            xzoff_done_mat(j) = true;
            break;
        end
    end
    
    if any(xzoff_done_mat) && ~Overwrite
        continue;
    end
    
    % clear old tmp file. 
    chunk_lock_clean(DSPath, 5 * numel(data_i.framePathsDS{1}) * nCh);
    tmp_filename = sprintf('%s/xz_offset.tmp', DSPath);
    if exist(tmp_filename, 'file')
        continue;
    else
        fclose(fopen(tmp_filename, 'w'));
    end    

    for k = 1 : numel(Fchosen)
        FInd = Fchosen(k);
        if ip.Results.DS
            prim_Path = data_i.framePathsDS{1}{FInd};
        else
            prim_Path = data_i.framePaths{1}{FInd};
        end
        
        frame_p = double(readtiff(prim_Path));
        frame_p(frame_p == 0) = nan;
        if ip.Results.sCMOSCameraFlip
            frame_p = flip(frame_p, 3);
        end            
        
        for j = 1 : nCh - 1
            ch_j = j + 1;
            if ip.Results.DS
                Path_sj = data_i.framePathsDS{ch_j}{FInd};
            else
                Path_sj = data_i.framePaths{ch_j}{FInd};
            end
            
            frame_sj = double(readtiff(Path_sj));
            frame_sj(frame_sj == 0) = nan;
            if ip.Results.sCMOSCameraFlip
                frame_sj = flip(frame_sj, 3);
            end            
    
            tic
            [xoff, zoff] = XZoffsetEstimation3D(frame_p, frame_sj);
            toc
            xzoff_mat(k, j, :) = [xoff, zoff];
        end
    end
    
    xoff_move_mat = [0, fix(median(xzoff_mat(:, :, 1)))];
    zoff_move_mat = [0, fix(median(xzoff_mat(:, :, 2)))];
    % [bottom, top] pad number
    xpad_mat = abs(xoff_move_mat' - [min(xoff_move_mat), max(xoff_move_mat)]);
    zpad_mat = abs(zoff_move_mat' - [min(zoff_move_mat), max(zoff_move_mat)]);
    channels = data_i.channels;
    
    pad_mat = cat(3, xpad_mat, zpad_mat);
    
    % perform the correction
    % reverse correction, last for Ch1
    for j = nCh : -1 : 1     
        DSPath = fileparts(data_i.framePathsDS{j}{1});
        xzoff_flag_filename = sprintf('%s/xz_offset.mat', DSPath);        
        if exist(xzoff_flag_filename, 'file') && ~ip.Results.Overwrite
            continue;
        end
        
        cur_pad_mat = squeeze(pad_mat(j, :, :))';
        if all(cur_pad_mat == 0)
            save('-v7.3', xzoff_flag_filename, 'xzoff_mat', 'pad_mat', 'channels');
            continue;
        end
        
        nPad = sum(cur_pad_mat, 2);
        DSCorrectedPath = [DSPath, '_tmp_xz_correct/']; 
        mkdir(DSCorrectedPath);
        for k = 1 : nF
            Path = data_i.framePathsDS{j}{k};
            im = readtiff(Path);
            im2 = zeros(size(im)+[0 nPad']);
            im2(:, cur_pad_mat(1, 1) + 1 : end - cur_pad_mat(1, 2), cur_pad_mat(2, 1) + 1 : end - cur_pad_mat(2, 2)) = im;
            
            [~, fname] = fileparts(Path);
            CorrectedFilePath = sprintf('%s/%s.tif', DSCorrectedPath, fname);
            
            writetiff(single(im2), CorrectedFilePath);           
        end
        
        % remove old folder, and move the results to DS folder
        rmdir(DSPath, 's');
        movefile(DSCorrectedPath, DSPath);
        
        save('-v7.3', xzoff_flag_filename, 'xzoff_mat', 'pad_mat', 'channels');
    end
    if exist(tmp_filename, 'file')
        delete(tmp_filename);
    end
    
    fprintf('Done!');
end


end


