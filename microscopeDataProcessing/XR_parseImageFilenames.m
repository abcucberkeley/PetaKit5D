function [fnames, fdinds, gfnames, partialvols, dataSizes, flipZstack_mat, latest_modify_times, FTP_inds, maskFullpaths] = ...
    XR_parseImageFilenames(dataPaths, ChannelPatterns, parseSettingFile, flipZstack, Decon, deconPaths, Streaming, minModifyTime, zarrFile)
% move the filename parsing code in microscope pipeline as an independent
% function to simpolify the microscope pipeline.
% Support both non-streaming and streaming modes
%
%
% Author: Xiongtao Ruan (07/01/2021)
% also include folder name for channel patterns
% xruan (08/25/2021): add support for zarr file
% xruan (04/05/2022): change to not include last number of channel patterns instead of last one
% xruan (04/27/2022): in streaming mode, multiply minModifyTime to the number of partial files 
% xruan (05/06/2022): in streaming mode, wait at lease minModifyTime

if nargin < 8
    zarrFile = false;
end
if zarrFile
    ext = '.zarr';
else
    ext = '.tif';
end

nd = numel(dataPaths);
% cast dataPaths to column cell array (in case of row arrays)
dataPaths = dataPaths(:);

% check existing files and parse channels
fnames_cell = cell(nd, 1);
gfnames_cell = cell(nd, 1); % for grouped partial volume files
partialvol_cell = cell(nd, 1);
datesize_cell = cell(nd, 1);
latest_modify_times = zeros(nd, 1);
for d = 1 : nd
    dataPath = dataPaths{d};
    % dir_info = dir([dataPath, '*.tif']);
    % fnames_d = {dir_info.name}';
    [containPartialVolume, groupedFnames_d, groupedDatenum, groupedDatasize] = groupPartialVolumeFiles(dataPath, ...
        'ext', ext, 'ChannelPatterns', ChannelPatterns);
    if any(containPartialVolume)
        fnames_d = cellfun(@(x) x{1}, groupedFnames_d, 'unif', 0);
        datenum_d = cellfun(@(x) max(x), groupedDatenum);
        % datesize_d = cellfun(@(x) max(x), groupedDatasize);
        datesize_d = groupedDatasize;
    else
        fnames_d = groupedFnames_d;
        datenum_d = groupedDatenum;
        datesize_d = groupedDatasize;
    end
    
    if isempty(fnames_d)
        continue;
    end    
    
    if Streaming
        last_modify_time = (datenum(clock) - datenum_d) * 24 * 60;
        % exclude last two frames
        latest_modify_time = max(mink(last_modify_time, numel(ChannelPatterns)));
        latest_modify_times(d) = latest_modify_time;
        
        % medium of number of partial files in group files
        ngf = median(cellfun(@numel, groupedFnames_d));

        % not include the lastest file if it is very recent
        if latest_modify_time < minModifyTime * ngf
            latest_modify_time = max(latest_modify_time, minModifyTime);
            fnames_d(last_modify_time <= latest_modify_time) = [];
            if any(containPartialVolume)
                groupedFnames_d(last_modify_time <= latest_modify_time) = [];
            end
        end
    end
    fnames_cell{d} = fnames_d;
    gfnames_cell{d} = groupedFnames_d;
    if any(containPartialVolume)
        partialvol_cell{d} = cellfun(@(x) numel(x) > 1, groupedFnames_d);
        datesize_cell{d} = cellfun(@(x) sum(x), datesize_d);
    else
        partialvol_cell{d} = false(numel(fnames_d), 1);
        datesize_cell{d} = datesize_d;
    end
end

fdinds = arrayfun(@(x) ones(numel(fnames_cell{x}), 1) * x, 1 : nd, 'unif', 0);
fnames = cat(1, fnames_cell{:});
fdinds = cat(1, fdinds{:});
gfnames = cat(1, gfnames_cell{:});
partialvols = cat(1, partialvol_cell{:});
dataSizes = cat(1, datesize_cell{:});

% filter filenames by channel patterns 
% 07/13/2021 also include folder names for channel pattern filtering
if isempty(fnames)
    warning('There is no image files in the dataPaths, please check if dataPaths are correct!');
    fnames = {};
    gfnames = {};
    flipZstack_mat = logical([]);
    FTP_inds = [];
    maskFullpaths = {};
    return;
end

fullnames = cellfun(@(x, y) [x, y], dataPaths(fdinds), fnames, 'unif', 0);
include_flag = false(numel(fnames), 1);
for c = 1 : numel(ChannelPatterns)
    include_flag = include_flag | contains(fullnames, ChannelPatterns{c}) | contains(fullnames, regexpPattern(ChannelPatterns{c}));
end
fnames = fnames(include_flag);
fdinds = fdinds(include_flag);
gfnames = gfnames(include_flag);
partialvols = partialvols(include_flag);
dataSizes = dataSizes(include_flag);

nF = numel(fnames);

% parse setting file 
flipZstack_mat = repmat(flipZstack, nF, 1);
if parseSettingFile
    frameFullpaths = arrayfun(@(x) [dataPaths{fdinds(x)}, fnames{x}], 1 : nF, 'unif', 0);
    settingInfo = XR_parseSettingFiles_wrapper(frameFullpaths);
    flipZstack_mat = [settingInfo.StageInterval] < 0;
end

% for ErodeByFTP, set up the frame numbers as first time point for each
% data
FTP_inds = zeros(nd, 1);
maskFullpaths = cell(nd, 1);
for d = 1 : nd
    c = 1;
    FTPfname = '';
    if isempty(fnames_cell{d})
        continue;
    end
    while isempty(FTPfname)
        fullnames_d = cellfun(@(x) [dataPaths{d}, x], fnames_cell{d}, 'unif', 0);
        all_inds = contains(fullnames_d, ChannelPatterns{c}) | contains(fullnames_d, regexpPattern(ChannelPatterns{c}));
        if ~isempty(all_inds) && ~isempty(find(all_inds, 1, 'first'))
            FTPfname = fnames_cell{d}{find(all_inds, 1, 'first')};
        end
        c = c + 1;
        if c > numel(ChannelPatterns)
            break;
        end
    end
    if ~isempty(FTPfname)
        ind_d = find(strcmp(fnames, FTPfname) & fdinds == d);
        FTP_inds(d) = ind_d;
        [~, FTPfsname] = fileparts(FTPfname);        
        if Decon
            maskFullpaths{d} = sprintf('%s/Masks/%s_eroded.tif', deconPaths{d}, FTPfsname);
        end
    end
end    


end

