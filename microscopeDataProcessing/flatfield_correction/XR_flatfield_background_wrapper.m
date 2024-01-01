function [] = XR_flatfield_background_wrapper(dataPaths, varargin)
% wrapper for flatfield background estimation. 
% 
% 
% Inputs :   
%        imageDirName : full path for the directory of the images to be stitched
%   imageListFileName : full path for the coordinate information csv file
%
% Options (as 'specifier'-value pairs): 
%
%         'axisOrder' : Axis mapping of the coordinate system. Default: 'xyz'.
%         'ffcorrect' : flat field correction. Default: false.
%        'Resolution' : Image resolution in um, 2D (same for xy) or 3D vector.  
%         'resultDir' : name of stitching result directory. The directory is child directory of imageDirname
%       'BlendMethod' : Blend method for overlap regions. Available: none, mean, max, median. Default: mean
%           'padSize' : Pad or crop the stitched image, empty (default) or 
%                       a 1X3 vector of integers (y, x, z). Postive for pad and negative for crop. 
%      'boundboxCrop' : Crop the stitched image by some bounding box, empty (default, no crop) 
%                       or a 3X2 vector for start and end indices of the bounding box (y, x, z).      
%        'zNormalize' : normalize background along z-axis by the median. Default: false.
%      'parseCluster' : Use slurm-based cluster computing. Default: true. 
%         'jobLogDir' : Log directory for the slurm jobs. 
%       'cpusPerTask' : Number of cpus for a job. Default: 12
%         'Save16bit' : true|{false}. Save final results as 16bit or single. 
%              'uuid' : unique string for a job for saving files. 
%       'maxTrialNum' : Max number of times to rerun failure cases. 
%      'unitWaitTime' : For computing without cluster, the wait time per      
%                       file in minutes, in order to check whether the computing is done. 
%
%
% Author: Xiongtao Ruan (10/21/2020)
%
% xruan (09/23/2021): refactor code to support per-tile estimation (for a tile with lots z slices)
% xruan (08/01/2023): add options to save MIPs (MIP z) for ff estimation
% (to avoid compute MIPs in a separate step); also add support for multiple
% datasets


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
% ip.addParameter('Overwrite', true, @islogical);
ip.addParameter('resultDirStr', 'matlab_flat_field_estimation', @ischar);
ip.addParameter('ChannelPatterns', {'CamA_ch0', 'CamB_ch0'}, @iscell); % only compute first time point (for deciding cropping bouding box)
ip.addParameter('estimationMode', 'per_image', @ischar); % 'cam_ch', 'came_ch_iter', 'came_ch_iter_stack', or 'per_image'
ip.addParameter('onlyFirstTP', false, @islogical); % only compute first time point (for deciding cropping bouding box)
ip.addParameter('Save16bit', false, @islogical);
ip.addParameter('saveMIP', false, @islogical); % save MIPs
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 24, @isnumeric);
ip.addParameter('mccMode', false, @ischar);
ip.addParameter('ConfigFile', '', @ischar);
ip.addParameter('uuid', '', @ischar);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
% Overwrite = pr.Overwrite;
ChannelPatterns = pr.ChannelPatterns;
estimationMode = pr.estimationMode;
resultDirStr = pr.resultDirStr;
saveMIP = pr.saveMIP;
jobLogDir = pr.jobLogDir;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;
uuid = pr.uuid;


% make root directory
if ischar(dataPaths)
    dataPaths = {dataPaths};
end

nd = numel(dataPaths);

for d = 1 : nd
    dataPath = dataPaths{d};
    if ~exist(dataPath, 'dir')
        error('%s does not exist, please double check the path to make sure it is correct!', dataPath);
    end

    bg_rt = [dataPath, filesep, resultDirStr];
    if ~exist(bg_rt, 'dir')
        mkdir(bg_rt);
        fileattrib(bg_rt, '+w', 'g');            
    end
    
    % save parameters 
    save('-v7.3', [bg_rt, '/parameters.mat'], 'pr');
    writetable(struct2table(pr, 'AsArray', true), [bg_rt, '/parameters.txt'])
    
    % temporary directory for intermediate results
    cdr_tmp = [bg_rt, filesep, 'tmp'];
    if ~exist(cdr_tmp, 'dir')
        mkdir(cdr_tmp);
        fileattrib(cdr_tmp, '+w', 'g');            
    end

    if saveMIP
        MIPPath = [bg_rt, filesep, 'MIPs'];
        if ~exist(MIPPath, 'dir')
            mkdir(MIPPath);
            fileattrib(MIPPath, '+w', 'g');            
        end
    end
end

% check if a slurm-based computing cluster exists
if parseCluster
    [parseCluster, job_log_fname, job_log_error_fname, slurm_constraint_str, jobLogDir] = checkSlurmCluster(dataPaths{1}, jobLogDir);
end

% uuid for the job
if isempty(uuid)
    uuid = get_uuid();
end

parseSettingFile = false;
flipZstack = false;
Decon = false;
deconPaths = '';
Streaming = false;

[fnames, fdinds, gfnames, partialvols, dataSizes, flipZstack_mat, FTP_inds, maskFullpaths] = ...
    XR_parseImageFilenames(dataPaths, ChannelPatterns, parseSettingFile, flipZstack, Decon, deconPaths, Streaming);

nF = numel(fnames);


%% parse image list information

% dir_info = dir([dataPath, '/*.tif']);
fn = fnames;
dbytes = dataSizes;

specifyCam = true;
if all(~cellfun(@isempty, regexp(fn, '_Cam\w_ch', 'match')))
    expression = '(?<prefix>\w*)Scan_Iter_(?<Iter>\d+)(?<subIter>_?\d+?_?\d+?)_Cam(?<Cam>\w+)_ch(?<ch>\d+)_CAM1_stack(?<stack>\d+)_(?<laser>\d+)nm_(?<abstime>\d+)msec_(?<fpgatime>\d+)msecAbs_(?<x>-?\d+)x_(?<y>-?\d+)y_(?<z>-?\d+)z_(?<t>\d+)t.tif';
elseif all(~cellfun(@isempty, regexp(fn, '_ch[0-9]_', 'match')))
    expression = '(?<prefix>\w*)Scan_Iter_(?<Iter>\d+)(?<subIter>_?\d+?_?\d+?)_ch(?<ch>\d+)_CAM1_stack(?<stack>\d+)_(?<laser>\d+)nm_(?<abstime>\d+)msec_(?<fpgatime>\d+)msecAbs_(?<x>-?\d+)x_(?<y>-?\d+)y_(?<z>-?\d+)z_(?<t>\d+)t.tif';
    specifyCam = false;
end

tmp = regexpi(fn, expression, 'names');

matched_inds = true(numel(tmp), 1);
% t = cell2table(cell(0, 15), 'VariableName', {'Filename', 'prefix', 'Iter', 'subIter', 'fullIter', 'camera', 'ch', 'stack', 'laser', 'abstime', 'fpgatime', 'x', 'y', 'z', 't'});
t = table();

for f = 1:numel(tmp)
    if isempty(tmp{f})
        matched_inds(f) = false;
        continue;
    end
    t.Filename{f} = fn{f};
    t.prefix{f} = tmp{f}.prefix;
    t.Iter(f) = str2double(tmp{f}.Iter);
    t.subIter{f} = tmp{f}.subIter;
    t.fullIter{f} = [tmp{f}.Iter, tmp{f}.subIter];
    if specifyCam
        t.camera(f) = (tmp{f}.Cam);
    else
        % use N to represent cam if it is not contained in the filename
        t.camera(f) = 'N';
    end
    t.ch(f) = str2double(tmp{f}.ch);
    t.stack(f) = str2double(tmp{f}.stack);
    t.laser(f) = str2double(tmp{f}.laser);
    t.abstime(f) = str2double(tmp{f}.abstime);
    t.fpgatime(f) = str2double(tmp{f}.fpgatime);
    t.dbytes(f) = dbytes(f);
    t.x(f) = str2double(tmp{f}.x);
    t.y(f) = str2double(tmp{f}.y);
    t.z(f) = str2double(tmp{f}.z);
    t.t(f) = str2double(tmp{f}.t);
    t.did(f) = fdinds(f);
end

t = t(matched_inds, :);

prefix = unique(t.prefix);
if ~isempty(prefix)
    prefix = prefix{1};
else
    prefix = '';
end
Iter = unique(t.Iter);
Ch = unique(t.ch);
Cam = unique(t.camera);
stackn = unique(t.stack);

if pr.onlyFirstTP
    Iter = Iter(1);
end

%% computing

switch estimationMode
    case 'cam_ch'
        inputFullpaths = cell(nd, numel(Cam), numel(Ch));
        input_str_filenames = cell(nd, numel(Cam), numel(Ch));
        outputFullpaths = cell(nd, numel(Cam), numel(Ch));
    case 'came_ch_iter'
        inputFullpaths = cell(nd, numel(Cam), numel(Ch), numel(Iter));
        input_str_filenames = cell(nd, numel(Cam), numel(Ch), numel(Iter));
        outputFullpaths = cell(nd, numel(Cam), numel(Ch), numel(Iter));
    case 'cam_ch_iter_stack'
        inputFullpaths = cell(nd, numel(Cam), numel(Ch), numel(Iter), nume(stackn));
        input_str_filenames = cell(nd, numel(Cam), numel(Ch), numel(Iter), nume(stackn));
        outputFullpaths = cell(nd, numel(Cam), numel(Ch), numel(Iter), nume(stackn));        
end

if ~strcmp(estimationMode, 'per_image')
    for d = 1 : nd
        dataPath = dataPaths{d};
        for c = 1:numel(Ch)   
            for ncam = 1:numel(Cam)
                cur_inds = t.did == d & t.ch == Ch(c) & t.camera == Cam(ncam);
                cur_t = t(cur_inds, :);
                if ~any(cur_inds)
                    continue;
                end
    
                switch estimationMode
                    case 'cam_ch'
                        inputFullpaths{c, ncam} = sprintf('%s/%s', dataPath, cur_t.Filename{1});
                        input_str_filenames{c, ncam} = cellfun(@(x) sprintf('%s/%s', dataPath, x), cat(1, gfnames{cur_inds}), 'unif', 0);
                        outputFullpaths{c, ncam} = sprintf('%s/%s.mat', bg_rt, cur_t.Filename{1}(1 : end - 4));
                    case 'came_ch_iter'
                        for n = 1:numel(Iter)                    
                            cur_inds = t.ch == Ch(c) & t.camera == Cam(ncam) & t.Iter == Iter(n);
                            cur_t = t(cur_inds, :);
                            if ~any(cur_inds)
                                continue;
                            end
    
                            inputFullpaths{c, ncam, n} = sprintf('%s/%s', dataPath, cur_t.Filename{1});
                            input_str_filenames{c, ncam, n} = cellfun(@(x) sprintf('%s/%s', dataPath, x), cat(1, gfnames{cur_inds}), 'unif', 0);
                            outputFullpaths{c, ncam, n} = sprintf('%s/%s/%s.mat', dataPath, resultDirStr, cur_t.Filename{1}(1 : end - 4));
                        end
                    case 'cam_ch_iter_stack'
                        % to be implemented
                end
            end
        end
    end

    empty_inds = cellfun(@isempty, inputFullpaths);
    inputFullpaths(empty_inds) = [];
    input_str_filenames(empty_inds) = [];
    outputFullpaths(empty_inds) = [];
elseif strcmp(estimationMode, 'per_image')
    inputFullpaths = arrayfun(@(x) sprintf('%s/%s', dataPaths{fdinds(x)}, fnames{x}), 1 : nF, 'unif', 0);
    input_str_filenames = arrayfun(@(x) cellfun(@(y) sprintf('%s/%s', dataPaths{fdinds(x)}, y), gfnames{x}, 'unif', 0), 1 : nF, 'unif', 0); 
    outputFullpaths = arrayfun(@(x) sprintf('%s/%s/%s.mat', dataPaths{fdinds(x)}, resultDirStr, fnames{x}(1 : end - 4)), 1 : nF, 'unif', 0);
end

funcStrs = arrayfun(@(x) sprintf('XR_flatfield_background_estimation(%s,''%s'',''saveMIP'',%s)', ...
    sprintf('{''%s''}', strjoin(input_str_filenames{x}, ''',''')), outputFullpaths{x}, string(saveMIP)), ...
    1 : numel(inputFullpaths), 'unif', 0);

imSizes = cellfun(@(y) getImageSize(sprintf('%s/%s', dataPaths{fdinds(1)}, y)), gfnames{1}, 'unif', 0);
imSizes = cat(1, imSizes{:});
imSize = [imSizes(1, 1 : 2), sum(imSizes(:, 3))];
dtype = getImageDataType([dataPaths{fdinds(1)}, '/', gfnames{1}{1}]);
byteNum = dataTypeToByteNumber(dtype);
memAllocate = prod(imSize) * byteNum / 2^30 * 5;

is_done_flag= generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, ...
    funcStrs, cpusPerTask=cpusPerTask, memAllocate=memAllocate, parseCluster=parseCluster, ...
    masterCompute=masterCompute, mccMode=mccMode, ConfigFile=ConfigFile);

if ~all(is_done_flag)
    generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, funcStrs, ...
        cpusPerTask=cpusPerTask, memAllocate=memAllocate*2, parseCluster=parseCluster, ...
        masterCompute=masterCompute, mccMode=mccMode, ConfigFile=ConfigFile);
end

%% collect results and take average
if ~all(is_done_flag, 'all')
    error('Some output files are missing!');
end

% collect results for specific channel and specific dataset
for d = 1 : nd
    dataPath = dataPaths{d};
    bg_rt = sprintf('%s/%s/', dataPath, resultDirStr);

    for c = 1:numel(Ch)   
        for ncam = 1:numel(Cam)
            bg_fname = sprintf('%s/%sScan_Cam%s_ch%d_CAM1_background.tif', ...
                bg_rt, prefix, Cam(ncam), Ch(c));
            
            if exist(bg_fname, 'file')
                is_done_flag(:, ncam, :, c) = true;
                continue;
            end
            
            cur_inds = find(contains(outputFullpaths, sprintf('Cam%s_ch%d', Cam(ncam), Ch(c))));
            if isempty(cur_inds)
                continue;
            end
            
            z_cell = cell(numel(cur_inds), 1);
            sz_mat = zeros(numel(cur_inds), 1);
            for i = 1 : numel(cur_inds)
                a = load(outputFullpaths{cur_inds(i)}, 'mu', 'numSlice');
                z_cell{i} = a.mu;
                sz_mat(i) = a.numSlice;
            end
            z_cell(sz_mat == 0) = [];
            sz_mat(sz_mat == 0) = [];
            
            z_nc = cat(3, z_cell{:});
            w_nc = sz_mat(:) ./ sum(sz_mat, 'all');
            w_nc = permute(w_nc(:), [2, 3, 1]);
            z_im = sum(z_nc .* w_nc, 3);

            bg_tmpname = sprintf('%s/%sScan_Cam%s_ch%d_CAM1_background_%s.tif', ...
                bg_rt, prefix, Cam(ncam), Ch(c), uuid);

            writetiff(z_im, bg_tmpname);
            movefile(bg_tmpname, bg_fname);
        end
    end
end

% take the 95% intensity of all background slices (test)
if false
    percent = 95;
    for c = 1:numel(Ch)   
        for ncam = 1:numel(Cam)
            bg_fname = sprintf('%s/%sScan_Cam%s_ch%d_CAM1_background_percentile_%d.tif', ...
                bg_rt, prefix, Cam(ncam), Ch(c), percent);
            
            if exist(bg_fname, 'file')
                is_done_flag(:, ncam, :, c) = true;
                continue;
            end
            
            cur_inds = find(contains(outputFullpaths, sprintf('Cam%s_ch%d', Cam(ncam), Ch(c))));
            if isempty(cur_inds)
                continue;
            end
            
            z_cell = cell(numel(cur_inds), 1);
            sz_mat = zeros(numel(cur_inds), 1);
            for i = 1 : numel(cur_inds)
                tic
                a = load(outputFullpaths{cur_inds(i)}, 'im', 'mu', 'numSlice');
                sz_mat(i) = a.numSlice;
                if a.numSlice < 100
                    continue;
                end
                z_cell{i} = prctile(a.im, percent, 3);
                toc
            end
            
            z_nc = cat(3, z_cell{:});
            z_im = median(z_nc, 3);

            bg_tmpname = sprintf('%s/%sScan_Cam%s_ch%d_CAM1_background_percentile_%d_%s.tif', ...
                bg_rt, prefix, Cam(ncam), Ch(c), percent, uuid);

            writetiff(z_im, bg_tmpname);
            movefile(bg_tmpname, bg_fname);
        end
    end
end



end

