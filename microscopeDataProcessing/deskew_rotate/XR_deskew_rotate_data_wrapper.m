function XR_deskew_rotate_data_wrapper(dataPaths, varargin)
% dataset level deskew/rotate wrapper, support small and large scale deskew/rotate
% 
%
% Required inputs:
%           dataPaths : char or cell array. Directory paths for the datasets. Either a string for a single dataset or a cell array of paths for several datasets with same experimental settings.
% 
% Parameters (as 'specifier'-value pairs): 
%           DSDirName : char (default: 'DS'). Deskew result directory name under dataPaths.
%          DSRDirName : char (default: 'DSR'). Deskew/rotation result directory name under dataPaths.
%              deskew : true|false (default: true). Deskew the data.
%              rotate : true|false (default: true). Rotate deskewed data.
%           overwrite : true|false (default: false). Overwrite existing results
%     channelPatterns : a cell array (default: {'CamA_ch0', 'CamA_ch1', 'CamB_ch0', 'CamB_ch1'}).  Channel identifiers for included channels. 
%                  dz : a number (default: 0.5). Scan interval in um.
%         xyPixelSize : a number (default: 0.108). Pixel size in um.
%           skewAngle : a number (default: 32.45). Skew angle (in degree) of the stage.
%       objectiveScan : true|false (default: false). Objective scan.
%          zStageScan : true|false (default: false). Z stage scan (orthogonal to objective scan).
%             reverse : true|false (default: false). Inverse direction of z axis. 
%          flipZstack : true|false (default: false). Flip z stacks.
%    parseSettingFile : true|false (default: false). Use the setting file to decide whether filp z stacks or not.
%                crop : true|false (default: false). Crop the data to smaller for rotation (removing extra empty regions).
%         DSRCombined : true|false (default: true). Use combined processing for deskew and rotation.
%        FFCorrection : true|false (default: false). Flat-field correction.
%           BKRemoval : true|false (default: false). Remove background during flat-field correction.
%          lowerLimit : a number between 0 and 1 (default: 0.4). Lower limit to cap the flat field image. 
%         constOffset : empty or a number (default: []). If empty, add the background; if not, add the const offset after flat field correction.
%        FFImagePaths : empty or a cell array of paths for corresponding channels (default: {'', '', ''}). Flat field image paths.
%     backgroundPaths : empty or a cell array of paths for corresponding channels (default: {'', '', ''}). Background image paths.
%           save16bit : true|false (default: true). Save 16bit result for deskew/rotate. 
%         save3DStack : true|false (default: true). Save 3D stack for DS and DSR. 
%             saveMIP : true|false (default: true). Save MIPs after deskew/rotate.
%           largeFile : true|false (default: false). Use large scale deskew/rotation strategy.
%            zarrFile : true|false (default: false). Use Zarr file as input.
%            saveZarr : true|false (default: false). Save results as Zarr files.
%           batchSize : 1x3 vector (default: [1024, 1024, 1024]). Batch size per stitching task.
%           blockSize : 1x3 vector (default: [256, 256, 256]). Block/chunk size for zarr output.
%           inputBbox : empty or 1x6 vector (default: []). Input bounding box for crop. Definiation: [ymin, xmin, zmin, ymax, xmax, zmax].
%            taskSize : empty or a number (default: []). Number of tasks per job for large-scale deskew/rotation.
%      resampleFactor : empty or 1x1, 1x2 or 1x3 vector (default: []). Resampling factor after rotation. Empty: no resampling; axis order yxz.
%        interpMethod : 'linear'|'cubic' (default: 'linear'). Interpolation method for geometric transformations in deskew and rotation.
%       maskFullpaths : empty or a cell array (default: {}). 2D projective masks to filter regions for large scale deskew/rotation, in xy, xz, and yz order.
%        parseCluster : true|false (default: true). Use slurm cluster for the processing.
%         parseParfor : true|false (default: false). Use matlab parfor for paralle processing.
%       masterCompute : true|false (default: true). Master job node is involved in the processing.
%           jobLogDir : char (default: '../job_logs'). Path for the slurm job logs.
%         cpusPerTask : a number (default: 1). The number of cpu cores per task for slurm job submission.
%                uuid : empty or a uuid string (default: ''). uuid string as part of the temporate result paths.
%               debug : true|false (default: false). Debug mode. Not actually used in this function. Reserved for future use.
%             mccMode : true|false (default: false). Use mcc mode.
%          configFile : empty or char (default: ''). Path for the config file for job submission.
%
%
% Author: Xiongtao Ruan (11/15/2023)


ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addParameter('DSDirName', 'DS', @ischar);
ip.addParameter('DSRDirName', 'DSR', @ischar);
ip.addParameter('deskew', true, @islogical);
ip.addParameter('rotate', true, @islogical);
ip.addParameter('overwrite', false, @islogical);
ip.addParameter('channelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0', 'CamB_ch1'}, @iscell);
ip.addParameter('dz', 0.5, @isscalar);
ip.addParameter('xyPixelSize', 0.108, @(x) isvector(x) && numel(x) <= 2);
ip.addParameter('skewAngle', 32.45, @isscalar);
ip.addParameter('objectiveScan', false, @islogical);
ip.addParameter('zStageScan', false, @islogical);
ip.addParameter('reverse', false, @islogical);
ip.addParameter('flipZstack', false, @islogical);
ip.addParameter('parseSettingFile', false, @islogical);
ip.addParameter('crop', false, @islogical);
ip.addParameter('DSRCombined', true, @(x) islogical(x));
ip.addParameter('FFCorrection', false, @islogical);
ip.addParameter('BKRemoval', false, @islogical);
ip.addParameter('lowerLimit', 0.4, @isnumeric);
ip.addParameter('constOffset', [], @(x) isnumeric(x));
ip.addParameter('FFImagePaths', {'','',''}, @iscell);
ip.addParameter('backgroundPaths', {'','',''}, @iscell);
ip.addParameter('save16bit', true , @islogical);
ip.addParameter('save3DStack', true , @islogical);
ip.addParameter('saveMIP', true , @islogical);
ip.addParameter('largeFile', false, @islogical);
ip.addParameter('zarrFile', false, @islogical);
ip.addParameter('saveZarr', false , @islogical);
ip.addParameter('batchSize', [1024, 1024, 1024] , @isvector);
ip.addParameter('blockSize', [256, 256, 256], @isvector);
ip.addParameter('inputBbox', [], @(x) isempty(x) || isvector(x));
ip.addParameter('taskSize', [], @isnumeric);
ip.addParameter('resampleFactor', [], @(x) isempty(x) || isnumeric(x));
ip.addParameter('interpMethod', 'linear', @(x) any(strcmpi(x, {'cubic', 'linear'})));
ip.addParameter('maskFullpaths', {}, @iscell);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('parseParfor', false, @islogical);
ip.addParameter('masterCompute', true, @islogical);
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 1, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('configFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

% make sure the function is in the root of XR_Repository or PetaKit5D. 
mpath = fileparts(which(mfilename));
repo_rt = [mpath, '/../../../'];
cd(repo_rt);
if ~exist([repo_rt, 'setup.m'], 'file')
    repo_rt = [mpath, '/../../'];
    cd(repo_rt);
end

pr = ip.Results;
DSDirName = pr.DSDirName;
DSRDirName = pr.DSRDirName;
deskew = pr.deskew;
rotate = pr.rotate;
overwrite = pr.overwrite;
% image related parameters
channelPatterns = pr.channelPatterns;
dz = pr.dz;
xyPixelSize = pr.xyPixelSize;
skewAngle = pr.skewAngle;
objectiveScan = pr.objectiveScan;
zStageScan = pr.zStageScan;
reverse = pr.reverse;
flipZstack = pr.flipZstack;
parseSettingFile = pr.parseSettingFile;
crop = pr.crop;
DSRCombined = pr.DSRCombined;
% flat field parameters
FFCorrection = pr.FFCorrection;
BKRemoval = pr.BKRemoval;
lowerLimit = pr.lowerLimit;
constOffset = pr.constOffset;
FFImagePaths = pr.FFImagePaths;
backgroundPaths = pr.backgroundPaths;
% input and output parameters
save16bit = pr.save16bit;
save3DStack = pr.save3DStack;
saveMIP = pr.saveMIP;
largeFile = pr.largeFile;
zarrFile = pr.zarrFile;
saveZarr = pr.saveZarr;
batchSize = pr.batchSize;
blockSize = pr.blockSize;
inputBbox = pr.inputBbox;
taskSize = pr.taskSize;
resampleFactor = pr.resampleFactor;
interpMethod = pr.interpMethod;
maskFullpaths = pr.maskFullpaths;
% job related parameters
parseCluster = pr.parseCluster;
parseParfor = pr.parseParfor;
masterCompute = pr.masterCompute;
jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
debug = pr.debug;
mccMode = pr.mccMode;
configFile = pr.configFile;

% suppress directory exists warning
warning('off', 'MATLAB:MKDIR:DirectoryExists');

if ischar(dataPaths)
    dataPaths = {dataPaths};
end

nd = numel(dataPaths);
for d = 1 : nd
    dataPath = dataPaths{d};
    if ~strcmp(dataPath(end), filesep)
        dataPaths{d} = [dataPath, filesep];
    end
end

if isempty(uuid)
    uuid = get_uuid();
end

% check if a slurm-based computing cluster exists
if parseCluster
    [parseCluster, job_log_fname, job_log_error_fname] = checkSlurmCluster(dataPath, jobLogDir);
end

if objectiveScan
    deskew = false;
    if strcmp(DSRDirName, 'DSR/')
        DSRDirName = 'rotated/';
    end
end

if largeFile
    if ~rotate
        error('Currently only combined DSR is supported for large scale deskew/rotation!\n')
    end
    DSRCombined = true;
    saveZarr = true;
end

% save zarr 
ext = '.tif';
if saveZarr
    ext = '.zarr';
end

saveDS = deskew && (~rotate || (rotate && ~DSRCombined));
if ~saveDS && ~rotate
    fprintf('Both deskew and rotation are disabled, skip the processing!\n');
    return;
end

dsPaths = cell(nd, 1);
dsrPaths = cell(nd, 1);
for d = 1 : nd
    dataPath = dataPaths{d};
    
    if saveDS
        dsPath = [dataPath, DSDirName, filesep];
        dsPaths{d, 1} = dsPath;        
    end
    if rotate
        dsrPath = [dataPath, DSRDirName, filesep];
        dsrPaths{d, 1} = dsrPath;        
    end

    if overwrite
        if saveDS && exist(dsPath, 'dir')
            rmdir(dsPath, 's');            
        end
        if rotate && exist(dsrPath, 'dir')
            rmdir(dsrPath, 's');
        end
    end
    
    if saveDS && ~exist(dsPath, 'dir')
        mkdir(dsPath);
        fileattrib(dsPath, '+w', 'g');
    end

    if rotate && ~exist(dsrPath, 'dir')
        mkdir(dsrPath);
        fileattrib(dsrPath, '+w', 'g');
    end

    % save ds/dsr parameters
    if saveDS && ~rotate
        save('-v7.3', [dsPath, '/parameters.mat'], 'pr');
        s = jsonencode(pr, PrettyPrint=true);
        fid = fopen([dsPath, '/parameters.json'], 'w');
        fprintf(fid, s);
        fclose(fid);
    else
        save('-v7.3', [dsrPath, '/parameters.mat'], 'pr');
        s = jsonencode(pr, PrettyPrint=true);
        fid = fopen([dsrPath, '/parameters.json'], 'w');
        fprintf(fid, s);
        fclose(fid);
    end
end


%% check existing files and parse channels
Decon = false;
deconPaths = {};
Streaming = false;
minModifyTime = 1;
[fnames, fdinds, gfnames, partialvols, dataSizes, flipZstack_mat, latest_modify_times, FTP_inds] = ...
    XR_parseImageFilenames(dataPaths, channelPatterns, parseSettingFile, flipZstack, Decon, deconPaths, Streaming, minModifyTime, zarrFile);

nF = numel(fnames);

% construct function strings and input and output files
inputFullpaths = cell(nF, 1);
outputFullpaths = cell(nF, 1);
funcStrs = cell(nF, 1);

for f = 1 : nF

    % first deskew and rotation
    fname = fnames{f};
    [~, fsname] = fileparts(fname);
    fdind = fdinds(f);
    partialvol = partialvols(f);
    gfname = gfnames{f};
    dataPath = dataPaths{fdind};
    dz_f = dz;

    frameFullpath = [dataPath, fname];
    inputFullpaths{f} = frameFullpath;

    if saveDS && ~rotate
        dsPath = dsPaths{fdind};
        dsFullpath = [dsPath, fsname, ext];
        outputFullpaths{f} = dsFullpath;
    end
    if rotate
        dsrPath = dsrPaths{fdind};
        dsrFullpath = [dsrPath, fsname, ext];
        outputFullpaths{f} = dsrFullpath;              
    end

    if largeFile
        maskFullpaths_str =  sprintf('{''%s''}', strjoin(maskFullpaths, ''','''));
        func_str = sprintf(['XR_deskewRotateZarr(''%s'',%.20d,%.20d,''resultDirStr'',''%s'',', ...
            '''overwrite'',%s,''objectiveScan'',%s,''skewAngle'',%d,''reverse'',%s,''flipZstack'',%s,', ...
            '''crop'',%s,''DSRCombined'',%s,''save16bit'',%s,''saveMIP'',%s,''saveZarr'',%s,', ...
            '''batchSize'',%s,''blockSize'',%s,''inputBbox'',%s,''taskSize'',%s,', ...
            '''resampleFactor'',%s,''interpMethod'',''%s'',''maskFullpaths'',%s,''parseCluster'',%s,', ...
            '''parseParfor'',%s,''masterCompute'',%s,''jobLogDir'',''%s'',''cpusPerTask'',%d,',...
            '''uuid'',''%s'',''debug'',%s,''mccMode'',%s,''configFile'',''%s'')'], ...
            frameFullpath, xyPixelSize, dz, DSRDirName, string(overwrite), string(objectiveScan), ...
            skewAngle, string(reverse), string(flipZstack), string(crop), string(DSRCombined), ...
            string(save16bit), string(saveMIP), string(saveZarr), strrep(mat2str(batchSize), ' ', ','), ...
            strrep(mat2str(blockSize), ' ', ','), strrep(mat2str(inputBbox), ' ', ','), ...
            strrep(mat2str(taskSize), ' ', ','), strrep(mat2str(resampleFactor), ' ', ','), ...
            interpMethod, maskFullpaths_str, string(parseCluster), string(parseParfor), ...
            string(masterCompute), jobLogDir, cpusPerTask, uuid, string(debug), ...
            string(mccMode), configFile);
    else
        if FFCorrection
            % LLFFMapping =  ~cellfun(@isempty, regexpi(fname, channelPatterns));
            % change to contains.m to unify the matching
            LLFFMapping =  cellfun(@(x) contains(frameFullpath, x), channelPatterns);
            FFImage = FFImagePaths{LLFFMapping};
            BackgroundImage = backgroundPaths{LLFFMapping};
        else
            FFImage = '';
            BackgroundImage = '';
        end
        
        flipZstack = flipZstack_mat(f);
        
        % set up input file for either single volume file or a
        % group of files
        ds_input_path = {frameFullpath};
        if partialvol
            ds_input_path = cellfun(@(x) [dataPath, x], gfname, 'unif', 0);
        end
        ds_input_str = sprintf('{''%s''}', strjoin(ds_input_path, ''','''));

        func_str = sprintf(['XR_deskewRotateFrame(%s,%.20d,%.20d,''DSDirName'',''%s'',', ...
            '''DSRDirName'',''%s'',''skewAngle'',%.20d,''objectiveScan'',%s,''zStageScan'',%s,', ...
            '''reverse'',%s,''FFCorrection'',%s,''BKRemoval'',%s,''lowerLimit'',%.20d,', ...
            '''constOffset'',[%s],''FFImage'',''%s'',''BackgroundImage'',''%s'',''rotate'',%s,', ...
            '''resampleFactor'',[%s],''DSRCombined'',%s,''inputBbox'',%s,''flipZstack'',%s,', ...
            '''save16bit'',%s,''save3DStack'',%s,''saveZarr'',%s,''interpMethod'',''%s'')'], ...
            ds_input_str, xyPixelSize, dz_f, DSDirName, DSRDirName, skewAngle, ...
            string(objectiveScan), string(zStageScan), string(reverse),  string(FFCorrection), ...
            string(BKRemoval), lowerLimit, num2str(constOffset, '%0.10f'), FFImage, ...
            BackgroundImage, string(rotate), strrep(mat2str(resampleFactor), ' ', ','), ...
            string(DSRCombined), strrep(mat2str(inputBbox), ' ', ','), string(flipZstack), ...
            string(save16bit(1)), string(save3DStack), string(saveZarr), interpMethod);
    end
    
    funcStrs{f} = func_str;
end

% estimate memory usuage
dsz = getImageSize(inputFullpaths{1});
if largeFile && ~masterCompute
    memAllocate = prod([min(batchSize(1), dsz(1)), dsz(2), dsz(3)]) * 4 / 2^30 * 10;
    cpusPerTask = 1;
else
    memAllocate = prod(dsz) * 4 / 2^30 * 10;                        
end

% cluster setting
is_done_flag = false;
maxTrialNum = 2;

% retry with more resources
for i = 1 : 3
    if  ~all(is_done_flag)
        is_done_flag = generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, ...
            funcStrs, 'cpusPerTask', cpusPerTask * 2^(i - 1),  'memAllocate', memAllocate * 2^(i - 1), ...
            'maxTrialNum', maxTrialNum, 'masterCompute', masterCompute, 'parseCluster', parseCluster, ...
            'mccMode', mccMode, 'configFile', configFile);
    end
end

if ~all(is_done_flag)
    error('%d / %d failed', sum(~is_done_flag), numel(is_done_flag));
end

end

