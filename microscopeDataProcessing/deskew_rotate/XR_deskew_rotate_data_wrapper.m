function XR_deskew_rotate_data_wrapper(dataPaths, varargin)
% data-level deskew/rotate wrapper, support small and large scale deskew/rotate
% 
% Author: Xiongtao Ruan (11/15/2023)


ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('dataPaths');
ip.addParameter('DSDirStr', 'DS/', @ischar);
ip.addParameter('DSRDirStr', 'DSR/', @ischar);
ip.addParameter('Deskew', true, @islogical);
ip.addParameter('Rotate', true, @islogical);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('ChannelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0', 'CamB_ch1'}, @iscell);
ip.addParameter('dz', 0.5, @isscalar);
ip.addParameter('xyPixelSize', 0.108, @isscalar);
ip.addParameter('SkewAngle', 32.45, @isscalar);
ip.addParameter('ObjectiveScan', false, @islogical);
ip.addParameter('ZstageScan', false, @islogical);
ip.addParameter('Reverse', false, @islogical);
ip.addParameter('flipZstack', false, @islogical);
ip.addParameter('parseSettingFile', false, @islogical); % use setting file to decide whether filp Z stack or not, it is  poirier over flipZstack
ip.addParameter('Crop', false, @islogical);
ip.addParameter('DSRCombined', true, @(x) islogical(x)); % combined processing 
ip.addParameter('LLFFCorrection', false, @islogical);
ip.addParameter('BKRemoval', false, @islogical);
ip.addParameter('LowerLimit', 0.4, @isnumeric); % this value is the lowest
ip.addParameter('constOffset', [], @(x) isnumeric(x)); % If it is set, use constant background, instead of background from the camera.
ip.addParameter('LSImagePaths', {'','',''}, @iscell);
ip.addParameter('BackgroundPaths', {'','',''}, @iscell);
ip.addParameter('Save16bit', false , @islogical); % saves deskewed data as 16 bit -- not for quantification
ip.addParameter('save3DStack', true , @islogical); % option to save 3D stack or not
ip.addParameter('SaveMIP', true , @islogical); % save MIP-z for ds and dsr. 
ip.addParameter('largeFile', false, @islogical);
ip.addParameter('zarrFile', false, @islogical); % use zarr file as input
ip.addParameter('saveZarr', false , @islogical); % save as zarr
ip.addParameter('BatchSize', [1024, 1024, 1024] , @isvector); % in y, x, z
ip.addParameter('BlockSize', [256, 256, 256], @isvector); % in y, x, z
ip.addParameter('zarrSubSize', [20, 20, 20], @isnumeric); % zarr subfolder size
ip.addParameter('inputBbox', [], @(x) isempty(x) || isvector(x));
ip.addParameter('taskSize', [], @isnumeric);
ip.addParameter('resample', [], @(x) isempty(x) || isnumeric(x)); % resampling after rotation 
ip.addParameter('Interp', 'linear', @(x) any(strcmpi(x, {'cubic', 'linear'})));
ip.addParameter('maskFns', {}, @iscell); % 2d masks to filter regions to deskew and rotate, in xy, xz, yz order
ip.addParameter('suffix', '', @ischar); % suffix for the folder
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('parseParfor', false, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 1, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('ConfigFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

% make sure the function is in the root of XR_Repository or LLSM5DTools. 
mpath = fileparts(which(mfilename));
repo_rt = [mpath, '/../../../'];
cd(repo_rt);
if ~exist([repo_rt, 'setup.m'], 'file')
    repo_rt = [mpath, '/../../'];
    cd(repo_rt);
end

pr = ip.Results;
DSDirStr = pr.DSDirStr;
DSRDirStr = pr.DSRDirStr;
Deskew = pr.Deskew;
Rotate = pr.Rotate;
Overwrite = pr.Overwrite;
% image related parameters
ChannelPatterns = pr.ChannelPatterns;
dz = pr.dz;
xyPixelSize = pr.xyPixelSize;
SkewAngle = pr.SkewAngle;
ObjectiveScan = pr.ObjectiveScan;
ZstageScan = pr.ZstageScan;
Reverse = pr.Reverse;
flipZstack = pr.flipZstack;
parseSettingFile = pr.parseSettingFile;
Crop = pr.Crop;
DSRCombined = pr.DSRCombined;
% flat field parameters
LLFFCorrection = pr.LLFFCorrection;
BKRemoval = pr.BKRemoval;
LowerLimit = pr.LowerLimit;
constOffset = pr.constOffset;
LSImagePaths = pr.LSImagePaths;
BackgroundPaths = pr.BackgroundPaths;
% input and output parameters
Save16bit = pr.Save16bit;
save3DStack = pr.save3DStack;
SaveMIP = pr.SaveMIP;
largeFile = pr.largeFile;
zarrFile = pr.zarrFile;
saveZarr = pr.saveZarr;
BatchSize = pr.BatchSize;
BlockSize = pr.BlockSize;
zarrSubSize = pr.zarrSubSize;
inputBbox = pr.inputBbox;
taskSize = pr.taskSize;
resample = pr.resample;
Interp = pr.Interp;
maskFns = pr.maskFns;
suffix = pr.suffix;
% job related parameters
parseCluster = pr.parseCluster;
parseParfor = pr.parseParfor;
masterCompute = pr.masterCompute;
jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
debug = pr.debug;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;

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

if ObjectiveScan
    Deskew = false;
    if strcmp(DSRDirStr, 'DSR/')
        DSRDirStr = 'Rotated/';
    end
end

if largeFile
    if ~Rotate
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

saveDS = Deskew && (~Rotate || (Rotate && ~DSRCombined));
if ~saveDS && ~Rotate
    fprintf('Both deskew and rotation are disabled, skip the processing!\n');
    return;
end

dsPaths = cell(nd, 1);
dsrPaths = cell(nd, 1);
for d = 1 : nd
    dataPath = dataPaths{d};
    
    if saveDS
        dsPath = [dataPath, DSDirStr, filesep];
        dsPaths{d, 1} = dsPath;        
    end
    if Rotate
        dsrPath = [dataPath, DSRDirStr, filesep];
        dsrPaths{d, 1} = dsrPath;        
    end

    if Overwrite
        if saveDS && exist(dsPath, 'dir')
            rmdir(dsPath, 's');            
        end
        if Rotate && exist(dsrPath, 'dir')
            rmdir(dsrPath, 's');
        end
    end
    
    if saveDS && ~exist(dsPath, 'dir')
        mkdir(dsPath);
        fileattrib(dsPath, '+w', 'g');
    end

    if Rotate && ~exist(dsrPath, 'dir')
        mkdir(dsrPath);
        fileattrib(dsrPath, '+w', 'g');
    end

    % save ds/dsr parameters
    if saveDS && ~Rotate
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
[fnames, fdinds, gfnames, partialvols, dataSizes, flipZstack_mat, latest_modify_times, FTP_inds, maskFullpaths] = ...
    XR_parseImageFilenames(dataPaths, ChannelPatterns, parseSettingFile, flipZstack, Decon, deconPaths, Streaming, minModifyTime, zarrFile);

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

    if saveDS && ~Rotate
        dsPath = dsPaths{fdind};
        dsFullpath = [dsPath, fsname, ext];
        outputFullpaths{f} = dsFullpath;
    end
    if Rotate
        dsrPath = dsrPaths{fdind};
        dsrFullpath = [dsrPath, fsname, ext];
        outputFullpaths{f} = dsrFullpath;              
    end

    if largeFile
        maskFns_str =  sprintf('{''%s''}', strjoin(maskFns, ''','''));
        func_str = sprintf(['XR_deskewRotateZarr(''%s'',%.20d,%.20d,''resultDirStr'',''%s'',', ...
            '''Overwrite'',%s,''ObjectiveScan'',%s,''SkewAngle'',%d,''Reverse'',%s,''flipZstack'',%s,', ...
            '''Crop'',%s,''DSRCombined'',%s,''Save16bit'',%s,''SaveMIP'',%s,''saveZarr'',%s,', ...
            '''BatchSize'',%s,''BlockSize'',%s,''zarrSubSize'',%s,''inputBbox'',%s,''taskSize'',%s,', ...
            '''resample'',%s,''Interp'',''%s'',''maskFns'',%s,''suffix'',''%s'',''parseCluster'',%s,', ...
            '''parseParfor'',%s,''masterCompute'',%s,''jobLogDir'',''%s'',''cpusPerTask'',%d,',...
            '''uuid'',''%s'',''debug'',%s,''mccMode'',%s,''ConfigFile'',''%s'')'], ...
            frameFullpath, xyPixelSize, dz, DSRDirStr, string(Overwrite), string(ObjectiveScan), ...
            SkewAngle, string(Reverse), string(flipZstack), string(Crop), string(DSRCombined), ...
            string(Save16bit), string(SaveMIP), string(saveZarr), strrep(mat2str(BatchSize), ' ', ','), ...
            strrep(mat2str(BlockSize), ' ', ','), strrep(mat2str(zarrSubSize), ' ', ','), ...
            strrep(mat2str(inputBbox), ' ', ','), strrep(mat2str(taskSize), ' ', ','), ...
            strrep(mat2str(resample), ' ', ','), Interp, maskFns_str, suffix, string(parseCluster), ...
            string(parseParfor), string(masterCompute), jobLogDir, cpusPerTask, uuid, string(debug), ...
            string(mccMode), ConfigFile);
    else
        if LLFFCorrection
            % LLFFMapping =  ~cellfun(@isempty, regexpi(fname, ChannelPatterns));
            % change to contains.m to unify the matching
            LLFFMapping =  cellfun(@(x) contains(frameFullpath, x), ChannelPatterns);
            LSImage = LSImagePaths{LLFFMapping};
            BackgroundImage = BackgroundPaths{LLFFMapping};
        else
            LSImage = '';
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

        func_str = sprintf(['XR_deskewRotateFrame(%s,%.20d,%.20d,''SkewAngle'',%.20d,', ...
            '''ObjectiveScan'',%s,''ZstageScan'',%s,''Reverse'',%s,''LLFFCorrection'',%s,', ...
            '''BKRemoval'',%s,''LowerLimit'',%.20d,''constOffset'',[%s],''LSImage'',''%s'',', ...
            '''BackgroundImage'',''%s'',''Rotate'',%s,''resample'',[%s],''DSRCombined'',%s,', ...
            '''inputBbox'',%s,''flipZstack'',%s,''Save16bit'',%s,''save3DStack'',%s,''saveZarr'',%s)'], ...
            ds_input_str, xyPixelSize, dz_f, SkewAngle, string(ObjectiveScan), string(ZstageScan), ...
            string(Reverse),  string(LLFFCorrection), string(BKRemoval), LowerLimit, ...
            num2str(constOffset, '%0.10f'),  LSImage, BackgroundImage, string(Rotate), ...
            strrep(mat2str(resample), ' ', ','), string(DSRCombined), strrep(mat2str(inputBbox), ' ', ','), ...
            string(flipZstack), string(Save16bit(1)), string(save3DStack), string(saveZarr));
    end
    
    funcStrs{f} = func_str;
end

% estimate memory usuage
dsz = getImageSize(inputFullpaths{1});
if largeFile && ~masterCompute
    memAllocate = prod([min(BatchSize(1), dsz(1)), dsz(2), dsz(3)]) * 4 / 2^30 * 10;
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
            'mccMode', mccMode, 'ConfigFile', ConfigFile);
    end
end

if ~all(is_done_flag)
    error('%d / %d failed', sum(~is_done_flag), numel(is_done_flag));
end

end

