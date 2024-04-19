function [] = XR_MIP_wrapper(dataPaths, varargin)   
% generate MIP for dataset

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addParameter('axis', [0, 0, 1], @isnumeric); % y, x, z
ip.addParameter('ChannelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0', 'CamB_ch1'}, @iscell);
ip.addParameter('zarrFile', false, @islogical); % use zarr file as input
ip.addParameter('largeZarr', false, @islogical); % use zarr file as input
ip.addParameter('BatchSize', [2048, 2048, 2048] , @isvector); % in y, x, z
ip.addParameter('Save16bit', true, @islogical);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('parseParfor', false, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('cpusPerTask', 3, @isscalar);
ip.addParameter('jobLogDir', '../job_logs/', @ischar);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('ConfigFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
axis = pr.axis;
ChannelPatterns =  pr.ChannelPatterns;
zarrFile = pr.zarrFile;
largeZarr = pr.largeZarr;
BatchSize = pr.BatchSize;
Save16bit = pr.Save16bit;
parseCluster = pr.parseCluster;
parseParfor = pr.parseParfor;
jobLogDir = pr.jobLogDir;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;

uuid = pr.uuid;
% uuid for the job
if isempty(uuid)
    uuid = get_uuid();
end
debug = pr.debug;

if ischar(dataPaths)
    dataPaths = {dataPaths};
end

if Save16bit
    dtype = 'uint16';
else
    dtype = 'single';
end

nd = numel(dataPaths);
resultPaths = cell(nd, 1);
for d = 1 : nd
    dataPath = dataPaths{d};
    resultPath = [dataPath, '/MIPs/'];
    mkdir(resultPath);
    
    dataPaths{d} = simplifyPath(dataPath);
    resultPaths{d} = simplifyPath(resultPath);    
    fileattrib(resultPath, '+w', 'g');

    save('-v7.3', [resultPath, '/parameters.mat'], 'pr');
    s = jsonencode(pr, PrettyPrint=true);
    fid = fopen([resultPath, '/parameters.json'], 'w');
    fprintf(fid, s);
    fclose(fid);        
end

% check if a slurm-based computing cluster exists
if parseCluster
    [parseCluster, job_log_fname, job_log_error_fname] = checkSlurmCluster(dataPaths{1}, jobLogDir);
end

% parse image filenames
[~, fsns, fd_inds, filepaths] = parseImageFilenames(dataPaths, zarrFile, ChannelPatterns);
nF = numel(fsns);
sz = getImageSize(filepaths{1});

%% use generic framework for the MIP computing

frameFullpaths = filepaths;
MIPFullpaths = cell(nF, 1);
func_strs = cell(nF, 1);

% use the last select axis as flag if more than one axis for MIP
axis_strs = {'y', 'x', 'z'};
fidx = find(axis, 1, 'last');
axis_str = axis_strs{fidx};

for f = 1 : nF
    frameFullpath = frameFullpaths{f};

    d = fd_inds(f);
    resultPath = resultPaths{d};
    MIPFullpath = sprintf('%s/%s_MIP_%s.tif', resultPath, fsns{f}, axis_str);
    MIPFullpaths{f} = MIPFullpath;
    
    if zarrFile
        if largeZarr
            func_strs{f} = sprintf(['XR_MIP_zarr(''%s'',''axis'',%s,''BatchSize'',%s,', ...
                '''parseCluster'',%s,''parseParfor'',%s,''jobLogDir'',''%s'',', ...
                '''mccMode'',%s,''ConfigFile'',''%s'')'], frameFullpath, mat2str_comma(axis), ...
                mat2str_comma(BatchSize), string(parseCluster), ...
                string(parseParfor), jobLogDir, string(mccMode), ConfigFile);
        else
            func_strs{f} = sprintf(['saveMIP_zarr(''%s'',''%s'',''%s'',%s)'], frameFullpath, ...
                MIPFullpath, dtype, mat2str_comma(axis));
        end
    else
        func_strs{f} = sprintf(['saveMIP_tiff(''%s'',''%s'',''dtype'',''%s'',''axis'',%s)'], ...
            frameFullpath, MIPFullpath, dtype, mat2str_comma(axis));
    end
end

if parseParfor && ~largeZarr
    matlab_parfor_generic_computing_wrapper(frameFullpaths, MIPFullpaths, func_strs, 'GPUJob', false, 'nworker', 12)
end
taskBatchNum = max(1, round(nF /1000));
memAllocate = prod(sz) * 4 / 1024^3 * 4;
generic_computing_frameworks_wrapper(frameFullpaths, MIPFullpaths, func_strs, ...
    parseCluster=parseCluster, jobLogDir=jobLogDir, masterCompute=masterCompute, ...
    taskBatchNum=taskBatchNum, cpusPerTask=cpusPerTask, memAllocate=memAllocate, ...
    mccMode=mccMode, ConfigFile=ConfigFile);

end

