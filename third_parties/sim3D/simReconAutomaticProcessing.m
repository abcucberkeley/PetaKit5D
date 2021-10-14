function [] = simReconAutomaticProcessing(dataPaths, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));


ip.addParameter('PSFs',{''},@(x) ischar(x) || iscell(x));

ip.addParameter('Deskew',true,@islogical);
ip.addParameter('Recon',true,@islogical);
ip.addParameter('Streaming',false,@islogical);

ip.addParameter('xyPixelSize',.108,@isnumeric); % typical value: 0.1
ip.addParameter('dz',.5,@isnumeric); % typical value: 0.2-0.5
ip.addOptional('SkewAngle', 32.45, @isscalar);
ip.addOptional('Reverse', true, @islogical);

ip.addParameter('Rotate',false,@islogical); % Rotate after deskew

ip.addParameter('ChannelPatterns',{'CamA','CamB'},@iscell);

ip.addParameter('islattice', true, @islogical); %Flag to indicate if this is light sheet data
ip.addParameter('NA_det', 1.0, @isnumeric);
ip.addParameter('NA_ext', 0.55, @isnumeric);
ip.addParameter('nimm', 1.33, @isnumeric);
ip.addParameter('wvl_em', .605, @isnumeric);
ip.addParameter('wvl_ext', .560, @isnumeric);
ip.addParameter('w', 5e-3, @isnumeric); %Wiener coefficient for regularization
ip.addParameter('apodize', true, @islogical); %Flag to indicate whether or not to apodize the final data

ip.addParameter('DS', true, @islogical);
ip.addParameter('nphases', 5, @isnumeric);
ip.addParameter('norders', 5, @isnumeric);
ip.addParameter('norientations', 1, @isnumeric);
ip.addParameter('lattice_period', 1.2021, @isnumeric); %Lattice period in microns - this is the coarsest period
ip.addParameter('lattice_angle', [pi/2], @isnumeric); %Angle parellel to pattern modulation (assuming horizontal is zero)
ip.addParameter('phase_step', .232, @isnumeric); %Phase step in microns
ip.addParameter('pxl_dim_data', [0.11,0.11,0.3*sind(32.4)], @isnumeric); %Voxel dimensions of the image in microns - note, stored as [dy,dx,dz]
ip.addParameter('pxl_dim_PSF', [0.11,0.11,0.2*sind(32.4)], @isnumeric); %Voxel dimensions of the PSF in microns - note, stored as [dy,dx,dz]
ip.addParameter('Background', 105, @isnumeric);

ip.addParameter('normalize_orientations', false, @islogical); %Flag to indicate whether or not to normalize total intensity for each orientation
ip.addParameter('perdecomp', false, @islogical); %Flag to indicate whether or not to use periodic/smooth decomposition to reduce edge effects
ip.addParameter('edgeTaper', false, @islogical); %Flag to indicate whether or not to window the data to reduce edge effects
ip.addParameter('edgeTaperVal', 0.1, @isnumeric); %Roll-off parameter for Tukey windowing

ip.addParameter('useGPU', true, @islogical);

ip.addParameter('Overwrite', false , @islogical);
ip.addParameter('Save16bit', false , @islogical);
ip.addParameter('flipZstack', false, @islogical);
ip.addParameter('EdgeSoften', 5, @isnumeric); % # ofxy px to soften
ip.addParameter('zEdgeSoften', 2, @isnumeric); % # ofxy px to soften
ip.addParameter('Crop', [], @isnumeric); % requires lower and higher values for cropping
ip.addParameter('zFlip', false, @islogical);
ip.addParameter('GenMaxZproj', [0,0,1] , @isnumeric);
ip.addParameter('ResizeImages', [] , @isnumeric);
ip.addParameter('EdgeErosion', 0 , @isnumeric); % erode edges for certain size.
ip.addParameter('ErodeMaskfile', '', @ischar); % erode edges file
ip.addParameter('SaveMaskfile', false, @islogical); % save mask file for common eroded mask
ip.addParameter('ChunkSize', [250, 250, 250] , @isvector); % in y, x, z
ip.addParameter('Overlap', 10, @isnumeric); % block overlap
ip.addParameter('maxSubVolume', 3e8, @isnumeric);
ip.addParameter('CPUMaxMem', 500, @isnumeric); % CPU Memory in Gb

ip.addParameter('splitJobsByChannelPattern', false, @islogical);

ip.addParameter('largeFile', false, @islogical);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', false, @islogical); % master node participate in the task computing.
ip.addParameter('masterCPU', false, @islogical); % master node is a cpu node, which is just for large file deconvolution.
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 5, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxModifyTime', 10, @isnumeric); % the maximum during of last modify time of a file, in minute.
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 2, @isnumeric);

ip.parse(dataPaths, varargin{:});

if ischar(dataPaths)
    dataPaths = {dataPaths};
end

pr = ip.Results;

PSFs = pr.PSFs;

if ischar(PSFs)
    PSFs = {PSFs};
end

Deskew = pr.Deskew;
Recon = pr.Recon;
Streaming = pr.Streaming;

xyPixelSize = pr.xyPixelSize;
dz = pr.dz;
SkewAngle = pr.SkewAngle;
Reverse = pr.Reverse;

Rotate = pr.Rotate;

ChannelPatterns = pr.ChannelPatterns;

islattice = pr.islattice;
NA_det = pr.NA_det;
NA_ext = pr.NA_ext;
nimm = pr.nimm;
wvl_em = pr.wvl_em;
wvl_ext = pr.wvl_ext;
w = pr.w;
apodize = pr.apodize;

DS = pr.DS;
nphases = pr.nphases;
norders = pr.norders;
norientations = pr.norientations;
lattice_period = pr.lattice_period;
lattice_angle = pr.lattice_angle;
phase_step = pr.phase_step;
pxl_dim_data = pr.pxl_dim_data;
pxl_dim_PSF = pr.pxl_dim_PSF;
Background = pr.Background;

normalize_orientations = pr.normalize_orientations;
perdecomp = pr.perdecomp;
edgeTaper = pr.edgeTaper;
edgeTaperVal = pr.edgeTaperVal;

useGPU = pr.useGPU;

flipZstack = pr.flipZstack;
Save16bit = pr.Save16bit;

EdgeErosion = pr.EdgeErosion;

ErodeMaskfile = pr.ErodeMaskfile;
if ~isempty(ErodeMaskfile) && exist(ErodeMaskfile, 'file')
    EdgeErosion = 0; % set EdgeErosion length as 0.
end
SaveMaskfile = pr.SaveMaskfile;

OL = pr.Overlap;
ChunkSize = pr.ChunkSize;
maxSubVolume = pr.maxSubVolume;

splitJobsByChannelPattern = pr.splitJobsByChannelPattern;

largeFile = pr.largeFile;
parseCluster = pr.parseCluster;
jobLogDir = pr.jobLogDir;
masterCompute = pr.masterCompute;
maxModifyTime = pr.maxModifyTime;
maxTrialNum = pr.maxTrialNum;
unitWaitTime = pr.unitWaitTime;
uuid = pr.uuid;
if isempty(uuid)
    uuid = get_uuid();
end

if(Rotate)
    folStr = 'DSR';
else
    folStr = 'DS';
end

firstTime = true;
latest_modify_time = 0;
allFullPaths = {''};

% Guess number of workers
%workers = cell(1,length(ChannelPatterns)*length(dataPaths));
workers = {};
cWorker = 1;
cStates = 'finished';

while(firstTime || (~isempty(workers) && ~all(strcmp(cStates,'finished'))) || (Streaming && (latest_modify_time < maxModifyTime)))
    if(firstTime) 
        firstTime = false;
    else
        % Sleep main thread before checking again
        pause(10);
    end
for i = 1:numel(dataPaths)
    for cPatt = 1:numel(ChannelPatterns)
        fnames = dir([dataPaths{i} filesep '*' ChannelPatterns{cPatt} '*.tif']);
        
        if(isempty(fnames))
            fprintf('No Files found for channel pattern: ''%s''. Skipping to next pattern.\n',ChannelPatterns{cPatt});
        end
        
        if(Streaming)
            allTimes = [fnames.datenum];
            allTimes = (datenum(clock)-allTimes)*24*60;
            if(min(allTimes) < unitWaitTime)
                newestIndex = allTimes == min(allTimes);
                fnames(newestIndex) = [];
            end
            % If there are no more files. Continue to next pattern 
            if(isempty(fnames))
                continue;
            end
        end
        
        % Check that jobs have not yet been submitted for these files
        match = ismember(strcat([fnames(1).folder filesep],{fnames.name}),allFullPaths);
        fnames(match) = [];
        
        % If there are no more files. Continue to next pattern 
        if(isempty(fnames))
            continue;
        end
        
        % Add files to list
        allFullPaths = horzcat(allFullPaths,strcat([fnames(1).folder filesep],{fnames.name}));
        
        fnames = {fnames.name};
        inputFullpaths = cell(numel(fnames), 1);
        outputFullpaths = cell(numel(fnames), 1);
        funcStrs = cell(numel(fnames), 1);
        
        if parseCluster
            if  ~exist(jobLogDir, 'dir')
                warning('The job log directory does not exist, use %s/job_logs as job log directory.', dataPaths{i})
                jobLogDir = sprintf('%s/job_logs', dataPaths{i});
                mkdir(jobLogDir);
            end
            job_log_fname = [jobLogDir, '/job_%A_%a.out'];
            job_log_error_fname = [jobLogDir, '/job_%A_%a.err'];
        end
        
        alreadyFinished = true;
        for j = 1: numel(fnames)
            [pathstr, fsname, ext] = fileparts(fnames{j});
            dataFullpath = [dataPaths{i} filesep fnames{j}];
            dataDSFullpath = [dataPaths{i} filesep folStr filesep fsname ext];
            
            if alreadyFinished
                if ~exist(dataDSFullpath,'file')
                    alreadyFinished = false;
                end
            end
            
            inputFullpaths{j} = dataFullpath;
            outputFullpaths{j} = dataDSFullpath;
            
            
            funcStrs{j} =  sprintf(['deskewPhasesFrame(''%s'',%.10f,%.10f,''SkewAngle'',%.10f,''Reverse'',%s,''nphases''', ...
                ',%.10f,''Rotate'',%s)'], dataFullpath, xyPixelSize, dz, SkewAngle, string(Reverse), nphases, string(Rotate));
        end
        
        if alreadyFinished
            fprintf('Skipping Deskew for pattern ''%s'' in folder ''%s'' because it already exists.\n',ChannelPatterns{cPatt},dataPaths{i});
            continue;
        end
        
        %if useGPU
        %   maxJobNum = inf;
        %  cpusPerTask = 5;
        % cpuOnlyNodes = false;
        %taskBatchNum = 5;
        %SlurmParam = '-p abc --qos abc_normal -n1 --mem=167G --gres=gpu:1';
        %else
        maxJobNum = inf;
        cpusPerTask = 24;
        cpuOnlyNodes = true;
        taskBatchNum = 1;
        SlurmParam = '-p abc --qos abc_normal -n1 --mem-per-cpu=21418M';
        %end
        
        
        
        if(~isempty(inputFullpaths))
            fprintf('Attempting Deskew on %d files for pattern ''%s'' in folder ''%s''\n',length(inputFullpaths),ChannelPatterns{cPatt},dataPaths{i})
            [workers{1,cWorker}] = parfeval(@slurm_cluster_generic_computing_wrapper,1,inputFullpaths, outputFullpaths, ...
                funcStrs, 'cpusPerTask', cpusPerTask, 'cpuOnlyNodes', cpuOnlyNodes, 'SlurmParam', SlurmParam, ...
                'maxJobNum', maxJobNum, 'taskBatchNum', taskBatchNum, 'masterCompute', masterCompute, 'parseCluster', parseCluster);
            workers{2,cWorker} = sprintf('Finished Deskew on %d files for pattern ''%s'' in folder ''%s''\n',length(inputFullpaths),ChannelPatterns{cPatt},dataPaths{i});
            cWorker = cWorker+1;
            %{
            is_done_flag= slurm_cluster_generic_computing_wrapper(inputFullpaths, outputFullpaths, ...
                funcStrs, 'cpusPerTask', cpusPerTask, 'cpuOnlyNodes', cpuOnlyNodes, 'SlurmParam', SlurmParam, ...
                'maxJobNum', maxJobNum, 'taskBatchNum', taskBatchNum, 'masterCompute', masterCompute, 'parseCluster', parseCluster);

            if ~all(is_done_flag)
                slurm_cluster_generic_computing_wrapper(inputFullpaths, outputFullpaths, ...
                    funcStrs, 'cpusPerTask', cpusPerTask, 'cpuOnlyNodes', cpuOnlyNodes, 'SlurmParam', SlurmParam, ...
                    'maxJobNum', maxJobNum, 'taskBatchNum', taskBatchNum, 'masterCompute', masterCompute, 'parseCluster', parseCluster);
            end
            %}
        end
    end
    
    
end

% Find the time of the most recently modified file
if Streaming
    latest_modify_time = inf;
    for i = 1:numel(dataPaths)
        for cPatt = 1:numel(ChannelPatterns)
            dir_info = dir([dataPaths{i} filesep '*' ChannelPatterns{cPatt} '*.tif']);

            last_modify_time = (datenum(clock) - [dir_info.datenum]) * 24 * 60;
            lowestTimeI = min(last_modify_time);
            latest_modify_time = min(latest_modify_time,lowestTimeI);
        end
    end
end

indices = {};
for i = 1:cWorker-1
    if strcmp(workers{1,i}.State,'finished')
       fprintf(workers{2,i});
       indices{end+1} = i;
       %workers(:,i) = [];
       %cWorker = cWorker-1;
       %i = i-1;
    end
end

workers(:,[indices{:}]) = [];
cWorker = cWorker-length(indices);


if(cWorker > 1)
   cStates = [workers{1,:}];
   cStates = {cStates.State};
end

end

% One last check for job output
for i = 1:cWorker-1
    if strcmp(workers{1,i}.State,'finished')
       fprintf(workers{2,i});
    end
end

end