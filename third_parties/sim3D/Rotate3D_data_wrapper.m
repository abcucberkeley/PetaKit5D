function [] = Rotate3D_data_wrapper(dataPaths, xyPixelSize, dz, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));

ip.addRequired('xyPixelSize'); % typical value: 0.1
ip.addRequired('dz'); % typical value: 0.2-0.5
ip.addOptional('SkewAngle', 32.45, @isscalar);
ip.addOptional('Reverse', true, @islogical);

ip.addParameter('ChannelPatterns',{'CamA','CamB'},@iscell);

ip.addParameter('ObjectiveScan', false, @islogical);
ip.addParameter('sCMOSCameraFlip', false, @islogical);
ip.addParameter('Save16bit', false , @islogical); % saves deskewed data as 16 bit -- not for quantification

ip.addParameter('Overwrite', false , @islogical);
ip.addParameter('flipZstack', false, @islogical);
ip.addParameter('EdgeSoften', 5, @isnumeric); % # of xy px to soften
ip.addParameter('zEdgeSoften', 2, @isnumeric); % # of xy px to soften
ip.addParameter('Crop', [], @isnumeric); % requires lower and higher values for cropping
ip.addParameter('maxSubVolume', 3e8, @isnumeric);
ip.addParameter('CPUMaxMem', 500, @isnumeric); % CPU Memory in Gb
ip.addParameter('largeFile', false, @islogical);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', false, @islogical); % master node participate in the task computing.
ip.addParameter('masterCPU', false, @islogical); % master node is a cpu node, which is just for large file deconvolution.
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 5, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 2, @isnumeric);

ip.parse(dataPaths, xyPixelSize, dz, varargin{:});

if ischar(dataPaths)
    dataPaths = {dataPaths};
end

pr = ip.Results;

SkewAngle = pr.SkewAngle;
Reverse = pr.Reverse;

ChannelPatterns = pr.ChannelPatterns;



Save16bit = pr.Save16bit;

ObjectiveScan = pr.ObjectiveScan;

largeFile = pr.largeFile;
parseCluster = pr.parseCluster;
jobLogDir = pr.jobLogDir;
masterCompute = pr.masterCompute;
maxTrialNum = pr.maxTrialNum;
uuid = pr.uuid;
if isempty(uuid)
    uuid = get_uuid();
end


for i = 1:numel(dataPaths)
    for cPatt = 1:numel(ChannelPatterns)
        fnames = dir([dataPaths{i} filesep '*' ChannelPatterns{cPatt} '*.tif']);
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
        for j = 1: numel(fnames)
            [pathstr, fsname, ext] = fileparts(fnames{j});
            dataFullpath = [dataPaths{i} filesep fnames{j}];
            dataDSFullpath = [dataPaths{i} filesep 'Rotated' filesep fsname ext];
            inputFullpaths{j} = dataFullpath;
            outputFullpaths{j} = dataDSFullpath;
            
            funcStrs{j} =  sprintf(['XR_RotateFrame3D(''%s'', %.10f, %.10f,''SkewAngle'',%.10f,''Reverse'',%s,''ObjectiveScan''', ...
                ',%s,''Save16bit'',%s)'], dataFullpath, xyPixelSize, dz, SkewAngle, string(Reverse), string(ObjectiveScan),string(Save16bit));
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
        
        is_done_flag= slurm_cluster_generic_computing_wrapper(inputFullpaths, outputFullpaths, ...
            funcStrs, 'cpusPerTask', cpusPerTask, 'cpuOnlyNodes', cpuOnlyNodes, 'SlurmParam', SlurmParam, ...
            'maxJobNum', maxJobNum, 'taskBatchNum', taskBatchNum, 'masterCompute', masterCompute, 'parseCluster', parseCluster);
        
        if ~all(is_done_flag)
            slurm_cluster_generic_computing_wrapper(inputFullpaths, outputFullpaths, ...
                funcStrs, 'cpusPerTask', cpusPerTask, 'cpuOnlyNodes', cpuOnlyNodes, 'SlurmParam', SlurmParam, ...
                'maxJobNum', maxJobNum, 'taskBatchNum', taskBatchNum, 'masterCompute', masterCompute, 'parseCluster', parseCluster);
        end
    end
    
    
end

end

