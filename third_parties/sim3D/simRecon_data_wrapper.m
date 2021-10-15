function [] = simRecon_data_wrapper(dataPaths, PSFs, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addRequired('PSFs');

ip.addParameter('ChannelPatterns',{'CamA','CamB'},@iscell);

ip.addParameter('islattice', true, @islogical); %Flag to indicate if this is light sheet data
ip.addParameter('NA_det', 1.0, @isnumeric);
ip.addParameter('NA_ext', 0.55, @isnumeric);
ip.addParameter('nimm', 1.33, @isnumeric);
ip.addParameter('wvl_em', .605, @isnumeric);
ip.addParameter('wvl_ext', .560, @isnumeric);
ip.addParameter('w', 5e-3, @isnumeric); %Wiener coefficient for regularization
ip.addParameter('apodize', true, @islogical); %Flag to indicate whether or not to apodize the final data

ip.addParameter('DS', false, @islogical);
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
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 2, @isnumeric);

ip.parse(dataPaths, PSFs, varargin{:});

if ischar(dataPaths)
    dataPaths = {dataPaths};
end

if ischar(PSFs)
    PSFs = {PSFs};
end

pr = ip.Results;

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
            dataReconFullpath = [dataPaths{i} filesep 'sim_recon' filesep fsname '_recon' ext];
            inputFullpaths{j} = dataFullpath;
            outputFullpaths{j} = dataReconFullpath;

            funcStrs{j} =  sprintf(['simReconFrame(''%s'',''%s'',''lattice_period'',%.10f,''phase_step'',%.10f,''norders''', ...
                ',%.10f,''nphases'',%.10f,''Overlap'',%.10f,''ChunkSize'',[%.10f,%.10f,%.10f],''edgeTaper'',%s,''edgeTaperVal'',%.10f,',...
                '''perdecomp'',%s,''useGPU'',%s,''DS'',%s,''Background'',%.10f)'], dataFullpath, PSFs{cPatt}, ...
                lattice_period, phase_step, norders, nphases, OL, ChunkSize, string(edgeTaper), edgeTaperVal, string(perdecomp), ...
                string(useGPU),  string(DS), Background);
        end
            
    if useGPU
        maxJobNum = inf;
        cpusPerTask = 5;
        cpuOnlyNodes = false;
        taskBatchNum = 5;
        SlurmParam = '-p abc --qos abc_normal -n1 --mem=167G --gres=gpu:1';
    else
        maxJobNum = inf;
        cpusPerTask = 24;
        cpuOnlyNodes = true;
        taskBatchNum = 1;
        SlurmParam = '-p abc --qos abc_normal -n1 --mem-per-cpu=21418M';
    end
    
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