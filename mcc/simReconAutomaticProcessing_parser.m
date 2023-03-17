function [] = simReconAutomaticProcessing_parser(dataPaths, varargin)

%#function deskewPhasesFrame
%#function simReconFrame

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));


ip.addParameter('PSFs',{''},@(x) ischar(x) || iscell(x));

ip.addParameter('Deskew',true,@(x) islogical(x) || ischar(x));
ip.addParameter('Recon',true,@(x) islogical(x) || ischar(x));
ip.addParameter('Streaming',false,@(x) islogical(x) || ischar(x));
ip.addParameter('resultsDirName', 'sim_recon', @ischar);

ip.addParameter('reconBatchNum', 5, @(x) isnumeric(x) || ischar(x));
ip.addParameter('parPoolSize', 24, @(x) isnumeric(x) || ischar(x));

ip.addParameter('xyPixelSize',.108,@(x) isnumeric(x) || ischar(x)); % typical value: 0.1
ip.addParameter('dz',.5,@(x) isnumeric(x) || ischar(x)); % typical value: 0.2-0.5
ip.addOptional('SkewAngle', 32.45, @(x) isscalar(x) || ischar(x));
ip.addOptional('Reverse', true, @(x) islogical(x) || ischar(x));

ip.addParameter('Rotate',false,@(x) islogical(x) || ischar(x)); % Rotate after deskew

ip.addParameter('ChannelPatterns',{'CamA','CamB'},@(x) iscell(x) || ischar(x));

ip.addParameter('islattice', true, @(x) islogical(x) || ischar(x)); %Flag to indicate if this is light sheet data
ip.addParameter('NA_det', 1.0, @(x) isnumeric(x) || ischar(x));
ip.addParameter('NA_ext', 0.55, @(x) isnumeric(x) || ischar(x));
ip.addParameter('nimm', 1.33, @(x) isnumeric(x) || ischar(x));
ip.addParameter('wvl_em', .605, @(x) isnumeric(x) || ischar(x));
ip.addParameter('wvl_ext', .560, @(x) isnumeric(x) || ischar(x));
ip.addParameter('w', 5e-3, @(x) isnumeric(x) || ischar(x)); %Wiener coefficient for regularization
ip.addParameter('apodize', true, @(x) islogical(x) || ischar(x)); %Flag to indicate whether or not to apodize the final data

ip.addParameter('DS', true, @(x) islogical(x) || ischar(x));
ip.addParameter('nphases', 5, @(x) isnumeric(x) || ischar(x));
ip.addParameter('norders', 5, @(x) isnumeric(x) || ischar(x));
ip.addParameter('norientations', 1, @(x) isnumeric(x) || ischar(x));
ip.addParameter('lattice_period', 1.2021, @(x) isnumeric(x) || ischar(x)); %Lattice period in microns - this is the coarsest period
ip.addParameter('lattice_angle', [pi/2], @(x) isnumeric(x) || ischar(x)); %Angle parellel to pattern modulation (assuming horizontal is zero)
ip.addParameter('phase_step', .232, @(x) isnumeric(x) || ischar(x)); %Phase step in microns
ip.addParameter('pxl_dim_data', [0.11,0.11,0.3*sind(32.4)], @(x) isnumeric(x) || ischar(x)); %Voxel dimensions of the image in microns - note, stored as [dy,dx,dz]
ip.addParameter('pxl_dim_PSF', [0.11,0.11,0.2*sind(32.4)], @(x) isnumeric(x) || ischar(x)); %Voxel dimensions of the PSF in microns - note, stored as [dy,dx,dz]
ip.addParameter('Background', 105, @(x) isnumeric(x) || ischar(x));

ip.addParameter('normalize_orientations', false, @(x) islogical(x) || ischar); %Flag to indicate whether or not to normalize total intensity for each orientation
ip.addParameter('perdecomp', false, @(x) islogical(x) || ischar(x)); %Flag to indicate whether or not to use periodic/smooth decomposition to reduce edge effects
ip.addParameter('edgeTaper', false, @(x) islogical(x) || ischar(x)); %Flag to indicate whether or not to window the data to reduce edge effects
ip.addParameter('edgeTaperVal', 0.1, @(x) isnumeric(x) || ischar(x)); %Roll-off parameter for Tukey windowing

ip.addParameter('useGPU', true, @(x) islogical(x) || ischar(x));

ip.addParameter('Overwrite', false , @(x) islogical(x) || ischar(x));
ip.addParameter('Save16bit', false , @(x) islogical(x) || ischar(x));
ip.addParameter('gpuPrecision', 'single', @ischar);

ip.addParameter('flipZstack', false, @(x) islogical(x) || ischar(x));
ip.addParameter('EdgeSoften', 5, @(x) isnumeric(x) || ischar(x)); % # ofxy px to soften
ip.addParameter('zEdgeSoften', 2, @(x) isnumeric(x) || ischar(x)); % # ofxy px to soften
ip.addParameter('Crop', [], @(x) isnumeric(x) || ischar(x)); % requires lower and higher values for cropping
ip.addParameter('zFlip', false, @(x) islogical(x) || ischar(x));
ip.addParameter('GenMaxZproj', [0,0,1] , @(x) isnumeric(x) || ischar(x));
ip.addParameter('ResizeImages', [] , @(x) isnumeric(x) || ischar(x));
ip.addParameter('EdgeErosion', 0 , @(x) isnumeric(x) || ischar(x)); % erode edges for certain size.
ip.addParameter('ErodeBefore', false, @(x) islogical(x) || ischar(x));
ip.addParameter('ErodeAfter', true,@(x) islogical(x) || ischar(x));
ip.addParameter('ErodeMaskfile', '', @ischar); % erode edges file
ip.addParameter('SaveMaskfile', false, @(x) islogical(x) || ischar(x)); % save mask file for common eroded mask
ip.addParameter('ChunkSize', [250, 250, 250] , @(x) isvector(x) || ischar(x)); % in y, x, z
ip.addParameter('Overlap', 10, @(x) isnumeric(x) || ischar(x)); % block overlap
ip.addParameter('maxSubVolume', 3e8, @(x) isnumeric(x) || ischar(x));
ip.addParameter('CPUMaxMem', 500, @(x) isnumeric(x) || ischar(x)); % CPU Memory in Gb

ip.addParameter('splitJobsByChannelPattern', false, @(x) islogical(x) || ischar(x));

ip.addParameter('largeFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('masterCompute', false, @(x) islogical(x) || ischar(x)); % master node participate in the task computing.
ip.addParameter('masterCPU', false, @(x) islogical(x) || ischar(x)); % master node is a cpu node, which is just for large file deconvolution.
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 5, @(x) isnumeric(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxModifyTime', 10, @(x) isnumeric(x) || ischar(x)); % the maximum during of last modify time of a file, in minute.
ip.addParameter('maxTrialNum', 3, @(x) isnumeric(x) || ischar(x));
ip.addParameter('unitWaitTime', 2, @(x) isnumeric(x) || ischar(x));
ip.addParameter('MatlabLaunchStr', 'module load matlab/r2021a; matlab -nodisplay -nosplash -nodesktop -nojvm -r', @ischar);
ip.addParameter('SlurmParam', '-p abc --qos abc_normal -n1 --mem-per-cpu=21418M', @ischar);
ip.addParameter('intThresh', 1, @(x) isnumeric(x) || ischar(x));
ip.addParameter('occThresh', 0.8, @(x) isnumeric(x) || ischar(x));

ip.parse(dataPaths, varargin{:});

pr = ip.Results;

PSFs = pr.PSFs;

Deskew = pr.Deskew;
Recon = pr.Recon;
Streaming = pr.Streaming;
resultsDirName = pr.resultsDirName;

reconBatchNum = pr.reconBatchNum;
parPoolSize = pr.parPoolSize;

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
intThresh = pr.intThresh;
occThresh = pr.occThresh;

useGPU = pr.useGPU;


Overwrite = pr.Overwrite;
flipZstack = pr.flipZstack;
Save16bit = pr.Save16bit;
gpuPrecision = pr.gpuPrecision;

EdgeErosion = pr.EdgeErosion;
ErodeBefore = pr.ErodeBefore;
ErodeAfter = pr.ErodeAfter;
ErodeMaskfile = pr.ErodeMaskfile;

SaveMaskfile = pr.SaveMaskfile;

Overlap = pr.Overlap;
ChunkSize = pr.ChunkSize;
maxSubVolume = pr.maxSubVolume;

splitJobsByChannelPattern = pr.splitJobsByChannelPattern;

%largeFile = pr.largeFile;
parseCluster = pr.parseCluster;
jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
masterCompute = pr.masterCompute;
maxModifyTime = pr.maxModifyTime;
maxTrialNum = pr.maxTrialNum;
unitWaitTime = pr.unitWaitTime;
MatlabLaunchStr = pr.MatlabLaunchStr;
SlurmParam = pr.SlurmParam;
uuid = pr.uuid;

if ischar(PSFs)
    PSFs = eval(PSFs);
end
if ischar(Deskew)
    Deskew = strcmp(Deskew,'true');
end
if ischar(Recon)
    Recon = strcmp(Recon,'true');
end
if ischar(Streaming)
    Streaming = strcmp(Streaming,'true');
end
if ischar(reconBatchNum)
    reconBatchNum = str2num(reconBatchNum);
end
if ischar(parPoolSize)
    parPoolSize = str2num(parPoolSize);
end
if ischar(xyPixelSize)
    xyPixelSize = str2num(xyPixelSize);
end
if ischar(dz)
    dz = str2num(dz);
end
if ischar(SkewAngle)
    SkewAngle = str2num(SkewAngle);
end
if ischar(Reverse)
    Reverse = strcmp(Reverse,'true');
end
if ischar(Rotate)
    Rotate = strcmp(Rotate,'true');
end
if ischar(ChannelPatterns)
    ChannelPatterns = eval(ChannelPatterns);
end
if ischar(islattice)
    islattice = strcmp(islattice,'true');
end
if ischar(NA_det)
    NA_det = str2num(NA_det);
end
if ischar(NA_ext)
    NA_ext = str2num(NA_ext);
end
if ischar(nimm)
    nimm = str2num(nimm);
end
if ischar(wvl_em)
    wvl_em = str2num(wvl_em);
end
if ischar(wvl_ext)
    wvl_ext = str2num(wvl_ext);
end
if ischar(w)
    w = str2num(w);
end
if ischar(apodize)
    apodize = strcmp(apodize,'true');
end
if ischar(DS)
    DS = strcmp(DS,'true');
end
if ischar(nphases)
    nphases = str2num(nphases);
end
if ischar(norders)
    norders = str2num(norders);
end
if ischar(norientations)
    norientations = str2num(norientations);
end
if ischar(lattice_period)
    lattice_period = str2num(lattice_period);
end
if ischar(lattice_angle)
    lattice_angle = str2num(lattice_angle);
end
if ischar(phase_step)
    phase_step = str2num(phase_step);
end
if ischar(pxl_dim_data)
    pxl_dim_data = str2num(pxl_dim_data);
end
if ischar(pxl_dim_PSF)
    pxl_dim_PSF = str2num(pxl_dim_PSF);
end
if ischar(Background)
    Background = str2num(Background);
end
if ischar(normalize_orientations)
    normalize_orientations = strcmp(normalize_orientations,'true');
end
if ischar(perdecomp)
    perdecomp = strcmp(perdecomp,'true');
end
if ischar(edgeTaper)
    edgeTaper = strcmp(edgeTaper,'true');
end
if ischar(edgeTaperVal)
    edgeTaperVal = str2num(edgeTaperVal);
end
if ischar(useGPU)
    useGPU = strcmp(useGPU,'true');
end
if ischar(Overwrite)
    Overwrite = strcmp(Overwrite,'true');
end
if ischar(Save16bit)
    Save16bit = strcmp(Save16bit,'true');
end
if ischar(flipZstack)
    flipZstack = strcmp(flipZstack,'true');
end
if ischar(EdgeErosion)
    EdgeErosion = str2num(EdgeErosion);
end
if ischar(ErodeAfter)
    ErodeAfter = strcmp(ErodeAfter,'true');
end
if ischar(SaveMaskfile)
    SaveMaskfile = strcmp(SaveMaskfile,'true');
end
if ischar(ChunkSize)
    ChunkSize = str2num(ChunkSize);
end
if ischar(Overlap)
    Overlap = str2num(Overlap);
end
if ischar(maxSubVolume)
    maxSubVolume = str2num(maxSubVolume);
end
if ischar(splitJobsByChannelPattern)
    splitJobsByChannelPattern = strcmp(splitJobsByChannelPattern,'true');
end
if ischar(parseCluster)
    parseCluster = strcmp(parseCluster,'true');
end
if ischar(masterCompute)
    masterCompute = strcmp(masterCompute,'true');
end
if ischar(cpusPerTask)
    cpusPerTask = str2num(cpusPerTask);
end
if ischar(maxModifyTime)
    maxModifyTime = str2num(maxModifyTime);
end
if ischar(maxTrialNum)
    maxTrialNum = str2num(maxTrialNum);
end
if ischar(unitWaitTime)
    unitWaitTime = str2num(unitWaitTime);
end
if ischar(intThresh)
    intThresh = str2num(intThresh);
end
if ischar(occThresh)
    occThresh = str2num(occThresh);
end




simReconAutomaticProcessing(dataPaths,'PSFs',PSFs,'Deskew',Deskew,'Recon',Recon,...
    'Streaming',Streaming,'resultsDirName',resultsDirName,'reconBatchNum',reconBatchNum,...
    'parPoolSize',parPoolSize,'xyPixelSize',xyPixelSize,'dz',dz,...
    'SkewAngle',SkewAngle,'Reverse',Reverse,'Rotate',Rotate,'ChannelPatterns',ChannelPatterns,...
    'islattice',islattice,'NA_det',NA_det,'NA_ext',NA_ext,'nimm',nimm,'wvl_em',wvl_em,...
    'wvl_ext',wvl_ext,'w',w,'apodize',apodize,'DS',DS,'nphases',nphases,...
    'norders',norders,'norientations',norientations,'lattice_period',lattice_period,...
    'lattice_angle',lattice_angle,'phase_step',phase_step,...
    'pxl_dim_data',pxl_dim_data,'pxl_dim_PSF',pxl_dim_PSF,'Background',Background,...
    'normalize_orientations',normalize_orientations,'perdecomp',perdecomp,...
    'edgeTaper',edgeTaper,'edgeTaperVal',edgeTaperVal,'useGPU',useGPU,...
    'Overwrite',Overwrite,'Save16bit',Save16bit,'gpuPrecision',gpuPrecision,...
    'flipZstack',flipZstack,'EdgeErosion',EdgeErosion,'ErodeBefore',ErodeBefore,...
    'ErodeAfter',ErodeAfter,'ErodeMaskfile',ErodeMaskfile,...
    'SaveMaskfile',SaveMaskfile,'ChunkSize',ChunkSize,'Overlap',Overlap,...
    'maxSubVolume',maxSubVolume,'splitJobsByChannelPattern',splitJobsByChannelPattern,...
    'parseCluster',parseCluster,'masterCompute',masterCompute,'jobLogDir',jobLogDir,...
    'cpusPerTask',cpusPerTask,'uuid',uuid,'maxModifyTime',maxModifyTime,...
    'maxTrialNum',maxTrialNum,'unitWaitTime',unitWaitTime,...
    'MatlabLaunchStr',MatlabLaunchStr,'SlurmParam',SlurmParam,...
    'intThresh',intThresh,'occThresh',occThresh);


end