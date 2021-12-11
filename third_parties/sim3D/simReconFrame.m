function [im] = simReconFrame(frameFullpaths, otf, varargin)

ip = inputParser;
ip.CaseSensitive = false;


ip.addRequired('frameFullpaths', @(x) ischar(x) || iscell(x));


ip.addRequired('otf');

ip.addParameter('resultsDirName', 'sim_recon', @ischar);

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
ip.addParameter('gpuPrecision', 'single', @ischar);

ip.addParameter('flipZstack', false, @islogical);
ip.addParameter('EdgeSoften', 5, @isnumeric); % # ofxy px to soften
ip.addParameter('zEdgeSoften', 2, @isnumeric); % # ofxy px to soften
ip.addParameter('Crop', [], @isnumeric); % requires lower and higher values for cropping
ip.addParameter('zFlip', false, @islogical);
ip.addParameter('GenMaxZproj', [0,0,1] , @isnumeric);
ip.addParameter('ResizeImages', [] , @isnumeric);
ip.addParameter('EdgeErosion', 0 , @isnumeric); % erode edges for certain size.
ip.addParameter('ErodeBefore', false, @islogical);
ip.addParameter('ErodeAfter', true,@islogical);
ip.addParameter('ErodeMaskfile', '', @ischar); % erode edges file
ip.addParameter('SaveMaskfile', false, @islogical); % save mask file for common eroded mask
ip.addParameter('ChunkSize', [250, 250, 250] , @isvector); % in y, x, z
ip.addParameter('Overlap', 10, @isnumeric); % block overlap
ip.addParameter('maxSubVolume', 3e8, @isnumeric);
ip.addParameter('CPUMaxMem', 500, @isnumeric); % CPU Memory in Gb
ip.addParameter('largeFile', false, @islogical);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing.
ip.addParameter('masterCPU', false, @islogical); % master node is a cpu node, which is just for large file deconvolution.
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 5, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 2, @isnumeric);
ip.addParameter('intThresh', 1, @isnumeric);
ip.addParameter('occThresh', 0.8, @isnumeric);

ip.parse(frameFullpaths, otf, varargin{:});

if ischar(frameFullpaths)
    frameFullpaths = {frameFullpaths};
end

pr = ip.Results;

resultsDirName = pr.resultsDirName;

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

Save16bit = pr.Save16bit;
gpuPrecision = pr.gpuPrecision;

useGPU = pr.useGPU;

useGPU = useGPU & gpuDeviceCount > 0;



if isempty(otf)
    error('You should provide an otf file for the frame...\n');
end

if ischar(otf)
    otf = sim_PSFtoOTF_gen(otf,'nphases',nphases,'norders',norders,'norientations',norientations, 'lattice_period',lattice_period, ...
        'phase_step',phase_step,'pxl_dim_PSF',pxl_dim_PSF,'Background',Background, 'useGPU', useGPU, 'Save16bit', Save16bit, 'gpuPrecision', gpuPrecision);
end

% parameters

flipZstack = pr.flipZstack;

EdgeErosion = pr.EdgeErosion;
ErodeBefore = pr.ErodeBefore;
ErodeAfter = pr.ErodeAfter;
ErodeMaskfile = pr.ErodeMaskfile;
if ~isempty(ErodeMaskfile) && exist(ErodeMaskfile, 'file')
    EdgeErosion = 0; % set EdgeErosion length as 0.
end
SaveMaskfile = pr.SaveMaskfile;

OL = pr.Overlap;
ChunkSize = pr.ChunkSize;
maxSubVolume = pr.maxSubVolume;
%CPUMaxMem = pr.CPUMaxMem;
%parseCluster = pr.parseCluster;
%jobLogDir = pr.jobLogDir;
%masterCompute = pr.masterCompute;
%maxTrialNum = pr.maxTrialNum;
uuid = pr.uuid;
if isempty(uuid)
    uuid = get_uuid();
end


% if master node is cpu node, msterCompute must be false, it is only for
% large file cluster computing
% if masterCPU
%     masterCompute = false;
%     largeFile = true;
%     parseCluster = true;
% end

%if parseCluster
%    [status, ~] = system('sinfo');
%    if status ~= 0
%        warning('A slurm-based computing cluster is not exist. Set parseCluster as false.')
%        parseCluster = false;
%    end
%end

nF = numel(frameFullpaths);
for f = 1 : nF
    frameFullpath = frameFullpaths{f};
    if ~exist(frameFullpath, 'file')
        warning('%s does not exist! Skip it!', frameFullpath);
        continue;
    end
    [pathstr, fsname, ext] = fileparts(frameFullpath);
    reconPath = '';
    if isempty(reconPath)
        reconPath = [pathstr '/' resultsDirName '/'];
        mkdir(reconPath);
        fileattrib(reconPath, '+w', 'g');
    end
    
    
    % check file size
    %[estMem] = XR_estimateComputingMemory(frameFullpath, {'deconvolution'}, 'cudaDecon', false);
    
    % The memory size is 500G
    %if estMem > CPUMaxMem
    %    largeFile = true;
    %end
    
    
    
    % first try single GPU deconvolution, if it fails split into multiple chunks
    reconFullPath = [reconPath '/' fsname '_recon.tif'];
    if exist(reconFullPath, 'file') && ~pr.Overwrite
        disp('Reconstruction results already exist, skip it!');
        continue;
    end
    switch ext
        case {'.tif', '.tiff'}
            im_raw = readtiff(frameFullpath);
        case '.zarr'
            bim = blockedImage(frameFullpath, 'Adapter', ZarrAdapter);
            im_raw = gather(bim);
    end
    
    if flipZstack
        im_raw = flip(im_raw, 3);
    end
    
    if EdgeErosion > 0
        if ~isempty(ErodeMaskfile) && exist(ErodeMaskfile, 'file')
            fprintf('Eroding edges for data...\n');
            im_bw_erode = readtiff(ErodeMaskfile);
            if ErodeBefore
                im_raw = im_raw .* cast(im_bw_erode, class(im_raw));
            end
        else                
            fprintf('Create eroded masks using raw data...\n');
            im_bw = im_raw > 0;
            % pad to avoid not erosion if a pixel touching the boundary
            im_bw_pad = false(size(im_bw) + 2);
            im_bw_pad(2 : end - 1, 2 : end - 1, 2 : end - 1) = im_bw;
            im_bw_erode = imerode(im_bw_pad, strel('disk', EdgeErosion));
            im_bw_erode = im_bw_erode(2 : end - 1, 2 : end - 1, 2 : end - 1);
            
            fprintf('Eroding edges for data...\n');
            if ErodeBefore
                im_raw = im_raw .* cast(im_bw_erode, class(im_raw));
            end
            
            clear im_bw im_bw_pad
            
            % save mask file as common one for other time points/channels
            if SaveMaskfile
                maskPath = [reconPath, '/', 'Masks'];
                if ~exist(maskPath, 'dir')
                    mkdir(maskPath);
                end
                maskFullPath = sprintf('%s/%s_eroded.tif', maskPath, fsname);
                maskTmpPath = sprintf('%s/%s_eroded_%s.tif', maskPath, fsname, uuid);
                writetiff(uint8(im_bw_erode), maskTmpPath);
                movefile(maskTmpPath, maskFullPath);
            end
        end                  
        if ErodeAfter
           im_bw_erode =  imresize3(im_bw_erode,size(im_bw_erode)./[.5,.5,nphases*norientations],'nearest'); 
        end
        %clear im_bw_erode;
    end
    
    
    reconTmpPath = sprintf('%s_%s_recon.tif', reconFullPath(1:end-10), uuid);
    imSize = size(im_raw);
    if Save16bit
        dtype = 'uint16';
    else
        dtype = 'single';
    end
    
    % if the image is skewed, don't chunk in x
    if~DS && useGPU
        ChunkSize(1) = round(ChunkSize(1) / imSize(2) * ChunkSize(2));
        ChunkSize(2) = imSize(2);
    end
    
    [xmin,xmax,ymin,ymax,zmin,zmax,nn] = XR_subVolumeCoordinatesExtraction(imSize ./ [1, 1, nphases*norientations], 'ChunkSize', ChunkSize, 'overlapSize', OL, 'maxSubVolume', maxSubVolume);
    
    % make sure the chunk size in x and y are even numbers, if not pad one
    % pixel in the left side (or right side if left side is the border)
    if any(rem(xmax - xmin + 1, 2) == 1)
        xmin = xmin - 1;
        xmax = xmax + (xmin == 0);
        xmin = xmin + (xmin == 0);
    end
    
    if any(rem(ymax - ymin + 1, 2) == 1)
        ymin = ymin - 1;
        ymax = ymax + (ymin == 0);
        ymin = ymin + (ymin == 0);
    end
    
    % create a folder for the file and write out the chunks
    fprintf('Processing image chunks...\n')
    im = zeros([imSize(1)*2,imSize(2)*2,imSize(3)/(nphases*norientations)], dtype);
    
    imSize = size(im);
    lol = floor(OL / 2);
    rol = ceil(OL / 2);
    chunkTimer = tic;
    for ck = 1:nn
        zmin_ck = (zmin(ck) - 1)* nphases*norientations + 1;
        zmax_ck = zmax(ck) * nphases*norientations;
        if(useGPU)
            im_chunk = gpuArray(im_raw(ymin(ck):ymax(ck), xmin(ck):xmax(ck),  zmin_ck : zmax_ck));
        else
            im_chunk = im_raw(ymin(ck):ymax(ck), xmin(ck):xmax(ck),  zmin_ck : zmax_ck);
        end
        
        % for blank region, just skip it
        if sum(im_chunk(:)) > intThresh && nnz(im_chunk(:))/ numel(im_chunk)>= occThresh
            fprintf('processing chunk:%d of %d \n',ck, nn)
            try
                im_chunk = simRecon(im_chunk, otf, 'islattice', islattice, 'NA_det', NA_det, 'NA_ext', NA_ext, 'nimm', nimm, ...
                    'wvl_em', wvl_em, 'wvl_ext', wvl_ext, 'w', w, 'apodize', apodize, 'nphases', nphases, 'norders', norders, ...
                    'norientations', norientations, 'lattice_period', lattice_period, 'lattice_angle', lattice_angle, 'phase_step', phase_step, ...
                    'pxl_dim_data', pxl_dim_data, 'pxl_dim_PSF', pxl_dim_PSF, 'Background', Background, 'useGPU', useGPU, ...
                    'normalize_orientations', normalize_orientations, 'perdecomp', perdecomp, 'edgeTaper', edgeTaper, 'edgeTaperVal', edgeTaperVal,'Save16bit',Save16bit,'gpuPrecision',gpuPrecision);
                
                
                tim = im_chunk;
                
                [tsy, tsx, tsz] = size(tim);
                
                
                ymin_ck = ymin(ck) * 2 - 1;
                xmin_ck = xmin(ck) * 2 - 1;
                
                
                yrange = ymin_ck + (ymin_ck ~= 1) * lol * 2 : ymax(ck) * 2 - (ymax(ck) * 2 ~= imSize(1)) * rol * 2;
                xrange = xmin_ck + (xmin_ck ~= 1) * lol * 2 : xmax(ck) * 2 - (xmax(ck) * 2 ~= imSize(2)) * rol * 2;
                zrange = zmin(ck) + (zmin(ck) ~= 1) * lol : zmax(ck) - (zmax(ck) ~= imSize(3)) * rol;
                yrange = yrange(yrange - ymin_ck + 1 <= tsy);
                xrange = xrange(xrange - xmin_ck + 1 <= tsx);
                zrange = zrange(zrange - zmin(ck) + 1 <= tsz);
                tyrange = yrange - ymin_ck + 1;
                txrange = xrange - xmin_ck + 1;
                tzrange = zrange - zmin(ck) + 1;
                im(yrange, xrange, zrange) = tim(tyrange, txrange, tzrange);
                clear tim;
            catch
                fprintf('chunk %d failed!!!!!!!!!!!!\n',ck)
            end
        end
        
    end
    toc(chunkTimer);
    %{
    if EdgeErosion > 0
        fprintf('Erode edges of deconvolved data w.r.t. raw data...\n');
        im = im .* cast(im_bw_erode, class(im));
    end
    
    if ~isempty(ErodeMaskfile) && exist(ErodeMaskfile, 'file')
        fprintf('Erode edges of deconvolved data using a predefined mask...\n');
        im_bw_erode = readtiff(ErodeMaskfile);
        im = im .* cast(im_bw_erode, class(im));
    end
    
    %}
    if EdgeErosion > 0 && ErodeAfter
       im = im .* cast(im_bw_erode, class(im));
    end
    
    reconTmpPath_eroded = sprintf('%s_%s_eroded.tif', reconFullPath(1:end-4), uuid);
    
    if Save16bit
        im = uint16(im);
        writetiff(im.*uint16((im>=0)), reconTmpPath_eroded);
    else
        im = single(im);
        writetiff(im.*(im>=0), reconTmpPath_eroded);
    end
    
    movefile(reconTmpPath_eroded, reconFullPath);
    delete(reconTmpPath);
    reconTmpMIPPath = sprintf('%s/%s_%s_MIP_z.tif', reconPath, fsname, uuid);
    delete(reconTmpMIPPath);
    
    reconMIPPath = sprintf('%s/MIPs/', reconPath);
    if ~exist(reconMIPPath, 'dir')
        mkdir(reconMIPPath);
        fileattrib(reconMIPPath, '+w', 'g');
    end
    reconMIPFullPath = sprintf('%s%s_MIP_z.tif', reconMIPPath, fsname);
    writetiff(max(im,[],3), reconMIPFullPath);
    
    if exist(reconFullPath, 'file')
        fprintf('Completed sim reconstruction of %s.\n', frameFullpath);
        continue;
    end
end

end
