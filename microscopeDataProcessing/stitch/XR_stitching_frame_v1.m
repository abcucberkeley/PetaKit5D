function XR_stitching_frame_v1(tileFullpaths, coordinates, varargin)
% , px, SkewAngle, dz, xyPixelSize, reversed, Reverse, Crop, dfx, dfy, dfz, Save16bit, EdgeArtifacts ,PSFpath, Background, nIter, dzPSF, Deconvolve)
% function XR_stitching_frame(rt, sfnl, st, px, SkewAngle, dz, xyPixelSize, reversed, Reverse, Crop, dfx, dfy, dfz, Save16bit, EdgeArtifacts ,PSFpath, Background, nIter, dzPSF, Deconvolve)
%
% assumes all tiles are the same size
% based on GU_stitching_frame.m and java_stitching_frame_wrapper.m
%
% Author: Xiongtao Ruan (03/20/2020)
% 
% implement xcorr-based stitching
% xruan (04/2020): add bound box crop of the stitched result
% xruan (05/06/2020): add option for enable/disable xcorr-based stitching
% xruan (05/08/2020): add 2 (instead of 1) as fixed padding to avoid outside of nv.
% xruan (05/11/2020): change halfOrder to [1, 2, 3] for all blend methods.
% xruan (06/19/2020): add option for support of primary/all xcorr shift modes.
%                     include 'Abs' in filename
% xruan (06/25/2020): change [xyz]ridx upper bound calculation (not add 1) to avoid out of boundary.
% xruan (07/03/2020): also save final stitched size in primary/primaryFirst
%                     options, to make sure the sizes for different channels/time points are the same.
% xruan (07/05/2020): test different resampling resolution
% xruan (07/10/2020): set default resampling method as (1, 1, r_z / r_x),
%                     and save the output voxel factor
% xruan (07/13/2020): add option for DSR dir, and also indirectly write result with uuid.
% xruan (07/22/2020): change to save stitching information for all primary
% channels (whether xcorrshift true or not). Force to use Scan_ in the beginning.
% xruan (07/26/2020): change psz to pImSz to avoid name conflict. Just use
% bbox if it is set. If it is not set, use primary channel image size as
% standard size for all time points/channels. Inf in bbox(4:6) mean to
% image size, a value larger than imSz (not Inf) means pad to that size. 
% xruan (07/31/2020): add MIP for stitching. Fix potential blank regions
% (especially for tall tiles). 
% xruan (08/03/2020): add support for user-defined axis orders.
% xruan (08/08/2020): not load all tiles in the begining to save time and space
% xruan (08/12/2020): add support for tiles with variable size.
% xruan (08/14/2020): add support for primary channel for noxcorr shifts
% (keep size across channels the same). Add support for distributed
% computing for xcorr. Change downsample in xcorr as [2, 2, 2]
% xruan (08/17/2020): add support for stitching of DSR decon (only for
% existing DSR decon files).
% xruan (08/20/2020): add support for objective scan
% xruan (08/23/2020): add option for overlap type (full)
% xruan (09/12/2020): add saving of pixel size in the result folder
% xruan (09/13/2020): cast dsr to specific dtype if not the same as requirement


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('tileFullpaths', @iscell);
ip.addRequired('coordinates', @isnumeric);
ip.addParameter('ResultDir', 'matlab_stitch', @ischar);
ip.addParameter('stitchInfoDir', 'stitchInfo', @ischar);
ip.addParameter('stitchInfoFullpath', '', @ischar); % filename that contain stitch information for secondrary channels
ip.addParameter('DSRDirstr', '', @ischar); % path for DSRDirstr, if it is not true
ip.addParameter('DSRDeconDirstr', '', @ischar); % path for DSRDir decon str, if it is not true
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('px', 0.108, @isscalar);
ip.addParameter('SkewAngle', 32.45, @isscalar);
ip.addParameter('axisOrder', 'x,y,z', @ischar);
ip.addParameter('dz', 0.5, @isscalar);
ip.addParameter('xyPixelSize', 0.108, @isscalar);
ip.addParameter('ObjectiveScan', false, @islogical);
ip.addParameter('sCMOSCameraFlip', false, @islogical);
ip.addParameter('Reverse', false, @islogical);
ip.addParameter('Crop', false, @islogical);
ip.addParameter('df', [], @isnumeric);
ip.addParameter('Save16bit', false , @islogical); % saves deskewed data as 16 bit -- not for quantification
ip.addParameter('EdgeArtifacts', 2, @isnumeric);
ip.addParameter('Decon', false, @islogical);
ip.addParameter('cudaDecon', ~false, @islogical);
ip.addParameter('DS', false, @islogical);
ip.addParameter('DSR', false, @islogical);
ip.addParameter('resampleType', 'xy_isotropic', @isstr);
ip.addParameter('Background', [], @isnumeric);
ip.addParameter('deconRotate', false, @islogical);
ip.addParameter('PSFpath', '', @isstr);
ip.addParameter('dzPSF', 0.1, @isnumeric);
ip.addParameter('BlendMethod', 'mean', @isstr);
ip.addParameter('halfOrder', [1,2,3], @isnumeric);
ip.addParameter('overlapType', '', @isstr); % '', 'none', 'half', or 'full'
ip.addParameter('xcorrShift', true, @islogical);
ip.addParameter('isPrimaryCh', true, @islogical);
ip.addParameter('stitchPadSize', [4, 4, 2], @(x) isnumeric(x) && numel(x) == 3);
ip.addParameter('padSize', [], @(x) isnumeric(x) && (isempty(x) || numel(x) == 3));
ip.addParameter('boundboxCrop', [], @(x) isnumeric(x) && (isempty(x) || all(size(x) == [3, 2]) || numel(x) == 6));
ip.addParameter('zNormalize', false, @islogical);
ip.addParameter('xcorrDownsample', [2, 2, 2], @isnumeric);
ip.addParameter('SaveMIP', true , @islogical); % save MIP-z for stitch. 
ip.addParameter('tileNum', [] , @isnumeric); % tile number in each axis (TODO, not now, for consistency with previous analysis).
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @isstr);
ip.addParameter('cpusPerTask', 2, @isnumeric);
ip.addParameter('cpuOnlyNodes', true, @islogical);
ip.addParameter('uuid', '', @isstr);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 30, @isnumeric);

ip.parse(tileFullpaths, coordinates, varargin{:});

pr = ip.Results;
 % Resolution = pr.Resolution;
ResultDir = pr.ResultDir;
stitchInfoDir = pr.stitchInfoDir;
stitchInfoFullpath = pr.stitchInfoFullpath;
px = pr.px;
SkewAngle = pr.SkewAngle;
axisOrder = pr.axisOrder;
dz = pr.dz;
xyPixelSize = pr.xyPixelSize;
ObjectiveScan = pr.ObjectiveScan;
Reverse = pr.Reverse;
Crop = pr.Crop;
% Deskew = pr.Deskew;
% Rotate = pr.Rotate;
Decon = pr.Decon;
cudaDecon = pr.cudaDecon;
DS = pr.DS;
DSR = pr.DSR;
Background = pr.Background;
PSFpath = pr.PSFpath;
BlendMethod = pr.BlendMethod;
halfOrder = pr.halfOrder;
overlapType = pr.overlapType;
xcorrShift = pr.xcorrShift;
isPrimaryCh = pr.isPrimaryCh;
stitchPadSize = pr.stitchPadSize;
padSize = pr.padSize;
boundboxCrop = pr.boundboxCrop;
zNormalize = pr.zNormalize;
xcorrDownsample = pr.xcorrDownsample;
% RotateAfterDecon = pr.RotateAfterDecon;
% ChannelPatterns = pr.ChannelPatterns;
jobLogDir = pr.jobLogDir;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
cpuOnlyNodes = pr.cpuOnlyNodes;
maxTrialNum = pr.maxTrialNum;
unitWaitTime = pr.unitWaitTime;
Save16bit = pr.Save16bit;
EdgeArtifacts = pr.EdgeArtifacts;
SaveMIP = pr.SaveMIP;
uuid = pr.uuid;

[dataPath, ~] = fileparts(tileFullpaths{1});

% check if a slurm-based computing cluster exist
if parseCluster 
    [status, ~] = system('sinfo');
    if status ~= 0
        warning('A slurm-based computing cluster is not exist. Set parseCluster as false.')
        parseCluster = false;
    end
    if parseCluster && ~exist(jobLogDir, 'dir')
        warning('The job log directory does not exist, use ${stitching_tmp}/job_logs as job log directory')
        jobLogDir = sprintf('%s/job_logs', dataPath);
        if ~exist(jobLogDir, 'dir')
            mkdir(jobLogDir);
            fileattrib(jobLogDir, '+w', 'g');
        end
    end
    job_log_fname = [jobLogDir, '/job_%A_%a.out'];
    job_log_error_fname = [jobLogDir, '/job_%A_%a.err'];
    
    if cpuOnlyNodes
        slurm_constraint_str = ' --constraint=c24 ';
    else
        slurm_constraint_str = '';
    end
end

if isempty(uuid)
    uuid = get_uuid();
end

fprintf('Stitch for tiles:\n%s\n', strjoin(tileFullpaths, '\n'));

% check if axis order is valid
axisOrder = strrep(axisOrder, ' ', '');
pattern = '^(-?x,?-?y,?-?z|-?y,?-?x,?-?z|-?z,?-?y,?-?x|-?x,?-?z,?-?y|-?x,?-?z,?-?y|-?y,?-?z,?-?x)$';
if ~regexpi(axisOrder, pattern)
    error("The axisOrder is not right, it must has the form like 'y,x,z' or '-x,y,z' (flipped in x-axis)!");
end

order_sign_mat = zeros(2, 3); % first row: order for permute; second row: sign
axisOrder_pure = strrep(strrep(axisOrder, ',', ''), '-', '');
xyz_str = 'xyz';
for i = 1 : 3
    [~, order_sign_mat(1, :)] = sort(axisOrder_pure);
    order_sign_mat(2, i) = 1 - 2 * contains(axisOrder, ['-', xyz_str(i)]);
end
xyz = coordinates(:, order_sign_mat(1, :)) .* order_sign_mat(2, :);

if ObjectiveScan
    zAniso = dz/px;    
else
    zAniso = sind(SkewAngle)*dz/px;
end
theta = SkewAngle * pi/180;
% dx = cos(theta)*dz/xyPixelSize;
resample_type = pr.resampleType;
switch resample_type
    case 'method_1'  %% old one
        zf = cot(abs(theta));
        yf = 1;
        % xf = cos(abs(theta)) + tan(abs(theta))*sin(abs(theta));
        xf = 1 / cos(abs(theta));
    case 'isotropic'  %% isotropic
        zf = 1;
        yf = 1;
        xf = 1;
    case 'xy_isotropic' %% x, y isotropic and z: r_z / r_x
        zf = sqrt((sin(theta) ^ 2 + zAniso ^ 2 * cos(theta) ^ 2) / (cos(theta) ^ 2 + zAniso ^ 2 * sin(theta) ^ 2));
        yf = 1;
        xf = 1;
end

% create a text file to indicate pixel size
pixelInfoFname = sprintf('px%0.5g_py%0.5g_pz%0.5g', px*xf, px*yf, px*zf);
pixelInfoFullpath = sprintf('%s/%s/%s', dataPath, ResultDir, pixelInfoFname);
% check if there is some other pixel size file
dir_info = dir(sprintf('%s/%s/px*_py*_pz*', dataPath, ResultDir));
pixelFnames = {dir_info.name}';
for i = 1 : numel(pixelFnames)
    if ~strcmp(pixelFnames{i}, pixelInfoFname)
        delete([dir_info(i).folder, '/', pixelFnames{i}])
    end
end

if ~exist(pixelInfoFullpath, 'file')
    fclose(fopen(pixelInfoFullpath, 'w'));
end

% create DSR or Rotated folder
if ~isempty(pr.DSRDirstr)
    dsrpath = [dataPath '/' pr.DSRDirstr];
elseif ~isempty(pr.DSRDeconDirstr)
    dsrpath = [dataPath '/' pr.DSRDeconDirstr];
    Decon = true;
else
    if ObjectiveScan
        dsrpath = [dataPath '/' 'Rotated_dx' num2str(px*xf) '_dz' num2str(px*zf)];        
    else
        dsrpath = [dataPath '/' 'DSR_dx' num2str(px*xf) '_dz' num2str(px*zf)];
    end
end

% normalize coordinates to zero
% for i = 1:3
%     xyz(:,i) = xyz(:,i) - min(xyz(:,i));
% end
xyz = xyz - min(xyz, [], 1);
nF = numel(tileFullpaths);

if zNormalize
    zmed_cell = cell(nF, 1);
end

for k = 1:nF
    tic
    tileFullpath = tileFullpaths{k};
    [dataPath, fsname] = fileparts(tileFullpath);
    
    dsrFullpath = [dsrpath '/' fsname '.tif'];
    if Decon && ~isempty(pr.DSRDeconDirstr)
        dsrFullpath = [dsrpath '/' fsname '_decon.tif'];
    end
    
    if ~exist(dsrFullpath, 'file')
        im = readtiff(tileFullpath);
        
        if ~ObjectiveScan
            fprintf('deskewing Data...')
            tic
            im = deskewFrame3D(im, SkewAngle, dz, xyPixelSize, Reverse,...
                'Crop', Crop);
            toc
        end
        
        % deconvovle
        if Decon && isempty(pr.DSRDeconDirstr)
            deconpath = [rt '/' 'matlab_decon'];
            if ~exist(deconpath, 'dir')
                mkdir(deconpath);
            end
            
            if ~exist([deconpath '/' fsname '.tif'], 'file')
                fprintf('Deconvolving Data...')
                tic
                volpath = [rt fsname];
                mask = logical(im);
                im = RLdecon_for_ExMPipeline(im, volpath, PSFpath, Background, nIter, dzPSF, zAniso, 0, [], 0, ...
                    px, 0, 0, 0, 0, [0,0,1],0, []) ;toc
                
                % remove edge artifacts
                fprintf('Removing Edges...'); tic
                r = EdgeArtifacts;
                if r > 0
                    se = strel('disk', r);
                    mask = imerode(mask,se);
                    mask(:,:,1:r) = 0;
                    mask(:,:,end-r:end) = 0;
                    im(~mask) = 0;
                    clear mask
                    toc
                end

                if Save16bit
                    writetiff(uint16(im), [deconpath '/' fsname '_' uuid '.tif']);
                else
                    writetiff(im, [deconpath '/' fsname '_' uuid '.tif']);
                end
                movefile([deconpath '/' fsname '_' uuid '.tif'], [deconpath '/' fsname '.tif']);
            end
        end
        fprintf('Rotating Data...')
        tic
        dsr = rotateFrame3D(im, SkewAngle, zAniso, Reverse,...
            'Crop', true, 'ObjectiveScan', ObjectiveScan);
        toc
        
        fprintf('resampling Rotated Data...')
        if any([xf, yf, zf] ~= 1)
            tic
            dsr = GU_resampleStack3D(dsr, xf, yf, zf,'Interp', 'linear');toc
        end
        
        if ~exist(dsrpath, 'dir')
            mkdir(dsrpath);
        end
        
        fprintf('saving resampled + Rotated Data...')
        tic
        writetiff(uint16(dsr), dsrFullpath); toc
    else
        % only load data when zNormalize is true
        if zNormalize
            fprintf('loading resampling Rotated Data...')
            tic
            dsr = readtiff(dsrFullpath); toc
        end
    end
    
    if zNormalize
        zmed = zeros(size(dsr, 3), 1);
        for z = 1 : size(dsr, 3)
            dsr_z = dsr(:, :, z);
            if any(dsr_z ~= 0, 'all')
                zmed(z) = median(dsr_z(dsr_z ~= 0));
            end
        end
        zmed_cell{k} = zmed;
    end    
end

if zNormalize
    zmed_mat = cat(2, zmed_cell{:});
    overall_z_median = median(zmed_mat(zmed_mat ~= 0));
end

if ~exist('dsr', 'var')
    tileFullpath = tileFullpaths{1};
    [dataPath, fsname] = fileparts(tileFullpath);
    dsrFullpath = [dsrpath '/' fsname '.tif'];  
    if Decon && ~isempty(pr.DSRDeconDirstr)
        dsrFullpath = [dsrpath '/' fsname '_decon.tif'];
    end    
    dsr = readtiff(dsrFullpath);
end
dtype = class(dsr);
imSizes = zeros(nF, 3);
for i = 1 : nF
    tileFullpath = tileFullpaths{i};
    [dataPath, fsname] = fileparts(tileFullpath);
    dsrFullpath = [dsrpath '/' fsname '.tif'];  
    if Decon && ~isempty(pr.DSRDeconDirstr)
        dsrFullpath = [dsrpath '/' fsname '_decon.tif'];
    end    
    imSizes(i, :) = getImageSize(dsrFullpath);
end

if Save16bit
    dtype = 'uint16';
end 

% identify pairs between pairs
% first slight shift xyz such that the distance between pairs are in
% integer fold of resolution. 
xyz_orig = xyz;
xyz = round((xyz - min(xyz, [], 1)) ./ ([xf, yf, zf] * px)) .* ([xf, yf, zf] * px);

overlap_matrix = false(nF);
overlap_regions = zeros(nF * (nF - 1) / 2, 6);
for i = 1 : nF - 1
    for j = i + 1 : nF
        xyz_i = xyz(i, :);
        xyz_j = xyz(j, :);
        cuboid_i = [xyz_i; xyz_i + (imSizes(i, [2, 1, 3]) - 1) .* [xf, yf, zf] * px]';
        cuboid_j = [xyz_j; xyz_j + (imSizes(j, [2, 1, 3]) - 1) .* [xf, yf, zf] * px]';
        [is_overlap, cuboid_overlap] = cuboids_overlaps(cuboid_i, cuboid_j);
        if is_overlap 
            overlap_matrix(i, j) = true;
            % ind = 0.5 * i * (2 * j - i - 1);
            ind = (i - 1) * nF - i * (i + 1) / 2 + j;
            overlap_regions(ind, :) = cuboid_overlap(:);
        end
    end
end

% calculate relative shifts between tiles
relative_shift_mat = zeros(nF * (nF - 1) / 2, 3);
max_xcorr_mat = zeros(nF * (nF - 1) / 2, 3);
if xcorrShift && isPrimaryCh
    fprintf('Compute cross-correlation based registration between overlap tiles...\n')
    
    xcorrDir = [dataPath '/' ResultDir '/' 'xcorr' '/'];
    if ~exist(xcorrDir, 'dir')
        mkdir(xcorrDir);
        fileattrib(xcorrDir, '+w', 'g');
    end
    nPair = nF * (nF - 1) / 2;
    pair_ind_mat = tril(ones(nF), -1);
    pair_ind_mat(pair_ind_mat == 1) = 1 : nPair;
    pair_ind_mat = pair_ind_mat';
    
    is_done_flag = false(nPair, 1);
    trial_counter = zeros(nPair, 1);
    if parseCluster
        job_ids = -ones(nPair, 1);
    end

    while ~all(is_done_flag | trial_counter >= maxTrialNum, 'all')
        lastP = find(~is_done_flag & trial_counter < maxTrialNum, 1, 'last');

        for i = 1 : nF - 1
            tileFullpath_i = tileFullpaths{i};
            [dataPath, fsname_i] = fileparts(tileFullpath_i);
            dsrFullpath_i = [dsrpath '/' fsname_i '.tif'];
            if Decon && ~isempty(pr.DSRDeconDirstr)
                dsrFullpath_i = [dsrpath '/' fsname_i '_decon.tif'];
            end

            for j = i + 1 : nF
                task_id = pair_ind_mat(i, j);
                ind = (i - 1) * nF - i * (i + 1) / 2 + j;

                if is_done_flag(task_id) || trial_counter(task_id) >= maxTrialNum 
                    continue;
                end
                
                if ~overlap_matrix(i, j)
                    is_done_flag(task_id) = true;
                    continue;
                end
                
                xcorrFullpath = sprintf('%s/xcorr_tile_%d_tile_%d.mat', xcorrDir, i, j);
                tmpFullpath = sprintf('%s.tmp', xcorrFullpath(1 : end - 4));                
                if exist(xcorrFullpath, 'file')
                    if all(max_xcorr_mat(ind, :) == 0)
                        a = load(xcorrFullpath);
                        relative_shift_mat(ind, :) = a.relative_shift;
                        max_xcorr_mat(ind, :) = [i, j, a.max_xcorr];
                    end
                    
                    is_done_flag(task_id) = true;
                    if exist(tmpFullpath, 'file')
                        delete(tmpFullpath);
                    end
                    continue;
                end
                
                tileFullpath_j = tileFullpaths{j};
                [dataPath, fsname_j] = fileparts(tileFullpath_j);
                dsrFullpath_j = [dsrpath '/' fsname_j '.tif'];
                if Decon && ~isempty(pr.DSRDeconDirstr)
                    dsrFullpath_j = [dsrpath '/' fsname_j '_decon.tif'];
                end

                xyz_i = xyz(i, :);
                xyz_j = xyz(j, :);
                cuboid_i = [xyz_i; xyz_i + (imSizes(i, [2, 1, 3]) - 1) .* [xf, yf, zf] * px]';
                cuboid_j = [xyz_j; xyz_j + (imSizes(j, [2, 1, 3]) - 1) .* [xf, yf, zf] * px]';

                % ind = 0.5 * i * (2 * j - i - 1);
                cuboid_overlap_ij = reshape(overlap_regions(ind, :), 3, 2);

                % tic
                if exist(tmpFullpath, 'file') || parseCluster
                    if parseCluster
                        job_status = check_slurm_job_status(job_ids(task_id), rem(task_id, 1000));

                        % kill the first pending job and use master node do the computing.
                        if job_status == 0.5 && (masterCompute && task_id == lastP)
                            system(sprintf('scancel %d_%d', job_ids(task_id), task_id), '-echo');
                            trial_counter(task_id) = trial_counter(task_id) - 1;
                        end

                        % if the job is still running, skip it. 
                        if job_status == 1 
                            continue;
                        end

                        % If there is no job, submit a job
                        if job_status == -1 && ~(masterCompute && task_id == lastP)
                            [estMem, estGPUMem, rawImageSize] = XR_estimateComputingMemory(dsrFullpath_i, {'deconvolution'}, 'cudaDecon', false);
                            if cpusPerTask * 20 < rawImageSize * 5
                                cpusPerTask = min(24, ceil(rawImageSize * 5 / 20));
                            end

                            matlab_cmd = sprintf(['addpath(genpath(pwd));tic;cross_correlation_stitching_matching(''%s'',''%s'',''%s'',', ...
                                '[%s],[%s],[%s],%0.20d,[%s],''downSample'',[%s]);toc;'], dsrFullpath_i, dsrFullpath_j, xcorrFullpath, ...
                                strrep(mat2str(cuboid_i), ' ', ','), strrep(mat2str(cuboid_j), ' ', ','), strrep(mat2str(cuboid_overlap_ij), ' ', ','), ...
                                px, sprintf('%.20d,%.20d,%.20d', xf, yf, zf), strrep(num2str(xcorrDownsample, '%.20d,'), ' ', ''));
                            xcorr_cmd = sprintf('module load matlab/r2020a; matlab -nodisplay -nosplash -nodesktop -nojvm -r \\"%s\\"', matlab_cmd);
                            cmd = sprintf('sbatch --array=%d %s -o %s -e %s -p abc --qos abc_normal -n1 --mem-per-cpu=21418M --cpus-per-task=%d --wrap="%s"', ...
                                rem(task_id, 5000), slurm_constraint_str, job_log_fname, job_log_error_fname, cpusPerTask, xcorr_cmd);
                            [status, cmdout] = system(cmd, '-echo');

                            job_id = regexp(cmdout, 'Submitted batch job (\d+)\n', 'tokens');
                            job_id = str2double(job_id{1}{1});
                            job_ids(task_id) = job_id;
                            trial_counter(task_id) = trial_counter(task_id) + 1;                                
                        end
                    else
                        temp_file_info = dir(tmpFullpath);
                        if (datenum(clock) - [temp_file_info.datenum]) * 24 * 60 < unitWaitTime
                            continue; 
                        else
                            fclose(fopen(tmpFullpath, 'w'));
                        end
                    end
                else
                    fclose(fopen(tmpFullpath, 'w'));
                end

                if ~parseCluster || (parseCluster && masterCompute && task_id == lastP)
                    [relative_shift, max_xcorr] = cross_correlation_stitching_matching( ...
                        dsrFullpath_i, dsrFullpath_j, xcorrFullpath, cuboid_i, cuboid_j, ...
                        cuboid_overlap_ij, px, [xf, yf, zf]', 'downSample', xcorrDownsample);
                    trial_counter(task_id) = trial_counter(task_id) + 1;    
                end
                % toc
                if exist(xcorrFullpath, 'file')
                    if all(max_xcorr_mat(ind, :) == 0)
                        a = load(xcorrFullpath);
                        relative_shift_mat(ind, :) = a.relative_shift;
                        max_xcorr_mat(ind, :) = [i, j, a.max_xcorr];
                    end
                    is_done_flag(task_id) = true;
                    if exist(tmpFullpath, 'file')
                        delete(tmpFullpath);
                    end
                end
            end
        end
        if ~all(is_done_flag | trial_counter >= maxTrialNum, 'all') 
            pause(30);
        end
    end
    
    % only keep primary pairs by mst 
    if ~false && nF > 2
        max_xcorr_mat(max_xcorr_mat(:, 1) == 0 | max_xcorr_mat(:, 2) == 0, :) = [];
        G = graph(max_xcorr_mat(:, 1), max_xcorr_mat(:, 2), -max_xcorr_mat(:, 3));
        T = minspantree(G);
        aj = full(adjacency(T));
        overlap_matrix = overlap_matrix .* aj;
    end
elseif ~isPrimaryCh
    if ~exist(stitchInfoFullpath, 'file')
        error('The stitch information filename %s does not exist!', stitchInfoFullpaths);
    end
    
    a = load(stitchInfoFullpath, 'relative_shift_mat', 'overlap_matrix', 'pImSz');
    overlap_matrix = a.overlap_matrix;
    relative_shift_mat = a.relative_shift_mat;
    pImSz = a.pImSz;
end

% calculate absolute shifts of tiles
xyz_shift = xyz;

for i = 1 : nF - 1
    for j = i + 1 : nF
        if ~overlap_matrix(i, j)
            continue;
        end
        % ind = 0.5 * i * (2 * j - i - 1);
        ind = (i - 1) * nF - i * (i + 1) / 2 + j;

        xyz_shift(j, :) = (xyz_shift(i, :) - xyz(i, :)) + (xyz_shift(j, :) + relative_shift_mat(ind, :) .* [xf, yf, zf] * px);
    end
end

% use half of overlap region for pairs of tiles with overlap
if isempty(overlapType)
    if strcmp(BlendMethod, 'none')
        overlapType = 'none';
        % halfOrder = [1, 2, 3];
    else
        overlapType = 'half';
        % halfOrder = [2, 1, 3];    
        % halfOrder = [1, 2, 3];    
    end
end

half_ol_region_cell = cell(nF);
ol_region_cell = cell(nF);
for i = 1 : nF - 1
    for j = i + 1 : nF
        if ~overlap_matrix(i, j)
            continue;
        end
        
        xyz_i = xyz_shift(i, :);
        xyz_j = xyz_shift(j, :);
        cuboid_i = [xyz_i; xyz_i + (imSizes(i, [2, 1, 3]) - 1) .* [xf, yf, zf] * px]';
        cuboid_j = [xyz_j; xyz_j + (imSizes(j, [2, 1, 3]) - 1) .* [xf, yf, zf] * px]';

        [mregion_1, mregion_2] = compute_half_of_overlap_region(cuboid_i, cuboid_j, ...
            px, [xf, yf, zf]', 'overlapType', overlapType, 'halfOrder', halfOrder);
        half_ol_region_cell{i, j} = {mregion_1, mregion_2};
        
        [mregion_11, mregion_21] = compute_half_of_overlap_region(cuboid_i, cuboid_j, ...
            px, [xf, yf, zf]', 'overlapType', 'full', 'halfOrder', halfOrder);
        ol_region_cell{i, j} = {mregion_11, mregion_21};        
    end
end

% put tiles into the stitched image
% also pad the tiles such that the images of all time points and different
% channels are the same
if nF > 1 && xcorrShift
    pad_size = [xf, yf, zf] .* stitchPadSize * (px * nF);
else
    pad_size = [0, 0, 0];
end
origin = min(xyz, [], 1) - pad_size;
dxyz = max(xyz, [], 1) + pad_size - origin;
xyz_shift = xyz_shift - origin;
dfx = dxyz(1);
dfy = dxyz(2);
dfz = dxyz(3);

sx = max(imSizes(:, 2));
sy = max(imSizes(:, 1));
sz = max(imSizes(:, 3));
nxs = round(dfx/(px*xf) + sx + 2);
nys = round(dfy/(px*yf) + sy + 2);
nzs = round(dfz/(px*zf) + sz + 2);

nv = zeros(nys, nxs, nzs, dtype);
nv_f = zeros(nys, nxs, nzs, dtype);
if strcmp(BlendMethod, 'mean') || strcmp(BlendMethod, 'median')
    nv_count = zeros(nys, nxs, nzs, 'uint16');    
    nv_f_count = zeros(nys, nxs, nzs, 'uint16');
end

for i = 1 : nF
    tileFullpath = tileFullpaths{i};
    [dataPath, fsname] = fileparts(tileFullpath);
    dsrFullpath = [dsrpath '/' fsname '.tif'];
    if Decon && ~isempty(pr.DSRDeconDirstr)
        dsrFullpath = [dsrpath '/' fsname '_decon.tif'];
    end
    
    dsr = readtiff(dsrFullpath);
    if ~strcmp(class(dsr), dtype)
        dsr = cast(dsr, dtype);
    end
    
    if zNormalize
        for z = 1 : size(zmed_mat, 1)
            if zmed_mat(z, i) == 0, continue; end
            dsr_z = dsr(:, :, z);
            dsr_z(dsr_z ~= 0) = dsr_z(dsr_z ~= 0) + uint16(overall_z_median - zmed_mat(z, i));
            dsr(:, :, z) = dsr_z;
        end
    end
    
    st_idx = round(xyz_shift(i, :) ./ ([xf, yf, zf] * px));
    % if any(st_idx < 1) || any(st_idx > [nxs, nys, nzs])
    % bound idx to the positions of the stitched image
    sx = imSizes(i, 2);
    sy = imSizes(i, 1);
    sz = imSizes(i, 3);
    xridx = max(1, 1 - st_idx(1)) : min(sx, nxs - st_idx(1));
    yridx = max(1, 1 - st_idx(2)) : min(sy, nys - st_idx(2));
    zridx = max(1, 1 - st_idx(3)) : min(sz, nzs - st_idx(3));
        
    xidx = st_idx(1) + xridx;
    yidx = st_idx(2) + yridx;
    zidx = st_idx(3) + zridx;

    fprintf('Removing Edges...\n');
    r = EdgeArtifacts;
    if r > 0
        se = strel('disk', r);
        mask = dsr > 0;
        mask = imerode(mask, se);
        dsr(~mask) = 0;
    end
    dsr_raw = dsr;

    % remove discard region in half overlap setting
    for j = 1 : nF
        if ~overlap_matrix(i, j) && ~overlap_matrix(j, i)
            continue;
        end
        if i < j
            mregion = half_ol_region_cell{i, j}{1};
        else
            mregion = half_ol_region_cell{j, i}{2};
        end
        
        midx = mregion;
        % dsr_ol = dsr(mregion_f(2, 1) : mregion_f(2, 2), mregion_f(1, 1) : mregion_f(1, 2), mregion_f(3, 1) : mregion_f(3, 2));
        dsr(midx(2, 1) : midx(2, 2), midx(1, 1) : midx(1, 2), midx(3, 1) : midx(3, 2)) = 0;
    end
    
    switch BlendMethod
        case 'none'
            tmp = nv(yidx, xidx, zidx);
            tmp = dsr(yridx, xridx, zridx) .* cast(tmp == 0, dtype);
            nv(yidx, xidx, zidx) = nv(yidx, xidx, zidx) + tmp;
            nv_f(yidx, xidx, zidx) = nv_f(yidx, xidx, zidx) + dsr_raw(yridx, xridx, zridx) .* cast(nv_f(yidx, xidx, zidx) == 0, dtype);
        case 'mean'
            nv(yidx, xidx, zidx) = nv(yidx, xidx, zidx) + dsr(yridx, xridx, zridx);
            nv_count(yidx, xidx, zidx) = nv_count(yidx, xidx, zidx) + cast(dsr(yridx, xridx, zridx) ~= 0, dtype);
            
            nv_f(yidx, xidx, zidx) = nv_f(yidx, xidx, zidx) + dsr_raw(yridx, xridx, zridx);
            nv_f_count(yidx, xidx, zidx) = nv_f(yidx, xidx, zidx) + cast(dsr_raw(yridx, xridx, zridx) ~= 0, dtype);
        case 'max'
            nv(yidx, xidx, zidx) = max(nv(yidx, xidx, zidx), dsr(yridx, xridx, zridx));
            nv_f(yidx, xidx, zidx) = max(nv_f(yidx, xidx, zidx), dsr_raw(yridx, xridx, zridx));
        case 'median'
            nv_count(yidx, xidx, zidx) = nv_count(yidx, xidx, zidx) + cast(dsr(yridx, xridx, zridx) ~= 0, dtype);
            nv(yidx, xidx, zidx) = nv(yidx, xidx, zidx) + dsr(yridx, xridx, zridx);
            
            nv_f_count(yidx, xidx, zidx) = nv_f_count(yidx, xidx, zidx) + cast(dsr_raw(yridx, xridx, zridx) ~= 0, dtype);
            nv_f(yidx, xidx, zidx) = nv_f(yidx, xidx, zidx) + dsr_raw(yridx, xridx, zridx);
    end
end

% apply complement nv to nv if any voxel in nv is zero while nv_comp is not
clear dsr mask
nv_zero_inds = nv == 0 & nv_f ~= 0;
nv(nv_zero_inds) = nv_f(nv_zero_inds);
clear nv_f;
if strcmp(BlendMethod, 'mean') || strcmp(BlendMethod, 'median')
    nv_count(nv_zero_inds) = nv_f_count(nv_zero_inds);
    clear nv_f_count;
end
clear nv_zero_inds;

% for median option, first check the indices with overlap tiles
% to do: update median for filling blank regions in overlap regions.
if strcmp(BlendMethod, 'median') 
    if any(nv_count > 2, 'all')
        % nv_count_orig = nv_count;
        
        % nv_med, (:, :, 1) for half overlap, (:, :, 2) for full overlap
        % (:, 1, 1): ind, (:, 2, 1): nv_count, (:, 3, :) for current counts.
        inds = find(nv_count(:) > 2);
        nv_med = zeros(numel(inds), max(nv_count(:)) + 3, 2);
        nv_med(:, 1, 1) = inds;
        nv_med(:, 2, 1) = nv_count(inds);
        
        % nv_count = zeros(nys, nxs, nzs, 'uint16');
        for i = 1 : nF
            st_idx = round(xyz_shift(i, :) ./ ([xf, yf, zf] * px));

            % bound idx to the positions of the stitched image
            sx = imSizes(i, 2);
            sy = imSizes(i, 1);
            sz = imSizes(i, 3);            
            xridx = max(1, 1 - st_idx(1)) : min(sx, nxs - st_idx(1));
            yridx = max(1, 1 - st_idx(2)) : min(sy, nys - st_idx(2));
            zridx = max(1, 1 - st_idx(3)) : min(sz, nzs - st_idx(3));

            xidx = st_idx(1) + xridx;
            yidx = st_idx(2) + yridx;
            zidx = st_idx(3) + zridx;
                        
            if any(nv_count(yidx, xidx, zidx) > 2, 'all')
                tileFullpath = tileFullpaths{i};
                [dataPath, fsname] = fileparts(tileFullpath);

                dsrFullpath = [dsrpath '/' fsname '.tif'];
                if Decon && ~isempty(pr.DSRDeconDirstr)
                    dsrFullpath = [dsrpath '/' fsname '_decon.tif'];
                end
                dsr = readtiff(dsrFullpath);
                if ~strcmp(class(dsr), dtype)
                    dsr = cast(dsr, dtype);
                end

                if zNormalize
                    for z = 1 : size(zmed_mat, 1)
                        if zmed_mat(z, i) == 0, continue; end
                        dsr_z = dsr(:, :, z);
                        dsr_z(dsr_z ~= 0) = dsr_z(dsr_z ~= 0) + uint16(overall_z_median - zmed_mat(z, i));
                        dsr(:, :, z) = dsr_z;
                    end
                end
                
                r = EdgeArtifacts;
                if r > 0
                    se = strel('disk', r);
                    mask = dsr > 0;
                    mask = imerode(mask, se);
                    dsr(~mask) = 0;
                end
                dsr_raw = dsr;
                
                % remove discard region in half overlap setting
                for j = 1 : nF
                    if ~overlap_matrix(i, j) && ~overlap_matrix(j, i)
                        continue;
                    end
                    if i < j
                        mregion = half_ol_region_cell{i, j}{1};
                    else
                        mregion = half_ol_region_cell{j, i}{2};
                    end

                    midx = mregion;
                    dsr(midx(2, 1) : midx(2, 2), midx(1, 1) : midx(1, 2), midx(3, 1) : midx(3, 2)) = 0;        
                end
                
                nv_tmp = zeros(nys, nxs, nzs, dtype);
                nv_tmp(yidx, xidx, zidx) = dsr(yridx, xridx, zridx);
                
                nz_els = nv_tmp(inds);
                nz_inds = find(nz_els ~= 0);
                nv_med(nz_inds, 3, 1) = nv_med(nz_inds, 3, 1) + 1;
                nv_med_inds = sub2ind(size(nv_med), nz_inds, 3 + nv_med(nz_inds, 3, 1), 1 * ones(numel(nz_inds), 1));
                nv_med(nv_med_inds) = nz_els(nz_inds);
                
                % use dsr raw to decide filling elements
                nv_tmp(yidx, xidx, zidx) = dsr_raw(yridx, xridx, zridx);
                nz_els = nv_tmp(inds);
                nz_inds = find(nz_els ~= 0);
                nv_med(nz_inds, 3, 2) = nv_med(nz_inds, 3, 2) + 1;
                nv_med_inds = sub2ind(size(nv_med), nz_inds, 3 + nv_med(nz_inds, 3, 2), 2 * ones(numel(nz_inds), 1));
                nv_med(nv_med_inds) = nz_els(nz_inds);
            end
        end
        nv = nv ./ nv_count;
        % median for more than 2 tile overlap regions
        % med_inds = nv_med(:, 2) > 2;
        % inds_1 = inds(med_inds);
        % nv_med = nv_med(med_inds, :);
        nv_med(nv_med(:, 3, 1) == 0, 3 : end, 1) = nv_med(nv_med(:, 3, 1) == 0, 3 : end, 2);
        nv_med(:, :, 2) = [];
        uniq_counts = unique(nv_med(:, 3));
        nv_med_final = zeros(size(nv_med, 1), 1);
        for c = 1 : numel(uniq_counts)
            inds_c = nv_med(:, 3) == uniq_counts(c);
            nv_med_final(inds_c) = median(nv_med(inds_c, 4 : 3 + uniq_counts(c)), 2);
        end
        nv(inds) = nv_med_final;        
    else
        nv = nv ./ nv_count;
    end
end

switch BlendMethod
    case 'none'
    case 'mean'
        nv = nv ./ nv_count;
    case 'max'
    case 'median'
end

clear nv_count;

% crop or pad file, if pad if padSize positive, otherwise crop
if ~isempty(padSize)
    psz = padSize;
    if psz > 0
        nv_orig = nv;
        nv = zeros(nys + psz(1) * 2, nxs + psz(2) * 2, nzs + psz(3) * 2, 'uint16');
        nv(psz(1) + 1 : psz(1) + nys, psz(2) + 1 : psz(2) + nxs, psz(3) + 1 : psz(3) + nzs) = nv_orig;
        clear nv_orig;
    else
        nv(-psz(1) : nys + psz(1), -psz(2) : nxs + psz(2), -psz(3) : nzs + psz(3)) = [];
    end
end

% first check the size in secondary channels are the same as primary
% channel after crop
if isempty(boundboxCrop) && ~isPrimaryCh
    out_sz = [nys, nxs, nzs];
    % handle size difference:
    % 1. if out_sz smaller than primary sz, pad it.
    % 2. if out_sz larger than primary sz, crop it in the next step. 
    if any(out_sz < pImSz)
        disp('The size of stitched data in some axis is smaller than the primary channel, pad it...');
        pd_sz = pImSz - out_sz;
        pd_sz = pd_sz .* (pd_sz > 0);
        nv = padarray(nv, pd_sz, 0, 'post');
    end
    if any(size(nv) > pImSz)
        disp('The size of stitched data in some axis is larger than the primary channel, crop it...');        
        nv = nv(1 : pImSz(1), 1 : pImSz(2), 1 : pImSz(3));
    end
end

% crop to specific size by given bounding box
if ~isempty(boundboxCrop)
    bbox = boundboxCrop;
    if any(isinf(bbox(4 : 6)))
        bbox(4 : 6) = isinf(bbox(4 : 6)) .* size(nv);
    end
    
    out_sz = size(nv);
    if any(size(nv) < bbox(4 : 6))
        disp('The size of stitched data in some axis is smaller than the bournding box, pad it...');
        pd_sz = bbox(4 : 6) - out_sz;
        pd_sz = pd_sz .* (pd_sz > 0);
        nv = padarray(nv, pd_sz, 0, 'post');
    end
    
    nv = nv(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6));
end

% save stitch information for primary channels if there is xcorr shift
[~, fsname] = fileparts(tileFullpaths{1});
% out = regexp(fsname, '_?(Scan_.*)', 'tokens');
% if ~isempty(out)
%     fsname = out{1}{1};
% end

if isPrimaryCh 
    stichInfoPath = [dataPath, '/', ResultDir, '/', stitchInfoDir];
    if ~exist(stichInfoPath, 'dir')
        mkdir(stichInfoPath);
        fileattrib(stichInfoPath, '+w', 'g');
    end
    stitch_info_tmp_fullname = sprintf('%s/%s/%s/%s_%s.mat', dataPath, ResultDir, stitchInfoDir, fsname(1:end-21), uuid);
    pImSz = size(nv);
    save('-v7.3', stitch_info_tmp_fullname, 'ip', 'overlap_regions', ...
        'overlap_matrix', 'relative_shift_mat', 'pImSz', 'xf', 'yf', 'zf');
    movefile(stitch_info_tmp_fullname, sprintf('%s/%s/%s/%s.mat', dataPath, ResultDir, stitchInfoDir, fsname(1:end-21)));
end

% save MIP
if SaveMIP
    stcMIPPath = sprintf('%s/%s/MIPs/', dataPath, ResultDir);
    if ~exist(stcMIPPath, 'dir')
        mkdir(stcMIPPath);
        fileattrib(stcMIPPath, '+w', 'g');
    end
    stcMIPname = sprintf('%s%s_MIP_z.tif', stcMIPPath, fsname(1:end-21));
    writetiff(uint16(max(nv, [], 3)), stcMIPname);
end

% save stitched result
stitch_tmp_fullname = [dataPath '/' ResultDir '/' fsname(1:end-21) '_' uuid '.tif'];
stitch_fullname = [dataPath '/' ResultDir '/' fsname(1:end-21) '.tif'];
writetiff(nv, stitch_tmp_fullname);
movefile(stitch_tmp_fullname, stitch_fullname);

if false
    % write to hdf5 
    stitch_tmp_hdf5_fullname = [dataPath '/' ResultDir '/' fsname(1:end-21) '_' uuid '.h5'];
    writehdf5(nv, stitch_tmp_hdf5_fullname)
    movefile(stitch_tmp_hdf5_fullname, [dataPath '/' ResultDir '/' fsname(1:end-21) '.h5']);
end


if false
    % write to n5 
    stitch_tmp_n5_fullname = [dataPath '/' ResultDir '/' fsname(1:end-21) '_' uuid '_n5'];
    tiff2n5_cmd = 'python /global/home/groups/software/sl-7.x86_64/modules/stitching-spark/startup-scripts/spark-local/convert-tiff-tiles-n5.py';
    cmd = sprintf('%s -i %s -o %s', tiff2n5_cmd, stitch_fullname, stitch_tmp_n5_fullname);
    system(cmd);
    movefile(stitch_tmp_n5_fullname, [dataPath '/' ResultDir '/' fsname(1:end-21) '_n5']);
end


end
