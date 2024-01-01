%runDetection3D(data) detects CCPs using a combination of model-based (PSF) fitting and statistical tests
%
% Inputs:
%      data : required, structure returned by loadConditionData()
%
% Options ('specifier', value):
%                    'Sigma' : sigma for the psf, N X 2 for N channels with form [sigma_xy, sigma_z] in each row. 
%                     'Mode' : estimate x, y, z, A (intensity) and c (background).  
%              'FitMixtures' : fit signle Gaussian or Mixture of Gaussian, default: false
%              'MaxMixtures' : number of mixture if choosing to fit mixtures. 
%          'RemoveRedundant' : remove redundant points if multiple points are very close in master channel. Default: true
%         'RedundancyRadius' : distance to define redundant points, scalar (total distance) or 1X2 vector (distances in xy-plane and z direction).
% 'RemoveSecondaryRedundant' : remove redundant points in secondary channels.
%                'Overwrite' : overwrite if the result exists, defalt: false
%                   'Master' : master channel. Default: [] (first channel). 
%                    'Alpha' : threshold of p-value for significant level. 
%                 'CellMask' : cell mask.
%                 'FileName' : detection result filename. 
%               'WindowSize' : window size for detection, default: [], and it is ceiling of 4 folds of sigma. 
%              'ResultsPath' : Result directory. 
%          'DetectionMethod' : detection method chosen from 'original', 'speedup', 'lowSNR', 'adaptive'. 
%                               'original': original method in cmeAnalysis3D package
%                               'speedup': the speedup of the original method
%                               'lowSNR' (default): the low-SNR detection methods, 
%                               'adaptive' is for the adative method between speedup and lowSNR 
%               'WriteAmira' : write detection results to amira file for visualization. Default: true
%     'BackgroundCorrection' : adjust illumination across z-stacks. Default: false
%    'FilterByResidualSigma' : remove noisy detected points by background std. Default: true
%
% Notes:
%  - 3D coordinates are in pixel space, except for the zCoord field passed
%    to the tracker, which is corrected for z-anisotropy
%  - All input data sets must have the same channels.

% Francois Aguet, 08/2013
% modifed Gokul U, 03/2017 -- added mixture model fitting option
% Xiongtao Ruan, 10/09/2019 add option to write to amira 
% Xiongtao Ruan 11/04/2019 add option for background correction.
% Xiongtao Ruan 11/06/2019 skip it if write amira file fails
% Xiongtao Ruan 11/08/2019 add option for filtering detection by residual
% sigma compared to DC standard deviation
% Xiongtao Ruan 11/13/2019 update the detection for slave channels
% Xiongtao Ruan 11/20/2019 add option for removing redundant points in
% slave channels
% Xiongtao Ruan 12/18/2019 change to use flag-based parallelize for
% frame-level parallelization 


function XR_runDetection3D_frame_parallel(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
% ip.addOptional('destPath', [], @ischar);
ip.addParameter('Sigma', []);%, @(x) numel(x)==length(data(1).channels));
% ip.addParamValue('SigmaSource', 'model', @(x) any(strcmpi(x, {'data', 'model'})));
ip.addParameter('Mode', 'xyzAc', @ischar);
ip.addParameter('FitMixtures', false, @islogical);
ip.addParameter('MaxMixtures', 3, @(x) numel(x)==1 && x>0 && round(x)==x);
ip.addParameter('RemoveRedundant', true, @islogical);
ip.addParameter('RedundancyRadius', 0.5, @(x) numel(x) <= 2 && all(x > 0));
ip.addParameter('RemoveSecondaryRedundant', true, @islogical);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('Master', [], @isnumeric);
ip.addParameter('Alpha', 0.05, @isnumeric);
ip.addParameter('CellMask', [],  @isnumeric);
ip.addParameter('FileName', 'Detection3D.mat');
ip.addParameter('WindowSize', []);
ip.addParameter('ResultsPath', arrayfun(@(i) [i.source 'Analysis' filesep], data, 'unif', 0), @iscell);
ip.addParameter('DetectionMethod', 'update_1', @ischar);
ip.addParameter('WriteAmira', true, @islogical);
ip.addParameter('BackgroundCorrection', false, @islogical); % background correction
ip.addParameter('FilterByResidualSigma', ~false, @islogical); % filter detection by residual sigma w.r.t DC std
ip.parse(data, varargin{:});
overwrite = ip.Results.Overwrite;
mCh = ip.Results.Master;
if isempty(mCh)
    mCh = unique(arrayfun(@(i) find(strcmpi(i.channels, i.source)), data));
    if numel(mCh)>1
        error('Master channel index is inconsistent accross data sets.');
    end
end

sigma = ip.Results.Sigma;

if any(sigma<1.1)
    fprintf(2, 'Sigma values < 1.1 were rounded to 1.1 to avoid poor localization performance.\n');
    sigma(sigma<1.1) = 1.1;
end

for i = 1:length(data)
    if ~(exist([ip.Results.ResultsPath{i} ip.Results.FileName], 'file') == 2) || overwrite
        fprintf('Running detection for %s ...', getShortPath(data(i)));

        is_done = main(data(i), sigma, mCh, ip.Results.ResultsPath{i}, ip.Results);
        if is_done
            fprintf(' done.\n');
        end
    else
        fprintf('Detection has already been run for %s\n', getShortPath(data(i)));
    end
end


function [is_done] = main(data, sigma, mCh, resultsPath, opts)

% xruan 11/08/2019 if need to filter by residual sigma, get DC std from ex setting file
FilterByResidualSigma = opts.FilterByResidualSigma;
if opts.FilterByResidualSigma
    ResidualSigmaThresh_device = XR_get_camera_standard_deviation_threshold(data);
end

% master channel
nCh = length(data.channels);
if opts.FitMixtures
    frameInfo(1:data.movieLength) = struct('x', [], 'y', [], 'z', [], 'A', [], 's', [], 'c', [],...
        'x_pstd', [], 'y_pstd', [], 'z_pstd', [], 'A_pstd', [], 'c_pstd', [],...
        'x_init', [], 'y_init', [], 'z_init', [],...
        'sigma_r', [], 'SE_sigma_r', [], 'RSS', [], 'pval_Ar', [], 'hval_Ar', [],  'hval_AD', [], 'isPSF', [],'mixtureIndex', [],...
        'xCoord', [], 'yCoord', [], 'zCoord', [], 'amp', [], 'dRange', []);
else
    frameInfo(1:data.movieLength) = struct('x', [], 'y', [], 'z', [], 'A', [], 's', [], 'c', [],...
        'x_pstd', [], 'y_pstd', [], 'z_pstd', [], 'A_pstd', [], 'c_pstd', [],...
        'x_init', [], 'y_init', [], 'z_init', [],...
        'sigma_r', [], 'SE_sigma_r', [], 'RSS', [], 'pval_Ar', [], 'hval_Ar', [],  'hval_AD', [], 'isPSF', [],...
        'xCoord', [], 'yCoord', [], 'zCoord', [], 'amp', [], 'dRange', []);
end

fmt = ['%.' num2str(ceil(log10(data.movieLength+1))) 'd'];
[~,~] = mkdir(resultsPath);
detectionResultDir = [resultsPath 'Detections'];
[~,~] = mkdir([resultsPath 'Detections']);
[~,~] = mkdir([resultsPath 'Masks']);
DetectionFullName = [resultsPath opts.FileName];

% double fields, multi-channel
dfields = {'x', 'y', 'z', 'A', 'c', 'x_pstd', 'y_pstd', 'z_pstd', 'A_pstd', 'c_pstd', 'sigma_r', 'SE_sigma_r', 'RSS', 'pval_Ar'};
% logical fields, multi-channel
lfields = {'hval_Ar', 'hval_AD', 'isPSF'};
% slave channel fields
sfields = [dfields {'hval_Ar', 'hval_AD'}]; % 'isPSF' is added later

rmfields = [dfields lfields {'x_init', 'y_init', 'z_init', 'mask_Ar'}];

DetectionMethod = opts.DetectionMethod;
if strcmp(DetectionMethod, 'update_3')
    kLevel = norminv(1 - opts.Alpha, 0, 1);
end

is_done = false;
% disable overwrite
if exist(DetectionFullName, 'file')
    is_done = true;
    return;
end

all_frame_computed = false(data.movieLength, 1);
% uuid for the job
temp = java.util.UUID.randomUUID;
myuuid = temp.toString;

% parfor k = 1:data.movieLength % parfor
% loop several times, so that all the frames for current data can be done. 
n_loops = 3;
for i = 1 : n_loops
    chunk_lock_clean(detectionResultDir, 2);    
    for k = 1:data.movieLength % parfor
        % in case the final detection file is generated in the middle of
        % computing
        if exist(DetectionFullName, 'file')
            is_done = true;            
            return;
        end

        cur_detection_filename = sprintf('%s/Detection_F_%d.mat', detectionResultDir, k);
        temp_filename = sprintf('%s/Detection_F_%d.tmp', detectionResultDir, k);
        if exist(cur_detection_filename, 'file')
            all_frame_computed(k) = true;
            continue;
        end
        if exist(temp_filename, 'file')
            continue;
        else
            fclose(fopen(temp_filename, 'w'));
        end
        fprintf('Compute Detection for Frame %d...\n', k);

        frame = double(readtiff(data.framePathsDS{mCh}{k})); %#ok<PFBNS>
        frame(frame==0) = NaN; % to avoid border effects

        [ny,nx,nz] = size(frame);
        switch DetectionMethod
            case 'original'
                [pstruct, mask] = pointSourceDetection3D(frame, sigma(mCh,:), 'Alpha', opts.Alpha,...
                    'Mask', opts.CellMask, 'RemoveRedundant', opts.RemoveRedundant, 'RedundancyRadius', opts.RedundancyRadius, ...
                    'RefineMaskLoG', false, 'WindowSize', opts.WindowSize, 'Mode', opts.Mode,...
                    'FitMixtures', opts.FitMixtures,'MaxMixtures', opts.MaxMixtures); %#ok<PFBNS>
            case 'speedup'
                [pstruct, mask] = XR_pointSourceDetection3D(frame, sigma(mCh,:), 'Alpha', opts.Alpha,...
                    'Mask', opts.CellMask, 'RemoveRedundant', opts.RemoveRedundant, 'RedundancyRadius', opts.RedundancyRadius, ...
                    'RefineMaskLoG', false, 'WindowSize', opts.WindowSize, 'Mode', opts.Mode,...
                    'FitMixtures', opts.FitMixtures,'MaxMixtures', opts.MaxMixtures, ...
                    'BackgroundCorrection', opts.BackgroundCorrection);
            case 'update_1'
                [pstruct, mask] = XR_pointSourceDetection3D_1(frame, sigma(mCh,:), 'Alpha', opts.Alpha,...
                    'Mask', opts.CellMask, 'RemoveRedundant', opts.RemoveRedundant, 'RedundancyRadius', opts.RedundancyRadius, ...
                    'RefineMaskLoG', false, 'WindowSize', opts.WindowSize, 'Mode', opts.Mode,...
                    'FitMixtures', opts.FitMixtures,'MaxMixtures', opts.MaxMixtures);
            case 'lowSNR'    % rigorous update  
                % xruan 11/08/2019
                if FilterByResidualSigma
                    [backgroundInfo] = XR_estimate_overall_background_information(frame, 'randseed', k);
                    bg_std_factor = 1.0;                
                    if isnan(ResidualSigmaThresh_device) || ResidualSigmaThresh_device < backgroundInfo(3) * bg_std_factor
                        % xruan 12/06/2019 if the estimate sigma is too large, which means the
                        % backgroud of the image may be uneven, in this case,
                        % we use the sigma from device. 
                        if ResidualSigmaThresh_device * 1.5 > (backgroundInfo(3) * bg_std_factor)
                            ResidualSigmaThresh =  backgroundInfo(3) * bg_std_factor;
                        else
                            ResidualSigmaThresh = ResidualSigmaThresh_device;
                        end
                    else
                        % If the estimated background is very small, use the
                        % maximun between relaxed device r std, and the
                        % estimated std. 
                        if ResidualSigmaThresh_device > 1.5 * (backgroundInfo(3) * bg_std_factor)
                            ResidualSigmaThresh = max(ResidualSigmaThresh_device * 0.9, backgroundInfo(3) * bg_std_factor);
                        else
                            ResidualSigmaThresh = ResidualSigmaThresh_device;
                        end
                    end
                    backgroundInfo = [backgroundInfo, ResidualSigmaThresh_device];
                end
                [pstruct, mask] = XR_pointSourceDetection3D_3(frame, sigma(mCh,:), 'Alpha', opts.Alpha,...
                    'Mask', opts.CellMask, 'RemoveRedundant', opts.RemoveRedundant, 'RedundancyRadius', opts.RedundancyRadius, ...
                    'RefineMaskLoG', false, 'WindowSize', opts.WindowSize, 'Mode', opts.Mode,...
                    'FitMixtures', opts.FitMixtures,'MaxMixtures', opts.MaxMixtures, ...
                    'BackgroundCorrection', opts.BackgroundCorrection, ...
                    'FilterByResidualSigma', FilterByResidualSigma, 'ResidualSigmaThresh', ResidualSigmaThresh, ...
                    'BackgroundInfo', backgroundInfo);
            case 'adaptive'    % adative method
                [pstruct, mask] = XR_pointSourceDetection3D_4(frame, sigma(mCh,:), 'Alpha', opts.Alpha,...
                    'Mask', opts.CellMask, 'RemoveRedundant', opts.RemoveRedundant, 'RedundancyRadius', opts.RedundancyRadius, ...
                    'RefineMaskLoG', false, 'WindowSize', opts.WindowSize, 'Mode', opts.Mode,...
                    'FitMixtures', opts.FitMixtures,'MaxMixtures', opts.MaxMixtures);
        end

        if ~isempty(pstruct)
            pstruct.s = sigma;
            pstruct = rmfield(pstruct, 's_pstd');

            pstruct.dRange{mCh} = [min(frame(:)) max(frame(:))];
            np = numel(pstruct.x);

            % expand structure for slave channels
            for f = 1:length(dfields)
                tmp = NaN(nCh, np);
                tmp(mCh,:) = pstruct.(dfields{f});
                pstruct.(dfields{f}) = tmp;
            end
            for f = 1:length(lfields)
                tmp = false(nCh, np);
                tmp(mCh,:) = pstruct.(lfields{f});
                pstruct.(lfields{f}) = tmp;
            end

            % retain only mask regions containing localizations
            CC = bwconncomp(mask);
            labels = labelmatrix(CC);
            loclabels = labels(sub2ind([ny nx nz], pstruct.y_init, pstruct.x_init, pstruct.z_init));
            idx = setdiff(1:CC.NumObjects, loclabels);
            CC.PixelIdxList(idx) = [];
            CC.NumObjects = length(CC.PixelIdxList);

            % clean mask
            labels = labelmatrix(CC);
            mask = labels~=0;

            for ci = setdiff(1:nCh, mCh)
                frame = double(readtiff(data.framePathsDS{ci}{k}));
                % xruan 11/12/2019 set 0 as NaN
                frame(frame == 0) = nan;            
                if FilterByResidualSigma
                    [backgroundInfo] = XR_estimate_overall_background_information(frame, 'randseed', k);
                    bg_std_factor = 1.3;
                    if isnan(ResidualSigmaThresh_device) || ResidualSigmaThresh_device < backgroundInfo(3) * bg_std_factor
                        ResidualSigmaThresh =  backgroundInfo(3) * bg_std_factor;
                    else
                        ResidualSigmaThresh = ResidualSigmaThresh_device;
                    end
                end                                

                pstruct.dRange{ci} = [min(frame(:)) max(frame(:))];
                X = [pstruct.x(mCh,:)' pstruct.y(mCh,:)' pstruct.z(mCh,:)'];

                if strcmp(DetectionMethod, 'original')
                    % linear index of positions            
                    linIdx = sub2ind(size(frame), roundConstr(X(:,2),ny), roundConstr(X(:,1),nx), roundConstr(X(:,3),nz));

                    % localize, and compare intensities & (x,y)-coordinates. Use localization result if it yields better contrast
                    [A_est, c_est] = estGaussianAmplitude3D(frame, sigma(ci,:));
                    pstructSlave = fitGaussians3D(frame, X, A_est(linIdx), sigma(ci,:), c_est(linIdx), 'Ac');

                    pstructSlaveLoc = fitGaussians3D(frame, X, pstructSlave.A, sigma(ci,:), pstructSlave.c, 'xyzAc');
                else
                    A_est_init = nan(size(X, 1), 1);
                    c_est_init = nan(size(X, 1), 1);

                    pstructSlave = XR_fitGaussians3D(frame, X, A_est_init, sigma(ci,:), c_est_init, 'Ac', 'FitGaussianMethod', 'CD_Newton');

                    pstructSlaveLoc = XR_fitGaussians3D(frame, X, pstructSlave.A, sigma(ci,:), pstructSlave.c, 'xyzAc');
                end

                % xruan 11/20/2019 change threshold to 2 sigma from 3 sigma
                idx = sqrt((pstruct.x(mCh,:)-pstructSlaveLoc.x).^2 ...
                    + (pstruct.y(mCh,:)-pstructSlaveLoc.y).^2) < 2*sigma(mCh,1) & pstructSlaveLoc.A > pstructSlave.A;

                % xruan 11/20/2019
                % check whether different points in master channel detect same points in slave channel, 
                % if so, only keep the nearest pair, and other points use pstructSlave
                if opts.RemoveSecondaryRedundant
                    pM = [pstructSlaveLoc.x' pstructSlaveLoc.y' pstructSlaveLoc.z'];
                    % only consider chosen points
                    pM = pM(idx, :);
                    pIdx = find(idx);
                    % xruan: change to allow radius for xy and z
                    RedundancyRadius = opts.RedundancyRadius;
                    if numel(RedundancyRadius) == 1
                        idxKD = KDTreeBallQuery(pM, pM, RedundancyRadius*ones(size(pM, 1),1));
                    elseif numel(RedundancyRadius) == 2
                        idxKD_xy = KDTreeBallQuery(pM(:, 1 : 2), pM(:, 1 : 2), RedundancyRadius(1)*ones(size(pM, 1),1));
                        idxKD = arrayfun(@(x) idxKD_xy{x}((abs(pM(x, 3) - pM(idxKD_xy{x}, 3)) < RedundancyRadius(2))), 1 : size(pM, 1), 'unif', false)';                    
                    end

                    % pick out the one with minimum RSS in secondary channel,
                    % and then find the point with closest distances in the
                    % primary channel, and pair them, and set other points in
                    % pstructSlaveLoc as removed. 
                    idxKD = idxKD(cellfun(@numel, idxKD)>1);
                    reassign_mat = zeros(length(idxKD), 2);
                    for iKD = 1:length(idxKD)
                        cur_idxKD = idxKD{iKD};
                        RSS = pstructSlaveLoc.RSS(pIdx(cur_idxKD));
                        % idx(idxKD{iKD}(RSS ~= min(RSS))) = 0;
                        [~, min_ind] = min(RSS);
                        pM_s = pM(cur_idxKD(min_ind), :);
                        [~, primary_min_ind] = min(sum((pM_s - X(pIdx(cur_idxKD), :)) .^ 2, 2));

                        if primary_min_ind ~= min_ind
                            reassign_mat(iKD, :) = [pIdx(cur_idxKD(primary_min_ind)), pIdx(cur_idxKD(min_ind))];
                        end
                        idx(pIdx(cur_idxKD([1 : numel(cur_idxKD)] ~= primary_min_ind))) = 0;
                    end

                    % reassign points in Secondary channel for better detection
                    % performance
                    reassign_mat(reassign_mat(:, 1) == 0, :) = [];
                    for f = 1:length(sfields)
                        pstructSlaveLoc.(sfields{f})(reassign_mat(:, 1)) = pstructSlaveLoc.(sfields{f})(reassign_mat(:, 2));
                    end
                end

                % fill slave channel information
                for f = 1:length(sfields)
                    pstruct.(sfields{f})(ci,~idx) = pstructSlave.(sfields{f})(~idx);
                    pstruct.(sfields{f})(ci,idx) = pstructSlaveLoc.(sfields{f})(idx);
                end

                nanIdx = isnan(pstructSlave.x); % points within slave channel border, remove from detection results
                for f = 1:length(rmfields)-1 %changed by GU
                    pstruct.(rmfields{f})(:,nanIdx) = []; % GU - error
                end

                pstruct.isPSF(ci,:) = ~pstruct.hval_AD(ci,:);
            end

            % add fields for tracker
            pstruct.xCoord = [pstruct.x(mCh,:)' pstruct.x_pstd(mCh,:)'];
            pstruct.yCoord = [pstruct.y(mCh,:)' pstruct.y_pstd(mCh,:)'];
            pstruct.zCoord = [pstruct.z(mCh,:)' pstruct.z_pstd(mCh,:)'] * data.zAniso;
            pstruct.amp =    [pstruct.A(mCh,:)' pstruct.A_pstd(mCh,:)'];
            % frameInfo(k) = orderfields(pstruct, fieldnames(frameInfo(k))); %#ok<PFOUS>
            frameInfo_k = orderfields(pstruct, fieldnames(frameInfo(k))); 
        else
            % frameInfo(k).dRange{mCh} = [min(frame(:)) max(frame(:))];
            frameInfo_k.dRange{mCh} = [min(frame(:)) max(frame(:))];
            for ci = setdiff(1:nCh, mCh)
                frame = double(imread(data.framePathsDS{ci}{k}));
                % frameInfo(k).dRange{ci} = [min(frame(:)) max(frame(:))];
                frameInfo_k.dRange{ci} = [min(frame(:)) max(frame(:))];
            end
            % frameInfo(k).s = sigma;
            frameInfo_k.s = sigma;
        end

        maskPath = [resultsPath filesep 'Masks' filesep 'dmask_' num2str(k, fmt) '.tif'];
        if ~exist(maskPath, 'file')
            writetiff(uint8(255*mask), maskPath);
        end

        if ~exist(cur_detection_filename, 'file')
            cur_detection_filename_tmp = sprintf('%s/Detection_F_%d_%s.mat', detectionResultDir, k, myuuid);
            
            save('-v7.3', cur_detection_filename_tmp, 'frameInfo_k');
            % then move the temp result file to the result file
            if ~exist(cur_detection_filename, 'file')
                movefile(cur_detection_filename_tmp, cur_detection_filename);
            else
                delete(cur_detection_filename_tmp);
            end
        end
        if exist(cur_detection_filename, 'file') && exist(temp_filename, 'file')
            all_frame_computed(k) = true;                                    
            delete(temp_filename);
        end    
    end
    if all(all_frame_computed)
        break;
    end
    pause(5);
end

if ~all(all_frame_computed) 
    is_done = false;
    return;
end

% collect all results
if ~exist(DetectionFullName, 'file')
    for k = 1 : data.movieLength
        cur_detection_filename = sprintf('%s/Detection_F_%d.mat', detectionResultDir, k);
        try
            load(cur_detection_filename, 'frameInfo_k');
        catch 
            delete(cur_detection_filename);
            return;
        end
        
        if all(contains(dfields, fieldnames(frameInfo_k)))
            frameInfo(k) = frameInfo_k;
        else
            fieldnames_k = fieldnames(frameInfo_k);
            for nf = 1 : numel(fieldnames_k)
                frameInfo(k).(fieldnames_k{nf}) = frameInfo_k.(fieldnames_k{nf});               
            end            
        end
    end
    % save('-v7.3', DetectionFullName, 'frameInfo');
end

% xruan 10/09/2019
if opts.WriteAmira
    amiraRootPath = [resultsPath filesep 'Amira'];
    mkdir(amiraRootPath);
    amiraPath = [amiraRootPath filesep 'Detection'];    
    mkdir(amiraPath);
    
    % xruan 11/06/2019 use try catch if it fails to write amira file
    % (contains empty detection)
    try
        XR_amira_write_detection_spots([amiraPath filesep 'detection'], frameInfo, data, 'FlipZ', false)
    catch
        warning('Write to Amira failed, Skip it!')
    end
end
   
% save('-v7.3', [resultsPath opts.FileName], 'frameInfo');
if ~exist(DetectionFullName, 'file')
    % save('-v7.3', DetectionFullName, 'frameInfo');
    DetectionFullName_tmp = sprintf('%s_%s%s', DetectionFullName(1 : end - 5), myuuid, DetectionFullName(end - 4 : end));

    save('-v7.3', DetectionFullName_tmp, 'frameInfo');
    % then move the temp result file to the result file
    if ~exist(DetectionFullName, 'file')
        movefile(DetectionFullName_tmp, DetectionFullName);
    else
        delete(DetectionFullName_tmp);
    end
end


function y = roundConstr(x, nx)
y = round(x);
y(y<1) = 1;
y(y>nx) = nx;
