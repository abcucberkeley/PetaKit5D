function [vol] = XR_ExM_PunctaRemovalPointDetection(vol, varargin)
% puncta removal based on point detection
% 
% xruan (02/24/2022): update thresholding method, and interpolate puncta
% region with surrounding region.
% xruan (07/14/2023): add gaussian smoothing for initial detection; add watershed to 
% define the region to segment the bright spot; using fillmissing2 to fill removed 
% regions, also fill them after segmet all bright spots

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('vol', @isnumeric);
ip.addParameter('Sigma', 2.5, @isnumeric); % for 3D Gauss filtering
ip.addParameter('OTSUMaxPer', 99.9, @isnumeric); % Max percentile of data for OTSU calculation
ip.addParameter('localWinSize', [41, 41, 41], @isvector); % local window size for cropping local region
ip.addParameter('SigmaThrsh', 4, @isvector); % sigma threshold for point detection
ip.addParameter('intThrsh', 5000, @isvector); % intensity threshold for the peak to be removed
ip.addParameter('initDetect', '2d', @ischar); % initial detection method, 2d mip or 3d stack
ip.addParameter('detVolThrsh', 5000, @isnumeric); % volume threshold for the detetion
ip.addParameter('debug', false, @islogical);

ip.parse(vol, varargin{:});

pr = ip.Results;
Sigma = pr.Sigma;
OTSUMaxPer = pr.OTSUMaxPer;
localWinSize = pr.localWinSize;
SigmaThrsh = pr.SigmaThrsh;
intThrsh = pr.intThrsh;
initDetect = pr.initDetect;
detVolThrsh = pr.detVolThrsh;

hbbox = round((localWinSize - 1) / 2);

vol_sm = imgaussfilt3(vol, Sigma);

switch initDetect
    case '2d'
        [mip, zinds] = max(vol_sm, [], 3);

        BW = imregionalmax(mip);

        [py, px] = find(BW);
        pz = zinds(BW(:));

        % use local thresholding to keep detected peaks
        std_thrsh = SigmaThrsh;
        kept_flag = false(numel(py), 1);
        im_std = std(mip(:));
        for i = 1 : numel(py)
            yi = max(1, py(i) - hbbox(1)) : min(size(mip, 1), py(i) + hbbox(1));
            xi = max(1, px(i) - hbbox(2)) : min(size(mip, 2), px(i) + hbbox(2));

            im_i = mip(yi, xi);

            peak_i = mip(py(i), px(i));

            bg_i = median(im_i(:));
            std_i = std(im_i(:));

            if peak_i - bg_i > std_thrsh * (std_i * 0.95 + im_std * .05)
                kept_flag(i) = true;
            end
        end

        if false
            figure, imshow(mip, [])
            hold on, plot(px(kept_flag), py(kept_flag), 'o')
        end        
    case '3d'
        BW = XR_locmax3d(vol_sm, [5, 5, 5]);

        pinds = find(BW > 0);
        [py, px, pz] = ind2sub(size(BW), pinds);

        % hbbox = [7, 7, 7];
        std_thrsh = SigmaThrsh;
        kept_flag = false(numel(py), 1);
        im_std = std(vol_sm(vol_sm > 0));
        sz = size(vol_sm);
        for i = 1 : numel(py)
            yi = max(1, py(i) - hbbox(1)) : min(sz(1), py(i) + hbbox(1));
            xi = max(1, px(i) - hbbox(2)) : min(sz(2), px(i) + hbbox(2));
            zi = max(1, pz(i) - hbbox(3)) : min(sz(3), pz(i) + hbbox(3));

            im_i = vol_sm(yi, xi, zi);

            peak_i = vol_sm(py(i), px(i), pz(i));

            bg_i = median(im_i(im_i > 0));
            std_i = std(im_i(im_i > 0));

            if peak_i - bg_i > std_thrsh * (std_i * 0.95 + im_std * .05)
                kept_flag(i) = true;
            end
        end
end

if ~any(kept_flag)
    return;
end

%% segment the region with bright spot

pyc = py(kept_flag);
pxc = px(kept_flag);
pzc = pz(kept_flag);

pcoords = [pyc, pxc, pzc];
% hbbox = [7, 7, 11];

rm_inds = cell(size(pcoords, 1), 1);
% fvals = cell(size(pcoords, 1), 1);
sz = size(vol);
for i = 1 : size(pcoords, 1)
    tic
    peak_i = max(vol_sm(pyc(i), pxc(i), pzc(i)), vol(pyc(i), pxc(i), pzc(i)));
    if peak_i < intThrsh
        continue;
    end
    
    yi = max(1, pyc(i) - hbbox(1)) : min(sz(1), pyc(i) + hbbox(1));
    xi = max(1, pxc(i) - hbbox(2)) : min(sz(2), pxc(i) + hbbox(2));
    zi = max(1, pzc(i) - hbbox(3)) : min(sz(3), pzc(i) + hbbox(3));
        
    pyi = pyc(i) - yi(1) + 1;
    pxi = pxc(i) - xi(1) + 1;
    pzi = pzc(i) - zi(1) + 1;
    
    im_i = vol_sm(yi, xi, zi);
    % im_gi = imgaussfilt3(im_i, Sigma);
    
    % T_i = thresholdOtsu(im_gi(im_gi > 0 & im_gi < prctile(im_gi(:), OTSUMaxPer)));
    % im_Ti = im_gi > T_i;
    % T_i = thresholdOtsu(im_i(im_i > median(im_i(im_i > 0))));
    % im_Ti = im_i > min(T_i, intThrsh);
    T_i = thresholdOtsu(im_i(im_i > 0));
    T_i_s = (min(T_i, intThrsh) * 0.2 + median(im_i(:)) * 0.8);
    im_Ti = im_i > T_i_s;

    L_i = bwlabeln(im_Ti);
    
    if L_i(pyi, pxi, pzi) == 0
        continue;
    end
    
    BW_i = L_i == L_i(pyi, pxi, pzi);

    alpha = 1.1;
    while sum(BW_i(:)) > detVolThrsh
        T_i_s = T_i_s * alpha;
        im_Ti = im_i > T_i_s;
        
        L_i = bwlabeln(im_Ti);
        if L_i(pyi, pxi, pzi) == 0
            break;
        end
        BW_i = L_i == L_i(pyi, pxi, pzi);
    end
    if L_i(pyi, pxi, pzi) == 0
        continue;
    end    

    binds = find(BW_i);

    % watershed segmentation
    im_i_1 = max(im_i(:)) + 1 - im_i;
    im_i_1(1 : end, 1 : end, [1, end]) = 1;
    im_i_1(1 : end, [1, end], 1 : end ) = 1;
    im_i_1([1, end], 1 : end, 1 : end) = 1;
    L_w = watershed(im_i_1);
    BW_w = L_w ==  L_w(pyi, pxi, pzi);
    
    BW_i = BW_i & imdilate(BW_w, strel('disk', 1));
    BW_i = imopen(BW_i, strel('disk', 1));
    
    [byi, bxi, bzi] = ind2sub(size(BW_i), binds);
    byi = byi + yi(1) - 1;
    bxi = bxi + xi(1) - 1;
    bzi = bzi + zi(1) - 1;

    % use linear interpolation to reset the values in puncta
    % im_i(BW_i) = NaN;
    % [im_i, TF] = fillmissing(im_i, 'linear', 2, 'EndValues', 'nearest');

    rm_inds{i} = sub2ind(sz, byi, bxi, bzi);
    % fvals{i} = im_i(BW_i);    
    toc
end

rinds = cat(1, rm_inds{:});
if isempty(rinds)
    return;
end
clear vol_sm;

%% removing bright spots slice by slice

borderSize = [3, 3];

BW_s = false(sz);
BW_s(rinds) = true;

% exclude large volumes to remove the bright spots
% BW_s_rm = bwareaopen(BW_s, detVolThrsh);
stats = regionprops3(BW_s, 'voxelIdxList', 'volume', 'PrincipalAxisLength');
stats = stats(stats.Volume > detVolThrsh * 2 | stats.PrincipalAxisLength(:, 1) > 35, :);
VoxelIdxList = cat(1, cat(1, stats.VoxelIdxList{:}));
BW_s(VoxelIdxList) = false;

% vol(rinds) = NaN;
for z = 1 : sz(3)
    tic
    BW_z = BW_s(:, :, z);
    
    stats = regionprops(BW_z, 'BoundingBox');
    if isempty(stats)
        continue;
    end

    BoundingBox = cat(1, stats.BoundingBox);
    bboxes = [ceil(BoundingBox(:, 2)), ceil(BoundingBox(:, 1)), ceil(BoundingBox(:, 4)), ceil(BoundingBox(:, 3))];
    bboxes(:, 3 : 4) = bboxes(:, 1: 2) + bboxes(:, 3: 4) - 1;
    
    rm_inds_z = cell(size(bboxes, 1), 1);
    fvals_z = cell(size(bboxes, 1), 1);

    for i = 1 : size(bboxes, 1)
        bbox_i = bboxes(i, :);
        yi = max(1, bbox_i(1) - borderSize(1)) : min(sz(1), bbox_i(3) + borderSize(1));
        xi = max(1, bbox_i(2) - borderSize(2)) : min(sz(2), bbox_i(4) + borderSize(2));
        zi = z;
                
        im_i = vol(yi, xi, zi);
        BW_zi = BW_z(yi, xi);
        im_i = fillmissing2(im_i, 'linear', MissingLocations=BW_zi);

        binds = find(BW_zi);
        [byi, bxi] = ind2sub(size(BW_zi), binds);
        byi = byi + yi(1) - 1;
        bxi = bxi + xi(1) - 1;

        rm_inds_z{i} = sub2ind(sz(1 : 2), byi, bxi);
        fvals_z{i} = im_i(BW_zi);    
    end
    
    vol_z = vol(:, :, z);
    rinds_z = cat(1, rm_inds_z{:});
    fvals_z = cat(1, fvals_z{:});
    vol_z(rinds_z) = fvals_z;
    vol(:, :, z) = vol_z;
    toc
end


end

