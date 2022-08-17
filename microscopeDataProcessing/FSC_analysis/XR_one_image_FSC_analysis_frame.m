function [fsc_mu, res_mu, fsc_cell, res_cell] = XR_one_image_FSC_analysis_frame(fn, fnout, varargin)
% perform FSC analysis for a single image for given image filename and save
% results to disk
% 
% xruan (10/19/2021): add support for bounding box for the region to calculate FSC
% xruan (04/14/2022): add option for skip the cone regions along axis
% xruan (05/12/2022): add support for FSC for multiple bboxes if bbox has multiple rows
% xruan (06/09/2022): add support for clipping very bright spots

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('fn', @(x) ischar(x) || isnumeric(x));
ip.addParameter('xyPixelSize', 0.108, @isnumeric);
ip.addParameter('dz', 0.108, @isnumeric);
ip.addParameter('dr', 1 , @isnumeric);
ip.addParameter('dtheta', pi / 12 , @isnumeric);
ip.addParameter('resThreshMethod', 'fixed', @ischar);
ip.addParameter('resThresh', 0.2, @isnumeric);
ip.addParameter('N', [501, 501, 501], @isnumeric);
ip.addParameter('bbox', [], @isnumeric);
ip.addParameter('resAxis', 'xz', @ischar);
ip.addParameter('skipConeRegion', true, @islogical);
ip.addParameter('clipPer', [], @isnumeric); % clip intensity higher than the given percentile
ip.addParameter('debug', false, @islogical);

ip.parse(fn,  varargin{:});

pr = ip.Results;
xyPixelSize = pr.xyPixelSize;
dz = pr.dz;
dr = pr.dr;
dtheta = pr.dtheta;
resThreshMethod = pr.resThreshMethod;
resThresh = pr.resThresh;
N = pr.N;
bbox = pr.bbox;
resAxis = pr.resAxis;
skipConeRegion = pr.skipConeRegion;
clipPer = pr.clipPer;
debug = pr.debug;

if exist(fnout, 'file')
    fprintf('Result file %s already exists!\n', fnout);
    return;
end

% run analysis
if ~isempty(bbox)
    centers = round((bbox(:, 4 : 6) + bbox(:, 1 : 3)) / 2);
    n = size(centers, 1);
else
    centers = getImageSize(fn);
    n = 1;    
end

fsc_mu = cell(n, 1);
res_mu = cell(n, 1);
fsc_cell = cell(n, 1);
res_cell = cell(n, 1);

for i = 1 : n
    if isempty(bbox)
        bbox_i = bbox;
    else
        bbox_i = bbox(i, :);
    end
    
    [fsc_mu_i, res_mu_i, fsc_cell_i, res_cell_i] = XR_one_image_FSC_analysis(fn, ...
        'xyPixelSize', xyPixelSize, 'dz', dz, 'dr', dr, 'dtheta', dtheta, 'N', N, ...
        'resThreshMethod', resThreshMethod, 'resThresh', resThresh, 'resAxis', resAxis, ...
        'skipConeRegion', skipConeRegion, 'bbox', bbox_i, 'clipPer', clipPer);

    fsc_mu{i} = fsc_mu_i;
    res_mu{i} = res_mu_i;
    fsc_cell{i} = fsc_cell_i;
    res_cell{i} = res_cell_i;    
end

if n == 1
    fsc_mu = fsc_mu{1};
    res_mu = res_mu{1};
    fsc_cell = fsc_cell{1};
    res_cell = res_cell{1};
end

save('-v7.3', fnout, 'fsc_mu', 'res_mu', 'fsc_cell', 'res_cell', 'centers', 'bbox');

end
            
            