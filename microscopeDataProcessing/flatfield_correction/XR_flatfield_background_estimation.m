function [] = XR_flatfield_background_estimation(frameFullnames, resultFullname, varargin)
% flat-field background estimation for scatter light with variant
% background in the field
% 
% Author: Xiongtao Ruan (09/24/2021)
%
% xruan (07/13/2022): change threshold of removing slices with intensity to 0.1% (rather than 0).
% xruan (08/01/2023): add options to save MIPs (MIP z) for ff estimation


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('frameFullnames'); 
ip.addRequired('resultFullname'); 
ip.addParameter('objectiveScan', false, @islogical);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('MovieSelector', 'cell', @ischar);
ip.addParameter('sigmaThresh', 8, @isnumeric); % camera sigma threshold for background image, for sCMOS, std should be around 4, set a higher one
ip.addParameter('sigmaFactor', 4, @isnumeric); % std factor to decide background image
ip.addParameter('probThresh', 0.001, @isnumeric); % probability threshold for high count images, 1% for 4sigma
ip.addParameter('save16bit', true , @islogical); % saves deskewed data as 16 bit -- not for quantification
ip.addParameter('saveMIP', false, @islogical); % save MIPs
ip.addParameter('uuid', '', @ischar);

ip.parse(frameFullnames, resultFullname, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
sigmaThresh = pr.sigmaThresh;
sigmaFactor = pr.sigmaFactor;
probThresh = pr.probThresh;
saveMIP = pr.saveMIP;

uuid = pr.uuid;
if isempty(uuid)
    uuid = get_uuid();
end

if exist(resultFullname, 'file') && ~Overwrite
    fprintf('The result %s already exists, skip it!\n', resultFullname);
    return;
end

im = cell(numel(frameFullnames), 1);
for i = 1 : numel(frameFullnames)
    try
        im{i} = parallelReadTiff(frameFullnames{i});
    catch ME
        im{i} = readtiff(frameFullnames{i});
    end
end
im = single(cat(3, im{:}));
% numSlice = size(source, 3);

if saveMIP
    MIP = max(im, [], 3);
    
    [outPath, fsn] = fileparts(resultFullname);
    MIPFullname = sprintf('%s/MIPs/%s_MIP_z.tif', outPath, fsn);
    MIPTmpname = sprintf('%s/MIPs/%s_MIP_z_%s.tif', outPath, fsn, uuid);
    writetiff(MIP, MIPTmpname);
    movefile(MIPTmpname, MIPFullname)
end

mu = mean(im, 3);
sigma_z = squeeze(std(im, 0, [1, 2]));
im = im(:, :, sigma_z < min(sigmaThresh * 5, median(sigma_z)));
% clear source;

% im_1 = single(im);
mu = mean(im, 3);
sigma = std(im, 0, 3);
sigma_z = squeeze(std(im - mu, 0, [1, 2]));
sigma_z_orig = sigma_z;

% removing slices with intensity
nIter = 100;
for i = 1 : nIter
    high_count = im > mu + sigmaFactor * sigma;
    inds_remove = squeeze(mean(high_count, [1, 2])) > probThresh ...
        & (sigma_z > min(sigmaThresh, median(sigma_z_orig)) ...
        |  squeeze(max(im - mu, [], [1, 2]) - min(im - mu,  [], [1, 2])) > sigma_z * 20);
    
    if mean(inds_remove) < 0.001
        break;
    end
    
    im = im(:, :, ~inds_remove);
    mu = mean(im, 3);
    sigma = std(im, 0, 3);
    sigma_z = squeeze(std(im - mu, 0, [1, 2]));
end

med = median(im, 3);
numSlice = size(im, 3);

% compute percentiles for each pixel in xy plane
pmat = 0 : 5 : 100;
imp = [];
if numSlice > 200
    imp = prctile(im, pmat, 3);
end

[pathstr, fname, ext] = fileparts(resultFullname);
resultTempname = [pathstr, fname, '_', uuid, ext];
% save('-v7.3', resultTempname, 'im', 'mu', 'sigma', 'med', 'numSlice');
save('-v7.3', resultTempname, 'imp', 'mu', 'sigma', 'med', 'numSlice', 'pmat');
movefile(resultTempname, resultFullname);

end

