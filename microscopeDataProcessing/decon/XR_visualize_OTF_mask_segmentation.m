function [] = XR_visualize_OTF_mask_segmentation(psfFn, OTFCumThresh, skewed, varargin)
% visualze OTF mask outline segmented with given threshold overlay with OTF
% 
% Author: Xiongtao Ruan (11/13/2023)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('psfFn', @ischar);
ip.addOptional('OTFCumThresh', 0.85, @(x) isscalar(x));
ip.addOptional('skewed', [], @(x) isempty(x) || islogical(x));
ip.addParameter('minIntThrsh', 1e-3, @(x) isnumeric(x));
ip.addParameter('visible',  true, @(x) islogical(x));
ip.addParameter('saveFig',  false, @(x) islogical(x));

ip.parse(psfFn, OTFCumThresh, skewed, varargin{:});

pr = ip.Results;
OTFCumThresh = pr.OTFCumThresh;
skewed = pr.skewed;
minIntThrsh = pr.minIntThrsh;
visible = pr.visible;
saveFig = pr.saveFig;

% load PSF and clear up the PSF
if ~exist(psfFn, 'file') && ~exist(psfFn, 'dir')
    error('The PSF file %s does not exist, please double check to make sure the path is correct!\n', psfFn);
end
[~, ~, ext] = fileparts(psfFn);
switch ext
    case {'.tif', '.tiff'}
        psf = single(readtiff(psfFn));
    case '.zarr'
        psf = single(readzarr(psfFn));
    otherwise
        error('Unknown format for PSF file %s\n!', psfFn);
end
dz_data = 0.5;
dz_psf = 0.5;
medFactor = 1.5;
PSFGenMethod = 'masked';
psf = psf_gen_new(psf, dz_psf, dz_data, medFactor, PSFGenMethod);
psf = psf ./ sum(psf, 'all');

% generate OTF mask and OMW backward projector
fprintf('OTFCumThresh: %f, skewed: %s\n', OTFCumThresh, string(skewed));
alpha = 0.01;
hanWinBounds = [0.8, 1.0];
[b_omw, OTF_bp_omw, abs_OTF_c, OTF_mask] = omw_backprojector_generation(psf, alpha, skewed, ...
    'OTFCumThresh', OTFCumThresh, 'hanWinBounds', hanWinBounds, minIntThrsh=minIntThrsh);

% visualize OTF mask on top of OTF 
fig = visualize_OTF_and_mask_outline(abs_OTF_c, OTF_mask, visible);

title_str = sprintf('OTFCumThresh: %f, skewed: %s', OTFCumThresh, string(skewed));
 % - Build title axes and title.
axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
text( 0.5, 0, title_str, 'FontSize', 14', 'FontWeight', 'Bold', ...
  'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;

if saveFig
    [pstr, fsn, ext] = fileparts(psfFn);
    figPath = sprintf('%s/figures/', pstr);
    if ~exist(figPath, 'dir')
        mkdir(figPath);
    end
    figFullname = sprintf('%s/%s.png', figPath, fsn);
    print(fig, '-painters','-dpng', '-loose', figFullname);
end

if ~visible
    close(fig);
end

end

