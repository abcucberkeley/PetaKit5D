function [psf_b, fsns] = omw_backprojector_files(psf, psfgenPath, psfFsn, wienerAlpha, OTFCumThresh, skewed, hannWinBounds, PSFfn, saveFigure)
% load the OMW back projector of a PSF from the psfgen folder, or generate it
% and save it there together with the OTF, the OTF mask and, if Java is
% available, a figure of the OTF mask outline.
%
% The back projector's file name holds wienerAlpha, OTFCumThresh and, unless
% they are the default [0.8, 1], hannWinBounds, each with 6 decimals (back
% projectors of the default bounds keep the name they had).
%
% Required inputs:
%                psf : 3d array. The (generated) PSF the back projector is made from.
%         psfgenPath : char. The psfgen folder of the deconvolution results.
%             psfFsn : char. PSF file name without extension, used in the file names.
%        wienerAlpha : a number. Wiener parameter.
%       OTFCumThresh : a number. Accumulated energy quantile to segment OTF mask.
%             skewed : empty or true|false. PSF in the skewed space, empty for automatic.
%      hannWinBounds : 1x2 vector. Apodization range for distance matrix with Hann window function.
%              PSFfn : char. PSF file path, only used in the messages.
%
% Optional input:
%         saveFigure : true|false (default: true). Save the figure of the OTF mask outline (needs Java).
%
% Outputs:
%              psf_b : the back projector.
%               fsns : names of the back projector, OTF and OTF mask files in psfgenPath.
%
% Moved here from RLdecon so that XR_decon_data_wrapper can make the back
% projector once before submitting the per-file jobs.


if nargin < 9
    saveFigure = true;
end

bpStr = sprintf('alpha_%0.6f_otf_cum_thresh_%0.6f', wienerAlpha, OTFCumThresh);
boundStr = sprintf('_hann_win_bounds_%0.6f_%0.6f', hannWinBounds(1), hannWinBounds(2));
if ~strcmp(boundStr, '_hann_win_bounds_0.800000_1.000000')
    bpStr = [bpStr, boundStr];
end
fsns = {sprintf('%s_back_projector_%s.tif', psfFsn, bpStr), sprintf('%s_OTF.tif', psfFsn), ...
    sprintf('%s_OTF_mask_otf_cum_thresh_%0.6f.tif', psfFsn, OTFCumThresh)};

bpFn = sprintf('%s/%s', psfgenPath, fsns{1});
if exist(bpFn, 'file')
    fprintf('Load existing OMW back projector %s for %s ...\n', bpFn, PSFfn);
    % kept in memory for the next batches or files of the session
    psf_b = decon_readtiff_cached(bpFn);
    return;
end

fprintf('OMW back projector generation for %s ...\n', PSFfn);
loc_uuid = get_uuid();
bpTmpFn = sprintf('%s/%s_back_projector_%s_%s.tif', psfgenPath, psfFsn, bpStr, loc_uuid);
[psf_b, ~, abs_OTF_c, OTF_mask] = omw_backprojector_generation(psf, wienerAlpha, skewed, 'OTFCumThresh', OTFCumThresh, 'hannWinBounds', hannWinBounds);

if saveFigure && usejava('jvm')
    visible = false;
    fig = visualize_OTF_and_mask_outline(abs_OTF_c, OTF_mask, visible);
    figFn = sprintf('%s/%s_back_projector_%s_figure.png', psfgenPath, psfFsn, bpStr);
    print(fig, figFn, '-dpng', '-r0');
    close(fig);
end

% write back projector
writetiff(psf_b, bpTmpFn);
movefile(bpTmpFn, bpFn);

% write OTF
otfFn = sprintf('%s/%s_OTF.tif', psfgenPath, psfFsn);
otfTmpFn = sprintf('%s/%s_OTF_%s.tif', psfgenPath, psfFsn, loc_uuid);
writetiff(abs_OTF_c, otfTmpFn);
movefile(otfTmpFn, otfFn);

% write OTF mask
maskFn = sprintf('%s/%s_OTF_mask_otf_cum_thresh_%0.6f.tif', psfgenPath, psfFsn, OTFCumThresh);
maskTmpFn = sprintf('%s/%s_OTF_mask_otf_cum_thresh_%0.6f_%s.tif', psfgenPath, psfFsn, OTFCumThresh, loc_uuid);
writetiff(uint8(OTF_mask), maskTmpFn);
movefile(maskTmpFn, maskFn);

end
