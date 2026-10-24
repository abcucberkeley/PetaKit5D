function [psf_b] = omw_backprojector_files(psf, psfgenPath, psfFsn, wienerAlpha, OTFCumThresh, skewed, hannWinBounds, PSFfn, saveFigure)
% load the OMW back projector of a PSF from the psfgen folder, or generate it
% and save it there together with the OTF, the OTF mask and, if Java is
% available, a figure of the OTF mask outline.
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
% Moved here from RLdecon so that XR_decon_data_wrapper can make the back
% projector once before submitting the per-file jobs.


if nargin < 9
    saveFigure = true;
end

bpFn = sprintf('%s/%s_back_projector_alpha_%0.6f_otf_cum_thresh_%0.6f.tif', psfgenPath, psfFsn, wienerAlpha, OTFCumThresh);
if exist(bpFn, 'file')
    fprintf('Load existing OMW back projector %s for %s ...\n', bpFn, PSFfn);
    % kept in memory for the next batches or files of the session
    psf_b = decon_readtiff_cached(bpFn);
    return;
end

fprintf('OMW back projector generation for %s ...\n', PSFfn);
loc_uuid = get_uuid();
bpTmpFn = sprintf('%s/%s_back_projector_alpha_%0.6f_otf_cum_thresh_%0.6f_%s.tif', psfgenPath, psfFsn, wienerAlpha, OTFCumThresh, loc_uuid);
[psf_b, ~, abs_OTF_c, OTF_mask] = omw_backprojector_generation(psf, wienerAlpha, skewed, 'OTFCumThresh', OTFCumThresh, 'hannWinBounds', hannWinBounds);

if saveFigure && usejava('jvm')
    visible = false;
    fig = visualize_OTF_and_mask_outline(abs_OTF_c, OTF_mask, visible);
    figFn = sprintf('%s/%s_back_projector_alpha_%0.6f_otf_cum_thresh_%0.6f_figure.png', psfgenPath, psfFsn, wienerAlpha, OTFCumThresh);
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
