function mip_masks_b10_ml(folder, what, threads, repeats)
% MIP mask benchmark (B10), the MATLAB side, with maxNumCompThreads set to each thread count:
% 'masks': XR_generate_MIP_masks of folder/mips (defaults, parseCluster false), the masks
% folder removed before each run (untimed); 'seg': XR_visualize_OTF_mask_segmentation of
% folder/psf/560_c.tif (0.88, skewed, saved, not visible), as the demo calls it (needs the JVM).

nt0 = maxNumCompThreads;
mips = [folder, '/mips'];
out = [mips, '/masks'];
psf = [folder, '/psf/560_c.tif'];
figs = [folder, '/psf/figures'];
for t = threads
    maxNumCompThreads(t);
    times = zeros(1, repeats);
    for r = 1 : repeats
        if strcmp(what, 'masks')
            if exist(out, 'dir')
                rmdir(out, 's');
            end
            t0 = tic;
            XR_generate_MIP_masks(mips, 'parseCluster', false);
        else
            if exist(figs, 'dir')
                rmdir(figs, 's');
            end
            t0 = tic;
            XR_visualize_OTF_mask_segmentation(psf, 0.88, true, 'visible', false, 'saveFig', true);
        end
        times(r) = toc(t0);
    end
    if strcmp(what, 'masks')
        pk5d_bench_record('generate_mip_masks', t, times);
    else
        pk5d_bench_record('otf_mask_segmentation', t, times);
    end
end
maxNumCompThreads(nt0);
end
