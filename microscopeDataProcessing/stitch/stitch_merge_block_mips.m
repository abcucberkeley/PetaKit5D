function [mip] = stitch_merge_block_mips(flagFullpaths, imSize, dtype)
% merge the z MIPs that the blending tasks (processStitchBlock with saveMIP
% true) save in their flag files into the z MIP of the stitched image, so
% that the image doesn't have to be read again. Returns [] if a flag has no
% MIPs (e.g., one from a run without them), then the caller computes the MIP
% from the image.
%
% flagFullpaths: flag files of all blending tasks
% imSize: size of the stitched image
% dtype: data type of the stitched image


% the smallest value of dtype, or NaN, which max ignores
if any(strcmp(dtype, {'single', 'double'}))
    mip = nan(imSize(1 : 2), dtype);
else
    mip = repmat(intmin(dtype), imSize(1 : 2));
end

for t = 1 : numel(flagFullpaths)
    try
        a = load(flagFullpaths{t});
    catch
        a = struct();
    end
    if ~isfield(a, 'mips') || ~isfield(a, 'mipBboxes')
        mip = [];
        return;
    end
    for b = 1 : numel(a.mips)
        bb = a.mipBboxes(b, :);
        mip(bb(1) : bb(3), bb(2) : bb(4)) = max(mip(bb(1) : bb(3), bb(2) : bb(4)), a.mips{b});
    end
end

end
