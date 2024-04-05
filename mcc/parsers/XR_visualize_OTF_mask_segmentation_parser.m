function [] = XR_visualize_OTF_mask_segmentation_parser(psfFn, OTFCumThresh, skewed)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('psfFn', @ischar);
ip.addOptional('OTFCumThresh', 0.85, @(x) isscalar(x) || ischar(x));
ip.addOptional('skewed', [], @(x) isempty(x) || islogical(x) || ischar(x));

ip.parse(psfFn, OTFCumThresh, skewed);

pr = ip.Results;

if ischar(OTFCumThresh)
    OTFCumThresh = str2num(OTFCumThresh);
end
if ischar(skewed)
    skewed = str2num(skewed);
end

XR_visualize_OTF_mask_segmentation(psfFn, OTFCumThresh, skewed);

end

