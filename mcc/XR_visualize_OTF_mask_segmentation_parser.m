function [] = XR_visualize_OTF_mask_segmentation_parser(psfFn, OTFCumThresh, skewed)

%#function XR_visualize_OTF_mask_segmentation

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('psfFn', @ischar);
ip.addRequired('OTFCumThresh', @(x) isnumeric(x) || ischar(x));
ip.addRequired('skewed', @(x) isempty(x) || isnumeric(x) || ischar(x));

ip.parse(psfFn, OTFCumThresh, skewed);

pr = ip.Results;
psfFn = pr.psfFn;
OTFCumThresh = pr.OTFCumThresh;
skewed = pr.skewed;

if ischar(OTFCumThresh)
    OTFCumThresh = str2num(OTFCumThresh);
end

if ischar(skewed)
    if strcmp(skewed, '[]')
        skewed = [];
    else
        skewed = strcmp(skewed, 'true');
    end
end

XR_visualize_OTF_mask_segmentation(psfFn, OTFCumThresh, skewed);

end
