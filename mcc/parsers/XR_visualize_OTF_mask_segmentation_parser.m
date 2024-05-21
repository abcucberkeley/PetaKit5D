function [] = XR_visualize_OTF_mask_segmentation_parser(psfFn, OTFCumThresh, skewed, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('psfFn', @ischar);
ip.addOptional('OTFCumThresh', 0.85, @(x) isscalar(x) || ischar(x));
ip.addOptional('skewed', [], @(x) isempty(x) || islogical(x) || ischar(x));
ip.addParameter('minIntThrsh', 2.5e-3, @(x) isnumeric(x) || ischar(x));

ip.parse(psfFn, OTFCumThresh, skewed, varargin{:});

pr = ip.Results;
minIntThrsh = pr.minIntThrsh;

if ischar(OTFCumThresh)
    OTFCumThresh = str2num(OTFCumThresh);
end
if ischar(skewed)
    skewed = str2num(skewed);
end
if ischar(minIntThrsh)
    minIntThrsh = str2num(minIntThrsh);
end

XR_visualize_OTF_mask_segmentation(psfFn, OTFCumThresh, skewed, minIntThrsh=minIntThrsh);

end

