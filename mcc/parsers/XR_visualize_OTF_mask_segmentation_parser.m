function [] = XR_visualize_OTF_mask_segmentation_parser(psfFn, OTFCumThresh, skewed, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('psfFn', @ischar);
ip.addOptional('OTFCumThresh', 0.85, @(x) isscalar(x) || ischar(x));
ip.addOptional('skewed', [], @(x) isempty(x) || islogical(x) || ischar(x));
ip.addParameter('minIntThrsh', 1e-3, @(x) isnumeric(x) || ischar(x));
ip.addParameter('visible',  true, @(x) islogical(x) || ischar(x));
ip.addParameter('saveFig',  false, @(x) islogical(x) || ischar(x));

ip.parse(psfFn, OTFCumThresh, skewed, varargin{:});

pr = ip.Results;
minIntThrsh = pr.minIntThrsh;
visible = pr.visible;
saveFig = pr.saveFig;

if ischar(OTFCumThresh)
    OTFCumThresh = str2num(OTFCumThresh);
end
if ischar(skewed)
    skewed = str2num(skewed);
end
if ischar(minIntThrsh)
    minIntThrsh = str2num(minIntThrsh);
end
if ischar(visible)
    visible = str2num(visible);
end
if ischar(saveFig)
    saveFig = str2num(saveFig);
end

XR_visualize_OTF_mask_segmentation(psfFn, OTFCumThresh, skewed, minIntThrsh=minIntThrsh, ...
    visible=visible, saveFig=saveFig);

end

