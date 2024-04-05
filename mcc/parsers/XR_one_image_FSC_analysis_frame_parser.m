function [fsc_mu, res_mu, fsc_cell, res_cell] = XR_one_image_FSC_analysis_frame_parser(fn, fnout, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('fn', @(x) ischar(x) || isnumeric(x));
ip.addParameter('xyPixelSize', 0.108, @(x) isnumeric(x) || ischar(x));
ip.addParameter('dz', 0.108, @(x) isnumeric(x) || ischar(x));
ip.addParameter('dr', 1 , @(x) isnumeric(x) || ischar(x));
ip.addParameter('dtheta', pi / 12 , @(x) isnumeric(x) || ischar(x));
ip.addParameter('resThreshMethod', 'fixed', @ischar);
ip.addParameter('resThresh', 0.2, @(x) isnumeric(x) || ischar(x));
ip.addParameter('N', [501, 501, 501], @(x) isnumeric(x) || ischar(x));
ip.addParameter('bbox', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('resAxis', 'xz', @ischar);
ip.addParameter('skipConeRegion', true, @(x) islogical(x) || ischar(x));
ip.addParameter('clipPer', [], @(x) isnumeric(x) || ischar(x)); % clip intensity higher than the given percentile
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));

ip.parse(fn,  varargin{:});

pr = ip.Results;
xyPixelSize = pr.xyPixelSize;
dz = pr.dz;
dr = pr.dr;
dtheta = pr.dtheta;
resThreshMethod = pr.resThreshMethod;
resThresh = pr.resThresh;
N = pr.N;
bbox = pr.bbox;
resAxis = pr.resAxis;
skipConeRegion = pr.skipConeRegion;
clipPer = pr.clipPer;
debug = pr.debug;

if ischar(xyPixelSize)
    xyPixelSize = str2num(xyPixelSize);
end
if ischar(dz)
    dz = str2num(dz);
end
if ischar(dr)
    dr = str2num(dr);
end
if ischar(dtheta)
    dtheta = str2num(dtheta);
end
if ischar(resThresh)
    resThresh = str2num(resThresh);
end
if ischar(N)
    N = str2num(N);
end
if ischar(bbox)
    bbox = str2num(bbox);
end
if ischar(skipConeRegion)
    skipConeRegion = str2num(skipConeRegion);
end
if ischar(clipPer)
    clipPer = str2num(clipPer);
end
if ischar(debug)
    debug = str2num(debug);
end

XR_one_image_FSC_analysis_frame(fn, xyPixelSize=xyPixelSize, dz=dz, dr=dr, ...
    dtheta=dtheta, resThreshMethod=resThreshMethod, resThresh=resThresh, N=N, ...
    bbox=bbox, resAxis=resAxis, skipConeRegion=skipConeRegion, clipPer=clipPer, ...
    debug=debug);

end

