function [] = XR_psf_analysis_plot_parser(frameFullname, figureFullname, RW_info_Fullname, ...
    ch_ind, source_descrip, xypixsize, zpixsize, NAdet, index, exc_lambda, det_lambda, PSFsubpix, gamma, bgFactor)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('frameFullname', @ischar);
ip.addRequired('figureFullname', @ischar);
ip.addRequired('RW_info_Fullname', @ischar);
ip.addRequired('ch_ind', @(x) isscalar(x) || ischar(x));
ip.addRequired('source_descrip', @ischar);
ip.addRequired('xypixsize', @(x) isscalar(x) || ischar(x));
ip.addRequired('zpixsize', @(x) isscalar(x) || ischar(x));
ip.addRequired('NAdet', @(x) isscalar(x) || ischar(x));
ip.addRequired('index', @(x) isscalar(x) || ischar(x));
ip.addRequired('exc_lambda', @(x) isscalar(x) || ischar(x));
ip.addRequired('det_lambda', @(x) isscalar(x) || ischar(x));
ip.addRequired('PSFsubpix', @ischar);
ip.addRequired('gamma', @(x) isscalar(x) || ischar(x));
ip.addRequired('bgFactor', @(x) isscalar(x) || ischar(x));

ip.parse(frameFullname, figureFullname, RW_info_Fullname, ch_ind, source_descrip, ...
    xypixsize, zpixsize, NAdet, index, exc_lambda, det_lambda, PSFsubpix, gamma, bgFactor);

pr = ip.Results;

if ischar(ch_ind)
    ch_ind = str2num(ch_ind);
end
if ischar(xypixsize)
    xypixsize = str2num(xypixsize);
end
if ischar(zpixsize)
    zpixsize = str2num(zpixsize);
end
if ischar(NAdet)
    NAdet = str2num(NAdet);
end
if ischar(index)
    index = str2num(index);
end
if ischar(exc_lambda)
    exc_lambda = str2num(exc_lambda);
end
if ischar(det_lambda)
    det_lambda = str2num(det_lambda);
end
if ischar(gamma)
    gamma = str2num(gamma);
end
if ischar(bgFactor)
    bgFactor = str2num(bgFactor);
end

XR_psf_analysis_plot(frameFullname, figureFullname, RW_info_Fullname, ch_ind, ...
    source_descrip, xypixsize, zpixsize, NAdet, index, exc_lambda, det_lambda, ...
    PSFsubpix, gamma, bgFactor);

end

