function [] = XR_resampleFrame_parser(fn, fnout, rsfactor, varargin)
% resample image for a single frame
%
% 
% Author: Xiongtao Ruan (01/15/2021)
% 
% xruan (05/11/2022): add support for 2d image

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('fn', @ischar);
ip.addRequired('fnout', @ischar);
ip.addRequired('rsfactor', @(x) isnumeric(x) || ischar(x));
ip.addParameter('Interp', 'linear', @ischar);
ip.addParameter('Save16bit', true ,@(x) islogical(x) || ischar(x)); % saves 16bit, else single
ip.addParameter('zarrFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('saveZarr', false, @(x) islogical(x) || ischar(x)); % use zarr file as output
ip.parse(fn, fnout, rsfactor, varargin{:});

ip.parse(fn, fnout, rsfactor, varargin{:});
pr = ip.Results;

Interp = pr.Interp;
Save16bit = pr.Save16bit;
zarrFile = pr.zarrFile;
saveZarr = pr.saveZarr;

if ischar(rsfactor)
    rsfactor = str2num(rsfactor);
end
if ischar(Save16bit)
    Save16bit = strcmp(Save16bit,'true');
end
if ischar(zarrFile)
    zarrFile = strcmp(zarrFile,'true');
end
if ischar(saveZarr)
    saveZarr = strcmp(saveZarr,'true');
end

XR_resampleFrame(fn, fnout, rsfactor, 'Interp', Interp, 'Save16bit', Save16bit, ...
    'zarrFile', zarrFile, 'saveZarr', saveZarr)
