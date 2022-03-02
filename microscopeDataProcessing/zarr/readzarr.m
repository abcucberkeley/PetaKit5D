function [data, bim] = readzarr(filepath, varargin)
% wrapper for zarr reader
% bbox provide the region bounding box [ymin, xmin, zmin, ymax, xmax, zmax]
% 
% Author: Xiongtao Ruan (01/25/2022)
% 
% xruan (02/01/2022): add support for parallel read zarr


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('filepath', @(x) ischar(x));
ip.addParameter('bbox', [], @isnumeric);
ip.parse(filepath, varargin{:});

pr = ip.Results;
bbox = pr.bbox;

try 
    if isempty(bbox)
        data = parallelReadZarr(filepath);
    else
        data = parallelReadZarr(filepath, bbox);        
    end
catch ME
    disp(ME);
    bim = blockedImage(filepath, "Adapter", ZarrAdapter);    
    if isempty(bbox)
        data = bim.Adapter.getIORegion([1, 1, 1], bim.Size);    
    else
        data = bim.Adapter.getIORegion(bbox(1 : 3), bbox(4 : 6));
    end
end

end