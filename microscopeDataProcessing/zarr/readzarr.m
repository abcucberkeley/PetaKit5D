function [data, bim] = readzarr(filepath, options)
% wrapper for zarr reader
% bbox provide the region bounding box [ymin, xmin, zmin, ymax, xmax, zmax]
% 
% Author: Xiongtao Ruan (01/25/2022)
% 
% xruan (02/01/2022): add support for parallel read zarr


arguments
    filepath char 
    options.bbox (1, :) {mustBeNumeric} = []
end

bbox = options.bbox;

try 
    if isempty(bbox)
        data = parallelReadZarr(filepath);
    else
        bbox = bbox(:)';
        data = parallelReadZarr(filepath, bbox);        
    end
    if nargout == 2
        bim = blockedImage(filepath, "Adapter", CZarrAdapter);
    end
catch ME
    disp(ME);
    disp('Use the alternative zarr reader (ZarrAdapter)...');   
    
    bim = blockedImage(filepath, "Adapter", ZarrAdapter);    
    if isempty(bbox)
        data = bim.Adapter.getIORegion([1, 1, 1], bim.Size);    
    else
        data = bim.Adapter.getIORegion(bbox(1 : 3), bbox(4 : 6));
    end
end

end
