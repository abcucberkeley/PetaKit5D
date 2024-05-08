function [data, bim] = readzarr(filepath, options)
% wrapper for zarr reader
% bbox provide the region bounding box [ymin, xmin, zmin, ymax, xmax, zmax]
% 
% Author: Xiongtao Ruan (01/25/2022)
% 
% xruan (02/01/2022): add support for parallel read zarr


arguments
    filepath char 
    options.inputBbox (1, :) {mustBeNumeric} = []
    options.sparseData (1, :) {mustBeNumericOrLogical} = true
end

inputBbox = options.inputBbox;
sparseData = options.sparseData;

try 
    if isempty(inputBbox)
        data = parallelReadZarr(filepath);
    else
        inputBbox = inputBbox(:)';
        data = parallelReadZarr(filepath, 'bbox', inputBbox, 'sparse', sparseData); 
    end
    
    if nargout == 2
        bim = blockedImage(filepath, "Adapter", CZarrAdapter);
    end
catch ME
    if ~exist(filepath, 'dir')
        error('zarr file %s does not exist!', filepath);
    end
    pe = pyenv;
    if pe.Status ~= "Loaded"
        throw(ME);
    end

    disp(ME);
    disp('Use the alternative zarr reader (ZarrAdapter)...');   
    
    bim = blockedImage(filepath, "Adapter", ZarrAdapter);    
    if isempty(inputBbox)
        data = bim.Adapter.getIORegion([1, 1, 1], bim.Size);    
    else
        data = bim.Adapter.getIORegion(inputBbox(1 : 3), inputBbox(4 : 6));
    end
end

end
