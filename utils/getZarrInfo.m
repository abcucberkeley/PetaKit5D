function [zInfo] = getZarrInfo(filePath)
% Get the info of .zarray
% 
% 
% Author: Xiongtao Ruan (06/2024)


warning ('off','all');

zInfo = struct();
if strcmp(filePath(end - 3 : end), 'zarr') || exist([filePath, '/.zarray'], 'file')
    zarrayFn = [filePath, '/.zarray'];
    fid = fopen(zarrayFn);
    raw = fread(fid);
    fclose(fid);
    str = char(raw');
    jdata = jsondecode(str);

    zInfo.size = jdata.shape(:)';
    zInfo.blockSize = jdata.chunks(:)';
    zInfo.dtype = jdata.dtype;
    zInfo.order = jdata.order;
    zInfo.compressor = jdata.compressor;
    zInfo.fill_value = jdata.fill_value;
    zInfo.filters = jdata.filters;
    zInfo.zarr_format = jdata.zarr_format;
else
    error('The input zarr file %s does not exist or the .zarray file is missing!', filePath)
end

end

