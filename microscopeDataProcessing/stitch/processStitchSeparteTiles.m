function [] = processStitchSeparteTiles(zarrFullpaths, stitchPath, stIndices, imSizes, pImSz, varargin)
% process the tiles in the corresponding location in the stitched image as
% separte images. 
% 
% 
% Author: Xiongtao Ruan (11/19/2020)
% distance map for the primary channel.

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('zarrFullpaths', @iscell);
ip.addRequired('stitchPath', @(x) ischar(x));
ip.addRequired('stIndices', @isnumeric);
ip.addRequired('imSizes', @isnumeric);
ip.addRequired('pImSz', @isnumeric);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('convertToTiff', true, @islogical);
ip.addParameter('SaveMIP', true, @islogical);

ip.parse(zarrFullpaths, stitchPath, stIndices, imSizes, pImSz, varargin{:});

Overwrite = ip.Results.Overwrite;
convertToTiff = ip.Results.convertToTiff;
SaveMIP = ip.Results.SaveMIP;

uuid = get_uuid();
nys = pImSz(1);
nxs = pImSz(2);
nzs = pImSz(3);

nF = numel(zarrFullpaths);
% incase of error of saving tiff after calling dask
writetiff(uint8(rand(10) > 0.5), sprintf('/tmp/%s.tif', uuid));
for i = 1 : nF
    zarrFullpath = zarrFullpaths{i};
    [~, fsname, ext] = fileparts(zarrFullpath);
    OutputFullname = [stitchPath, filesep, fsname, ext];
    OutputTempname = [stitchPath, filesep, fsname, uuid, ext];
    
    if exist(OutputFullname, 'dir')
        continue;
    end
    
    fprintf('Process tile %d : %s ...\n', i, fsname);
    st_idx = stIndices(i, :);

    sx = imSizes(i, 2);
    sy = imSizes(i, 1);
    sz = imSizes(i, 3);
    xridx = max(1, 1 - st_idx(1)) : min(sx, nxs - st_idx(1));
    yridx = max(1, 1 - st_idx(2)) : min(sy, nys - st_idx(2));
    zridx = max(1, 1 - st_idx(3)) : min(sz, nzs - st_idx(3));
        
    xidx = st_idx(1) + xridx;
    yidx = st_idx(2) + yridx;
    zidx = st_idx(3) + zridx;
    
    pad_width = [yidx(1) - 1, nys - yidx(end); xidx(1) - 1, nxs - xidx(end); zidx(1) - 1, nzs - zidx(end)];
    pad_width = cast(pad_width, 'int32');
    pad_width_cell = mat2cell(pad_width', 2, [1, 1, 1]);

    pypad_width = {py.tuple(pad_width_cell{1}), py.tuple(pad_width_cell{2}), py.tuple(pad_width_cell{3})};
    pypad_width = py.tuple(pypad_width);
    py.daskAPI.daskZarrPadArray(zarrFullpath, OutputTempname, pypad_width);
    
    if ~exist(OutputFullname, 'dir')
        movefile(OutputTempname, OutputFullname);
    else
        rmdir(OutputFullname, 's');
    end
    
    if convertToTiff
        bim = blockedImage(OutputFullname, 'Adapter', ZarrAdapter);
        tiffFullname = [stitchPath, filesep, fsname, '.tif'];
        tiffTempname = [stitchPath, filesep, fsname, uuid, '.tif'];
        write(bim, tiffTempname, 'Adapter', MPageTiffAdapter, 'BlockSize', [bim.Size(1), bim.Size(2), 5]);
        movefile(tiffTempname, tiffFullname);
    end
    
    if SaveMIP
        stcMIPPath = sprintf('%s/MIPs/', stitchPath);
        if ~exist(stcMIPPath, 'dir')
            mkdir(stcMIPPath);
            fileattrib(stcMIPPath, '+w', 'g');
        end
        stcMIPname = sprintf('%s%s_MIP_z.tif', stcMIPPath, fsname);
        saveMIP_zarr(OutputFullname, stcMIPname);
    end
    fprintf('Done!\n');
end

end
