function [] = zarrToTiff(zarrFilename, tiffFilename, varargin)
% use image block function to convert zarr to tiff
% adapted from tiffToZarr.m
% 
% Author: Xiongtao Ruan (01/19/2021)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('tifFilename', @ischar);
ip.addRequired('zarrFilename', @ischar);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('usrFcn', '', @(x) isempty(x) || isa(x,'function_handle') || ischar(x));
ip.addParameter('uuid', '', @isstr);

ip.parse(zarrFilename, tiffFilename, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
usrFcn = pr.usrFcn;
uuid = pr.uuid;

if exist(tiffFilename, 'file')
    if ~Overwrite
        fprintf('Tiff result %s already exists, skip it!\n', tiffFilename);
        return;
    else
        delete(tiffFilename);
    end
end

[pthstr, fsname] = fileparts(tiffFilename);
if ~exist(pthstr, 'dir')
    mkdir_recursive(pthstr);
end

if isempty(uuid)
    uuid = get_uuid();
end

fprintf('Convert Zarr to Tiff for %s ...\n', zarrFilename);    

% load tiff file as image block file
% Direct conversion should be faster. 
try 
    bim = blockedImage(zarrFilename, 'Adapter', CZarrAdapter);
catch ME
    disp(ME)
    bim = blockedImage(zarrFilename, 'Adapter', ZarrAdapter);    
end
sz = bim.Size;
if any(sz == 0)
    fprintf('The input image is empty, skip it!\n');
    return;
end    

if ~isempty(usrFcn)
    if ischar(usrFcn)
        usrFcn = str2func(usrFcn);
    end
    bim = apply(bim, @(bs) usrFcn(bs.Data), 'BlockSize', sz);
end

tmpFilename = sprintf('%s_%s.tif', tiffFilename(1 : end - 4), uuid);
% in case of too large block size, we reduce block size to 5
bsz_z = floor(10 * 1024^3 / prod(sz(1 : 2)) / 16);
try
    % BlockSize = [sz(1 : 2), min(sz(3), bsz_z)];
    % write(bim, tmpFilename, "BlockSize", BlockSize, "Adapter", MPageTiffAdapter);
    writetiff(bim.Adapter.getIORegion([1, 1, 1], bim.Size), tmpFilename);
catch ME
    disp(ME)
    BlockSize = [bim.BlockSize(1 : 2), min(sz(3), bsz_z)];
    write(bim, tmpFilename, "BlockSize", BlockSize, "Adapter", MPageTiffAdapter);
end

% mv tmp result folder to output folder
movefile(tmpFilename, tiffFilename);
fprintf('Done!\n');

end

