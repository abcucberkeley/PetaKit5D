function [] = XR_bioformats_to_tiff_or_zarr_wrapper(dataPaths, varargin)
% convert various microscopy data format (nd2 and czi by default) to separated 3d tiff or zarr files (xyz)
%
% Required inputs:
%           dataPaths : char or cell array. Directory paths for the datasets. Either a string for a single dataset or a cell array of paths for several datasets with same experimental settings.
%
% Parameters (as 'specifier'-value pairs):
%       resultDirName : char (default: 'tiffs'). Result directory under data paths for saving results.
%     channelPatterns : cell array (default: {'.nd2', '.czi'}). File name patterns or extensions for identifying the image channels to be included.
%           nChannels : numeric scalar (default: 2). Number of channels in each dataset.
%        dataFormats  : cell array (default: {'.nd2', '.czi'}). File formats to be processed. Other formats can be added as long as bioFormats can read it.
%           saveZarr  : true|false (default: false). Save output data in Zarr format.
%           blockSize : 1x3 numeric vector (default: [256, 256, 256]). Block/chunk size for Zarr output.
%          overWrite  : true|false (default: false). Overwrite existing output if it exists.
%               uuid  : char (default: ''). UUID string used as part of the temporary or result paths.
% 
% Author: Xiongtao Ruan (08/04/2025)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addParameter('resultDirName', 'tiffs', @ischar);
ip.addParameter('channelPatterns', {'.nd2', '.czi'}, @iscell);
ip.addParameter('nChannels', 2, @isnumeric);
ip.addParameter('dataFormats', {'.nd2', '.czi'}, @iscell);
ip.addParameter('saveZarr', false, @islogical);
ip.addParameter('blockSize', [256, 256, 256], @isnumeric);
ip.addParameter('overWrite', false, @islogical);
ip.addParameter('uuid', '', @ischar);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
resultDirName = pr.resultDirName;
channelPatterns = pr.channelPatterns;
nChannels = pr.nChannels;
dataFormats = pr.dataFormats;
saveZarr = pr.saveZarr;
blockSize = pr.blockSize;
overWrite = pr.overWrite;
uuid = pr.uuid;

if ischar(dataPaths)
    dataPaths = {dataPaths};
end

nd = numel(dataPaths);
resultPaths = cell(nd, 1);
for d = 1 : nd
    dataPath = dataPaths{d};
    if ~strcmp(dataPath(end), filesep)
        dataPaths{d} = [dataPath, filesep];
    end

    dataPaths{d} = [simplifyPath(dataPaths{d}), '/'];
    resultPath = [dataPaths{d}, '/' resultDirName, '/'];
    resultPaths{d} = resultPath;
    mkdir(resultPath);
    fileattrib(resultPath, '+w', 'g');
end

if isempty(uuid)
    uuid = get_uuid();
end

% parse image filenames
zarrFile = false;
[fnames, fsns, fd_inds, filepaths] = parseImageFilenames(dataPaths, zarrFile, channelPatterns, dataFormats=dataFormats);
nF = numel(fnames);

ext = '.tif';
if saveZarr
    ext = '.zarr';
end

for f = 1 : nF
    tic
    fname = fnames{f};
    fsn = fsns{f};
    fd_ind = fd_inds(f);
    resultPath = resultPaths{fd_ind};
    fn = sprintf('%s/%s', dataPaths{fd_ind}, fname);
    v = bfopen(fn);
    v = v(:, 1);
    
    M = size(v, 1);
    for m = 1 : M
        vol_m = cat(3, v{m}{:, 1});
        vol_m = reshape(vol_m, size(vol_m, 1), size(vol_m, 2), nChannels, []);
        
        for c = 1 : nChannels
            vc = squeeze(vol_m(:, :, c, :));
            
            if M == 1
                fnout = sprintf('%s/%s_ch%d%s', resultPath, fsn, c - 1, ext);
                tmpout = sprintf('%s/%s_ch%d_%s%s', resultPath, fsn, c - 1, uuid, ext);
            else
                fnout = sprintf('%s/%s_fov%04d_ch%d%s', resultPath, fsn, m, c - 1, ext);
                tmpout = sprintf('%s/%s_fov%04d_ch%d_%s%s', resultPath, fsn, m, c - 1, uuid, ext);
            end
            if ~overWrite && (exist(fnout, 'file') || exist(fnout, 'dir'))
                continue;
            end
            if saveZarr
                blockSize = min(size(vc), blockSize);
                dimSeparator = '.';
                if prod(size(vc) ./ blockSize) > 10000
                    dimSeparator = '/';
                end
                writezarr(vc, tmpout, blockSize=blockSize, dimSeparator=dimSeparator);                
            else
                writetiff(vc, tmpout);
            end
            movefile(tmpout, fnout);
        end
    end
    toc
end

end
