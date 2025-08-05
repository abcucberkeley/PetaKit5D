function [] = XR_bioformats_to_tiff_data_wrapper(dataPaths, varargin)
% convert various microscopy data format (nd2 and czi by default) to
% separated 3d tiff files (xyz)
%
% Author: Xiongtao Ruan (08/04/2025)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addParameter('resultDirName', 'tiffs', @ischar);
ip.addParameter('channelPatterns', {'.nd2', '.czi'}, @iscell);
ip.addParameter('nChannels', 2, @isnumeric);
ip.addParameter('dataFormats', {'.nd2', '.czi'}, @iscell);
ip.addParameter('overWrite', false, @islogical);
ip.addParameter('uuid', '', @ischar);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
resultDirName = pr.resultDirName;
channelPatterns = pr.channelPatterns;
nChannels = pr.nChannels;
dataFormats = pr.dataFormats;
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
                fnout = sprintf('%s/%s_ch%d.tif', resultPath, fsn, c - 1);
                tmpout = sprintf('%s/%s_ch%d_%s.tif', resultPath, fsn, c - 1, uuid);
            else
                fnout = sprintf('%s/%s_fov%04d_ch%d.tif', resultPath, fsn, m, c - 1);
                tmpout = sprintf('%s/%s_fov%04d_ch%d_%s.tif', resultPath, fsn, m, c - 1, uuid);
            end
            if ~overWrite && exist(fnout, 'file')
                continue;
            end
            writetiff(vc, tmpout);
            movefile(tmpout, fnout);
        end
    end
    toc
end

end
