function [containPartVolume, groupedFnames, groupedDatenum, groupedDatasize] = groupPartialVolumeFiles(dataPath, varargin)
% check and group filenames if they are parts of the same volume. The
% filenames are ordered as .tif (with no parts), _part0001.tif, _part0002.tif ...
% 
% Either dataPath or a list of file full path (in same folder) should be provided.
% 
% Author: Xiongtao Ruan (10/06/2020)
% xruan (10/15/2020): also add support for only processing first timepoint
% xruan (02/12/2022): add channel patterns to exclude unnecessary files.
% xruan (12/20/2022): use actual data size for tiff files instead of compressed size.


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPath', @(x) isempty(x) || ischar(x));
ip.addOptional('fileFullpathList', {}, @(x) isempty(x) || iscell(x));
ip.addParameter('ext', '.tif', @(x) ischar(x));
ip.addParameter('onlyFirstTP', false, @(x) ischar(x));
ip.addParameter('ChannelPatterns', {}, @(x) iscell(x));

ip.parse(dataPath, varargin{:});

fileFullpathList = ip.Results.fileFullpathList;
ext = ip.Results.ext;
onlyFirstTP = ip.Results.onlyFirstTP;
ChannelPatterns = ip.Results.ChannelPatterns;

fprintf('Check partial volume files... ');

if isempty(dataPath)
    if onlyFirstTP
        fileFullpathList = fileFullpathList(contains(fileFullpathList, 'Iter_0000_'));
    end
    [dataPath, ~, ext] = fileparts(fileFullpathList{1});
    nF = numel(fileFullpathList);    
    fnames = cell(nF, 1);
    datenum = zeros(nF, 1);
    datasize = zeros(nF, 1);
    for i = 1 : nF
        dir_info = dir(fileFullpathList{i});
        fnames{i} = dir_info.name;
        datenum(i) = dir_info.datenum;
        datasize(i) = dir_info.bytes;
    end
else
    if onlyFirstTP
        expression = '*_Iter_0000_*';
    else
        expression = '*';
    end
    dir_info = dir([dataPath, filesep, expression, ext]);
    fnames = {dir_info.name}';
    datenum = [dir_info.datenum];
    datasize = [dir_info.bytes]';
    
    if ~isempty(ChannelPatterns)
        include_flag = false(numel(fnames), 1);
        fns = cellfun(@(x) [dataPath, filesep, x], fnames, 'unif', 0);
        for c = 1 : numel(ChannelPatterns)
            include_flag = include_flag | contains(fns, ChannelPatterns{c}) | contains(fns, regexpPattern(ChannelPatterns{c}));
        end
        fnames = fnames(include_flag);
        datenum = datenum(include_flag);
        datasize = datasize(include_flag);
    end

    switch ext
        case {'.tif', '.tiff'}
            bim = blockedImage([dataPath, filesep, fnames{1}], 'Adapter', MPageTiffAdapter);            
        case '.zarr'
            bim = blockedImage([dataPath, filesep, fnames{1}], 'Adapter', ZarrAdapter);            
    end

    dtype = bim.ClassUnderlying;
    switch dtype
        case 'uint8'
            dbytes = 1;
        case 'uint16'
            dbytes = 2;
        case 'single'
            dbytes = 4;
        otherwise
            dbytes = 8;
    end
    
    % for tiff file, only check the first and last file and use the max one as data size
    switch ext
        case {'.tif', '.tiff'}
            for f = [1, numel(fnames)]
                fn = [dataPath, filesep, fnames{f}];
                datasize(f) = prod(getImageSize(fn)) * dbytes;
            end
            datasize(:) = max(datasize);
        case {'.zarr'}
            for f = 1 : numel(fnames)
                fn = [dataPath, filesep, fnames{f}];
                datasize(f) = prod(getImageSize(fn)) * dbytes;
            end
    end

    fileFullpathList = cellfun(@(x) [dataPath, filesep, x], fnames, 'unif', 0);
    nF = numel(fileFullpathList);
end

containPartVolume = ~cellfun(@isempty, regexpi(fileFullpathList, ['_part\d+', ext, '$']));
% if there is no partial volume file, just return the filename list and
% datanum in original format (different from the ones with partial volumes)
if ~any(containPartVolume)
    fprintf('\nThere is no partial volume files\n');
    groupedFnames = fnames;
    groupedDatenum = datenum;
    groupedDatasize = datasize;
    return;
end
    
% fsnames = cell(nF, 1);
% for i = 1 : nF
%     [~, fsnames{i}] = fileparts(fileFullpathList{i});
% end
[~, fsnames] = fileparts(fileFullpathList);

% find the main file without _part\d+ pattern
mstr_inds = cellfun(@isempty, regexpi(fsnames, '_part\d+$'));
mstr = fsnames(mstr_inds);
mdn = datenum(mstr_inds);
mds = datasize(mstr_inds);

% find the part files and group them in the defined order
groupedFnames = cell(numel(mstr), 1);
groupedDatenum = cell(numel(mstr), 1);
groupedDatasize = cell(numel(mstr), 1);
for i = 1 : numel(mstr)
    mstr_i = mstr{i};
    mdn_i = mdn(i);
    mds_i = mds(i);
    pstr_inds = ~cellfun(@isempty, regexpi(fsnames, [mstr_i, '_part\d+$']));
    
    if any(pstr_inds)
        pdn = datenum(pstr_inds);
        pds = datasize(pstr_inds);
        pstr_i = sort(fsnames(pstr_inds));
        gstr = [mstr_i; pstr_i];
        groupedFnames{i} = cellfun(@(x) [x, ext], gstr, 'unif', 0);
        groupedDatenum{i} = [mdn_i; pdn(:)];
        groupedDatasize{i} = [mds_i; pds(:)];
    else
        groupedFnames{i} = {[mstr_i, ext]};
        groupedDatenum{i} = mdn_i;
        groupedDatasize{i} = mds_i;
    end
end

fprintf('Done!\n');

end
