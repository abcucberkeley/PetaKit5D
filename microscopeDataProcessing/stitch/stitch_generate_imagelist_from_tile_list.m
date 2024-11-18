function [] = stitch_generate_imagelist_from_tile_list(dataPaths, tileFilenames, tileIndices, tileInterval, varargin)
% generate image list file from a given list of tile files.
% 
%
% Author: Xiongtao Ruan (11/18/2024)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addRequired('tileFilenames', @(x) ischar(x) || iscell(x));
ip.addRequired('tileIndices', @(x) isnumeric(x) && size(x, 2) == 5); % 5 x nF: tcxyz
ip.addRequired('tileInterval', @(x) isnumeric(x) && size(x, 2) == 3); % in um
ip.addParameter('uuid', '', @ischar);

ip.parse(dataPaths, tileFilenames, tileIndices, tileInterval, varargin{:});

pr = ip.Results;
uuid = pr.uuid;

if isempty(uuid)
    uuid = get_uuid();
end

if ischar(dataPaths)
    dataPaths = {dataPaths};
end

fnames = tileFilenames;
nF = numel(fnames);

filepaths = cellfun(@(x) [dataPaths{1}, x], fnames, 'UniformOutput', false);

tab = cell2table(cat(2, filepaths(:), fnames(:)), 'VariableNames', ["Filepath", "Filename"]);

for f = 1 : nF
    tab.prefix{f} = 'Scan';
    tab.Iter(f) = tileIndices(f, 1);
    tab.subIter{f} = '';
    tab.fullIter{f} = sprintf('%04d', tileIndices(f, 1));
    tab.camera(f) = 'A';

    tab.ch(f) = tileIndices(f, 2);
    tab.stack(f) = 0;
    tab.x(f) = tileIndices(f, 3);
    tab.y(f) = tileIndices(f, 4);
    tab.z(f) = tileIndices(f, 5);
    tab.t(f) = tileIndices(f, 1);
    tab.mappedFilename(f) = {sprintf(['Scan_Iter_%s_Cam%s_ch%d_CAM1_stack%04d_', ...
        '488nm_000000msec_000000msecAbs_%03dx_%03dy_%03dz_%04dt.tif'], ...
        tab.fullIter{f}, tab.camera(f), tab.ch(f), tab.stack(f), tab.x(f), ...
        tab.y(f), tab.z(f), tab.t(f))};
end


% get image sizes for all tiles
nF = size(tab, 1);

% generate coordiantes 
tileIndices_xyz = [tab.x, tab.y, tab.z];
ts = min(tileIndices_xyz);
tt = max(tileIndices_xyz);

% calculate physical size in xyz order
fullIter = unique(tab.fullIter);
Ch = unique(tab.ch);
Cam = unique(tab.camera);
stackn = unique(tab.stack);

xyz = zeros(nF, 3);

for n = 1:numel(fullIter)
    for ncam = 1:numel(Cam)
        for s = 1:numel(stackn)
            for c = 1:numel(Ch)
                inds = tab.ch == Ch(c) & tab.camera == Cam(ncam) & strcmp(tab.fullIter, fullIter{n}) & tab.stack == stackn(s);
                cur_tab = tab(inds, :);                            

                tileIndices_xyz = [cur_tab.x, cur_tab.y, cur_tab.z];

                xyz(inds, :) = (tileIndices_xyz - ts) .* tileInterval;
            end
        end
    end
end

% visualize tile positions
if false
    figure, scatter3(xyz(:, 1), xyz(:, 2), xyz(:, 3), '*')
end

% save image list csv file
counter = 0;
sz = [nF, 9];
varTypes = ["string", "string", "string", "double", "double", "double", "double", "double", "double"];
varNames = ["Filepath","Filename","mappedFilename","StageX_um_","StageY_um_","StageZ_um_","ObjectiveX_um_","ObjectiveY_um_","ObjectiveZ_um_"];

t = table('Size', sz, 'VariableTypes', varTypes, 'VariableNames', varNames);

for f = 1 : nF
    t(counter + 1, :) = {tab.Filepath{f}, tab.Filename{f}, tab.mappedFilename{f}, xyz(f, 1), xyz(f, 2), xyz(f, 3), xyz(f, 1), xyz(f, 2), xyz(f, 3)};
    counter = counter + 1;
end
t(counter + 1 : end, :) = [];

tmpout = sprintf('%s/ImageList_from_tile_list_%s.csv', dataPaths{1}, uuid);
fnout = sprintf('%s/ImageList_from_tile_list.csv', dataPaths{1});
writetable(t, tmpout, 'Delimiter', ',');
fileattrib(tmpout, '+w', 'g');
movefile(tmpout, fnout);

fprintf('Image list file %s is generated!\n', fnout);

end

