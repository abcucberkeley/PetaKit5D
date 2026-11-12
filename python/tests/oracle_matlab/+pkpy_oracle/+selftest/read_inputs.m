function [names, tif, zr, matv] = read_inputs(dataPath, matFn)
% read the inputs the harness copied: the names under dataPath (hidden files
% included), a tiff, a zarr and a MAT-file


d = dir(fullfile(dataPath, '**'));
d = d(~[d.isdir]);
base = strrep(dataPath, '\', '/');
if ~endsWith(base, '/')
    base = [base, '/'];
end
names = sort(strrep(cellfun(@(f, n) [strrep(f, '\', '/'), '/', n], {d.folder}, {d.name}, 'unif', 0), base, ''));
tif = readtiff([dataPath, 'a.tif']);
zr = readzarr([dataPath, 'b.zarr']);
matv = load(matFn);

end
