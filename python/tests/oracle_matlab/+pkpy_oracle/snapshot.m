function [s] = snapshot(root)
% relative path, folder flag, size and modification time of everything
% under root (hidden files included)


d0 = dir(root);
base = strrep(d0(strcmp({d0.name}, '.')).folder, '\', '/');
d = dir(fullfile(root, '**'));
d = d(~ismember({d.name}, {'.', '..'}));
rel = strrep(cellfun(@(f, n) [f, '/', n], {d.folder}, {d.name}, 'unif', 0), '\', '/');
rel = cellfun(@(x) x(numel(base) + 2 : end), rel, 'unif', 0);
s = struct('rel', {reshape(rel, [], 1)}, 'isdir', reshape([d.isdir], [], 1), ...
    'bytes', reshape([d.bytes], [], 1), 'datenum', reshape([d.datenum], [], 1));

end
