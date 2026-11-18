function out = pkpy_n5_summary(n5)
% An N5 dataset as the oracle compares it: {attributes.json text, the chunk files (paths
% relative to the dataset, '/'-separated, sorted), their first 16 bytes (the header of a 3D
% chunk), the data read back by N5Adapter}.

txt = fileread([n5, '/attributes.json']);
d = dir([n5, '/**/*']);
d = d(~[d.isdir]);
n = numel(n5) + 2;
rel = cell(1, numel(d));
for i = 1 : numel(d)
    f = [d(i).folder, '/', d(i).name];
    rel{i} = strrep(f(n : end), '\', '/');
end
rel = sort(rel(~strcmp(rel, 'attributes.json')));
heads = zeros(numel(rel), 16, 'uint8');
for i = 1 : numel(rel)
    fid = fopen([n5, '/', rel{i}], 'r');
    h = fread(fid, 16, '*uint8');
    fclose(fid);
    heads(i, 1 : numel(h)) = h';
end
bim = blockedImage(n5, 'Adapter', N5Adapter);
data = bim.Adapter.getIORegion(ones(1, numel(bim.Size)), bim.Size);
out = {txt, rel, heads, data};

end
