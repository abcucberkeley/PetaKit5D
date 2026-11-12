function [names] = write_files(outPath, im)
% write one file of each kind the harness compares; return their names


mkdir([outPath, 'sub']);
writetiff(im, [outPath, 'sub/im.tif']);
writezarr(im, [outPath, 'im.zarr']);
s = struct('name', 'demo', 'sz', size(im), 'vals', {{1, 'a'}}, 't', 1.5);
save([outPath, 'info.mat'], '-v7.3', '-struct', 's');
save([outPath, 'info_v6.mat'], '-v6', '-struct', 's');
fid = fopen([outPath, 'params.json'], 'w');
fprintf(fid, '%s', jsonencode(struct('a', 1, 'b', 'x', 'c', [1, 2])));
fclose(fid);
fid = fopen([outPath, 'log.txt'], 'w');
fprintf(fid, 'line 1\nwork %s\n', outPath);
fclose(fid);
mkdir([outPath, 'empty']);
flag = ['flag_', char(java.util.UUID.randomUUID()), '.txt'];
fid = fopen([outPath, flag], 'w');
fclose(fid);
names = {'sub/im.tif', 'im.zarr', 'info.mat', 'info_v6.mat', 'params.json', 'log.txt', flag};

end
