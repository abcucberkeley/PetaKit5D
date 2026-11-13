function bench_io_mex(inputFile, workdir, threads, repeats, outFile, label)
% MATLAB side of bench_io.py: readtiff, writetiff, readzarr and writezarr (the
% cpp-tiff and cpp-zarr mex PetaKit5D uses) with maxNumCompThreads set to each
% thread count; seconds per call, after one warm-up call (second: the first timed
% call, the second in the process, as on the Python side), and the minor page
% faults of each call (/proc/self/stat, on Linux). A read's output is cleared
% after toc. Run after setup.

tStart = posixtime(datetime('now', 'TimeZone', 'local'));
if ~exist(workdir, 'dir')
    mkdir(workdir);
end
tifOut = fullfile(workdir, 'mex_out.tif');
zarrOut = fullfile(workdir, 'mex_out.zarr');
zarrIn = fullfile(workdir, 'mex_in.zarr');
nt0 = maxNumCompThreads;
im = readtiff(inputFile);
if exist(zarrIn, 'dir')
    rmdir(zarrIn, 's');
end
writezarr(im, zarrIn, 'blockSize', [256, 256, 256]);
ops = {'tiff_read', 'tiff_write', 'zarr_read', 'zarr_write'};
rows = {};
for t = threads(:)'
    maxNumCompThreads(t);
    for k = 1 : numel(ops)
        op = ops{k};
        ts = zeros(1, repeats);
        flt = zeros(1, repeats + 1);
        for r = 0 : repeats
            f0 = minflt();
            switch op
                case 'tiff_read'
                    t0 = tic; out = readtiff(inputFile); dt = toc(t0); %#ok<NASGU>
                case 'tiff_write'
                    t0 = tic; writetiff(im, tifOut); dt = toc(t0);
                case 'zarr_read'
                    t0 = tic; out = readzarr(zarrIn); dt = toc(t0); %#ok<NASGU>
                case 'zarr_write'
                    if exist(zarrOut, 'dir')
                        rmdir(zarrOut, 's');
                    end
                    t0 = tic; writezarr(im, zarrOut, 'blockSize', [256, 256, 256]); dt = toc(t0);
            end
            flt(r + 1) = minflt() - f0;
            clear out
            if r > 0
                ts(r) = dt;
            end
        end
        rows{end + 1} = struct('impl', 'matlab', 'op', op, 'threads', t, 'times', ts, ...
            'median', median(ts), 'second', ts(1), 'minflt', flt); %#ok<AGROW>
        fprintf('matlab %-18s t=%2d median %.3f s  second %.3f  %s  minflt %s\n', op, t, ...
            median(ts), ts(1), mat2str(ts, 4), mat2str(flt));
    end
end
maxNumCompThreads(nt0);
if exist(zarrOut, 'dir')
    rmdir(zarrOut, 's');
end
rmdir(zarrIn, 's');
if exist(tifOut, 'file')
    delete(tifOut);
end
[~, host] = system('hostname');
res = struct('label', label, 'impl', 'matlab', 'host', strtrim(host), 'version', version, ...
    'shape', size(im), 'dtype', class(im), 't_start', tStart, 't_end', posixtime(datetime('now', 'TimeZone', 'local')), 'rows', {rows});
fid = fopen(outFile, 'w');
fwrite(fid, jsonencode(res));
fclose(fid);
end


function n = minflt()
% this process's minor page faults so far (NaN where /proc is missing)
n = NaN;
try
    s = fileread('/proc/self/stat');
    f = strsplit(strtrim(s(find(s == ')', 1, 'last') + 2 : end)), ' ');
    n = str2double(f{8});
catch
end
end
