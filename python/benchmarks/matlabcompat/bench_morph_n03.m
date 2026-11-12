function bench_morph_n03(threads, repeats)
% The N03 benchmark operations, MATLAB side (the Python side is bench_morph_n03.py; the same
% data, made in each process in its own layout):
%   imerode_disk11     imerode of a 512 x 1800 x 500 mask (a band as in skewed data) by
%                      strel('disk', 11)
%   mex_erode_disk11   the same erosion by decon_mask_edge_erosion_mex (the disk's decomposition
%                      with morphop's padding, no border clearing, no z line)
%   bwdist_2048        bwdist of a 2048 x 2048 slice with sparse features
%   medfilt3_psf       medfilt3 of a 128 x 128 x 101 single volume, [3 3 3]
%   mex_medfilt3_psf   the same by omw_backprojector_mex('medfilt3')
%   bwlabeln_vol       bwlabeln of a 512 x 512 x 200 mask, connectivity 26
%   bwconncomp_vol     bwconncomp of the same mask
% Each is timed repeats times after one untimed call, at each thread count (maxNumCompThreads).

x = 1 : 1800;
z = reshape(1 : 500, 1, 1, []);
A = repmat((x > 0.6 * z + 90) & (x < 0.6 * z + 900), 512, 1, 1);
se = strel('disk', 11);
seq = decompose(se);
nh = arrayfun(@getnhood, seq, 'unif', 0);
[pu, pl] = getpadsize(se);
padUL = pl(1 : 2);
padLR = pu(1 : 2);

[yy, xx] = ndgrid(1 : 2048, 1 : 2048);
S = mod(yy * 7 + xx * 13, 997) == 0;

y = (1 : 128)';
x2 = 1 : 128;
z2 = reshape(1 : 101, 1, 1, []);
P = single(mod(y * 7 + x2 * 13 + z2 * 17, 4096));

y3 = (1 : 512)';
x3 = 1 : 512;
z3 = reshape(1 : 200, 1, 1, []);
V = mod(y3 * 7 + x3 * 13 + z3 * 17, 31) < 12;

pk5d_bench_mark('ready');
for t = threads
    maxNumCompThreads(t);
    pk5d_bench_record('imerode_disk11', t, timed(@() imerode(A, se), repeats));
    pk5d_bench_record('mex_erode_disk11', t, ...
        timed(@() decon_mask_edge_erosion_mex(A, nh, padUL, padLR, 0, false, t), repeats));
    pk5d_bench_record('bwdist_2048', t, timed(@() bwdist(S), repeats));
    pk5d_bench_record('medfilt3_psf', t, timed(@() medfilt3(P, [3, 3, 3]), repeats));
    pk5d_bench_record('mex_medfilt3_psf', t, timed(@() omw_backprojector_mex('medfilt3', P, t), repeats));
    pk5d_bench_record('bwlabeln_vol', t, timed(@() bwlabeln(V, 26), repeats));
    pk5d_bench_record('bwconncomp_vol', t, timed(@() bwconncomp(V, 26), repeats));
end

end


function [ts] = timed(f, repeats)
% seconds per call of f after one untimed call
o = f(); %#ok<NASGU>
clear o
ts = zeros(1, repeats);
for i = 1 : repeats
    t0 = tic;
    o = f(); %#ok<NASGU>
    ts(i) = toc(t0);
    clear o
end

end
