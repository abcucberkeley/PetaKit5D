function bench_resize_n01(threads, repeats)
% The N01 benchmark operations, MATLAB side (the Python side is bench_resize_n01.py; the same
% data, made in each process in its own layout): a uint16 512 x 1800 x 500 volume
%   rs3_half_cubic, rs3_half_linear   imresize3 to round(size .* [0.5 0.5 0.5])
%   rs3_z2_cubic, rs3_z2_linear       imresize3 to round(size .* [1 1 2])
%   rs3_zhalf_linear                  imresize3 to round(size ./ [1 1 2]) (XR_resampleFrame's
%                                     resample factor [1 1 2])
%   fast_half_cubic, fast_z2_cubic    imresize3_fast to the same sizes
% Each is timed repeats times after one untimed call, at each thread count (maxNumCompThreads).

y = single(1 : 512)';
x = single(1 : 1800);
z = single(reshape(1 : 500, 1, 1, []));
V = uint16(mod(y * 7 + x * 13 + z * 17, single(4000)));
sz = size(V);
half = round(sz .* [0.5, 0.5, 0.5]);
z2 = round(sz .* [1, 1, 2]);
zh = round(sz ./ [1, 1, 2]);

pk5d_bench_mark('ready');
for t = threads
    maxNumCompThreads(t);
    pk5d_bench_record('rs3_half_cubic', t, timed(@() imresize3(V, half), repeats));
    pk5d_bench_record('rs3_half_linear', t, timed(@() imresize3(V, half, 'linear'), repeats));
    pk5d_bench_record('rs3_z2_cubic', t, timed(@() imresize3(V, z2), repeats));
    pk5d_bench_record('rs3_z2_linear', t, timed(@() imresize3(V, z2, 'linear'), repeats));
    pk5d_bench_record('rs3_zhalf_linear', t, timed(@() imresize3(V, zh, 'linear'), repeats));
    pk5d_bench_record('fast_half_cubic', t, timed(@() imresize3_fast(V, half), repeats));
    pk5d_bench_record('fast_z2_cubic', t, timed(@() imresize3_fast(V, z2), repeats));
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
