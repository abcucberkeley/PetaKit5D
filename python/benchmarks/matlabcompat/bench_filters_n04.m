function bench_filters_n04(threads, repeats)
% The N04 benchmark operations, MATLAB side (the Python side is bench_filters_n04.py; the same
% data, made in each process in its own layout):
%   ig3_unmix_single     imgaussfilt3 of a single 512 x 512 x 400 block, sigma [1 1 1] (the
%                        Gaussian unmixing)
%   ig3_puncta_uint16    imgaussfilt3 of a uint16 512 x 512 x 400 block, sigma 2.5 (puncta removal)
%   ig2_phaseff_2048     imgaussfilt of a single 2048 x 2048 image, sigma 100 (the phase flat
%                        field; the frequency domain)
%   ig2_phaseff_1024     the same on a 1024 x 1024 image (the spatial domain)
% Each is timed repeats times after one untimed call, at each thread count (maxNumCompThreads).

y = single(1 : 512)';
x = single(1 : 512);
z = single(reshape(1 : 400, 1, 1, []));
V = mod(y * 7 + x * 13 + z * 17, single(1000)) / single(1000);
U = uint16(mod(y * 7 + x * 13 + z * 17, single(4000)));
y2 = single(1 : 2048)';
x2 = single(1 : 2048);
P2 = mod(y2 * 7 + x2 * 13, single(1000)) / single(1000);
P1 = P2(1 : 1024, 1 : 1024);

pk5d_bench_mark('ready');
for t = threads
    maxNumCompThreads(t);
    pk5d_bench_record('ig3_unmix_single', t, timed(@() imgaussfilt3(V, [1, 1, 1]), repeats));
    pk5d_bench_record('ig3_puncta_uint16', t, timed(@() imgaussfilt3(U, 2.5), repeats));
    pk5d_bench_record('ig2_phaseff_2048', t, timed(@() imgaussfilt(P2, 100), repeats));
    pk5d_bench_record('ig2_phaseff_1024', t, timed(@() imgaussfilt(P1, 100), repeats));
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
