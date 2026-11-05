function [tasks, scenarios] = isolation_prepare(rootDir, varargin)
% build the synthetic datasets of the isolation harness in rootDir and capture
% the tasks that the data wrappers hand to the job framework (see
% isolation_capture): each scenario folder is left as it is just before its
% tasks run, so running the tasks (in any order, in any process) writes their
% outputs there. The tasks keep the absolute paths of rootDir, so they must
% run with the scenario folders at the same place (restore a copy of rootDir
% there before each run).
%
% Parameters (as 'specifier'-value pairs):
%    scenarios : cell array (default: {}, all). Names of the scenarios to build (see scenario_table below).
%          gpu : true|false|[] (default: []). Include the GPU scenarios ([]: if a GPU is available).
%        heavy : true|false (default: false). Include the scenarios that need much memory (PSF analysis, ~22 GB).
%     demoData : char (default: ''). Folder that holds PetaKit5D_demo_cell_image_dataset: include the
%                scenarios on two of its tiles (copied to rootDir).
%      explore : true|false (default: false). Run all framework calls of each scenario and print them (no tasks).
%
% Output:
%   tasks : struct array (also saved in rootDir/tasks.mat) with fields
%              id : char. <scenario>_<number>.
%        scenario : char. Scenario name (its folder in rootDir).
%            func : char. Task function.
%         funcStr : char. Task function string, as the job framework runs it.
%           input : char or cell array. Input the framework checks for the task.
%          output : char. Output the framework checks for the task ('' if none).
%             gpu : true|false. Needs a GPU.
%         threads : [] or a number. maxNumCompThreads for the task ([]: as the process has).
%   scenarios : struct array of the scenarios built (name, gpu, heavy, knownIssue, description).


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('rootDir', @ischar);
ip.addParameter('scenarios', {}, @iscell);
ip.addParameter('gpu', [], @(x) isempty(x) || islogical(x));
ip.addParameter('heavy', false, @islogical);
ip.addParameter('demoData', '', @ischar);
ip.addParameter('explore', false, @islogical);
ip.parse(rootDir, varargin{:});

pr = ip.Results;
rootDir = strrep(rootDir, '\', '/');
if ~endsWith(rootDir, '/')
    rootDir = [rootDir, '/'];
end
useGPU = pr.gpu;
if isempty(useGPU)
    useGPU = false;
    try
        useGPU = usable_gpu_num() > 0;
    catch
    end
end

demoDir = '';
if ~isempty(pr.demoData)
    demoDir = [strrep(pr.demoData, '\', '/'), '/PetaKit5D_demo_cell_image_dataset/'];
end
scTable = scenario_table(demoDir);
sel = true(numel(scTable), 1);
if ~isempty(pr.scenarios)
    sel = ismember({scTable.name}', pr.scenarios);
end
sel = sel & (~[scTable.gpu]' | useGPU) & (~[scTable.heavy]' | pr.heavy) & (~[scTable.demo]' | ~isempty(demoDir));
scTable = scTable(sel);

if ~exist(rootDir, 'dir')
    mkdir(rootDir);
end

tasks = struct('id', {}, 'scenario', {}, 'func', {}, 'funcStr', {}, 'input', {}, 'output', {}, 'gpu', {}, 'threads', {});
for s = 1 : numel(scTable)
    sc = scTable(s);
    dataDir = [rootDir, sc.name, '/'];
    if exist(dataDir, 'dir')
        rmdir(dataDir, 's');
    end
    mkdir(dataDir);
    stopAtCall = sc.stopAtCall;
    if pr.explore
        stopAtCall = inf;
    end
    fprintf('Scenario %s: %s\n', sc.name, sc.description);
    rs = RandStream('twister', 'Seed', s);
    if pr.explore
        try
            calls = isolation_capture(@() sc.build(dataDir, rs), 'stopAtCall', stopAtCall, 'singleJob', sc.singleJob);
        catch ME
            fprintf('  failed: %s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
            continue;
        end
        for c = 1 : numel(calls)
            funcs = unique(cellfun(@(x) regexp(x, '^\w+', 'match', 'once'), calls{c}.funcStrs, 'unif', 0));
            fprintf('  call %d (%s): %d tasks of %s\n', c, calls{c}.type, numel(calls{c}.funcStrs), strjoin(funcs, ', '));
        end
        continue;
    end
    calls = isolation_capture(@() sc.build(dataDir, rs), 'stopAtCall', stopAtCall, 'singleJob', sc.singleJob);

    % the tasks of the stopped call, or all single job submissions
    if sc.singleJob
        calls = calls(cellfun(@(x) strcmp(x.type, 'singleJob'), calls));
        funcStrs = unique(cellfun(@(x) x.funcStrs{1}, calls, 'unif', 0), 'stable');
        [outputs, inputs] = cellfun(@decon_output, funcStrs, 'unif', 0);
    else
        if isempty(calls) || numel(calls) < sc.stopAtCall
            error('Scenario %s made %d framework calls, fewer than %d.', sc.name, numel(calls), sc.stopAtCall);
        end
        funcStrs = calls{sc.stopAtCall}.funcStrs(:);
        outputs = calls{sc.stopAtCall}.outputFullpaths(:);
        inputs = calls{sc.stopAtCall}.inputFullpaths(:);
    end
    for i = 1 : numel(funcStrs)
        t.id = sprintf('%s_%02d', sc.name, i);
        t.scenario = sc.name;
        t.func = regexp(funcStrs{i}, '^\w+', 'match', 'once');
        t.funcStr = funcStrs{i};
        t.input = inputs{i};
        t.output = outputs{i};
        t.gpu = sc.gpu;
        t.threads = [];
        tasks(end + 1) = t;
    end
    fprintf('  %d tasks of %s\n', numel(funcStrs), strjoin(unique({tasks(strcmp({tasks.scenario}, sc.name)).func}), ', '));
end

scenarios = rmfield(scTable, {'build', 'stopAtCall', 'singleJob'});
if ~pr.explore
    save([rootDir, 'tasks.mat'], 'tasks', 'scenarios');
end

end


function [scTable] = scenario_table(demoDir)
% name, builder, framework call whose tasks are captured, single job
% submissions instead (decon), GPU, heavy, known issue, description; the
% scenarios on the demo dataset (demoDir) come last

c = {
    'dsr', @build_dsr, 1, false, false, false, false, 'deskew/rotate frames, tiff and MIPs'
    'dsr_zarr_ff', @build_dsr_zarr_ff, 1, false, false, false, false, 'deskew/rotate frames with flat field, zarr'
    'dsr_large', @build_dsr_large, 1, false, false, false, false, 'large-file deskew/rotate (XR_deskewRotateZarr tasks)'
    'dsr_blocks', @build_dsr_blocks, 1, false, false, false, false, 'large-file deskew/rotate blocks'
    'decon_omw', @(d, rs) build_decon(d, rs, 'omw', false), 1, true, false, false, false, 'OMW decon frames, CPU'
    'decon_rl', @(d, rs) build_decon(d, rs, 'simplified', false), 1, true, false, false, false, 'RL decon frames, CPU'
    'decon_omw_gpu', @(d, rs) build_decon(d, rs, 'omw', true), 1, true, true, false, false, 'OMW decon frames, GPU (one multi-file task)'
    'decon_rl_gpu', @(d, rs) build_decon(d, rs, 'simplified', true), 1, true, true, false, false, 'RL decon frames, GPU (one multi-file task)'
    'decon_blocks', @build_decon_blocks, 1, false, false, false, false, 'large-file OMW decon blocks'
    'decon_psf_a', @(d, rs) build_decon_psf_pair(d, rs, false), 1, true, false, false, true, 'RL decon with PSF A'
    'decon_psf_b', @(d, rs) build_decon_psf_pair(d, rs, true), 1, true, false, false, true, 'RL decon with PSF B (A mirrored in z: same size and energy)'
    'tiff2zarr', @build_tiff2zarr, 1, false, false, false, false, 'tiff to zarr'
    'zarr2tiff', @build_zarr2tiff, 1, false, false, false, false, 'zarr to tiff'
    'stitch_convert', @(d, rs) build_stitch(d, rs, false), 1, false, false, false, false, 'stitching: tile conversion with flat field'
    'stitch_xcorr', @(d, rs) build_stitch(d, rs, false), 2, false, false, false, false, 'stitching: tile registration'
    'stitch_blockinfo', @(d, rs) build_stitch(d, rs, false), 3, false, false, false, false, 'stitching: block info'
    'stitch_bwdist', @(d, rs) build_stitch(d, rs, false), 4, false, false, false, false, 'stitching: distance maps'
    'stitch_blocks', @(d, rs) build_stitch(d, rs, false), 5, false, false, false, false, 'stitching: feather blending blocks'
    'stitch_large_mips', @(d, rs) build_stitch(d, rs, true), 2, false, false, false, false, 'large-file stitching: MIP slabs of the tiles'
    'stitch_large_bwdist', @(d, rs) build_stitch(d, rs, true), 14, false, false, false, false, 'large-file stitching: distance maps from MIP slabs'
    'stitch_large_blocks', @(d, rs) build_stitch(d, rs, true), 15, false, false, false, false, 'large-file stitching: blending blocks'
    'resample', @build_resample, 1, false, false, false, false, 'resample frames'
    'resample_blocks', @build_resample_blocks, 1, false, false, false, false, 'resample zarr blocks'
    'crop', @build_crop, 1, false, false, false, false, 'crop frames'
    'crop_large', @build_crop_blocks, 1, false, false, false, false, 'large-file crop (XR_crop_frame of a zarr)'
    'crop_blocks', @build_crop_blocks, 2, false, false, false, false, 'crop zarr blocks'
    'mip', @build_mip, 1, false, false, false, false, 'MIPs of tiffs'
    'mip_blocks', @build_mip_blocks, 1, false, false, false, false, 'MIP blocks of a zarr'
    'unmix', @build_unmix, 1, false, false, false, false, 'channel unmixing frames'
    'unmix_blocks', @build_unmix_blocks, 2, false, false, false, false, 'channel unmixing zarr blocks'
    'fsc', @build_fsc, 1, false, false, false, false, 'FSC analysis'
    'fft', @build_fft, 1, false, false, false, false, 'FFT spectrum'
    'resave_blocks', @build_resave_blocks, 1, false, false, false, false, 'zarr resave blocks'
    'chromatic', @build_chromatic, 1, false, false, false, false, 'chromatic shift correction'
    'demo_fw', @build_demo_fw, 1, false, false, false, false, 'the framework demo''s crop and deskew/rotate function'
    'psf_analysis', @build_psf_analysis, 1, false, false, true, false, 'PSF analysis plots'
    };
c(:, end + 1) = {false};
d = {
    'demo_dsr', @(d, rs) build_demo_dsr(d, demoDir), 1, false, false, false, false, 'demo tiles: deskew/rotate', true
    'demo_decon', @(d, rs) build_demo_decon(d, demoDir, false), 1, true, false, false, false, 'demo tiles: OMW decon, CPU', true
    'demo_decon_gpu', @(d, rs) build_demo_decon(d, demoDir, true), 1, true, true, false, false, 'demo tiles: OMW decon, GPU', true
    };
scTable = cell2struct([c; d], {'name', 'build', 'stopAtCall', 'singleJob', 'gpu', 'heavy', 'knownIssue', 'description', 'demo'}, 2);

end


function [output, input] = decon_output(funcStr)
% the result of the last file of an XR_RLdeconFrame3D task string (its
% output for the job framework), and its first file

output = '';
input = '';
tok = regexp(funcStr, '^XR_RLdeconFrame3D\((\{[^}]*\}|''[^'']*''),[^,]+,[^,]+,''([^'']*)''', 'tokens', 'once');
if isempty(tok)
    return;
end
frames = regexp(tok{1}, '''([^'']*)''', 'tokens');
input = frames{1}{1};
[~, fsn] = fileparts(frames{end}{1});
ext = '.tif';
if contains(funcStr, '''saveZarr'',true')
    ext = '.zarr';
end
deconPath = tok{2};
if ~endsWith(deconPath, '/')
    deconPath = [deconPath, '/'];
end
output = [deconPath, fsn, ext];

end


%% data

function [fsn] = llsm_name(cam, t)
% LLSM-style file name of camera cam (e.g. 'A') and time point t

fsn = sprintf('Scan_Iter_%04d_0000_Cam%s_ch0_CAM1_stack0000_488nm_%07dmsec_0106060251msecAbs_000x_000y_000z_%04dt', ...
    t, cam, t * 1000, t);

end


function [im] = bead_frame(rs, sz, nBeads, sigma, background)
% uint16 frame of random blurred beads on a background

pos = (0.1 + 0.8 * rand(rs, nBeads, 3)) .* sz;
amp = 500 + 1500 * rand(rs, nBeads, 1);
im = uint16(pktest.beads(sz, pos, sigma, amp, background) + 5 * rand(rs, sz));

end


function [fns] = write_frames(dataDir, rs, sz, names, ext)
% bead frames written as tiff or zarr

fns = cell(numel(names), 1);
for i = 1 : numel(names)
    im = bead_frame(rs, sz, 30, [1.5, 2], 100);
    fns{i} = [dataDir, names{i}, ext];
    if strcmp(ext, '.zarr')
        writezarr(im, fns{i}, 'blockSize', [32, 32, 32]);
    else
        writetiff(im, fns{i});
    end
end

end


function [ffFns, bgFns] = write_flat_fields(dataDir, sz)
% flat field and background images for cameras A and B

mkdir([dataDir, 'ff']);
bg = zeros(sz(1 : 2), 'single') + 100;
ffA = bg + single(linspace(500, 1000, sz(2)));
ffB = bg + single(linspace(1000, 500, sz(1)))';
ffFns = {[dataDir, 'ff/ffA.tif'], [dataDir, 'ff/ffB.tif']};
bgFns = {[dataDir, 'ff/bg.tif'], [dataDir, 'ff/bg.tif']};
writetiff(ffA, ffFns{1});
writetiff(ffB, ffFns{2});
writetiff(bg, bgFns{1});

end


%% deskew/rotate

function [] = build_dsr(dataDir, rs)

write_frames(dataDir, rs, [128, 256, 101], {llsm_name('A', 0), llsm_name('B', 0), llsm_name('A', 1)}, '.tif');
XR_deskew_rotate_data_wrapper({dataDir}, 'deskew', true, 'rotate', true, 'DSRCombined', true, ...
    'xyPixelSize', 0.108, 'dz', 0.2, 'skewAngle', 32.45, 'reverse', true, 'channelPatterns', {'CamA', 'CamB'}, ...
    'save16bit', true, 'saveMIP', true, 'parseCluster', false, 'masterCompute', true);

end


function [] = build_dsr_zarr_ff(dataDir, rs)

sz = [128, 256, 101];
write_frames(dataDir, rs, sz, {llsm_name('A', 0), llsm_name('B', 0)}, '.tif');
[ffFns, bgFns] = write_flat_fields(dataDir, sz);
XR_deskew_rotate_data_wrapper({dataDir}, 'deskew', true, 'rotate', true, 'DSRCombined', true, ...
    'xyPixelSize', 0.108, 'dz', 0.3, 'skewAngle', 32.45, 'reverse', true, 'channelPatterns', {'CamA', 'CamB'}, ...
    'FFCorrection', true, 'lowerLimit', 0.4, 'FFImagePaths', ffFns, 'backgroundPaths', bgFns, ...
    'save16bit', true, 'saveZarr', true, 'blockSize', [32, 64, 32], 'parseCluster', false, 'masterCompute', true);

end


function [] = build_dsr_large(dataDir, rs)

write_frames(dataDir, rs, [96, 192, 81], {llsm_name('A', 0), llsm_name('A', 1)}, '.zarr');
XR_deskew_rotate_data_wrapper({dataDir}, 'deskew', true, 'rotate', true, 'DSRCombined', true, ...
    'xyPixelSize', 0.108, 'dz', 0.2, 'skewAngle', 32.45, 'reverse', true, 'channelPatterns', {'CamA'}, ...
    'zarrFile', true, 'largeFile', true, 'batchSize', [32, 1024, 1024], 'blockSize', [32, 64, 64], ...
    'save16bit', true, 'parseCluster', false, 'masterCompute', true);

end


function [] = build_dsr_blocks(dataDir, rs)

fns = write_frames(dataDir, rs, [96, 192, 81], {llsm_name('A', 0)}, '.zarr');
XR_deskewRotateZarr(fns{1}, 0.108, 0.2, 'skewAngle', 32.45, 'reverse', true, 'batchSize', [24, 1024, 1024], ...
    'blockSize', [24, 64, 64], 'saveMIP', true, 'parseCluster', false);

end


%% deconvolution

function [psf] = decon_psf(sz, sigma, mirrored)
% Gaussian PSF with a weaker second lobe 3 slices below its peak (above it
% when mirrored), zero in its first and last 5 slices

c = floor(sz / 2) + 1;
psf = pktest.beads(sz, [c; c + [0, 0, 3]], sigma, [10000; 2500], 0);
psf(:, :, [1 : 5, end - 4 : end]) = 0;
if mirrored
    psf = flip(psf, 3);
end
psf = uint16(psf);

end


function [] = write_decon_data(dataDir, rs, nT, mirrored)

mkdir([dataDir, 'data']);
mkdir([dataDir, 'psf']);
writetiff(decon_psf([31, 31, 31], [1.5, 2], mirrored), [dataDir, 'psf/PSF_CamA.tif']);
for t = 1 : nT
    writetiff(bead_frame(rs, [96, 96, 64], 20, sqrt(0.8 ^ 2 + [1.5, 2] .^ 2), 100), ...
        [dataDir, 'data/', llsm_name('A', t - 1), '.tif']);
end

end


function [] = build_decon(dataDir, rs, method, gpu)
% the jobs the decon wrapper submits on a cluster (the OMW back projector is
% made by the wrapper then); with edge erosion of each file's own mask (the
% mask of the first time point would come from its job)

write_decon_data(dataDir, rs, 3, false);
iter = 10;
if strcmp(method, 'omw')
    iter = 2;
end
% GPU jobs of two and one files
taskBatchNum = [];
if gpu
    taskBatchNum = 2;
end
XR_decon_data_wrapper({[dataDir, 'data/']}, 'channelPatterns', {'CamA_ch0'}, ...
    'psfFullpaths', {[dataDir, 'psf/PSF_CamA.tif']}, 'xyPixelSize', 0.108, 'dz', 0.2, 'dzPSF', 0.2, ...
    'background', 100, 'RLMethod', method, 'deconIter', iter, 'wienerAlpha', 0.005, 'OTFCumThresh', 0.88, ...
    'skewed', true, 'edgeErosion', 4, 'erodeByFTP', false, 'save16bit', true, 'parseCluster', true, ...
    'GPUJob', gpu, 'taskBatchNum', taskBatchNum, 'jobLogDir', [dataDir, 'logs'], 'maxTrialNum', 1);

end


function [] = build_decon_psf_pair(dataDir, rs, mirrored)
% the same frames deconvolved with PSF A or B (A mirrored in z); in one
% process, the second one reused the first one's OTF while the cache was keyed
% on the size and sum of squares of the PSF (the same for both;
% DeconFunctionTest/changedPSFWithSameEnergyIsNotIgnored)

% the same frames for A and B
rs = RandStream('twister', 'Seed', 100);
write_decon_data(dataDir, rs, 2, mirrored);
% integer values that sum to 2^16 and whose squares sum to less than 2^24,
% so that the normalized PSF and its sum of squares are exact in single
% precision in any summation order (the old keys of A and B are equal)
sz = [31, 31, 31];
c = floor(sz / 2) + 1;
psf = round(pktest.beads(sz, [c; c + [0, 0, 3]], [2.2, 2.87], [240; 60], 0));
psf(:, :, [1 : 5, end - 4 : end]) = 0;
psf(c(1), c(2), c(3)) = psf(c(1), c(2), c(3)) + 2 ^ 16 - sum(psf(:));
if sum(psf(:)) ~= 2 ^ 16 || min(psf(:)) < 0 || sum(psf(:) .^ 2) >= 2 ^ 24
    error('The PSF for the same-energy pair is not exact (sum %g, min %g, sum of squares %g).', ...
        sum(psf(:)), min(psf(:)), sum(psf(:) .^ 2));
end
if mirrored
    psf = flip(psf, 3);
end
writetiff(uint16(psf), [dataDir, 'psf/PSF_CamA.tif']);
XR_decon_data_wrapper({[dataDir, 'data/']}, 'channelPatterns', {'CamA_ch0'}, ...
    'psfFullpaths', {[dataDir, 'psf/PSF_CamA.tif']}, 'xyPixelSize', 0.108, 'dz', 0.2, 'dzPSF', 0.2, ...
    'background', 100, 'RLMethod', 'simplified', 'deconIter', 5, 'psfGen', false, 'edgeErosion', 0, ...
    'save16bit', true, 'parseCluster', true, 'GPUJob', false, 'jobLogDir', [dataDir, 'logs'], 'maxTrialNum', 1);

end


function [] = build_decon_blocks(dataDir, rs)

write_decon_data(dataDir, rs, 1, false);
fn = [dataDir, 'data/', llsm_name('A', 0)];
writezarr(readtiff([fn, '.tif']), [fn, '.zarr'], 'blockSize', [16, 16, 16]);
delete([fn, '.tif']);
XR_decon_data_wrapper({[dataDir, 'data/']}, 'channelPatterns', {'CamA_ch0'}, ...
    'psfFullpaths', {[dataDir, 'psf/PSF_CamA.tif']}, 'xyPixelSize', 0.108, 'dz', 0.2, 'dzPSF', 0.2, ...
    'background', 100, 'RLMethod', 'omw', 'deconIter', 2, 'wienerAlpha', 0.005, 'OTFCumThresh', 0.88, ...
    'skewed', true, 'edgeErosion', 0, 'save16bit', true, 'zarrFile', true, 'largeFile', true, ...
    'largeMethod', 'inplace', 'batchSize', [48, 48, 32], 'blockSize', [16, 16, 16], ...
    'parseCluster', false, 'GPUJob', false);

end


%% conversion

function [] = build_tiff2zarr(dataDir, rs)

write_frames(dataDir, rs, [128, 200, 60], {llsm_name('A', 0), llsm_name('B', 0), llsm_name('A', 1)}, '.tif');
XR_tiffToZarr_wrapper(dataDir, 'channelPatterns', {'CamA', 'CamB'}, 'blockSize', [32, 32, 16], ...
    'parseCluster', false, 'masterCompute', true);

end


function [] = build_zarr2tiff(dataDir, rs)

write_frames(dataDir, rs, [128, 200, 60], {llsm_name('A', 0), llsm_name('B', 0)}, '.zarr');
XR_zarrToTiff_wrapper(dataDir, 'channelPatterns', {'CamA', 'CamB'}, 'parseCluster', false, 'masterCompute', true);

end


%% stitching

function [csvFn] = write_stitch_tiles(dataDir, rs, gt, starts, tileSize, stageXYZ)
% tiles cut from gt at 0-based offsets starts ([y, x, z]); stage positions
% (um) from the offsets unless given

px = 0.108;
nT = size(starts, 1);
if nargin < 6
    stageXYZ = starts(:, [2, 1, 3]) * px;
end
xyzStarts = starts(:, [2, 1, 3]);
tileIdx = zeros(nT, 3);
for i = 1 : 3
    [~, ~, tileIdx(:, i)] = unique(xyzStarts(:, i));
end
tiles = struct('im', cell(nT, 1), 'idx', [], 'xyz', []);
for i = 1 : nT
    s = starts(i, :);
    im = gt(s(1) + 1 : s(1) + tileSize(1), s(2) + 1 : s(2) + tileSize(2), s(3) + 1 : s(3) + tileSize(3));
    % some noise per tile, so that overlaps differ a little
    tiles(i).im = im + uint16(10 * rand(rs, size(im)));
    tiles(i).idx = tileIdx(i, :) - 1;
    tiles(i).xyz = stageXYZ(i, :);
end
[~, csvFn] = pktest.writeStitchTiles(dataDir, tiles);

end


function [] = build_stitch(dataDir, rs, largeFile)
% 2 x 2 tiles, the last one listed 2 voxels off in x, registered and
% feather blended; without largeFile, the tiles are flat field corrected
% when they are converted (a process function, as in
% demo_phase_and_2photon_stitching.m)

px = 0.108;
tileSize = [80, 84, 40];
gt = uint16(1000 + 4000 * imgaussfilt3(rand(rs, 120, 140, 40), 1.5));
starts = [0, 0, 0; 0, 56, 0; 40, 0, 0; 40, 56, 0];
stageXYZ = starts(:, [2, 1, 3]) * px;
stageXYZ(4, 1) = stageXYZ(4, 1) + 2 * px;
csvFn = write_stitch_tiles(dataDir, rs, gt, starts, tileSize, stageXYZ);
args = {'DSR', true, 'resampleType', 'isotropic', 'xyPixelSize', px, 'axisOrder', 'x,y,z', ...
    'channelPatterns', {'CamA_ch0'}, 'primaryCh', 'CamA_ch0', 'xcorrShift', true, 'xcorrDownsample', [1, 1, 1], ...
    'xyMaxOffset', 20, 'zMaxOffset', 10, 'blendMethod', 'feather', 'parseCluster', false};
if largeFile
    args = [args, {'largeFile', true, 'poolSize', [2, 2, 2], 'batchSize', [32, 32, 16], 'blockSize', [16, 16, 16]}];
else
    mkdir([dataDir, 'flatfield']);
    ffFn = [dataDir, 'flatfield/flatfield.tif'];
    bgFn = [dataDir, 'flatfield/background.tif'];
    writetiff(uint16(900 + 200 * rand(rs, tileSize(1 : 2))), ffFn);
    writetiff(100 * ones(tileSize(1 : 2), 'uint16'), bgFn);
    usrFun = sprintf("@(x)erodeVolumeBy2DProjection(processFFCorrectionFrame(x,'%s','%s',%d,%d,%s),%0.10d)", ...
        ffFn, bgFn, 1, 0.4, string(false), 0);
    funFn = [dataDir, 'processFunction.mat'];
    save('-v7.3', funFn, 'usrFun');
    args = [args, {'batchSize', [64, 64, 32], 'blockSize', [32, 32, 16], 'processFunPath', {funFn}}];
end
XR_matlab_stitching_wrapper(dataDir, csvFn, args{:});

end


%% tools

function [] = build_resample(dataDir, rs)

write_frames(dataDir, rs, [100, 120, 60], {llsm_name('A', 0), llsm_name('B', 0)}, '.tif');
XR_resample_dataset(dataDir, [2, 2, 1.5], 'channelPatterns', {'CamA', 'CamB'}, 'parseCluster', false, 'masterCompute', true);

end


function [] = build_resample_blocks(dataDir, rs)

fns = write_frames(dataDir, rs, [100, 120, 60], {llsm_name('A', 0)}, '.zarr');
XR_resampleSingleZarr(fns{1}, [dataDir, 'resampled/', llsm_name('A', 0), '.zarr'], [2, 2, 2], ...
    'batchSize', [32, 32, 16], 'blockSize', [16, 16, 8], 'borderSize', [4, 4, 4], 'parseCluster', false);

end


function [] = build_crop(dataDir, rs)

write_frames(dataDir, rs, [100, 120, 60], {llsm_name('A', 0), llsm_name('B', 0)}, '.tif');
XR_crop_dataset(dataDir, [5, 7, 3, 90, 101, 50], 'channelPatterns', {'CamA', 'CamB'}, 'parseCluster', false, 'masterCompute', true);

end


function [] = build_crop_blocks(dataDir, rs)

write_frames(dataDir, rs, [100, 120, 60], {llsm_name('A', 0)}, '.zarr');
XR_crop_dataset(dataDir, [5, 7, 3, 90, 101, 50], 'channelPatterns', {'CamA'}, 'zarrFile', true, 'largeFile', true, ...
    'batchSize', [32, 32, 16], 'blockSize', [16, 16, 8], 'parseCluster', false, 'masterCompute', true);

end


function [] = build_mip(dataDir, rs)

write_frames(dataDir, rs, [100, 120, 60], {llsm_name('A', 0), llsm_name('B', 0)}, '.tif');
XR_MIP_wrapper(dataDir, 'axis', [1, 1, 1], 'channelPatterns', {'CamA', 'CamB'}, 'parseCluster', false, 'masterCompute', true);

end


function [] = build_mip_blocks(dataDir, rs)

fns = write_frames(dataDir, rs, [100, 120, 60], {llsm_name('A', 0)}, '.zarr');
XR_MIP_zarr(fns{1}, 'batchSize', [32, 32, 16], 'parseCluster', false);

end


function [] = build_unmix(dataDir, rs)

write_frames(dataDir, rs, [100, 120, 40], {llsm_name('A', 0), llsm_name('B', 0), llsm_name('A', 1), llsm_name('B', 1)}, '.tif');
XR_unmix_channels_data_wrapper(dataDir, 'unmixFactors', [1, -0.1], 'mode', 'linear', 'channelPatterns', {'CamA', 'CamB'}, ...
    'channelInd', 1, 'constBackground', [100, 100], 'constOffset', 1, 'save16bit', true, 'parseCluster', false, ...
    'masterCompute', true);

end


function [] = build_unmix_blocks(dataDir, rs)

write_frames(dataDir, rs, [100, 120, 40], {llsm_name('A', 0), llsm_name('B', 0)}, '.zarr');
XR_unmix_channels_data_wrapper(dataDir, 'unmixFactors', [1, -0.1], 'mode', 'linear', 'channelPatterns', {'CamA', 'CamB'}, ...
    'channelInd', 1, 'constBackground', [100, 100], 'zarrFile', true, 'largeFile', true, 'saveZarr', true, ...
    'batchSize', [64, 64, 32], 'blockSize', [32, 32, 16], 'parseCluster', false, 'masterCompute', true);

end


function [] = build_fsc(dataDir, rs)

b = imgaussfilt3(randn(rs, 64, 64, 64, 'single'), 2);
signal = 1000 + 300 * b / std(b, 0, 'all');
for i = 1 : 3
    writetiff(signal + 200 * i * randn(rs, 64, 64, 64, 'single'), sprintf('%sIter_%04d.tif', dataDir, 2 * i));
end
XR_FSC_analysis_wrapper({dataDir}, 'resultDirName', 'FSCs', 'xyPixelSize', 0.1, 'dz', 0.1, 'dr', 2, ...
    'dtheta', pi / 6, 'resThreshMethod', 'fixed', 'resThresh', 0.2, 'halfSize', [24, 24, 24], ...
    'channelPatterns', {'Iter'}, 'channels', 488, 'iterInterval', 2, 'parseCluster', false, 'masterCompute', true);

end


function [] = build_fft(dataDir, rs)

write_frames(dataDir, rs, [62, 60, 24], {llsm_name('A', 0), llsm_name('A', 1)}, '.tif');
XR_fftSpectrumComputingWrapper({dataDir}, 'resultDirName', 'FFT', 'xyPixelSize', 0.108, 'dz', 0.3, ...
    'outSize', [81, 81, 81], 'outPixelSize', 0.108, 'save3DStack', true, 'parseCluster', false, 'masterCompute', true);

end


function [] = build_resave_blocks(dataDir, rs)

fns = write_frames(dataDir, rs, [100, 120, 60], {llsm_name('A', 0)}, '.zarr');
XR_resaveSingleZarr(fns{1}, [dataDir, 'resaved/', llsm_name('A', 0), '.zarr'], 'blockSize', [16, 16, 8], ...
    'batchSize', [64, 64, 32], 'parseCluster', false);

end


function [] = build_chromatic(dataDir, rs)

write_frames(dataDir, rs, [80, 90, 30], {llsm_name('A', 0), llsm_name('B', 0)}, '.tif');
XR_chromatic_shift_correction_data_wrapper({dataDir}, 'mode', 'same', 'channelPatterns', {'CamA', 'CamB'}, ...
    'chromaticOffset', [0, 0, 0; 1, -2, 1], 'parseCluster', false, 'masterCompute', true);

end


function [] = build_demo_fw(dataDir, rs)
% task strings built as in demos/demo_generic_computing_framework.m

outPath = [dataDir, 'Cropped_DSR/'];
mkdir(outPath);
names = {llsm_name('A', 0), llsm_name('B', 0)};
inputFullnames = strcat(dataDir, names, '.tif');
outputFullnames = strcat(outPath, names, '.tif');
funcStrs = cell(2, 1);
for f = 1 : 2
    writetiff(bead_frame(rs, [128, 256, 60], 30, [1.5, 2], 100), inputFullnames{f});
    funcStrs{f} = sprintf(['crop_deskew_rotate_demo_function(''%s'',''%s'',%d,', ...
        '%d,%d,%s)'], inputFullnames{f}, outputFullnames{f}, 32.45, 0.108, 0.3, strrep(mat2str([1, 1, 1, 128, 200, 60]), ' ', ','));
end
generic_computing_frameworks_wrapper(inputFullnames, outputFullnames, funcStrs, parseCluster=false, ...
    masterCompute=true, memAllocate=20, configFile='', mccMode=false);

end


function [] = copy_demo_tiles(dataDir, demoDir)
% the CamA and CamB tiles of the first position (100 slices each) and the
% two PSFs of the demo dataset

pos = '000x_000y_000z_0000t';
for cam = {'CamA', 'CamB'}
    d = dir(sprintf('%sScan_Iter_0000_0000_%s_ch0_*%s.tif', demoDir, cam{1}, pos));
    copyfile([d(1).folder, '/', d(1).name], dataDir);
end
mkdir([dataDir, 'PSF']);
copyfile([demoDir, 'PSF/488_2_c.tif'], [dataDir, 'PSF/']);
copyfile([demoDir, 'PSF/560_c.tif'], [dataDir, 'PSF/']);

end


function [] = build_demo_dsr(dataDir, demoDir)
% the demos' deskew/rotate settings

copy_demo_tiles(dataDir, demoDir);
XR_deskew_rotate_data_wrapper({dataDir}, 'deskew', true, 'rotate', true, 'DSRCombined', true, ...
    'xyPixelSize', 0.108, 'dz', 0.3, 'skewAngle', 32.45, 'reverse', true, 'channelPatterns', {'CamA_ch0', 'CamB_ch0'}, ...
    'save16bit', true, 'saveMIP', true, 'parseCluster', false, 'masterCompute', true);

end


function [] = build_demo_decon(dataDir, demoDir, gpu)
% demo_RL_deconvolution.m's OMW settings, as jobs on a cluster

copy_demo_tiles(dataDir, demoDir);
XR_decon_data_wrapper({dataDir}, 'channelPatterns', {'CamA_ch0', 'CamB_ch0'}, ...
    'psfFullpaths', {[dataDir, 'PSF/488_2_c.tif'], [dataDir, 'PSF/560_c.tif']}, 'xyPixelSize', 0.108, ...
    'dz', 0.3, 'dzPSF', 0.3, 'reverse', true, 'background', 100, 'RLMethod', 'omw', 'deconIter', 2, ...
    'wienerAlpha', 0.005, 'OTFCumThresh', 0.87, 'skewed', true, 'edgeErosion', 0, 'save16bit', true, ...
    'parseCluster', true, 'GPUJob', gpu, 'jobLogDir', [dataDir, 'logs'], 'maxTrialNum', 1);

end


function [] = build_psf_analysis(dataDir, rs)

psfPath = [dataDir, 'PSF/'];
mkdir([psfPath, 'RW_PSFs/']);
writetiff(uint16(pktest.beads([64, 64, 101], [32, 32, 51], [1.6, 3], 2000, 100) + 5 * rand(rs, [64, 64, 101])), ...
    [psfPath, '488_CamA_ch0.tif']);
rwFn = [psfPath, 'RW_PSFs/rw.tif'];
writetiff(single(pktest.beads([128, 128, 101], [64, 64, 51], [2, 6], 1, 0)), rwFn);
XR_psf_analysis_wrapper({psfPath}, 'xyPixelSize', 0.108, 'dz', 0.3, 'skewAngle', 32.45, ...
    'channelPatterns', {'CamA_ch0'}, 'channels', 488, 'deskew', true, 'RWFn', {rwFn}, ...
    'sourceStr', 'test', 'visible', false, 'parseCluster', false, 'masterCompute', true);

end
