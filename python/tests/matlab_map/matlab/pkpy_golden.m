function [varargout] = pkpy_golden(what, scenario, fields, work)
% the golden scenarios of the MATLAB suite (tests/golden/GoldenOutputTest.m) for the
% oracle cases of the port (tests/matlab_map/test_golden.py). The outputs are the fields of
% the scenario's struct, in the order of fields.
%
%       what : 'stored'  the reference outputs in tests/golden/data/<scenario>.mat
%              'compute' GoldenOutputTest.computeScenario(scenario, work) of this tree
%              'inputs'  the arrays computeScenario makes its outputs from, generated the
%                        same way (rng(1, 'twister'), then the same calls in the same order)
%   scenario : char. A scenario of GoldenOutputTest.
%     fields : cell array of char. All field names of the scenario's struct (an error otherwise).
%       work : char. Folder for the files of a scenario (with trailing slash).

root = fileparts(which('setup'));
addpath(fullfile(root, 'tests'), fullfile(root, 'tests', 'golden'));
if nargin < 4 || isempty(work)
    work = [strrep(tempname, '\', '/'), '/'];
    mkdir(work);
end

switch what
    case 'stored'
        ref = load(fullfile(root, 'tests', 'golden', 'data', [scenario, '.mat']));
        s = ref.expected;
    case 'compute'
        s = GoldenOutputTest.computeScenario(scenario, work);
    case 'inputs'
        s = scenario_inputs(scenario);
    otherwise
        error('pkpy_golden:what', 'Unknown request %s', what);
end
fns = cellstr(fields);
if ~isequal(sort(fieldnames(s)), sort(fns(:)))
    error('pkpy_golden:fields', 'Fields of %s %s: %s', what, scenario, strjoin(fieldnames(s)', ', '));
end
varargout = cellfun(@(x) s.(x), fns, 'unif', 0);

end


function [s] = scenario_inputs(name)
% as GoldenOutputTest.computeScenario and its private helpers make them

rng(1, 'twister');
switch name
    case {'deskew_rotate_sample_scan', 'deskew_rotate_skewed_interp', 'deskew_rotate_resampled'}
        s.vol = single(skewed_volume());
    case {'decon_rl_simplified', 'decon_omw'}
        [s.im, s.psf] = blurred_beads();
    case 'psf_gen_masked'
        psf = pktest.beads([63, 63, 41], [32, 32, 21], [1.6, 2.2], 2000, 100);
        psf = psf + pktest.beads([63, 63, 41], [32, 32, 21], [12, 1.2], 50, 0);
        s.psf = single(round(psf + 3 * randn(size(psf))));
    case 'stitch_feather_xcorr'
        s.gt = pktest.texture([32, 104, 32], 3, [1000, 5000]);
    case 'point_detection'
        pos = [16.3, 20.7, 18.2; 16, 48, 32; 48, 24, 40; 40, 40, 20; 32, 32, 48];
        s.vol = double(pktest.beads([64, 64, 64], pos, 1.5, 200, 100)) + 5 * randn(64, 64, 64);
    otherwise
        error('pkpy_golden:scenario', 'Unknown golden output scenario %s', name);
end

end


function [vol] = skewed_volume()
[xr, zr] = ndgrid([12, 24, 36], [8, 20, 32]);
raw = [4 + 2.6 * mod((0 : 8)' * 4, 9), xr(:), zr(:)];
vol = pktest.beads([32, 48, 40], raw, [1.5, 2], 1000, 100) + pktest.texture([32, 48, 40], 3, [0, 200]);
vol = uint16(vol);
end


function [im, psf] = blurred_beads()
psf = pktest.beads([21, 21, 21], [11, 11, 11], [1.5, 2.5], 1, 0);
psf = psf / sum(psf(:));
obj = pktest.beads([32, 32, 32], [10, 12, 9; 22, 20, 24; 16, 8, 20], 0.8, 5000, 0) + ...
    pktest.texture([32, 32, 32], 2, [0, 300]);
im = imfilter(obj, double(psf), 'circular', 'conv') + 100;
im = single(round(im + 2 * randn(size(im))));
end
