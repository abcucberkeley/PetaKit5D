function [J, err_mat, overlapOut] = decon_lucy_omw_function(I, PSF_f, PSF_b, NUMIT, varargin)
% RL decon with alternative backprojector OTF masked wiener (OMW) method
%   save16bit: for GPU compute, convert to 16bit before transfer from GPU to CPU 
%   bbox: for GPU compute, crop the data first before transfer from GPU to CPU
%   useMex: run the iterations with half-spectrum FFTs in decon_lucy_rfft_mex
%       (CPU, the same results as the MATLAB code, also the debug outputs) or
%       decon_lucy_rfft_gpu_mex (GPU); empty: as decon_use_rfft_mex decides
%   overlapFcn: function the CPU runs while the GPU computes (after the
%       iterations are queued, before the result is copied back; on the CPU
%       after the iterations); its output is overlapOut
%   useGPU: true runs on a GPU and errors (decon:noUsableGPU) if none is
%       usable, false on the CPU, [] (the default) on a GPU if one is usable,
%       else on the CPU with a warning
%
% Author: Xiongtao Ruan (11/10/2022)
%
% xruan (08/31/2023): add scaleFactor, deconOffset, dampFactor, and
%   EdgeErosion, and use input parser for parameters.


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('I', @isnumeric);
ip.addRequired('PSF_f', @isnumeric);
ip.addRequired('PSF_b', @isnumeric);
ip.addRequired('NUMIT', @isnumeric);
ip.addParameter('Background', [], @isnumeric);
ip.addParameter('useGPU', [], @(x) isempty(x) || islogical(x)); % GPU processing: [] if one is usable, true requires one
ip.addParameter('save16bit', true, @islogical);
ip.addParameter('dampFactor', 1, @isnumeric); % damp factor for decon result
ip.addParameter('scaleFactor', 1, @isnumeric); % scale factor for result
ip.addParameter('deconOffset', 0, @isnumeric); % offset for decon result
ip.addParameter('EdgeErosion', 0, @isnumeric); % edge erosion for decon result
ip.addParameter('bbox', [], @isnumeric); % bounding box to crop data after decon
ip.addParameter('debug', false, @islogical);
ip.addParameter('debug_folder', './debug/', @ischar);
ip.addParameter('saveStep', 1, @isnumeric); % save intermediate results every given iterations
ip.addParameter('useMex', [], @(x) isempty(x) || islogical(x)); % half-spectrum FFT mex for the iterations
ip.addParameter('overlapFcn', [], @(x) isempty(x) || isa(x, 'function_handle')); % CPU work during the GPU iterations

ip.parse(I, PSF_f, PSF_b, NUMIT, varargin{:});

pr = ip.Results;
Background = pr.Background;
useGPU = pr.useGPU;
% [] (auto): on a GPU if one is usable; true: a usable GPU is required
requireGPU = isequal(useGPU, true);
useGPU = isempty(useGPU) || useGPU;
save16bit = pr.save16bit;
dampFactor = pr.dampFactor;
scaleFactor = pr.scaleFactor;
deconOffset = pr.deconOffset;
EdgeErosion = pr.EdgeErosion;
bbox = pr.bbox;
debug = pr.debug;
debug_folder = pr.debug_folder;
saveStep = pr.saveStep;

if isempty(Background)
    Background = 0;
end

PSF_f = single(PSF_f);
PSF_b = single(PSF_b);
sizeI = size(I);

% 1. Prepare PSFs: the cached OTFs are used for the same PSFs (bit for bit),
% image size, GPU and mex use, and number of threads; cached CPU OTFs of a
% call without a GPU are computed again when a GPU is required
persistent OTF_f OTF_b psz last_PSF_f last_PSF_b last_nThreads gpuNum last_useGPU last_useMex
nThreads = maxNumCompThreads;
newOTF = isempty(psz) || isempty(OTF_f) || ~(numel(psz) == numel(sizeI) && all(psz == sizeI) && last_useGPU == useGPU ...
        && last_nThreads == nThreads && same_psf(PSF_f, last_PSF_f) && same_psf(PSF_b, last_PSF_b)) ...
        || (useGPU && ~isempty(OTF_f) && gpuNum > 0 && ~existsOnGPU(OTF_f)) || (requireGPU && isequal(gpuNum, 0));
if newOTF
    % (no cached OTF if the check errors)
    gpuNum = 0;
    psz = [];
    if useGPU
        % (GPUs newer than this MATLAB get forward compatibility; errors if
        % a GPU is required and none is usable, else warns with MATLAB's reason)
        gpuNum = decon_gpu_num(requireGPU);
    end
end
% the mex uses half OTFs, the MATLAB code full OTFs
useMex = decon_use_rfft_mex(pr.useMex, 'omw', I, useGPU && gpuNum > 0, debug, Background);
if newOTF || ~isequal(last_useMex, useMex)
    psz = sizeI;
    last_useGPU = useGPU;
    last_useMex = useMex;
    last_nThreads = nThreads;
    last_PSF_f = gather(PSF_f);
    last_PSF_b = gather(PSF_b);

    if useGPU && gpuNum > 0
        PSF_f = gpuArray(PSF_f);
        PSF_b = gpuArray(PSF_b);
        PSF_f = PSF_f ./ sum(PSF_f(:), 'double');
        PSF_b = PSF_b ./ sum(PSF_b(:), 'double');
    else
        % summed in an order that doesn't depend on the number of threads
        PSF_f = PSF_f ./ decon_sum_double(PSF_f);
        PSF_b = PSF_b ./ decon_sum_double(PSF_b);
    end
    if useMex
        OTF_f = [];
        OTF_b = [];
        OTF_f = decon_lucy_rfft(useGPU && gpuNum > 0, 'otf', PSF_f, double(sizeI));
        OTF_b = decon_lucy_rfft(useGPU && gpuNum > 0, 'otf', PSF_b, double(sizeI));
    else
        OTF_f = decon_psf2otf(PSF_f, double(sizeI));
        OTF_b = decon_psf2otf(PSF_b, double(sizeI));
    end
end
clear PSF_f PSF_b;

if gpuNum == 0
    useGPU = false;
end

if useMex
    % iterations in the mex, which also subtracts the background
    if useGPU
        I = gpuArray(I);
    elseif ~any(strcmp(class(I), {'uint16', 'single'}))
        I = single(I);
    end
    if deconOffset ~= 0 || EdgeErosion > 0
        I_mask = I ~= 0;
    end
    if debug
        % (CPU only) the iterations in pieces up to the debug outputs
        [J, err_mat] = omw_mex_debug(I, OTF_f, OTF_b, NUMIT, Background, saveStep, debug_folder);
    else
        J = decon_lucy_rfft(useGPU, 'omw', I, OTF_f, OTF_b, NUMIT, Background);
        err_mat = zeros(NUMIT, 4);
    end
    % with 0 iterations, the background-subtracted input
    if dampFactor > 1
        I = decon_lucy_rfft(useGPU, 'omw', I, OTF_f, OTF_b, 0, Background);
    end
else
    if useGPU
        I = gpuArray(I);
        if deconOffset ~= 0 || EdgeErosion > 0
            I_mask = I ~= 0;
        end
        if strcmp(underlyingType(I), 'uint16')
            % keep uint16 input as uint16 on the GPU (2 instead of 4 bytes per
            % voxel); it is cast to the class of J where it is used (exact)
            I = I - uint16(Background);            
        else
            I = max(I - Background, 0);    
        end
    else
        if ~isa(I, 'double')
            I = single(I);
        end
        if deconOffset ~= 0 || EdgeErosion > 0
            I_mask = I ~= 0;
        end    
        I = max(I - Background, 0);
    end
    J = I;
    if strcmp(underlyingType(J), 'uint16')
        J = single(J);
    end

    % 2. L_R Iterations
    estep = saveStep;
    err_mat = zeros(NUMIT, 4);

    for k = 1 : NUMIT
        % valina RL update step
        CX = real(ifftn(fftn(J) .* OTF_f));
        % CX = CX + (CX < eps) * eps;
        CX = max(CX, eps);
        % one FFT step per statement, and CX freed at the end, so that at most two
        % complex volumes are on the GPU at a time (the same operations)
        CX = fftn(cast(I, 'like', J) ./ CX);
        CX = CX .* OTF_b;
        CX = ifftn(CX);
        J = real(CX) .* J;
        clear CX;
        % J = J .* (J > 0);
        J = max(J, 0);
            
        if debug
            err_mat(k, 1:2) = [k, mean((double(J) - double(I)) .^ 2, 'all')];
            err_mat(k, 3) = err_mat(k, 2) ./ err_mat(1, 2); 
            if rem(k, estep) == 0 && k > estep * 2
                istp = k/estep;
                err_mat(k, 4) = min(abs(err_mat(istp, 3) - err_mat(istp-1, 3)) / 10, abs(err_mat(istp, 3) + err_mat(istp-2, 3) - 2 * err_mat(istp-1, 3)) / 2);
            else
                if k < 3 * estep
                    err_mat(k, 4) = 1;
                else
                    err_mat(k, 4) = err_mat(istp * estep, 4);
                end
            end
        
            if k == 1
                mkdir(debug_folder);
            end
            if rem(k, estep) == 0
                if useGPU 
                    writetiff(uint16(gather(J)), sprintf('%s/Iter_%04d.tif', debug_folder, k));
                else
                    writetiff(uint16(J), sprintf('%s/Iter_%04d.tif', debug_folder, k));
                end
            end
        end
    end

    if debug
        err_mat = err_mat(1 : floor(NUMIT / estep), :);
    end
end

% the caller's CPU work while the GPU runs the queued iterations
overlapOut = [];
if ~isempty(pr.overlapFcn)
    overlapOut = pr.overlapFcn();
end

% post-processing of deconvolved result
if dampFactor > 1
    I = cast(I, 'like', J);
    J = J .* (J <= dampFactor * I) + min(J, dampFactor * I .* (J >= dampFactor * I));
end
clear I CX;

if EdgeErosion > 0
    I_mask = decon_mask_edge_erosion(I_mask, EdgeErosion);
end

if scaleFactor ~= 1 || deconOffset ~= 0 || EdgeErosion > 0
    if deconOffset ~= 0 || EdgeErosion > 0
        J = (J * scaleFactor + deconOffset) .* I_mask;
    else
        J = J * scaleFactor;
    end
end

if save16bit
    J = uint16(J);
end

if useGPU
    if ~isempty(bbox)
        J = J(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6));
    end    
    J = gather(J);
else
    if ~isempty(bbox)
        J = crop3d(J, bbox);
    end
end

end


function [J, err_mat] = omw_mex_debug(I, OTF_f, OTF_b, NUMIT, Background, estep, debug_folder)
% the debug outputs of the MATLAB loop (err_mat, the result every estep
% iterations) with the iterations in the CPU mex, run in pieces from the last
% result to the iterations the outputs need: the error metric of the rows of
% err_mat that are kept (the first floor(NUMIT / estep) iterations) and the
% saved results. The rest of err_mat is computed for every iteration as in the
% loop (it is cheap and errors in the same cases).

% the background-subtracted input of the MATLAB loop
Is = max(single(I) - Background, 0);
nErr = floor(NUMIT / estep);
kLast = floor(NUMIT);
err_mat = zeros(NUMIT, 4);
J = [];
kd = 0;
for k = 1 : NUMIT
    if k <= nErr || rem(k, estep) == 0 || k == kLast
        if isempty(J)
            J = decon_lucy_rfft(false, 'omw', I, OTF_f, OTF_b, k, Background);
        else
            J = decon_lucy_rfft(false, 'omw', I, OTF_f, OTF_b, k - kd, Background, J);
        end
        kd = k;
    end

    err_mat(k, 1) = k;
    if k <= nErr
        err_mat(k, 2) = mean((double(J) - double(Is)) .^ 2, 'all');
    end
    err_mat(k, 3) = err_mat(k, 2) ./ err_mat(1, 2); 
    if rem(k, estep) == 0 && k > estep * 2
        istp = k/estep;
        err_mat(k, 4) = min(abs(err_mat(istp, 3) - err_mat(istp-1, 3)) / 10, abs(err_mat(istp, 3) + err_mat(istp-2, 3) - 2 * err_mat(istp-1, 3)) / 2);
    else
        if k < 3 * estep
            err_mat(k, 4) = 1;
        else
            err_mat(k, 4) = err_mat(istp * estep, 4);
        end
    end

    if k == 1
        mkdir(debug_folder);
    end
    if rem(k, estep) == 0
        writetiff(uint16(J), sprintf('%s/Iter_%04d.tif', debug_folder, k));
    end
end
if isempty(J)
    J = decon_lucy_rfft(false, 'omw', I, OTF_f, OTF_b, 0, Background);
end
err_mat = err_mat(1 : floor(NUMIT / estep), :);

end


function [tf] = same_psf(a, b)
% whether single arrays a and b have the same size and bits (unlike isequal,
% -0 differs from 0 and a NaN equals itself); b on the CPU

tf = isequal(size(a), size(b)) && isequal(typecast(gather(a(:)), 'uint32'), typecast(b(:), 'uint32'));

end
