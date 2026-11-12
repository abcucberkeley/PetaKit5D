function [J_2, err_mat, k, overlapOut] = decon_lucy_function(I, PSF, NUMIT, varargin)
% adapted from matlab deconvlucy.m
% 
% xruan (05/18/2021): add support for early stop with stop criteria
% xruan (05/19/2021): add support for gpu computing
% xruan (06/02/2021): add support for fix iteration decon
% xruan (06/03/2021): add support for debug (evaluate error and save
% intermediate result).
% xruan (07/09/2021): add output for the number of iterations run; also
% output err mat for fixed iterations.
% xruan (11/11/2021): when using GPU, use single precision unless the data
% is in double
% xruan (11/12/2021): refactor code to save memory (J_4), and disable
% weighting or subsampling
% xruan (11/17/2022): add save16bit, bbox, background option for faster processing 
% on GPU, also remove err_thrsh and fixIter parameters
% xruan (08/31/2023): add scaleFactor, deconOffset, damper, and
%   EdgeErosion, and use input parser for parameters.
%
% useMex: run the FFT steps of the iterations in decon_lucy_rfft_mex (CPU,
%   the same results as the MATLAB code, also the debug outputs) or
%   decon_lucy_rfft_gpu_mex (GPU); empty: as decon_use_rfft_mex decides
% overlapFcn: function the CPU runs while the GPU computes (after the
%   iterations are queued, before the result is copied back; on the CPU after
%   the iterations); its output is overlapOut
% useGPU: true runs on a GPU and errors (decon:noUsableGPU) if none is
%   usable, false on the CPU, [] (the default) on a GPU if one is usable, else
%   on the CPU with a warning
%
% On the CPU, the PSF normalization and the acceleration factor are summed in
% an order that doesn't depend on the number of threads (decon_sum_double,
% decon_lucy_accel), and the acceleration step runs in decon_lucy_accel_mex
% for arrays of 89000 or more elements, with the results of the MATLAB code.


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('I', @isnumeric);
ip.addRequired('PSF', @isnumeric);
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
ip.addParameter('useMex', [], @(x) isempty(x) || islogical(x)); % FFT steps in a mex
ip.addParameter('overlapFcn', [], @(x) isempty(x) || isa(x, 'function_handle')); % CPU work during the GPU iterations

ip.parse(I, PSF, NUMIT, varargin{:});

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

PSF = single(PSF);
sizeI = size(I);

% 1. Prepare PSFs: the cached OTF is used for the same PSF (bit for bit),
% image size, GPU and mex use, and number of threads; a cached CPU OTF of a
% call without a GPU is computed again when a GPU is required
persistent H psz last_PSF last_nThreads gpuNum last_useGPU last_useMex
nThreads = maxNumCompThreads;
newOTF = isempty(psz) || ~(numel(psz) == numel(sizeI) && all(psz == sizeI) && last_useGPU == useGPU ...
        && last_nThreads == nThreads && same_psf(PSF, last_PSF)) || (useGPU && ~isempty(H) && gpuNum > 0 && ~existsOnGPU(H)) ...
        || (requireGPU && isequal(gpuNum, 0));
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
% the mex uses a half OTF, the MATLAB code the full OTF
useMex = decon_use_rfft_mex(pr.useMex, 'simplified', I, useGPU && gpuNum > 0, debug, Background);
if newOTF || ~isequal(last_useMex, useMex)
    psz = sizeI;
    last_useGPU = useGPU;
    last_useMex = useMex;
    last_nThreads = nThreads;
    last_PSF = gather(PSF);

    if useGPU && gpuNum > 0
        PSF = gpuArray(PSF);
        PSF = PSF ./ sum(PSF(:), 'double');
    else
        % summed in an order that doesn't depend on the number of threads
        PSF = PSF ./ decon_sum_double(PSF);
    end
    if useMex
        H = [];
        H = decon_lucy_rfft(useGPU && gpuNum > 0, 'otf', PSF, double(sizeI));
    else
        H = decon_psf2otf(PSF, double(sizeI));
    end
end
clear PSF;

if gpuNum == 0
    useGPU = false;
end

% 2. Prepare parameters for iterations
if useGPU
    I = gpuArray(I);
    if deconOffset ~= 0 || EdgeErosion > 0
        I_mask = I ~= 0;
    end    
    if strcmp(underlyingType(I), 'uint16')
        % keep uint16 input as uint16 on the GPU (2 instead of 4 bytes per
        % voxel); it is cast to the class of J_2 where it is used (exact)
        I = I - uint16(Background);            
        classI = 'single';
        J_2 = single(I);
    else
        I = max(I - Background, 0);    
        classI = underlyingType(I);
        J_2 = I;
    end

    J_3 = zeros(sizeI, classI, 'gpuArray');
    J_4 = zeros(prod(sizeI), 1, classI, 'gpuArray');
else
    if ~isa(I, 'double')
        I = single(I);
    end
    if deconOffset ~= 0 || EdgeErosion > 0
        I_mask = I ~= 0;
    end    
    I = max(I - Background, 0);    
    classI = underlyingType(I);

    J_2 = I;
    J_3 = zeros(sizeI, classI);
    J_4 = zeros(prod(sizeI), 1, classI);
end

% 3. L_R Iterations
estep = saveStep;
err_mat = zeros(NUMIT, 4);

lambda = 2*any(J_4(:)~=0);
for k = 1 : NUMIT
    
    % 3.a Make an image predictions for the next iteration
    if k > 2
        % lambda = (J_4(:,1).'*J_4(:,2))/(J_4(:,2).'*J_4(:,2) +eps);        
        if useGPU
            lambda = sum((J_2(:) - Y(:)) .* J_4) / (sum(J_4 .* J_4) + eps);
            lambda = max(min(lambda,1),0);% stability enforcement
            J_4 = J_2(:) - Y(:);
        else
            % the same, summed in double in an order that doesn't depend on
            % the number of threads
            [lambda, J_4] = decon_lucy_accel('lambda', J_2, Y, J_4);
        end
    elseif k == 2
        J_4 = J_2(:) - Y(:);
    end
    if useGPU
        Y = max(J_2 + lambda*(J_2 - J_3), 0);% plus positivity constraint
    else
        Y = decon_lucy_accel('predict', J_2, J_3, lambda);
    end
        
    % 3.b  Make core for the LR estimation
    % 3.c Determine next iteration image & apply positivity constraint
    if useMex
        % both steps in the mex, with the half OTF
        J_3 = J_2;
        J_2 = decon_lucy_rfft(useGPU, 'rlstep', Y, I, H);
    else
        % directly compute CC within the same function to reduce overhead 
        ReBlurred = max(real(ifftn(H.*fftn(Y))), eps);
        % ReBlurred = ReBlurred + (ReBlurred == 0) * eps;
        ReBlurred = cast(I, 'like', J_2)./ReBlurred;

        J_3 = J_2;
        % one FFT step per statement, and ReBlurred freed at the end, so that at
        % most two complex volumes are on the GPU at a time (the same operations)
        ReBlurred = fftn(ReBlurred);
        ReBlurred = conj(H) .* ReBlurred;
        ReBlurred = ifftn(ReBlurred);
        J_2 = max(Y .* real(ReBlurred), 0);
        clear ReBlurred;
    end
        
    if debug
        err_mat(k, 1:2) = [k, mean((double(J_2) - double(I)) .^ 2, 'all')];
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
                writetiff(uint16(gather(J_2)), sprintf('%s/Iter_%04d.tif', debug_folder, k));
            else
                writetiff(uint16(J_2), sprintf('%s/Iter_%04d.tif', debug_folder, k));
            end
        end
    end
end

if debug
    err_mat = err_mat(1 : istp, :);
end

clear J_3 J_4 Y;

% the caller's CPU work while the GPU runs the queued iterations
overlapOut = [];
if ~isempty(pr.overlapFcn)
    overlapOut = pr.overlapFcn();
end

% post-processing of deconvolved result
if dampFactor > 1
    I = cast(I, 'like', J_2);
    J_2 = J_2 .* (J_2 <= dampFactor * I) + min(J_2, dampFactor * I .* (J_2 >= dampFactor * I));
end
clear I;

if EdgeErosion > 0
    I_mask = decon_mask_edge_erosion(I_mask, EdgeErosion);
end

if scaleFactor ~= 1 || deconOffset ~= 0 || EdgeErosion > 0
    if deconOffset ~= 0 || EdgeErosion > 0
        J_2 = (J_2 * scaleFactor + deconOffset) .* I_mask;
    else
        J_2 = J_2 * scaleFactor;
    end
end

if save16bit
    J_2 = uint16(J_2);
end

if useGPU
    if ~isempty(bbox)
        J_2 = J_2(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6));
    end
    J_2 = gather(J_2);
else
    if ~isempty(bbox)
        try 
            J_2 = crop3d_mex(J_2, bbox);
        catch ME
            mex_fallback('crop3d_mex', ME, {'crop:dataTypeError'});
            J_2 = J_2(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6));
        end
    end
end

end


function [tf] = same_psf(a, b)
% whether single arrays a and b have the same size and bits (unlike isequal,
% -0 differs from 0 and a NaN equals itself); b on the CPU

tf = isequal(size(a), size(b)) && isequal(typecast(gather(a(:)), 'uint32'), typecast(b(:), 'uint32'));

end
