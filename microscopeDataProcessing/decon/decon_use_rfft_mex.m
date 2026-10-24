function [useMex, canUseMex] = decon_use_rfft_mex(useMex, method, I, useGPU, debug, Background)
% whether the RL iterations of decon_lucy_function (method 'simplified') and
% decon_lucy_omw_function ('omw') run in the mex with half-spectrum
% (real-to-complex) FFTs (decon_lucy_rfft_mex on the CPU,
% decon_lucy_rfft_gpu_mex on the GPU) instead of the MATLAB code. true or
% false forces either path (true errors if the mex can't be used). Empty
% (the default) uses the mex for the methods in defaultMethodsCPU or
% defaultMethodsGPU if it can: canUseMex is true if the mex exists for the
% platform and supports the call (real 2D or 3D data with at least 2 elements
% along y and x on the CPU, any class but double on the CPU, uint16 or single
% on the GPU, a scalar background, debug output only on the CPU).
%
% The CPU mex gives the same results as the MATLAB code on Linux, bit for bit
% (it runs the same FFTW transforms as MATLAB's fftn and ifftn there, only on
% the half spectrum). The GPU mex uses other cuFFT transforms than MATLAB's gpuArray
% FFTs, which changes 1 in a small fraction of the uint16 voxels (demo tile:
% OMW-2 0.004%, RL-30 0.01%; as much as the MATLAB code changes between CPU
% and GPU).


% methods that use the mex by default: on the CPU both on Linux, where it
% gives the same results as the MATLAB code (on macOS MATLAB's FFTs differ
% from it in the last bits, and Windows isn't checked), on the GPU none
defaultMethodsCPU = {};
if isunix && ~ismac
    defaultMethodsCPU = {'omw', 'simplified'};
end
defaultMethodsGPU = {};

if useGPU
    mexName = 'decon_lucy_rfft_gpu_mex';
    supported = any(strcmp(underlyingType(I), {'uint16', 'single'}));
    defaultMethods = defaultMethodsGPU;
else
    mexName = 'decon_lucy_rfft_mex';
    supported = ~isa(I, 'double') && size(I, 1) > 1 && size(I, 2) > 1;
    defaultMethods = defaultMethodsCPU;
end
supported = supported && isreal(I) && ~isempty(I) && ndims(I) <= 3 && isscalar(Background) ...
    && isreal(Background) && ~(debug && useGPU);
available = exist(mexName, 'file') == 3;
canUseMex = supported && available;

if isempty(useMex)
    useMex = canUseMex && any(strcmp(method, defaultMethods));
    return;
end
if useMex && ~canUseMex
    if ~available
        error('%s is not available for this platform.', mexName);
    end
    error('%s does not support the data (class %s, size %s) or the options (debug on the GPU or a non-scalar background).', ...
        mexName, underlyingType(I), mat2str(size(I)));
end

end
