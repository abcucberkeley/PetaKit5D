function [gpuNum, reason, gpuAll] = usable_gpu_num()
% the number of GPUs MATLAB can use (gpuDeviceCount("available")). If GPUs are
% found but none can be used because they are newer than this MATLAB release
% supports natively (e.g., B200 GPUs, compute capability 10.0, in R2024b),
% forward compatibility is turned on for the session, unless the environment
% variable MW_CUDA_FORWARD_COMPATIBILITY is 0: the CUDA driver then recompiles
% MATLAB's GPU libraries for them (once; the CUDA cache keeps them for later
% sessions). MathWorks cautions that GPU results may change with forward
% compatibility; the deconvolution's on B200 GPUs match the CPU's as on
% natively supported GPUs.
%
% reason: why no GPU can be used, empty if one can or none is found (MATLAB's
% error of selecting GPU 1; nothing is selected when that fails).
% gpuAll: the number of GPUs found (gpuDeviceCount("all")).


gpuNum = 0;
reason = '';
gpuAll = 0;
try
    gpuNum = gpuDeviceCount("available");
    gpuAll = gpuDeviceCount("all");
catch ME
    reason = sprintf('Counting the GPUs failed: %s', ME.message);
    return;
end
if gpuNum > 0 || gpuAll == 0
    return;
end

try
    gpuDevice(1);
catch ME
    reason = sprintf('Selecting GPU 1 failed: %s', strtrim(regexprep(ME.message, '\s+', ' ')));
    if ~strcmp(ME.identifier, 'parallel:gpu:device:DeviceTooNew') || strcmp(strtrim(getenv('MW_CUDA_FORWARD_COMPATIBILITY')), '0')
        return;
    end
    try
        parallel.gpu.enableCUDAForwardCompatibility(true);
        gpuNum = gpuDeviceCount("available");
    catch ME
        gpuNum = 0;
        reason = sprintf('%s Turning on forward compatibility failed: %s', reason, ME.message);
    end
    if gpuNum > 0
        reason = '';
        fprintf(['The GPUs (%d) are newer than MATLAB R%s supports natively, so forward compatibility is ' ...
            'turned on for them (MW_CUDA_FORWARD_COMPATIBILITY=0 turns this off).\n'], gpuNum, version('-release'));
        return;
    end
    % it didn't make a GPU usable: as it was
    try
        parallel.gpu.enableCUDAForwardCompatibility(false);
    catch
    end
end

end
