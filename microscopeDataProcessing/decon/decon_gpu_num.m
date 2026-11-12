function [gpuNum] = decon_gpu_num(required)
% the number of GPUs the deconvolution can use (usable_gpu_num: GPUs newer
% than this MATLAB supports get forward compatibility). With none, an error
% (decon:noUsableGPU) if a GPU is required, else a warning that the
% deconvolution runs on the CPU. Both give the number of GPUs found,
% CUDA_VISIBLE_DEVICES and, when GPUs are found but none can be used (broken,
% in use, or not supported by this MATLAB), MATLAB's reason.
%
% required: true|false (default: false).


if nargin < 1
    required = false;
end

[gpuNum, reason, gpuAll] = usable_gpu_num();
if gpuNum > 0
    return;
end

if ~isempty(reason)
    reason = [' ', reason];
end
msg = sprintf('No GPU is available (%d found, CUDA_VISIBLE_DEVICES=''%s'')', gpuAll, getenv('CUDA_VISIBLE_DEVICES'));
if required
    error('decon:noUsableGPU', '%s, but the deconvolution is to run on a GPU (useGPU true).%s', msg, reason);
end
warning('%s, so the deconvolution runs on the CPU.%s', msg, reason);

end
