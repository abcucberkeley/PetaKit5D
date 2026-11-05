function [gpuNum] = decon_gpu_num()
% the number of GPUs the deconvolution can use (usable_gpu_num: GPUs newer
% than this MATLAB supports get forward compatibility). With none, a warning
% that the deconvolution runs on the CPU, with the number of GPUs found,
% CUDA_VISIBLE_DEVICES and, when GPUs are found but none can be used (broken,
% in use, or not supported by this MATLAB), MATLAB's reason.


[gpuNum, reason, gpuAll] = usable_gpu_num();
if gpuNum > 0
    return;
end

if ~isempty(reason)
    reason = [' ', reason];
end
warning('No GPU is available (%d found, CUDA_VISIBLE_DEVICES=''%s''), so the deconvolution runs on the CPU.%s', ...
    gpuAll, getenv('CUDA_VISIBLE_DEVICES'), reason);

end
