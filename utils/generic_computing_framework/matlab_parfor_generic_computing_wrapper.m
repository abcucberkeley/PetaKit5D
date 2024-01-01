function [is_done_flag] = matlab_parfor_generic_computing_wrapper(inputFullpaths, outputFullpaths, funcStrs, varargin)
% generic computing framework that use a function handle/string as input
% for the computing based on matlab parfor
% 
% 
% Author: Xiongtao Ruan (10/30/2021)
% 
% xruan (11/03/2022): limit workers threads to maxNumCompThreads / nworker

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('inputFullpaths', @(x) iscell(x) || ischar(x));
ip.addRequired('outputFullpaths', @(x) iscell(x) || ischar(x));
ip.addRequired('functionStrs', @(x) iscell(x) || ischar(x));
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('GPUJob', true, @islogical); % master node participate in the task computing. 
ip.addParameter('jobsPerGPU', 1, @isnumeric); % number of jobs per GPU
ip.addParameter('nworker', 24, @isnumeric); % number of workers 
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('tmpDir', '', @ischar);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('taskBatchNum', 1, @isnumeric); % aggragate several tasks together
ip.addParameter('language', 'matlab', @ischar); % support matlab, bash

ip.parse(inputFullpaths, outputFullpaths, funcStrs, varargin{:});

% move to the root path
paths = split(which('matlab_parfor_generic_computing_wrapper'), 'XiongtaoScripts');
cd(paths{1});

pr = ip.Results;
GPUJob = pr.GPUJob;
jobsPerGPU = pr.jobsPerGPU;
nworker = pr.nworker;
taskBatchNum = pr.taskBatchNum;
uuid = pr.uuid;
language = pr.language;

if isempty(uuid)
    uuid = 'tmp';
end

[dataPath, ~] = fileparts(inputFullpaths{1});

% if the job is on GPU, count ther number of GPU available and override
% number of batches to the number of GPUs
if GPUJob
    nGPU = gpuDeviceCount("available");
else 
    nGPU = 0;
end

if GPUJob 
    nworker = nGPU * jobsPerGPU;
end

persistent c;
NumThreads = maxNumCompThreads;
if isempty(gcp('nocreate'))
    c = parcluster('local');
    c.NumThreads = max(1, round(maxNumCompThreads / nworker));
    parpool(c, nworker);
    if GPUJob 
        spmd
            if gpuDeviceCount("available") > 0
                gd = gpuDevice(rem(labindex - 1, nGPU) + 1);                
                idx = gd.Index;
                disp(['Using GPU ',num2str(idx)]);                
            end
            setup([],true);
        end
    else
        spmd
            setup([],true);
            maxNumCompThreads(max(1, round(NumThreads / nworker)));
            maxNumCompThreads
        end
    end
end

loop_counter = 0;
max_loop = 3;

output_exist_mat = batch_file_exist(outputFullpaths, [], true);
is_done_flag = output_exist_mat;

tmpFullpaths = cellfun(@(x) [x, '_', uuid], outputFullpaths, 'unif', 0);
while ~all(is_done_flag) && loop_counter < max_loop
    % check output files
    if loop_counter > 0
        output_exist_mat = batch_file_exist(outputFullpaths, [], true);
        is_done_flag(~is_done_flag) = output_exist_mat;
    end
    if all(is_done_flag)
        break;
    end
    outputFullpaths = outputFullpaths(~output_exist_mat);
    tmpFullpaths = tmpFullpaths(~output_exist_mat);
    funcStrs = funcStrs(~output_exist_mat);
    nF = sum(~output_exist_mat);

    parfor f = 1 : nF
        if exist(outputFullpaths{f}, 'file') 
            continue;
        end
        if exist(tmpFullpaths{f}, 'file')
            continue;
        else
            fclose(fopen(tmpFullpaths{f}, 'w'));
        end
        func_str = funcStrs{f};  
        switch language
            case 'matlab'
                try
                    feval(str2func(['@()', func_str]));
                    % eval(func_str);
                catch ME
                    pause(1);
                    try 
                        feval(str2func(['@()', func_str]));
                    catch ME_1
                        pause(2);
                        try
                            feval(str2func(['@()', func_str]));  
                        catch ME_2
                            disp(ME_2);
                        end
                        disp(ME_1)
                    end
                    disp(ME);
                end
            case 'bash'
                gpu_str = '';
                if GPUJob
                    t = getCurrentTask();              
                    idx = rem(t.ID - 1, nGPU) + 1;
                    gpu_str = sprintf('CUDA_VISIBLE_DEVICES=%s', num2str(idx-1));
                end
                fprintf(gpu_str);
                t0=tic; system([gpu_str, ' ', func_str], '-echo'); toc(t0)
        end
        if exist(outputFullpaths{f}, 'file')
            delete(tmpFullpaths{f});
        end
    end

    loop_counter = loop_counter + 1;
end


% maxQueue = nworker + 1;
% 
% out = cell(nF, 1);
% s = 1;
% for f =1 : nF
%     if exist(outputFullpaths{f}, 'file')
%         s = f + 1;
%         continue;
%     end
%     
%     func_str = funcStrs{f};
%     % if GPUJob
%     %     gpuDevice(rem(f, nGPU) + 1);
%     % end
%     try
%         out{f} = parfeval(str2func(['@()', func_str]), 0);
%     catch ME
%     end
% %     feval(str2func(['@()', func_str]));
%     numQueued = cellfun(@(x) strcmp(x.State, 'queued'), out(s : f));
%     while sum(numQueued) > maxQueue
%         pause(5);
%         numQueued = cellfun(@(x) strcmp(x.State, 'queued'), out(s : f));
%     end
%     if any(numQueued == 0)
%         s = s + find(numQueued == 0, 1, 'last') - 1;
%     end
% end

% wait(out);

end

