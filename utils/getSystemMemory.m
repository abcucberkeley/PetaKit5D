function [totalMem] = getSystemMemory()
% get system memory in GB
% For Linux, also check SLURM environmental variable for actual available memory
%
% Author: Xiongtao Ruan (08/06/2024)


% NaN if it can't be found, so that callers take the low memory path
totalMem = NaN;
if ismac
    [status, output] = system('sysctl -n hw.memsize');
    totalMem = str2double(strip(output)) / 1024^3;
elseif isunix
    job_id = getenv('SLURM_JOB_ID');
    slurmJob = ~isempty(job_id);

    if slurmJob
        % the limit of the job's cgroup: a job submitted from another job
        % (e.g., one that coordinates jobs) inherits its SLURM_MEM_PER_NODE
        totalMem = job_cgroup_memory(job_id);
        if isnan(totalMem) && isempty(getenv('SLURM_MEM_PER_CPU'))
            % the job's memory per node in MB (unset with --mem-per-cpu; with
            % both set, one is another job's, so scontrol tells)
            totalMem = str2double(getenv("SLURM_MEM_PER_NODE")) / 1024;
        end
        if isnan(totalMem)
            array_id = getenv('SLURM_ARRAY_TASK_ID');
            if isempty(array_id)
                job_array_id = sprintf('%s', job_id);
            else
                job_array_id = sprintf('%s_%s', job_id, array_id);
            end
            [status, output] = system(sprintf('scontrol show job %s', job_array_id));
            pattern = 'mem=([\d.]+[TGMK])';
            matches = regexp(output, pattern, 'tokens');
            if ~isempty(matches)
                totalMem = matches{1}{1};
                switch totalMem(end)
                    case 'T'
                        unit_denominator = 1 / 1024;
                    case 'G'
                        unit_denominator = 1;
                    case 'M'
                        unit_denominator = 1024;                        
                    case 'K'
                        unit_denominator = 1024^2;
                    otherwise
                        error('Unknown unit for %s', totalMem);
                end
                totalMem = str2double(totalMem(1 : end - 1)) / unit_denominator;
            else
                disp('No match found');
            end
        end
    else
        [status, output] = system('cat /proc/meminfo | grep MemTotal');
        pattern = 'MemTotal:\s+(\d+)\s+kB';
        matches = regexp(output, pattern, 'tokens');
        if ~isempty(matches)
            totalMem = str2double(matches{1}{1});
            totalMem = totalMem / 1024^2;
        else
            disp('No match found');
        end
    end
elseif ispc
    a = memory;
    totalMem = a.MaxPossibleArrayBytes / 1024^3;
end

end


function [mem] = job_cgroup_memory(job_id)
% the smallest memory limit (GB) of the cgroups on this process's cgroup path
% (cgroup v1: memory.limit_in_bytes of the memory controller; v2: memory.max),
% if the path is that of the Slurm job job_id (job_<id>); NaN otherwise

mem = NaN;
try
    txt = fileread('/proc/self/cgroup');
catch
    return;
end
p = regexp(txt, '^\d+:[^:\n]*\<memory\>[^:\n]*:(\S+)$', 'tokens', 'once', 'lineanchors');
if ~isempty(p)
    root = '/sys/fs/cgroup/memory';
    limitName = 'memory.limit_in_bytes';
else
    p = regexp(txt, '^0::(\S+)$', 'tokens', 'once', 'lineanchors');
    root = '/sys/fs/cgroup';
    limitName = 'memory.max';
end
if isempty(p) || isempty(regexp(p{1}, ['/job_', job_id, '(/|$)'], 'once'))
    return;
end
cgPath = p{1};
limit = inf;
while ~isempty(cgPath) && ~strcmp(cgPath, '/')
    try
        v = str2double(strtrim(fileread([root, cgPath, '/', limitName])));
        if v > 0
            limit = min(limit, v);
        end
    catch
    end
    cgPath = fileparts(cgPath);
end
% (cgroup v1 shows no limit as a huge number)
if limit < 2^60
    mem = limit / 2^30;
end

end

