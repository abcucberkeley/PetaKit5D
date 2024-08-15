function [totalMem] = getSystemMemory()
% get system memory in GB
% For Linux, also check SLURM environmental variable for actual available memory
%
% Author: Xiongtao Ruan (08/06/2024)


if ismac
    [status, output] = system('sysctl -n hw.memsize');
    totalMem = str2double(strip(output)) / 1024^3;
elseif isunix
    job_id = getenv('SLURM_JOB_ID');
    slurmJob = ~isempty(job_id);

    if slurmJob
        totalMem = getenv("SLURM_MEM_PER_NODE");
        if isempty(totalMem)
            array_id = getenv('SLURM_ARRAY_TASK_ID');
            if isempty(array_id)
                job_array_id = sprintf('%s', job_id);
            else
                job_array_id = sprintf('%s_%s', job_id, array_id);
            end
            [status, output] = system(sprintf('scontrol show job %s', job_array_id));
            pattern = 'mem=(\d+[GMK]),';
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

