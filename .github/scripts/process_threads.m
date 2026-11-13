function [n, cpus] = process_threads()
% the threads of this MATLAB process now (NaN if unknown) and the CPUs it may
% run on (Windows: the affinity mask; Linux: Cpus_allowed_list), for the
% diagnostics of the CI job exit-probe (mex_exit_probe, suite_exit_probe)

n = NaN;
cpus = '';
pid = feature('getpid');
try
    if ispc
        [s, out] = system(sprintf(['powershell -NoProfile -NonInteractive -Command "$p = Get-Process -Id %d; ', ...
            '''{0} {1}'' -f $p.Threads.Count, [int64]$p.ProcessorAffinity"'], pid));
        if s == 0
            v = sscanf(strtrim(out), '%d %d');
            n = v(1);
            cpus = sprintf('mask 0x%X', v(2));
        end
    elseif ismac
        [s, out] = system(sprintf('ps -M -p %d | tail -n +2 | wc -l', pid));
        if s == 0
            n = str2double(strtrim(out));
        end
    else
        st = fileread(sprintf('/proc/%d/status', pid));
        n = str2double(regexp(st, 'Threads:\s*(\d+)', 'tokens', 'once'));
        cpus = char(regexp(st, 'Cpus_allowed_list:\s*(\S+)', 'tokens', 'once'));
    end
catch
end

end
