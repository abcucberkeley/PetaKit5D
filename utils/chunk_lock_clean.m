function [out_stat, deleted_files] = chunk_lock_clean(path, duration, initial_wait)
% This function is for cleaning tmp files whose jobs are not running in the folder.
%

% Author: Xiongtao Ruan
% Date: Feb. 24, 2016
% change datetime('now') to clock to be compatible with 8.1
if nargin < 3
    initial_wait = false;
end
if nargin < 2
    delete_by_time = false;
else
    delete_by_time = true;
end
hard_duration = 60;
deleted_files = {};  
out_stat = false;
keep_flag = true;

if initial_wait
    job_id = regexprep(getenv('SLURM_JOB_ID'), '[\r\n]', '');    
    wait_time_str = fliplr(job_id);
    if ~isempty(wait_time_str)        
        wait_time_str = [wait_time_str(1), '.', wait_time_str(2:end)];
    end
    wait_time = str2double(wait_time_str);
    if ~isnumeric(wait_time) || isnan(wait_time)
        wait_time = rand() * 5;
    end
    pause(wait_time + rand());
end

[status, output] = system(sprintf('find %s -type f -name "*.tmp"', path));
if ~isempty(output)
    % temp_filenames = strsplit(strtrim(output), '\n');
    % use more robust way
    temp_filenames = strsplit(strtrim(output), 'tmp'); 
    
    if numel(temp_filenames) > 0 
        for i = 1 : numel(temp_filenames)
            current_temp_filename = [strtrim(temp_filenames{i}), 'tmp'];
            if delete_by_time
                temp_file_info = dir(current_temp_filename);
                %if (datenum(datetime('now')) - [temp_file_info.datenum]) * 24 * 60 > duration
                if (datenum(clock) - [temp_file_info.datenum]) * 24 * 60 > duration
                    deleted_files{end + 1} = current_temp_filename;
                    delete(current_temp_filename);
                end
            elseif exist(current_temp_filename, 'file')
                % current_temp_filename
                try
                    job_info = textread(current_temp_filename, '%s', 'delimiter', '\n');
                catch ME
                    disp(ME);
                    if ~exist(current_temp_filename, 'file')
                        continue;
                    end
                end
                if numel(job_info) ~=2
                    deleted_files{end + 1} = current_temp_filename;
                    delete(current_temp_filename);
                    keep_flag = false;
                    continue
                end
                id_type = job_info{1};
                id = job_info{2};
                 
                % for the cluster to check job id
                if strcmp(id_type, 'jid')
                    [status, output] = system(sprintf('sacct -j %s --format=State', id));           
                
                    if isempty(regexp(output, 'RUNNING'))  
                        deleted_files{end + 1} = current_temp_filename;
                        delete(current_temp_filename);
                        keep_flag = false;
                    end
                elseif strcmp(id_type, 'pid')
                    [status, output] = system(sprintf('ps -o pid= -p %s', id));  
                    if isempty(output)
                        deleted_files{end + 1} = current_temp_filename;
                        delete(current_temp_filename);
                        keep_flag = false;
                    end
                else
                    deleted_files{end + 1} = current_temp_filename;
                    delete(current_temp_filename);
                    keep_flag = false;			
                end 
		% add hard clean condition, that is, if the tmp file exist for a certain time(here set as 1hr)
                if keep_flag
                    temp_file_info = dir(current_temp_filename);
                    if (datenum(clock) - [temp_file_info.datenum]) * 24 * 60 > hard_duration
                            deleted_files{end + 1} = current_temp_filename;
                            delete(current_temp_filename);
                    end
                end
            end
        end
    end
end

if ~isempty(deleted_files)
    out_stat = true;
end

end
    

