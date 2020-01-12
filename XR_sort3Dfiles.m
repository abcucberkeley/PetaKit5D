function [] = XR_sort3Dfiles(data, channel_names, channel_patterns)
% The function is to reorganize images, mainly create channel
% subdirectories under an experiment directory, if the images are in the
% experiment directory. 
% 
% example: XR_sort3Dfiles(data, {'ch1_560', 'ch2_642'}, {'*560nm*.tif', '*642nm*.tif'})
% 
% 09/16/2019 xruan: copied from GU_sort3Dfiles.m and modify it as a
% function.



for k = 1:numel(data)
    tic
    rt = [data(k).source];
    % cd(rt)
    
    for i = 1 : numel(channel_names)
        cur_channel_name = channel_names{i};
        mkdir([rt cur_channel_name]);       
        
        cur_channel_pattern = channel_patterns{i};
        
        stackFiles = dir([rt cur_channel_pattern]);  
        if numel(stackFiles) >0
            stackFiles = {stackFiles.name};

            i = regexpi(stackFiles, '(?<=msec_)\d+', 'match');
            [~,i] = sort(cellfun(@str2double, i));
            stackFiles = stackFiles(i);
            stackFiles = cellfun(@(i) [i], stackFiles, 'unif', 0);
            % cd(rt)

            for j = 1:numel(stackFiles)
                movefile([rt stackFiles{j}], [rt cur_channel_name filesep stackFiles{j}]) 
            end
        end
    end
    toc
    
end


end