function e2e_cluster_workflow(workflow, dataPath, configFile, mccMode, timingFile)
% One end-to-end workflow in cluster mode (or locally with configFile ''), as MATLAB runs it,
% for the speed comparison with the port (python/benchmarks/e2e/e2e_cluster_workflow.py does the
% same steps with the same parameters). Appends "<workflow> <seconds>" to timingFile.
%
% workflow:
%   generic  demo_generic_computing_framework's 16 tasks (crop and deskew/rotate each tile)
%   dsr      demo_geometric_transformation's dataset deskew/rotation (16 tiles, zarr results)
%   dsr_large  the large-file deskew/rotation of the zarr volumes in dataPath (two channels
%            CamA_ch0 and CamB_ch0; xyPixelSize 0.098, dz 0.5)
% mccMode: true for compiled jobs (the config names the compiled master)

parseCluster = ~isempty(configFile);
if ~parseCluster
    configFile = '';
end
if dataPath(end) ~= '/'
    dataPath = [dataPath, '/'];
end
t0 = tic;
switch workflow
    case 'generic'
        dir_info = dir([dataPath, '*.tif']);
        fsns = {dir_info.name};
        inputFullnames = cellfun(@(x) [dataPath, x], fsns, 'UniformOutput', false);
        outPath = [dataPath, 'Cropped_DSR/'];
        mkdir(outPath);
        outputFullnames = cellfun(@(x) [outPath, x], fsns, 'UniformOutput', false);
        inBbox = [101, 1, 1, 1700, 512, 500];
        inBbox(4 : 6) = min(inBbox(4 : 6), getImageSize(inputFullnames{1}));
        nF = numel(inputFullnames);
        funcStrs = cell(nF, 1);
        for f = 1 : nF
            funcStrs{f} = sprintf(['crop_deskew_rotate_demo_function(''%s'',''%s'',%d,', ...
                '%d,%d,%s)'], inputFullnames{f}, outputFullnames{f}, 32.45, 0.108, 0.3, ...
                strrep(mat2str(inBbox), ' ', ','));
        end
        generic_computing_frameworks_wrapper(inputFullnames, outputFullnames, funcStrs, ...
            parseCluster=parseCluster, masterCompute=true, memAllocate=20, ...
            configFile=configFile, mccMode=mccMode);
    case 'dsr'
        XR_deskew_rotate_data_wrapper({dataPath}, Deskew=true, Rotate=true, ...
            DSRCombined=true, xyPixelSize=0.108, dz=0.3, SkewAngle=32.45, ...
            ObjectiveScan=false, Reverse=true, ChannelPatterns={'CamA', 'CamB'}, ...
            largeFile=false, zarrFile=false, saveZarr=true, blockSize=[256, 256, 256], ...
            Save16bit=true, parseCluster=parseCluster, masterCompute=true, ...
            configFile=configFile, mccMode=mccMode);
    case 'dsr_large'
        XR_deskew_rotate_data_wrapper({dataPath}, xyPixelSize=0.098, dz=0.5, ...
            Reverse=true, ChannelPatterns={'CamA_ch0', 'CamB_ch0'}, largeFile=true, ...
            zarrFile=true, saveZarr=true, Save16bit=true, parseCluster=parseCluster, ...
            masterCompute=true, configFile=configFile, mccMode=mccMode);
    otherwise
        error('e2e:workflow', 'unknown workflow %s', workflow);
end
t = toc(t0);
fid = fopen(timingFile, 'a');
fprintf(fid, '%s %.3f\n', workflow, t);
fclose(fid);
end
