%% script to compile and configure mccMaster.m

[fpath, fname] = fileparts(which(mfilename));
cd(fpath);

if ismac
    tmp = matlab.desktop.editor.getActive;
    cd(fileparts(tmp.Filename));
    if ~exist('/Applications/LLSM5DTools', 'dir')
        mkdir('/Applications/LLSM5DTools');
    end
    mcc -v -R -nodisplay -C -d /Applications/LLSM5DTools -m mccMaster.m

    cd('mac');
    system("sed -i '' '20i\'$'\n''  DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${exe_dir}/../../microscopeDataProcessing/io/c-zarr/parallelWriteZarr/mac; '$'\n''' run_mccMaster.sh");
    system("sed -i '' '20i\'$'\n''  # add custom library paths'$'\n''' run_mccMaster.sh");
elseif isunix
    if ~exist('linux', 'dir')
        mkdir('linux');
    end
    mcc -v -R -nodisplay -d linux -m mccMaster.m

    cd('linux');
    system('sed -i ''21i\  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${exe_dir}/../../microscopeDataProcessing/io/c-zarr/parallelWriteZarr/linux;'' run_mccMaster.sh');
    system('sed -i ''21i\  # add custom library paths'' run_mccMaster.sh');
elseif ispc
    tmp = matlab.desktop.editor.getActive;
    cd(fileparts(tmp.Filename));

    if ~exist('windows', 'dir')
        mkdir('windows');
    end
    mcc -v -R -nodisplay -d windows -m mccMaster.m
    
    cd('windows');
end

