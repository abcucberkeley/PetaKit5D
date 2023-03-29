%% script to compile and configure mccMaster.m

[fpath, fname] = fileparts(which('compile_mccMaster.m'));
cd(fpath);

if ismac
    if ~exist('/Applications/LLSM5DToolsMCC', 'dir')
        mkdir('/Applications/LLSM5DToolsMCC');
    end
    mcc -v -R -nodisplay -C -d /Applications/LLSM5DToolsMCC -m mccMaster.m

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
    if ~exist('windows', 'dir')
        mkdir('windows');
    end
    mcc -v -d windows -m mccMaster.m
    
    cd('windows');
    copyfile('../../microscopeDataProcessing/io/c-zarr/parallelWriteZarr/windows/*dll', './');
end

