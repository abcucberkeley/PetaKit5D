%% script to compile and configure mccMaster.m

cpath = pwd;

[fpath, fname] = fileparts(which('compile_mccMaster.m'));
cd(fpath);

if ismac
    if ~exist('/Applications/LLSM5DToolsMCC', 'dir')
        mkdir('/Applications/LLSM5DToolsMCC');
    end
    mcc -v -R -nodisplay -C -d /Applications/LLSM5DToolsMCC -m mccMaster.m

    cd('/Applications/LLSM5DToolsMCC');
    system("sed -i '' '20i\'$'\n''  DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${exe_dir}/../../microscopeDataProcessing/io/cpp-tiff/mac; '$'\n''' run_mccMaster.sh");
    system("sed -i '' '20i\'$'\n''  DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${exe_dir}/../../microscopeDataProcessing/io/cpp-zarr/max; '$'\n''' run_mccMaster.sh");
    system("sed -i '' '20i\'$'\n''  # add custom library paths'$'\n''' run_mccMaster.sh");
    % create a zip file in the repo folder
    cd(fpath);
    if ~exist('mac', 'dir')
        mkdir('mac');
    end
    zip('mac/LLSM5DToolsMCC.zip', '/Applications/LLSM5DToolsMCC/*')
elseif isunix
    if ~exist('linux', 'dir')
        mkdir('linux');
    end
%     mcc -v -R -nodisplay -R -singleCompThread -d linux -m mccMaster.m
    mcc -v -R -nodisplay -R -nojvm -d linux -m mccMaster.m

    cd('linux');
    system('sed -i ''21i\  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${exe_dir}/../../microscopeDataProcessing/io/cpp-tiff/linux;'' run_mccMaster.sh');
    system('sed -i ''21i\  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${exe_dir}/../../microscopeDataProcessing/io/cpp-zarr/linux;'' run_mccMaster.sh');
    system('sed -i ''21i\  # add custom library paths'' run_mccMaster.sh');
elseif ispc
    if ~exist('windows', 'dir')
        mkdir('windows');
    end
    mcc -v -d windows -m mccMaster.m
    
    cd('windows');
    copyfile('../../microscopeDataProcessing/io/cpp-tiff/windows/*dll', './');
    copyfile('../../microscopeDataProcessing/io/cpp-zarr/windows/*dll', './');
end

cd(cpath);
