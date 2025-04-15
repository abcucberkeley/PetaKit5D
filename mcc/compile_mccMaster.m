function [] = compile_mccMaster(nojvm, createMCRInstaller, zipMac)
% script to compile and configure mccMaster.m

if nargin < 1
    nojvm = true;
end
if nargin < 2
    createMCRInstaller = false;
end
if nargin < 3
    zipMac = true;
end

cpath = pwd;

[fpath, fname] = fileparts(which('compile_mccMaster.m'));
cd(fpath);

if ismac
    if nojvm
        mdir = 'PetaKit5DMCC/mac';
    else
        mdir = 'PetaKit5DMCC/mac_with_jvm';
    end
    
    mccPath = ['/Applications/', mdir];
    if ~exist(mccPath, 'dir')
        mkdir_recursive(mccPath);
    end
    ioReleaseFolder = 'mac';
    if computer == "MACA64"
        ioReleaseFolder = 'macArm';
    end
    % mcc -v -R -nodisplay -C -d /Applications/PetaKit5DMCC -m mccMaster.m
    copyfile([fpath '/../microscopeDataProcessing/io/cpp-tiff/' ioReleaseFolder '/*.dylib'], ['/Applications/', mdir, '/']);
    copyfile([fpath '/../microscopeDataProcessing/io/cpp-zarr/' ioReleaseFolder '/*.dylib'], ['/Applications/', mdir, '/']);
    if nojvm
        mcc -v -R -nodisplay -R -nojvm -C -d /Applications/PetaKit5DMCC/mac -m mccMaster.m
    else
        mcc -v -C -d /Applications/PetaKit5DMCC/mac_with_jvm -m mccMaster.m
    end
    cd(mccPath);
    if nojvm
        system("sed -i '' '20i\'$'\n''  DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:/Applications/PetaKit5DMCC/mac; '$'\n''' run_mccMaster.sh");
    else
        system("sed -i '' '20i\'$'\n''  DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:/Applications/PetaKit5DMCC/mac_with_jvm; '$'\n''' run_mccMaster.sh");
    end
    system("sed -i '' '20i\'$'\n''  # add custom library paths'$'\n''' run_mccMaster.sh");
    % disable printout of environment variables
    system('sed -i '''' ''9s/^/#/'' run_mccMaster.sh');
    system('sed -i '''' ''14s/^/#/'' run_mccMaster.sh');
    system('sed -i '''' ''16s/^/#/'' run_mccMaster.sh');
    system('sed -i '''' ''23s/^/#/'' run_mccMaster.sh');
    % create a zip file in the repo folder
    cd(fpath);
    if ~exist('mac', 'dir')
        mkdir('mac');
    end
    if zipMac
        zip('mac/PetaKit5DMCC.zip', '/Applications/PetaKit5DMCC/*');
    end
elseif isunix
    if nojvm
        mdir = 'linux';
    else
        mdir = 'linux_with_jvm';
    end
    if ~exist(mdir, 'dir')
        mkdir(mdir);
    end
    % mcc -v -R -nodisplay -R -singleCompThread -d linux -m mccMaster.m
    if nojvm
        mcc -v -R -nodisplay -R -nojvm -d linux -m mccMaster.m
    else
        mcc -v -d linux_with_jvm -m mccMaster.m
    end
    
    cd(mdir);
    system('sed -i ''21i\  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${exe_dir}/../../microscopeDataProcessing/io/cpp-tiff/linux;'' run_mccMaster.sh');
    system('sed -i ''21i\  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${exe_dir}/../../microscopeDataProcessing/io/cpp-zarr/linux;'' run_mccMaster.sh');
    system('sed -i ''21i\  # add custom library paths'' run_mccMaster.sh');
    % disable printout of environment variables
    system('sed -i ''9s/^/#/'' run_mccMaster.sh');
    system('sed -i ''14s/^/#/'' run_mccMaster.sh');
    system('sed -i ''16s/^/#/'' run_mccMaster.sh');
    system('sed -i ''25s/^/#/'' run_mccMaster.sh');
elseif ispc
    if nojvm
        mdir = 'windows';
    else
        mdir = 'windows_with_jvm';
    end
    if ~exist(mdir, 'dir')
        mkdir(mdir);
    end
    if nojvm
        mcc -v -R -nojvm -d windows -m mccMaster.m
    else
        mcc -v -d windows_with_jvm -m mccMaster.m
    end
    cd(mdir);
    copyfile('../../microscopeDataProcessing/io/cpp-tiff/windows/*dll', './');
    copyfile('../../microscopeDataProcessing/io/cpp-zarr/windows/*dll', './');
end

if createMCRInstaller
    cd(fpath);
    if ismac
        mdir = ['/Applications/', mdir];
    end
    compiler.runtime.customInstaller("installMCR", [mdir '/requiredMCRProducts.txt'], RuntimeDelivery="web");
end

cd(cpath);

end
