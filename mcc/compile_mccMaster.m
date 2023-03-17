%% script to compile and configure mccMaster.m

[fpath, fname] = fileparts(which(mfilename));
cd(fpath);

mcc -v -R -nodisplay -m mccMaster.m

%% add custom library paths
system('sed -i ''21i\  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${exe_dir}/../LLSM5DTools/microscopeDataProcessing/io/c-zarr/parallelWriteZarr/linux;'' run_mccMaster.sh');
system('sed -i ''21i\  # add custom library paths'' run_mccMaster.sh');
