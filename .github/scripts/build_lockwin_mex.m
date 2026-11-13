function [] = build_lockwin_mex(outDir, mingwDir)
% for the CI job exit-probe: build the PetaKit5D OpenMP mex files of Windows
% into outDir from copies of their sources whose mexLock block also applies on
% Windows. Kept as the name of the locked build; build_variant_mex is the one
% builder (it also builds the mexAtExit variants).

build_variant_mex(outDir, mingwDir, 'lock');

end
