%PSF_folder = ['/clusterfs/fiona/matthewmueller/oldSim3D/data/'];
%PSF_file = ['PSF_5phase.tif'];

%Load the PSF data
PSF=double(loadtiff('/clusterfs/fiona/Data/20210918_latticeSIM/PSFs/phasePSF/cropped/DS/RAW_3phaseWFPSF_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0177990936msecAbs_000x_000y_000z_0000t_combined.tif'));

O = sim_PSFtoOTF_gen(PSF,'lattice_period', 1.4, 'phase_step', .4667, 'norders', 3, 'nphases', 3, 'useGPU', false);

data = '/clusterfs/fiona/Data/20210918_latticeSIM/exp01_3phase/DS/RAW_exp01_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0013904967msecAbs_000x_000y_000z_0000t_combined.tif';
output = simReconFrame(data, O, 'lattice_period', 1.4, 'phase_step', .4667, 'norders', 3, 'nphases', 3, 'Overlap', 100, 'ChunkSize', [250,250,250], 'edgeTaper', true, 'perdecomp', true);
%dataFull = loadtiff(data);
%output = simRecon(dataFull, O, 'lattice_period', 1.4, 'phase_step', .4667, 'norders', 3, 'nphases', 3, 'edgeTaper', true, 'perdecomp', true, 'useGPU', false);

writetiff(single(output.*(output>=0)),'/clusterfs/fiona/Data/20210918_latticeSIM/exp01_3phase/DS/sim_recon/RAW_exp01_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0013904967msecAbs_000x_000y_000z_0000t_OL100_CS250.tif');








%%
output = simReconFrame(data, O, 'lattice_period', 1.4, 'phase_step', .4667, 'norders', 3, 'nphases', 3, 'Overlap', 32, 'ChunkSize', [128,128,128], 'edgeTaper', false, 'perdecomp', false);
writetiff(single(output.*(output>=0)),'/clusterfs/fiona/Data/20210918_latticeSIM/exp01_3phase/sim_recon/RAW_exp01_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0013904967msecAbs_000x_000y_000z_0000t_OL32_CS128_noedge.tif');



output = simRecon(dataFull, O, 'lattice_period', 1.4, 'phase_step', .4667, 'norders', 3, 'nphases', 3, 'edgeTaper', false, 'perdecomp', false, 'useGPU', false);
writetiff(single(output.*(output>=0)),'/clusterfs/fiona/Data/20210918_latticeSIM/exp01_3phase/sim_recon/RAW_exp01_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0013904967msecAbs_000x_000y_000z_0000t_CPU_Nochunk_noedge.tif');
%%
PSF=double(loadtiff('/clusterfs/fiona/Data/20210923_Aang_latticeSIM/phasePSF/isolated/DS/cropped/RAW_488_slow_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0013199486msecAbs_000x_000y_000z_0000t.tif'));

O = sim_PSFtoOTF_gen(PSF,'lattice_period', 1.4256, 'phase_step', .2851, 'norders', 5, 'nphases', 5, 'useGPU', true);

data = '/clusterfs/fiona/Data/20210923_latticeSIM/data06_100perc/DS/RAW_exp08_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0006674367msecAbs_000x_000y_000z_0000t_cropped.tif';
output = simReconFrame(data, O, 'lattice_period', 1.4256, 'phase_step', .2851, 'norders', 5, 'nphases', 5, 'Overlap', 64, 'ChunkSize', [256,256,256], 'edgeTaper', true, 'perdecomp', true);
%dataFull = loadtiff(data);
%output = simRecon(dataFull, O, 'lattice_period', 1.4, 'phase_step', .4667, 'norders', 3, 'nphases', 3, 'edgeTaper', true, 'perdecomp', true, 'useGPU', false);

writetiff(single(output.*(output>=0)),'/clusterfs/fiona/Data/20210923_latticeSIM/data06_100perc/DS/sim_recon/RAW_exp08_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0006674367msecAbs_000x_000y_000z_0000t_cropped.tif');

%%
PSF=double(loadtiff('/clusterfs/fiona/Data/20211005_latticeSIM/PSFs/488_NA0p4_sig0p1_highSN/DS/RAW_exp01_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0000558259msecAbs_000x_000y_000z_0000t.tif'));

O = sim_PSFtoOTF_gen(PSF,'lattice_period', 1.4, 'phase_step', 1.4./5, 'norders', 5, 'nphases', 5, 'useGPU', true);

data = '/clusterfs/fiona/Data/20211005_latticeSIM/wait4ms_2p6msInt_10perc/DS/RAW_exp01_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0005858093msecAbs_000x_000y_000z_0000t.tif';
output = simReconFrame(data, O, 'lattice_period', 1.4, 'phase_step', 1.4./5, 'norders', 5, 'nphases', 5, 'Overlap', 32, 'ChunkSize', [32,32,32], 'edgeTaper', true, 'perdecomp', true, 'DS', true);
%dataFull = loadtiff(data);
%output = simRecon(dataFull, O, 'lattice_period', 1.4, 'phase_step', .4667, 'norders', 3, 'nphases', 3, 'edgeTaper', true, 'perdecomp', true, 'useGPU', false);

%writetiff(single(output.*(output>=0)),'/clusterfs/fiona/Data/20210923_latticeSIM/data06_100perc/DS/sim_recon/RAW_exp08_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0006674367msecAbs_000x_000y_000z_0000t_cropped.tif');
