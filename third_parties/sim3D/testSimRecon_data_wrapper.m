

dataF = {'/clusterfs/fiona/Data/20211006_Aang_latticeSIM_conSweep_interleave/continuesSweep_150Stacks_FOV02/DS/',...
    '/clusterfs/fiona/Data/20211006_Aang_latticeSIM_conSweep_interleave/continuesSweep_150Stacks/DS/', ...
    '/clusterfs/fiona/Data/20211005_latticeSIM/wait4ms_2p6msInt_10perc_150stacks_1800x512/DS/',...
    '/clusterfs/fiona/Data/20211005_latticeSIM/wait4ms_2p6msInt_10perc_150stacks_1800x512_FOV2/DS/'};
rt_psf = '/clusterfs/fiona/Data/20211006_Aang_latticeSIM_conSweep_interleave/PSFs/';
fn = '488_NA0p4_sig0p1_highSN/DS/RAW_exp01_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0000558259msecAbs_000x_000y_000z_0000t.tif';
fn2 = '560_NA0p46_sig0p1_highSN/DS/RAW_exp01_CamB_ch0_CAM1_stack0000_488nm_0000000msec_0000751255msecAbs_000x_000y_000z_0000t.tif';

lattice_period = 1.4;
norders = 5;
nphases = 5;
phase_step = lattice_period./nphases; % in um
useGPU = true;

ChunkSize = [32, 32, 32];
Overlap = 32;
edgeTaperVal = 0.1;
ds = true;
Background488 = 105;
Background560 = 105;

simRecon_data_wrapper(dataF,{[rt_psf fn], [rt_psf fn2]},'lattice_period', lattice_period,...
                        'phase_step', phase_step, 'norders', norders, 'nphases', nphases,...
                        'Overlap', Overlap, 'ChunkSize', ChunkSize, 'edgeTaper', true, 'edgeTaperVal', edgeTaperVal,...
                        'perdecomp', true, 'useGPU', useGPU, 'DS', ds,'Background', Background488);