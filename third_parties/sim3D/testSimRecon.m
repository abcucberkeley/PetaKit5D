tStart = tic;

data = '/clusterfs/fiona/matthewmueller/20210830SimRecon3D/2020_11_12(Code_for_Kitware)/2020_11_12(Code_for_Kitware)/data/Cell_5phase_small.tif';
otf = '/clusterfs/fiona/matthewmueller/20210830SimRecon3D/2020_11_12(Code_for_Kitware)/2020_11_12(Code_for_Kitware)/data/PSF_5phas_OTF_normalized.mat';
dataL = '/clusterfs/fiona/matthewmueller/20210830SimRecon3D/2020_11_12(Code_for_Kitware)/2020_11_12(Code_for_Kitware)/data/Cell_5phase_large.tif';

otf = load(otf);
fns = fieldnames(otf);
otf=otf.(fns{1});

%dataSR = readtiff(data);
%dataSRL = readtiff(dataL);
%otfSR = load(otf);
%fns = fieldnames(otfSR);
%otfSR=otfSR.(fns{1});
simReconFrame(dataL,otf,'useGPU',true, 'perdecomp', true, 'edgeTaper', true, 'DS', true);
%simRecon(dataSRL,otfSR,'useGPU',false);
toc(tStart)

