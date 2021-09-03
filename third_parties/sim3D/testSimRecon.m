tStart = tic;

data = loadtiff('/clusterfs/fiona/matthewmueller/20210830SimRecon3D/2020_11_12(Code_for_Kitware)/2020_11_12(Code_for_Kitware)/data/Cell_5phase_small.tif');
otf = load('/clusterfs/fiona/matthewmueller/20210830SimRecon3D/2020_11_12(Code_for_Kitware)/2020_11_12(Code_for_Kitware)/data/PSF_5phas_OTF_normalized.mat');


test = simRecon_Frame(data,otf);
toc(tStart)

