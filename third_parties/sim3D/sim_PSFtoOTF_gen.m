function [O] = sim_PSFtoOTF_gen(PSF, varargin)
%Generates and saves a centered SIM OTF for each orientation

tStart = tic;


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('PSF');
ip.addParameter('nphases', 5, @isnumeric);
ip.addParameter('norders', 5, @isnumeric);
ip.addParameter('norientations', 1, @isnumeric);
ip.addParameter('lattice_period', 1.2021, @isnumeric); %Lattice period in microns - this is the coarsest period
ip.addParameter('phase_step', .232, @isnumeric); %Phase step in microns
ip.addParameter('pxl_dim_PSF', [0.11,0.11,0.2*sind(32.4)], @isnumeric); %Voxel dimensions of the PSF in microns - note, stored as [dy,dx,dz]
ip.addParameter('Background', 105, @isnumeric); 

ip.addParameter('displayFit', false, @islogical);

ip.addParameter('Save16bit', false , @islogical);
ip.addParameter('gpuPrecision', 'double', @ischar);

ip.addParameter('useGPU', true, @islogical);

ip.addParameter('saveOTF', false, @islogical);
ip.addParameter('Overwrite', false, @islogical);

ip.parse(PSF, varargin{:});

pr = ip.Results;
nphases = pr.nphases;
norders = pr.norders;
norientations = pr.norientations;
lattice_period = pr.lattice_period;
phase_step = pr.phase_step;
pxl_dim_PSF = pr.pxl_dim_PSF;
Background = pr.Background;

displayFit = pr.displayFit;

Save16bit = pr.Save16bit;
gpuPrecision = pr.gpuPrecision;

useGPU = pr.useGPU;
useGPU = useGPU & gpuDeviceCount > 0;

saveOTF = pr.saveOTF;
Overwrite = pr.Overwrite;

%Load the PSF data

if(ischar(PSF))
    PSF = readtiff(PSF);
end

if useGPU
    PSF = gpuArray(PSF);
end

PSF=PSF-Background;
PSF(PSF<0)=0;
[ny_PSF,nx_PSF,nimgs_PSF] = size(PSF);
nz_PSF=nimgs_PSF/(nphases*norientations);

if useGPU
    dk_PSF=cast(gpuArray(1./([ny_PSF,nx_PSF,nz_PSF].*pxl_dim_PSF)),gpuPrecision);
    wfPSF=gpuArray(zeros(ny_PSF,nx_PSF,nz_PSF,gpuPrecision));
else
    dk_PSF=1./([ny_PSF,nx_PSF,nz_PSF].*pxl_dim_PSF);
    wfPSF=zeros(ny_PSF,nx_PSF,nz_PSF);
end
%First - center the bead in the volume to get rid of any phase ramps in the
%OTF. We do this by fitting the widefield PSF which is the sum of all
%phase stepped images per plane
tt=1;
for ii=1:nz_PSF
    wfPSF(:,:,ii)=sum(PSF(:,:,tt:tt+(nphases*norientations)-1),3);
    tt=tt+nphases*norientations;
end

%Compute the centroid of the wfPSF by fitting a 3D Gaussian
if useGPU
    wfPSF = gather(wfPSF);
end
fitParams=Lsq_GaussianFit_3D(double(wfPSF),displayFit);

%Separate the images for each phase and generate the OTFs

if useGPU
    O=gpuArray(zeros(ny_PSF,nx_PSF,nz_PSF,norders,norientations,gpuPrecision));
else
    O=zeros(ny_PSF,nx_PSF,nz_PSF,norders,norientations);
end

%['Separating information components']
%Make the forward separation matrix
[sep_matrix]=make_forward_separation_matrix(nphases,norders,lattice_period,phase_step);

%Make the inverse separation matrix
inv_sep_matrix=cast(pinv(sep_matrix),gpuPrecision);
    
if useGPU
    inv_sep_matrix = gpuArray(inv_sep_matrix);
end

for jj=1:norientations
    
    if useGPU
        Dr=gpuArray(zeros(ny_PSF,nx_PSF,nz_PSF,nphases,gpuPrecision));
        Dr_shift=gpuArray(zeros(ny_PSF,nx_PSF,nz_PSF,nphases,gpuPrecision));
        Dk=gpuArray(zeros(ny_PSF,nx_PSF,nz_PSF,nphases,gpuPrecision));
    else
        Dr=zeros(ny_PSF,nx_PSF,nz_PSF,nphases);
        Dr_shift=zeros(ny_PSF,nx_PSF,nz_PSF,nphases);
        Dk=zeros(ny_PSF,nx_PSF,nz_PSF,nphases);
    end
    
    for ii=1:nphases
        Dr(:,:,:,ii)=PSF(:,:,ii+(jj-1)*nphases:nphases*norientations:end);
        %Note - the ceil((nx+1)/2) is there to account for how fftshift and
        %ifftshift deal with even vs. odd sized datasets. Centering the bead at
        %ceil((nx+1)/2)... gives zero phase ramp in the OTF.
        Dr_shift(:,:,:,ii)=fourierShift3D(Dr(:,:,:,ii),[ceil((ny_PSF+1)/2)-fitParams(2),ceil((nx_PSF+1)/2)-fitParams(3),ceil((nz_PSF+1)/2)-fitParams(4)],useGPU,gpuPrecision);
        Dk(:,:,:,ii)=fftshift(ifftn(ifftshift(Dr_shift(:,:,:,ii)))).*1/prod(dk_PSF);
    end
    
    %Separate OTF orders by solving the linear system of equations for each pixel
    % for ii=1:nphases
    %    for kk=1:nphases
    %         O(:,:,:,ii,jj)=O(:,:,:,ii,jj)+inv_sep_matrix(ii,kk)*Dk(:,:,:,kk);
    %     end
    % end
    O(:, :, :, :, jj) = reshape(reshape(Dk, [], nphases) * inv_sep_matrix.', size(O, 1 : 4));

end

%Normalize OTF's so that 0th order of each orientation has unit energy
for jj=1:norientations
    OTF0=O(:,:,:,ceil(norders/2),jj);
    energy=sum(OTF0(:).*conj(OTF0(:))*prod(dk_PSF));
    O(:,:,:,:,jj)=O(:,:,:,:,jj)/sqrt(energy);
end

if useGPU
    O = gather(O);
end

if saveOTF
    [path, fn, ext] = fileparts(PSF);
    if ~exist([path '/' 'otf' '/' fn '.mat'],'file')
        mkdir([path '/' 'otf' '/']);
        save([path '/' 'otf' '/' fn '.mat'],'O');
    end
end
%save([PSF_folder,PSF_file(1:end-5),'_OTF_normalized.mat'],'O')
toc(tStart)
end
