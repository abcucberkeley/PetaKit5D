function [Data_r_space] = simRecon_Frame(data, otf, varargin)
% This code will generate a assemble a super-resolved SIM image from a deskewed 
% dataset and OTF acquired with 5-phase lattice SIM data.

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data');
ip.addRequired('otf');
ip.addParameter('islattice', true, @islogical); %Flag to indicate if this is light sheet data
ip.addParameter('NA_det', 1.0, @isnumeric);
ip.addParameter('NA_ext', 0.55, @isnumeric);
ip.addParameter('nimm', 1.33, @isnumeric);
ip.addParameter('wvl_em', .605, @isnumeric);
ip.addParameter('wvl_ext', .560, @isnumeric);
ip.addParameter('w', 5e-3, @isnumeric); %Wiener coefficient for regularization
ip.addParameter('apodize', true, @islogical); %Flag to indicate whether or not to apodize the final data

ip.addParameter('nphases', 5, @isnumeric);
ip.addParameter('norders', 5, @isnumeric);
ip.addParameter('norientations', 1, @isnumeric);
ip.addParameter('lattice_period', 1.2021, @isnumeric); %Lattice period in microns - this is the coarsest period
ip.addParameter('lattice_angle', [pi/2], @isnumeric); %Angle parellel to pattern modulation (assuming horizontal is zero)
ip.addParameter('phase_step', .232, @isnumeric); %Phase step in microns
ip.addParameter('pxl_dim_data', [0.11,0.11,0.3*sind(32.4)], @isnumeric); %Voxel dimensions of the image in microns - note, stored as [dy,dx,dz]
ip.addParameter('pxl_dim_PSF', [0.11,0.11,0.2*sind(32.4)], @isnumeric); %Voxel dimensions of the PSF in microns - note, stored as [dy,dx,dz]
ip.addParameter('background', 105, @isnumeric); 

ip.parse(data, otf, varargin{:});

pr = ip.Results;
islattice = pr.islattice;
NA_det = pr.NA_det;
NA_ext = pr.NA_ext;
nimm = pr.nimm;
wvl_em = pr.wvl_em;
wvl_ext = pr.wvl_ext;
w = pr.w;
apodize = pr.apodize;

nphases = pr.nphases;
norders = pr.norders;
norientations = pr.norientations;
lattice_period = pr.lattice_period;
lattice_angle = pr.lattice_angle;
phase_step = pr.phase_step;
pxl_dim_data = pr.pxl_dim_data;
pxl_dim_PSF = pr.pxl_dim_PSF;
background = pr.background;



%Setup initial parameters:
%islattice=1; %Flag to indicate if this is light sheet data
%NA_det=1.0;
%NA_ext=0.55;
%nimm=1.33;
%wvl_em=.605;
%wvl_ext=.560;
%w = 5e-3; %Wiener coefficient for regularization
%apodize=1; %Flag to indicate whether or not to apodize the final data
%kvec_search_range=.5; %The +- range of k-space to search for information overlap (represented at microns in real space around estimated pattern period)

%{
nphases = 5;
norders = 5;
norientations = 1;
lattice_period = 1.2021; %Lattice period in microns - this is the coarsest period
lattice_angle = [pi/2]; %Angle parellel to pattern modulation (assuming horizontal is zero)
phase_step = .232; %Phase step in microns
pxl_dim_data = [0.11,0.11,0.3*sind(32.4)]; %Voxel dimensions of the image in microns - note, stored as [dy,dx,dz]
pxl_dim_PSF = [0.11,0.11,0.2*sind(32.4)];  %Voxel dimensions of the PSF in microns - note, stored as [dy,dx,dz]
background=105; %Background to subtract from the images - necessary to remove DC from fourier space
%}

%Folders and paths to the data and OTF
%data_folder = ['/clusterfs/fiona/matthewmueller/20210830SimRecon3D/2020_11_12(Code_for_Kitware)/2020_11_12(Code_for_Kitware)/data/'];
%data_file = ['Cell_5phase_small.tif'];
%OTF_folder = ['/clusterfs/fiona/matthewmueller/20210830SimRecon3D/2020_11_12(Code_for_Kitware)/2020_11_12(Code_for_Kitware)/data/'];
%OTF_file = ['PSF_5phas_OTF_normalized.mat'];

%Load the OTF
%disp("SIM_processdata_withShiftRefinement_5_phase_lattice_small.mlx")
disp("Load OTF")
tic
%O=load([OTF_folder,'/',OTF_file]);
fns = fieldnames(otf);
otf=otf.(fns{1});
[ny_PSF,nx_PSF,nz_PSF,~,~] = size(otf);
dk_PSF=gpuArray(1./([ny_PSF,nx_PSF,nz_PSF].*pxl_dim_PSF));
toc

%Load the data
disp("Load data")
tic
%data=loadtiff([data_folder,data_file]);
data=data-background;
data(data<0)=0;

[ny_data,nx_data,nimgs_data] = size(data);
nz_data=nimgs_data/(nphases*norientations);
dk_data=gpuArray(1./([ny_data,nx_data,nz_data].*pxl_dim_data));
toc

%Rescale OTF to be the same size as the image data - this is done by
%interpolating the complex OTF onto a grid that's the same size as the
%k-space representation of the data
disp("Rescale OTF")
tic
[O_scaled]=resampleOTF(otf,pxl_dim_PSF,data,pxl_dim_data,nphases,norders,norientations);
toc

%Normalize resampled OTF's so that 0th order of each orientation has unit energy
disp("Normalize OTF")
%dk_PSF = gpuArray(dk_PSF);
tic
%O_scaled = gather(O_scaled);
for jj=1:norientations
    OTF0=O_scaled(:,:,:,ceil(norders/2),jj);
    energy=sum(OTF0(:).*conj(OTF0(:))*prod(dk_PSF));
    O_scaled(:,:,:,:,jj)=O_scaled(:,:,:,:,jj)/sqrt(energy);
end
clear OTF0
clear O
toc

%Make a cylindrical mask based on the theoretical lateral and axial OTF extent
disp("cylMask")
tic
cylmask=cylMask(NA_det,wvl_em,NA_ext,wvl_ext,nimm,norders,ny_data,nx_data,nz_data,dk_data,islattice);
toc

% Notes: right now, this masking is just used when determining the starting 
% phase and modulation depth. There are likely other ways to mask the OTFs (e.g. 
% via SNR).

%Separate data information orders by solving the linear system of equations for each pixel
disp("Separate data information orders by solving the linear system of equations for each pixel (perdecomp_3D/make_forward_seperation_matrix)")
tic
Dk_sep=gpuArray(zeros(ny_data,nx_data,nz_data,norders,norientations));
for jj=1:norientations
    
    %Separate the images for each phase and generate the separated Dk orders
    Dr = gpuArray(zeros(ny_data,nx_data,nz_data,nphases));
    Dk = gpuArray(zeros(ny_data,nx_data,nz_data,nphases));
    
    for ii=1:nphases
        Dr(:,:,:,ii)=double(data(:,:,ii+(jj-1)*nphases:nphases*norientations:end));
        Dk(:,:,:,ii)=double(fftshift(ifftn(ifftshift(Dr(:,:,:,ii))))*1/prod(dk_data));
    end
    
    %['Separating information components']
    %Make the forward separation matrix
    [sep_matrix]=make_forward_separation_matrix(nphases,norders,lattice_period,phase_step);
    
    %Make the inverse separation matrix
    inv_sep_matrix=pinv(sep_matrix);
    
    for ii=1:nphases
        %ii
        for kk=1:nphases
            Dk_sep(:,:,:,ii,jj)=Dk_sep(:,:,:,ii,jj)+inv_sep_matrix(ii,kk)*Dk(:,:,:,kk);
        end
    end
end
clear Dk
clear Dr
clear data
toc

%Now we need to correct shift vector and starting lateral phase for each
%information component Dk_sep
disp("correct shift vector and starting lateral phase for each and scale shifted copies (fourierShift3D/cylmask/MaskedTranslationRegistration2D_fit)")
tic

for jj=1:norientations
    ai=[]; bi=[];
    for kk=1:floor(norders/2) %Only consider negative orders - positive orders a just the complex conjugate and opposite direction
        transform_tot=[0,0,0];
        order=(kk-ceil(norders/2)); %m'th order information component
        k_shift_scaling = 1/lattice_period./dk_data; %Convert pattern frequency into k-space pixel units
        
        %Estimated shift vector
        p_vec_guess(kk,:,jj)=order*k_shift_scaling.*[sin(lattice_angle(jj)),cos(lattice_angle(jj)),0];
        
        for qq=1:2
        %Shift image frequency information - D˜m(k+mp)
        shift_Dk_sep=fourierShift3D(Dk_sep(:,:,:,kk,jj),p_vec_guess(kk,:,jj));
        
        %Shift transfer function - O_m(k+mp)
        shift_Om=fourierShift3D(O_scaled(:,:,:,kk,jj),p_vec_guess(kk,:,jj));
        
        %Shift mask - we use imtranslate here because for small images,
        %fourier shifting can wrap around the image
        %p_vec_guess = gather(p_vec_guess);
        cylmask_shift=imtranslate_function(double(cylmask(:,:,:,kk)),[p_vec_guess(kk,2,jj),p_vec_guess(kk,1,jj),p_vec_guess(kk,3,jj)],'FillValues',0);
        cylmask_shift=abs(cylmask_shift)>.5; %Get rid of non-logical values due to interpolation
        
        %Overlap mask with the zero-information component
        overlap_mask= cylmask(:,:,:,3)&cylmask_shift;
        overlap_mask = gather(overlap_mask);
        

% Next, scale the shifted copies by either the shifted or non-shifted OTF as 
% per equations 9* and 9**. 

        DmO0=shift_Dk_sep.*O_scaled(:,:,:,ceil(norders/2),jj); %D˜m(k+mp)O_0(k) - eqn (9*)
        D0Om=Dk_sep(:,:,:,ceil(norders/2),jj).*shift_Om; %D˜0(k)O_m(k+mp) - eqn (9**)
        D0Om = gather(D0Om);
        DmO0 = gather(DmO0);
        %Right now, we only do 1 round of shift vector refinement. In
        %practice, we could iterate.
        if qq==1
                [transform,maxC,C,numberOfOverlapMaskedPixels] = MaskedTranslationRegistration2D_fit(abs(sum(D0Om,3)),abs(sum(DmO0,3)),max(overlap_mask,[],3),max(overlap_mask,[],3),.5);
                transform=[transform(2),transform(1),0];
                p_vec_guess(kk,:,jj)=p_vec_guess(kk,:,jj)-transform;
                transform_tot=transform_tot+transform;
                norm(transform);
                maxC;
        end
        end
        
        %['lattice period estimate from order ', num2str(order), ' = ', num2str(abs(order)/norm(p_vec_guess(kk,:).*dk_data)), ' microns'];
        %['lattice angle estimate from order ', num2str(order), ' = ', num2str(acosd(p_vec_guess(kk,2)/p_vec_guess(kk,1))),' degrees'];
        
        %B is a complex scale factor that incorporates both the difference in
        %modulation depth and starting phase from when the OTF and data are
        %acquired.
        ai=D0Om(overlap_mask(:));
        bi=DmO0(overlap_mask(:));
        
        B_TLS(kk,jj) = tls(ai,bi); %Total least squares regression on two complex vectors
        B_PLS(kk,jj) = sum(conj(ai).*bi)/sum(abs(ai).^2); %Partial least squares regression on the two complex vectors
    end
    
    %Positive order shifts are just the opposite of negative order shifts and
    %zero order is 0
    p_vec_guess(ceil(norders/2),:,jj)=0;
    B_TLS(ceil(norders/2),jj)=1+1i*0;
    B_PLS(ceil(norders/2),jj)=1+1i*0;
    
    p_vec_guess(ceil(norders/2)+1:norders,:,jj)=-p_vec_guess(floor(norders/2):-1:1,:,jj);
    B_TLS(ceil(norders/2)+1:norders,jj)=conj(B_TLS(floor(norders/2):-1:1,jj));
    B_PLS(ceil(norders/2)+1:norders,jj)=conj(B_PLS(floor(norders/2):-1:1,jj));
end
clear shift_O
clear shift_Om
clear shift_Dk_sep
clear D0Om
clear DmO0
clear cylmask_shift
toc


disp("Generate Weiner Filter (fourierShift3D)")
tic
%'Generating Weiner Filter'
supersample=[2,2,1];
padrange=floor((supersample.*[ny_data,nx_data,nz_data]-[ny_data,nx_data,nz_data])/2);

B_PLS = gpuArray(B_PLS);
Dk_sep_scaled=gpuArray(padarray(zeros(ny_data,nx_data,nz_data),padrange,'both'));

%Assemble denominator of Weiner filter - this will then be shifted opposite to the p-vector direction to
%normalize each of the separated information components prior to shifting
%to their true locations in k-space

%denom=padarray(zeros(ny_data,nx_data,nz_data),padrange,'both');
denom = gpuArray(padarray(zeros(ny_data,nx_data,nz_data),padrange,'both'));
%p_vec_guess = gpuArray(p_vec_guess);

ii=ceil(norders/2);
for qq=1:norientations
    for kk=1:norders
        big_Ospace=padarray(O_scaled(:,:,:,kk,qq).*B_PLS(kk,qq),padrange,'both');
        p_vec_shift=p_vec_guess(kk,:,qq);
        shift_O=fourierShift3D(big_Ospace,p_vec_shift);
        denom=denom+abs(shift_O).^2;
    end
end

clear big_Ospace
clear shift_O
toc

%TESTING
%disp("Filtering information components")
%tic
%Dk_sep_scaled = filterInformationComponents(ny_data,nx_data,nz_data,norientations,norders,denom,p_vec_guess,Dk_sep,O_scaled,B_PLS,padrange,w);
%toc

%'filtering information components'
disp("Filtering information components")
tic
for jj=1:norientations
    for ii=1:norders
        shift_denom=imtranslate_function(denom,-[p_vec_guess(ii,2,jj),p_vec_guess(ii,1,jj),p_vec_guess(ii,3,jj)]);
        shift_denom = gpuArray(shift_denom);
        numerator=padarray((Dk_sep(:,:,:,ii,jj).*conj(O_scaled(:,:,:,ii,jj).*B_PLS(ii,jj))),padrange,'both');
        Dk_sep_scaled(:,:,:,ii,jj)=numerator./(shift_denom+w^2);
    end
end
clear O_scaled
clear Dk_sep
clear numerator
clear denom
clear shift_denom
toc


% The next step is to shift these scaled information components back into their 
% correct locations in frequency space and then add them together into the final 
% dataset. This is pretty straightforward since we already know the correct shift 
% vectors. To account for this shift, we need a larger frequency-space (which 
% is basically super-sampling with smaller pixels in real space) compared to our 
% original image. For now, we'll use a super-sampling factor of 2.

%'assembling final dataset'
disp("assembling final dataset")
tic
Data_k_space = gpuArray(padarray(zeros(ny_data,nx_data,nz_data),padrange,'both'));
p_vec_guess = gpuArray(p_vec_guess);
%Shift information components and assemble final image
for jj=1:norientations
%Assemble final dataset
for ii=1:norders
    big_kspace_shifted=fourierShift3D(Dk_sep_scaled(:,:,:,ii,jj),p_vec_guess(ii,:,jj));
    Data_k_space=Data_k_space+big_kspace_shifted;
    clear big_kspace_shifted
end
end

clear Dk_sep_scaled
toc

%Apodize final frequency space data using a triangular apodization function
disp("apodizeEllipse")
tic
if apodize
[Data_k_space_apodized,apo_ellipse] = apodizeEllipse(Data_k_space,dk_data,p_vec_guess,lattice_angle,NA_det,nimm,NA_ext,wvl_em,wvl_ext,islattice);
end
toc

%Transform final dataset back into real space
disp("transform data back to real space")
tic
Data_r_space=real(fftshift(fftn(ifftshift(Data_k_space_apodized))))*prod(dk_data);
toc

%disp("write to disk (write3Dtiff)")
tic
Data_r_space = gather(Data_r_space);
%Write the data to disk
write3Dtiff(single(abs(Data_r_space)),['/clusterfs/fiona/matthewmueller/20210830SimRecon3D/2020_11_12(Code_for_Kitware)/2020_11_12(Code_for_Kitware)/data/siRecon_5phase_small.tif'])
toc

end