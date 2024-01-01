function [Data_r_space] = simRecon(data, otf, varargin)
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
ip.addParameter('wVer', 1, @isnumeric); %Wiener coefficient version
ip.addParameter('apodize', true, @islogical); %Flag to indicate whether or not to apodize the final data

ip.addParameter('normalize_orientations', false, @islogical); %Flag to indicate whether or not to normalize total intensity for each orientation
ip.addParameter('perdecomp', true, @islogical); %Flag to indicate whether or not to use periodic/smooth decomposition to reduce edge effects
ip.addParameter('edgeTaper', true, @islogical); %Flag to indicate whether or not to window the data to reduce edge effects
ip.addParameter('edgeTaperVal', 0.1, @isnumeric); %Roll-off parameter for Tukey windowing

ip.addParameter('nphases', 5, @isnumeric);
ip.addParameter('norders', 5, @isnumeric);
ip.addParameter('norientations', 1, @isnumeric);
ip.addParameter('lattice_period', 1.2021, @isnumeric); %Lattice period in microns - this is the coarsest period
ip.addParameter('lattice_angle', [pi/2], @isnumeric); %Angle parellel to pattern modulation (assuming horizontal is zero)
%ip.addParameter('starting_angle', 35, @isnumeric);
ip.addParameter('phase_step', .232, @isnumeric); %Phase step in microns
ip.addParameter('pxl_dim_data', [0.11,0.11,0.3*sind(32.4)], @isnumeric); %Voxel dimensions of the image in microns - note, stored as [dy,dx,dz]
ip.addParameter('pxl_dim_PSF', [0.11,0.11,0.2*sind(32.4)], @isnumeric); %Voxel dimensions of the PSF in microns - note, stored as [dy,dx,dz]
ip.addParameter('Background', 105, @isnumeric);

ip.addParameter('Save16bit', false , @islogical);
ip.addParameter('gpuPrecision', 'single', @ischar);

ip.addParameter('useGPU', true, @islogical);

ip.parse(data, otf, varargin{:});

pr = ip.Results;
islattice = pr.islattice;
NA_det = pr.NA_det;
NA_ext = pr.NA_ext;
nimm = pr.nimm;
wvl_em = pr.wvl_em;
wvl_ext = pr.wvl_ext;
w = pr.w;
wVer = pr.wVer;
apodize = pr.apodize;

normalize_orientations = pr.normalize_orientations;
perdecomp = pr.perdecomp;
edgeTaper = pr.edgeTaper;
edgeTaperVal = pr.edgeTaperVal;

nphases = pr.nphases;
norders = pr.norders;
norientations = pr.norientations;
lattice_period = pr.lattice_period;
lattice_angle = pr.lattice_angle;
%starting_angle = pr.starting_angle;

phase_step = pr.phase_step;
pxl_dim_data = pr.pxl_dim_data;
pxl_dim_PSF = pr.pxl_dim_PSF;
Background = pr.Background;

Save16bit = pr.Save16bit;
gpuPrecision = pr.gpuPrecision;

useGPU = pr.useGPU;
useGPU = useGPU & gpuDeviceCount > 0;

%kvec_search_range=.5; %The +- range of k-space to search for information overlap (represented at microns in real space around estimated pattern period)

%Load the OTF
%disp("Load OTF")
%tic
[ny_PSF,nx_PSF,nz_PSF,~,~] = size(otf);

if(useGPU)
    dk_PSF=cast(gpuArray(1./([ny_PSF,nx_PSF,nz_PSF].*pxl_dim_PSF)),gpuPrecision);
else
    dk_PSF=1./([ny_PSF,nx_PSF,nz_PSF].*pxl_dim_PSF);
end

%toc

%Load the data
%disp("Load data")
%tic
data=data-Background;
data(data<0)=0;

[ny_data,nx_data,nimgs_data] = size(data);
nz_data=nimgs_data/(nphases*norientations);

if(useGPU)
    dk_data=cast(gpuArray(1./([ny_data,nx_data,nz_data].*pxl_dim_data)),gpuPrecision);
else
    dk_data=1./([ny_data,nx_data,nz_data].*pxl_dim_data);
end

if normalize_orientations
    index=zeros(nz_data*nphases,norientations);
    for jj=1:norientations
        tt=1;
        for kk=1:nz_data
            index(tt:tt+nphases-1,jj)=(((jj-1)*nphases)+(kk-1)*nphases*norientations)+1:((jj-1)*nphases)+((kk-1)*nphases*norientations)+nphases;
            tt=tt+nphases;
        end
        wf_orientation_intensity(jj)=sum(data(:,:,index(:,jj)),"all");
    end
    for jj=1:norientations
        data(:,:,index(:,jj))=data(:,:,index(:,jj))./(wf_orientation_intensity(jj)/max(wf_orientation_intensity));
    end
end
%toc

%Rescale OTF to be the same size as the image data - this is done by
%interpolating the complex OTF onto a grid that's the same size as the
%k-space representation of the data
%disp("Rescale OTF")
%tic
[O_scaled]=resampleOTF(otf,pxl_dim_PSF,data,pxl_dim_data,nphases,norders,norientations,useGPU,gpuPrecision);
%toc

%Normalize resampled OTF's so that 0th order of each orientation has unit energy
%disp("Normalize OTF")
%tic
for jj=1:norientations
    OTF0=O_scaled(:,:,:,ceil(norders/2),jj);
    energy=sum(OTF0(:).*conj(OTF0(:))*prod(dk_PSF));
    O_scaled(:,:,:,:,jj)=O_scaled(:,:,:,:,jj)/sqrt(energy);
end
clear OTF0
clear O
%toc

%Make a cylindrical mask based on the theoretical lateral and axial OTF extent
%disp("cylMask")
%tic
cylmask=cylMask(NA_det,wvl_em,NA_ext,wvl_ext,nimm,norders,ny_data,nx_data,nz_data,dk_data,islattice,useGPU);
%toc
sz = size(cylmask);
sz = sz(1:3);


% Notes: right now, this masking is just used when determining the starting
% phase and modulation depth. There are likely other ways to mask the OTFs (e.g.
% via SNR).

%Separate data information orders by solving the linear system of equations for each pixel
%disp("Separate data information orders by solving the linear system of equations for each pixel (perdecomp_3D/make_forward_seperation_matrix)")
%tic
if(useGPU)
    Dk_sep=gpuArray(zeros(ny_data,nx_data,nz_data,norders,norientations,gpuPrecision));
else
    Dk_sep=zeros(ny_data,nx_data,nz_data,norders,norientations);
end

%['Separating information components']
%Make the forward separation matrix
[sep_matrix]=make_forward_separation_matrix(nphases,norders,lattice_period,phase_step);

%Make the inverse separation matrix
inv_sep_matrix=cast(pinv(sep_matrix),gpuPrecision);
if(useGPU)
    inv_sep_matrix = gpuArray(inv_sep_matrix);
end

for jj=1:norientations

    %Separate the images for each phase and generate the separated Dk orders
    if(useGPU)
        Dr = gpuArray(zeros(ny_data,nx_data,nz_data,nphases,gpuPrecision));
        Dk = gpuArray(zeros(ny_data,nx_data,nz_data,nphases,gpuPrecision));
    else
        Dr = zeros(ny_data,nx_data,nz_data,nphases);
        Dk = zeros(ny_data,nx_data,nz_data,nphases);
    end

    if edgeTaper
        [window] = tukwin(Dr(:,:,:,1),edgeTaperVal,useGPU);
    end

    for ii=1:nphases
        Dr(:,:,:,ii)=data(:,:,ii+(jj-1)*nphases:nphases*norientations:end);
        if edgeTaper
            Dr(:,:,:,ii) = Dr(:,:,:,ii).*window;
        end
        if perdecomp
            Dr(:,:,:,ii) = perdecomp_3D(Dr(:,:,:,ii),useGPU,gpuPrecision);
        end
        Dk(:,:,:,ii)=fftshift(ifftn(ifftshift(Dr(:,:,:,ii))))*1/prod(dk_data);
    end

   % for ii=1:nphases
    %    for kk=1:norders
     %       Dk_sep(:,:,:,kk,jj)=Dk_sep(:,:,:,kk,jj)+inv_sep_matrix(kk,ii)*Dk(:,:,:,ii);
      %  end
    %end
    Dk_sep(:, :, :, :, jj) = reshape(reshape(Dk, [], nphases) * inv_sep_matrix.', size(Dk_sep, 1 : 4));
end
clear Dk
clear Dr
clear data
%toc

%Now we need to correct shift vector and starting lateral phase for each
%information component Dk_sep
%disp("correct shift vector and starting lateral phase for each and scale shifted copies (fourierShift3D/cylmask/MaskedTranslationRegistration2D_fit)")
%tic

%starting_angle=35;
p_vec_guess = zeros(floor(norders/2), 3, norientations);
B_TLS = zeros(floor(norders/2), norientations);
B_PLS = zeros(floor(norders/2), norientations);

for jj=1:norientations
    %lattice_angle = [-pi/2+(starting_angle+(jj-1)*180/norientations)*pi/180];
    for kk=1:floor(norders/2) %Only consider negative orders - positive orders a just the complex conjugate and opposite direction
        transform_tot=[0,0,0];
        order=(kk-ceil(norders/2)); %m'th order information component
        k_shift_scaling = 1/lattice_period./dk_data; %Convert pattern frequency into k-space pixel units

        %Estimated shift vector
        % REPLACE WITH
        % p_vec_guess(kk,:,jj)=order*k_shift_scaling.*[sin(lattice_angle),cos(lattice_angle),0];
        % when adding starting angle
        p_vec_guess(kk,:,jj)=order*k_shift_scaling.*[sin(lattice_angle(jj)),cos(lattice_angle(jj)),0];

        for qq=1:2
            %Shift image frequency information - D˜m(k+mp)
            shift_Dk_sep=fourierShift3D(Dk_sep(:,:,:,kk,jj),p_vec_guess(kk,:,jj),useGPU,gpuPrecision);

            %Shift transfer function - O_m(k+mp)
            shift_Om=fourierShift3D(O_scaled(:,:,:,kk,jj),p_vec_guess(kk,:,jj),useGPU,gpuPrecision);

            %Shift mask - we use imtranslate here because for small images,
            %fourier shifting can wrap around the image
            % cylmask_shift=imtranslate_function(double(cylmask(:,:,:,kk)),[p_vec_guess(kk,2,jj),p_vec_guess(kk,1,jj),p_vec_guess(kk,3,jj)],'FillValues',0);
            % cylmask_shift_0=abs(cylmask_shift)>.5; %Get rid of non-logical values due to interpolation

            %sz = size(cylmask, [1, 2, 3]);
            cylmask_shift = false(sz);
            shift = round([p_vec_guess(kk,1,jj),p_vec_guess(kk,2,jj),p_vec_guess(kk,3,jj)]);
            so = max(1, 1 - shift);
            to = min(sz, sz - shift);
            sn = max(1, shift + 1);
            tn = min(sz, sz + shift);
            cylmask_shift(sn(1) : tn(1), sn(2) : tn(2), sn(3) : tn(3)) = cylmask(so(1) : to(1), so(2) : to(2), so(3) : to(3), kk);

            %Overlap mask with the zero-information component
            overlap_mask= cylmask(:,:,:,3)&cylmask_shift;


            % Next, scale the shifted copies by either the shifted or non-shifted OTF as
            % per equations 9* and 9**.

            DmO0=shift_Dk_sep.*O_scaled(:,:,:,ceil(norders/2),jj); %D˜m(k+mp)O_0(k) - eqn (9*)
            D0Om=Dk_sep(:,:,:,ceil(norders/2),jj).*shift_Om; %D˜0(k)O_m(k+mp) - eqn (9**)

            %Right now, we only do 1 round of shift vector refinement. In
            %practice, we could iterate.
            if qq==1
                [transform,~,~,~] = MaskedTranslationRegistration2D_fit(abs(sum(D0Om,3)),abs(sum(DmO0,3)),max(overlap_mask,[],3),max(overlap_mask,[],3),.5,useGPU,gpuPrecision);
                transform=[transform(2),transform(1),0];
                p_vec_guess(kk,:,jj)=p_vec_guess(kk,:,jj)-transform;
                transform_tot=transform_tot+transform;
                norm(transform);
                %maxC;
            end
        end

        %['lattice period estimate from order ', num2str(order), ' = ', num2str(abs(order)/norm(p_vec_guess(kk,:).*dk_data)), ' microns'];
        %['lattice angle estimate from order ', num2str(order), ' = ', num2str(acosd(p_vec_guess(kk,2)/p_vec_guess(kk,1))),' degrees'];

        %B is a complex scale factor that incorporates both the difference in
        %modulation depth and starting phase from when the OTF and data are
        %acquired.
        ai=D0Om(overlap_mask(:));
        bi=DmO0(overlap_mask(:));

        B_TLS(kk,jj) = tls(ai,bi,useGPU); %Total least squares regression on two complex vectors
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
clear cylmask
clear cylmask_shift
clear overlap_mask
clear ai
clear bi
%toc


%disp("Generate Weiner Filter (fourierShift3D)")
%tic
%'Generating Weiner Filter'
supersample=[2,2,1];
padrange=floor((supersample.*[ny_data,nx_data,nz_data]-[ny_data,nx_data,nz_data])/2);

if(useGPU)
    B_PLS = gpuArray(B_PLS);
    Dk_sep_scaled = gpuArray(zeros(ny_data + padrange(1) * 2, nx_data + padrange(2) * 2,nz_data + padrange(3) * 2,gpuPrecision));
    denom = gpuArray(zeros(ny_data + padrange(1) * 2, nx_data + padrange(2) * 2,nz_data + padrange(3) * 2,gpuPrecision));
else
    Dk_sep_scaled = zeros(ny_data + padrange(1) * 2, nx_data + padrange(2) * 2,nz_data + padrange(3) * 2);
    denom = zeros(ny_data + padrange(1) * 2, nx_data + padrange(2) * 2,nz_data + padrange(3) * 2);
end

%Assemble denominator of Weiner filter - this will then be shifted opposite to the p-vector direction to
%normalize each of the separated information components prior to shifting
%to their true locations in k-space
if wVer == 0 || wVer == 2
    for qq=1:norientations
        for kk=1:norders
            big_Ospace=padarray(O_scaled(:,:,:,kk,qq).*B_PLS(kk,qq),padrange,'both');
            p_vec_shift=p_vec_guess(kk,:,qq);
            shift_O=fourierShift3D(big_Ospace,p_vec_shift,useGPU,gpuPrecision);
            denom=denom+abs(shift_O).^2;
        end
    end
end

clear big_Ospace
clear shift_O
%toc

%'filtering information components'
%disp("Filtering information components")
%tic
for jj=1:norientations
    for ii=1:norders
        numerator=padarray((Dk_sep(:,:,:,ii,jj).*conj(O_scaled(:,:,:,ii,jj).*B_PLS(ii,jj))),padrange,'both');
        if wVer == 0
            shift_denom=imtranslate_function(denom,-[p_vec_guess(ii,2,jj),p_vec_guess(ii,1,jj),p_vec_guess(ii,3,jj)]);
            Dk_sep_scaled(:,:,:,ii,jj)=numerator./(shift_denom+w^2);
        elseif wVer == 1
             Dk_sep_scaled(:,:,:,ii,jj)=numerator./(abs(padarray(O_scaled(:,:,:,ii,jj).*B_PLS(ii,jj),padrange,'both')).^2+w^2);
        elseif wVer == 2
            Dk_sep_scaled(:,:,:,ii,jj)=numerator;
        end
        
    end
end
clear O_scaled
clear Dk_sep
clear numerator
clear denom
clear shift_denom
%toc


% The next step is to shift these scaled information components back into their
% correct locations in frequency space and then add them together into the final
% dataset. This is pretty straightforward since we already know the correct shift
% vectors. To account for this shift, we need a larger frequency-space (which
% is basically super-sampling with smaller pixels in real space) compared to our
% original image. For now, we'll use a super-sampling factor of 2.

%'assembling final dataset'
%disp("assembling final dataset")
%tic

if(useGPU)
    % Data_k_space = gpuArray(padarray(zeros(ny_data,nx_data,nz_data),padrange,'both'));
    Data_k_space = gpuArray(zeros(ny_data + padrange(1) * 2, nx_data + padrange(2) * 2,nz_data + padrange(3) * 2,gpuPrecision));
    p_vec_guess = gpuArray(p_vec_guess);
else
    % Data_k_space = padarray(zeros(ny_data,nx_data,nz_data),padrange,'both');
    Data_k_space = zeros(ny_data + padrange(1) * 2, nx_data + padrange(2) * 2,nz_data + padrange(3) * 2);
end
%Shift information components and assemble final image
for jj=1:norientations
    %Assemble final dataset
    for ii=1:norders
        Data_k_space=Data_k_space + fourierShift3D(Dk_sep_scaled(:,:,:,ii,jj),p_vec_guess(ii,:,jj),useGPU,gpuPrecision);
        %big_kspace_shifted=fourierShift3D(Dk_sep_scaled(:,:,:,ii,jj),p_vec_guess(ii,:,jj),useGPU,gpuPrecision);
        %Data_k_space=Data_k_space+big_kspace_shifted;
        %clear big_kspace_shifted
    end
end

if wVer == 2
    Data_k_space = Data_k_space./(denom+w^2);
end

clear Dk_sep_scaled
%toc

%Apodize final frequency space data using a triangular apodization function
%disp("apodizeEllipse")
%tic
if apodize
    [Data_k_space_apodized,~] = apodizeEllipse(Data_k_space,dk_data,p_vec_guess,lattice_angle,NA_det,nimm,NA_ext,wvl_em,wvl_ext,islattice,useGPU);
    %clear Data_k_space
    %Transform final dataset back into real space
    Data_r_space=real(fftshift(fftn(ifftshift(Data_k_space_apodized))))*prod(dk_data);
else
    Data_r_space=real(fftshift(fftn(ifftshift(Data_k_space))))*prod(dk_data);
end
%toc

%TESTING
%if(useGPU)
%    Data_r_space = gather(Data_r_space);
%end

end