function [mask] = cylMask(NA_det,wvl_em,NA_ext,wvl_ext,nimm,norders,ny_data,nx_data,nz_data,dk_data,islattice,useGPU)
%Make a cylindrical mask based on the theoretical lateral and axial OTF
%extents. Works for both 3D SIM and lattice SIM
alpha_det=asin(NA_det/nimm); %Half angle for detection NA
kmag_det=nimm/wvl_em; %k-vector magnitude for detection
alpha_ext=asin(NA_ext/nimm); %Half angle for detection NA
kmag_ext=nimm/wvl_ext; %k-vector magnitude for detection
radial_extent=(2*NA_det/wvl_em);

if(useGPU)
    mask = gpuArray(zeros(ny_data,nx_data,nz_data,norders,'logical'));
else
    mask = zeros(ny_data,nx_data,nz_data,norders,'logical');
end

%Set up vectors based on data image size
kxxdata=[-ceil((nx_data-1)/2):floor((nx_data-1)/2)]*dk_data(2);
kyydata=[-ceil((ny_data-1)/2):floor((ny_data-1)/2)]*dk_data(1);
kzzdata=[-ceil((nz_data-1)/2):floor((nz_data-1)/2)]*dk_data(3);

if(useGPU)
    kxxdata = gpuArray(kxxdata);
    kyydata = gpuArray(kyydata);
    kzzdata = gpuArray(kzzdata);
end

%Set up arrays for data frequency space pixels
[kxxdata_arr,kyydata_arr,kzzdata_arr]=meshgrid(kxxdata,kyydata,kzzdata);

if(useGPU)
    axial_extent = gpuArray(zeros(norders,1));
else
    axial_extent=zeros(norders,1);
end
for ii=1:norders
    if islattice
        axial_extent(ii)=kmag_det*(1-cos(alpha_det))+(floor(norders/2)-abs((ceil(norders/2)-ii)*.5))*kmag_ext*(sin(alpha_ext));
    else
        axial_extent(ii)=kmag_det*(1-cos(alpha_det))+(1-mod(ii,2))*kmag_ext*(1-cos(alpha_ext));
    end
    
    mask(:,:,:,ii)=(sqrt(kxxdata_arr.^2+kyydata_arr.^2)<=radial_extent)&(sqrt(kzzdata_arr.^2)<=axial_extent(ii));
end
end

