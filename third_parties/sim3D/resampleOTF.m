function [OTF_scaled]=resampleOTF(inputOTF,pxl_dim_PSF,inputData,pxl_dim_data,nphases,norders,norientations,useGPU,gpuPrecision)
[ny_PSF,nx_PSF,nz_PSF,~,~] = size(inputOTF);
dk_PSF=1./([ny_PSF,nx_PSF,nz_PSF].*pxl_dim_PSF);

[ny_data,nx_data,nimgs_data] = size(inputData);
nz_data=nimgs_data/(nphases*norientations);
dk_data=1./([ny_data,nx_data,nz_data].*pxl_dim_data);

OTF_yy=[-ceil((ny_PSF-1)/2):floor((ny_PSF-1)/2)]*dk_PSF(1);
OTF_xx=[-ceil((nx_PSF-1)/2):floor((nx_PSF-1)/2)]*dk_PSF(2);
OTF_zz=[-ceil((nz_PSF-1)/2):floor((nz_PSF-1)/2)]*dk_PSF(3);

map_yy=[-ceil((ny_data-1)/2):floor((ny_data-1)/2)]*dk_data(1);
map_xx=[-ceil((nx_data-1)/2):floor((nx_data-1)/2)]*dk_data(2);
map_zz=[-ceil((nz_data-1)/2):floor((nz_data-1)/2)]*dk_data(3);

if(useGPU)
    OTF_yy = gpuArray(OTF_yy);
    OTF_xx = gpuArray(OTF_xx);
    OTF_zz = gpuArray(OTF_zz);

    map_yy = gpuArray(map_yy);
    map_xx = gpuArray(map_xx);
    map_zz = gpuArray(map_zz);
end

[OTF_xx_arr,OTF_yy_arr,OTF_zz_arr]=meshgrid(OTF_xx,OTF_yy,OTF_zz);
[map_xx_arr,map_yy_arr,map_zz_arr]=meshgrid(map_xx,map_yy,map_zz);

if(useGPU)
    OTF_scaled=gpuArray(zeros(size(map_xx_arr),gpuPrecision));
    inputOTF = gpuArray(inputOTF);
else
    OTF_scaled=zeros(size(map_xx_arr));
end

for jj=1:norientations
    for ii=1:norders
        OTF_scaled(:,:,:,ii,jj) = interp3(OTF_xx_arr,OTF_yy_arr,OTF_zz_arr,inputOTF(:,:,:,ii,jj),map_xx_arr,map_yy_arr,map_zz_arr,'linear',0+0i);
    end
end
end