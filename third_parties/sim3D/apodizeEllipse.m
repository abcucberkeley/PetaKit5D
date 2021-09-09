function [Data_k_space_apodized,apo_ellipse] = apodizeEllipse(Data_k_space,dk_data,p_vec_guess,angle,NA_det,nimm,NA_ext_max,wvl_em,wvl_ext,islattice,useGPU)

alpha_det=asin(NA_det/nimm); %Half angle for detection NA
alpha_ext=asin(NA_ext_max/nimm); %Half angle for excitation NA
kmag_det=nimm/wvl_em; %k-vector magnitude for detection
kmag_ext=nimm/wvl_ext; %k-vector magnitude for excitation

%Fill out k-space array for super-sampled data
[ny_datass,nx_datass,nz_datass]=size(Data_k_space);
kyy_ss=[-ceil((ny_datass-1)/2):floor((ny_datass-1)/2)]*dk_data(1);
kxx_ss=[-ceil((nx_datass-1)/2):floor((nx_datass-1)/2)]*dk_data(2);
kzz_ss=[-ceil((nz_datass-1)/2):floor((nz_datass-1)/2)]*dk_data(3);

if(useGPU)
    kyy_ss = gpuArray(kyy_ss);
    kxx_ss = gpuArray(kxx_ss);
    kzz_ss = gpuArray(kzz_ss);
end

[kxx_ss_arr,kyy_ss_arr,kzz_ss_arr]=meshgrid(kxx_ss,kyy_ss,kzz_ss);
    
%Generate apodization filter
%Horizontal major axis of the ellipse
a=max(sqrt((p_vec_guess(:,1)*dk_data(1)).^2+(p_vec_guess(:,2)*dk_data(2)).^2))+kmag_det*sin(alpha_det)*2;

%Horizontal minor axis of the ellipse
b=kmag_det*sin(alpha_det)*2;

%Vertical axis of the ellipse
if islattice
c=kmag_ext*sin(alpha_ext)*2+kmag_det*(1-cos(alpha_det));
else
%3beam SIM
c=kmag_ext*(1-cos(alpha_ext))*2+kmag_det*(1-cos(alpha_det));
end

%Rotate apodization filter to match the pattern orientation
rot_coords=[cos(angle),-sin(angle),0;sin(angle),cos(angle),0;0,0,1]*[kxx_ss_arr(:)';kyy_ss_arr(:)';kzz_ss_arr(:)'];
apo_ellipse=-((rot_coords(1,:).^2/a^2+rot_coords(2,:).^2/b^2+rot_coords(3,:).^2/c^2))+1;
apo_ellipse(apo_ellipse<0)=0;
apo_ellipse=reshape(apo_ellipse,ny_datass,nx_datass,nz_datass);

%Apply apodization filter
Data_k_space_apodized=Data_k_space.*apo_ellipse;

end

