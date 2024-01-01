function [shiftedImage]=fourierShift3D(inputImage,shifts,useGPU,gpuPrecision)
%Performs sub-pixel real-space shifts of a 3D image by applying a phase
%ramp in Fourier space

%Inputs:
%inputImage is a 3D array of either real or complex values
%shifts is a 1 x 3 vector of shift values in real-space pixels

%Set up vectors based on image size
[ny,nx,nz]=size(inputImage);
kxx=cast([-ceil((nx-1)/2):floor((nx-1)/2)],gpuPrecision);
kyy=cast([-ceil((ny-1)/2):floor((ny-1)/2)],gpuPrecision);
kzz=cast([-ceil((nz-1)/2):floor((nz-1)/2)],gpuPrecision);

if(useGPU)
    kxx = gpuArray(kxx);
    kyy = gpuArray(kyy);
    kzz = gpuArray(kzz);
    inputImage = gpuArray(inputImage);
end

%Set up arrays for frequency space pixels
[kxx_arr,kyy_arr,kzz_arr]=meshgrid(kxx/nx,kyy/ny,kzz/nz);

%FFT image, apply phase ramp, and inverse transform
Fimage=fftshift(fftn((inputImage)));
Fimage_shift=Fimage.*exp(-1i*2*pi*(kyy_arr.*shifts(1)+kxx_arr.*shifts(2)+kzz_arr.*shifts(3)));
shiftedImage=ifftn(ifftshift(Fimage_shift));

end