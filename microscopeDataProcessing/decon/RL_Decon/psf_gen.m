function [ psf ] = psf_gen(psf_tiff_file, dz_psf, dz_data, cropToSize)
%psf_gen resample and crop raw PSF
%   

%load TIFF
psftiff=Tiff(psf_tiff_file,'r');

z=1;
while ~psftiff.lastDirectory
    psf_raw(:,:,z)=psftiff.read();
    z = z+1;
    psftiff.nextDirectory()
end
psf_raw(:,:,z)=psftiff.read();
nz = z;

% subtract background estimated from the last Z section
psf_raw = double(psf_raw) - mean(mean(psf_raw(:,:,nz)));
% convert all negative pixels to 0
psf_raw(psf_raw<0) = 0.0;

% locate the peak pixel
[peak, peakInd] = max(psf_raw(:));
[peakx,peaky,peakz] = ind2sub(size(psf_raw), peakInd);
% crop x-y around the peak
psf_cropped = psf_raw(peakx-cropToSize/2:peakx+cropToSize/2, ...
    peaky-cropToSize/2:peaky+cropToSize/2, :);

% center the PSF in z; otherwise RL decon results are shifted in z
psf_cropped=circshift(psf_cropped, [0, 0, round(nz/2-peakz)]);

[nx,ny,nz]=size(psf_cropped);
% resample PSF to match axial pixel size of raw image
dz_ratio = dz_data / dz_psf;
if dz_ratio > 1
    psf_fft=fftn(psf_cropped);
    new_nz = uint16(round(nz / dz_ratio));
    psf_fft_trunc=complex(zeros(nx,ny,new_nz));
    psf_fft_trunc(:,:,1:new_nz/2)=psf_fft(:,:,1:new_nz/2);
    psf_fft_trunc(:,:,new_nz-new_nz/2+1:new_nz)=psf_fft(:,:,nz-new_nz/2+1:nz);
    psf=real(ifftn(psf_fft_trunc));
else
    psf = psf_cropped;
end

psf(psf<0) = 0.0;
end

