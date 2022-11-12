function psf = decon_otf2psf(otf, outSize)
% decon_otf2psf Convert optical transfer function to point-spread function.

otf = double(otf);
otfSize = size(otf);
if isempty(otfSize)
  outSize = otfSize;
end

if ~all(otf(:)==0)
   
   % Pad the PSF to outSize
   cropSize = double(otfSize - outSize);
   % psf     = padarray(psf, padSize, 'post');
   psf = real(ifftn(otf));

   % Circularly shift otf so that the "center" of the PSF is at the
   % (1,1) element of the array.
   psf    = circshift(psf,-floor(cropSize/2));
   psf = psf(1 : outSize(1), 1 : outSize(2), 1 : outSize(3));
   
else
   psf = zeros(outSize, class(otf));
end

end

