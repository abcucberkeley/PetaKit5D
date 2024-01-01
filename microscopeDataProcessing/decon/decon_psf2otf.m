function otf = decon_psf2otf(psf, outSize)
%PSF2OTF Convert point-spread function to optical transfer function.

psfSize = size(psf);
if nargin < 2 
    outSize = psfSize;
end

if ~isempty(psf) % empty arrays are treated similar as in the fftn
  [psfSize, outSize] = decon_padlength(psfSize, outSize(:).');
  if any(outSize < psfSize)
    error('outSizeIsSmallerThanPsfSize')
  end
end

if ~all(psf(:)==0)
   
   % Pad the PSF to outSize
   padSize = double(outSize - psfSize);
   psf     = padarray(psf, padSize, 'post');

   % Circularly shift otf so that the "center" of the PSF is at the
   % (1,1) element of the array.
   psf    = circshift(psf,-floor(psfSize/2));

   % Compute the OTF
   otf = fftn(psf);

   % Estimate the rough number of operations involved in the 
   % computation of the FFT.
%    nElem = prod(psfSize);
%    nOps  = 0;
%    for k=1:ndims(psf)
%       nffts = nElem/psfSize(k);
%       nOps  = nOps + psfSize(k)*log2(psfSize(k))*nffts; 
%    end

   % Discard the imaginary part of the psf if it's within roundoff error.
%    if max(abs(imag(otf(:))))/max(abs(otf(:))) <= nOps*eps
%       otf = real(otf);
%    end
else
   otf = zeros(outSize, class(psf));
end

end


function [psz_1, psz_2] = decon_padlength(sz_1, sz_2)
%PADLENGTH Pad input vectors with ones to give them equal lengths.

numDims = 3;

psz_1 = [sz_1 ones(1,numDims-length(sz_1))];
psz_2 = [sz_2 ones(1,numDims-length(sz_2))];

end

