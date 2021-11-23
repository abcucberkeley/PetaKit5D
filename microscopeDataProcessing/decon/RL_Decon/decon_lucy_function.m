function [J_2, err_mat, k] = decon_lucy_function(I, PSF, NUMIT, fixIter, err_thrsh, debug, debug_folder, useGPU)
% adapted from matlab deconvlucy.m
% 
% xruan (05/18/2021): add support for early stop with stop criteria
% xruan (05/19/2021): add support for gpu computing
% xruan (06/02/2021): add support for fix iteration decon
% xruan (06/03/2021): add support for debug (evaluate error and save
% intermediate result).
% xruan (07/09/2021): add output for the number of iterations run; also
% output err mat for fixed iterations.
% xruan (11/11/2021): when using GPU, use single precision unless the data
% is in double
% xruan (11/12/2021): refactor code to save memory (J_4), and disable
% weighting or subsampling

% switch nargin
%     case 3 %                 deconvlucy(I,PSF,NUMIT)
%         NUMIT = varargin{1};
%     case 4 %                 deconvlucy(I,PSF,NUMIT,DAMPAR)
%         NUMIT = varargin{1};
%         DAMPAR = varargin{2};
%     case 5 %                 deconvlucy(I,PSF,NUMIT,DAMPAR,WEIGHT)
%         NUMIT = varargin{1};
%         DAMPAR = varargin{2};
%     case 6 %                 deconvlucy(I,PSF,NUMIT,DAMPAR,WEIGHT,READOUT)
%         NUMIT = varargin{1};
%         DAMPAR = varargin{2};
%         READOUT = varargin{3};
%     case 7 %                 deconvlucy(I,PSF,NUMIT,DAMPAR,WEIGHT,READOUT,SUBSMPL)
%         NUMIT = varargin{1};
%         DAMPAR = varargin{2};
%         READOUT = varargin{3};
%         SUBSMPL = varargin{4};
% end

if nargin < 4
    fixIter = false;
end

if nargin < 5 || isempty(err_thrsh)
    err_thrsh = 5e-11;
end

if nargin < 6
    debug = false;
end

if debug
    fixIter = true;
end

if nargin < 7
    debug_folder = './debug/';
end

% option to use gpu
if nargin < 8
    useGPU = true;
end
useGPU = useGPU & gpuDeviceCount > 0;

if useGPU && ~isa(I, 'double')
    I = single(I);
else
    I = double(I);
end
classI = class(I);
sizeI = size(I);
sizePSF = size(PSF);

J_1 = I;
J_2 = J_1;
J_3 = zeros(sizeI, classI);
% J{4} = 0;

NUMIT = cast(NUMIT, classI);
DAMPAR = 0;
% WEIGHT = ones(sizeI);
WEIGHT = 1;
READOUT = 0;
SUBSMPL = 1;
% [sizeI, sizePSF] = padlength(size(I), size(PSF));
numNSdim = 3;

% if length(J)==3 % assign the 4-th element of input cell now
%     J{4}(prod(sizeI)*SUBSMPL^length(numNSdim),2) = 0;
% end;

% J_4 = J_1;
J_4 = zeros(prod(sizeI)*SUBSMPL^length(numNSdim),1,classI);

if useGPU
    PSF = gpuArray(PSF);
    J_1 = gpuArray(J_1);
    J_2 = gpuArray(J_2);
    J_3 = gpuArray(J_3);
    J_4 = gpuArray(J_4);
end

% 1. Prepare PSF. If PSF is known at a higher sampling rate, it has to be
% padded with zeros up to sizeI(numNSdim)*SUBSMPL in all non-singleton
% dimensions. Or its OTF could take care of it:
% sizeOTF = sizeI;
% sizeOTF(numNSdim) = SUBSMPL*sizeI(numNSdim);
% test normalization of PSF first
PSF = PSF ./ sum(PSF(:));
H = decon_psf2otf(PSF,double(sizeI));
scale = sum(PSF(:));
clear PSF

% 2. Prepare parameters for iterations
%
% Create indexes for image according to the sampling rate
% xruan: redefine subsample factors as separate ones
% idx = repmat({':'},[1 length(sizeI)]);
% for k = numNSdim % index replicates for non-singleton PSF sizes only
%     idx{k} = reshape(repmat(1:sizeI(k),[SUBSMPL 1]),[SUBSMPL*sizeI(k) 1]);
% end
% idx_y = reshape(repmat(1:sizeI(1),[SUBSMPL 1]),[SUBSMPL*sizeI(1) 1]);
% idx_x = reshape(repmat(1:sizeI(2),[SUBSMPL 1]),[SUBSMPL*sizeI(2) 1]);
% idx_z = reshape(repmat(1:sizeI(3),[SUBSMPL 1]),[SUBSMPL*sizeI(3) 1]);


wI = max(WEIGHT.*(READOUT + J_1),0);% at this point  - positivity constraint
% J_1 = J_1 .* (J_1 > 0);
% J_2 = J_2(idx_y, idx_x, idx_z);
% scale = real(ifftn(conj(H).*fftn(WEIGHT(idx_y, idx_x, idx_z)))) + sqrt(eps);
clear J_1 WEIGHT;

% DAMPAR22 = (DAMPAR.^2)/2;

if SUBSMPL~=1 % prepare vector of dimensions to facilitate the reshaping
    % when the matrix is binned within the iterations.
    vec = zeros(1, numNSdim * 2);
    vec(2:2:2*length(sizeI)) = sizeI;
    vec(2*numNSdim-1) = -1;
    vec(vec==0) = [];
    num = fliplr(find(vec==-1));
    vec(num) = SUBSMPL;
else
    vec = [];
    num = [];
end

% 3. L_R Iterations
%
estep = 5;
err_mat = zeros(NUMIT, 4);

lambda = 2*any(J_4(:)~=0);
for k = lambda + 1 : lambda + NUMIT
    
    % 3.a Make an image predictions for the next iteration
    if k > 2
        % lambda = (J_4(:,1).'*J_4(:,2))/(J_4(:,2).'*J_4(:,2) +eps);        
        lambda = ((J_2(:) - Y(:)).'*J_4)/(J_4.'*J_4 +eps);
        lambda = max(min(lambda,1),0);% stability enforcement
        J_4 = J_2(:) - Y(:);
    elseif k == 2
        J_4 = J_2(:) - Y(:);
    end
    Y = max(J_2 + lambda*(J_2 - J_3),0);% plus positivity constraint
        
    % 3.b  Make core for the LR estimation
    % CC = corelucy(Y,H,DAMPAR22,wI,READOUT,SUBSMPL,vec,num);
    % directly compute CC within the same function to reduce overhead 
    ReBlurred = real(ifftn(H.*fftn(Y)));
    % ReBlurred = ReBlurred + READOUT;
    % ReBlurred(ReBlurred == 0) = eps;
    ReBlurred = ReBlurred + (ReBlurred == 0) * eps;
    ReBlurred = wI./ReBlurred + eps;
    
    % 3.c Determine next iteration image & apply positivity constraint
    J_3 = J_2;
    J_2 = max(Y.*real(ifftn(conj(H).*fftn(ReBlurred)))./scale,0);
    
    % clear CC;
    % CC = 0;
    % J_4 = [J_2(:)-Y(:) J_4(:,1)];
    % J_4 = flip(J_4, 2);
    % J_4(:, 1) = J_2(:)-Y(:);
    
    if ~debug && rem(k, estep) == 0
        istp = k/estep;
        if ~useGPU
            % err_mat(istp) = sum((J_2 - I) .^ 2, 'all') / numel(J_2);
            err_mat(istp, 1:2) = [k, sum((J_2 - I) .^ 2, 'all')];
            err_mat(istp, 3) = err_mat(istp, 2) ./ err_mat(1, 2);
            if istp > 2
                err_mat(istp, 4) = min(abs(err_mat(istp, 3) - err_mat(istp-1, 3)) / 10, ...
                    abs(err_mat(istp, 3) + err_mat(istp-2, 3) - 2 * err_mat(istp-1, 3)) / 2);
            else
                err_mat(istp, 4) = 1;
            end

            if ~fixIter && k > estep * 2 && err_mat(istp, 4) < err_thrsh
                break;
            end
        end
    end
    
    if debug
        err_mat(k, 1:2) = [k, mean((double(J_2) - double(I)) .^ 2, 'all')];
        err_mat(k, 3) = err_mat(k, 2) ./ err_mat(1, 2); 
        if rem(k, estep) == 0 && k > estep * 2
            istp = k/estep;
            err_mat(k, 4) = min(abs(err_mat(istp, 3) - err_mat(istp-1, 3)) / 10, abs(err_mat(istp, 3) + err_mat(istp-2, 3) - 2 * err_mat(istp-1, 3)) / 2);
        else
            if k < 3 * estep
                err_mat(k, 4) = 1;
            else
                err_mat(k, 4) = err_mat(istp * estep, 4);
            end
        end
        
        if k == 1
            mkdir(debug_folder);
        end
        if rem(k, 5) == 0
            if useGPU 
                writetiff(single(gather(J_2)), sprintf('%s/Iter_%04d.tif', debug_folder, k));
            else
                writetiff(single(J_2), sprintf('%s/Iter_%04d.tif', debug_folder, k));
            end
        end
    end
end

if ~fixIter
    err_mat = err_mat(1 : istp, :);
end

if useGPU
    J_2 = gather(J_2);
end

% clear wI H scale Y;

% 4. Convert the right array (for cell it is first array, for notcell it is
% second array) to the original image class & output whole thing
% num = 1 + strcmp(classI{1},'notcell');
% if ~strcmp(classI{2},'double')
%     J{num} = images.internal.changeClass(classI{2},J{num});
% end

% J = images.internal.changeClass(classI, J_2);


end


function f = corelucy(Y,H,DAMPAR22,wI,READOUT,SUBSMPL,vec,num)
%CORELUCY Accelerated Damped Lucy-Richarson Operator.
%  Calculates function that when used with the scaled projected array 
%  produces the next iteration array that maximizes the likelihood that 
%  the entire suite satisfies the Poisson statistics. 
%
% See also DECONVLUCY and DECONVBLIND.

%  Copyright 1993-2003 The MathWorks, Inc.  

%   References
%   ----------
%   "Acceleration of iterative image restoration algorithms, by D.S.C. Biggs 
%   and M. Andrews, Applied Optics, Vol. 36, No. 8, 1997.
%   "Deconvolutions of Hubble Space Telescope Images and Spectra",
%   R.J. Hanisch, R.L. White, and R.L. Gilliland. in "Deconvolution of Images 
%   and Spectra", Ed. P.A. Jansson, 2nd ed., Academic Press, CA, 1997.

ReBlurred = real(ifftn(H.*fftn(Y)));

% 1. Resampling if needed
if SUBSMPL ~= 1 % Bin ReBlurred back to the sizeI for non-singleton dims
  
  %1.Reshape so that the-to-binned dimension separates into two
  %dimensions, with one of them consisting of elements of a single bin.
  ReBlurred = reshape(ReBlurred,vec);

  %2. Bin (==calculate mean) along the first of the-to-binned dimension,
  %that dimension consists of the bin elements. Reshape to get rid off
  for k = num % new appeared singleton.
    vec(k) = [];
    ReBlurred = reshape(mean(ReBlurred,k),vec);
  end
  
end

% 2. An Estimate for the next step
ReBlurred = ReBlurred + READOUT;
% ReBlurred(ReBlurred == 0) = eps;
ReBlurred = ReBlurred + (ReBlurred == 0) * eps;
ReBlurred = wI./ReBlurred + eps;

% 3. Damping if needed
if DAMPAR22 == 0 % No Damping
  % ImRatio = AnEstim(idx_y, idx_x, idx_z);
  % ImRatio = AnEstim;
else % Damping of the image relative to DAMPAR22 = (N*sigma)^2
%   gm = 10;
%   g = (wI.*log(AnEstim)+ ReBlurred - wI)./DAMPAR22;
%   g = min(g,1);
%   G = (g.^(gm-1)).*(gm-(gm-1)*g);
%   % ImRatio = 1 + G(idx_y, idx_x, idx_z).*(AnEstim(idx_y, idx_x, idx_z) - 1);
%   ImRatio = 1 + G.*(AnEstim - 1);
end

f = fftn(ReBlurred);

end

function otf = decon_psf2otf(psf, outSize)

%PSF2OTF Convert point-spread function to optical transfer function.
psfSize = size(psf);
if isempty(outSize)
  outSize = psfSize;
elseif ~isempty(psf) % empty arrays are treated similar as in the fftn
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
   nElem = prod(psfSize);
   nOps  = 0;
   for k=1:ndims(psf)
      nffts = nElem/psfSize(k);
      nOps  = nOps + psfSize(k)*log2(psfSize(k))*nffts; 
   end

   % Discard the imaginary part of the psf if it's within roundoff error.
   if max(abs(imag(otf(:))))/max(abs(otf(:))) <= nOps*eps
      otf = real(otf);
   end
else
   otf = zeros(outSize, class(psf));
end

end


function [psz_1, psz_2] = decon_padlength(sz_1, sz_2)
%PADLENGTH Pad input vectors with ones to give them equal lengths.
%
%   Example
%   -------
%       [a,b,c] = padlength([1 2],[1 2 3 4],[1 2 3 4 5])
%       a = 1 2 1 1 1
%       b = 1 2 3 4 1
%       c = 1 2 3 4 5

%   Copyright 1993-2003 The MathWorks, Inc.  

% Find longest size vector.  Call its length "numDims".
% numDims = zeros(nargin, 1);
% for k = 1:nargin
%     numDims(k) = length(varargin{k});
% end
% numDims = max(numDims);
numDims = 3;

% Append ones to input vectors so that they all have the same length;
% assign the results to the output arguments.
% limit = max(1,nargout);
% varargout = cell(1,limit);
% for k = 1 : limit
%     pSize = [varargin{k} ones(1,numDims-length(varargin{k}))];
% end

psz_1 = [sz_1 ones(1,numDims-length(sz_1))];
psz_2 = [sz_2 ones(1,numDims-length(sz_2))];

end
