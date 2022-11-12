function [J_2, err_mat, k] = decon_lucy_function(I, PSF, NUMIT, fixIter, err_thrsh, debug, debug_folder, saveStep, useGPU)
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

if nargin < 8
    saveStep = 5;
end

% option to use gpu
if nargin < 9
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
estep = saveStep;
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
        if rem(k, estep) == 0
            if useGPU 
                writetiff(uint16(gather(J_2)), sprintf('%s/Iter_%04d.tif', debug_folder, k));
            else
                writetiff(uint16(J_2), sprintf('%s/Iter_%04d.tif', debug_folder, k));
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


end

