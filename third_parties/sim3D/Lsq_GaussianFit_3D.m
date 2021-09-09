function fitParams = Lsq_GaussianFit_3D(PeakImg,displayFit)
%Performs least square fitting to a 3D stack of a bead
%Determines centroid position in units of pixels

[ny,nx,nz]=size(PeakImg);
%Generate pixel coordinates for Gaussian fitting
[xx,yy,zz]=meshgrid(1:nx,1:ny,1:nz);
PeakData=[PeakImg(:),yy(:),xx(:),zz(:)];

%Use the value and location of the max voxel as initial guesses
[maxval,ind]=max(PeakImg(:));
[cy,cx,cz]=ind2sub([ny,nx,nz],ind);

params_guess(1) = maxval;  %max value within the window
params_guess(2) = cy;  %y pixel value of peak pixel in the window
params_guess(3) = cx;  %x pixel value of peak pixel in the window
params_guess(4) = cz;  %zo(guess) = z pixel value of peak pixel in the window
params_guess(5) = 1.5; %sigy(guess) = y sigma value
params_guess(6) = 1.5; %sigx(guess) = x sigma value
params_guess(7) = 1.5; %sigz(guess) = z sigma value

%find the initial best fit with unadjusted weights:
options = optimset('Display', 'off', 'MaxIter', 10000, 'MaxFunEvals', 10000, 'TolX', 1e-6, 'TolFun', 1e-6);
[fitParams, resnorm] = fminsearch(@(params) objFun_GaussianFit_3D(params,PeakData), params_guess,options);

%% display results----------------------------------------------------------------
if displayFit
    F=plot3DGaussian(fitParams,PeakData,[ny,nx,nz]);
    figure
    subplot(2,2,1)
    imagesc(squeeze(max(reshape(PeakData(:,1),[ny,nx,nz]),[],3)))
    axis equal
    title('XY projection of the input image')
    subplot(2,2,2)
    imagesc(squeeze(max(reshape(PeakData(:,1),[ny,nx,nz]),[],1))')
    axis equal
    title('XZ projection of the input image')
    subplot(2,2,3)
    imagesc(squeeze(max(F,[],3)))
    axis equal
    title('XY projection of the fit image')
    subplot(2,2,4)
    imagesc(squeeze(max(F,[],1))')
    axis equal
    title('XZ projection of the fit image')
end
end

%% --------------------------------------------------------------------------------------------
%Helper functions below
%Calculate the objective function for the least-square fit
function F = objFun_GaussianFit_3D(params,PeakData)
F = sum(((PeakData(:,1) - (params(1) .* exp(-((PeakData(:,2) - params(2)).^2./(2*params(5)^2) + (PeakData(:,3) - params(3)).^2./(2*params(6)^2) + (PeakData(:,4) - params(4)).^2./(2*params(7)^2))))).^2));
end

%Reshape fitdata for display
function im = plot3DGaussian(x,PeakData,dims)
im = (x(1) .* exp(-((PeakData(:,2) - x(2)).^2./(2*x(5)^2) + (PeakData(:,3) - x(3)).^2./(2*x(6)^2) + (PeakData(:,4) - x(4)).^2./(2*x(7)^2))));
im=reshape(im,dims(1),dims(2),dims(3));
end