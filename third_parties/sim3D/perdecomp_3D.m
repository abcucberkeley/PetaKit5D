

%------------------------------ PERDECOMP ------------------------------


%       Periodic plus Smooth Image Decomposition

%

%               author: Lionel Moisan

%

%   This program is freely available on the web page

%

%   http://www.mi.parisdescartes.fr/~moisan/p+s

%

%   I hope that you will find it useful.

%   If you use it for a publication, please mention 

%   this web page and the paper it makes reference to.

%   If you modify this program, please indicate in the

%   code that you did so and leave this message.

%   You can also report bugs or suggestions to 

%   lionel.moisan [AT] parisdescartes.fr

%

% This function computes the periodic (p) and smooth (s) components

% of an image (3D array) u

%

% usage:    p = perdecomp(u)    or    [p,s] = perdecomp(u)

%

% note: this function also works for 1D signals (line or column vectors)

%

% v1.0 (01/2014): initial Matlab version from perdecomp.sci v1.2


function [p,s] = perdecomp_3D(u)


[ny,nx,nz] = size(u); 

u = double(u);

X = 1:nx; Y = 1:ny; Z = 1:nz;

v1 = zeros(ny,nx,nz);

v2 = zeros(ny,nx,nz);

v3 = zeros(ny,nx,nz);

v1(1,:,:) = u(1,:,:) - u(ny,:,:);

v1(ny,:,:) = -v1(1,:,:);

v2(:,1,:) = u(:,1,:) - u(:,nx,:);

v2(:,nx,:) = -v2(:,1,:);

v3(:,:,1) = u(:,:,1) - u(:,:,nz);

v3(:,:,nz) = -v3(:,:,1);

v = v1+v2+v3;

% v(1,X)  = u(1,X)-u(ny,X);

% v(ny,X) = -v(1,X);

% v(Y,1 ) = v(Y,1 )+u(Y,1)-u(Y,nx);

% v(Y,nx) = v(Y,nx)-u(Y,1)+u(Y,nx);

fxx = cos(2.*pi*(X -1)/nx);

fyy = cos(2.*pi*(Y-1)/ny);

fzz = cos(2.*pi*(Z-1)/nz);

[fx,fy,fz] = meshgrid(fxx,fyy,fzz);

fx(1,1,1)=0.;   % avoid division by 0 in the line below

s = real(ifftn(fftn(v)*0.5./(3.-fx-fy-fz)));

p = u-s;

