function [xy_exp_PSF, xz_exp_PSF, yz_exp_PSF, xy_exp_OTF, xz_exp_OTF, yz_exp_OTF, xOTF_linecut, yOTF_linecut, zOTF_linecut, zOTF_bowtie_linecut, zOTF_bowtie_linecut_yz] = ...
    Load_and_Plot_Exp_Overall_xzPSF_xzOTF_update(filenm, source_descrip, xypixsize, zpixsize, NAdet, index, exc_lambda, det_lambda, PSFsubpix, gamma, bgFactor)
%
%LOAD_AND_PLOT_EXP_OVERALL_xzPSF_xzOTF  Loads a 3D TIFF stack of an experimentally
%measured overall PSF, and plots it along with the OTF determined by its
%Fourier transform
%
%   filenm is a string with the full path to the TIFF file
%   source_descrip is a string giving a description of the source file data
%   xypixsize is the xy pixel size of the data in nm
%   zpixsize is the z pixel size of the data in nm
%   NA = numerical aperture of the objective
%   index = refractive index of the medium
%   exc_lambda = excitation wavelength in free space
%   det_lambda = detection wavelength in free space
%   0< gamma < 1 is the non-linear contrast factor for secondary plots of
%       the PSF and OTF
%   
%   xz_exp_PSF is a 2D slice through the measured PSF, rescaled to the
%       pixel size and FOV used in the overall PSF/OTF simulations
%   xz_exp_OTF is a 2D slice through the measured OTF, rescaled to the
%       pixel size and FOV used in the overall PSF/OTF simulations
%   xOTF_linecut is a 1D linecut of the kx variation in the OTF along the 
%       line where ky=kz=0, at the pixel size used in the PSF/OTF simulations
%   zOTF_bowtie linecut is a 1D linecut of the kz variation in the OTF along the 
%       line where the kz oTF bowite is the fattest, i.e., kx=kabbe/2 and
%       ky = 0, at the pixel size used in the PSF/OTF simulations
%
% xruan: unify normalization factor, fix crop issue, add user input psf sub
% size. 
% xuran: make the peak to the center (after subtracting background). 
% xruan (07/23/2021): correct xz axis as yz for figures; add actual xz
% pannels.
% xruan (04/09/2022): output all central PSFs/OTFs and linecuts

%plot parameters from simulation code:
sim_pixsize = 0.1;   %plot xz pix size in media excitation wavelengths
sim_half_FOV = 50;   %half width of the xz plot FOV
MaxCalcK = 0.25./0.1; %calc the OTF from -MaxCalcK to +MaxCalcK in excitation 
                      %media wavelngths
% PSFsubpix = [128 128 101];  %region to be extracted from the raw PSF data for
                            %calculation of the OTF -- thus making the
                            %DC scaling at the center peak the same for all data sets 
% xruan: change PSFsubpix as an input
if nargin < 9
    PSFsubpix = [128, 128, 101];
end
if nargin < 10
    gamma = 0.5;
end
if nargin < 11
    % bgFactor = 1.15;
    bgFactor = 1.5;
end
%
%hardwire the filename location for now:
%filenm = 'C:\Users\betzige\Dropbox (HHMI)\HexLLS_CF0p02_FF1_complete_benchmark\complete_benchmark\Hex\totalPSF\560nm_Hex_CFp02_FF1_p55p40_c-7p5um_PSF.tif';
%filenm = 'C:\Users\betzige\Dropbox (HHMI)\For Gokul and Gaoxiang\20210223_WF-PSF\488nm\488nm_WF-PSF_bead1.tif';
%filenm = 'C:\Users\betzige\Dropbox (HHMI)\For Gokul and Gaoxiang\20210223_WF-PSF\488nm\488nm_WF-PSF_RotatedDO_bead1.tif';
% filenm = 'C:\AO-SwAK\LLS calcs Jan 2021\detection PSF and OTF\NA1p0\experimental OTFs\PSF 488_450ms from Gaoxiang.tif';
%
%load all images in the TIF stack to a single 3D matrix:
warning('off','all') % Suppress all the tiff warnings
tstack  = Tiff(filenm);
[I,J] = size(tstack.read());
K = length(imfinfo(filenm));
PSF3D = zeros(I,J,K);
PSF3D(:,:,1)  = tstack.read();
for n = 2:K
    tstack.nextDirectory()
    PSF3D(:,:,n) = tstack.read();
end
%%convert from integer type to double:
%PSF3D = woof;
PSF3Dexp = cast(PSF3D, 'double');
A = size(PSF3Dexp);
%
%use the first 10 x/y rows/columns to estimate and subtract the dark background;
% backgnd = sum(sum(sum(PSF3Dexp(1:10,1:10,:))))./(100.*A(3));
% backgnd = 1.5*mode(PSF3Dexp(:));
backgnd = bgFactor * median(PSF3Dexp(PSF3Dexp > 0));
PSF3Dexp = max(PSF3Dexp-backgnd,0);

% xruan: pad the xy to make them same size, and also to PSFsubpix size if
% they are too small
size_xy = max([PSFsubpix(1), PSFsubpix(2), A(1), A(2)]);
hsz = [ceil((size_xy - A(1) - 1) / 2), ceil((size_xy - A(2) - 1) / 2)];
psize_1 = [hsz(1), hsz(2), 0];
psize_2 = [size_xy - A(1) - hsz(1), size_xy - A(2) - hsz(2), 0];

% pad z if the number of slices is fewer than the requirement
if A(3) < PSFsubpix(3)
    hsz = ceil((PSFsubpix(3) - A(3) -1) / 2);
    psz = [hsz, abs(PSFsubpix(3) - A(3)) - hsz];
    psize_1(3) = psz(1);
    psize_2(3) = psz(2);
end    
    
PSF3Dexp = padarray(PSF3Dexp, psize_1, 0, 'pre');
PSF3Dexp = padarray(PSF3Dexp, psize_2, 0, 'post');
A = size(PSF3Dexp);

% xruan: shift peak to the center
[peak, peakInd] = max(PSF3Dexp(:));
[peaky,peakx,peakz] = ind2sub(A, peakInd);
PSF3Dexp=circshift(PSF3Dexp, round((A + 1) / 2 - [peaky,peakx,peakz]));

%
%resize to the pixel sizes used in the PSF/OTF simulations:
%find the measured x and z pixel sizes in media excitation wavelengths:
sim_xzpixsize = sim_pixsize.*(exc_lambda./index);  %simulation pix size in nm
magxy = xypixsize./sim_xzpixsize; 
magz = zpixsize./sim_xzpixsize;
%find the plot dimensions after resizing the data to simulation pixel size:
B(1) = round(magxy.*A(1));
B(2) = round(magxy.*A(2));
B(3) = round(magz.*A(3));
% [x,y,z] = ndgrid(1:A(1),1:A(2),1:A(3));
% [xq,yq,zq] = ndgrid(1:(A(1)./B(1)):A(1),1:(A(2)./B(2)):A(2),1:(A(3)./B(3)):A(3));
% PSF3D = interpn(x,y,z, PSF3Dexp, xq,yq,zq);
% PSF3D = PSF3D./max(max(max(PSF3D)));
PSF3D = imresize3(PSF3Dexp, B, 'linear');
PSF3D = PSF3D./max(PSF3D(:));
B = size(PSF3D);
%
%put this data into an array of size equal to that used for the simulations,
%padding with zeros as necessary:
array_size = 2.*sim_half_FOV./sim_pixsize + 1;
PSF_array = zeros(array_size, array_size, array_size);
xyi = max(round(sim_half_FOV./sim_pixsize-B(1)./2),1);
xyf = xyi+min(B(1), array_size)-1;
% xyf = xyi+array_size-1;
zi = max(round(sim_half_FOV./sim_pixsize-B(3)./2),1);
zf = zi+min(B(3), array_size)-1;
% zf = zi+array_size-1;
Bc = round((B + 1) / 2);
lhsz = ceil((min(B, array_size) - 1) / 2);
rhsz = min(B, array_size) - lhsz - 1;

PSF_array(xyi:xyf,xyi:xyf,zi:zf) = PSF3D(Bc(1) - lhsz(1) : Bc(1) + rhsz(1), Bc(2) - lhsz(2) : Bc(2) + rhsz(2), Bc(3) - lhsz(3) : Bc(3) + rhsz(3));
% PSF_array(xyi:xyf,xyi:xyf,zi:zf) = PSF3D;
PSF3D = PSF_array;

A = size(PSF3D);

% xruan: fine adjust the center by the centroid
% BW = PSF3D > prctile(PSF3D(:), 99.99);
BW = PSF3D > multithresh(PSF3D(PSF3D > 0), 1);
% BW = imclose(BW, strel('sphere', 1));
% L = bwlabeln(BW);
py = (A(1) + 1) / 2;
px = (A(2) + 1) / 2;
pz = (A(3) + 1) / 2;
% BW = L == L(round(py), round(px), round(pz));
% % calculate centroid
% [y, x, z] = ind2sub(A, find(BW));
% f = PSF3D(BW);
% cy = y' * f / sum(f);
% cx = x' * f / sum(f);
% cz = z' * f / sum(f);
% sigma = [1.3 * magxy, 3 * magz];
% 
% [pstruct, mask] = XR_pointSourceDetection3D(PSF3D, sigma, 'Mode', 'xyzAc', 'Alpha', .05,...
%     'Mask', BW, 'RemoveRedundant', true, 'RefineMaskLoG', false, ...
%     'FitGaussianMethod', 'original', 'BackgroundCorrection', false);
% idx = find(pstruct.A==max([pstruct.A]));
% % sigmaXY = pstruct.s(1,idx);
% % sigmaZ = pstruct.s(2,idx);
% cx = pstruct.x(idx);
% cy = pstruct.y(idx);
% cz = pstruct.z(idx);
% change to use peak
PSF3D_1 = imgaussfilt3(PSF3D, 2);
[peak, peakInd] = max(PSF3D_1(:));
clear PSF3D_1;
[interp_peaky,interp_peakx,interp_peakz] = ind2sub(A, peakInd);
cx = interp_peakx;
cy = interp_peaky;
cz = interp_peakz;
% make centroid within [-0.5, 0.5] to the center
if abs(py - cy) > 0.5 || abs(px - cx) > 0.5 || abs(pz - cz) > 0.5
    PSF3D=circshift(PSF3D, round((A + 1) / 2 - [cy,cx,cz]));
    [peak, peakInd] = max(PSF3D(:));
    [interp_peaky,interp_peakx,interp_peakz] = ind2sub(A, peakInd);
end

%
%find the location of the maximum of the PSF:
%assume the maximum occurs at the center z plane:
midpix = (A(3)+1)./2;
PSFxy = squeeze(PSF3D(:,:,midpix));
PSFxy = PSFxy./max(max(PSFxy));
PSFy = sum(PSFxy,1);
[maxy, peakypix] = max(PSFy);  %peakypix is the y pixel value at whic hte PSF is maximum
%
%calculate the 3D OTF:
%pick out the subset of pixels in the PSF used for comparisons across
%different measurements:
C = size(PSF3Dexp);
pixi = max(round((C-PSFsubpix)./2),1);
PSFsubset = PSF3Dexp(pixi(1):(pixi(1)+PSFsubpix(1)-1),pixi(2):(pixi(2)+PSFsubpix(2)-1),pixi(3):(pixi(3)+PSFsubpix(3)-1));
OTF3Dexp = abs(fftshift(fftn(PSFsubset)));
maxval = max(max(max(OTF3Dexp)));
OTF3Dexp = OTF3Dexp ./ maxval;
C = size(OTF3Dexp);
%
%resize to the pixel sizes used in the OTF simulations:
%find the measured x and z half FOV in media excitation wavelengths:
xy_half_FOV = C(1).*xypixsize./2;  %xy half FOV in nm
z_half_FOV = C(3).*zpixsize./2;    %z half FOV in nm
sim_half_FOVnm = sim_half_FOV.*(exc_lambda./index);  %simulation half FOV in nm
xyOTFpixsize = 1./xy_half_FOV;  %OTF pixel size inversely proportional to FOV
zOTFpixsize = 1./z_half_FOV;
simOTFpixsize = 1./sim_half_FOVnm;
magxy = xyOTFpixsize./simOTFpixsize; 
magz = zOTFpixsize./simOTFpixsize;
%find the plot dimensions after resizing the data to simulation pixel size:
D(1) = round(magxy.*C(1));
D(2) = round(magxy.*C(2));
D(3) = round(magz.*C(3));
[x,y,z] = ndgrid(1:C(1),1:C(2),1:C(3));
[xq,yq,zq] = ndgrid(1:(C(1)./D(1)):C(1),1:(C(2)./D(2)):C(2),1:(C(3)./D(3)):C(3));
OTF3D = interpn(x,y,z, OTF3Dexp, xq,yq,zq);
OTF3D = OTF3D./max(max(max(OTF3D)));
D = size(OTF3D);
%
%put this data into an array of size equal to that used for the simulations,
%padding with zeros as necessary:
OTF_array = zeros(array_size, array_size, array_size);
xyi = max(round(sim_half_FOV./sim_pixsize-D(1)./2),1);
% xyf = xyi+D(1)-1;
xyf = xyi+min(D(1), array_size)-1;
zi = max(round(sim_half_FOV./sim_pixsize-D(3)./2),1);
% zf = zi+D(3)-1;
zf = zi+min(D(3), array_size)-1;
Dc = round((D + 1) / 2);
lhsz = ceil((min(D, array_size) - 1) / 2);
rhsz = min(D, array_size) - lhsz - 1;
% OTF_array(xyi:xyf,xyi:xyf,zi:zf) = OTF3D;
OTF_array(xyi:xyf,xyi:xyf,zi:zf) = OTF3D(Dc(1) - lhsz(1) : Dc(1) + rhsz(1), Dc(2) - lhsz(2) : Dc(2) + rhsz(2), Dc(3) - lhsz(3) : Dc(3) + rhsz(3));
OTF3D = OTF_array;
B = size(OTF3D);
%
%loop through the various y slices of the 3D PSF to find the brightest one
% that will be deemed the central y slice form which to extract the xz overall PSF:
% note: xz should be yz actually here
maxval = zeros(1,11);
for yoffset = -5:1:5
    yplane = peakypix + yoffset;
    %extract the xz PSF at the current y plane:
    yz_exp_PSF = squeeze(PSF3D(:,yplane,:))';
    maxval(yoffset+6) = max(max(yz_exp_PSF));
end
[totmax,maxplane] = max(maxval);
maxyplane = peakypix + maxplane - 6;
yz_exp_PSF = squeeze(PSF3D(:,maxyplane,:))'; 
hsz = round((size(yz_exp_PSF) + 1) / 2);
% xruan: narrow down the range to [-50, 50]
yz_exp_PSF = yz_exp_PSF(hsz(1) - 50 : hsz(1) + 50, hsz(2) - 50 : hsz(2) + 50);
A = size(yz_exp_PSF);
%

fig = figure('Renderer', 'painters', 'Position', [10 10 2500 1100]);
%calc and plot the experimental overall PSF:
% figure  %create a new figure window for the plots
% actually yz psf
subplot(2,5,4)
% set(gcf, 'Position', [550 100 600 600]);
% axes_h = axes('Position', [0.13, 0.1, 0.8, 0.8]);
image(256 .* yz_exp_PSF);
colormap hot(256);
axis([1 A(1) 1 A(2)]);
axis square;
set(gca, 'XTick', [1:(A(1)-1)./4:A(1)]);
% set(gca, 'XTickLabel', sim_half_FOV .* [-1:0.5:1]);
set(gca, 'XTickLabel', (A(1)-1) / 2 * sim_pixsize .* [-1:0.5:1]);
xlabel(['y / (\lambda_{exc}/n)'], 'FontSize', 14);
set(gca, 'YTick', [1:(A(1)-1)./4:A(1)]);
% set(gca, 'YTickLabel', sim_half_FOV .* [-1:0.5:1]);
set(gca, 'YTickLabel', (A(1)-1) / 2 * sim_pixsize .* [-1:0.5:1]);
ylabel(['z / (\lambda_{exc}/n)'], 'FontSize', 14);
text(0.05 .*A(1), -0.03 .* A(2), ['Overall PSF From ', source_descrip], 'FontSize', 12);
if ~false
    hold on
    plot((A(2) + 1) / 2, (A(1) + 1) / 2, 'o')
    hold on
    plot(interp_peaky - (hsz - 50) + 1, interp_peakz - (hsz - 50) + 1, '*')
end

%
%calc and plot the experimental overall PSF, gamma adjusted:
% figure  %create a new figure window for the plots
subplot(2,5,9)
% actually yz psf
% set(gcf, 'Position', [550 100 600 600]);
% axes_h = axes('Position', [0.13, 0.1, 0.8, 0.8]);
image(256 .* (yz_exp_PSF.^gamma));
colormap hot(256);
axis([1 A(1) 1 A(2)]);
axis square;
set(gca, 'XTick', [1:(A(1)-1)./4:A(1)]);
% set(gca, 'XTickLabel', sim_half_FOV .* [-1:0.5:1]);
set(gca, 'XTickLabel', (A(1)-1) / 2 * sim_pixsize .* [-1:0.5:1]);
xlabel(['y / (\lambda_{exc}/n)'], 'FontSize', 14);
set(gca, 'YTick', [1:(A(1)-1)./4:A(1)]);
% set(gca, 'YTickLabel', sim_half_FOV .* [-1:0.5:1]);
set(gca, 'YTickLabel', (A(1)-1) / 2 * sim_pixsize .* [-1:0.5:1]);
ylabel(['z / (\lambda_{exc}/n)'], 'FontSize', 14);
text(-0.05.*A(1), -0.03.*A(2), ['Overall PSF From ', source_descrip, ', gamma = ', num2str(gamma, '%1.2f')], 'FontSize', 12);
%

%extract the xz plane at the center of the OTF:
midpt = (B(2)+1)./2;
%extract the xz PSF at the current y plane:
yz_exp_OTF = squeeze(OTF3D(:,midpt,:))';
A = size(yz_exp_OTF);
%
%truncate the OTF data to the phyiscally possible range of kmax = +/-4*pi*n/lambda_exc
%   given by the maximum spatial freq possible (defined by two counterpropagating 
%   plane waves):
minkpix = round((1-1./MaxCalcK).*(A(1)-1)./2 + 1);
maxkpix = round((1+1./MaxCalcK).*(A(1)-1)./2 + 1);
plot_yz_exp_OTF = yz_exp_OTF(minkpix:maxkpix, minkpix:maxkpix);
A = size(plot_yz_exp_OTF);
%
%plot the experimental overall xz OTF
% figure  %create a new figure window for the OTF plot
% xruan not include this panel

% subplot(3,3,2)
% % set(gcf, 'Position', [100 100 600 600]);
% % axes_h = axes('Position', [0.13, 0.12, 0.8, 0.8]);
% image(256 .* plot_xz_exp_OTF);
% colormap hot(256);
% axis([1 A(1) 1 A(2)]);
% axis square;
% set(gca, 'XTick', [1:(A(1)-1)./4:A(1)]);
% set(gca, 'XTickLabel', [-1:0.5:1], 'FontSize', 12);
% xlabel(['k_{y} / (4\pin/\lambda_{exc})'], 'FontSize', 14);
% set(gca, 'YTick', [1:(A(2)-1)./4:A(2)]);
% set(gca, 'YTickLabel', [-1:0.5:1], 'FontSize', 12);
% ylabel(['k_{z} / (4\pin/\lambda_{exc})'], 'FontSize', 14);
% text(0.0 .*A(1), -0.03 .* A(2), ['Overall OTF From ', source_descrip], 'FontSize', 12);
%
%plot the gamma adjusted overall xz OTF
% figure  %create a new figure window for the OTF plot
% actually yz OTF
subplot(2,5,5)
% set(gcf, 'Position', [100 100 600 600]);
% axes_h = axes('Position', [0.13, 0.12, 0.8, 0.8]);
image(256 .* (plot_yz_exp_OTF.^gamma));
colormap hot(256);
axis([1 A(1) 1 A(2)]);
axis square;
midpt = (A(1)+1)./2;
line([1 A(2)], [midpt midpt], 'LineWidth', 1, 'Color', [0 0 1]);
line([midpt midpt], [1 A(1)], 'LineWidth', 1, 'Color', [1 0 0]);
FattestColumn = round(midpt + A(1).*NAdet./1.33./4);
line([FattestColumn FattestColumn], [1 A(1)], 'LineWidth', 1, 'Color', [0 1 0]);
set(gca, 'XTick', [1:(A(1)-1)./4:A(1)]);
set(gca, 'XTickLabel', [-1:0.5:1], 'FontSize', 12);
xlabel(['k_{y} / (4\pin/\lambda_{exc})'], 'FontSize', 14);
set(gca, 'YTick', [1:(A(2)-1)./4:A(2)]);
set(gca, 'YTickLabel', [-1:0.5:1], 'FontSize', 12);
ylabel(['k_{z} / (4\pin/\lambda_{exc})'], 'FontSize', 14);
text(-0.05 .*A(1), -0.03 .* A(2), ['Overall OTF From ', source_descrip, ', gamma = ', num2str(gamma, '%1.2f')], 'FontSize', 12);
%
%plot orthogonal linecuts through the lattice overall OTF, as well as  
%    an axial linecut at the fattest part of the xz OTF:
LateralOTFCrossSection = squeeze(plot_yz_exp_OTF(midpt,:));
max_LateralOTFCrossSection = max(LateralOTFCrossSection);
LateralOTFCrossSection = LateralOTFCrossSection ./ max(LateralOTFCrossSection);
yOTF_linecut = LateralOTFCrossSection;

midpt = (A(1)+1)./2;
AxialOTFCrossSection = squeeze(plot_yz_exp_OTF(:,midpt)');
AxialOTFCrossSection = AxialOTFCrossSection ./ max_LateralOTFCrossSection;
zOTF_linecut = AxialOTFCrossSection;
OffsetAxialLinecut = squeeze(plot_yz_exp_OTF(:,FattestColumn)');
OffsetAxialLinecut = OffsetAxialLinecut./max_LateralOTFCrossSection;
zOTF_bowtie_linecut_yz = OffsetAxialLinecut;


% figure  %create a new figure window for the plots
subplot(2,5,10)
% set(gcf, 'Position', [950 100 600 600]);
% axes_h = axes('Position', [0.13, 0.1, 0.8, 0.77]);
D = size(AxialOTFCrossSection);
plot(log10(AxialOTFCrossSection), 'r', 'LineWidth', 2);
hold on
plot(log10(LateralOTFCrossSection), 'b', 'LineWidth', 2);
plot(log10(OffsetAxialLinecut), 'g', 'LineWidth', 2);
axis([1 D(2) -3 0]);
axis square;
grid on;
set(gca, 'XTick', [1:(D(2)-1)./10:D(2)]);
set(gca, 'XTickLabel', [-1:0.2:1]);
xlabel(['k / (4\pi/\lambda)'], 'FontSize', 14);
set(gca, 'YTick', [-3:1:0]);
set(gca, 'YTickLabel', 10.^[-3:1:0]);
ylabel(['OTF Strength'], 'FontSize', 14);
text(-0.1 .*A(2), 0.15, ['Overall OTF linecuts From ', source_descrip], 'FontSize', 12);
text(0.6.*A(2), -0.15, 'OTF along ky', 'Color', [0 0 1], 'FontSize', 14);
text(0.6.*A(2), -0.3, 'OTF along kz', 'Color', [1 0 0], 'FontSize', 14);
text(0.6.*A(2), -0.45, 'Bowtie OTF along kz', 'Color', [0 0.75 0], 'FontSize', 14);
%
%extract the xy plane at the center of the OTF:
midpt = (B(2)+1)./2;
%extract the xz PSF at the current y plane:
xy_exp_OTF = squeeze(OTF3D(:,:,midpt));
A = size(xy_exp_OTF);
%
%truncate the OTF data to the phyiscally possible range of kmax = +/-4*pi*n/lambda_exc
%   given by the maximum spatial freq possible (defined by two counterpropagating 
%   plane waves):
minkpix = (1-1./MaxCalcK).*(A(1)-1)./2 + 1;
maxkpix = (1+1./MaxCalcK).*(A(1)-1)./2 + 1;
plot_xy_exp_OTF = xy_exp_OTF(minkpix:maxkpix, minkpix:maxkpix);
A = size(plot_xy_exp_OTF);
%
%plot the gamma adjusted overall xy OTF
% figure  %create a new figure window for the OTF plot
subplot(2,5,6)
% set(gcf, 'Position', [100 100 600 600]);
% axes_h = axes('Position', [0.13, 0.12, 0.8, 0.8]);
image(256 .* (plot_xy_exp_OTF.^gamma));
colormap hot(256);
axis([1 A(1) 1 A(2)]);
axis square;
midpt = (A(1)+1)./2;
%line([1 A(2)], [midpt midpt], 'LineWidth', 1, 'Color', [0 0 1]);
set(gca, 'XTick', [1:(A(1)-1)./4:A(1)]);
set(gca, 'XTickLabel', [-1:0.5:1], 'FontSize', 12);
xlabel(['k_{x} / (4\pin/\lambda_{exc})'], 'FontSize', 14);
set(gca, 'YTick', [1:(A(2)-1)./4:A(2)]);
set(gca, 'YTickLabel', [-1:0.5:1], 'FontSize', 12);
ylabel(['k_{y} / (4\pin/\lambda_{exc})'], 'FontSize', 14);
text(-0.05 .*A(1), -0.03 .* A(2), ['Overall OTF From ', source_descrip, ', gamma = ', num2str(gamma, '%1.2f')], 'FontSize', 12);


% xruan: add new pannels
% plot xy PSF max plane slice with gamma
% find the brightest one slice
[~, peakInd] = max(PSF3D(:));
[maxyplane, maxxplane, maxzplane] = ind2sub(size(PSF3D), peakInd);
% maxzplane = peakypix + maxplane - 6;
xy_exp_PSF = squeeze(PSF3D(:,:,maxzplane)); 
hsz = round((size(xy_exp_PSF) + 1) / 2);
% xruan: narrow down the range to [-50, 50]
xy_exp_PSF = xy_exp_PSF(hsz(1) - 50 : hsz(1) + 50, hsz(2) - 50 : hsz(2) + 50);
A = size(xy_exp_PSF);

subplot(2, 5, 1)
% set(gcf, 'Position', [550 100 600 600]);
% axes_h = axes('Position', [0.13, 0.1, 0.8, 0.8]);
image(256 .* (xy_exp_PSF));
colormap hot(256);
axis([1 A(1) 1 A(2)]);
axis square;
set(gca, 'XTick', [1:(A(1)-1)./4:A(1)]);
% set(gca, 'XTickLabel', sim_half_FOV .* [-1:0.5:1]);
set(gca, 'XTickLabel', (A(1)-1) / 2 * sim_pixsize .* [-1:0.5:1]);
xlabel(['x / (\lambda_{exc}/n)'], 'FontSize', 14);
set(gca, 'YTick', [1:(A(1)-1)./4:A(1)]);
% set(gca, 'YTickLabel', sim_half_FOV .* [-1:0.5:1]);
set(gca, 'YTickLabel', (A(1)-1) / 2 * sim_pixsize .* [-1:0.5:1]);
ylabel(['y / (\lambda_{exc}/n)'], 'FontSize', 14);
text(-0.05.*A(1), -0.03.*A(2), ['Overall PSF From ', source_descrip], 'FontSize', 12);
if ~false
    hold on
    plot((A(2) + 1) / 2, (A(1) + 1) / 2, 'o')
    hold on
    plot(interp_peakx - (hsz - 50) + 1, interp_peaky - (hsz - 50) + 1, '*')
end


% plot xz PSF max plane slice withpout gamma
xz_exp_PSF = squeeze(PSF3D(maxyplane,:,:))'; 
hsz = round((size(xz_exp_PSF) + 1) / 2);
% xruan: narrow down the range to [-50, 50]
xz_exp_PSF = xz_exp_PSF(hsz(1) - 50 : hsz(1) + 50, hsz(2) - 50 : hsz(2) + 50);
A = size(xz_exp_PSF);
subplot(2, 5, 2)
% set(gcf, 'Position', [550 100 600 600]);
% axes_h = axes('Position', [0.13, 0.1, 0.8, 0.8]);
image(256 .* (xz_exp_PSF));
colormap hot(256);
axis([1 A(1) 1 A(2)]);
axis square;
set(gca, 'XTick', [1:(A(1)-1)./4:A(1)]);
% set(gca, 'XTickLabel', sim_half_FOV .* [-1:0.5:1]);
set(gca, 'XTickLabel', (A(1)-1) / 2 * sim_pixsize .* [-1:0.5:1]);
xlabel(['x / (\lambda_{exc}/n)'], 'FontSize', 14);
set(gca, 'YTick', [1:(A(1)-1)./4:A(1)]);
% set(gca, 'YTickLabel', sim_half_FOV .* [-1:0.5:1]);
set(gca, 'YTickLabel', (A(1)-1) / 2 * sim_pixsize .* [-1:0.5:1]);
ylabel(['z / (\lambda_{exc}/n)'], 'FontSize', 14);
text(-0.05.*A(1), -0.03.*A(2), ['Overall PSF From ', source_descrip,], 'FontSize', 12);
if ~false
    hold on
    plot((A(2) + 1) / 2, (A(1) + 1) / 2, 'o')
    hold on
    plot(interp_peakx - (hsz - 50) + 1, interp_peakx - (hsz - 50) + 1, '*')
end


%plot the gamma adjusted overall xz PSF
% figure  %create a new figure window for the PSF plot
subplot(2,5,7)
% actually yz psf
% set(gcf, 'Position', [550 100 600 600]);
% axes_h = axes('Position', [0.13, 0.1, 0.8, 0.8]);
image(256 .* (xz_exp_PSF .^ gamma));
colormap hot(256);
axis([1 A(1) 1 A(2)]);
axis square;
set(gca, 'XTick', [1:(A(1)-1)./4:A(1)]);
% set(gca, 'XTickLabel', sim_half_FOV .* [-1:0.5:1]);
set(gca, 'XTickLabel', (A(1)-1) / 2 * sim_pixsize .* [-1:0.5:1]);
xlabel(['x / (\lambda_{exc}/n)'], 'FontSize', 14);
set(gca, 'YTick', [1:(A(1)-1)./4:A(1)]);
% set(gca, 'YTickLabel', sim_half_FOV .* [-1:0.5:1]);
set(gca, 'YTickLabel', (A(1)-1) / 2 * sim_pixsize .* [-1:0.5:1]);
ylabel(['z / (\lambda_{exc}/n)'], 'FontSize', 14);
text(-0.05.*A(1), -0.03.*A(2), ['Overall PSF From ', source_descrip, ', gamma = ', num2str(gamma, '%1.2f')], 'FontSize', 12);


% plot xz OTF max plane slice with gamma
%extract the xz plane at the center of the OTF:
midpt = (B(1)+1)./2;
%extract the xz PSF at the current y plane:
xz_exp_OTF = squeeze(OTF3D(midpt,:,:))';
A = size(xz_exp_OTF);
%
%truncate the OTF data to the phyiscally possible range of kmax = +/-4*pi*n/lambda_exc
%   given by the maximum spatial freq possible (defined by two counterpropagating 
%   plane waves):
minkpix = round((1-1./MaxCalcK).*(A(1)-1)./2 + 1);
maxkpix = round((1+1./MaxCalcK).*(A(1)-1)./2 + 1);
plot_xz_exp_OTF = xz_exp_OTF(minkpix:maxkpix, minkpix:maxkpix);
A = size(plot_xz_exp_OTF);


%plot the gamma adjusted overall xz OTF
% figure  %create a new figure window for the OTF plot
% actually yz OTF
subplot(2,5,3)
% set(gcf, 'Position', [100 100 600 600]);
% axes_h = axes('Position', [0.13, 0.12, 0.8, 0.8]);
image(256 .* (plot_xz_exp_OTF.^gamma));
colormap hot(256);
axis([1 A(1) 1 A(2)]);
axis square;
midpt = (A(1)+1)./2;
line([1 A(2)], [midpt midpt], 'LineWidth', 1, 'Color', [0 0 1]);
line([midpt midpt], [1 A(1)], 'LineWidth', 1, 'Color', [1 0 0]);
FattestColumn = round(midpt + A(1).*NAdet./1.33./4);
line([FattestColumn FattestColumn], [1 A(1)], 'LineWidth', 1, 'Color', [0 1 0]);
set(gca, 'XTick', [1:(A(1)-1)./4:A(1)]);
set(gca, 'XTickLabel', [-1:0.5:1], 'FontSize', 12);
xlabel(['k_{x} / (4\pin/\lambda_{exc})'], 'FontSize', 14);
set(gca, 'YTick', [1:(A(2)-1)./4:A(2)]);
set(gca, 'YTickLabel', [-1:0.5:1], 'FontSize', 12);
ylabel(['k_{z} / (4\pin/\lambda_{exc})'], 'FontSize', 14);
text(-0.05 .*A(1), -0.03 .* A(2), ['Overall OTF From ', source_descrip, ', gamma = ', num2str(gamma, '%1.2f')], 'FontSize', 12);
%
%plot orthogonal linecuts through the lattice overall OTF, as well as  
%    an axial linecut at the fattest part of the xz OTF:
LateralOTFCrossSection = squeeze(plot_xz_exp_OTF(midpt,:));
max_LateralOTFCrossSection = max(LateralOTFCrossSection);
LateralOTFCrossSection = LateralOTFCrossSection ./ max(LateralOTFCrossSection);
xOTF_linecut = LateralOTFCrossSection;

midpt = (A(1)+1)./2;
AxialOTFCrossSection = squeeze(plot_xz_exp_OTF(:,midpt)');
AxialOTFCrossSection = AxialOTFCrossSection ./ max_LateralOTFCrossSection;
zOTF_linecut = AxialOTFCrossSection;
OffsetAxialLinecut = squeeze(plot_xz_exp_OTF(:,FattestColumn)');
OffsetAxialLinecut = OffsetAxialLinecut./max_LateralOTFCrossSection;
zOTF_bowtie_linecut = OffsetAxialLinecut;


% figure  %create a new figure window for the plots
subplot(2,5,8)
% set(gcf, 'Position', [950 100 600 600]);
% axes_h = axes('Position', [0.13, 0.1, 0.8, 0.77]);
D = size(AxialOTFCrossSection);
plot(log10(AxialOTFCrossSection), 'r', 'LineWidth', 2);
hold on
plot(log10(LateralOTFCrossSection), 'b', 'LineWidth', 2);
plot(log10(OffsetAxialLinecut), 'g', 'LineWidth', 2);
axis([1 D(2) -3 0]);
axis square;
grid on;
set(gca, 'XTick', [1:(D(2)-1)./10:D(2)]);
set(gca, 'XTickLabel', [-1:0.2:1]);
xlabel(['k / (4\pi/\lambda)'], 'FontSize', 14);
set(gca, 'YTick', [-3:1:0]);
set(gca, 'YTickLabel', 10.^[-3:1:0]);
ylabel(['OTF Strength'], 'FontSize', 14);
text(-0.1 .*A(2), 0.15, ['Overall OTF linecuts From ', source_descrip], 'FontSize', 12);
text(0.6.*A(2), -0.15, 'OTF along kx', 'Color', [0 0 1], 'FontSize', 14);
text(0.6.*A(2), -0.3, 'OTF along kz', 'Color', [1 0 0], 'FontSize', 14);
text(0.6.*A(2), -0.45, 'Bowtie OTF along kz', 'Color', [0 0.75 0], 'FontSize', 14);


% add text for the annotation of the coordinates
text(-A(1) * 0.25, -3.75, 'Image y = x galvo dither direction & y-stage scan direction', 'fontSize', 14)
text(-A(1) * 0.25, -4.05, 'In deskewed data image x = LLS(y) propogation direction', 'fontSize', 14);
text(-A(1) * 0.1, -4.3, sprintf('Peak coordinate (xyz): (%d, %d, %d)', peakx, peaky, peakz), 'fontSize', 14)


end





