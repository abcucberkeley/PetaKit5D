function RLdecon_CBC(input_tiff, psf, background, nIter, dz_psf, dz_data, varargin)
%UNTITLED2 Summary of this function goes here
%   input_tiff: input TIFF file name
%   psf: psf array in 'double'
%   background: backgroud to subtract
%   nIter: number of iterations

%Putty command line to complie: mcc -m -R -nojvm -v RLdecon.m

if ischar(dz_psf)
    dz_psf=str2double(dz_psf);
end
if ischar(dz_data)
    dz_data=str2double(dz_data);
end

if ischar(psf)
    [a,b,suffix]=fileparts(psf);
    if strcmp(suffix, '.mat')
        load(psf, 'psf');
    else
        if strcmp(suffix, '.tif')
            psf=psf_gen(psf, dz_psf, dz_data, 48);
        end
    end     
end
if ischar(background)
    background=str2double(background);
end
if ischar(nIter)
    nIter=str2num(nIter);
end

nTapering = 0;

for k = 1 : length(varargin);
    switch k
        case 1
            % number of pixel for x-y tapering
            if ischar(varargin{k})
                nTapering = str2num(varargin{k});
            else
                nTapering = varargin{k};
            end
        otherwise
            disp('Unknown varargin index')
    end
end

rawtiff=Tiff(input_tiff,'r');
z=1;
while ~rawtiff.lastDirectory
    rawdata(:,:,z)=single(rawtiff.read());
    z=z+1;
    rawtiff.nextDirectory()
end
rawdata(:,:,z)=single(rawtiff.read());
nz = z;

rawtiff.close()

rawdata = rawdata - background;
tic
rawdata(rawdata<0.0) = 0.0;

if nTapering > 0
    taperKernel = fspecial('gaussian', nTapering+30, nTapering);
    for iz=1:nz
        rawdata(:,:,iz) = edgetaper(rawdata(:,:,iz), taperKernel);
    end
end

deconvolved = deconvlucy(rawdata, psf, nIter) * numel(rawdata);
toc
% construct output file name
output_tiff = strrep(input_tiff,'.tif','_decon.tif');
inds=findstr(output_tiff, '/');
if any(inds)
    l0=inds(length(inds));
    l1=length(output_tiff);
    output_tiff = output_tiff(l0+1:l1);
end

write3Dtiff(deconvolved, output_tiff);


end

