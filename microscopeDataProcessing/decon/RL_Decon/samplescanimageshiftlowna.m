
tagstruct.ImageLength=416;
tagstruct.ImageWidth=416;
tagstruct.Photometric=Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample=32;
tagstruct.SampleFormat=3;
tagstruct.SamplesPerPixel=1;
tagstruct.PlanarConfiguration =1;
tagstruct.Software='MATLAB';
% Using YResolution tag to indicate z step
tagstruct.YResolution=0.25*cos(45.6*pi/180); %0.2*cos(58.5*pi/180);
nz=301;%zplane
nphases=5;%phase
xdata=[71 486];% number higher, shift toward left
ydata=[1 416];

[FileName,PathName,Filterindex] = uigetfile('*.tif');
infile=[PathName FileName];
outfile=[PathName 'translated' FileName];
tifin=Tiff(infile, 'r');
tifout=Tiff(outfile, 'w');
a=102.4633; %mean value 
b=0.3053; %stad value 
c=a+b*randn(416,416);

for z=0:nz-1
%   T=maketform('affine', [1 0 0; 0 1 0; 0.8261*(nz-1-z) 0 1]); %[1 0 0; 0 1 0; 1.5982*z 0 1]);[1 0 0; 0 1 0; 0.74297*z 0 1] zstep*cos(theta)./pixel size
T=maketform('affine', [1 0 0; 0 1 0; 1.4029*(nz-1-z) 0 1]); 
% T=maketform('affine', [1 0 0; 0 1 0; 0.9186*z 0 1]);
  for ph=0:nphases-1
    rawsec=z*nphases+ph+1;
    tifin.setDirectory(rawsec);
    img=single(tifin.read());
    imgT = imtransform(img, T, 'bicubic','xdata', xdata, 'ydata', ydata, 'fill', 102.5);
    C=imgT;
    C(find(imgT==102.5))=c(find(imgT==102.5));
    imgT=C;
    tifout.setTag(tagstruct);
    tifout.write(imgT);
    if rawsec<nz*nphases
      tifout.writeDirectory();
    end
  end
end

tifout.close();
tifin.close();