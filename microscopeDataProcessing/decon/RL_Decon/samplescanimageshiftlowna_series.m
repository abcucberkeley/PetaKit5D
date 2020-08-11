clear all
tagstruct.ImageLength=672
tagstruct.ImageWidth=672;
tagstruct.Photometric=Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample=32;
tagstruct.SampleFormat=3;
tagstruct.SamplesPerPixel=1;
tagstruct.PlanarConfiguration =1;
tagstruct.Software='MATLAB';
% Using YResolution tag to indicate z step
tagstruct.YResolution=0.4*cos(32.6*pi/180); %0.2*cos(58.5*pi/180); not necessary parameter
nz=251;%zplane
nphases=5;%phase
xdata=[141 812]; % number higher, shift toward left 
ydata=[1 672];

[FileName,PathName,Filterindex] = uigetfile('*.tif',  'MultiSelect', 'on');
N=size(FileName,2);
for n=1:N
    
infile=[PathName char(FileName(n))];
% infile=strcat(PathName, FileName(1));
outfile=[PathName 'translated' char(FileName(n))];
tifin=Tiff(infile, 'r');
tifout=Tiff(outfile, 'w');

a=110.4633; %mean value 
b=0.3053; %stad value 
c=a+b*randn(672,672);

for z=0:nz-1
  T=maketform('affine', [1 0 0; 0 1 0; 2.4383*(nz-1-z) 0 1]); %[1 0 0; 0 1 0; 1.5982*z 0 1]);[1 0 0; 0 1 0; 0.74297*z 0 1] zstep*cos(theta)./pixel size
% T=maketform('affine', [1 0 0; 0 1 0; 2.3065*(nz-1-z) 0 1]); 
% T=maketform('affine', [1 0 0; 0 1 0; 2.4383*z 0 1]);
  for ph=0:nphases-1
    rawsec=z*nphases+ph+1;
    tifin.setDirectory(rawsec);
    img=single(tifin.read());
    imgT = imtransform(img, T, 'bicubic','xdata', xdata, 'ydata', ydata, 'fill', 110);
    C=imgT;
    C(find(imgT==110))=c(find(imgT==110));
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
end 
