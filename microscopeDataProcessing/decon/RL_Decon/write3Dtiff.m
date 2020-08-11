function write3Dtiff(array, filename)

outtiff=Tiff(filename, 'w');

dims = size(array);
if length(dims)==2
    dims(3)=1;
end

tagstruct.ImageLength = dims(1);
tagstruct.ImageWidth = dims(2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
if isa(array, 'single')
    tagstruct.BitsPerSample = 32;
    tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
else if isa(array, 'uint16')
    tagstruct.BitsPerSample = 16;
    tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
    end
end
tagstruct.SamplesPerPixel = 1;
%tagstruct.RowsPerStrip = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;

%outtiff.setTag(tagstruct)

for z=1:dims(3)
    %outtiff.currentDirectory
    outtiff.setTag(tagstruct)
%     outtiff.setTag('ImageWidth', dims(1))
%     outtiff.setTag('ImageLength', dims(2))
%     outtiff.setTag('Photometric',rawtiff.getTag('Photometric'))
%     outtiff.setTag('PlanarConfiguration',rawtiff.getTag('PlanarConfiguration'))
%     outtiff.setTag('BitsPerSample',32)
%     outtiff.setTag('SampleFormat', Tiff.SampleFormat.IEEEFP);
    outtiff.write(array(:,:,z))
    outtiff.writeDirectory()
end
outtiff.close()
