
%load raw images
background=60;
tiff_files = dir('*.tif');
nfiles = length(tiff_files);

for n=1:nfiles
    
    if any(regexp(tiff_files(n).name, '[0-9]ch_[0-9]{4}.tif'))
        rawtiff=Tiff(tiff_files(n).name,'r');
        z=1;
        while ~rawtiff.lastDirectory
            rawdata(:,:,z)=single(rawtiff.read());
            z=z+1;
            rawtiff.nextDirectory()
        end
        rawdata(:,:,z)=single(rawtiff.read());
        
        deconvolved = deconvlucy(rawdata-background, psf_resampled, 10);
        write3Dtiff(deconvolved, strrep(tiff_files(n).name,'.tif','_decon10Itr_bg60_nonSquare.tif'));
        rawtiff.close()
    end
end