im = rand(100,100,100,'single');
parallelWriteTiff('test.tif',im);

%%
imT = parallelReadTiff('test.tif');

%%
imSize = getImageSizeMex('test.tif');