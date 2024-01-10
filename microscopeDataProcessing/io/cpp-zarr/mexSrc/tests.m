im = rand(100,100,100,'single');
parallelWriteZarr('test.zarr',im,1);

%%
imT = parallelReadZarr('test.zarr');

%%
createZarrFile('testT.zarr');