im = rand(100,100,100,'single');
parallelWriteZarr('test.zarr',im);

%%
imT = parallelReadZarr('test.zarr');

%%
createZarrFile('testT.zarr');