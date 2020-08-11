function resampleTiffs(filestem,scale_factor,threshold)
tiffs=dir([filestem,'\*.tif']);
mkdir(filestem,'downsampled_data');
for i=1:15
eval(['data=loadtiff(''',filestem,'\',tiffs(i).name,''');']);
temp=sort(data(:));
data2=data./temp(round(length(temp)*threshold));
data3=image_resize(data2,round(size(data2)/scale_factor));
eval(['saveastiff(uint16(data3.*2^16),''',filestem,'\downsampled_data\',tiffs(i).name(1:end-4),'_rescaled.tif'');']);
end