function cropSeries (foldername, zSubRanges)
%cropSeries: z crops for all TIFFs under foldername
%   zSubRanges: organize lower/upper limit pairs into row vectors; for
%   example, [20 30; 40 70; 80 100] -- doing three crops between (20, 30),
%   (40, 70) and (80, 100)
%   Results wil be saved in subfolders 'zcropped xx-yy' named 

oldfolder = cd(foldername);
tiffs=dir('*.tif');

% make subfolders named after subranges
nCrops = size(zSubRanges, 1);
for i=1:nCrops
    subdirNames{i} = sprintf('zcropped %d-%d', zSubRanges(i, 1), zSubRanges(i, 2));
    mkdir(subdirNames{i})
end

for i=1:length(tiffs)
    img=loadtiff(tiffs(i).name);
    for j=1:nCrops
        if zSubRanges(j, 1) < zSubRanges(j, 2)
            write3Dtiff(img(:, :, zSubRanges(j, 1) : zSubRanges(j, 2)), ...
                strcat(subdirNames{j},'/',tiffs(i).name))
        else
            write3Dtiff(img(:, :, zSubRanges(j, 2) : zSubRanges(j, 1)), ...
                strcat(subdirNames{j},'/',tiffs(i).name))
        end
    end
end

cd(oldfolder)