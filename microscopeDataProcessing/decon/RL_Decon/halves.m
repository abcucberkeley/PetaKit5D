function halves( tiffs, centline, midline )
%HALVES make left-right or top-bottom halves (in total 4) of a series of
%TIFF files
%   "centline" defines the left-right divide line
%   "midline" defines the top-bottom divide line

[s, mess, messid] = mkdir('righthalf');
[s, mess, messid] = mkdir('lefthalf');
[s, mess, messid] = mkdir('bothalf');
[s, mess, messid] = mkdir('tophalf');

if s
    for i=1:length(tiffs)
        img=loadtiff(tiffs(i).name);
        [nx, ny, nz] = size(img);
        write3Dtiff(img(:, centline:ny,:), strcat('righthalf/', strrep(tiffs(i).name, '_decon.tif', '_decon_righthalf.tif')))
        write3Dtiff(img(:, 1:centline,:), strcat('lefthalf/', strrep(tiffs(i).name, '_decon.tif', '_decon_lefthalf.tif')))
        write3Dtiff(img(1:midline, :,:), strcat('tophalf/', strrep(tiffs(i).name, '_decon.tif', '_decon_tophalf.tif')))
        write3Dtiff(img(midline:nx, :,:), strcat('bothalf/', strrep(tiffs(i).name, '_decon.tif', '_decon_bothalf.tif')))
    end
else
    disp(mess)
end

end

