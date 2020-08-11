function timeEqualize( tiffs )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
for i=1:length(tiffs)
    img=loadtiff(tiffs(i).name);
    meanIntensity = mean(img(:));
    mean_above_background = mean(img(img>meanIntensity));
    if i==1
        mean_above_background0=mean_above_background;
    else
        ratio1 = mean_above_background0 / mean_above_background
        img = img * ratio1;
    end
   write3Dtiff(img, strcat('../rescaled/',regexprep(tiffs(i).name, '_[0-9]*nm_[0-9]*msec_decon','_decon')))
end
