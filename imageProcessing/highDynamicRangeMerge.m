function [combinedMystery, logMystery, logImage] = highDynamicRangeMerge(image_low, image_high, varargin)
% this function takes two images of the same field taken at different
% exposure levels, removes the saturated region from the higher exposed
% image, and replaces it with a rescaled version of the same region 
% from the lower exposed image. This should result in a smooth transition
% that forms one high dynamic range image. Currently, this function also
% flattens the image to allow viewing of all features easily despite
% initial brightness and then converts it to a 16 bit image and writes it
% to a file
%- Jessica Tytell October 25, 2011

% parse input
ip = inputParser;
ip.addRequired('image_low', @isnumeric); 
ip.addRequired('image_high', @isnumeric);
ip.addOptional('mysteryOffsetFactor', 1.25, @isnumeric);
ip.parse(image_low, image_high, varargin{:});
mysteryOffsetFactor=ip.Results.mysteryOffsetFactor;

%show that files are uploaded _ for troubleshooting only
% figure('Name', 'Dim'), imshow(image_low,[]);
% figure('Name', 'bright'), imshow(image_high,[]);

% Make mask of saturated pixels and display. Saturated pixels =1
satMask = image_high == 4095;
% figure('Name', 'satMask'), imshow(satMask,[]);

%make masked image of bright that gets rid of saturated pixels
brightMasked = image_high.*~satMask;
%figure('Name', 'BrightMasked'), imshow(brightMasked,[]);

%make image of only regions from dim image that are saturated in bright
%image
dimMasked = image_low.*satMask;
% figure('Name', 'dimMasked'), imshow(dimMasked,[]);

% Make a mask for background region (coarse, for merging purpose only)
segment_level = thresholdRosin(image_high);
backgroundMask = image_high < segment_level;

% A combined mask for both the saturated area and the background area.
bg_st_Mask = or(satMask,backgroundMask);

%find polynomial fit (linear) that says relationship between two images
% figure, plot(image_high(~satMask(:)), image_low(~satMask(:)));
lineValues = polyfit(image_high(~satMask(:)), image_low(~satMask(:)), 1);
coeff=lineValues(1);
offset = lineValues(2);
% hold on; plot(image_high(~satMask(:)), coeff*image_high(~satMask(:))+offset,'-r')

% find a linear relationship between two images, linear, not affine.
High_data = image_high(~bg_st_Mask(:));
Low_data = image_low(~bg_st_Mask(:));
bright2dim_scaler = (High_data'*Low_data)/(High_data'*High_data);

%NOTE various methods for combining images are listed below. Currently some
%are used and some are not. None are being deleted until I figure out how
%to make a smoother image

%combine images using correction
combined = (brightMasked.*coeff + offset) + dimMasked;
% figure('Name', 'Combined'), imshow(combined,[]);
% maxNum = max(combined(:));
% scaledCombined = scaleContrast(combined, [0 maxNum], [0 255]);
% imwrite(uint8(scaledCombined), 'scaledCombined.tif', 'tiff');

%flatten combined image by taking the log of the image
combinedTwentyFivePercent = (brightMasked.*0.25) + dimMasked;
logImage = log(combinedTwentyFivePercent);
% figure('Name', 'logimage'), imshow(logImage,[]);

% combine images using linear correction
combined_linear = (brightMasked.*bright2dim_scaler) + dimMasked;
% flatten combined image by taking the log of the image
% here, this temporary overwrites the output logImage, will check in further
% for how correct it is
logImage = log(combined_linear);

% scale image using mystery factor
combinedMystery= (brightMasked.*coeff.*mysteryOffsetFactor + offset) + dimMasked;
% figure('Name', 'CombinedMystery'), imshow(combinedMystery,[]);
% imwrite(scaleContrast(combinedMystery, [], [0 65535]));

logMystery = log(combinedMystery);
%figure('Name', 'log Mystery'), imshow(logMystery,[]);
% scaledMystery = (logMystery/10) * 2^15;
% intlogMystery = uint16(scaledMystery);
% imwrite(intlogMystery, 'test2.tif', 'tiff');


