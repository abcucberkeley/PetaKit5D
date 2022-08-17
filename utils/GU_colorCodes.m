function [ceR, ceB, ceG, ceP, ceO, cfR, cfB, cfG, cfP, cfO,cfK,ceK] = GU_colorCodes

% Gokul Upadhyayula, Jan 2015

ceR = hsv2rgb([1.00 0.85 0.85]);
ceB = hsv2rgb([0.55 0.85 0.85]);
ceG = hsv2rgb([0.30 0.85 0.85]);
ceP = hsv2rgb([0.80 0.85 .85]); %purple
ceO = hsv2rgb([0.10 0.85 .85]); % orange
cfR = hsv2rgb([1.00 0.3 1]);
cfB = hsv2rgb([0.55 0.3 1]);
cfG = hsv2rgb([0.30 0.3 1.0]);
cfP = hsv2rgb([0.80 0.3 1.0]); %purple
cfO = hsv2rgb([0.10 0.8 1]); % orange

cfK = hsv2rgb([0.5 0.1 0.5]); 
ceK = hsv2rgb([0 0 0]); 