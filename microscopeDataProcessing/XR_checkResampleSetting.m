function [resample, zAniso] = XR_checkResampleSetting(resample_type, resample, ObjectiveScan, SkewAngle, xyPixelSize, dz);
% check resample settings
%
% Author: Xiongtao Ruan (01/13/2021)


if ObjectiveScan
    zAniso = dz/xyPixelSize;    
else
    zAniso = sind(SkewAngle)*dz/xyPixelSize;
end

switch resample_type
    case 'given'  %% given resample 
        if isempty(resample)
            error('For resample Type "given", the parameter resample must not be empty!');
        else
            if numel(resample) == 1
                resample = ones(1, 3) * resample;
            elseif numel(resample) == 2
                resample = [ones(1, 2) * resample(1), resample(2)];
            end
        end
    case 'isotropic'  %% isotropic
        resample = [1,1,1];
    case 'xy_isotropic' %% x, y isotropic and z: r_z / r_x
        theta = SkewAngle * pi/180;
        
        zf = sqrt((sin(theta) ^ 2 + zAniso ^ 2 * cos(theta) ^ 2) / (cos(theta) ^ 2 + zAniso ^ 2 * sin(theta) ^ 2));
        resample = [1, 1, zf];
end

end


