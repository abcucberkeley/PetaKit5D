function [sigmaXY, sigmaZ, XYZ, pstruct, mask] = GU_estimateSigma3D(rt,FileName)
% estimates the position of the bead, and fits the sigma as well
% filtered based on the highest intensity

% Gokul Upadhyayula, 2015
% rt = '/Users/Gokul/Desktop/2015-April-13_p505_p6_nup133counting/488totalPSF.tif';
%  rt = '/Users/Gokul/Desktop/2015-April-13_p505_p6_nup133counting/560totalPSF.tif';
% 


if ~isempty(FileName)
    im = [rt FileName];
else
    im = rt;
end

opts.Alpha = 0.05;
opts.CellMask = [];
opts.RemoveRedundant = true;
opts.WindowSize = [];
opts.Mode = 'xyzAsrc';
for i = 1:3
    if i == 1
        sigma = [1.2, 2.5]; % estimate
    else
        sigma = [sigmaXY, sigmaZ];
    end
    
    frame = double(readtiff(im)); %#ok<PFBNS>
    frame(frame==0) = NaN; % to avoid border effects

%         [ny,nx,nz] = size(frame);

    % xruan 12/31/2019
    % [pstruct, mask] = pointSourceDetection3D(frame, sigma, 'Mode', opts.Mode, 'Alpha', opts.Alpha,...
    %     'Mask', opts.CellMask, 'RemoveRedundant', opts.RemoveRedundant,...
    %     'RefineMaskLoG', false, 'WindowSize', opts.WindowSize); %#ok<PFBNS>
    [pstruct, mask] = XR_pointSourceDetection3D(frame, sigma, 'Mode', opts.Mode, 'Alpha', opts.Alpha,...
        'Mask', opts.CellMask, 'RemoveRedundant', opts.RemoveRedundant,...
        'RefineMaskLoG', false, 'WindowSize', opts.WindowSize, ...
        'FitGaussianMethod', 'original', 'BackgroundCorrection', false);
    
    idx = find(pstruct.A==max([pstruct.A]));
    sigmaXY = pstruct.s(1,idx);
    sigmaZ = pstruct.s(2,idx);
    XYZ.x = pstruct.x(idx);
    XYZ.y = pstruct.y(idx);
    XYZ.z = pstruct.z(idx);
    % cd(rt);
    % mkdir(['SigmaMask']);
    % writetiff(uint8(255*mask), [rt 'SigmaMask' filesep FileName '_mask.tif']);
    
end