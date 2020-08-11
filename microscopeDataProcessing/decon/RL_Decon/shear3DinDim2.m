function [ output ] = shear3DinDim2(inArray, angle, bReverse, dz, xypixelsize, fillVal, varargin)
%shear3DinDim2 -- deskew inArray about the 1st dimension (y in image terms)
%x in sheet-scan terms)
%   inArray -- input 3D array
%   angle -- skewing angle
%   bReverse -- is skewing in reverse direction ? (like in Bi-Chang's data)
%   dz -- sample-scan step size (in microns)
%   xypixelsize -- in microns
%   the right)
%   fillVal -- use this value to fill the new void

trans = 0;
nphases = 1;
for k = 1 : length(varargin);
    switch k
        case 1
            % user-specified final output's width
            if ~isempty(varargin{k})
                outputWidth = varargin{k};
            end
        case 2
            % extra left-right translation (positive number translates the result toward
            if ~isempty(varargin{k})
                trans = varargin{k};
            end
        case 3
            % number of phases in raw data; default to 1 in non-SIM mode
            if ~isempty(varargin{k})
                nphases = varargin{k};
            end
        otherwise
            disp('Unknown varargin index')
    end
end
center = (size(inArray)./[1 1 nphases]+1)/2;

trans1 = [1 0 0 0
    0 1 0 0
    0 0 1 0
    -center 1];

rot = [1 0 0 0
    0 1 0 0
    0 cos(angle*pi/180)*dz/xypixelsize 1 0
    0 0 0 1];

if bReverse
    rot(3, 2) = -rot(3, 2);
end

if exist('outputWidth', 'var') && outputWidth > 0
    widenBy = outputWidth - size(inArray, 2);
else
    widenBy = ceil(size(inArray, 3)/nphases *dz * (cos(angle*pi/180)) / xypixelsize/2);
end

trans2 = [1 0 0 0
    0 1 0 0
    0 0 1 0
    center+[0 widenBy/2+trans 0] 1];

T = trans1*rot*trans2;

tform = maketform('affine', T);
R = makeresampler('linear', 'fill');

if nphases > 1
    output = zeros(size(inArray)+[0 widenBy 0], 'single');
    for p = 1:nphases
         output(:,:,p:nphases:size(inArray,3)) = ... 
            tformarray(inArray(:,:,p:nphases:size(inArray,3)), ...
            tform, R, [1 2 3], [1 2 3], (size(inArray)+[0 widenBy 0])./[1 1 nphases] , [], fillVal);
    end
else
    output = tformarray(inArray, tform, R, [1 2 3], [1 2 3], ...
        size(inArray)+[0 widenBy 0] , [], fillVal);
end

end
