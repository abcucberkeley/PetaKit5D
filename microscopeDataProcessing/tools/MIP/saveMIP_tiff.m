function [] = saveMIP_tiff(frameFullname, MIPFullname, varargin)
% get the mip for a given image with image path
% 
% Author: Xiongtao Ruan (03/30/2021)
% 
% xruan (05/05/2022): add support for user defined axis MIPs

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('frameFullpath'); 
ip.addRequired('mipFullpath'); 
ip.addParameter('dtype', 'uint16', @ischar); % suffix for the folder
ip.addParameter('axis', [0, 0, 1], @isnumeric); % suffix for the folder

ip.parse(frameFullname, MIPFullname, varargin{:});

if exist(MIPFullname, 'file')
    fprintf('The mip file %s already exists, skip it!\n', MIPFullname)
    return;
end

dtype = ip.Results.dtype;
axis = ip.Results.axis;

im = readtiff(frameFullname);

for i = 1 : 3
    if axis(i) == 0
        continue;
    end
    switch i
        case 1
            MIP = project3DImageto2D(im, 'mip_xz');
            MIPFullname = sprintf('%s_MIP_y.tif', MIPFullname(1 : end - 10));
        case 2
            MIP = project3DImageto2D(im, 'mip_yz');     
            MIPFullname = sprintf('%s_MIP_x.tif', MIPFullname(1 : end - 10));
        case 3
            MIP = project3DImageto2D(im, 'mip_xy');     
            MIPFullname = sprintf('%s_MIP_z.tif', MIPFullname(1 : end - 10));
    end
    MIP = cast(MIP, dtype);
    writetiff(MIP, MIPFullname);
end

% if ~exist(MIPFullname, 'file')
%     imwrite(MIP, MIPFullname);
% end

end

