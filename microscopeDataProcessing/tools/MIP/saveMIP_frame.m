function [] = saveMIP_frame(im, MIPFullname, varargin)
% get the mip for a given image array
% 
% Author: Xiongtao Ruan (11/10/2022)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('im'); 
ip.addRequired('MIPFullname'); 
ip.addParameter('dtype', 'uint16', @ischar); % suffix for the folder
ip.addParameter('axis', [0, 0, 1], @isnumeric); % suffix for the folder

ip.parse(im, MIPFullname, varargin{:});

dtype = ip.Results.dtype;
axis = ip.Results.axis;

if exist(MIPFullname, 'file')
    fprintf('The mip file %s already exists, skip it!\n', MIPFullname)
    return;
end

if all(axis == 0)
    return;
end

MIPPath = fileparts(MIPFullname);
if ~exist(MIPPath, 'dir')
    mkdir(MIPPath);
    if ~ispc                
        fileattrib(MIPPath, '+w', 'g');
    end
end

im = cast(im, dtype);
axislabel = 'yxz';
uuid = get_uuid();
for d=1:3
    if axis(d) ~= 0
        MIPTmpname = sprintf('%s_MIP_%s_%s.tif', MIPFullname(1 : end - 10), axislabel(d), uuid);
        MIPFn = sprintf('%s_MIP_%s.tif', MIPFullname(1 : end - 10), axislabel(d));
        writetiff(max(im, [], d), MIPTmpname);
        movefile(MIPTmpname, MIPFn);                
    end
end

end
