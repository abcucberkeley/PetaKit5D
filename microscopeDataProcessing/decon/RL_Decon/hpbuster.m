function [ output ] = hpbuster(inArr, background, harshness, border_width)
%UNTITLED detect hot pixels and replace them with the mean of its
%neighbours
%   Detailed explanation goes here

[nx, ny, nz] = size(inArr);
output = inArr;

factor = 4 + harshness;
numPixels = (2*border_width+1) * (2*border_width+1) - 1;

%     function result=hotpixelReplace(a)
% 
%         dim = int32(size(a));
%         
%         center = idivide(dim, 2)+1;
%         centerval = a(center(1), center(2));
%         if centerval > background
%             stddev = sqrt(double(centerval)-background);
%             diff = (centerval - a) > factor*stddev;
%             how_many_different = sum(diff(:));
%             % If all neighbouring pixel are significantly different,
%             % replace it with the mean
%             if how_many_different == numPixels
%                 result = (sum(a(:)) - centerval) / numPixels;
%             else
%                 result = centerval;
%             end
%         else
%             result = centerval;
%         end
%     end


for z=1:nz
%     % nlfilter: General sliding-neighborhood operations
%     output(:,:,z) = nlfilter(inArr(:,:,z), ... 
%         [2*border_width+1 2*border_width+1], ...
%         @hotpixelReplace);

    for y=border_width+1:ny-border_width
        for x=border_width+1:nx-border_width
            if inArr(x, y, z)>background
                stddev = sqrt(single(inArr(x, y, z))-background);
            else
                continue
            end
            diff = (inArr(x, y, z) - inArr(x-border_width:x+border_width, ...
                y-border_width:y+border_width, z)) > factor*stddev;
            how_many_different = sum(diff(:));
%             If all neighbouring pixel are significantly different,
%             replace it with the mean
            if how_many_different == numPixels
                output(x, y, z) = (sum(sum(inArr(x-border_width:x+border_width, ...
                    y-border_width:y+border_width, z))) - inArr(x,y,z)) / numPixels;
            end
        end
    end
end

end


