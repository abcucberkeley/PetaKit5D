function [A] = indexing4d(A, R, bbox, rbbox)
% assin cropped region of R with rbbox (if provided) to A according to bbox
%
% Author: Xiongtao Ruan (04/18/2024)


if nargin < 4
    try
        % indexing4d_crop_mex(A, bbox, R);
        indexing4d_mex(A, bbox, R);
    catch ME
        disp(ME);
        disp(ME.stack);
        A(bbox(1) : bbox(5), bbox(2) : bbox(6), bbox(3) : bbox(7), bbox(4) : bbox(8)) = R;
    end
else
    try
        indexing4d_crop_mex(A, bbox, R, rbbox);
    catch ME
        disp(ME);
        disp(ME.stack);
        A(bbox(1) : bbox(5), bbox(2) : bbox(6), bbox(3) : bbox(7), bbox(4) : bbox(8)) ...
            = R(rbbox(1) : rbbox(5), rbbox(2) : rbbox(6), rbbox(3) : rbbox(7), rbbox(4) : rbbox(8));
    end 
end

end
