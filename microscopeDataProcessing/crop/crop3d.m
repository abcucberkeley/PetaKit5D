function [A] = crop3d(A, bbox)
% crop A according to bbox
%
% Author: Xiongtao Ruan (04/18/2024)


try
    A = crop3d_mex(A, bbox);
catch ME
    disp(ME)
    disp(ME.stack);
    A = A(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6));            
end

end
