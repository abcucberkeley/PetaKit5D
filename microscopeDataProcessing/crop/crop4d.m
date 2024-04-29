function [A] = crop4d(A, bbox)
% crop A according to bbox
%
% Author: Xiongtao Ruan (04/24/2024)


try
    A = crop4d_mex(A, bbox);
catch ME
    disp(ME)
    disp(ME.stack);
    A = A(bbox(1) : bbox(5), bbox(2) : bbox(6), bbox(3) : bbox(7), bbox(4) : bbox(8));
end

end
