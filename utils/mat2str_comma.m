function [S] = mat2str_comma(A, sn)
% convert matrix to string with comma separation for given significant number
% A: input vector; sn: significant number, assume integers if not provided
% S: output char string


if nargin == 1 || isempty(sn)
    S = strrep(mat2str(A), ' ', ',');
else
    S = sprintf('[%s]', strip(num2str(A, sprintf('%%0.%df,', sn)), ','));
end

end