function [output] = get_hostname()
% show hostname
%
% Author: Xiongtao Ruan (03/31/2023)

if ispc
    [~, output] = system('hostname');
else
    [~, output] = system('echo $HOSTNAME');
end
fprintf('Hostname: %s \n', strip(output));

end