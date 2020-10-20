function [uuid] = get_uuid()
% generate uuid. First use java-based method, if fails, use system uuidgen
% 
% Author: Xiongtao Ruan (02/22/2020)

try 
    uuid = char(java.util.UUID.randomUUID.toString);
catch
    [status, cmdout] = system('uuidgen');
    uuid = strip(cmdout);
end

end