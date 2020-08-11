function [uuid] = get_uuid()
% generate uuid. First use java-based method, if fails, use system uuidgen
% 
% Author: Xiongtao Ruan (02/22/2020)

try 
    temp = java.util.UUID.randomUUID;
    uuid = char(temp.toString);
catch
    [status, cmdout] = system('uuidgen');
    uuid = strip(cmdout);
end

end