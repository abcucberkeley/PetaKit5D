function [uuid] = get_uuid()
% generate uuid. First use java-based method, if fails, use system uuidgen
% 
% Author: Xiongtao Ruan (02/22/2020)
%
% xruan (11/08/2022): for none-java method on windows, use a random num string instead. 

if usejava('jvm') 
    uuid = char(java.util.UUID.randomUUID.toString);
else
    if ispc
        rng('shuffle');
        uuid = num2str(randi(2^53-1));
    else
        [~, cmdout] = system('uuidgen');
        uuid = strip(cmdout);
    end
end

end