function [uuid] = get_uuid()
% generate uuid. First use java-based method, if fails, use system uuidgen
% 
% Author: Xiongtao Ruan (02/22/2020)
%
% xruan (11/08/2022): for none-java method on windows, use a random num string instead. 
% xruan (1/16/2023): for Windows, only keep first 4 letters to avoid too long file path issue

if usejava('jvm') 
    uuid = char(java.util.UUID.randomUUID.toString);
else
    if ispc
        rng('shuffle');
        uuid = num2str(randi(2^53-1));
    else
        [status, cmdout] = system('uuidgen');
        if status == 0
            uuid = strip(cmdout);
        else
            uuid = num2str(randi(2^53-1));
        end
    end
end

if ispc
    if numel(uuid) > 4
        uuid = uuid(1 : 4);
    end
end

end