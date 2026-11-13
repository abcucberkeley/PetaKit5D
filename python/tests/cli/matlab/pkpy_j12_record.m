function out = pkpy_j12_record(action, value)
% The arguments the recorder of pkpy_j12_master got: ('set', args) keeps them, ('get') returns them.

persistent got
switch action
    case 'set'
        got = value;
        out = [];
    case 'get'
        out = got;
end
end
