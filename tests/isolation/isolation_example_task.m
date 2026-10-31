function [] = isolation_example_task(outFn, leak)
% a task for testing the isolation harness: it saves a value to outFn that
% depends on the process state it leaks into, as a task with that leak would
%
%   leak : 'rng' (saves random numbers of the global stream), 'persistent'
%          (counts its calls in a persistent variable), 'warning' (turns a
%          warning off), 'cd' (changes the current folder), 'file' (leaves
%          a file open) or 'none'.


persistent nCalls
if isempty(nCalls)
    nCalls = 0;
end
nCalls = nCalls + 1;

value = 0;
switch leak
    case 'rng'
        value = rand(1, 3);
    case 'persistent'
        value = nCalls;
    case 'warning'
        warning('off', 'isolation:example');
    case 'cd'
        cd(tempdir);
end
save(outFn, 'value');
if strcmp(leak, 'file')
    fopen(outFn, 'r');
end

end
