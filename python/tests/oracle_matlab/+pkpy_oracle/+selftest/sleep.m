function [x] = sleep(t)
% wait t seconds, then return t (for tests that stop a session mid-case)


pause(t);
x = t;

end
