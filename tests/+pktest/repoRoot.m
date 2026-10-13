function [root] = repoRoot()
% absolute path of the PetaKit5D repository root

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));

end
