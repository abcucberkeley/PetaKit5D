function setup(codeRt, addPython, pythonPath, addPrivate)
% automatically load all code repositories to path
% 
% Author: Xiongtao Ruan
% Date: 10/05/2020


if nargin < 1 || isempty(codeRt)
    codeRt = fileparts(which(mfilename));
end

if nargin < 2
    addPython = false;
end

% show hostname
if ispc
    [~, output] = system('hostname');
else
    [~, output] = system('echo $HOSTNAME');
end
fprintf('Hostname: %s \n', strip(output));

if nargin < 3
    pythonPath = pyenv().Executable;
    % https://www.mathworks.com/help/matlab/matlab_external/install-supported-python-implementation.html
end

if addPython && isempty(pythonPath)
    if ispc
        error('Please provide your python install path!');
    else
        pythonPath = '~/anaconda3/bin/python';
        fprintf('Use default python: %s .\n', pythonPath);        
    end
end

if nargin < 4
    addPrivate = false;
end

fprintf('Add matlab libraries to path...\n')
addpath(genpath([codeRt]));
rmpath(genpath([codeRt, '/mcc/mac']));

% also add python libary
if addPython
    fprintf('Add python library to path...\n')
    try
        dir_info = dir(pythonPath);
        pe = pyenv;
        if pe.Status ~= "Loaded" && ~contains(pe.Executable, 'miniconda3', 'IgnoreCase', true) && ...
                ~contains(pe.Executable, 'anaconda3', 'IgnoreCase', true) && ...
                ~strcmp(fileparts(pe.Executable), dir_info.folder)
            pyenv('Version', pythonPath);
        end

        % resolve license issue when using multiprocessing on Windows
        if ispc
            py.multiprocessing.spawn.set_executable(pe.Executable)
        end
            
        pymod_rel_path = [codeRt, '/microscopeDataProcessing/python/'];
        dir_info = dir(pymod_rel_path);
        pymod_path = dir_info.folder;
        insert(py.sys.path, int64(0), pymod_path);
        py.importlib.import_module('zarrAPI');
        py.importlib.import_module('daskAPI');
    catch ME
        disp(ME);
    end
end

end

