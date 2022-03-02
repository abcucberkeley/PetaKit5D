function llsm_setup(codeRt, addPython, pythonPath, addPrivate)
% automatically load all code repositories to path
% 
% Author: Xiongtao Ruan
% Date: 10/05/2020


if nargin < 1 || isempty(codeRt)
    % code_rt = '/clusterfs/fiona/ABCCode/';
    codeRt = fileparts(which(mfilename));
end

if nargin < 2
    addPython = false;
end

if nargin < 3 || isempty(pythonPath)
    pythonPath = '~/anaconda3/bin/python';
end

if nargin < 4
    addPrivate = false;
end

if ~ispc
    [~, output] = system('echo $HOSTNAME');
end
fprintf('Hostname: %s \n', strip(output))
fprintf('Add matlab libraries to path...\n')
addpath(genpath([codeRt]));

% also add python libary
if addPython
    fprintf('Add python library to path...\n')
    try
        dir_info = dir(pythonPath);
        pe = pyenv;
        if ~contains(pe.Executable, 'miniconda3') && ~contains(pe.Executable, 'anaconda3') ...
                && ~strcmp(fileparts(pe.Executable), dir_info.folder)
            pyenv('Version', pythonPath);
        end
        % addpath(genpath('../third_parties/blockedDemo'));

        % add some components from matlab R2020b for older versions.
        if verLessThan('matlab', '9.9') 
            addpath(genpath([codeRt, '/../third_parties/matlab_R2020b']));
        end

        pymod_rel_path = [codeRt, '/microscopeDataProcessing/python/'];
        dir_info = dir(pymod_rel_path);
        pymod_path = dir_info.folder;
        insert(py.sys.path, int64(0), pymod_path);
        py.importlib.import_module('zarrAPI');
        py.importlib.import_module('daskAPI');
        % py.importlib.import_module('setZarrData');

        % insert dask functions
        % pymod_path = [codeRt, '/XiongtaoScripts/stitch/'];
        % insert(py.sys.path, int64(0), pymod_path);
        % py.importlib.import_module('daskZarrMaxProjection');
        % py.importlib.import_module('daskZarrPadArray');
    catch ME
        disp(ME);
    end
end

end

