function pkpy_n5_setup()
% Python for MATLAB's N5Adapter in an oracle session: the interpreter of the active virtual
% environment (it needs zarr-python 2 and numcodecs) and the oracle tree's zarrAPI module.

pe = pyenv;
if pe.Status ~= "Loaded"
    venv = getenv('VIRTUAL_ENV');
    if ~isempty(venv)
        pyenv('Version', fullfile(venv, 'bin', 'python'));
    end
end
p = fullfile(fileparts(which('setup')), 'microscopeDataProcessing', 'python');
if count(py.sys.path, p) == 0
    insert(py.sys.path, int64(0), p);
end
py.importlib.import_module('zarrAPI');

end
