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
% (an environment without zarr-python 2 failed in the call, MATLAB:undefinedVarOrClass of
% py.zarr.N5Store; the switch to 45b2bb1 ran the cases in one first)
if ~py.hasattr(py.importlib.import_module('zarr'), 'N5Store')
    error('pkpy:n5:noZarr2', ['MATLAB''s N5Adapter needs zarr-python 2 (N5Store) in the ', ...
        'active environment: run the N5 cases with orun.sh -e <an environment with zarr 2>']);
end
p = fullfile(fileparts(which('setup')), 'microscopeDataProcessing', 'python');
if count(py.sys.path, p) == 0
    insert(py.sys.path, int64(0), p);
end
py.importlib.import_module('zarrAPI');

end
