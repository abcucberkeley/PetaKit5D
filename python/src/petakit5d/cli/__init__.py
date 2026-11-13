"""The command lines of PetaKit5D in Python (docs/cli.md).

- ``pk5d`` (``petakit5d.cli.main``): the snake_case functions with JSON values, ``pk5d list``,
  ``pk5d doc <function>``, and the helpers ``pk5d parse``, ``pk5d master`` and ``pk5d mcc-shim``.
- ``pk5d-master`` (``petakit5d.cli.master``): the compiled MATLAB master's command line
  (``mccMaster <function> <MATLAB-style strings>``), for the GUI and PyPetaKit5D, with a
  ``run_mccMaster.sh`` shim that ignores the MATLAB Runtime argument.

Importing this package loads nothing else, so that the master starts fast.
"""
