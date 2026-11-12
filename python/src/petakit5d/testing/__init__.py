"""Test infrastructure of the port.

- ``cases``: oracle cases (the MATLAB call, the Python call, inputs, comparison)
- ``oracle``: golden store, MATLAB packs, repo manifests, the pytest check, the command line
  (``python -m petakit5d.testing.oracle``)
- ``compare``: MATLAB values in Python, value files, comparison with diffs, file trees

Nothing is imported here; import the submodules.
"""
