# 0004: tiff and zarr through the published cpp-tiff and cpp-zarr packages

Status: proposed (task IO1). Checked with cpp-tiff 1.7.3 and cpp-zarr 1.7.3 from PyPI,
NumPy 2.5.4, Python 3.12, Linux x86-64, against the mex files PetaKit5D ships (built from
the same 1.7.3 releases).

## Decision

`petakit5d.io` reads and writes tiff and zarr through the `cpptiff` and `cppzarr` packages,
with no reader or writer of its own. The files are the mex's files:

- **Tiff**: `write_tiff` writes the same bytes as `parallelWriteTiff` (and so as `writetiff`
  with LZW or zstd) for all ten numeric types, 2D and 3D, single rows and columns, multi-strip
  pages, more than 65,535 pages and files over 4 GB, with one thread or many. The layout is
  BigTIFF, `RowsPerStrip` 512, eleven tags, nothing else.
- **Zarr**: with PetaKit5D's settings (below) `write_zarr` writes the same `.zarray` text and
  the same chunk bytes as `writezarr`/`createzarr`, for every compressor, F and C order, both
  separators, with one thread or many. Each package reads every file the other side writes,
  whole and by region, including the mex's v2 sharded arrays.

## How PetaKit5D calls them

- Volumes are `(y, x, z)` Fortran order (PLAN 3.2). Zarr needs no conversion: `read_zarr`
  returns Fortran-ordered arrays for order-F files and `write_zarr` takes them as they are.
- Tiff arrays are `(z, y, x)` C order in the package. A read is followed by a copy to
  `(y, x, z)` Fortran order, split over the call's threads. A write passes the Fortran memory
  to the extension's `transpose=True` entry (`pybind11_write_tiff(name, a.transpose(1, 0, 2),
  True, compression, 0, threads)`), which is what the mex does and writes identical bytes
  without a copy; `write_tiff` itself would copy. Both are stand-ins until cpp-tiff offers
  `(y, x, z)` reads and writes in its API (requested).
- `write_zarr(name, a, order='F', chunks=..., dimension_separator='.', cname='zstd', clevel=1)`,
  with chunks capped at the array size (as `writezarr` caps `blockSize`) and 2D data given a
  third axis of 1 (PetaKit5D stores `[ny, nx, 1]`; the package keeps 2D arrays 2D). The
  package's own defaults (256 chunks, `/`) are not PetaKit5D's.
- Regions: MATLAB's `[ymin, xmin, zmin, ymax, xmax, zmax]` (1-based, inclusive) is
  `start = [ymin-1, xmin-1, zmin-1]`, `end = [ymax, xmax, zmax]`; tiff `range [a, b]` is
  `z_range = [a-1, b]`. One page or one slice comes back 2D, as MATLAB's size says.

## What the port does itself until the packages close the gaps

| Gap (IO1 report has the reproducers) | Until fixed |
|---|---|
| No metadata-only creation (`createZarrFile`) | write the `.zarray` text, identical to the mex's (tested), then fill by regions |
| Writes always skip all-zero chunks; `writezarr` writes them | a different set of chunk files for the same data (reads are identical); needs a `sparse` option |
| Region writes convert with NumPy's `astype` (truncate, wrap) | convert with MATLAB's rounding (`petakit5d.core`) before writing |
| `uuid=False` not exposed | none (speed on network storage only; no production caller uses it) |
| Shards only as Zarr v3; `createzarr` with `shardSize` writes a v2 sharded `.zarray` | no Python path for `shardSize` yet |
| `read_tiff` cannot fill `(y, x, z)`; `write_tiff` copies | the copy and the extension entry above |
| `get_image_shape` reports 1 page for single-IFD ImageJ stacks | size such files by reading them, or the mex's rule |
| float16, empty arrays, unknown compression names, big-endian input, `os.PathLike` | checked by `petakit5d.io` before calling |

## Not covered by the packages, and kept out of the port

`writetiff(..., 'Compression', 'none')` goes through MATLAB's `saveastiff` (classic TIFF,
64-bit integers stored as doubles). Python writes uncompressed BigTIFF through cpp-tiff
instead; no production caller writes uncompressed tiff.
