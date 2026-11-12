# 0004: tiff and zarr through the published cpp-tiff and cpp-zarr packages

Status: proposed (task IO1); amended 2026-11-12 (task IO3) for the user's decisions on the
array convention and the zarr layout (`(z, y, x)` C-ordered arrays, native zarr outputs).
Checked with cpp-tiff 1.7.3 and cpp-zarr 1.7.3 from PyPI, NumPy 2.5.4, Python 3.12, Linux
x86-64, against the mex files PetaKit5D ships (built from the same 1.7.3 releases).

## Decision

`petakit5d.io` reads and writes tiff and zarr through the `cpptiff` and `cppzarr` packages,
with no reader or writer of its own. Arrays are NumPy `(z, y, x)` in C order (images `(y, x)`);
boxes are `(z0, y0, x0, z1, y1, x1)` and ranges `(start, stop)`, 0-based, end-exclusive; sizes
and per-axis values (block, shard and subfolder sizes) are `(z, y, x)`.

- **Tiff is unchanged on disk**: `write_tiff(vol)` writes the same bytes as `parallelWriteTiff`
  (and so as `writetiff` with LZW or zstd) for MATLAB's array of the same data,
  `to_matlab(vol)`, for all ten numeric types, 2D and 3D, single rows and columns, multi-strip
  pages, more than 65,535 pages and files over 4 GB, with one thread or many. The layout is
  BigTIFF, `RowsPerStrip` 512, eleven tags, nothing else. MATLAB's `readtiff` reads them as
  its `(y, x, z)` arrays.
- **Zarr written by the port is native**: Zarr v2, shape `[nz, ny, nx]`, `"order": "C"`, chunks
  `(cz, cy, cx)`, chunk files named by `(z, y, x)` chunk index (`iz.iy.ix`), blosc zstd level 1
  and the `.` separator by default; an image is stored as `[1, ny, nx]`. The `.zarray` text is
  the one cpp-zarr and the mex write for those values; zarr-python, napari and the OME-Zarr
  readers read the same `(z, y, x)` array. No OME-Zarr metadata is written (open with the
  user).
- **Zarr written by MATLAB is still read**: an array with `"order": "F"` is taken as MATLAB's
  `[ny, nx, nz]` (or `[ny, nx]`) and returned as `(z, y, x)`. `read_zarr(..., axes=...)`
  (`'auto'`, `'zyx'`, `'yxz'`) overrides the detection. Region writes into such arrays are
  refused; `write_zarr(..., overwrite=True)` replaces one with a native array.

## What changed on disk (2026-11-12)

| | before (IO2) | now (IO3) |
|---|---|---|
| tiff | the mex's files | the same files |
| zarr written by the port | MATLAB's layout: `[ny, nx, nz]`, order F, chunks `[cy, cx, cz]`, chunk files `iy.ix.iz` | native: `[nz, ny, nx]`, order C, chunks `[cz, cy, cx]`, chunk files `iz.iy.ix` |
| an image in zarr | `[ny, nx, 1]` | `[1, ny, nx]` |
| sharded zarr (`shardSize`) | MATLAB's v2 `.zarray` with a `sharding_indexed` codec | standard Zarr v3 with sharding (`zarr.json`; chunk files `c.iz.iy.ix`); a listed difference in files, the user's decision of 2026-11-12 |
| zarr MATLAB wrote | read as `(y, x, z)` without a copy | read as `(z, y, x)`: one threaded transposing pass |

**MATLAB PetaKit5D reads the port's new zarr files with their axes in another order**:
`readzarr` returns a `nz x ny x nx` array (its `getImageSize` says `[nz ny nx]`), the data
otherwise intact (each MATLAB dimension is one of ours). Mixing MATLAB and Python stages on
zarr data needs a conversion step; tiff data is exchanged as before.

## How PetaKit5D calls them

- **Tiff**: cpp-tiff's `(z, y, x)` arrays are used as they are, with no copy on either side
  (`read_tiff` returns the package's array; `write_tiff` passes C-ordered memory to the public
  `write_tiff`). A slice range `(start, stop)` is the package's `z_range`. One page comes back
  2D. MATLAB's tiff reader filled its layout while decoding, so a file-level task has one pass
  to spare.
- **Native zarr**: `write_zarr(name, a, start, stop)` into the array `create_zarr` made (as
  MATLAB's `writezarr` creates and then writes a region, also for a whole array), and
  `read_zarr(name, start, stop, order='C')`. Chunks are capped at the array size (as `writezarr`
  caps `blockSize`); MATLAB's defaults (`blockSize [500 500 500]`, `[y x z]`) go through
  `core.axes.from_yxz`. The package's own defaults (256 chunks, `/`) are not PetaKit5D's.
- **MATLAB-layout zarr**: `read_zarr(name, start, stop, order='F')` with the box in the file's
  `(y, x, z)` order, then each z plane transposed in the buffer the package returned
  (`petakit5d._ext.io.transpose_pages_inplace`: tiled, SSE2, OpenMP; one volume in memory),
  or `core.arrays.from_matlab` without the extension.
- **Sharded zarr**: `create_zarr(..., shard_size=...)` writes the `zarr.json` cpp-zarr writes for
  a v3 array with shards of `block_size` holding chunks of `shard_size` (tested against the
  package's own file, text and shard bytes), and the data goes in as regions. MATLAB's
  v2-sharded files (order F) are read through the package like its other arrays.
- **Conversions**: a region of another class written into an array is converted first with
  `core.dtypes.cast` (MATLAB's rounding and saturation, which the mex applies; the package
  would use NumPy's `astype`).

## What the port does itself until the packages close the gaps

| Gap (IO1 report has the reproducers) | Until fixed |
|---|---|
| No metadata-only creation (`createZarrFile`) | write the `.zarray` text (tested against the package's and the mex's), then fill by regions |
| Writes always skip all-zero chunks; `writezarr` writes them | a different set of chunk files for the same data (reads are identical); `sparse` is passed once the package takes it |
| Region writes convert with NumPy's `astype` (truncate, wrap) | `core.dtypes.cast` before writing |
| `uuid=False` not exposed | passed once the package takes it (speed on network storage only) |
| `get_image_shape` reports 1 page for single-IFD ImageJ stacks | the mex's rule from the first directory's description |
| float16, empty arrays and zero dimensions (a read past the buffer), compression names other than `lzw`/`zstd`/`none` (any other writes LZW), big-endian input, `os.PathLike`, ranges past the end | checked by `petakit5d.io` before calling |

## Not covered by the packages, and kept out of the port

`writetiff(..., 'Compression', 'none')` goes through MATLAB's `saveastiff` (classic TIFF,
64-bit integers stored as doubles). Python writes uncompressed BigTIFF through cpp-tiff
instead; no production caller writes uncompressed tiff.
