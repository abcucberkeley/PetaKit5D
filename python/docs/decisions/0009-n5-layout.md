# 0009: N5 written by the port is standard N5 (`[nx, ny, nz]`)

Status: accepted by the user, 2026-11-17 ("Standard N5 is fine"; task B09). A deliberate
difference from MATLAB's files, as the native zarr layout of record 0004.

## Context

- MATLAB writes N5 through its `N5Adapter`, over zarr-python 2's `N5Store`. The store keeps
  the C-order array it is given and lists its shape reversed as the N5 `dimensions`. MATLAB's
  array is `(y, x, z)`, so its N5 datasets have `dimensions` `[nz, nx, ny]`, `blockSize`
  `[bz, bx, by]` and chunk files `iz/ix/iy`.
- Standard N5, as Fiji, BigStitcher, BigDataViewer, n5-java, tensorstore and zarr-python's
  N5 drivers read it, is x fastest: `dimensions` `[nx, ny, nz]`. Tools that read MATLAB's
  files show them with their axes permuted. MATLAB's BigStitcher XML (`constructXmlFile` with
  `dataOrder` `'yxz'`) compensates, by giving sizes, voxel sizes and translations in the
  permuted order.

## Decision

- **`zarr_to_n5` and `write_n5` write standard N5.** The port's `(z, y, x)` C array is the
  dataset's C-order array: `dimensions` `[nx, ny, nz]`, `blockSize` `[bx, by, bz]`, chunk files
  `ix/iy/iz`.
- **Everything else is as MATLAB writes it:**
  - whole chunks (edge chunks padded with zeros), big-endian;
  - blosc zstd level 1 with byte shuffle;
  - `attributes.json` as zarr-python 2's text (`"clevel": 1.0`, `"n5": "2.0.0"`).
- **Reading MATLAB's files needs `axes='yxz'`.** `open_n5`, `read_n5`, `n5_to_zarr` and
  `n5_to_zarr_block` take it (Python-only); the default is `'zyx'`.
- **The layout can't be detected.** N5 has no order flag, and MATLAB's `attributes.json` is
  what zarr-python 2 writes for any C-order array. The one trace of MATLAB is `"clevel": 1.0`, a
  float from MATLAB's double, but any zarr-python user who passes a float writes the same.
- **The BigStitcher XML defaults to the port's files.** `construct_xml_file`'s `data_order`
  default is `'zyx'`: sizes, voxel sizes and translations are listed in `[x, y, z]`. `'yxz'`
  (MATLAB's default) still describes MATLAB's files, and gives MATLAB's XML. The MATLAB-name
  layer passes `dataOrder` unchanged (kind `same`, `differs=`).

## What it means

| | MATLAB | the port |
|---|---|---|
| `dimensions` of a `(z, y, x)` volume | `[nz, nx, ny]` | `[nx, ny, nz]` |
| chunk files | `iz/ix/iy` | `ix/iy/iz` |
| `constructXmlFile` default `dataOrder` | `'yxz'` | `'zyx'` |
| N5 written by MATLAB | read as is | read with `axes='yxz'` |
| N5 from Fiji or BigStitcher (`[x, y, z]`) | read with its axes permuted | read as the image |

- **MATLAB reads the port's N5 permuted**, just as it reads the port's native zarr. Mixing
  MATLAB and Python stages on N5 data needs the `axes` option.
- **The oracle cases compare after the permutation.** The zarrToN5 cases compare the decoded
  voxels, and the attributes, chunk names and headers permuted from `[x, y, z]` to
  `[z, x, y]`. The N5ToZarr cases read MATLAB-layout inputs with `axes='yxz'`.
- **Checked with an independent reader.** zarr-python 2's `N5Store` reads the port's file as the
  image `(z, y, x)`, and the port reads its files, in both layouts.
- **To undo:** make `"yxz"` the default of `create_n5`/`open_n5` and `data_order`. Nothing else
  depends on it.
