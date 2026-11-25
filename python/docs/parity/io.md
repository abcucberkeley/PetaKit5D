# Differences from MATLAB: io

The port reads and writes MATLAB's files (tiff and zarr through cpp-tiff and cpp-zarr); the file
format choices are rows F1, F2 and F8 of the [index](index.md) (decision record
[0004](../decisions/0004-io-bindings.md)). This page holds the io functions' rows beyond those.

## The arguments blocks of readtiff, writetiff, readzarr, writezarr and createzarr (A1; proposed 2026-11-22, for the user)

The five wrappers declare their inputs in `arguments` blocks (`range (1, :) {mustBeNumeric,
mustBeFinite}`, `groupWrite (1, 1) logical`, `filepath char`, ...: 11 declarations of 36
parameters). MATLAB converts a value to the declared class and size before the wrapper runs, then
validates it. The MATLAB-name layer (`petakit5d.matlab`, and the command line) does the same
since 2026-11-22 (`params/declarations.py`): a vector of any orientation (a column, 1x1xN) and
`[]` are made rows, numbers are made logical values (0 false, other numbers true), and MATLAB's
errors are raised with its identifiers and messages (a flag of two elements
`MATLAB:validation:IncompatibleSize`, NaN for a flag `MATLAB:validation:UnableToConvert`, a
matrix as a range, a box or a size `IncompatibleSize`, NaN or Inf in a range or a box
`MATLAB:validators:mustBeFinite`, `'LZW'` as the compression `MATLAB:validators:mustBeMember`).
The snake_case functions are unchanged. Of MATLAB's conversions, the user decided on 2026-11-22 to
leave out those that only hide mistakes; the port refuses those values with an error of its own:

| Row | Where | MATLAB | The port | Evidence |
|---|---|---|---|---|
| A1 | the MATLAB-name layer of `readtiff`, `writetiff`, `readzarr`, `writezarr` and `createzarr` | for a `char` declaration (the file paths, `dtype`, `order`, `compressor`, `dimSeparator`, `Compression`, `Mode`): numbers are taken as the characters of their codes (truncated, saturated at 0 and 65535, NaN 0), so `writetiff(img, 65)` writes a file named `A`, `'dimSeparator', 47` makes an array with the separator `'/'` and `[]` is `''`; a cell of char arrays is made a char matrix of its texts (`{'c.tif'}` is `'c.tif'`, `{}` is `''`); a function handle is taken as its name (`writetiff(img, @sin)` writes `sin`). For a `logical` declaration (the flags, writezarr's `sparseData` and `create`): characters are taken as true, all but `char(0)` (`'a'` and `'0'` are true). Where the converted value fails a later check, MATLAB raises that check's error (`IncompatibleSize`, `mustBeMember`) | refuses each of these values with `ConversionRefusedError` (`petakit5d:params:conversionRefused`), before the wrapper writes anything; every other value is taken, converted or refused as MATLAB's arguments block does, with MATLAB's identifiers and messages (`MatlabValidationError`) | `tests/params/test_declarations.py` (oracle job 4642229): per declaration, the values the port converts as MATLAB does (`decl.<key>`, 13 cases of 47 to 106 values, the same outcome for each in every form the registry takes) and the values it refuses (`decl.<key>.refused`, 8 cases of 21 to 39 values: MATLAB's outcomes recorded, none `UnableToConvert`, the port's refusal checked); `tests/io/test_io_arguments.py`: the wrappers themselves, 20 calls whose conversions or errors the port copies, and 8 it refuses (`args.refused.*`: the files MATLAB writes recorded, nothing written by the port) |
