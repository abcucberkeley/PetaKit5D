# Benchmarks

Benchmarks of the port against MATLAB PetaKit5D, run by the benchmark harness
(`petakit5d.testing.bench`; method, spec format and gates in `../docs/benchmarks.md`).

| Folder | What |
|---|---|
| `common/` | the harness's MATLAB helpers (`pk5d_bench_mark.m`, `pk5d_bench_record.m`, `pk5d_bench_info.m`), the process-cost specs and probe, the self-check, the harness's tests |
| `io/` | reads and writes of the tiff and zarr packages against the mex (task IO1) |

One spec per comparison, in the folder of its area (`benchmarks/<area>/<name>.json`): MATLAB
and Python variants that run the same work on the same inputs, timed inside the processes with
`record` / `pk5d_bench_record` under the same op names and thread counts. Run one with

```
python -m petakit5d.testing.bench run benchmarks/<area>/<name>.json --out <results folder>
```

from `python/`, with the oracle tree in `PETAKIT5D_ORACLE_DIR` for the MATLAB side. Numbers
that count come from one job per comparison on a pinned node type (the cluster helper does
that and keeps the results with their history).

The harness's own tests: `python -m pytest tests/testing/test_bench.py` (run by CI with the rest of the suite).
