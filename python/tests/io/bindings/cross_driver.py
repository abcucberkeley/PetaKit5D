"""Cross-read and byte matrix between the cpp-tiff / cpp-zarr Python packages and the
MATLAB mex files (task IO1). Stand-in for the oracle harness (F04) until it lands.

    python cross_driver.py prepare <dir> [--big]   inputs, manifest, Python-written files
    matlab -batch "setup([]); io1_cross_mex('<dir>')"   (from the oracle tree, matlab/ on the path)
    python cross_driver.py check <dir>             Python reads MATLAB's files, compares bytes

check writes <dir>/results.json (every check of every case, both sides) and prints a summary.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import time
import traceback
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import iobind_cases as ic  # noqa: E402


def _data(c, folder):
    return ic.read_raw(os.path.join(folder, 'inputs', c['name'] + '.bin'), c['dtype'], c['shape'])


def _zarr_kwargs(c, threads=None):
    shape = tuple(c['shape']) + (1,) * (3 - len(c['shape']))
    chunks = [min(s, b) for s, b in zip(shape, c['chunks'])]
    kw = dict(order=c['order'], chunks=chunks, dimension_separator=c['sep'], cname=c['cname'], clevel=1, num_threads=threads)
    if c.get('shard'):
        # MATLAB's blockSize is the shard, shardSize the inner chunk
        kw.update(zarr_format=3, shards=chunks, chunks=[min(b, s) for b, s in zip(chunks, c['shard'])])
    return kw


def _as3d(a):
    return a.reshape(a.shape + (1,) * (3 - a.ndim), order='F')


def _region_inputs(c, folder):
    regs = []
    pats = c.get('region_patterns') or ['random'] * len(c['regions'])
    for i, (b, pat) in enumerate(zip(c['regions'], pats)):
        sz = [b[3] - b[0] + 1, b[4] - b[1] + 1, b[5] - b[2] + 1]
        if c['writer'] == 'convert':
            d = ic.convert_values(sz)
        else:
            d = ic.make_data(c['dtype'], sz, pat, ic.seed_of(f"{c['name']}_r{i}"))
        regs.append((b, d))
    return regs


def _expected_regions(c, regs):
    exp = np.zeros(c['shape'], dtype=c['dtype'], order='F')
    for b, d in regs:
        exp[b[0] - 1:b[3], b[1] - 1:b[4], b[2] - 1:b[5]] = ic.matlab_uint16(d) if c['writer'] == 'convert' else d
    return exp


def _py_create(fn, c, threads=None):
    """Stand-in for createzarr: an all-zero write, which leaves only .zarray (the
    Python package has no metadata-only creation and skips all-zero chunks)."""
    zeros = np.zeros(tuple(c['shape']) + (1,) * (3 - len(c['shape'])), dtype=c['dtype'], order='F')
    import cppzarr
    cppzarr.write_zarr(str(fn), zeros, **_zarr_kwargs(c, threads))


def _py_regions(fn, regs):
    import cppzarr
    for b, d in regs:
        cppzarr.write_zarr(str(fn), d, [b[0] - 1, b[1] - 1, b[2] - 1], [b[3], b[4], b[5]])


def prepare(folder, big=False):
    import cppzarr
    folder = Path(folder)
    cases = ic.tiff_cases(big) + ic.zarr_cases(big)
    if not big:
        cases = [c for c in cases if not c.get('big')]
    else:
        cases = [c for c in cases if c.get('big')]
    ic.manifest_write(folder, cases)
    timing = {}
    for c in cases:
        t0 = time.perf_counter()
        name = c['name']
        data = ic.make_data(c['dtype'], c['shape'], c['pattern'], ic.seed_of(name))
        ic.write_raw(folder / 'inputs' / (name + '.bin'), data)
        if c['kind'] == 'tiff':
            for sub, fn_write, th in (('py', ic.tiff_write_yxz, None), ('py1', ic.tiff_write_yxz, 1),
                                      ('pydirect', ic.tiff_write_yxz_direct, None)):
                if c.get('big') and sub == 'py1':
                    continue
                (folder / sub).mkdir(exist_ok=True)
                fn_write(folder / sub / (name + '.tif'), data, c['compression'], th)
        else:
            for sub, th in (('py', None), ('py1', 1)):
                if c.get('big') and sub == 'py1':
                    continue
                (folder / sub).mkdir(exist_ok=True)
                fn = folder / sub / (name + '.zarr')
                if c['writer'] in ('regions', 'convert'):
                    regs = _region_inputs(c, folder)
                    if sub == 'py':
                        for i, (b, d) in enumerate(regs):
                            ic.write_raw(folder / 'inputs' / f'{name}_r{i}.bin', d)
                        ic.write_raw(folder / 'inputs' / f'{name}_expected.bin', _expected_regions(c, regs))
                        _py_create(folder / sub / (name + '_base.zarr'), c)
                    _py_create(fn, c, th)
                    _py_regions(fn, regs)
                elif c['writer'] == 'createonly':
                    _py_create(fn, c, th)
                else:
                    cppzarr.write_zarr(str(fn), _as3d(data), **_zarr_kwargs(c, th))
                    if c['shape'] == [21, 17] or tuple(c['shape']) == (21, 17):
                        # what the package does with 2D data as it is
                        cppzarr.write_zarr(str(folder / sub / (name + '_as2d.zarr')), data, **dict(_zarr_kwargs(c, th), chunks=list(c['shape'])))
        timing[name] = time.perf_counter() - t0
    with open(folder / 'prepare_timing.json', 'w') as f:
        json.dump(timing, f)
    print(f'prepared {len(cases)} cases in {folder}')


# --- check -------------------------------------------------------------------------------------

def _file_bytes_equal(a, b):
    if not (os.path.isfile(a) and os.path.isfile(b)):
        return None
    if os.path.getsize(a) != os.path.getsize(b):
        return False
    with open(a, 'rb') as fa, open(b, 'rb') as fb:
        while True:
            x, y = fa.read(1 << 24), fb.read(1 << 24)
            if x != y:
                return False
            if not x:
                return True


def _chunk_files(fn):
    out = []
    for root, _, files in os.walk(fn):
        for f in files:
            if not f.startswith('.') and f != 'zarr.json':
                out.append(os.path.relpath(os.path.join(root, f), fn).replace(os.sep, '/'))
    return sorted(out)


def _meta_file(fn):
    for m in ('.zarray', 'zarr.json'):
        if os.path.isfile(os.path.join(fn, m)):
            return os.path.join(fn, m)
    return None


def _zarr_compare(a, b):
    """(metadata text equal, chunk file sets equal, all chunk bytes equal, detail)"""
    ma, mb = _meta_file(a), _meta_file(b)
    meta = ma is not None and mb is not None and os.path.basename(ma) == os.path.basename(mb) and \
        open(ma, 'rb').read() == open(mb, 'rb').read()
    ca, cb = _chunk_files(a), _chunk_files(b)
    same_set = ca == cb
    differ = [f for f in sorted(set(ca) & set(cb)) if not _file_bytes_equal(os.path.join(a, f), os.path.join(b, f))]
    detail = ''
    if not same_set:
        detail += f'chunks only in first: {sorted(set(ca) - set(cb))[:5]} ({len(set(ca) - set(cb))}); ' \
                  f'only in second: {sorted(set(cb) - set(ca))[:5]} ({len(set(cb) - set(ca))}); '
    if differ:
        detail += f'{len(differ)} of {len(set(ca) & set(cb))} common chunks differ (e.g. {differ[:3]})'
    return meta, same_set, not differ, detail


def check(folder):
    import cpptiff
    import cppzarr
    folder = Path(folder)
    man = json.load(open(folder / 'manifest.json'))
    mex_results = json.load(open(folder / 'results_mex.json')) if (folder / 'results_mex.json').exists() else []
    if isinstance(mex_results, dict):
        mex_results = [mex_results]
    recs = []

    def rec(c, check, ok, detail=''):
        status = 'n/a' if ok is None else ('pass' if ok else 'fail')
        recs.append(dict(case=c['name'], kind=c['kind'], check=check, side='python', status=status, detail=str(detail)[:500]))

    def guarded(c, check, f):
        try:
            r = f()
            if isinstance(r, tuple):
                rec(c, check, r[0], r[1])
            else:
                rec(c, check, r)
        except Exception as e:  # a failing binding call is a result, not a crash
            recs.append(dict(case=c['name'], kind=c['kind'], check=check, side='python', status='error',
                             detail=f'{type(e).__name__}: {str(e)[:300]}'))

    for c in man['cases']:
        name = c['name']
        shape = tuple(c['shape'])
        data = _data(c, folder)
        if c['kind'] == 'tiff':
            mex, wt, mex1 = (folder / s / (name + '.tif') for s in ('mex', 'wt', 'mex1'))
            py, py1, pyd = (folder / s / (name + '.tif') for s in ('py', 'py1', 'pydirect'))
            guarded(c, 'py_reads_mex', lambda: ic.bits_equal_m(ic.tiff_read_yxz(mex), data))
            guarded(c, 'py_reads_writetiff', lambda: ic.bits_equal_m(ic.tiff_read_yxz(wt), data))
            if c['range']:
                r0, r1 = c['range']
                guarded(c, 'py_reads_mex_range', lambda: ic.bits_equal_m(ic.tiff_read_yxz(mex, [r0 - 1, r1]), data[:, :, r0 - 1:r1]))
            nz = shape[2] if len(shape) == 3 else 1
            guarded(c, 'py_shape_of_mex', lambda: (tuple(cpptiff.get_image_shape(str(mex))) == (nz, shape[0], shape[1]),
                                                    cpptiff.get_image_shape(str(mex))))
            guarded(c, 'bytes_py_eq_mex', lambda: _file_bytes_equal(py, mex))
            guarded(c, 'bytes_py1_eq_mex1', lambda: _file_bytes_equal(py1, mex1))
            guarded(c, 'bytes_pydirect_eq_py', lambda: _file_bytes_equal(pyd, py))
            guarded(c, 'bytes_writetiff_eq_mex', lambda: _file_bytes_equal(wt, mex))
            if not c.get('big'):
                guarded(c, 'tifffile_reads_py', lambda: _tifffile_check(py, data))
                guarded(c, 'tifffile_reads_mex', lambda: _tifffile_check(mex, data))
        else:
            mex, mex1 = folder / 'mex' / (name + '.zarr'), folder / 'mex1' / (name + '.zarr')
            py, py1 = folder / 'py' / (name + '.zarr'), folder / 'py1' / (name + '.zarr')
            w = c['writer']
            if w in ('regions', 'convert'):
                exp = ic.read_raw(folder / 'inputs' / f'{name}_expected.bin', c['dtype'], shape)
                regs = [(b, ic.read_raw(folder / 'inputs' / f'{name}_r{i}.bin', 'float64' if w == 'convert' else c['dtype'],
                                        [b[3] - b[0] + 1, b[4] - b[1] + 1, b[5] - b[2] + 1])) for i, b in enumerate(c['regions'])]
                guarded(c, 'py_regions_eq_expected', lambda: ic.bits_equal(cppzarr.read_zarr(str(py)), exp))
                guarded(c, 'py_reads_mex_regions', lambda: ic.bits_equal(cppzarr.read_zarr(str(mex)), exp))
                pyonmex = folder / 'pyonmex' / (name + '.zarr')

                def regions_on_mex_base():
                    if pyonmex.exists():
                        shutil.rmtree(pyonmex)
                    shutil.copytree(folder / 'mex' / (name + '_base.zarr'), pyonmex)
                    _py_regions(pyonmex, regs)
                    return ic.bits_equal(cppzarr.read_zarr(str(pyonmex)), exp)
                guarded(c, 'py_regions_on_mex_created', regions_on_mex_base)
                guarded(c, 'zarray_py_created_eq_mex_created',
                        lambda: _zarr_compare(folder / 'py' / (name + '_base.zarr'), folder / 'mex' / (name + '_base.zarr'))[0])
                guarded(c, 'chunk_set_py_eq_mex', lambda: (lambda r: (r[1], r[3]))(_zarr_compare(py, mex)))
                guarded(c, 'chunk_bytes_py1_eq_mex1', lambda: (lambda r: (r[2] if r[1] else None, r[3]))(_zarr_compare(py1, mex1)))
            elif w == 'createonly':
                guarded(c, 'zarray_py_eq_mex', lambda: (lambda r: (r[0], open(_meta_file(mex)).read() if _meta_file(mex) else 'no metadata'))(_zarr_compare(py, mex)))
                guarded(c, 'chunk_set_py_eq_mex', lambda: (lambda r: (r[1], r[3]))(_zarr_compare(py, mex)))
                guarded(c, 'py_reads_mex_created_as_zeros', lambda: ic.bits_equal(cppzarr.read_zarr(str(mex)), np.zeros(shape, c['dtype'], order='F')))
            else:
                exp = _as3d(data)
                guarded(c, 'py_reads_mex', lambda: ic.bits_equal(cppzarr.read_zarr(str(mex)), exp))
                guarded(c, 'py_reads_mex_is_F', lambda: cppzarr.read_zarr(str(mex)).flags['F_CONTIGUOUS'] or c['order'] == 'C')
                if c['bboxes']:
                    def boxes():
                        bad = [b for b in c['bboxes'] if not ic.bits_equal(
                            cppzarr.read_zarr(str(mex), [b[0] - 1, b[1] - 1, b[2] - 1], [b[3], b[4], b[5]]),
                            exp[b[0] - 1:b[3], b[1] - 1:b[4], b[2] - 1:b[5]])]
                        return not bad, bad
                    guarded(c, 'py_reads_mex_bbox', boxes)
                guarded(c, 'zarray_py_eq_mex', lambda: (lambda r: (r[0], '' if r[0] else _meta_diff(py, mex)))(_zarr_compare(py, mex)))
                guarded(c, 'chunk_set_py_eq_mex', lambda: (lambda r: (r[1], r[3]))(_zarr_compare(py, mex)))
                guarded(c, 'chunk_bytes_py_eq_mex', lambda: (lambda r: (r[2], r[3]))(_zarr_compare(py, mex)))
                if not c.get('big'):
                    guarded(c, 'chunk_bytes_py1_eq_mex1', lambda: (lambda r: (r[2], r[3]))(_zarr_compare(py1, mex1)))
                    guarded(c, 'chunk_bytes_py_eq_py1', lambda: (lambda r: (r[2], r[3]))(_zarr_compare(py, py1)))
                    guarded(c, 'zarr_python_reads_mex', lambda: _zarr_python_check(mex, exp))
                if (folder / 'py' / (name + '_as2d.zarr')).exists():
                    guarded(c, 'zarray_py2d_eq_mex', lambda: (_zarr_compare(folder / 'py' / (name + '_as2d.zarr'), mex)[0],
                                                               _meta_diff(folder / 'py' / (name + '_as2d.zarr'), mex)))
    for r in mex_results:
        r = dict(r)
        r['side'] = 'matlab'
        recs.append(r)
    with open(folder / 'results.json', 'w') as f:
        json.dump({'versions': ic.versions(), 'mex_versions': man.get('versions'), 'records': recs}, f, indent=1)
    summarize(recs)
    return recs


def _meta_diff(a, b):
    ma, mb = _meta_file(a), _meta_file(b)
    if not ma or not mb:
        return f'metadata missing: {ma} {mb}'
    ja, jb = json.load(open(ma)), json.load(open(mb))
    diff = {k: (ja.get(k), jb.get(k)) for k in sorted(set(ja) | set(jb)) if ja.get(k) != jb.get(k)}
    return f'fields differ: {diff}' if diff else 'same fields, text differs (formatting)'


def _tifffile_check(fn, data):
    import tifffile
    a = tifffile.imread(str(fn))
    return ic.bits_equal_m(ic.zyx_to_yxz(np.ascontiguousarray(a)) if a.ndim == 3 else np.asfortranarray(a), data)


def _zarr_python_check(fn, exp):
    import zarr
    a = zarr.open_array(str(fn), mode='r')[...]
    return ic.bits_equal(np.asfortranarray(a), exp)


def summarize(recs):
    from collections import Counter, defaultdict
    by = defaultdict(Counter)
    for r in recs:
        by[(r['side'], r['check'])][r['status']] += 1
    print(f"{'side':7s} {'check':38s} pass fail error n/a")
    for (side, chk), cnt in sorted(by.items()):
        print(f"{side:7s} {chk:38s} {cnt['pass']:4d} {cnt['fail']:4d} {cnt['error']:5d} {cnt['n/a']:3d}")
    bad = [r for r in recs if r['status'] in ('fail', 'error')]
    print(f'{len(recs)} checks, {len(bad)} fail or error')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('step', choices=['prepare', 'check'])
    p.add_argument('folder')
    p.add_argument('--big', action='store_true')
    a = p.parse_args()
    try:
        prepare(a.folder, a.big) if a.step == 'prepare' else check(a.folder)
    except Exception:
        traceback.print_exc()
        sys.exit(1)
