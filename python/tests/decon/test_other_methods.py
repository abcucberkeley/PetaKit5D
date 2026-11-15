"""RLdecon's 'original' method (MATLAB's deconvlucy, petakit5d.decon.other_methods) against
MATLAB's deconvlucy at the oracle commit, as a black box.

Oracle cases (inputs made in MATLAB's form; the roles give the Python calls the port's form):

- lucy.*: the call RLdecon makes, ``deconvlucy(max(single(raw) - Background, 0), psf, nIter)``
  with the PSF normalized by its double sum: 2D and 3D data, odd and even sizes, a PSF cropped to
  the data along z, 1, 2, 3, 5, 10 and 15 iterations.
- class.*: the classes deconvlucy takes (uint8, uint16, int16, single, double; single and double
  PSFs), the default NUMIT.
- opt.*: DAMPAR, WEIGHT, READOUT and their combination (deconvlucy's documented parameters; RLdecon
  passes none of them), [] for a default, and errors.
- cell.*: the cell-array form (the image, the last and the previous estimates and the internal
  array) and an interrupted run resumed from it.
- rldecon.*: RLdecon with RLMethod 'original' in memory.
- otf.*: MATLAB's psf2otf (deconvlucy's OTF): real for a symmetric PSF, and where it turns complex.
- Each case runs on both FFT backends: ``.exact`` (FFTW with MATLAB's driver; skipped without the
  add-on or where FFTW's plans differ from the reference platform's) and ``.permissive``
  (decision record 0006's tolerances). On the exact backend single, integer and logical results
  are compared bit for bit; double results within ``DOUBLE_TOL`` (MATLAB's acceleration factor
  comes from BLAS dot products, whose order is the CPU's: see ``decon.other_methods``).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from petakit5d import fft as P
from petakit5d.core.arrays import to_matlab
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, File, In, InPath, Value, Vol
from petakit5d.testing.compare import tol
from petakit5d.testing.fft_reference import digests_hold_here
from petakit5d.testing.synth import CounterRNG, beads, decon_psf, write_tiff

BACKENDS = ("exact", "permissive")
DOUBLE_TOL = tol(
    reason="proposed (B11): deconvlucy's double results on the exact FFT backend, whose "
    "acceleration factor MATLAB takes from BLAS dot products in the CPU's order (max <= 1e-13 of "
    "the largest value, relative L2 <= 1e-13)",
    rel_max=1e-13,
    rel_l2=1e-13,
)
PERMISSIVE_TOL = tol(
    reason="0006: deconvolution with the permissive FFT backend (relative L2 <= 1e-4)", rel_l2=1e-4
)
PERMISSIVE_INT = tol(
    reason="0006: integer deconvolution with the permissive FFT backend (at most 1 count in at "
    "most 0.1% of the voxels)",
    abs=1,
    max_frac=1e-3,
)

BLUR_SIGMA = (math.sqrt(0.8**2 + 2.0**2), math.sqrt(0.8**2 + 1.5**2))  # (z, yx)


def _py_size(msize) -> tuple[int, ...]:
    sz = tuple(int(v) for v in msize)
    return (sz[2], sz[0], sz[1]) if len(sz) == 3 else sz


def raw_beads(msize, seed: int) -> np.ndarray:
    """beads blurred by the test PSF plus some noise on a background of 100, uint16, (z, y, x)"""
    py = _py_size(msize)
    rng = CounterRNG(seed, "lucy")
    n = max(3, round(math.prod(py) / 4000))
    pos = rng.uniform((n, len(py))) * (np.array(py, dtype=np.float64) - 1)
    im = beads(py, pos, BLUR_SIGMA if len(py) == 3 else BLUR_SIGMA[1], 1000.0, 100.0)
    im = im + np.float32(30) * rng.uniform(py, np.float32)
    return np.floor(im + np.float32(0.5)).astype(np.uint16)  # im >= 0


def rldecon_data(msize, seed: int) -> np.ndarray:
    """what RLdecon passes: max(single(raw) - 100, 0) (exact in single); MATLAB's form"""
    a = raw_beads(msize, seed).astype(np.float32) - np.float32(100)
    return to_matlab(np.maximum(a, np.float32(0)))


def class_data(msize, seed: int, cls: str) -> np.ndarray:
    """the raw beads in a class deconvlucy takes; MATLAB's form"""
    u = raw_beads(msize, seed)
    if cls == "uint16":
        out = u
    elif cls == "uint8":
        out = np.floor(u.astype(np.float64) / 8 + 0.5).astype(np.uint8)  # MATLAB's uint8(u / 8)
    elif cls == "int16":
        out = u.astype(np.int16)  # < 2^15
    elif cls == "double":
        out = u.astype(np.float64) + 0.25
    else:
        out = (u.astype(np.float32) + np.float32(0.37)).astype(np.float32)
    return to_matlab(out)


def make_psf(msize, cls: str = "single", factor: float = 1.0) -> np.ndarray:
    """the test PSF cropped to the data as RLdecon crops it (a plane for 2D data), divided by its
    sum (correctly rounded double sum) as RLdecon divides it, times ``factor``; MATLAB's form"""
    psf = decon_psf((31, 31, 31), (2.0, 1.5))
    py = _py_size(msize)
    if len(py) == 2:
        psf = psf[15]
    s = [max(0, (p - n) // 2) for p, n in zip(psf.shape, py, strict=True)]
    t = [a + min(p, n) for a, p, n in zip(s, psf.shape, py, strict=True)]
    psf = np.ascontiguousarray(psf[tuple(slice(a, b) for a, b in zip(s, t, strict=True))])
    total = math.fsum(psf.astype(np.float64).ravel().tolist())
    d = psf.astype(np.float64) / total * factor
    out = d.astype(np.float32) if cls == "single" else d
    return to_matlab(np.ascontiguousarray(out))


def asym_psf(msize) -> np.ndarray:
    """a PSF without symmetry (its OTF is complex): the test PSF plus a shifted, narrower one"""
    a = from_m(make_psf(msize, "double"))
    b = np.roll(a, (1, 2), axis=(-2, -1)) * 0.5
    c = a + b
    total = math.fsum(c.ravel().tolist())
    return to_matlab(np.ascontiguousarray((c / total).astype(np.float32)))


def from_m(a):
    from petakit5d.core.arrays import from_matlab

    return from_matlab(a)


def weight_mask(msize) -> np.ndarray:
    """WEIGHT: ones, zeros in the first 5 rows (MATLAB's y) and 0.5 in a block; MATLAB's form"""
    w = np.ones(_py_size(msize), np.float64)
    w[..., :5, :] = 0
    w[..., 10:20, 8:14] = 0.5
    return to_matlab(w)


# ---------------------------------------------------------------------------------------------
# Python sides
# ---------------------------------------------------------------------------------------------


def py_lucy(*args, backend):
    from petakit5d.decon import other_methods as O

    return O.deconvlucy(*args, backend=backend)


# MATLAB: deconvlucy with the arguments, [] where a default is meant; an error gives its id
_TRY = "try, J = deconvlucy(in{:}); out = {'', J}; catch e, out = {e.identifier, []}; end"


def py_try(*args, backend):
    from petakit5d.decon import other_methods as O

    try:
        return ["", O.deconvlucy(*args, backend=backend)]
    except O.DeconvlucyError as e:
        return [e.identifier, np.zeros((0, 0))]


_CELL = (
    "J = deconvlucy({in{1}}, in{2}, in{3}); K = deconvlucy(J, in{2}, in{4}); "
    "L = deconvlucy({in{1}}, in{2}, in{3} + in{4}); "
    "out = {J{1}, J{2}, J{3}, J{4}, K{2}, K{3}, K{4}, L{2}};"
)


def py_cell(im, psf, n1, n2, *, backend):
    from petakit5d.decon import other_methods as O

    def cols(g):  # MATLAB's numel x 2 array of the changes
        return np.stack(
            [to_matlab(g[0]).reshape(-1, order="F"), to_matlab(g[1]).reshape(-1, order="F")], axis=1
        )

    j = O.deconvlucy_state(im, psf, n1, backend=backend)
    k = O.deconvlucy_state(j, psf, n2, backend=backend)
    m = O.deconvlucy_state(im, psf, n1 + n2, backend=backend)
    return [
        j.image,
        j.estimate,
        j.previous,
        cols(j.internal),
        k.estimate,
        k.previous,
        cols(k.internal),
        m.estimate,
    ]


_OTF = "o = psf2otf(in{1}, in{2}); out = {o, isreal(o)};"

_RLDECON = "r = RLdecon('', '', in{1}, 0.108, 0.2, 0.2, 'rawdata', in{2}, in{3 : end}); out = {r};"


def py_rldecon(psf_fn, im, *opts, backend):
    from petakit5d import matlab

    with P.using(backend):
        r = matlab.RLdecon("", "", psf_fn, 0.108, 0.2, 0.2, "rawdata", im, *opts)
    return r[0]


def _psf_file(path):
    p = np.floor(decon_psf((31, 31, 31), (2.0, 1.5)) + np.float32(0.5)).astype(np.uint16)
    write_tiff(path, p)


# ---------------------------------------------------------------------------------------------
# cases
# ---------------------------------------------------------------------------------------------


def _spec(kind: str, backend: str):
    """the comparison of an output: 'exact' (texts, flags, errors), else by the result's class"""
    if kind == "exact":
        return None
    if backend == "exact":
        return DOUBLE_TOL if kind == "double" else None
    return PERMISSIVE_INT if kind == "int" else PERMISSIVE_TOL


def cases(
    cid,
    *,
    inputs,
    args=None,
    code=None,
    python,
    nargout=1,
    out_roles=("vol",),
    kinds=("single",),
    backends=BACKENDS,
):
    res = []
    for be in backends:
        specs = {i: _spec(k, be) for i, k in enumerate(kinds)}
        compare = {i: v for i, v in specs.items() if v is not None}
        res.append(
            Case(
                f"{cid}.{be}",
                matlab=None if code else "deconvlucy",
                code=code,
                args=args,
                nargout=nargout,
                python=python,
                py_kwargs={"backend": be},
                inputs=inputs,
                out_roles=out_roles,
                py_error=Exception,
                files=False,
                **({"compare": compare} if compare else {}),
            )
        )
    return res


def V(make):
    return Value(make)


CASES: list[Case] = []

# RLdecon's call: sizes and iterations
_MS = [48, 46, 36]
SIZES = [[47, 45, 33], [40, 38, 25], [64, 50], [31, 45], [10, 64, 23]]
for _n in (1, 2, 3, 5, 15):
    CASES += cases(
        f"lucy.48x46x36.iter{_n}",
        inputs={"im": V(lambda: rldecon_data(_MS, 1)), "psf": V(lambda: make_psf(_MS))},
        args=[Vol(In("im")), Vol(In("psf")), _n],
        python=py_lucy,
    )
for _i, _ms in enumerate(SIZES):
    CASES += cases(
        f"lucy.{'x'.join(map(str, _ms))}.iter10",
        inputs={
            "im": V(lambda m=_ms, s=_i: rldecon_data(m, 10 + s)),
            "psf": V(lambda m=_ms: make_psf(m)),
        },
        args=[Vol(In("im")), Vol(In("psf")), 10],
        python=py_lucy,
    )
CASES += cases(
    "lucy.asym",
    inputs={"im": V(lambda: rldecon_data(_MS, 2)), "psf": V(lambda: asym_psf(_MS))},
    args=[Vol(In("im")), Vol(In("psf")), 5],
    python=py_lucy,
)
CASES += cases(
    "lucy.unnormalized",
    inputs={"im": V(lambda: rldecon_data(_MS, 3)), "psf": V(lambda: make_psf(_MS, factor=3.0))},
    args=[Vol(In("im")), Vol(In("psf")), 5],
    python=py_lucy,
)

# classes
_SM = [32, 30, 20]
for _cls in ("uint8", "uint16", "int16", "single", "double"):
    CASES += cases(
        f"class.{_cls}",
        inputs={"im": V(lambda c=_cls: class_data(_SM, 20, c)), "psf": V(lambda: make_psf(_SM))},
        args=[Vol(In("im")), Vol(In("psf")), 5],
        python=py_lucy,
        kinds=(_cls if _cls in ("single", "double") else "int",),
    )
CASES += cases(
    "class.double_psf_single",
    inputs={
        "im": V(lambda: rldecon_data(_SM, 21)),
        "psf": V(lambda: make_psf(_SM, "double")),
    },
    args=[Vol(In("im")), Vol(In("psf")), 5],
    python=py_lucy,
)
CASES += cases(
    "class.default_numit",
    inputs={"im": V(lambda: rldecon_data(_SM, 22)), "psf": V(lambda: make_psf(_SM))},
    args=[Vol(In("im")), Vol(In("psf"))],
    python=py_lucy,
)

# documented options and errors (MATLAB's error identifiers)
_OPT_IN = {
    "im": V(lambda: rldecon_data(_SM, 23)),
    "psf": V(lambda: make_psf(_SM)),
    "w": V(lambda: weight_mask(_SM)),
}
_E = np.zeros((0, 0))
OPTS = {
    "dampar": [5, 5.0],
    "dampar_single": [5, np.float32(5.0)],
    "weight": [4, _E, Vol(In("w"))],
    "readout": [4, _E, _E, 3.0],
    "all": [4, 5.0, Vol(In("w")), 3.0],
    "subsmpl1": [3, _E, _E, _E, 1],
    "numit0": [0],
    "numit_frac": [2.5],
}
for _k, _extra in OPTS.items():
    CASES += cases(
        f"opt.{_k}",
        inputs=_OPT_IN,
        code=_TRY,
        args=[Vol(In("im")), Vol(In("psf")), *_extra],
        python=py_try,
        nargout=2,
        out_roles=(None, "vol"),
        kinds=("exact", "single"),
    )
CASES += cases(
    "opt.zeros",
    inputs={
        "im": V(lambda: to_matlab(np.zeros((20, 32, 30), np.float32))),
        "psf": V(lambda: make_psf(_SM)),
    },
    code=_TRY,
    args=[Vol(In("im")), Vol(In("psf")), 3],
    python=py_try,
    nargout=2,
    out_roles=(None, "vol"),
    kinds=("exact", "single"),
)
CASES += cases(
    "opt.negative",
    inputs={
        "im": V(lambda: to_matlab(from_m(rldecon_data(_SM, 24)) - np.float32(50))),
        "psf": V(lambda: make_psf(_SM)),
    },
    code=_TRY,
    args=[Vol(In("im")), Vol(In("psf")), 3],
    python=py_try,
    nargout=2,
    out_roles=(None, "vol"),
    kinds=("exact", "single"),
)
CASES += cases(
    "opt.psf_larger",
    inputs={"im": V(lambda: rldecon_data([20, 22, 12], 25)), "psf": V(lambda: make_psf(_SM))},
    code=_TRY,
    args=[Vol(In("im")), Vol(In("psf")), 3],
    python=py_try,
    nargout=2,
    out_roles=(None, "vol"),
    kinds=("exact", "single"),
)
_SUB = [24, 22]
CASES += cases(
    "opt.subsmpl2",
    inputs={"im": V(lambda: rldecon_data(_SUB, 26)), "psf": V(lambda: make_psf([15, 15]))},
    code=_TRY,
    args=[Vol(In("im")), Vol(In("psf")), 3, _E, _E, _E, 2],
    python=py_try,
    nargout=2,
    out_roles=(None, "vol"),
    kinds=("exact", "single"),
)


def opt_cases(cid, inputs, extra, kind="single", backends=BACKENDS):
    return cases(
        cid,
        inputs=inputs,
        code=_TRY,
        args=[Vol(In("im")), Vol(In("psf")), *extra],
        python=py_try,
        nargout=2,
        out_roles=(None, "vol"),
        kinds=("exact", kind),
        backends=backends,
    )


_S = np.float32
OPTS2 = {
    "readout_single": [4, _E, _E, _S(3.0)],
    "readout_frac": [4, _E, _E, _S(2.7)],
    "dampar_frac": [4, _S(3.3)],
    "all_single": [4, _S(3.3), Vol(In("w")), _S(2.7)],
    "weight_single": [3, _E, Vol(In("ws"))],
    "weight_size": [3, _E, Vol(In("wbad"))],
    "numit_empty": [_E],
    "dampar_negative": [3, _S(-1.0)],
    "subsmpl0": [3, _E, _E, _E, 0],
}
_OPT_IN2 = dict(
    _OPT_IN,
    ws=V(lambda: weight_mask(_SM).astype(np.float32)),
    wbad=V(lambda: weight_mask([32, 30, 19])),
)
for _k, _extra in OPTS2.items():
    CASES += opt_cases(f"opt.{_k}", _OPT_IN2, _extra)
CASES += opt_cases(
    "opt.double_damp_readout",
    {"im": V(lambda: class_data(_SM, 28, "double")), "psf": V(lambda: make_psf(_SM))},
    [4, 3.3, _E, 1.5],
    "double",
)
CASES += opt_cases(
    "opt.uint16_damp_readout",
    {"im": V(lambda: class_data(_SM, 29, "uint16")), "psf": V(lambda: make_psf(_SM))},
    [4, np.uint16(20), _E, np.uint16(5)],
    "int",
)
CASES += opt_cases(
    "opt.psf_zero",
    {"im": V(lambda: rldecon_data(_SM, 30)), "psf": V(lambda: np.zeros((9, 9, 9), np.float32))},
    [3],
)
CASES += opt_cases(
    "opt.psf_negative",
    {
        "im": V(lambda: rldecon_data(_SM, 31)),
        "psf": V(lambda: make_psf(_SM) - np.float32(1e-4)),
    },
    [3],
    # negative PSF values amplify FFT round-off (relative L2 3e-4 on the permissive backend)
    backends=("exact",),
)
CASES += opt_cases(
    "opt.psf2d_data3d",
    {"im": V(lambda: rldecon_data(_SM, 32)), "psf": V(lambda: make_psf([15, 15]))},
    [3],
)
CASES += opt_cases(
    "opt.subsmpl2_3d",
    {
        "im": V(lambda: rldecon_data([16, 14, 10], 33)),
        "psf": V(lambda: make_psf([9, 9, 9])),
    },
    [3, _E, _E, _E, 2],
)
CASES += opt_cases(
    "opt.subsmpl3",
    {"im": V(lambda: rldecon_data([12, 10], 34)), "psf": V(lambda: make_psf([15, 15]))},
    [3, _E, _E, _E, 3],
)
CASES += opt_cases(
    "opt.subsmpl2_damp_weight",
    {
        "im": V(lambda: rldecon_data(_SUB, 35)),
        "psf": V(lambda: make_psf([15, 15])),
        "w": V(lambda: weight_mask(_SUB)),
    },
    [4, _S(3.3), Vol(In("w")), _S(2.7), 2],
)

# the cell-array form, and a run resumed from it
for _ms in ([64, 50], [32, 30, 20]):
    CASES += cases(
        f"cell.{'x'.join(map(str, _ms))}",
        inputs={"im": V(lambda m=_ms: rldecon_data(m, 27)), "psf": V(lambda m=_ms: make_psf(m))},
        code=_CELL,
        args=[Vol(In("im")), Vol(In("psf")), 3, 2],
        python=py_cell,
        nargout=8,
        out_roles=("vol", "vol", "vol", None, "vol", "vol", None, "vol"),
        kinds=("single", *("double",) * 7),
    )

# psf2otf (deconvlucy's OTF): real for a symmetric PSF?
for _name, _mk in (
    ("sym", lambda: make_psf(_MS)),
    ("asym", lambda: asym_psf(_MS)),
    ("sym_double", lambda: make_psf(_MS, "double")),
):
    CASES.append(
        Case(
            f"otf.{_name}",
            code=_OTF,
            args=[In("psf"), [48, 46, 36]],
            nargout=2,
            python=None,
            inputs={"psf": V(_mk)},
            files=False,
        )
    )

# where psf2otf's result turns complex: a symmetric PSF plus d * max(PSF) next to its center
_OTF_T = (
    "p = in{1}; d = in{3}; c = num2cell(floor(size(p) / 2) + 1); c{1} = c{1} + 1; "
    "k = sub2ind(size(p), c{:}); m = max(p(:)); f = false(1, numel(d)); "
    "for i = 1 : numel(d), q = p; q(k) = q(k) + d(i) * m; f(i) = isreal(psf2otf(q, in{2})); end; "
    "out = {f};"
)
DELTAS = [float(f"{m}e{e}") for e in range(-17, -8) for m in (1, 2, 5)]


def py_otf_flags(psf, msize, deltas, *, backend):
    from petakit5d.core.arrays import from_matlab, size_from_matlab
    from petakit5d.decon import other_methods as O

    p = from_matlab(np.asarray(psf)).astype(np.float64)
    c = [n // 2 for n in p.shape]
    c[-2] += 1
    m = float(p.max())
    flags = []
    for d in deltas:
        q = p.copy()
        q[tuple(c)] += d * m
        flags.append(np.isrealobj(O.psf2otf_real(q, size_from_matlab(msize), backend=backend)))
    return np.array([flags])


for _ms in ([48, 46, 36], [64, 50], [96, 80, 40]):
    CASES += cases(
        f"otf.threshold.{'x'.join(map(str, _ms))}",
        inputs={"psf": V(lambda m=_ms: make_psf(m, "double"))},
        code=_OTF_T,
        args=[In("psf"), _ms, DELTAS],
        python=py_otf_flags,
        out_roles=(None,),
        kinds=("exact",),
    )
# the bound, closer: relative to max|OTF| (a PSF of sum 1000), by the PSF's number of elements
for _name, _ms, _mk, _ds in (
    ("2d", [64, 50], lambda: make_psf([64, 50], "double"), range(200, 401, 5)),
    ("2d_x1000", [64, 50], lambda: make_psf([64, 50], "double", 1000.0), range(200, 401, 5)),
    ("2d_psf15", [64, 50], lambda: make_psf([15, 15], "double"), range(10, 201, 5)),
    ("3d", [48, 46, 36], lambda: make_psf([48, 46, 36], "double"), range(600, 801, 5)),
):
    _sc = "e-13" if _name != "3d" else "e-11"
    CASES += cases(
        f"otf.bound.{_name}",
        inputs={"psf": V(_mk)},
        code=_OTF_T,
        args=[In("psf"), _ms, [float(f"{k}{_sc}") for k in _ds]],
        python=py_otf_flags,
        out_roles=(None,),
        kinds=("exact",),
    )

# RLdecon with RLMethod 'original' in memory
_RLD = ["RLMethod", "original", "useGPU", False, "Background", 100, "DeconIter", 5]
for _k, _extra in {
    "single": ["psfGen", False, "save16bit", False],
    "psfgen": ["save16bit", False],
    "bbox_scale": ["psfGen", False, "deconBbox", [3, 4, 2, 40, 41, 30], "scaleFactor", 2],
}.items():
    CASES += cases(
        f"rldecon.{_k}",
        inputs={
            "root": File(_psf_file, "psf.tif"),
            "im": V(lambda: to_matlab(raw_beads(_MS, 30))),
        },
        code=_RLDECON,
        args=[InPath("root"), In("im"), *_RLD, *_extra],  # the MATLAB-name layer maps the data
        python=py_rldecon,
        out_roles=(None,),  # the MATLAB-name layer returns MATLAB's form
    )


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    if case.id.endswith(".exact"):
        reason = digests_hold_here()
        if reason:
            pytest.skip(reason)
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# the port's own rules
# ---------------------------------------------------------------------------------------------


def _vol(shape, seed):
    a = raw_beads((shape[1], shape[2], shape[0]), seed).astype(np.float32) - np.float32(100)
    return np.maximum(a, np.float32(0))


def _psf(shape):
    return from_m(make_psf((shape[1], shape[2], shape[0])))


def test_results_do_not_depend_on_the_thread_count():
    from petakit5d.decon import other_methods as O

    im, psf = _vol((40, 256, 256), 40), _psf((40, 256, 256))
    ref = O.deconvlucy(im, psf, 3, threads=1)
    for t in (3, 8):
        assert np.array_equal(O.deconvlucy(im, psf, 3, threads=t), ref), t


def test_a_resumed_run_continues_the_iterations():
    from petakit5d.decon import other_methods as O

    im, psf = _vol((20, 32, 30), 41), _psf((20, 32, 30))
    w = from_m(weight_mask([32, 30, 20]))
    a = O.deconvlucy_state(im, psf, 3, np.float32(2.0), w)
    b = O.deconvlucy_state(a, psf, 2, np.float32(2.0), w)
    c = O.deconvlucy_state(im, psf, 5, np.float32(2.0), w)
    assert np.array_equal(b.image, im)
    for x, y in zip(b[1:], c[1:], strict=True):
        assert np.array_equal(x, y)
    assert b.internal.shape == (2, 20, 32, 30) and b.estimate.dtype == np.float64


def test_the_data_is_not_changed():
    from petakit5d.decon import other_methods as O

    im = _vol((20, 32, 30), 44).astype(np.float64) - 50.0
    keep = im.copy()
    st = O.deconvlucy_state(im, _psf((20, 32, 30)), 2)
    assert np.array_equal(im, keep) and st.image is im
    assert not np.shares_memory(st.estimate, im) and not np.shares_memory(st.previous, im)


def test_subsampling_returns_the_fine_grid():
    from petakit5d.decon import other_methods as O

    im = _vol((10, 16, 14), 42)
    psf = from_m(make_psf([9, 9, 9]))
    st = O.deconvlucy_state(im, psf, 2, subsmpl=2)
    assert st.estimate.shape == (20, 32, 28) and st.internal.shape == (2, 20, 32, 28)
    assert O.deconvlucy(im, psf, 2, subsmpl=2).dtype == np.float32


_ONES = np.ones((4, 8, 8), np.float32)


@pytest.mark.parametrize(
    ("im", "kw", "identifier"),
    [
        (np.zeros((4, 8, 8), np.int32), {}, "MATLAB:deconvlucy:invalidType"),
        (_ONES, {"numit": 0}, "MATLAB:deconvlucy:expectedPositive"),
        (_ONES, {"numit": -1.5}, "MATLAB:deconvlucy:expectedPositive"),
        (_ONES, {"subsmpl": 1.5}, "MATLAB:deconvlucy:expectedInteger"),
        (_ONES, {"subsmpl": 0}, "MATLAB:deconvlucy:expectedPositive"),
        (_ONES, {"dampar": 2.0}, "images:deconvlucy:damparMustBeSameClassAsInputImage"),
        (_ONES, {"readout": 1.0}, "images:deconvlucy:readoutMustBeSameClassAsInputImage"),
        (_ONES, {"weight": _ONES}, "MATLAB:deconvlucy:invalidType"),
        (_ONES, {"weight": np.ones((4, 8, 7))}, "images:deconvlucy:weightMustBeSameSizeAsImage"),
    ],
)
def test_arguments_are_checked_as_matlab_checks_them(im, kw, identifier):
    from petakit5d.decon import other_methods as O

    with pytest.raises(O.DeconvlucyError) as e:
        O.deconvlucy(im, np.ones((3, 3, 3), np.float32), **{"numit": 2, **kw})
    assert e.value.identifier == identifier


def test_psf_rules():
    from petakit5d.decon import other_methods as O

    with pytest.raises(O.DeconvlucyError) as e:
        O.deconvlucy(_ONES, np.zeros((3, 3, 3)), 2)
    assert e.value.identifier == "images:deconvlucy:psfMustNotBeZeroEverywhere"
    with pytest.raises(O.DeconvlucyError) as e:
        O.deconvlucy(_ONES, np.ones((5, 3, 3)), 2)
    assert e.value.identifier == "images:deconvlucy:psfMustBeSmallerThanImage"
    # with subsampling the PSF may fill the fine grid; a 2D PSF serves 3D data
    assert O.deconvlucy(_ONES, np.ones((8, 16, 16)), 1, subsmpl=2).shape == (8, 16, 16)
    assert O.deconvlucy(_ONES, np.ones((3, 3)), 1).shape == (4, 8, 8)


def test_psf2otf_real_drops_round_off_only():
    from petakit5d.decon import other_methods as O

    p = _psf((20, 32, 30)).astype(np.float64)
    assert np.isrealobj(O.psf2otf_real(p, (20, 32, 30)))
    q = p.copy()
    q[10, 16, 14] += 1e-3
    assert np.iscomplexobj(O.psf2otf_real(q, (20, 32, 30)))


def test_rldecon_branch_crops_and_scales():
    from petakit5d.core import dtypes
    from petakit5d.decon import other_methods as O

    im, psf = _vol((20, 32, 30), 43), _psf((20, 32, 30))
    box = (1, 2, 3, 18, 30, 25)
    out = O.rl_decon_original(im, psf, 3, decon_bbox=box, scale_factor=1.5)
    ref = np.asarray(dtypes.times(O.deconvlucy(im, psf, 3)[1:18, 2:30, 3:25], 1.5))
    assert out.dtype == np.float32 and np.array_equal(out, ref)
