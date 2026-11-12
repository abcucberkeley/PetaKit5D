"""petakit5d.matlabcompat.warp against MATLAB's imwarp (3D affine), and its own tests.

Oracle cases: one family per transform shape PetaKit5D's code gives imwarp (deskew: a shear of
x by z; rotation in the x-z plane with z scaled, sample scan and objective scan; deskew and
rotation combined; resampling; output boxes; PSF rotation in double; generic matrices), over
classes, methods and sizes, plus edge and error cases. The MATLAB side calls
``imwarp(V, affine3d(T), method, 'FillValues', fill, 'OutputView', R)`` with ``R`` the
``imref3d`` of an output size or of a box (as the lab's deskewRotateFrame3D builds it).
Transforms are built here from the math in MATLAB's row-vector convention on (x, y, z).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from petakit5d.core import arrays, bbox
from petakit5d.matlabcompat import warp
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, In, Value
from petakit5d.testing.synth import CounterRNG, stable_sincos_2pi, texture

pytestmark = pytest.mark.skipif(
    not __import__("petakit5d._ext", fromlist=["available"]).available("warp"),
    reason="extension module warp not built",
)

_from_matlab = arrays.from_matlab
_to_matlab = arrays.to_matlab


def _same_bits(a, b) -> bool:
    a, b = np.asarray(a), np.asarray(b)
    if a.dtype != b.dtype or a.shape != b.shape:
        return False
    return bool(
        np.all(np.ascontiguousarray(a).view(np.uint8) == np.ascontiguousarray(b).view(np.uint8))
    )


# ------------------------------------------------------------------------------------------
# transforms (MATLAB's convention: [x y z 1] * T)
# ------------------------------------------------------------------------------------------


def _mul(*ms):
    """matrix product, left to right, each element summed in index order"""
    out = ms[0]
    for b in ms[1:]:
        out = [
            [
                ((out[i][0] * b[0][j] + out[i][1] * b[1][j]) + out[i][2] * b[2][j])
                + out[i][3] * b[3][j]
                for j in range(4)
            ]
            for i in range(4)
        ]
    return out


def _translate(x, y, z):
    return [[1.0, 0, 0, 0], [0, 1.0, 0, 0], [0, 0, 1.0, 0], [float(x), float(y), float(z), 1.0]]


def _scale(x, y, z):
    return [[float(x), 0, 0, 0], [0, float(y), 0, 0], [0, 0, float(z), 0], [0, 0, 0, 1.0]]


def _shear(s, tx):
    """x' = x + s * z + tx (a sample-scan deskew: plane z moves by s along x)"""
    return [[1.0, 0, 0, 0], [0, 1.0, 0, 0], [float(s), 0, 1.0, 0], [float(tx), 0, 0, 1.0]]


def _sincos(deg):
    s, c = stable_sincos_2pi((deg % 360.0) / 360.0)
    return float(s), float(c)


def _rot_xz(deg):
    """rotation in the x-z plane by deg degrees"""
    s, c = _sincos(deg)
    return [[c, 0, -s, 0], [0, 1.0, 0, 0], [s, 0, c, 0], [0, 0, 0, 1.0]]


def _center(n_yxz):
    """centers (x, y, z) of a volume [ny nx nz]: (n + 1) / 2"""
    ny, nx, nz = n_yxz
    return (nx + 1) / 2, (ny + 1) / 2, (nz + 1) / 2


def _rotation_chain(in_yxz, out_yxz, deg, zscale):
    """move the input center to the origin, scale z, rotate in x-z, move to the output center"""
    ci, co = _center(in_yxz), _center(out_yxz)
    return _mul(
        _translate(-ci[0], -ci[1], -ci[2]), _scale(1, 1, zscale), _rot_xz(deg), _translate(*co)
    )


def _rotated_size(ny, nx, nz, deg, zscale):
    s, c = _sincos(deg)
    w = nx * abs(c) + nz * zscale * abs(s)
    h = nx * abs(s) + nz * zscale * abs(c)
    return [ny, max(1, int(w + 0.5)), max(1, int(h + 0.5))]


def _resample(out_yxz, r_xyz):
    """resampling about the center by factors r (x, y, z): the output size shrinks by r"""
    new = [
        max(1, int(out_yxz[0] / r_xyz[1] + 0.5)),
        max(1, int(out_yxz[1] / r_xyz[0] + 0.5)),
        max(1, int(out_yxz[2] / r_xyz[2] + 0.5)),
    ]
    c0, c1 = _center(out_yxz), _center(new)
    m = _mul(
        _translate(-c0[0], -c0[1], -c0[2]),
        _scale(1 / r_xyz[0], 1 / r_xyz[1], 1 / r_xyz[2]),
        _translate(*c1),
    )
    return m, new


def _deskew_shift(dz, xy, deg):
    s, c = _sincos(deg)
    return c * dz / xy


# ------------------------------------------------------------------------------------------
# inputs
# ------------------------------------------------------------------------------------------

_CLS = {
    "uint8": np.uint8,
    "int16": np.int16,
    "uint16": np.uint16,
    "single": np.float32,
    "double": np.float64,
}


def _volume(size_yxz, cls, seed, kind="texture"):
    """MATLAB's form: (y, x, z), Fortran order"""

    def make():
        rng = CounterRNG(seed, f"warp/{cls}/{kind}/{size_yxz}")
        if kind == "texture":
            t = texture(size_yxz, 2.0, (100.0, 4000.0), rng).astype(np.float64)
            if cls == "single":
                v = t.astype(np.float32)
            elif cls == "double":
                v = t / 1000.0
            elif cls == "uint8":
                v = (t / 16.0).astype(np.uint8)
            elif cls == "int16":
                v = (t * 8.0 - 16000.0).astype(np.int16)
            else:
                v = (t * 16.0).astype(np.uint16)
        else:
            if cls in ("single", "double"):
                v = (rng.normal(size_yxz) * 1000.0).astype(_CLS[cls])
            else:
                info = np.iinfo(_CLS[cls])
                v = rng.integers(info.min, info.max, size_yxz, dtype=_CLS[cls])
        return np.asfortranarray(v)

    return Value(make)


MATLAB_CODE = (
    "w = in{3}; "
    "if numel(w) == 3, R = imref3d(w, 1, 1, 1); "
    "else, R = imref3d(w(4:6) - w(1:3) + 1, [w(2) - 0.5, w(5) + 0.5], [w(1) - 0.5, w(4) + 0.5], "
    "[w(3) - 0.5, w(6) + 0.5]); end; "
    "out = {imwarp(in{1}, affine3d(in{2}), in{4}, 'FillValues', in{5}, 'OutputView', R)};"
)


def py_imwarp(v, T, view, method, fill):
    """the Python side in MATLAB's terms: v (y, x, z), view [ny nx nz] or a MATLAB bbox (until
    F14's roles: Vol, PerAxis('yxz') or Box, out_roles=('vol',))"""
    view = [int(x) for x in np.ravel(view)]
    if len(view) == 3:
        kw = {"output_size": arrays.size_from_matlab(view)}
    else:
        kw = {"output_box": bbox.box_from_matlab(view)}
    out = warp.imwarp(
        _from_matlab(v), np.asarray(T, dtype=np.float64), str(method), fill_value=float(fill), **kw
    )
    return _to_matlab(out)


def _case(cid, size_yxz, cls, T, view, method="linear", fill=0.0, kind="texture", seed=1, **kw):
    return Case(
        cid,
        code=MATLAB_CODE,
        args=[
            In("v"),
            np.array(T, dtype=np.float64),
            [float(x) for x in view],
            method,
            float(fill),
        ],
        python=py_imwarp,
        inputs={"v": _volume(size_yxz, cls, seed, kind)},
        **kw,
    )


def _families():
    cases = []
    # deskew (shear of x by z), with and without a crop of the width, both directions
    for n, (size, deg, dz, xy, rev, crop) in enumerate(
        [
            ([20, 36, 17], 31.5, 0.4, 0.108, False, False),
            ([17, 29, 14], 32.45, 0.3, 0.104, True, False),
            ([23, 31, 12], 45.0, 0.25, 0.147, False, True),
            ([16, 40, 9], 32.8, 0.5, 0.1, True, True),
            ([19, 22, 11], 60.0, 0.2, 0.1, False, False),
        ]
    ):
        ny, nx, nz = size
        s = _deskew_shift(dz, xy, deg)
        width = math.ceil((nz - 1) * s)
        tx, step = (-s, s) if not rev else (s + width, -s)
        out_nx = nx + width
        if crop:
            cut = min(nx // 2, width // 2)
            tx -= cut
            out_nx -= 2 * cut
        T = _shear(step, tx)
        for cls, methods in (
            ("single", ("linear", "nearest", "cubic")),
            ("uint16", ("linear", "cubic")),
        ):
            for m in methods:
                cases.append(
                    _case(
                        f"deskew{n}_{cls}_{m}",
                        size,
                        cls,
                        T,
                        [ny, out_nx, nz],
                        m,
                        kind="texture" if n % 2 else "random",
                    )
                )
    # an exact integer shear and a dyadic one: coordinates on the grid
    for n, s in enumerate((2.0, 0.25)):
        T = _shear(s, -s)
        for cls in ("single", "uint16"):
            cases.append(
                _case(
                    f"deskew_exact{n}_{cls}", [12, 20, 10], cls, T, [12, 20 + math.ceil(9 * s), 10]
                )
            )
    # rotation (sample scan: z scaled by sin(angle) dz / xy; objective scan: dz / xy)
    for n, (size, deg, zs) in enumerate(
        [
            ([18, 60, 14], 31.5, abs(_sincos(31.5)[0]) * 0.4 / 0.108),
            ([15, 52, 12], -32.45, 0.3 / 0.104),
            ([21, 44, 16], 45.0, 1.0),
            ([14, 70, 10], 32.8, 2.2),
            ([17, 38, 20], 90.0 - 32.45, 1.6),
        ]
    ):
        out = _rotated_size(*size, deg, zs)
        T = _rotation_chain(size, out, deg, zs)
        for cls, methods in (
            ("single", ("linear", "nearest", "cubic")),
            ("uint16", ("linear", "nearest", "cubic")),
        ):
            for m in methods:
                cases.append(
                    _case(
                        f"rotate{n}_{cls}_{m}",
                        size,
                        cls,
                        T,
                        out,
                        m,
                        kind="random" if n % 2 else "texture",
                        seed=2,
                    )
                )
    # deskew and rotation combined, then resampled, sample scan
    for n, (size, deg, dz, xy, rev, r) in enumerate(
        [
            ([16, 30, 14], 31.5, 0.4, 0.108, False, None),
            ([13, 26, 12], 32.45, 0.3, 0.104, True, None),
            ([15, 28, 13], 31.5, 0.5, 0.108, False, (2.0, 2.0, 2.0)),
            ([14, 24, 15], 32.45, 0.3, 0.104, True, (1.5, 1.5, 1.5)),
            ([12, 30, 11], 45.0, 0.25, 0.147, False, (1.0, 1.0, 2.0)),
            ([10, 18, 9], 31.5, 0.4, 0.108, False, (0.5, 0.5, 0.5)),
        ]
    ):
        ny, nx, nz = size
        s = _deskew_shift(dz, xy, deg)
        width = math.ceil((nz - 1) * s)
        step, tx = (s, -s) if not rev else (-s, s + width)
        sn, _ = _sincos(deg)
        zs = abs(sn) * dz / xy
        dsk = [ny, nx + width, nz]
        out = _rotated_size(ny, nx, nz, deg, zs)
        T = _mul(_shear(step, tx), _rotation_chain(dsk, out, -deg if rev else deg, zs))
        if r is not None:
            R, out = _resample(out, r)
            T = _mul(T, R)
        for cls, methods in (
            ("single", ("linear", "cubic")),
            ("uint16", ("linear", "nearest", "cubic")),
        ):
            for m in methods:
                cases.append(
                    _case(
                        f"dsr{n}_{cls}_{m}",
                        size,
                        cls,
                        T,
                        out,
                        m,
                        kind="texture" if n % 2 else "random",
                        seed=3,
                    )
                )
        # output boxes of the same transform: inside, and reaching outside the full view
        oy, ox, oz = out
        for b, box in enumerate(
            [
                [2, 3, 2, oy - 1, ox - 2, oz - 1],
                [1, ox // 3, oz // 4, oy, ox // 3 + 7, oz // 4 + 5],
                [oy - 3, ox // 2 - 3, oz // 2 - 2, oy + 2, ox // 2 + 3, oz // 2 + 3],
            ]
        ):
            cases.append(
                _case(f"dsr{n}_box{b}_uint16", size, "uint16", T, box, kind="random", seed=3)
            )
            if b == 1:
                cases.append(
                    _case(f"dsr{n}_box{b}_single_cubic", size, "single", T, box, "cubic", seed=3)
                )
    # objective scan combined: rotation of the raw volume, resampled
    for n, (size, deg, zs, r) in enumerate(
        [
            ([14, 40, 12], 31.5, 0.4 / 0.108, None),
            ([12, 34, 14], 32.45, 0.3 / 0.104, (1.5, 1.5, 1.5)),
        ]
    ):
        out = _rotated_size(*size, deg, zs)
        T = _rotation_chain(size, out, deg, zs)
        if r is not None:
            R, out = _resample(out, r)
            T = _mul(T, R)
        for cls in ("single", "uint16"):
            for m in ("linear", "cubic"):
                cases.append(_case(f"dsr_os{n}_{cls}_{m}", size, cls, T, out, m, seed=4))
    # resampling alone (scale about the centers)
    for n, (size, r) in enumerate(
        [
            ([20, 24, 16], (2.0, 2.0, 2.0)),
            ([15, 18, 12], (0.75, 0.75, 1.5)),
            ([9, 11, 7], (1 / 3, 1 / 3, 1 / 3)),
        ]
    ):
        R, out = _resample(size, r)
        for cls in ("single", "uint16"):
            for m in ("linear", "nearest", "cubic"):
                cases.append(
                    _case(f"scale{n}_{cls}_{m}", size, cls, R, out, m, kind="random", seed=5)
                )
    # PSF rotation (double; the decon PSF path) and rotation after decon (single)
    for n, (size, deg, zs) in enumerate(
        [([31, 31, 41], 31.5, 0.25 / 0.108), ([25, 25, 33], -32.45, 0.2 / 0.104)]
    ):
        out = _rotated_size(*size, deg, zs)
        T = _rotation_chain(size, out, deg, zs)
        cases.append(_case(f"psf{n}_double_linear", size, "double", T, out, seed=6))
        cases.append(_case(f"psf{n}_double_cubic", size, "double", T, out, "cubic", seed=6))
    # generic affine matrices about the centers (the full LU path, y mixed in), every class
    rng = CounterRNG(7, "warp/generic")
    ci, co = _center([13, 15, 11]), _center([12, 16, 13])
    for n in range(6):
        lin = np.eye(4)
        lin[:3, :3] += rng.integers(-24, 24, (3, 3)) / 64.0
        shift = rng.integers(-2048, 2048, 3) / 1024.0
        T = _mul(
            _translate(-ci[0], -ci[1], -ci[2]),
            lin.tolist(),
            _rot_xz(10.0 + 7 * n),
            _translate(co[0] + shift[0], co[1] + shift[1], co[2] + shift[2]),
        )
        cls = ("uint8", "int16", "uint16", "single", "double", "single")[n]
        for m in ("linear", "nearest", "cubic"):
            cases.append(
                _case(
                    f"generic{n}_{cls}_{m}",
                    [13, 15, 11],
                    cls,
                    T,
                    [12, 16, 13],
                    m,
                    kind="random",
                    seed=7,
                )
            )
    return cases


def _edges():
    cases = []
    size = [8, 9, 7]
    # translations at the edges of the input (in-bounds rules of warp3d and interp3)
    for n, (ax, d) in enumerate(
        [
            (0, 2.0**-40),
            (0, -(2.0**-40)),
            (1, 2.0**-52),
            (1, -(2.0**-52)),
            (2, 0.5),
            (2, -0.5),
            (0, 8.0),
            (1, -8.0),
            (2, 6.0),
            (0, 1 - 2.0**-53),
        ]
    ):
        t = [0.0, 0.0, 0.0]
        t[ax] = d
        T = _translate(*t)
        for m in ("linear", "nearest", "cubic"):
            cases.append(
                _case(
                    f"edge_translate{n}_single_{m}",
                    size,
                    "single",
                    T,
                    size,
                    m,
                    kind="random",
                    seed=8,
                )
            )
    # identity: the input itself
    for cls in ("uint16", "single"):
        for m in ("linear", "nearest", "cubic"):
            cases.append(
                _case(
                    f"identity_{cls}_{m}",
                    size,
                    cls,
                    _translate(0, 0, 0),
                    size,
                    m,
                    kind="random",
                    seed=9,
                )
            )
    # fill values cast like the input (uint16: rounding and saturation)
    T = _mul(_shear(0.7, -0.7), _rot_xz(12.0))
    for n, (cls, fill) in enumerate(
        [
            ("uint16", 70000.6),
            ("uint16", 2.5),
            ("uint16", -3.0),
            ("single", 7.25),
            ("single", 0.1),
            ("double", -1.5),
        ]
    ):
        for m in ("linear", "cubic"):
            cases.append(
                _case(
                    f"fill{n}_{cls}_{m}", size, cls, T, [8, 14, 9], m, fill, kind="random", seed=10
                )
            )
    # method names as imwarp accepts them
    for n, m in enumerate(("bilinear", "bicubic", "Nearest", "LINEAR", "cub")):
        cases.append(
            _case(f"method_name{n}", size, "single", T, [8, 14, 9], m, kind="random", seed=11)
        )
    # the smallest volumes each method takes
    for n, sz in enumerate(([3, 3, 3], [2, 2, 2], [3, 5, 2], [1, 6, 4], [6, 1, 4])):
        for m in ("linear", "nearest", "cubic"):
            if m == "cubic" and min(sz) < 3:
                continue
            cases.append(
                _case(
                    f"small{n}_single_{m}",
                    sz,
                    "single",
                    _translate(
                        0.3 if sz[1] > 1 else 0.0,
                        -0.2 if sz[0] > 1 else 0.0,
                        0.1 if sz[2] > 1 else 0.0,
                    ),
                    sz,
                    m,
                    kind="random",
                    seed=12,
                )
            )
    # NaN in the input
    cases += [
        Case(
            f"nan_single_{m}",
            code=MATLAB_CODE,
            args=[In("v"), np.array(T), [8.0, 14.0, 9.0], m, 0.0],
            python=py_imwarp,
            inputs={"v": Value(_nan_volume)},
        )
        for m in ("linear", "nearest", "cubic")
    ]
    return cases


def _nan_volume():
    v = np.asfortranarray((CounterRNG(13, "warp/nan").normal((8, 9, 7)) * 100.0).astype(np.float32))
    v[3, 4, 2] = np.float32(np.nan)
    v[0, 0, 0] = -np.float32(np.nan)
    return v


def _errors():
    good = _translate(0.5, 0.5, 0.5)
    cases = []
    for cid, size, T, m, ident in [
        ("error_2d_input", [6, 7, 1], good, "linear", "images:imwarp:tformDoesNotAgreeWithSizeOfA"),
        (
            "error_last_column",
            [6, 7, 5],
            [[1.0, 0, 0, 0.5], [0, 1.0, 0, 0], [0, 0, 1.0, 0], [0, 0, 0, 1.0]],
            "linear",
            "images:geotrans:invalidAffineMatrix",
        ),
        (
            "error_singular",
            [6, 7, 5],
            [[1.0, 0, 0, 0], [2.0, 0, 0, 0], [0, 0, 1.0, 0], [0, 0, 0, 1.0]],
            "linear",
            "images:geotrans:singularTransformationMatrix",
        ),
        (
            "error_not_finite",
            [6, 7, 5],
            [[1.0, 0, 0, 0], [0, 1.0, 0, 0], [0, 0, 1.0, 0], [math.inf, 0, 0, 1.0]],
            "linear",
            True,
        ),
        (
            "error_cubic_one_sample",
            [6, 1, 5],
            good,
            "cubic",
            "MATLAB:griddedInterpolant:DegenerateGridErrId",
        ),
        ("error_method", [6, 7, 5], good, "bi", True),
    ]:
        cases.append(
            _case(
                cid,
                size,
                "single",
                T,
                size,
                m,
                kind="random",
                seed=14,
                expect_error=ident,
                py_error=(warp.ImwarpError,),
            )
        )
    return cases


def _halves():
    """values on exact halves: a half-voxel shift (linear: (a + b) / 2; cubic: the weights
    -1/16, 9/16, 9/16, -1/16), negative ones too (int16): MATLAB's two roundings"""
    cases = []
    for ax in range(3):
        t = [0.0, 0.0, 0.0]
        t[ax] = 0.5
        for cls in ("int16", "uint16"):
            for m in ("linear", "cubic"):
                cases.append(
                    _case(
                        f"halves{ax}_{cls}_{m}",
                        [9, 10, 8],
                        cls,
                        _translate(*t),
                        [9, 10, 8],
                        m,
                        kind="random",
                        seed=15,
                    )
                )
    return cases


CASES = _families() + _edges() + _errors() + _halves()


def _divide_cases():
    """MATLAB's a / T (the solve inside imwarp's composite matrix) for every transform above with
    its output view's shift, and matrices whose pivots tie"""
    seen = set()
    out = []
    extra = [
        _rot_xz(45.0),  # |cos| == |sin|, no translation: tied pivots
        _mul(_rot_xz(45.0), _translate(3.0, 1.0, -2.0)),
        _mul(_scale(2.0, 1.0, 2.0), _rot_xz(-45.0)),
        _rot_xz(90.0),
        _mul(_shear(-1.0, 1.0), _rot_xz(45.0)),
    ]
    items = [(c.args[1], c.args[2]) for c in CASES if not c.expect_error] + [
        (np.array(T), [9.0, 9.0, 9.0]) for T in extra
    ]
    for T, view in items:
        T = np.asarray(T, dtype=np.float64)
        a = np.eye(4)
        if len(view) == 6:
            a[3, :3] = [view[1] - 1, view[0] - 1, view[2] - 1]
        key = T.tobytes() + a.tobytes()
        if key in seen:
            continue
        seen.add(key)
        out.append(
            Case(
                f"divide{len(out)}",
                code="out = {in{1} / in{2}};",
                args=[a, T],
                python=lambda a, T: np.array(warp._mrdivide(a.tolist(), T.tolist())),
            )
        )
    return out


DIVIDE_CASES = _divide_cases()


@oracle.parametrize(CASES + DIVIDE_CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


def test_oracle_inputs_are_deterministic():
    for c in CASES[::17]:
        for x in c.inputs.values():
            a, b = x.make(), x.make()
            assert _same_bits(a, b)


# ------------------------------------------------------------------------------------------
# a NumPy model of MATLAB's arithmetic per output voxel (MATLAB's form inside), for tests
# without MATLAB
# ------------------------------------------------------------------------------------------


def _model(v_yxz, t, out_yxz, method, fill=0.0):
    V = np.asarray(v_yxz, dtype=np.float64)
    ny, nx, nz = out_yxz
    i = np.arange(1, ny + 1, dtype=np.float64)[:, None, None]
    j = np.arange(1, nx + 1, dtype=np.float64)[None, :, None]
    k = np.arange(1, nz + 1, dtype=np.float64)[None, None, :]
    shape = (ny, nx, nz)
    out = np.full(shape, fill)
    if method == "cubic":
        p = [
            np.broadcast_to(((t[0, a] * j + t[1, a] * i) + t[2, a] * k) + t[3, a], shape)
            for a in range(3)
        ]
        G = np.transpose(V, (1, 0, 2))  # (x, y, z)
        E = G
        for ax in range(3):
            E = np.moveaxis(E, ax, 0)
            E = np.concatenate(
                [(3 * (E[0] - E[1]) + E[2])[None], E, (3 * (E[-1] - E[-2]) + E[-3])[None]]
            )
            E = np.moveaxis(E, 0, ax)
        inb = np.ones(shape, bool)
        for a in range(3):
            inb &= (p[a] >= 1) & (p[a] <= G.shape[a])
        ms, ws = [], []
        for a in range(3):
            q = p[a][inb]
            m = np.minimum(np.floor(q), G.shape[a] - 1)
            tt = q - m
            sq = tt * tt
            ws.append(
                (
                    ((-0.5 * tt + 1) * tt - 0.5) * tt,
                    (1.5 * tt - 2.5) * sq + 1,
                    ((-1.5 * tt + 2) * tt + 0.5) * tt,
                    (0.5 * tt - 0.5) * sq,
                )
            )
            ms.append(m.astype(np.int64))

        def s4(w, e):
            return ((w[0] * e[0] + w[1] * e[1]) + w[2] * e[2]) + w[3] * e[3]

        lz = [
            s4(
                ws[1],
                [
                    s4(ws[0], [E[ms[0] - 1 + a, ms[1] - 1 + b, ms[2] - 1 + c] for a in range(4)])
                    for b in range(4)
                ],
            )
            for c in range(4)
        ]
        out[inb] = s4(ws[2], lz)
        return out
    q = [
        np.broadcast_to(i * t[1, a] + (j * t[0, a] + (k * t[2, a] + (t[3, a] - 1))), shape)
        for a in range(3)
    ]
    qs, qf, qt = q
    n1, n2, n3 = V.shape
    inb = (qf >= 0) & (qf <= n1 - 1) & (qs >= 0) & (qs <= n2 - 1) & (qt >= 0) & (qt <= n3 - 1)
    qf, qs, qt = qf[inb], qs[inb], qt[inb]
    if method == "nearest":
        out[inb] = V[
            np.floor(qf + 0.5).astype(int),
            np.floor(qs + 0.5).astype(int),
            np.floor(qt + 0.5).astype(int),
        ]
        return out
    f0, s0, t0 = np.floor(qf), np.floor(qs), np.floor(qt)
    fd, sd, td = qf - f0, qs - s0, qt - t0
    f0, s0, t0 = f0.astype(int), s0.astype(int), t0.astype(int)
    f1, s1, t1 = np.ceil(qf).astype(int), np.ceil(qs).astype(int), np.ceil(qt).astype(int)
    c00 = V[f0, s0, t0] * (1 - sd) + V[f0, s1, t0] * sd
    c10 = V[f1, s0, t0] * (1 - sd) + V[f1, s1, t0] * sd
    c01 = V[f0, s0, t1] * (1 - sd) + V[f0, s1, t1] * sd
    c11 = V[f1, s0, t1] * (1 - sd) + V[f1, s1, t1] * sd
    c0 = c00 * (1 - fd) + c10 * fd
    c1 = c01 * (1 - fd) + c11 * fd
    out[inb] = c0 * (1 - td) + c1 * td
    return out


def _model_class(r, cls, method):
    if method == "nearest":  # the sample itself
        return r.astype(cls)
    if np.dtype(cls).kind in "ui":
        info = np.iinfo(cls)
        if method == "cubic":  # MATLAB's cast
            r = np.where(np.isnan(r), 0, np.where(r < 0, np.ceil(r - 0.5), np.floor(r + 0.5)))
        else:  # warp3d: the truncation of v -+ 0.5
            r = np.trunc(np.where(r < 0, r - 0.5, r + 0.5))
        return np.clip(r, info.min, info.max).astype(cls)
    return r.astype(cls)


def _random_case(seed, cls):
    rng = np.random.default_rng(seed)
    size = tuple(int(s) for s in rng.integers(4, 14, 3))  # (y, x, z)
    T = np.eye(4)
    T[:3, :3] += rng.normal(size=(3, 3)) * 0.3
    T[3, :3] = rng.normal(size=3) * 2
    if np.dtype(cls).kind == "f":
        v = (rng.normal(size=size) * 100).astype(cls)
    else:
        info = np.iinfo(cls)
        v = rng.integers(info.min, info.max, size=size, endpoint=True).astype(cls)
    out = tuple(int(s) for s in rng.integers(4, 14, 3))
    return np.asfortranarray(v), T, out


@pytest.mark.parametrize("cls", [np.uint8, np.int16, np.uint16, np.float32, np.float64])
@pytest.mark.parametrize("method", ["linear", "nearest", "cubic"])
def test_kernel_matches_the_model(cls, method):
    for seed in range(12):
        v, T, out = _random_case(seed, cls)
        t = warp.composite_matrix(T, output_size=arrays.size_from_matlab(out))
        want = _model_class(_model(v, t, out, method), cls, method)
        got = warp.imwarp(_from_matlab(v), T, method, output_size=arrays.size_from_matlab(out))
        assert got.dtype == np.dtype(cls) and got.shape == (out[2], out[0], out[1])
        assert _same_bits(_to_matlab(got), want), (seed, method)


@pytest.mark.parametrize("method", ["linear", "nearest", "cubic"])
def test_results_do_not_depend_on_threads(method):
    v, T, out = _random_case(99, np.float32)
    vol = _from_matlab(v)
    r = [
        warp.imwarp(vol, T, method, output_size=arrays.size_from_matlab(out), num_threads=n)
        for n in (1, 2, 3)
    ]
    assert all(_same_bits(r[0], x) for x in r[1:])


def test_other_layouts_are_copied():
    v, T, out = _random_case(5, np.uint16)
    vol = _from_matlab(v)
    a = warp.imwarp(vol, T, output_size=arrays.size_from_matlab(out))
    b = warp.imwarp(np.asfortranarray(vol), T, output_size=arrays.size_from_matlab(out))
    assert np.array_equal(a, b) and b.flags.c_contiguous


def test_a_box_is_its_part_of_the_full_view():
    # a transform with exact arithmetic (dyadic): the box's voxels equal the full view's
    rng = np.random.default_rng(3)
    vol = (rng.random((9, 11, 13)) * 1000).astype(np.float32)
    T = _mul(_shear(0.25, -0.5), _scale(0.5, 2.0, 1.0), _translate(1.5, -2.25, 0.75))
    full = warp.imwarp(vol, T, output_size=(10, 12, 14))
    for box in ((1, 2, 3, 9, 10, 11), (0, 0, 0, 10, 12, 14), (3, 1, 6, 4, 7, 8)):
        part = warp.imwarp(vol, T, output_box=box)
        assert np.array_equal(part, full[box[0] : box[3], box[1] : box[4], box[2] : box[5]]), box


def test_axes_of_the_transform_are_matlabs():
    # a shift by +1 along MATLAB's x moves the content along the last NumPy axis
    vol = np.zeros((4, 5, 6), np.float32)
    vol[1, 2, 3] = 1
    out = warp.imwarp(vol, _translate(1, 0, 0), output_size=(4, 5, 6))
    assert out[1, 2, 4] == 1 and out.sum() == 1
    out = warp.imwarp(vol, _translate(0, 1, 0), output_size=(4, 5, 6))
    assert out[1, 3, 3] == 1 and out.sum() == 1
    out = warp.imwarp(vol, _translate(0, 0, 1), output_size=(4, 5, 6))
    assert out[2, 2, 3] == 1 and out.sum() == 1


def test_errors():
    vol = np.zeros((5, 6, 7), np.float32)
    T = _translate(0, 0, 0)
    with pytest.raises(warp.ImwarpError) as e:
        warp.imwarp(np.zeros((1, 6, 7), np.float32), T, output_size=(1, 6, 7))
    assert e.value.identifier == "images:imwarp:tformDoesNotAgreeWithSizeOfA"
    with pytest.raises(warp.ImwarpError):
        warp.imwarp(np.zeros((6, 7), np.float32), T, output_size=(1, 6, 7))
    with pytest.raises(TypeError):
        warp.imwarp(vol.astype(np.int32), T, output_size=(5, 6, 7))
    with pytest.raises(warp.ImwarpError) as e:
        warp.imwarp(vol, np.eye(3), output_size=(5, 6, 7))
    assert e.value.identifier == "MATLAB:affine3d:set:T:incorrectSize"
    with pytest.raises(ValueError):
        warp.imwarp(vol, T)
    with pytest.raises(ValueError):
        warp.imwarp(vol, T, output_size=(5, 6, 7), output_box=(0, 0, 0, 5, 6, 7))
    with pytest.raises(ValueError):
        warp.imwarp(vol, T, output_box=(0, 0, 3, 5, 6, 3))
    with pytest.raises(NotImplementedError):
        warp.imwarp(np.zeros((2, 6, 7), np.float32), T, "cubic", output_size=(2, 6, 7))
    with pytest.raises(warp.ImwarpError):
        warp.imwarp(vol, T, "bi", output_size=(5, 6, 7))
    with pytest.raises(warp.ImwarpError):
        warp.imwarp(vol, T, "spline", output_size=(5, 6, 7))
    assert warp.imwarp(vol, T, "bilinear", output_size=(5, 6, 7)).shape == (5, 6, 7)
