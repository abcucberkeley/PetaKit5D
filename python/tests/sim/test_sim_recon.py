"""SIM reconstruction (petakit5d.sim, task C06) against MATLAB R2024b's ``third_parties/sim3D``:
sim_PSFtoOTF_gen, simRecon, simReconFrame, deskewPhasesFrame, separatePhases and
combinePhases, on synthetic 5-phase lattice data with a known pattern (``sim.synthetic``: beads
under a lattice of period 1.2021 um along y, a bead image as the PSF), through the MATLAB-name
layer; the 2020 lattice SIM data of the demo data (``sim3D_demo_data``: 5 and 3 phases) with
the settings verified against the original 2020 workflow (C06's report); and GPU cases (MATLAB's
gpuArray path against CuPy's).

The dataset wrappers fail in today's MATLAB (``cpuOnlyNodes``, which the job framework no longer
takes): reference-only cases record that error; the port runs them, tested here in local mode.

Tolerances: docs/parity/sim.md (proposed rows R1 to R6).
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from petakit5d import matlab
from petakit5d.core.arrays import to_matlab
from petakit5d.sim import synthetic
from petakit5d.sim.recon import SimError
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, Data, In, InPath, Tree, Value, mcell
from petakit5d.testing.compare import tol, tree

# tolerances (docs/parity/sim.md, proposed)
R4_OTF = tol(
    reason="proposed (C06, R4): the OTF (single phase ramps, the fit's sums)",
    rel_max=1e-5,
    rel_l2=1e-5,
)
R5_RECON = tol(
    reason="proposed (C06, R5): the reconstruction on the CPU", rel_max=2e-4, rel_l2=1e-4
)
R6_GPU = tol(
    reason="proposed (C06, R6): the reconstruction on a GPU (CuPy against gpuArray)",
    rel_max=1e-4,
    rel_l2=5e-5,
)
ERRORS = (SimError, ValueError, IndexError, OSError)

SHAPE = (10, 36, 32)  # (z, y, x) of one phase
PSF_SHAPE = (11, 32, 30)
BEADS = [
    (2.5, 8.2, 7.1, 900.0), (4.0, 17.3, 21.6, 1400.0), (6.2, 25.4, 9.8, 1100.0),
    (7.7, 12.6, 26.3, 800.0), (3.3, 29.1, 17.2, 1250.0), (5.1, 5.9, 15.4, 700.0),
    (8.4, 21.7, 4.6, 1000.0), (1.6, 14.4, 12.9, 1300.0),
]  # fmt: skip


def data_m(dtype=np.uint16, shape=SHAPE, nphases=5, **kw) -> np.ndarray:
    """the raw phases in MATLAB's form (y, x, pages)"""
    return to_matlab(synthetic.phase_stack(shape, BEADS, nphases=nphases, dtype=dtype, **kw))


def psf_m(nphases=5, **kw) -> np.ndarray:
    return to_matlab(synthetic.phase_psf(PSF_SHAPE, nphases=nphases, **kw))


def write_tif(path, matlab_form):
    from petakit5d.core.arrays import from_matlab
    from petakit5d.testing.synth import write_tiff

    write_tiff(path, from_matlab(matlab_form))


def frame_tree(*, data=None, psf=True, existing=False):
    def make(folder):
        folder = Path(folder)
        write_tif(folder / "data.tif", data_m() if data is None else data())
        if psf:
            write_tif(folder / "psf.tif", psf_m())
        if existing:
            (folder / "sim_recon").mkdir()
            write_tif(
                folder / "sim_recon" / "data_recon.tif", np.ones((4, 4, 2), np.float32, order="F")
            )

    return make


# ---------------------------------------------------------------------------------------------
# sim_PSFtoOTF_gen and simRecon
# ---------------------------------------------------------------------------------------------


def py_otf(psf, *opts):
    return matlab.sim_PSFtoOTF_gen(psf, "useGPU", False, *opts)


_RECON = (
    "otf = sim_PSFtoOTF_gen(in{1}, 'useGPU', false, 'gpuPrecision', 'single', in{3}{:}); "
    "out = {simRecon(in{2}, otf, 'useGPU', false, in{3}{:}, in{4:end})};"
)


def py_recon(psf, data, common, *opts):
    otf = matlab.sim_PSFtoOTF_gen(psf, "useGPU", False, "gpuPrecision", "single", *common)
    return matlab.simRecon(data, otf, "useGPU", False, *common, *opts)


def _recon_cases() -> list[Case]:
    psf = {"psf": Value(psf_m)}
    cs = [
        Case("otf.single", matlab="sim_PSFtoOTF_gen",
             args=[In("psf"), "useGPU", False, "gpuPrecision", "single"],
             python=lambda psf, *a: py_otf(psf, *a[2:]), inputs=psf, compare=R4_OTF),
        Case("otf.double", matlab="sim_PSFtoOTF_gen", args=[In("psf"), "useGPU", False],
             python=lambda psf, *a: py_otf(psf), inputs=psf, compare=R4_OTF),
        Case("otf.save_fails", matlab="sim_PSFtoOTF_gen",
             args=[In("psf"), "useGPU", False, "saveOTF", True],
             python=lambda psf, *a: py_otf(psf, "saveOTF", True), inputs=psf,
             expect_error=True, py_error=ERRORS),
    ]  # fmt: skip
    three = ["nphases", 3, "norders", 3, "lattice_period", 0.52932, "phase_step", 0.170]
    for cid, common, opts, make in [
        ("recon.defaults", [], [], lambda: data_m()),
        ("recon.frame_defaults", [], ["perdecomp", False, "edgeTaper", False], lambda: data_m()),
        ("recon.wver0", [], ["wVer", 0], lambda: data_m()),
        ("recon.double", [], ["gpuPrecision", "double"], lambda: data_m()),
        ("recon.no_apodize", [], ["apodize", False], lambda: data_m()),
        ("recon.single_input", [], [], lambda: data_m(np.float32)),
        ("recon.3beam", [], ["islattice", False], lambda: data_m()),
        # 20 planes: the zero order's mask (the third order's, MATLAB's quirk) is axially
        # smaller than the first orders' on this grid
        ("recon.3beam.20planes", [], ["islattice", False], lambda: data_m(shape=(20, 36, 32))),
        ("recon.3phase", three, ["NA_ext", 0.3],
         lambda: data_m(nphases=3, lattice_period=0.52932, phase_step=0.170, order_amps=(0.45,),
                        order_sigmas=(0.9,))),
    ]:  # fmt: skip
        psf_make = (
            (
                lambda: psf_m(
                    3,
                    lattice_period=0.52932,
                    phase_step=0.170,
                    order_amps=(0.45,),
                    order_sigmas=(0.9,),
                )
            )
            if common
            else psf_m
        )
        cs.append(
            Case(
                cid,
                code=_RECON,
                args=[In("psf"), In("data"), mcell(*common), *opts],
                python=py_recon,
                py_args=[In("psf"), In("data"), list(common), *opts],
                inputs={"psf": Value(psf_make), "data": Value(make)},
                compare=R5_RECON,
            )
        )
    for cid, make, opts in [
        ("recon.nz_not_whole", lambda: data_m()[:, :, :-1], []),
        ("recon.two_orders", lambda: data_m(), ["norders", 2]),
        ("recon.wver2_fails", lambda: data_m(), ["wVer", 2, "perdecomp", False]),
    ]:
        cs.append(
            Case(
                cid,
                code=_RECON,
                args=[In("psf"), In("data"), mcell(), *opts],
                python=py_recon,
                py_args=[In("psf"), In("data"), [], *opts],
                inputs={"psf": Value(psf_m), "data": Value(make)},
                expect_error=True,
                py_error=ERRORS,
            )
        )
    return cs


# ---------------------------------------------------------------------------------------------
# simReconFrame, deskewPhasesFrame, separatePhases, combinePhases
# ---------------------------------------------------------------------------------------------

_FRAME = "simReconFrame([in{1}, in{2}], [in{1}, 'psf.tif'], 'useGPU', false, in{3:end}); out = {};"
_RECON_TREE = tree(per_file={"*.tif": R5_RECON})


def py_frame(d, name, *opts):
    matlab.simReconFrame(d + name, d + "psf.tif", "useGPU", False, *opts)
    return []


def py_phases(d):
    matlab.separatePhases(d + "data.tif")
    matlab.combinePhases(d + "data.tif")
    return []


def _frame_cases() -> list[Case]:
    cs = []
    for cid, opts, make in [
        ("frame.psf", [], frame_tree()),
        ("frame.16bit", ["Save16bit", True], frame_tree()),
        ("frame.test_call", ["perdecomp", True, "edgeTaper", True, "DS", True], frame_tree()),
        ("frame.erode_after", ["EdgeErosion", 2], frame_tree()),
        ("frame.erode_before_mask", ["EdgeErosion", 2, "ErodeBefore", True, "ErodeAfter", False,
                                     "SaveMaskfile", True], frame_tree()),
        ("frame.chunks", ["ChunkSize", [20, 18, 6], "Overlap", 4], frame_tree()),
        # odd chunks in y (15 voxels): every chunk widened, the first at its end, the others at
        # their start
        ("frame.odd_chunks", ["ChunkSize", [8, 11, 6], "Overlap", 4], frame_tree()),
        ("frame.flip", ["flipZstack", True], frame_tree()),
        ("frame.blank_chunks", ["occThresh", 1.5], frame_tree()),
        ("frame.existing", [], frame_tree(existing=True)),
    ]:  # fmt: skip
        cs.append(
            Case(
                cid,
                code=_FRAME,
                args=[InPath("d"), "data.tif", *opts],
                nargout=0,
                python=py_frame,
                inputs={"d": Tree(make, "d")},
                tree=_RECON_TREE,
            )
        )
    cs.append(
        Case(
            "frame.missing",
            code=_FRAME,
            args=[InPath("d"), "nothing.tif"],
            nargout=0,
            python=py_frame,
            inputs={"d": Tree(frame_tree(), "d")},
            tree=_RECON_TREE,
        )
    )
    cs.append(
        Case(
            "frame.no_otf",
            code="simReconFrame([in{1}, 'data.tif'], [], 'useGPU', false); out = {};",
            args=[InPath("d")],
            nargout=0,
            python=lambda d: matlab.simReconFrame(d + "data.tif", None, "useGPU", False),
            inputs={"d": Tree(frame_tree(), "d")},
            expect_error=True,
            py_error=ERRORS,
        )
    )
    for cid, opts in [
        ("deskew.ds", ["Rotate", False]),
        ("deskew.ds.16bit", ["Rotate", False, "Save16bit", True]),
        ("deskew.ds.forward", ["Reverse", False]),
        ("deskew.dsr", ["Rotate", True]),
    ]:
        cs.append(
            Case(
                cid,
                code="deskewPhasesFrame([in{1}, 'data.tif'], 0.108, 0.3, 'SkewAngle', 32.45, "
                "in{2:end}); out = {};",
                args=[InPath("d"), *opts],
                nargout=0,
                python=lambda d, *o: matlab.deskewPhasesFrame(
                    d + "data.tif", 0.108, 0.3, "SkewAngle", 32.45, *o
                ),
                inputs={"d": Tree(frame_tree(psf=False), "d")},
            )
        )
    cs.append(
        Case(
            "phases.separate_combine",
            code="separatePhases([in{1}, 'data.tif']); combinePhases([in{1}, 'data.tif']); "
            "out = {};",
            args=[InPath("d")],
            nargout=0,
            python=py_phases,
            inputs={"d": Tree(frame_tree(psf=False), "d")},
        )
    )
    return cs


# ---------------------------------------------------------------------------------------------
# the dataset wrappers in MATLAB: their error (reference only)
# ---------------------------------------------------------------------------------------------


def _wrapper_cases() -> list[Case]:
    return [
        Case(
            "wrapper.recon_fails_in_matlab",
            code="simRecon_data_wrapper(in{1}, {[in{1}, 'psf.tif']}, 'ChannelPatterns', {'data'}, "
            "'parseCluster', false, 'useGPU', false); out = {};",
            args=[InPath("d")],
            nargout=0,
            python=None,
            inputs={"d": Tree(frame_tree(), "d")},
            expect_error="MATLAB:InputParser:UnmatchedParameter",
        ),
        Case(
            "wrapper.deskew_fails_in_matlab",
            code="deskewPhases_data_wrapper(in{1}, 0.108, 0.3, 'ChannelPatterns', {'data'}, "
            "'parseCluster', false); out = {};",
            args=[InPath("d")],
            nargout=0,
            python=None,
            inputs={"d": Tree(frame_tree(psf=False), "d")},
            expect_error="MATLAB:InputParser:UnmatchedParameter",
        ),
        Case(
            "wrapper.rotate_fails_in_matlab",
            code="Rotate3D_data_wrapper(in{1}, 0.108, 0.3, 'ChannelPatterns', {'data'}, "
            "'parseCluster', false); out = {};",
            args=[InPath("d")],
            nargout=0,
            python=None,
            inputs={"d": Tree(frame_tree(psf=False), "d")},
            expect_error="MATLAB:InputParser:UnmatchedParameter",
        ),
    ]


# ---------------------------------------------------------------------------------------------
# the 2020 lattice SIM data (demo data), with the settings verified against the 2020 workflow
# ---------------------------------------------------------------------------------------------

_DEMO = "$PETAKIT5D_DEMO_DATA/sim3D_demo_data/"
_MARKS = ("slow",)


def py_real_frame(data, otf, *opts):
    matlab.simReconFrame(data, otf, *opts)
    return []


def py_real_frame_mat(data, otf_mat, *opts):
    from petakit5d.io.matfile import read_mat

    o = read_mat(otf_mat)
    otf = next(v for k, v in o.items() if not k.startswith("__"))
    matlab.simReconFrame(data, np.asarray(otf), *opts)
    return []


def _realdata_cases() -> list[Case]:
    d5 = {"data": Data(_DEMO + "Cell_5phase_small.tif", "Cell_5phase_small.tif"),
          "psf": Data(_DEMO + "PSF_5phase.tif", "PSF_5phase.tif")}  # fmt: skip
    d5m = {"data": Data(_DEMO + "Cell_5phase_small.tif", "Cell_5phase_small.tif"),
           "otf": Data(_DEMO + "PSF_5phas_OTF_normalized.mat", "otf.mat")}  # fmt: skip
    d3 = {"data": Data(_DEMO + "Cell_3phase.tif", "Cell_3phase.tif"),
          "psf": Data(_DEMO + "PSF_3phase.tif", "PSF_3phase.tif")}  # fmt: skip
    three = ["nphases", 3, "norders", 3, "lattice_period", 0.52932, "phase_step", 0.170,
             "pxl_dim_data", [0.11, 0.11, 0.2], "pxl_dim_PSF", [0.11, 0.11, 0.1],
             "NA_ext", 0.3]  # fmt: skip
    load = "o = load(in{2}); f = fieldnames(o); otf = o.(f{1}); "
    return [
        # the live script's values (sim3D's defaults), the OTF from the PSF
        Case("realdata.frame5", code="simReconFrame(in{1}, in{2}, in{3:end}); out = {};",
             args=[InPath("data"), InPath("psf"), "useGPU", False, "DS", True], nargout=0,
             python=py_real_frame, inputs=d5, tree=_RECON_TREE, marks=_MARKS),
        # testSimRecon.m's call, with the 2020 OTF
        Case("realdata.frame5_test_call",
             code=load + "simReconFrame(in{1}, otf, in{3:end}); out = {};",
             args=[InPath("data"), InPath("otf"), "useGPU", False, "perdecomp", True,
                   "edgeTaper", True, "DS", True], nargout=0, python=py_real_frame_mat,
             inputs=d5m, tree=_RECON_TREE, marks=_MARKS),
        # the 3-phase set with its live script's values
        Case("realdata.frame3", code="simReconFrame(in{1}, in{2}, in{3:end}); out = {};",
             args=[InPath("data"), InPath("psf"), "useGPU", False, "DS", True, "perdecomp", True,
                   *three],
             nargout=0, python=py_real_frame, inputs=d3, tree=_RECON_TREE, marks=_MARKS),
    ]  # fmt: skip


# ---------------------------------------------------------------------------------------------
# GPU: MATLAB's gpuArray path against CuPy's
# ---------------------------------------------------------------------------------------------


def py_gpu_otf(psf):
    return matlab.sim_PSFtoOTF_gen(psf, "useGPU", True, "gpuPrecision", "single")


def py_gpu_recon(psf, data):
    otf = matlab.sim_PSFtoOTF_gen(psf, "useGPU", True, "gpuPrecision", "single")
    return matlab.simRecon(data, otf, "useGPU", True)


def py_gpu_frame(d, name, *opts):
    matlab.simReconFrame(d + name, d + "psf.tif", "useGPU", True, *opts)
    return []


def _gpu_cases() -> list[Case]:
    gtree = tree(per_file={"*.tif": R6_GPU})
    return [
        Case("gpu.otf",
             code="out = {sim_PSFtoOTF_gen(in{1}, 'useGPU', true, 'gpuPrecision', 'single')};",
             args=[In("psf")], python=py_gpu_otf, inputs={"psf": Value(psf_m)}, device="gpu",
             compare=R4_OTF),
        Case("gpu.recon",
             code="otf = sim_PSFtoOTF_gen(in{1}, 'useGPU', true, 'gpuPrecision', 'single'); "
             "out = {gather(simRecon(gpuArray(in{2}), otf, 'useGPU', true))};",
             args=[In("psf"), In("data")], python=py_gpu_recon,
             inputs={"psf": Value(psf_m), "data": Value(data_m)}, device="gpu", compare=R6_GPU),
        Case("gpu.frame",
             code="simReconFrame([in{1}, in{2}], [in{1}, 'psf.tif'], 'useGPU', true, in{3:end}); "
             "out = {};",
             args=[InPath("d"), "data.tif", "DS", True], nargout=0, python=py_gpu_frame,
             inputs={"d": Tree(frame_tree(), "d")}, device="gpu", tree=gtree),
        Case("gpu.realdata.frame5_test_call",
             code="o = load(in{2}); f = fieldnames(o); otf = o.(f{1}); "
             "simReconFrame(in{1}, otf, in{3:end}); out = {};",
             args=[InPath("data"), InPath("otf"), "useGPU", True, "perdecomp", True,
                   "edgeTaper", True, "DS", True], nargout=0, python=py_real_frame_mat,
             inputs={"data": Data(_DEMO + "Cell_5phase_small.tif", "Cell_5phase_small.tif"),
                     "otf": Data(_DEMO + "PSF_5phas_OTF_normalized.mat", "otf.mat")},
             device="gpu", tree=gtree, marks=_MARKS),
    ]  # fmt: skip


CASES = _recon_cases() + _frame_cases() + _wrapper_cases() + _realdata_cases() + _gpu_cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# the port's own tests: physics on the synthetic data, the wrappers in local mode
# ---------------------------------------------------------------------------------------------


def _otf():
    from petakit5d.sim.otf import sim_psf_to_otf

    return sim_psf_to_otf(synthetic.phase_psf(PSF_SHAPE), use_gpu=False, gpu_precision="single")


def test_reconstruction_doubles_the_lateral_size_and_extends_the_spectrum():
    from petakit5d.sim.recon import sim_recon_volume

    raw = synthetic.phase_stack(SHAPE, BEADS)
    r = sim_recon_volume(raw, _otf(), use_gpu=False, perdecomp=False, edge_taper=False)
    assert r.shape == (SHAPE[0], 2 * SHAPE[1], 2 * SHAPE[2]) and r.dtype == np.float32
    # the brightest bead (z 4, y 17.3, x 21.6) is the brightest spot of its plane, at twice its
    # lateral position
    iy, ix = np.unravel_index(int(np.argmax(r[4])), r[4].shape)
    assert abs(iy - 2 * 17.3) <= 3 and abs(ix - 2 * 21.6) <= 3


def test_shift_vectors_start_at_the_lattice_and_refine_a_little(monkeypatch):
    """the guessed shift vectors are the lattice's frequency (orders -2, -1 along y), the masked
    NCC refines them by a fraction of their length, and the positive orders mirror them"""
    from petakit5d.sim import recon as R

    seen = {}
    real = R._assemble

    def spy(ops, dk_sep, o_scaled, b_pls, p_vec, shape, dk_data, k):
        seen["p"] = p_vec.copy()
        seen["dk"] = np.asarray(dk_data, dtype=np.float64)
        return real(ops, dk_sep, o_scaled, b_pls, p_vec, shape, dk_data, k)

    monkeypatch.setattr(R, "_assemble", spy)
    raw = synthetic.phase_stack(SHAPE, BEADS)
    R.sim_recon_volume(raw, _otf(), use_gpu=False, perdecomp=False, edge_taper=False)
    p, dk = seen["p"][:, :, 0], seen["dk"]
    for row, order in ((0, -2), (1, -1)):
        guess = order / 1.2021 / dk[0]
        assert abs(p[row, 0] - guess) < 0.2 * abs(guess) and abs(p[row, 1]) < 1.5
    assert np.array_equal(p[2], np.zeros(3))
    assert np.array_equal(p[3:], -p[1::-1])


def test_sim_recon_file_writes_the_result_and_its_mip(tmp_path):
    d = str(tmp_path) + "/"
    frame_tree()(tmp_path)
    from petakit5d.sim.file import sim_recon_file

    vol = sim_recon_file(d + "data.tif", d + "psf.tif", use_gpu=False)
    from petakit5d.io import read_tiff

    out = read_tiff(d + "sim_recon/data_recon.tif")
    mip = read_tiff(d + "sim_recon/MIPs/data_MIP_z.tif")
    assert out.shape == (SHAPE[0], 2 * SHAPE[1], 2 * SHAPE[2]) and out.dtype == np.float32
    assert np.all(out >= 0)
    assert np.array_equal(mip, vol.max(axis=0))


def test_sim_recon_dataset_runs_its_tasks_locally(tmp_path):
    from petakit5d.io import read_tiff
    from petakit5d.sim.dataset import sim_recon_dataset

    frame_tree()(tmp_path)
    os.rename(tmp_path / "data.tif", tmp_path / "cell_CamA.tif")
    done = sim_recon_dataset(str(tmp_path), [str(tmp_path / "psf.tif")], channel_patterns=["CamA"],
                             parse_cluster=False, use_gpu=False)  # fmt: skip
    assert done == [True]
    assert read_tiff(str(tmp_path / "sim_recon" / "cell_CamA_recon.tif")).shape[1:] == (72, 64)


def test_deskew_phases_dataset_and_rotate_dataset_run_locally(tmp_path):
    from petakit5d.sim.dataset import deskew_phases_dataset

    frame_tree(psf=False)(tmp_path)
    os.rename(tmp_path / "data.tif", tmp_path / "cell_CamA.tif")
    assert deskew_phases_dataset(str(tmp_path), 0.108, 0.3, channel_patterns=["CamA"],
                                 parse_cluster=False) == [True]  # fmt: skip
    assert (tmp_path / "DS" / "cell_CamA.tif").exists()


def test_automatic_processing_deskews_then_reconstructs(tmp_path, monkeypatch):
    from petakit5d.sim import automatic

    monkeypatch.setattr(automatic, "_PAUSE", 0.05)
    frame_tree()(tmp_path)
    os.rename(tmp_path / "data.tif", tmp_path / "cell_CamA.tif")
    automatic.sim_recon_automatic_processing(
        str(tmp_path), psfs=[str(tmp_path / "psf.tif")], channel_patterns=["CamA"],
        parse_cluster=False, use_gpu=False, par_pool_size=2, xy_pixel_size=0.108, dz=0.3,
    )  # fmt: skip
    assert (tmp_path / "DS" / "cell_CamA.tif").exists()
    assert (tmp_path / "DS" / "sim_recon" / "cell_CamA_recon.tif").exists()


@pytest.mark.gpu
def test_gpu_path_runs_on_the_gpu_and_agrees_with_the_cpu(monkeypatch):
    from petakit5d.sim import recon as R

    raw = synthetic.phase_stack(SHAPE, BEADS)
    otf = _otf()
    paths = []
    real = R._Ops.__init__

    def spy(self, gpu, *a):
        paths.append(gpu)
        real(self, gpu, *a)

    monkeypatch.setattr(R._Ops, "__init__", spy)
    g = R.sim_recon_volume(raw, otf, use_gpu=True, perdecomp=True, edge_taper=True)
    c = R.sim_recon_volume(raw, otf, use_gpu=False, perdecomp=True, edge_taper=True)
    assert paths == [True, False]
    assert isinstance(g, np.ndarray) and g.dtype == np.float32 and g.shape == c.shape
    assert np.linalg.norm(g - c) / np.linalg.norm(c) < 1e-3
