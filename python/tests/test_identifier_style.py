"""The port's own error and warning identifiers have one style, ``petakit5d:<area>:<name>``
(decided by the user on 2026-11-25; docs/guide/index.md, "Errors and warnings").

``<area>`` is the lowercase name of the petakit5d package or module that raises it (``params``,
``detection``, ``bioformats``), ``<name>`` lowerCamelCase. MATLAB's identifiers stay as MATLAB
has them: MathWorks' (``MATLAB:...`` and the toolboxes': ``images:``, ``stats:``, ...) and those
that MATLAB code raises and the port copies for parity (MATLAB_COPIED, each with the file that
raises it: PetaKit5D's MATLAB code beside python/, or the oracle harness's).

Scanned:
- the package's Python code (src/petakit5d, tools, examples, benchmarks): every string whose
  whole text has an identifier's form (names joined by colons), f-strings with the module
  constants they format filled in (one whose first or last part is another formatted value is
  not an identifier: ``f"custom:{name}"``);
- every text file of python/ (code, tests, docs) for mentions in the port's form (``petakit5d:``
  or ``pk5d:``, in any case): each one an identifier the package gives, in the style, one of
  MATLAB_COPIED or a figure's group id.
"""

import ast
import re
from pathlib import Path

import pytest

PYROOT = Path(__file__).resolve().parents[1]
REPO = PYROOT.parent  # PetaKit5D's repository, whose MATLAB code is beside python/
PACKAGE = PYROOT / "src" / "petakit5d"
CODE = ("src", "tools", "examples", "benchmarks")
TEXT = {".py", ".md", ".m", ".json", ".txt", ".toml", ".yml", ".cfg", ".cmake", ".c", ".cpp",
        ".h", ".hpp", ".cu"}  # fmt: skip
SKIP_DIRS = {"__pycache__", "build", "dist", "node_modules"}

NAME = "[A-Za-z][A-Za-z0-9_]*"
FORM = re.compile(rf"{NAME}(?::{NAME})+")
STYLE = re.compile(r"petakit5d:([a-z][a-z0-9_]*):([a-z][A-Za-z0-9]*)")
PORT_FORM = re.compile(rf"(?i)(?<![\w.:/\\-])(?:petakit5d|pk5d)(?::{NAME})+(?![\w:])")
HOLE = "\0"  # a formatted value of an f-string that is not a module constant

# the first parts of MathWorks' identifiers in the port (MATLAB's and its toolboxes': Image
# Processing, Statistics, Optimization, Parallel Computing, Database)
MATHWORKS = frozenset({"MATLAB", "images", "stats", "optimlib", "parallel", "database"})

# identifiers of MATLAB code that the port copies, by the file that raises them (from the
# repository's root; python/ for the oracle harness's MATLAB code)
MATLAB_COPIED = {
    "microscopeDataProcessing/crop/mex/crop3d_mex.c": ("crop:dataTypeError", "crop:inputError"),
    "microscopeDataProcessing/crop/mex/crop4d_mex.cpp": ("tiff:dataTypeError",),
    "microscopeDataProcessing/crop/mex/indexing4d_mex.c": ("indexing:inputError",),
    "microscopeDataProcessing/decon/decon_gpu_num.m": ("decon:noUsableGPU",),
    "microscopeDataProcessing/decon/mex/decon_lucy_accel_mex.cpp": (
        "decon_lucy_accel_mex:inputError",
    ),
    "microscopeDataProcessing/decon/mex/decon_lucy_rfft_gpu_mex.cu": (
        "decon_lucy_rfft_gpu_mex:inputError",
        "decon_lucy_rfft_gpu_mex:outSizeIsSmallerThanPsfSize",
    ),
    "microscopeDataProcessing/decon/mex/decon_lucy_rfft_mex.cpp": (
        "decon_lucy_rfft_mex:inputError",
        "decon_lucy_rfft_mex:outSizeIsSmallerThanPsfSize",
    ),
    "microscopeDataProcessing/decon/mex/decon_mask_edge_erosion_gpu_mex.cu": (
        "decon_mask_edge_erosion_gpu_mex:inputError",
    ),
    "microscopeDataProcessing/decon/mex/decon_mask_edge_erosion_mex.cpp": (
        "decon_mask_edge_erosion_mex:inputError",
    ),
    "microscopeDataProcessing/decon/mex/omw_backprojector_mex.cpp": (
        "omw_backprojector_mex:invalidInput",
    ),
    "microscopeDataProcessing/deskew_rotate/mex/volume_deskew_rotate_warp_mex.cpp": (
        "volume_warp:inputError",
    ),
    "microscopeDataProcessing/flatfield_correction/mex/flat_field_correction_3d_mex.cpp": (
        "flat_field_correction_3d_mex:dataTypeError",
        "flat_field_correction_3d_mex:inputError",
    ),
    "microscopeDataProcessing/chromatic_shift_correction/"
    "XR_chromatic_shift_correction_data_wrapper.m": (
        "XR_chromatic_shift_correction_data_wrapper:noOffsets",
        "XR_chromatic_shift_correction_data_wrapper:offsetSize",
        "XR_chromatic_shift_correction_data_wrapper:padValue",
    ),
    "microscopeDataProcessing/FSC_analysis/XR_FSC.m": ("XR_FSC:noSignal",),
    "microscopeDataProcessing/FSC_analysis/XR_FSC_resolution.m": (
        "XR_FSC_resolution:oddSectorNumber",
    ),
    "microscopeDataProcessing/io/cpp-zarr/mexSrc/parallelreadzarrmex.cpp": ("zarr:inputError",),
    "microscopeDataProcessing/stitch/mex/any_4th_dim_mex.cpp": (
        "any_4th_dim_mex:dataRangeError",
        "any_4th_dim_mex:dataTypeError",
        "any_4th_dim_mex:inputError",
    ),
    "microscopeDataProcessing/stitch/mex/feather_blending_3d_mex.cpp": (
        "feather_blending:dataTypeError",
        "feather_blending:inputError",
    ),
    "microscopeDataProcessing/stitch/mex/feather_distance_map_resize_3d_mex.c": (
        "distmap_resize:inputError",
    ),
    "microscopeDataProcessing/stitch/mex/integral_image_3d_mex.cpp": (
        "integral_image_3d:dataTypeError",
        "integral_image_3d:inputError",
    ),
    "microscopeDataProcessing/stitch/mex/normxcorr3_window_mex.cpp": (
        "normxcorr3_window_mex:inputError",
    ),
    "microscopeDataProcessing/stitch/mex/tile_blending_3d_mex.cpp": (
        "tile_blending:dataTypeError",
        "tile_blending:inputError",
    ),
    "microscopeDataProcessing/stitch/stitch_assignment_lsqlin.m": (
        "stitch_assignment_lsqlin:rowError",
    ),
    "microscopeDataProcessing/stitch/stitch_generate_imagelist_from_sqlite.m": (
        "stitch_generate_imagelist_from_sqlite:missingMetadataHashInView",
    ),
    "microscopeDataProcessing/stitch/XR_stitching_frame_zarr_dev_v1.m": (
        "XR_stitching_frame_zarr_dev_v1:axisOrder",
        "XR_stitching_frame_zarr_dev_v1:stitchedSize",
        "XR_stitching_frame_zarr_dev_v1:tileSizes",
    ),
    "microscopeDataProcessing/tools/MIP/mex/max_pooling_3d_mex.cpp": (
        "mat:dataTypeError",
        "pooling:inputError",
    ),
    "microscopeDataProcessing/tools/MIP/mex/min_bbox_3d_mex.cpp": ("minbbox:inputError",),
    "microscopeDataProcessing/tools/MIP/mex/range_max_3d_mex.cpp": ("range_max_3d:inputError",),
    "microscopeDataProcessing/unmix/mex/unmix_channels_sum_mex.cpp": (
        "unmix_channels:dataTypeError",
        "unmix_channels:inputError",
    ),
    "pointDetection/mex_functions/replace_nan_with_value.c": (
        "data:TypeError",
        "data:inputError",
    ),
    "utils/generic_computing_framework/read_task_metrics.m": ("read_task_metrics:sacct",),
    "utils/generic_computing_framework/slurm_alternatives.m": (
        "slurm_alternatives:entry",
        "slurm_alternatives:key",
        "slurm_alternatives:testOnly",
    ),
    "utils/generic_computing_framework/task_metrics.m": ("task_metrics:failed",),
    "utils/generic_computing_framework/task_metrics_run.m": ("task_metrics_run:failed",),
    "utils/generic_computing_framework/task_runner_run.m": ("task_runner_run:single",),
    "utils/mex_fallback.m": ("PetaKit5D:mexFallback",),
    # the oracle harness's MATLAB code: its mark for warnings it could not parse
    "python/tests/oracle_matlab/+pkpy_oracle/parse_warnings.m": ("pkpy:oracle:unparsedWarning",),
}
COPIED = {i: src for src, ids in MATLAB_COPIED.items() for i in ids}

# the figures' axes group ids (matplotlib set_gid: tools.figures, tools.psf_analysis)
FIGURE_GIDS = frozenset(
    f"pk5d:{g}"
    for g in ("colorbar", "linecuts_rw", "linecuts_xz", "linecuts_yz", "mask", "mip", "otf_xy",
              "otf_xz", "otf_yz", "psf_xy", "psf_xy_gamma", "psf_xz", "psf_xz_gamma", "psf_yz",
              "psf_yz_gamma", "tiles", "title")
)  # fmt: skip

# strings of the package's code in an identifier's form that are not identifiers
NOT_IDENTIFIERS = {
    **dict.fromkeys(FIGURE_GIDS, "a figure's axes group id"),
    **dict.fromkeys(
        ("port:string", "port:vector", "port:evenly", "port:whole", "port:poswhole"),
        "the registry's code of a validator word of the port (params.spec)",
    ),
    "no:cacheprovider": "pytest's -p option",
    "r:gz": "a tarfile mode",
    "crop_deskew_rotate_demo_function:crop_deskew_rotate_demo_function": "a task's function",
}


def kind(ident: str, areas: set[str]) -> str | None:
    """what a string in an identifier's form is: 'port' (the style, with an area of the
    package), 'MathWorks', 'copied' (MATLAB_COPIED), 'not an identifier', or None"""
    m = STYLE.fullmatch(ident)
    if m:
        return "port" if m.group(1) in areas else None
    if ident.split(":", 1)[0] in MATHWORKS:
        return "MathWorks"
    if ident in COPIED:
        return "copied"
    if ident in NOT_IDENTIFIERS:
        return "not an identifier"
    return None


def package_areas(package: Path = PACKAGE) -> set[str]:
    """the names of the package's subpackages and modules"""
    names = {p.name for p in package.rglob("*") if p.is_dir() and (p / "__init__.py").is_file()}
    names |= {p.stem for p in package.rglob("*.py")}
    return {n for n in names if re.fullmatch("[a-z][a-z0-9_]*", n)}


def _constants(tree: ast.Module) -> dict[str, str]:
    """the module's string constants assigned once at its top level"""
    seen: dict[str, list] = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1:
            target, value = node.targets[0], node.value
        elif isinstance(node, ast.AnnAssign) and node.value is not None:
            target, value = node.target, node.value
        else:
            continue
        if isinstance(target, ast.Name):
            seen.setdefault(target.id, []).append(value)
    return {
        k: v[0].value
        for k, v in seen.items()
        if len(v) == 1 and isinstance(v[0], ast.Constant) and isinstance(v[0].value, str)
    }


def _fstring(node: ast.JoinedStr, consts: dict[str, str]) -> str:
    out = []
    for v in node.values:
        if isinstance(v, ast.Constant):
            out.append(v.value)
        elif (
            isinstance(v.value, ast.Name)
            and v.value.id in consts
            and v.conversion == -1
            and v.format_spec is None
        ):
            out.append(consts[v.value.id])
        else:
            out.append(HOLE)
    return "".join(out)


def code_identifiers(source: str) -> list[tuple[str, int]]:
    """the strings of a module's code in an identifier's form, with their lines (HOLE for a
    formatted value that is not a module constant; none with a HOLE in the first or last part,
    unless it starts in the port's form)"""
    tree = ast.parse(source)
    consts = _constants(tree)
    inner = {
        id(n)
        for f in ast.walk(tree)
        if isinstance(f, ast.JoinedStr)
        for n in ast.walk(f)
        if n is not f
    }
    out = []
    for node in ast.walk(tree):
        if id(node) in inner:
            continue
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            text = node.value
        elif isinstance(node, ast.JoinedStr):
            text = _fstring(node, consts)
        else:
            continue
        if not FORM.fullmatch(text.replace(HOLE, "x")):
            continue
        parts = text.split(":")
        if HOLE in parts[0] + parts[-1] and not PORT_FORM.match(text.replace(HOLE, "x")):
            continue
        out.append((node.lineno, node.col_offset, text))
    return [(t, n) for n, _, t in sorted(out)]


def mentions(text: str) -> list[tuple[str, int]]:
    """the mentions in the port's form in a text, with their lines"""
    return [(m.group(), text.count("\n", 0, m.start()) + 1) for m in PORT_FORM.finditer(text)]


def _files(root: Path, suffixes: set[str]):
    for p in sorted(root.rglob("*")):
        rel = p.relative_to(root).parts
        if (
            p.suffix in suffixes
            and p.is_file()
            and not any(
                d in SKIP_DIRS or d.startswith(".") or d.endswith(".egg-info") for d in rel[:-1]
            )
        ):
            yield p


def _rel(p: Path) -> str:
    return p.relative_to(PYROOT).as_posix()


@pytest.fixture(scope="module")
def scan():
    """the package's code strings (file, line, text) and the mentions in python/'s text files"""
    code, files = [], []
    for d in CODE:
        for p in _files(PYROOT / d, {".py"}):
            files.append(_rel(p))
            code += [(_rel(p), n, t) for t, n in code_identifiers(p.read_text(encoding="utf-8"))]
    found = []
    this = Path(__file__).resolve()
    for p in _files(PYROOT, TEXT):
        if p.resolve() != this:
            text = p.read_text(encoding="utf-8", errors="replace")
            found += [(_rel(p), n, t) for t, n in mentions(text)]
    return code, found, files


def bad_code(code: list[tuple[str, int, str]], areas: set[str]) -> list[str]:
    """the code's strings (file, line, text) in an identifier's form that are none of kind's"""
    return [f"{f}:{n}: {t.replace(HOLE, '{}')!r}" for f, n, t in code if kind(t, areas) is None]


def bad_mentions(
    code: list[tuple[str, int, str]], found: list[tuple[str, int, str]], areas: set[str]
) -> list[str]:
    """the mentions (file, line, text) in the port's form that name no identifier the code gives
    in the style, none of MATLAB_COPIED and no figure's group id"""
    ok = {t for _, _, t in code if kind(t, areas) == "port"} | set(COPIED) | FIGURE_GIDS
    return [f"{f}:{n}: {t!r}" for f, n, t in found if t not in ok]


def test_the_package_gives_its_own_identifiers_in_one_style(scan):
    bad = bad_code(scan[0], package_areas())
    assert not bad, (
        "identifiers neither in the port's style (petakit5d:<area>:<name>, the area a module "
        "of petakit5d) nor MathWorks', MATLAB_COPIED's or NOT_IDENTIFIERS':\n" + "\n".join(bad)
    )
    # the scan reads every folder of code
    assert {f.split("/", 1)[0] for f in scan[2]} == {"src", "tools", "examples", "benchmarks"}


def test_mentions_of_the_ports_identifiers_name_ones_the_package_gives(scan):
    bad = bad_mentions(scan[0], scan[1], package_areas())
    assert not bad, (
        "mentions in the port's form of no identifier the package gives in the style:\n"
        + "\n".join(bad)
    )
    # the scan sees the code, the tests and the docs
    assert {f.split("/", 1)[0] for f, _, _ in scan[1]} >= {"src", "tests", "docs"}


def test_every_listed_string_is_in_the_package(scan):
    code = {t for _, _, t in scan[0]}
    texts = "\n".join(
        p.read_text(encoding="utf-8") for d in CODE for p in _files(PYROOT / d, {".py"})
    )
    stale = [i for i in COPIED if i not in code and i not in texts]
    stale += [s for s in NOT_IDENTIFIERS if s not in code and s not in texts]
    assert not stale, f"listed but not in the package's code: {stale}"


def test_the_copied_identifiers_come_from_their_files():
    matlab = (REPO / "setup.m").is_file() and (REPO / "microscopeDataProcessing").is_dir()
    bad = []
    for src, ids in MATLAB_COPIED.items():
        if src.startswith("python/"):
            path = PYROOT / src.removeprefix("python/")
        elif matlab:
            path = REPO / src
        else:
            continue
        text = path.read_text(encoding="utf-8", errors="replace") if path.is_file() else ""
        bad += [f"{i} is not in {src}" for i in ids if i not in text]
    assert not bad, "\n".join(bad)
    if not matlab:
        pytest.skip("PetaKit5D's MATLAB code is not beside python/ (checked: the harness's)")


AREAS = {"params", "detection", "bioformats", "sim"}


@pytest.mark.parametrize(
    ("ident", "expected"),
    [
        ("petakit5d:params:badValue", "port"),
        ("petakit5d:bioformats:unsupportedAxes", "port"),
        ("petakit5d:sim:naAboveIndex", "port"),
        ("petakit5d:detection:noUsableGPU", "port"),
        ("PetaKit5D:sim:naAboveIndex", None),
        ("PETAKIT5D:params:badValue", None),
        ("pk5d:detection:tooFewPoints", None),
        ("pk5d:regexp", None),
        ("petakit5d:Params:badValue", None),
        ("petakit5d:nosuchmodule:badValue", None),
        ("petakit5d:params:BadValue", None),
        ("petakit5d:params:bad_value", None),
        ("petakit5d:params", None),
        ("petakit5d:params:bad:value", None),
        ("MATLAB:InputParser:UnmatchedParameter", "MathWorks"),
        ("images:imresize3:badScaleOrSize", "MathWorks"),
        ("Images:imresize3:badScaleOrSize", None),
        ("crop:inputError", "copied"),
        ("crop:inputErrors", None),
        ("PetaKit5D:mexFallback", "copied"),
        ("pk5d:otf_xy", "not an identifier"),
        ("no:cacheprovider", "not an identifier"),
        ("generic_single_job_submit_wrapper:taskFile", None),
    ],
)
def test_kind(ident, expected):
    assert kind(ident, AREAS) == expected


def test_code_identifiers_fill_in_module_constants():
    src = (
        "_MEX = 'decon_lucy_rfft_mex'\n"
        "TWICE = 'a'\n"
        "TWICE = 'b'\n"
        "def f(name, x):\n"
        "    raise E(f'{_MEX}:inputError', f'{name} is wrong')\n"
        "    raise E(f'MATLAB:{name}:expectedReal')\n"
        "    raise E(f'custom:{name}', f'{TWICE}:inputError', f'{_MEX!r}:inputError')\n"
        "    raise E(f'petakit5d:{name}:fails', 'crop:inputError', 'a: b', 'not one')\n"
        "    raise E(f'{x:{name}}', 'tile_blending:' 'inputError')\n"
    )
    assert code_identifiers(src) == [
        ("decon_lucy_rfft_mex:inputError", 5),
        (f"MATLAB:{HOLE}:expectedReal", 6),
        (f"petakit5d:{HOLE}:fails", 8),
        ("crop:inputError", 8),
        ("tile_blending:inputError", 9),
    ]
    assert kind(f"petakit5d:{HOLE}:fails", AREAS) is None


def test_mentions_find_the_ports_form_in_any_case():
    text = (
        "``PETAKIT5D:params:badValue``\n"
        "'pk5d:regexp' and Pk5d:otf_xy, petakit5d:sim:naAboveIndex.\n"
        "petakit5d.cli.main:main, pk5d: error, PetaKit5D:mexFallback\n"
    )
    assert mentions(text) == [
        ("PETAKIT5D:params:badValue", 1),
        ("pk5d:regexp", 2),
        ("Pk5d:otf_xy", 2),
        ("petakit5d:sim:naAboveIndex", 2),
        ("PetaKit5D:mexFallback", 3),
    ]


def test_bad_code_and_mentions():
    code = [
        ("src/a.py", 1, "petakit5d:params:badValue"),
        ("src/a.py", 2, f"petakit5d:params:{HOLE}"),
        ("src/b.py", 3, "crop:inputError"),
        ("src/b.py", 4, "crop:outputError"),
        ("tools/c.py", 5, "pk5d:otf_xy"),
    ]
    assert bad_code(code, AREAS) == [
        "src/a.py:2: 'petakit5d:params:{}'",
        "src/b.py:4: 'crop:outputError'",
    ]
    found = [
        ("docs/a.md", 1, "petakit5d:params:badValue"),
        ("docs/a.md", 2, "petakit5d:params:badValues"),
        ("tests/t.py", 3, "PETAKIT5D:params:badValue"),
        ("src/d.py", 4, "PetaKit5D:mexFallback"),
        ("tests/t.py", 5, "pk5d:otf_xy"),
        ("tests/t.py", 6, "pk5d:regexp"),
        ("tests/t.py", 7, "petakit5d:sim:naAboveIndex"),
    ]
    assert bad_mentions(code, found, AREAS) == [
        "docs/a.md:2: 'petakit5d:params:badValues'",
        "tests/t.py:3: 'PETAKIT5D:params:badValue'",
        "tests/t.py:6: 'pk5d:regexp'",
        "tests/t.py:7: 'petakit5d:sim:naAboveIndex'",
    ]
