"""Oracle harness: golden store, MATLAB packs, repo manifests and the pytest check.

In a test module::

    from petakit5d.testing import oracle
    from petakit5d.testing.cases import Case, In, Value

    CASES = [Case('sum_dim3', matlab='sum', args=[In('x'), 3], python=my_sum,
                  inputs={'x': Value(lambda: ...)})]

    @oracle.parametrize(CASES)
    def test_matlab(case, tmp_path):
        oracle.check(case, tmp_path)

``check`` compares the Python result with MATLAB's: with the golden store
(``PETAKIT5D_GOLDEN_DIR``) against the full reference, else against the module's manifest
(``<module>.oracle.json``: hashes, sizes and small arrays), as CI does. References are made on
the cluster by ``python -m petakit5d.testing.oracle run`` (through ``orun.sh``), which runs every
case missing from the store in one MATLAB session and, with ``--manifest``, writes the
manifests. Details: ``python/docs/testing.md``.
"""

from __future__ import annotations

import argparse
import contextlib
import dataclasses
import inspect
import json
import os
import re
import shutil
import socket
import statistics
import subprocess
import sys
import tempfile
import time
import uuid
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from petakit5d.testing.cases import (
    FORMAT,
    Case,
    Data,
    File,
    PyResult,
    Tree,
    Value,
    role_hints,
    role_kind,
)
from petakit5d.testing.compare import (
    NO_MASKS,
    FileTree,
    Masks,
    MValue,
    Report,
    _cmp,
    compare_trees,
    describe,
    describe_tree,
    from_descriptor,
    layout,
    mask_value,
    read_pkv,
    scan_tree,
    sha256_file,
    to_mvalue,
    tree_from_descriptor,
    uuid_masks,
    write_pkv,
)

ENV_GOLDEN = "PETAKIT5D_GOLDEN_DIR"  # the golden store (cluster only)
ENV_ORACLE = "PETAKIT5D_ORACLE_DIR"  # the oracle tree (an immutable copy of the reference commit)
ENV_COMMIT = "PETAKIT5D_ORACLE_COMMIT"  # override of the oracle commit
ENV_MATLAB = "PETAKIT5D_MATLAB"  # MATLAB command (default: matlab)
ENV_MISSING = "PETAKIT5D_ORACLE_MISSING"  # 'skip': skip cases without a reference (default: fail)
ENV_REPORT = "PETAKIT5D_ORACLE_REPORT_DIR"  # write a JSON report per compared case here
ENV_HARNESS = "PETAKIT5D_ORACLE_HARNESS"  # the MATLAB side (default: python/tests/oracle_matlab)
ENV_WORK = "PETAKIT5D_ORACLE_WORK"  # MATLAB work folders (default: the temp folder)

_REGISTRY: dict[str, list[Case]] = {}  # module file -> its cases (filled by parametrize)
_ABS_PATH = re.compile(r'(?<![\w<>.~-])/(?:[^\s/\'":;,()\[\]{}<>]+/)+[^\s/\'":;,()\[\]{}<>]*')


# ---------------------------------------------------------------------------------------------
# configuration
# ---------------------------------------------------------------------------------------------


def harness_dir() -> Path:
    """python/tests/oracle_matlab: the MATLAB package and the repo's oracle commit"""
    p = os.environ.get(ENV_HARNESS)
    return Path(p) if p else Path(__file__).resolve().parents[3] / "tests" / "oracle_matlab"


def python_root() -> Path:
    """the repo's python/ folder, which the cases' helper folders and the manifests' module
    names are relative to: the harness's (python/tests/oracle_matlab), since an installed
    petakit5d (site-packages, as in CI) is not inside it"""
    return harness_dir().resolve().parents[1]


def expected_commit() -> str:
    """the oracle commit the repo's manifests are made with (tests/oracle_matlab/ORACLE_COMMIT)"""
    return (harness_dir() / "ORACLE_COMMIT").read_text().split()[0].strip()


def read_oracle_commit(tree: str | os.PathLike) -> str:
    """the commit of an oracle tree, from its .git folder (no git needed)"""
    g = Path(tree) / ".git"
    head = (g / "HEAD").read_text().strip()
    if head.startswith("ref:"):
        ref = head[4:].strip()
        p = g / ref
        if p.is_file():
            return p.read_text().strip()
        for line in (g / "packed-refs").read_text().splitlines():
            if line.endswith(" " + ref):
                return line.split()[0]
        raise ValueError(f"{tree}: cannot resolve {ref}")
    return head


def oracle_dir() -> Path | None:
    p = os.environ.get(ENV_ORACLE)
    return Path(p) if p else None


def current_commit() -> str:
    c = os.environ.get(ENV_COMMIT)
    if c:
        return c
    d = oracle_dir()
    if d is not None and (d / ".git").exists():
        return read_oracle_commit(d)
    return expected_commit()


def short(commit: str) -> str:
    return commit[:7]


def same_commit(a: str, b: str) -> bool:
    n = min(len(a), len(b))
    return n >= 7 and a[:n] == b[:n]


# ---------------------------------------------------------------------------------------------
# golden store
# ---------------------------------------------------------------------------------------------


def _replace_same(tmp: Path, p: Path) -> None:
    """os.replace(tmp, p) of a file named by its content: Windows refuses the replace while
    another writer replaces or reads p (PermissionError); p then holds the same bytes"""
    try:
        os.replace(tmp, p)
    except PermissionError:
        if not p.is_file():
            raise
        with contextlib.suppress(OSError):
            os.remove(tmp)


class Store:
    """the shared store of MATLAB results, written only by publishing folders by rename

    <root>/inputs/<hh>/<hash>.pkv|.file|.tree   inputs, by content hash
    <root>/<commit7>/<hh>/<key>/                one case's result (meta.json, outputs.pkv,
                                                files/, decoded/, stdout.txt)
    <root>/<commit7>/.staging/<pack>/<key>/     results being written (never read)
    """

    def __init__(self, root: str | os.PathLike, commit: str):
        self.root = Path(root)
        self.commit = commit
        self.dir = self.root / short(commit)

    @classmethod
    def from_env(cls) -> Store | None:
        root = os.environ.get(ENV_GOLDEN)
        if not root or not Path(root).is_dir():
            return None
        return cls(root, current_commit())

    def case_dir(self, key: str) -> Path:
        return self.dir / key[:2] / key

    def has(self, key: str) -> bool:
        return (self.case_dir(key) / "meta.json").is_file()

    def staging(self, pack_id: str) -> Path:
        return self.dir / ".staging" / pack_id

    def load(self, key: str) -> Reference:
        return Reference.from_store(self, key)

    # inputs

    def input_path(self, h: str, ext: str) -> Path:
        return self.root / "inputs" / h[:2] / f"{h}.{ext}"

    def _tmp(self) -> Path:
        t = self.root / "inputs" / ".tmp"
        t.mkdir(parents=True, exist_ok=True)
        return t / f"{uuid.uuid4().hex}"

    def put_value(self, v: MValue, h: str) -> Path:
        p = self.input_path(h, "pkv")
        if not p.is_file():
            p.parent.mkdir(parents=True, exist_ok=True)
            tmp = self._tmp()
            write_pkv(tmp, v)
            _replace_same(tmp, p)
        return p

    def put_file(self, src: Path, h: str) -> Path:
        p = self.input_path(h, "file")
        if not p.is_file():
            p.parent.mkdir(parents=True, exist_ok=True)
            tmp = self._tmp()
            shutil.copyfile(src, tmp)
            _replace_same(tmp, p)
        return p

    def put_tree(self, src: Path, h: str) -> Path:
        p = self.input_path(h, "tree")
        if not p.is_dir():
            p.parent.mkdir(parents=True, exist_ok=True)
            tmp = self._tmp()
            shutil.copytree(src, tmp)
            try:
                os.rename(tmp, p)
            except OSError:
                if not p.is_dir():
                    raise
                shutil.rmtree(tmp, ignore_errors=True)
        return p

    def put_input(self, x: Any) -> Path:
        if isinstance(x, Value):
            return self.put_value(x.mvalue(), x.hash())
        if isinstance(x, File):
            return self.put_file(x.local(), x.hash())
        if isinstance(x, Tree):
            return self.put_tree(x.local(), x.hash())
        raise TypeError(type(x).__name__)


_WARN_BLOCK = re.compile(r"\[\x08Warning: (.*?)\]\x08", re.S)
_WARN_BLOCK_PLAIN = re.compile(r"\[Warning: (.*?)\]\s*(?:\n|$)", re.S)
_WARN_TYPE = re.compile(
    r'^(.*?)\s*\(Type\s+"(?:<a href="[^"]*">)?warning\s+off\s+([^"<]*?)(?:</a>)?"\s+'
    r"to\s+suppress\s+this\s+warning\.\)\s*$",
    re.S,
)
UNKNOWN_IDS = ("", "pkpy:oracle:unparsedWarning")


def parse_matlab_warnings(
    text: str, last_id: str | None = None, last_msg: str | None = None
) -> list[dict]:
    """the warnings in MATLAB's captured command window text (as +pkpy_oracle/parse_warnings.m)

    MATLAB writes ``[<backspace>Warning: message]<backspace>``, with the line ``(Type "warning
    off <id>" to suppress this warning.)`` before the ``]`` unless the warning was issued with
    verbose off; long lines are wrapped at the window width, also inside that line. A warning
    without the line has an empty identifier. The last warning takes ``last_id`` and the exact
    ``last_msg`` (MATLAB's lastwarn) when its identifier is ``last_id`` or its message is
    ``last_msg`` but for line breaks.
    """
    blocks = _WARN_BLOCK.findall(text)
    if not blocks and "[Warning: " in text:
        blocks = _WARN_BLOCK_PLAIN.findall(text)
    out = []
    for body in blocks:
        m = _WARN_TYPE.match(body)
        if m:
            out.append({"identifier": m.group(2).strip(), "message": m.group(1)})
        else:
            out.append({"identifier": "", "message": body.strip()})
    if out and last_msg is not None:
        ident = out[-1]["identifier"]
        same = " ".join(out[-1]["message"].split()) == " ".join(last_msg.split())
        if (ident == (last_id or "") and (last_id or same)) or (not ident and same):
            out[-1] = {"identifier": last_id or ident, "message": last_msg}
    n = text.count("[Warning: ") + text.count("[\x08Warning: ")
    if n > len(blocks):
        out.append(
            {
                "identifier": "pkpy:oracle:unparsedWarning",
                "message": f"{n - len(blocks)} warnings could not be parsed; see stdout.txt",
            }
        )
    return out


def _seq(v: Any) -> list:
    if v is None:
        return []
    return v if isinstance(v, list) else [v]


@dataclasses.dataclass
class Reference:
    """MATLAB's result for one case, from the store (full) or from a manifest (descriptors)"""

    key: str
    status: str  # 'ok' or 'error'
    error: dict | None
    warnings: list[dict]
    outputs: list[MValue] | None
    files: FileTree | None
    meta: dict
    source: str  # 'store' or 'manifest'
    masks: Masks
    commit: str

    @classmethod
    def from_store(cls, store: Store, key: str) -> Reference:
        d = store.case_dir(key)
        meta = json.loads((d / "meta.json").read_text(encoding="utf-8"))
        outputs = None
        if meta["status"] == "ok":
            v = read_pkv(d / "outputs.pkv")
            outputs = list(v.items)
        masks = Masks().with_paths(work=meta.get("work"), host=meta.get("host"))
        files = None
        if meta.get("capture"):
            files = scan_tree(
                d / "files",
                _seq(meta.get("files")),
                decoded_root=d / "decoded",
                deleted=_seq(meta.get("deleted")),
            )
        err = meta.get("error") if meta["status"] == "error" else None
        warns = [w for w in _seq(meta.get("warnings")) if isinstance(w, dict)]
        so = d / "stdout.txt"
        if so.is_file():
            # parsed again from the command window text, so that results made by an older
            # runner get today's rules; the last warning's exact text is in the meta
            lw = meta.get("last_warning")
            lid = lmsg = None
            if isinstance(lw, dict):
                lid, lmsg = lw.get("identifier", ""), lw.get("message", "")
            elif warns and warns[-1].get("identifier") != UNKNOWN_IDS[1]:
                # made before last_warning was recorded: the runner had the last one exact
                lid, lmsg = warns[-1].get("identifier", ""), warns[-1].get("message", "")
            text = so.read_text(encoding="utf-8", errors="replace")
            warns = parse_matlab_warnings(text, lid, lmsg)
        return cls(
            key,
            meta["status"],
            err,
            warns,
            outputs,
            files,
            meta,
            "store",
            masks,
            meta.get("oracle_commit", store.commit),
        )

    @classmethod
    def from_manifest(cls, entry: dict) -> Reference:
        outs = (
            None if entry.get("outputs") is None else [from_descriptor(d) for d in entry["outputs"]]
        )
        files = None if entry.get("files") is None else tree_from_descriptor(entry["files"])
        return cls(
            entry["key"],
            entry["status"],
            entry.get("error"),
            entry.get("warnings") or [],
            outs,
            files,
            entry,
            "manifest",
            NO_MASKS,
            entry["oracle"],
        )


# ---------------------------------------------------------------------------------------------
# manifests (in the repo, next to each test module)
# ---------------------------------------------------------------------------------------------


def manifest_path(case: Case) -> Path:
    g = case._module_globals or {}
    if g.get("ORACLE_MANIFEST"):
        return Path(g["ORACLE_MANIFEST"])
    f = g.get("__file__")
    if not f:
        raise ValueError(f"case {case.id} was not declared with oracle.parametrize in a module")
    return Path(f).with_suffix(".oracle.json")


_MANIFESTS: dict[Path, tuple[bytes, dict]] = {}


def load_manifest(path: Path) -> dict:
    """{case id: entry}; parsed again only when the file's bytes change (not by time stamp:
    two writes within the clock's resolution keep it)"""
    path = Path(path)
    if not path.is_file():
        return {}
    data = path.read_bytes()
    c = _MANIFESTS.get(path)
    if c is None or c[0] != data:
        try:
            cases = json.loads(data.decode("utf-8")).get("cases", {})
        except ValueError as e:
            msg = f"{path}: not a valid manifest ({e}); rerun orun.sh --manifest"
            raise ValueError(msg) from None
        c = (data, cases)
        _MANIFESTS[path] = c
    return c[1]


def _clean_text(s: str, masks: Masks) -> str:
    """text for the repo: masks, then any remaining absolute path as <PATH>"""
    return _ABS_PATH.sub("<PATH>", masks.text(s or ""))


def manifest_entry(case: Case, ref: Reference) -> dict:
    """the manifest entry of a stored reference: no paths, hosts or times"""
    m = ref.masks.with_paths(oracle=str(oracle_dir() or ""), golden=os.environ.get(ENV_GOLDEN, ""))
    outs = None
    if ref.outputs is not None:
        outs = [describe(mask_value(v, ref.masks)) for v in ref.outputs]
    files = None
    if ref.files is not None:
        files = describe_tree(ref.files, case.tree, ref.masks.plus(*case.tree.masks))
    err = None
    if ref.error:
        err = {
            "identifier": ref.error.get("identifier", ""),
            "message": _clean_text(ref.error.get("message", ""), m),
        }
    return {
        "key": ref.key,
        "oracle": ref.commit,
        "format": FORMAT,
        "matlab": case.matlab or "(code)",
        "status": ref.status,
        "error": err,
        "warnings": [
            {"identifier": w.get("identifier", ""), "message": _clean_text(w.get("message", ""), m)}
            for w in ref.warnings
        ],
        "outputs": outs,
        "files": files,
    }


def write_manifest(path: Path, entries: dict[str, dict], module: str) -> bool:
    """write a manifest (sorted, stable; one compact line per case, so that a changed case is
    one changed line); True if it changed"""

    def dump(x: Any) -> str:
        return json.dumps(x, sort_keys=True, separators=(",", ":"), ensure_ascii=True)

    cases = [f"  {dump(cid)}: {dump(entries[cid])}" for cid in sorted(entries)]
    head = [
        f' "format": {FORMAT}',
        f' "module": {dump(module)}',
        f' "oracle": {dump(expected_commit())}',
    ]
    text = "{\n" + ",\n".join(head) + ',\n "cases": {\n' + ",\n".join(cases) + "\n }\n}\n"
    if path.is_file() and path.read_text(encoding="utf-8") == text:
        return False
    tmp = path.with_name(path.name + f".tmp{os.getpid()}")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)
    return True


def check_manifests(root: str | os.PathLike | None = None) -> list[str]:
    """problems of the manifests under root (default python/tests): entries made with another
    oracle commit or another case format"""
    root = Path(root) if root else harness_dir().parent
    exp = expected_commit()
    out = []
    for p in sorted(root.rglob("*.oracle.json")):
        try:
            doc = json.loads(p.read_text(encoding="utf-8"))
        except ValueError as e:
            out.append(f"{p}: not valid JSON ({e})")
            continue
        for cid, e in sorted(doc.get("cases", {}).items()):
            if not same_commit(e.get("oracle", ""), exp):
                out.append(
                    f"{p.name}: {cid}: reference made with oracle {short(e.get('oracle', '?'))}, "
                    f"the repo's oracle is {short(exp)}"
                )
            if e.get("format") != FORMAT:
                out.append(
                    f"{p.name}: {cid}: case format {e.get('format')}, the harness has {FORMAT}"
                )
    return out


# ---------------------------------------------------------------------------------------------
# comparing a Python run with a reference
# ---------------------------------------------------------------------------------------------


def compare_case(case: Case, ref: Reference, res: PyResult) -> Report:
    rep = Report()
    am = Masks().with_paths(work=str(res.work), host=socket.gethostname())
    am = am.plus(*uuid_masks(res.uuids))
    rm = ref.masks
    if ref.status == "error":
        rid = (ref.error or {}).get("identifier", "")
        rmsg = (ref.error or {}).get("message", "")
        if res.error is None:
            rep.add(
                "error",
                "error",
                f"MATLAB raised {rid or 'an error'}: {rmsg}",
                "Python returned normally",
            )
        else:
            pid = getattr(res.error, "identifier", None)
            if not isinstance(res.error, case.py_error):
                rep.add(
                    "error",
                    "error",
                    f"Python raised {type(res.error).__name__}, the case allows {case.py_error}",
                )
            if rid and pid is not None and pid != rid:
                rep.add("error", "identifier", f"MATLAB {rid}, Python {pid}")
            elif rid and pid is None:
                rep.note(
                    "error",
                    f"MATLAB raised {rid}; the Python exception has no identifier attribute",
                )
            if case.error_message and _clean_text(str(res.error), am) != _clean_text(rmsg, rm):
                rep.add("error", "message", f"MATLAB {rmsg!r}", f"Python {str(res.error)!r}")
    else:
        if res.error is not None:
            tail = res.error_text.strip().splitlines()[-6:]
            rep.add(
                "error",
                "error",
                f"Python raised {type(res.error).__name__}: {res.error}; MATLAB returned normally",
                *tail,
            )
        else:
            names = case.output_names()
            roles = list(case.out_roles) if case.out_roles is not None else [None] * len(names)
            for i, name in enumerate(names):
                r = mask_value(ref.outputs[i], rm) if ref.source == "store" else ref.outputs[i]
                try:
                    a = mask_value(to_mvalue(res.outputs[i], like=r), am)
                except (TypeError, ValueError) as e:
                    rep.add(name, "kind", f"the Python value cannot be read as a MATLAB value: {e}")
                    continue
                sub = Report()
                with layout("vol" if role_kind(roles[i]) == "vol" else None):
                    _cmp(r, a, case.spec_for(i), name, sub)
                if sub.problems:
                    try:
                        hints = role_hints(case, roles[i], r, a)
                    except Exception:  # a hint never stops a report
                        hints = []
                    sub.problems[0].lines.extend(f"hint: {h}" for h in hints)
                rep.extend(sub)
    if case.captures_files and ref.files is not None and res.files is not None:
        rep.extend(compare_trees(ref.files, res.files, case.tree, ref_masks=rm, act_masks=am))
    if case.warnings == "ids":
        rw = sorted(w.get("identifier", "") for w in ref.warnings)
        pw = sorted(w.get("identifier", "") for w in res.warnings)
        if rw != pw:
            rep.add("warnings", "warnings", f"warning identifiers: MATLAB {rw}, Python {pw}")
    return rep


def _title(case: Case, ref: Reference) -> str:
    what = case.matlab or "MATLAB code"
    return (
        f"case {case.id} ({what}; reference from the {ref.source}, oracle {short(ref.commit)}, key "
        f"{ref.key[:12]})"
    )


# ---------------------------------------------------------------------------------------------
# pytest
# ---------------------------------------------------------------------------------------------


def parametrize(cases: Sequence[Case]):
    """pytest decorator: one test per case (id = case id), marked ``oracle`` and with the
    case's marks; the cases are registered for ``oracle run`` and their manifest"""
    import pytest

    g = inspect.currentframe().f_back.f_globals
    ids = [c.id for c in cases]
    dup = sorted({i for i in ids if ids.count(i) > 1})
    if dup:
        raise ValueError(f"duplicate case ids: {dup}")
    reg = _REGISTRY.setdefault(g.get("__file__", "?"), [])
    for c in cases:
        c._module_globals = g
        if all(c is not x for x in reg):
            if any(x.id == c.id for x in reg):
                raise ValueError(f"two different cases have id {c.id!r} in {g.get('__file__')}")
            reg.append(c)
    params = [
        pytest.param(c, id=c.id, marks=[getattr(pytest.mark, m) for m in c.all_marks()])
        for c in cases
    ]

    def deco(fn):
        return pytest.mark.oracle(pytest.mark.parametrize("case", params)(fn))

    return deco


def reference(case: Case) -> Reference:
    """the reference of a case: the store if it has it, else the manifest (checked for the
    oracle commit and the case key); fails or skips the test when there is none"""
    import pytest

    exp = expected_commit()
    key = case.key()
    mpath = manifest_path(case)
    entry = load_manifest(mpath).get(case.id)
    missing = os.environ.get(ENV_MISSING, "fail")
    store = Store.from_env()
    if store is not None:
        if not same_commit(store.commit, exp):
            pytest.fail(
                f"{ENV_ORACLE} is oracle {short(store.commit)}; the repo expects {short(exp)} "
                f"(tests/oracle_matlab/ORACLE_COMMIT)",
                pytrace=False,
            )
        if store.has(key):
            if entry is None or entry.get("key") != key:
                import warnings

                warnings.warn(
                    f"case {case.id}: manifest {mpath.name} lacks this reference or is out of "
                    "date; "
                    "run orun.sh --manifest",
                    stacklevel=2,
                )
            return store.load(key)
    if entry is None:
        msg = (
            f"case {case.id}: no MATLAB reference (key {key[:12]}): run orun.sh, then orun.sh "
            "--manifest"
        )
        if missing == "skip":
            pytest.skip(msg)
        pytest.fail(msg, pytrace=False)
    if not same_commit(entry.get("oracle", ""), exp):
        pytest.fail(
            f"case {case.id}: reference made with oracle {short(entry.get('oracle', '?'))}, the "
            "repo's "
            f"oracle is {short(exp)}: rerun orun.sh --manifest and review what changed",
            pytrace=False,
        )
    if entry.get("key") != key:
        msg = (
            f"case {case.id}: the case or its inputs changed since its reference was made (key "
            f"{key[:12]}, "
            f"manifest {str(entry.get('key'))[:12]}): rerun orun.sh --manifest. If nothing "
            "changed, an input "
            "generator gives other bytes here (orun.sh checks that they are deterministic)"
        )
        if missing == "skip" and store is not None:
            pytest.skip(msg)
        pytest.fail(msg, pytrace=False)
    return Reference.from_manifest(entry)


def check(case: Case, tmp_path: str | os.PathLike) -> Report:
    """pytest check of one case: run the Python side and compare it with MATLAB's result"""
    import pytest

    ref = reference(case)
    if case.expect_error and ref.status != "error":
        pytest.fail(
            f"case {case.id}: expects an error, but MATLAB returned normally", pytrace=False
        )
    if (
        isinstance(case.expect_error, str)
        and ref.error
        and ref.error.get("identifier") != case.expect_error
    ):
        pytest.fail(
            f"case {case.id}: expects error {case.expect_error}, MATLAB raised "
            f"{ref.error.get('identifier')!r}: {ref.error.get('message')}",
            pytrace=False,
        )
    if not case.expect_error and ref.status == "error":
        e = ref.error or {}
        pytest.fail(
            f"case {case.id}: MATLAB raised {e.get('identifier')!r}: {e.get('message')} "
            "(set expect_error=True if that is what the case tests)",
            pytrace=False,
        )
    if case.python is None:
        return Report()
    res = case.run_python(tmp_path)
    rep = compare_case(case, ref, res)
    out = os.environ.get(ENV_REPORT)
    if out:
        Path(out).mkdir(parents=True, exist_ok=True)
        (Path(out) / f"{case.id}.json").write_text(
            json.dumps(
                {"case": case.id, "key": ref.key, "source": ref.source, **rep.to_json()}, indent=1
            )
        )
    if not rep.ok:
        pytest.fail(rep.text(_title(case, ref)), pytrace=False)
    if rep.undecided:
        pytest.skip(
            f"case {case.id}: not decided with the manifest alone: " + "; ".join(rep.undecided[:3])
        )
    return rep


# ---------------------------------------------------------------------------------------------
# packs: running cases in MATLAB
# ---------------------------------------------------------------------------------------------


def _q(p: Any) -> str:
    return str(p).replace("'", "''")


def write_pack(cases: Sequence[Case], store: Store, pack_dir: Path, work_root: Path) -> Path:
    """write a pack (pack.json, argument files, the MATLAB side); inputs go to the store"""
    pack_dir.mkdir(parents=True, exist_ok=True)
    (pack_dir / "args").mkdir(exist_ok=True)
    mat = pack_dir / "matlab"
    if not mat.exists():
        shutil.copytree(harness_dir() / "+pkpy_oracle", mat / "+pkpy_oracle")
    pack_id = f"{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}_{uuid.uuid4().hex[:8]}"
    helpers: list[str] = []
    entries = []
    for c in cases:
        key = c.key()
        paths: dict[str, str] = {}

        def input_path(name: str, x: Any, paths: dict = paths) -> str:
            if name not in paths:
                paths[name] = str(store.put_input(x))
            return paths[name]

        args_fn = pack_dir / "args" / f"{key}.pkv"
        write_pkv(args_fn, c.matlab_args(input_path))
        copies = []
        for _name, x in sorted(c.inputs.items()):
            if isinstance(x, (File, Tree)):
                copies.append({"src": str(store.put_input(x)), "dst": x.name, "kind": x.kind})
            elif isinstance(x, Data) and x.copy:
                copies.append(
                    {
                        "src": str(x.local()),
                        "dst": x.name,
                        "kind": "tree" if x.local().is_dir() else "file",
                    }
                )
        for h in c.helpers:
            src = python_root() / h
            dst = pack_dir / "helpers" / h.replace("/", "__")
            if not dst.exists():
                shutil.copytree(src, dst)
            if str(dst) not in helpers:
                helpers.append(str(dst))
        entries.append(
            {
                "key": key,
                "id": c.id,
                "func": c.matlab or "",
                "code": c.code or "",
                "nargout": c.nargout,
                "args": str(args_fn),
                "capture": bool(c.captures_files),
                "copy": copies,
                "device": c.device,
            }
        )
    pack = {
        "format": FORMAT,
        "pack": pack_id,
        "oracle_commit": store.commit,
        "oracle_dir": str(oracle_dir() or ""),
        "store": str(store.dir),
        "staging": str(store.staging(pack_id)),
        "work_root": str(work_root),
        "helpers": helpers,
        "summary": str(pack_dir / "summary.json"),
        "cases": entries,
    }
    fn = pack_dir / "pack.json"
    fn.write_text(json.dumps(pack, indent=1), encoding="utf-8")
    return fn


def matlab_command(pack_json: Path, oracle: Path) -> str:
    pack = json.loads(pack_json.read_text(encoding="utf-8"))
    s = f"cd('{_q(oracle)}'); setup([]); addpath('{_q(pack_json.parent / 'matlab')}'); "
    s += "".join(f"addpath('{_q(h)}'); " for h in pack["helpers"])
    return s + f"pkpy_oracle.run_pack('{_q(pack_json)}');"


def run_matlab(
    pack_json: Path,
    oracle: Path,
    *,
    matlab: str | None = None,
    timeout: float | None = None,
    echo: bool = True,
) -> tuple[int, float]:
    """run a pack in one MATLAB session; returns (exit code, wall time)"""
    matlab = matlab or os.environ.get(ENV_MATLAB, "matlab")
    env = dict(os.environ)
    env.pop("DISPLAY", None)
    log = pack_json.parent / "matlab.log"
    t0 = time.perf_counter()
    with open(log, "w", encoding="utf-8") as lf:
        p = subprocess.Popen(
            [*matlab.split(), "-batch", matlab_command(pack_json, oracle)],
            cwd=oracle,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            errors="replace",
        )
        try:
            for line in p.stdout:
                if "does not support locale" in line:
                    continue
                lf.write(line)
                lf.flush()
                if echo and (line.startswith("PKPY") or "Error" in line):
                    print(line.rstrip(), flush=True)
                if timeout is not None and time.perf_counter() - t0 > timeout:
                    p.terminate()
                    break
            rc = p.wait()
        finally:
            if p.poll() is None:
                p.terminate()
                p.wait()
    return rc, time.perf_counter() - t0


# ---------------------------------------------------------------------------------------------
# command line
# ---------------------------------------------------------------------------------------------


def collect(pytest_args: Sequence[str]) -> list[tuple[str, Case]]:
    """the oracle cases pytest selects with these arguments"""
    import pytest

    class _Plugin:
        def __init__(self):
            self.items: list[tuple[str, Case]] = []

        def pytest_collection_finish(self, session):
            for it in session.items:
                cs = getattr(it, "callspec", None)
                c = cs.params.get("case") if cs is not None else None
                if isinstance(c, Case):
                    self.items.append((it.nodeid, c))

    p = _Plugin()
    rc = pytest.main(["--collect-only", "-qq", "-p", "no:cacheprovider", *pytest_args], plugins=[p])
    if rc not in (0, 5):
        raise SystemExit(f"pytest collection failed (exit {rc})")
    return p.items


def _store_and_oracle(allow_other: bool) -> tuple[Store, Path]:
    root, od = os.environ.get(ENV_GOLDEN), oracle_dir()
    if not root or od is None:
        raise SystemExit(f"set {ENV_GOLDEN} and {ENV_ORACLE} (orun.sh does)")
    commit = current_commit()
    exp = expected_commit()
    if not same_commit(commit, exp) and not allow_other:
        raise SystemExit(
            f"the oracle tree is {short(commit)}, the repo expects {short(exp)} "
            "(tests/oracle_matlab/ORACLE_COMMIT); --allow-other-oracle to run anyway"
        )
    return Store(root, commit), od


def _device_now() -> str:
    return "gpu" if os.environ.get("CUDA_VISIBLE_DEVICES", "").strip() not in ("", "-1") else "cpu"


def _write_manifests(
    items: list[tuple[str, Case]], store: Store, pack_dir: Path | None
) -> list[str]:
    files = sorted({c.module_file for _, c in items if c.module_file})
    written = []
    py_root = python_root()
    for f in files:
        cases = _REGISTRY.get(f, [])
        if not cases:
            continue
        path = manifest_path(cases[0])
        try:
            old = load_manifest(path)
        except ValueError as e:
            print(f"{e}: writing it anew")
            old = {}
        entries, lacking = {}, []
        for c in cases:
            k = c.key()
            if store.has(k):
                entries[c.id] = manifest_entry(c, store.load(k))
            else:
                lacking.append(c.id)
                if c.id in old:
                    entries[c.id] = old[c.id]
        try:
            module = Path(f).resolve().relative_to(py_root).as_posix()
        except ValueError:
            module = Path(f).name
        if write_manifest(path, entries, module):
            try:
                written.append(path.resolve().relative_to(py_root).as_posix())
            except ValueError:
                written.append(str(path))
        print(
            f"manifest {path.name}: {len(entries)} cases"
            + (f"; no reference yet (old entries kept): {lacking}" if lacking else "")
        )
    if pack_dir is not None:
        pack_dir.mkdir(parents=True, exist_ok=True)
        # with each file's SHA-256: orun.sh checks its copy against it (network file systems
        # can show another machine a file before all its bytes)
        py = python_root()
        lines = [f"{sha256_file(py / w)}  {w}\n" for w in written]  # py / an absolute w is w
        (pack_dir / "manifests.txt").write_text("".join(lines))
    return written


def _summarize(pack_dir: Path, wall: float) -> None:
    fn = pack_dir / "summary.json"
    if not fn.is_file():
        print("MATLAB wrote no summary (it crashed or was stopped); see", pack_dir / "matlab.log")
        return
    s = json.loads(fn.read_text(encoding="utf-8"))
    rows = [r for r in _seq(s.get("cases")) if isinstance(r, dict)]
    ran = [r for r in rows if r.get("status") not in ("cached",)]
    over = [
        1000 * (r["t_case"] - r["t_call"])
        for r in ran
        if isinstance(r.get("t_call"), (int, float)) and isinstance(r.get("t_case"), (int, float))
    ]
    print(
        f"MATLAB {s.get('matlab_version', '?')} on {s.get('threads', '?')} threads: process "
        f"{wall:.1f} s, "
        f"pack {s.get('t_total', 0):.1f} s, start-up and setup {wall - s.get('t_total', 0):.1f} s"
    )
    if over:
        print(
            f"harness overhead per case beyond the call: median {statistics.median(over):.1f} ms, "
            f"mean {statistics.fmean(over):.1f} ms, max {max(over):.1f} ms ({len(over)} cases)"
        )


def cmd_run(ns: argparse.Namespace) -> int:
    store, od = _store_and_oracle(ns.allow_other_oracle)
    items = collect(ns.pytest_args)
    dev = ns.device or _device_now()
    seen: dict[str, Case] = {}
    skipped_dev = []
    for _, c in items:
        if c.device.split(":")[0] != dev:
            skipped_dev.append(c.id)
            continue
        seen.setdefault(c.key(), c)
    bad = [p for c in seen.values() for p in c.check_inputs()]
    if bad:
        print("\n".join(bad))
        print("no MATLAB run: fix the input generators first")
        return 2
    missing = [c for k, c in seen.items() if not store.has(k)]
    print(
        f"{len(items)} cases selected, {len(seen)} distinct for this {dev} job; "
        f"{len(seen) - len(missing)} "
        f"already in the store ({short(store.commit)}), {len(missing)} to run"
    )
    if skipped_dev:
        print(
            f"{len(skipped_dev)} cases need a {'GPU' if dev == 'cpu' else 'CPU'} pack: "
            f"{skipped_dev[:10]}"
        )
    pack_dir = Path(ns.pack_dir or tempfile.mkdtemp(prefix="pkpy_pack_"))
    rc = 0
    if missing and not ns.dry_run:
        work_root = (
            Path(os.environ.get(ENV_WORK) or tempfile.gettempdir())
            / f"pkpy_work_{os.getpid()}_{uuid.uuid4().hex[:6]}"
        )
        work_root.mkdir(parents=True)
        try:
            pj = write_pack(sorted(missing, key=lambda c: c.id), store, pack_dir, work_root)
            mrc, wall = run_matlab(pj, od, timeout=ns.timeout)
            print(f"MATLAB exit code {mrc}; log {pack_dir / 'matlab.log'}")
            _summarize(pack_dir, wall)
        finally:
            shutil.rmtree(work_root, ignore_errors=True)
    still = [c for k, c in seen.items() if not store.has(k)]
    errs = []
    for k, c in seen.items():
        if store.has(k):
            meta = json.loads((store.case_dir(k) / "meta.json").read_text(encoding="utf-8"))
            if meta["status"] == "error" and not c.expect_error:
                errs.append(
                    f"{c.id}: MATLAB raised {meta.get('error', {}).get('identifier')!r}: "
                    f"{meta.get('error', {}).get('message')}"
                )
    if still:
        rc = 1
        print(f"{len(still)} cases have no reference after the run: {[c.id for c in still][:20]}")
    for e in errs:
        print("unexpected MATLAB error:", e)
    if ns.manifest:
        _write_manifests(items, store, pack_dir)
    return rc


def cmd_manifest(ns: argparse.Namespace) -> int:
    store, _ = _store_and_oracle(ns.allow_other_oracle)
    items = collect(ns.pytest_args)
    _write_manifests(items, store, Path(ns.pack_dir) if ns.pack_dir else None)
    return 0


def cmd_status(ns: argparse.Namespace) -> int:
    store = Store.from_env()
    items = collect(ns.pytest_args)
    exp = expected_commit()
    n_bad = 0
    for _nodeid, c in items:
        k = c.key()
        st = (
            "in store"
            if store is not None and store.has(k)
            else ("-" if store is not None else "no store")
        )
        e = load_manifest(manifest_path(c)).get(c.id)
        if e is None:
            ms = "no manifest entry"
        elif e.get("key") != k:
            ms = "manifest out of date (case changed)"
        elif not same_commit(e.get("oracle", ""), exp):
            ms = f"manifest from oracle {short(e.get('oracle', '?'))}"
        else:
            ms = "manifest current"
        n_bad += ms != "manifest current"
        print(f"{c.id:40s} {k[:12]}  {st:9s}  {ms}")
    for p in check_manifests():
        print("manifest problem:", p)
    return 1 if n_bad else 0


def cmd_clean(ns: argparse.Namespace) -> int:
    """remove staging folders and temporary inputs older than --hours (left by stopped jobs)"""
    root = os.environ.get(ENV_GOLDEN)
    if not root:
        raise SystemExit(f"set {ENV_GOLDEN}")
    limit = time.time() - ns.hours * 3600
    n = 0
    for d in list(Path(root).glob("*/.staging/*")) + list(Path(root).glob("inputs/.tmp/*")):
        if d.stat().st_mtime < limit:
            shutil.rmtree(d, ignore_errors=True) if d.is_dir() else d.unlink(missing_ok=True)
            n += 1
    print(f"removed {n} leftovers older than {ns.hours} h")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="python -m petakit5d.testing.oracle", description=__doc__.split("\n\n")[0]
    )
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name, fn, hlp in (
        ("run", cmd_run, "run the selected cases missing from the store in one MATLAB session"),
        ("manifest", cmd_manifest, "write the manifests of the selected modules from the store"),
        ("status", cmd_status, "show store and manifest state of the selected cases"),
    ):
        p = sub.add_parser(name, help=hlp)
        p.set_defaults(fn=fn)
        p.add_argument("--pack-dir", help="folder for the pack, its MATLAB log and summary")
        p.add_argument(
            "--allow-other-oracle",
            action="store_true",
            help="run with an oracle tree the repo does not expect",
        )
        p.add_argument("--args-file", help="more pytest arguments, one per line")
        p.add_argument(
            "pytest_args", nargs="*", help="pytest selection (paths, node ids, -k/-m expressions)"
        )
        if name == "run":
            p.add_argument(
                "--manifest",
                action="store_true",
                help="also write the manifests of the selected modules",
            )
            p.add_argument("--dry-run", action="store_true", help="only report what would run")
            p.add_argument(
                "--device",
                choices=("cpu", "gpu"),
                help="device of this job (default: from CUDA_VISIBLE_DEVICES)",
            )
            p.add_argument("--timeout", type=float, help="stop MATLAB after this many seconds")
    p = sub.add_parser("clean-staging", help="remove leftovers of stopped jobs")
    p.set_defaults(fn=cmd_clean)
    p.add_argument("--hours", type=float, default=24.0)
    ns = ap.parse_args(argv)
    if getattr(ns, "args_file", None):
        ns.pytest_args = list(ns.pytest_args) + [
            ln.strip() for ln in Path(ns.args_file).read_text().splitlines() if ln.strip()
        ]
    return ns.fn(ns)


if __name__ == "__main__":
    # run in the imported module, whose registry the test modules fill (not in this __main__ copy)
    from petakit5d.testing import oracle as _oracle

    sys.exit(_oracle.main())
