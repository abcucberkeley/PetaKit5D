"""The publish-order spy (MATLAB's pktest.MovefileSpyFixture) and the duplicate-writer race.

PublishSpy
----------
Tasks publish their result by renaming it into place, and the job framework counts a task as
done (and may cancel its job) once the result exists, so everything else a task writes must
exist before that rename. ``PublishSpy`` wraps ``os.rename``, ``os.replace`` and ``os.link``
(the publish of a file on Linux and macOS) in this process: it records each call with the
files and folders next to its destination at that moment (``list_files``), then does the
call as usual::

    with PublishSpy() as spy:
        run_the_task(...)
    assert "MIPs/result_MIP_z.tif" in spy.files_when_moved(result)

A pytest fixture is two lines: ``with PublishSpy() as spy: yield spy``. Only calls of this
process are seen (tasks run locally); calls are recorded before they run, so a rename that
fails is recorded too (``ok`` says how it went).

The race
--------
``race_round`` lets ``n`` writer processes make the same result at the same moment and checks
what MATLAB's duplicate-writer rules promise: exactly one result, complete whenever it can be
seen, never removed or replaced once it was seen, exactly one writer that published it, no
temporary results or flag folders left. ``mode='publish'``: every writer publishes;
``mode='writers'``: every writer registers (``jobs.writers``), waits for the earlier ones and
only the first computes. ``semantics='windows'`` emulates Windows' rename (no replace) on
Linux and macOS. For writers on several machines on one network file system, run
``python -m petakit5d.testing.publish_spy node`` on each (see ``main``).
"""

from __future__ import annotations

import argparse
import errno
import hashlib
import json
import os
import random
import re
import socket
import sys
import threading
import time
import uuid as _uuid
from collections.abc import Callable, Sequence
from dataclasses import asdict, dataclass, field
from typing import Any

__all__ = [
    "PublishSpy",
    "RaceResult",
    "check_round",
    "emulate_windows",
    "list_files",
    "normalize",
    "race_round",
    "verify_result",
    "windows_rename",
]


# ------------------------------------------------------------------------------------------
# the spy
# ------------------------------------------------------------------------------------------


def normalize(p: str | os.PathLike[str]) -> str:
    """forward slashes, no repeated (except a leading //) or trailing separators: the code
    under test builds paths like 'dir//name'"""
    p = re.sub(r"(?<=.)/+", "/", os.fspath(p).replace("\\", "/"))
    return p.rstrip("/") or p[:1]


def list_files(folder: str | os.PathLike[str]) -> list[str]:
    """files and folders under folder, as sorted paths relative to it with '/', without the
    contents of zarr folders (entries whose folder path has '.zarr'); [] if it isn't a
    folder"""
    folder = os.fspath(folder)
    out: list[str] = []
    if not os.path.isdir(folder):
        return out
    for root, dirs, fns in os.walk(folder):
        rel = os.path.relpath(root, folder).replace("\\", "/")
        rel = "" if rel == "." else rel
        if ".zarr" in rel:
            dirs[:] = []
            continue
        for name in dirs + fns:
            out.append(f"{rel}/{name}" if rel else name)
    return sorted(out)


@dataclass
class SpyCall:
    op: str  # 'rename', 'replace' or 'link'
    src: str
    dst: str
    files: list[str]  # list_files of dst's folder before the call
    ok: bool | None = None  # whether the call succeeded (None while it runs)


class PublishSpy:
    """records os.rename, os.replace and os.link calls of this process (see the module
    docstring); use it as a context manager, or install() and uninstall()"""

    def __init__(self) -> None:
        self.calls: list[SpyCall] = []
        self._lock = threading.Lock()
        self._saved: dict[str, Callable[..., Any]] = {}

    def install(self) -> PublishSpy:
        if self._saved:
            return self
        for op in ("rename", "replace", "link"):
            real = getattr(os, op)
            self._saved[op] = real
            setattr(os, op, self._wrap(op, real))
        return self

    def uninstall(self) -> None:
        for op, real in self._saved.items():
            setattr(os, op, real)
        self._saved.clear()

    def __enter__(self) -> PublishSpy:
        return self.install()

    def __exit__(self, *exc: object) -> None:
        self.uninstall()

    def _wrap(self, op: str, real: Callable[..., Any]) -> Callable[..., Any]:
        def spy(src, dst, *args, **kwargs):
            dst_n = normalize(dst)
            folder = dst_n.rsplit("/", 1)[0] if "/" in dst_n else "."
            call = SpyCall(op, normalize(src), dst_n, list_files(folder))
            with self._lock:
                self.calls.append(call)
            try:
                out = real(src, dst, *args, **kwargs)
            except BaseException:
                call.ok = False
                raise
            call.ok = True
            return out

        spy.__wrapped__ = real  # type: ignore[attr-defined]
        return spy

    def moves_to(self, dst: str | os.PathLike[str]) -> list[SpyCall]:
        """the recorded calls whose destination is dst"""
        d = normalize(dst)
        return [c for c in self.calls if c.dst == d]

    def files_when_moved(self, dst: str | os.PathLike[str]) -> list[str]:
        """the files next to dst (relative to its folder) when something was last moved to
        dst; [] if nothing was"""
        calls = self.moves_to(dst)
        return list(calls[-1].files) if calls else []

    def was_moved(self, dst: str | os.PathLike[str]) -> bool:
        """whether something was moved (or linked) to dst"""
        return bool(self.moves_to(dst))


# ------------------------------------------------------------------------------------------
# the race
# ------------------------------------------------------------------------------------------

CHUNKS = 8  # files of a folder result
CHUNK_BYTES = 16384
MANIFEST = "manifest.txt"  # written last into a folder result: the chunks' hashes


@dataclass
class RaceResult:
    """what one round of the race found (``ok``: all checks held)"""

    result: str
    writers: int
    kind: str
    mode: str
    published: list[str] = field(default_factory=list)  # writers that published
    computed: list[str] = field(default_factory=list)  # writers that computed (mode writers)
    waited: list[str] = field(default_factory=list)  # writers that found the result made
    content_of: str = ""  # the writer whose data the result holds
    attempt_spread_ms: float = 0.0  # first to last publish attempt
    # (mode writers) writers that saw fewer records than writers before their wait
    late_records: list[str] = field(default_factory=list)
    nodes: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors


def _payload(writer: str, i: int) -> bytes:
    seed = hashlib.sha256(f"{writer}/{i}".encode()).digest()
    head = f"writer={writer}\nchunk={i}\n".encode()
    return head + (seed * (CHUNK_BYTES // len(seed) + 1))[: CHUNK_BYTES - len(head)]


def _make_temp(tmp: str, kind: str, writer: str) -> None:
    # the writer's temporary result: a folder of chunks and a manifest written last, or a
    # file with its hash at the end
    if kind == "folder":
        os.mkdir(tmp)
        lines = []
        for i in range(CHUNKS):
            data = _payload(writer, i)
            with open(os.path.join(tmp, f"c{i}"), "wb") as f:
                f.write(data)
            lines.append(f"c{i} {hashlib.sha256(data).hexdigest()}")
        with open(os.path.join(tmp, MANIFEST), "w", encoding="utf-8") as f:
            f.write(f"writer={writer}\n" + "\n".join(lines) + "\n")
    else:
        data = b"".join(_payload(writer, i) for i in range(CHUNKS))
        with open(tmp, "wb") as f:
            f.write(data + f"\nsha256={hashlib.sha256(data).hexdigest()}\n".encode())


def verify_result(path: str, kind: str) -> str:
    """the writer whose complete result path holds; raises ValueError if it is incomplete or
    mixed"""
    if kind == "folder":
        with open(os.path.join(path, MANIFEST), encoding="utf-8") as f:
            lines = f.read().splitlines()
        writer = lines[0].removeprefix("writer=")
        names = sorted(os.listdir(path))
        want = sorted([MANIFEST] + [f"c{i}" for i in range(CHUNKS)])
        if names != want:
            raise ValueError(f"{path} holds {names}")
        for line in lines[1:]:
            name, digest = line.split()
            with open(os.path.join(path, name), "rb") as f:
                data = f.read()
            if hashlib.sha256(data).hexdigest() != digest or data != _payload(
                writer, int(name[1:])
            ):
                raise ValueError(f"{path}/{name} is not {writer}'s")
        return writer
    with open(path, "rb") as f:
        data = f.read()
    body, _, tail = data.rpartition(b"\nsha256=")
    if hashlib.sha256(body).hexdigest().encode() + b"\n" != tail:
        raise ValueError(f"{path} is incomplete")
    writer = body.split(b"\n", 1)[0].decode().removeprefix("writer=")
    if body != b"".join(_payload(writer, i) for i in range(CHUNKS)):
        raise ValueError(f"{path} is not {writer}'s")
    return writer


_real_rename, _real_link, _real_unlink = os.rename, os.link, os.unlink


def windows_rename(src: str | os.PathLike[str], dst: str | os.PathLike[str]) -> None:
    """os.rename as Windows does it (MoveFileEx without replace), on Linux and macOS: it fails
    with FileExistsError if dst exists. A file is linked and unlinked (atomic, no replace);
    a folder renamed (that fails onto a folder with files in it)."""
    if os.path.isdir(src):
        if os.path.lexists(dst):
            raise FileExistsError(errno.EEXIST, "Cannot create a file when it exists", dst)
        try:
            _real_rename(src, dst)
        except OSError as err:
            if err.errno in (errno.EEXIST, errno.ENOTEMPTY):
                raise FileExistsError(
                    errno.EEXIST, "Cannot create a file when it exists", dst
                ) from None
            raise
        return
    _real_link(src, dst)
    _real_unlink(src)


def emulate_windows(monkeypatch: Any = None) -> None:
    """publish as on Windows from now on (with a pytest monkeypatch: for one test):
    os.rename is ``windows_rename`` and jobs.publish takes the Windows path"""
    from petakit5d.jobs import publish

    if monkeypatch is None:
        os.rename = windows_rename
        publish._WINDOWS = True
    else:
        monkeypatch.setattr(os, "rename", windows_rename)
        monkeypatch.setattr(publish, "_WINDOWS", True)


def _wait_until(t: float) -> None:
    while True:
        d = t - time.time()
        if d <= 0:
            return
        time.sleep(min(d, 0.05) if d > 0.002 else 0)


def _count_records(result: str) -> int:
    from petakit5d.jobs import files, writers

    folder, name, _ = files.fileparts(result)
    root = os.path.join(folder, writers.FLAG_FOLDER)
    if os.path.isdir(root):
        writers.refresh_listing(root)
    try:
        subs = [e.name for e in os.scandir(root) if e.name.startswith(f"{name}_")]
    except OSError:
        return 0
    return sum(os.path.isfile(os.path.join(root, s, writers.RECORD_NAME)) for s in subs)


def _writer_main(spec: dict[str, Any]) -> None:
    # one writer of the race (a process of its own)
    from petakit5d.jobs import files, publish, writers

    t_start = time.time()
    if spec["semantics"] == "windows":
        emulate_windows()
    if not spec.get("refresh", True):
        # (to measure what refresh_listing does)
        writers.refresh_listing = lambda folder: None
    result, kind, writer = spec["result"], spec["kind"], spec["writer"]
    uid = f"{writer}_{_uuid.uuid4().hex[:8]}"
    tmp = publish.temp_path(result, uid)
    folder, name, _ = files.fileparts(result)
    flag = os.path.join(folder, writers.FLAG_FOLDER, f"{name}_{uid}")
    os.makedirs(flag, exist_ok=True)
    out: dict[str, Any] = {"writer": writer, "node": socket.gethostname(), "pid": os.getpid()}
    out["t_start"] = t_start
    try:
        computed = True
        if spec["mode"] == "writers":
            writers.register(result, uid, flag)
            out["t_registered"] = time.time()
            deadline = time.time() + spec["barrier_timeout"]
            while _count_records(result) < spec["total"] and time.time() < deadline:
                time.sleep(0.02)
            out["records_seen"] = _count_records(result)
            out["t_barrier"] = time.time()
            out["registered_all"] = out["records_seen"] >= spec["total"]
            _wait_until(spec["go"])
            if writers.wait(result, uid):
                computed = False
                out["waited"] = True
                publish.remove_path(flag, folder_only=True)
        if computed:
            _make_temp(tmp, kind, writer)
            with open(os.path.join(flag, "blocks_1_1.flag"), "w"):
                pass
            if spec["mode"] == "writers":
                time.sleep(spec["compute_time"])
            else:
                _wait_until(spec["go"] + random.uniform(0, spec["jitter"]))
            out["t_attempt"] = time.time()
            out["published"] = publish.publish_result(tmp, result, flag)
            out["t_done"] = time.time()
        out["computed"] = computed
    except BaseException as err:
        out["error"] = f"{type(err).__name__}: {err}"
    tmp_out = f"{spec['outcome']}.{os.getpid()}.tmp"
    with open(tmp_out, "w", encoding="utf-8") as f:
        json.dump(out, f)
    os.replace(tmp_out, spec["outcome"])


class _Watcher(threading.Thread):
    # looks at the result every 2 ms: once seen it must stay, be the same file or
    # folder and be complete
    def __init__(self, result: str, kind: str, verify_every: float = 0.05) -> None:
        super().__init__(daemon=True)
        self.result, self.kind, self.verify_every = result, kind, verify_every
        self.stop = threading.Event()
        self.errors: list[str] = []
        self.identity: tuple[int, int] | None = None
        self.writer = ""
        self.looks = 0

    def run(self) -> None:
        parent = os.path.dirname(self.result)
        last_verify = 0.0
        denied = 0
        while not self.stop.is_set():
            self.looks += 1
            try:
                with os.scandir(parent):
                    pass
                st = os.lstat(self.result)
            except FileNotFoundError:
                if self.identity is not None:
                    self.errors.append(f"the result was gone after it was seen ({time.time():.3f})")
                    return
                time.sleep(0.002)
                continue
            ident = (st.st_dev, st.st_ino)
            if self.identity is None:
                self.identity = ident
            elif ident != self.identity:
                self.errors.append("the result was replaced after it was seen")
                return
            now = time.monotonic()
            if now - last_verify >= self.verify_every:
                last_verify = now
                try:
                    w = verify_result(self.result, self.kind)
                except PermissionError as err:
                    # Windows refuses a file that a scanner (or a closing handle) holds for a
                    # moment after its rename: looked at again, up to 2 s; elsewhere an error
                    if sys.platform.startswith("win") and denied < 100:
                        denied += 1
                        last_verify = 0.0
                        time.sleep(0.02)
                        continue
                    self.errors.append(f"an incomplete result was seen: {err}")
                    return
                except (OSError, ValueError) as err:
                    self.errors.append(f"an incomplete result was seen: {err}")
                    return
                denied = 0
                if self.writer and w != self.writer:
                    self.errors.append(f"the result changed from {self.writer}'s to {w}'s")
                    return
                self.writer = w
            time.sleep(0.002)


def _start_writers(specs: Sequence[dict[str, Any]]) -> list[Any]:
    import multiprocessing as mp

    ctx = mp.get_context("spawn")
    # (daemons: writers left by an interrupted round end with this process)
    procs = [ctx.Process(target=_writer_main, args=(s,), daemon=True) for s in specs]
    for p in procs:
        p.start()
    return procs


def run_writers(
    result: str,
    names: Sequence[str],
    outcome_dir: str,
    *,
    total: int,
    kind: str = "folder",
    mode: str = "publish",
    semantics: str = "native",
    go: float,
    jitter: float = 0.0,
    compute_time: float = 1.5,
    barrier_timeout: float = 60.0,
    timeout: float = 300.0,
    refresh: bool = True,
) -> tuple[list[str], str]:
    """Starts the writers ``names`` of one round in processes of their own and watches the
    result until they end; returns the watcher's errors and the writer whose data it saw.
    Their outcomes go to ``<outcome_dir>/<writer>.json``."""
    os.makedirs(outcome_dir, exist_ok=True)
    specs = [
        {
            "result": result,
            "writer": w,
            "kind": kind,
            "mode": mode,
            "semantics": semantics,
            "go": go,
            "jitter": jitter,
            "total": total,
            "compute_time": compute_time,
            "barrier_timeout": barrier_timeout,
            "refresh": refresh,
            "outcome": os.path.join(outcome_dir, f"{w}.json"),
        }
        for w in names
    ]
    procs = _start_writers(specs)
    watcher = _Watcher(result, kind)
    watcher.start()
    errors = []
    try:
        deadline = time.time() + timeout
        for p in procs:
            p.join(max(0.0, deadline - time.time()))
    finally:
        for p in procs:
            if p.is_alive():
                # (a test may end processes it started)
                p.kill()
                p.join()
                errors.append(f"writer process {p.pid} did not end in {timeout} s")
        time.sleep(0.05)
        watcher.stop.set()
        watcher.join()
    return errors + watcher.errors, watcher.writer


def check_round(
    result: str,
    outcome_dir: str,
    names: Sequence[str],
    kind: str,
    mode: str,
    errors: Sequence[str] = (),
) -> RaceResult:
    """Checks one round after its writers ended: their outcomes, the result and its
    folder (no temporary results or flag folders of the round left)."""
    from petakit5d.jobs import files, writers

    rr = RaceResult(result, len(names), kind, mode, errors=list(errors))
    attempts = []
    nodes = set()
    for w in names:
        fn = os.path.join(outcome_dir, f"{w}.json")
        try:
            with open(fn, encoding="utf-8") as f:
                out = json.load(f)
        except (OSError, ValueError) as err:
            rr.errors.append(f"no outcome of {w}: {err}")
            continue
        nodes.add(out.get("node", ""))
        if out.get("error"):
            rr.errors.append(f"{w}: {out['error']}")
        if out.get("published"):
            rr.published.append(w)
        if out.get("computed"):
            rr.computed.append(w)
        if out.get("waited"):
            rr.waited.append(w)
        if "t_attempt" in out:
            attempts.append(out["t_attempt"])
        if mode == "writers" and not out.get("registered_all", True):
            # (the race's own barrier timed out: reported, the checks below decide)
            rr.late_records.append(w)
    rr.nodes = sorted(nodes)
    if attempts:
        rr.attempt_spread_ms = 1000 * (max(attempts) - min(attempts))
    if len(rr.published) != 1:
        rr.errors.append(f"{len(rr.published)} writers published: {rr.published}")
    try:
        rr.content_of = verify_result(result, kind)
    except (OSError, ValueError) as err:
        rr.errors.append(f"the result is not complete: {err}")
    if rr.published and rr.content_of and rr.content_of != rr.published[0]:
        rr.errors.append(f"the result holds {rr.content_of}'s data, {rr.published[0]} published")
    if mode == "writers":
        if len(rr.computed) != 1:
            rr.errors.append(f"{len(rr.computed)} writers computed: {rr.computed}")
        if len(rr.waited) != len(names) - len(rr.computed):
            rr.errors.append(f"{len(rr.waited)} writers found the result made, not all others")
    folder, name, ext = files.fileparts(result)
    left = [x for x in os.listdir(folder) if x.startswith(f"{name}{ext}_")]
    if left:
        rr.errors.append(f"temporary results left: {left}")
    flag_root = os.path.join(folder, writers.FLAG_FOLDER)
    if os.path.isdir(flag_root):
        flags = [x for x in os.listdir(flag_root) if x.startswith(f"{name}_")]
        if flags:
            rr.errors.append(f"flag folders left: {flags}")
    return rr


def race_round(
    folder: str,
    n: int,
    kind: str = "folder",
    mode: str = "publish",
    semantics: str = "native",
    jitter: float = 0.0,
    prepare: float = 3.0,
    name: str = "result",
) -> RaceResult:
    """One round on one machine: ``n`` writers of ``<folder>/<name>.zarr`` (folder) or
    ``<name>.tif`` (file) start publishing at the same moment (``prepare`` s from now, each
    up to ``jitter`` s later); see the module docstring."""
    result = os.path.join(folder, f"{name}.zarr" if kind == "folder" else f"{name}.tif")
    outcome_dir = os.path.join(folder, f"{name}_outcomes")
    names = [f"w{i:02d}" for i in range(n)]
    go = time.time() + prepare
    errors, _ = run_writers(
        result,
        names,
        outcome_dir,
        total=n,
        kind=kind,
        mode=mode,
        semantics=semantics,
        go=go,
        jitter=jitter,
    )
    return check_round(result, outcome_dir, names, kind, mode, errors)


# ------------------------------------------------------------------------------------------
# several machines: python -m petakit5d.testing.publish_spy node ...
# ------------------------------------------------------------------------------------------

ROUND_KINDS = (
    # (writers in all, kind, mode, jitter in s)
    *(
        (n, k, "publish", j)
        for n in (2, 4, 8, 16, 32)
        for k in ("folder", "file")
        for j in (0.0, 0.005)
    ),
    *((n, "folder", "writers", 0.0) for n in (2, 8, 32)),
)


def _rounds(repeats: int, modes: Sequence[str] = ("publish", "writers")) -> list[tuple]:
    return [r for _ in range(repeats) for r in ROUND_KINDS if r[2] in modes]


def _split(n: int, tags: Sequence[str], tag: str) -> list[str]:
    # the writers of tag when n writers are dealt to the tags in turn
    return [f"{tags[i % len(tags)]}{i:02d}" for i in range(n) if tags[i % len(tags)] == tag]


def _wait_for(path: str, timeout: float) -> bool:
    deadline = time.time() + timeout
    parent = os.path.dirname(path)
    while time.time() < deadline:
        try:
            with os.scandir(parent):
                pass
        except OSError:
            pass
        if os.path.exists(path):
            return True
        time.sleep(0.2)
    return False


def node_main(args: argparse.Namespace) -> int:
    """One machine's part of a race on a shared folder: meets the other machines (files
    ``ready_<tag>``), runs its writers of every round at the agreed times, and the first tag
    checks all rounds at the end (``summary.json``; exit status 1 if a check failed)."""
    tags = args.tags.split(",")
    base = os.path.abspath(args.dir)
    os.makedirs(base, exist_ok=True)
    host = socket.gethostname()
    with open(os.path.join(base, f"ready_{args.tag}"), "w", encoding="utf-8") as f:
        f.write(host)
    for t in tags:
        if not _wait_for(os.path.join(base, f"ready_{t}"), args.meet_timeout):
            print(f"no peer {t} within {args.meet_timeout} s", file=sys.stderr)
            return 2
    hosts = {}
    for t in tags:
        with open(os.path.join(base, f"ready_{t}"), encoding="utf-8") as f:
            hosts[t] = f.read().strip()
    print(f"machines: {hosts}", flush=True)
    if len(tags) > 1 and len(set(hosts.values())) < len(tags):
        print("two tags run on the same machine", file=sys.stderr)
        return 2
    go_file = os.path.join(base, "go.json")
    if args.tag == tags[0]:
        with open(go_file + ".tmp", "w", encoding="utf-8") as f:
            json.dump({"t0": time.time() + 15}, f)
        os.replace(go_file + ".tmp", go_file)
    if not _wait_for(go_file, 120):
        return 2
    with open(go_file, encoding="utf-8") as f:
        t0 = json.load(f)["t0"]
    rounds = _rounds(args.repeats, args.modes.split(","))
    for r, (n, kind, mode, jitter) in enumerate(rounds):
        folder = os.path.join(base, f"round_{r:03d}")
        os.makedirs(folder, exist_ok=True)
        result = os.path.join(folder, "result.zarr" if kind == "folder" else "result.tif")
        mine = _split(n, tags, args.tag)
        start = t0 + r * args.period
        _wait_until(start)
        errors, _ = run_writers(
            result,
            mine,
            os.path.join(folder, "outcomes"),
            total=n,
            kind=kind,
            mode=mode,
            semantics=args.semantics,
            go=start + args.prepare,
            jitter=jitter,
            barrier_timeout=args.period / 2,
            refresh=not args.no_refresh,
            timeout=args.period * 3,
        )
        with open(os.path.join(folder, f"watch_{args.tag}.json"), "w", encoding="utf-8") as f:
            json.dump(errors, f)
        print(
            f"round {r} ({n} {kind} {mode} jitter {jitter}): {len(mine)} writers here, "
            f"watcher {'ok' if not errors else errors}",
            flush=True,
        )
    with open(os.path.join(base, f"done_{args.tag}"), "w", encoding="utf-8") as f:
        f.write(host)
    if args.tag != tags[0]:
        return 0
    for t in tags:
        _wait_for(os.path.join(base, f"done_{t}"), args.period * 3 + 600)
    time.sleep(2)
    summary = []
    bad = 0
    for r, (n, kind, mode, _) in enumerate(rounds):
        folder = os.path.join(base, f"round_{r:03d}")
        result = os.path.join(folder, "result.zarr" if kind == "folder" else "result.tif")
        errors = []
        for t in tags:
            try:
                with open(os.path.join(folder, f"watch_{t}.json"), encoding="utf-8") as f:
                    errors += [f"{t}: {e}" for e in json.load(f)]
            except OSError as err:
                errors.append(f"no watcher result of {t}: {err}")
        names = [x for t in tags for x in _split(n, tags, t)]
        rr = check_round(result, os.path.join(folder, "outcomes"), names, kind, mode, errors)
        bad += not rr.ok
        summary.append(asdict(rr) | {"round": r, "ok": rr.ok})
    with open(os.path.join(base, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=1)
    print(f"{len(rounds) - bad} of {len(rounds)} rounds ok", flush=True)
    for s in summary:
        if not s["ok"]:
            print(f"round {s['round']}: {s['errors']}", flush=True)
    return 1 if bad else 0


def main(argv: Sequence[str] | None = None) -> int:
    """``python -m petakit5d.testing.publish_spy node --dir D --tag A --tags A,B`` on each
    machine (``--tags A`` alone: one machine, e.g. a local disk)"""
    p = argparse.ArgumentParser(prog="python -m petakit5d.testing.publish_spy")
    sub = p.add_subparsers(dest="cmd", required=True)
    n = sub.add_parser("node", help="one machine's part of a race on a shared folder")
    n.add_argument("--dir", required=True, help="the shared folder (made if needed)")
    n.add_argument("--tag", required=True)
    n.add_argument("--tags", required=True, help="all machines' tags, the checker first")
    n.add_argument("--repeats", type=int, default=3, help="times each kind of round runs")
    n.add_argument("--period", type=float, default=12.0, help="s between round starts")
    n.add_argument("--prepare", type=float, default=5.0, help="s from a round's start to go")
    n.add_argument("--semantics", default="native", choices=("native", "windows"))
    n.add_argument("--modes", default="publish,writers", help="kinds of rounds to run")
    n.add_argument("--meet-timeout", type=float, default=3600.0)
    n.add_argument("--no-refresh", action="store_true", help="writers list without refresh")
    args = p.parse_args(argv)
    return node_main(args)


if __name__ == "__main__":
    sys.exit(main())
