"""Diagnostic for the CI job exit-probe where the runner has no cdb (thread_stacks.ps1): the threads of a minidump
(x64, as procdump writes it), one line each: its id, suspend count, priority, start address and instruction pointer,
and the first values on its stack that point into a module (return addresses, mostly), as module+offset.
Usage: python dump_threads.py <minidump>"""
import struct
import sys

d = open(sys.argv[1], "rb").read()
sig, _, count, at = struct.unpack_from("<4sIII", d, 0)
if sig != b"MDMP":
    sys.exit("not a minidump")
streams = {}
for k in range(count):
    kind, _, rva = struct.unpack_from("<III", d, at + 12 * k)
    streams[kind] = rva

modules = []
if 4 in streams:  # ModuleListStream
    at = streams[4]
    for k in range(struct.unpack_from("<I", d, at)[0]):
        base, size = struct.unpack_from("<QI", d, at + 4 + 108 * k)
        name = struct.unpack_from("<I", d, at + 4 + 108 * k + 20)[0]
        n = struct.unpack_from("<I", d, name)[0]
        modules.append((base, size, d[name + 4:name + 4 + n].decode("utf-16-le").split("\\")[-1]))


def where(a):
    for base, size, name in modules:
        if base <= a < base + size:
            return f"{name}+0x{a - base:x}"
    return None


starts = {}
if 17 in streams:  # ThreadInfoListStream
    at = streams[17]
    head, entry, n = struct.unpack_from("<III", d, at)
    for k in range(n):
        tid = struct.unpack_from("<I", d, at + head + entry * k)[0]
        starts[tid] = struct.unpack_from("<Q", d, at + head + entry * k + 48)[0]
at = streams[3]  # ThreadListStream
n = struct.unpack_from("<I", d, at)[0]
print(f"{n} threads, {len(modules)} modules")
for k in range(n):
    tid, suspend, _, priority, _, stack, size, rva, csize, crva = struct.unpack_from("<IIIIQQIIII", d, at + 4 + 48 * k)
    rip, rsp = (struct.unpack_from("<Q", d, crva + 0xF8)[0], struct.unpack_from("<Q", d, crva + 0x98)[0]) \
        if csize >= 0x100 else (0, stack)
    found = []
    for o in range(rva + max(rsp - stack, 0) // 8 * 8, rva + size - 7, 8):
        w = where(struct.unpack_from("<Q", d, o)[0])
        if w:
            found.append(w)
            if len(found) == 12:
                break
    s = starts.get(tid, 0)
    print(f"thread {tid}: suspend {suspend}, priority {priority}, start {where(s) or hex(s)}, "
          f"at {where(rip) or hex(rip)}; on the stack {', '.join(found) or '-'}")
