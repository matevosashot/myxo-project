import gc
import logging
import numexpr as ne
import os
import sys
import time
import tracemalloc
from datetime import datetime

import numpy as np
import psutil

tracemalloc.start()

_begin_time = time.perf_counter()
_last_log_time = time.perf_counter()
_process = psutil.Process(os.getpid())

_rss_gb_max = 0
_vms_gb_max = 0

def log(text, notime=False):
    global _last_log_time, _rss_gb_max, _vms_gb_max

    logger = logging.getLogger("main")

    t = time.perf_counter()
    current, peak = tracemalloc.get_traced_memory()
    # mem = _process.memory_info()
    # rss_gb = mem.rss / 1024**3
    # vms_gb = mem.vms / 1024**3
    # _rss_gb_max = max(_rss_gb_max, rss_gb)
    # _vms_gb_max = max(_vms_gb_max, vms_gb)

    dt = t - _last_log_time
    _last_log_time = t
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    if notime:
        logger.debug(f"{text}")
    else:
        logger.debug(f"{dt:4.2f}s | {text:60.60s} | {t - _begin_time:4.2f}s total | {current/1024**3:.2f}GB<{peak/1024**3:.2f}GB | {now}")
        # logger.debug(f"{now}" + " " * 40 + f"RSS {rss_gb:.2f}GB<{_rss_gb_max:.2f}GB | VMS {vms_gb:.2f}GB<{_vms_gb_max:.2f}GB")


def clear_numexpr_last_cache():        
    a = np.array([1, 2, 3])
    b = np.array([4, 5, 6])
    ne.evaluate("a+b", local_dict={"a": a, "b": b}, out=a)


def diag_large_ndarrays(min_gb=0.5, max_paths=3, max_depth=6, also_walk_frames=True):
    """
    Print live numpy arrays larger than ``min_gb`` GiB and the names that
    reach them. Intended for hunting hidden refs that survive a `del`.

    Note: numpy ndarrays are NOT tracked by Python's cyclic GC (they're
    non-container objects), so ``gc.get_objects()`` does not enumerate
    them. We instead walk the object graph from module globals and the
    calling frame stack, discovering ndarrays as we encounter them.

    For each found ndarray (>= ``min_gb`` GiB):
      - shape, dtype, nbytes,
      - up to ``max_paths`` qualified name paths reached during the walk.
    """
    import inspect

    rss_gb = vms_gb = 0.0
    try:
        with open("/proc/self/status") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    rss_gb = int(line.split()[1]) / 1024**2
                elif line.startswith("VmSize:"):
                    vms_gb = int(line.split()[1]) / 1024**2
    except OSError:
        pass
    tm = ""
    if tracemalloc.is_tracing():
        cur, peak = tracemalloc.get_traced_memory()
        tm = f"  tracemalloc cur={cur/1024**3:5.2f}GB peak={peak/1024**3:5.2f}GB"
    print(f"--- mem  RSS={rss_gb:6.2f}GB  VmSize={vms_gb:6.2f}GB{tm}")

    threshold = min_gb * 1024**3
    found = {}              # id(arr) -> (arr, [paths])
    visited = set()         # container ids already walked

    def _record(arr, path):
        key = id(arr)
        if key not in found:
            found[key] = (arr, [])
        if len(found[key][1]) < max_paths:
            found[key][1].append(path)

    def walk(node, prefix, depth):
        if depth > max_depth:
            return
        nid = id(node)
        if nid in visited:
            return
        visited.add(nid)

        if isinstance(node, np.ndarray):
            if node.nbytes >= threshold:
                _record(node, prefix)
            return

        if isinstance(node, dict):
            items = ((repr(k), v) for k, v in node.items())
            joiner = lambda p, k: f"{p}[{k}]"
        elif isinstance(node, (list, tuple)):
            items = ((str(i), v) for i, v in enumerate(node))
            joiner = lambda p, k: f"{p}[{k}]"
        else:
            d = getattr(node, "__dict__", None)
            if d is None or isinstance(node, type):
                return
            items = ((k, v) for k, v in d.items())
            joiner = lambda p, k: f"{p}.{k}"

        for k, v in items:
            try:
                walk(v, joiner(prefix, k), depth + 1)
            except (AttributeError, TypeError, RecursionError):
                pass

    # Walk every loaded module.
    for modname in list(sys.modules):
        mod = sys.modules.get(modname)
        if mod is None or not hasattr(mod, "__dict__"):
            continue
        try:
            walk(mod.__dict__, modname, 0)
        except Exception:
            pass

    # Walk the calling frame stack (locals from the caller, its caller, ...).
    # In Python 3.13+ frame.f_locals is a FrameLocalsProxy, not a dict
    # subclass — convert explicitly so the walk's isinstance(dict) check
    # picks it up.
    if also_walk_frames:
        import threading

        main_tid = threading.get_ident()
        # Current thread (this thread): walk from f_back up.
        frame = inspect.currentframe()
        if frame is not None:
            frame = frame.f_back
        while frame is not None:
            try:
                walk(dict(frame.f_locals), f"<frame {frame.f_code.co_name}>", 0)
            except Exception:
                pass
            frame = frame.f_back

        # All other threads (workers etc.): walk their current top frame
        # and the entire stack via f_back. This catches arrays held by
        # scipy.fft / concurrent.futures worker pools.
        for tid, top_frame in list(sys._current_frames().items()):
            if tid == main_tid:
                continue
            f = top_frame
            while f is not None:
                try:
                    walk(
                        dict(f.f_locals),
                        f"<thread {tid} frame {f.f_code.co_name}>",
                        0,
                    )
                except Exception:
                    pass
                f = f.f_back

    # gc-tracked containers: walk dicts/lists/tuples in gc.get_objects()
    # and inspect function closures. ndarrays themselves aren't gc-tracked,
    # but their containers usually are — this catches arrays held by
    # closures, cell vars, frame cells, or random module-detached dicts.
    # When a hit is found in a gc-tracked dict, also report the dict's other
    # keys and one level of its referrers so the caller can identify which
    # subsystem is holding the cache.
    suspect_dicts = []
    for obj in gc.get_objects():
        try:
            if isinstance(obj, dict):
                for k, v in obj.items():
                    if isinstance(v, np.ndarray) and v.nbytes >= threshold:
                        _record(v, f"<gc dict @{id(obj):x}>[{k!r}]")
                        suspect_dicts.append(obj)
            elif isinstance(obj, (list, tuple)):
                for i, v in enumerate(obj):
                    if isinstance(v, np.ndarray) and v.nbytes >= threshold:
                        _record(v, f"<gc {type(obj).__name__} @{id(obj):x}>[{i}]")
            elif callable(obj) and getattr(obj, "__closure__", None):
                qname = getattr(obj, "__qualname__", repr(obj))
                for cell in obj.__closure__:
                    try:
                        v = cell.cell_contents
                    except ValueError:
                        continue
                    if isinstance(v, np.ndarray) and v.nbytes >= threshold:
                        _record(v, f"<closure of {qname}>")
        except Exception:
            pass

    # For each suspect gc dict, print its keys and what's holding it,
    # so the caller can attribute the cache to a specific subsystem.
    for d in suspect_dicts:
        keys_summary = ", ".join(repr(k) for k in list(d.keys())[:10])
        if len(d) > 10:
            keys_summary += f", ... ({len(d)} keys)"
        print(f"  └ suspect dict @{id(d):x} keys: {{{keys_summary}}}")
        for r in gc.get_referrers(d)[:5]:
            rtype = type(r).__name__
            if isinstance(r, dict):
                # try to find the key that maps to d
                for rk, rv in r.items():
                    if rv is d:
                        print(f"     ↑ held by dict[{rk!r}] of type {rtype}")
                        break
                else:
                    print(f"     ↑ held by {rtype} @{id(r):x}")
            elif callable(r):
                qname = getattr(r, "__qualname__", repr(r))
                mod = getattr(r, "__module__", "?")
                print(f"     ↑ held by callable {mod}.{qname}")
            elif hasattr(r, "__class__"):
                cls = r.__class__
                print(f"     ↑ held by instance of {cls.__module__}.{cls.__qualname__}")
            else:
                print(f"     ↑ held by {rtype}")

    print(f"--- {len(found)} ndarrays >= {min_gb} GiB:")
    for _, (arr, paths) in sorted(found.items(), key=lambda kv: -kv[1][0].nbytes):
        size = arr.nbytes / 1024**3
        print(f"  {size:6.2f} GiB shape={arr.shape!s:30s} dtype={arr.dtype!s:10s}  "
              f"{' | '.join(paths[:max_paths]) or '<no path found>'}")
    print("---")
