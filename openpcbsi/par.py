"""Run independent pieces of one check in forked worker processes (the board and its rasters are
shared with the workers, only the items and the results are passed)."""
import multiprocessing as mp
import os

_STATE = {}


def workers(cfg):
    """[dfm] workers, or a quarter of the CPUs up to 8."""
    n = int(cfg["dfm"].get("workers") or 0)
    return n if n > 0 else min(8, max(1, (os.cpu_count() or 2) // 4))


def _call(item):
    return _STATE["fn"](item)


def pmap(fn, items, n):
    """[fn(i) for i in items], in up to n processes. `fn` may be a closure: it is inherited, not pickled."""
    items = list(items)
    if n <= 1 or len(items) <= 1 or "fn" in _STATE:
        return [fn(i) for i in items]
    _STATE["fn"] = fn
    try:
        with mp.get_context("fork").Pool(min(n, len(items))) as pool:
            return pool.map(_call, items, chunksize=1)
    except (OSError, ValueError):
        return [fn(i) for i in items]
    finally:
        _STATE.pop("fn", None)
