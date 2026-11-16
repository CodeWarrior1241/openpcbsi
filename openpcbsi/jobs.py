"""Job slots for the heavy solvers (#33): KiPI and openEMS runs from any number of openpcbsi processes
share the machine through a fixed number of slots (lock files in the cache) and a thread budget.

    with jobs.slot("kipi") as threads:      # blocks until a slot is free; threads = this job's share
        ...

Slots: OPENPCBSI_JOB_SLOTS (default: CPU count // 16, at least 1). Threads per job:
OPENPCBSI_JOB_THREADS (default: CPU count // slots).
"""
import contextlib
import fcntl
import os
import time

from .run import CACHE


def _n():
    cpu = os.cpu_count() or 4
    slots = int(os.environ.get("OPENPCBSI_JOB_SLOTS", max(1, cpu // 16)))
    threads = int(os.environ.get("OPENPCBSI_JOB_THREADS", max(1, cpu // slots)))
    return slots, threads


@contextlib.contextmanager
def slot(kind, log=None):
    slots, threads = _n()
    os.makedirs(CACHE, exist_ok=True)
    open(os.path.join(CACHE, "jobs.lock"), "a").close()      # marker: the cache is shared through job slots
    t0, fh = time.time(), None
    while fh is None:
        for k in range(slots):
            f = open(os.path.join(CACHE, f"jobs.lock.{k}"), "a+")
            try:
                fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
                f.seek(0)
                f.truncate()
                f.write(f"{os.getpid()} {kind} {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.flush()
                fh = f
                break
            except OSError:
                f.close()
        if fh is None:
            if log and int(time.time() - t0) % 60 == 0:
                log(f"  waiting for a free {kind} job slot ({slots} slots in use)")
            time.sleep(2)
    try:
        yield threads
    finally:
        fcntl.flock(fh, fcntl.LOCK_UN)
        fh.close()
