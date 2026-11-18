"""A private, headless KiCad PCB editor for the live-board test (#32).

Starts stock KiCad's pcbnew on a virtual display (xvfb-run), in its own process group, with a private
settings folder (API server on; the user's library tables copied so no first-run dialog blocks it), and
stops only that process group afterwards. The user's own KiCad and settings are not touched.
"""
import glob
import os
import shutil
import signal
import subprocess
import time

PCBNEW = os.environ.get("OPENPCBSI_TEST_PCBNEW", "/usr/bin/pcbnew")
KIPY_PY = os.environ.get("OPENPCBSI_KIPY_PYTHON", "")       # no default: the kicad-python (kipy) venv's Python


KICAD_GUIS = {"kicad", "pcbnew", "eeschema", "pl_editor", "gerbview", "bitmap2component", "pcb_calculator"}


def _other_kicad_running():
    out = subprocess.run(["ps", "-eo", "comm"], capture_output=True, text=True).stdout.split()
    return any(c in KICAD_GUIS for c in out)


def _clear_stale(since=None):
    """Remove /tmp/kicad API lock/socket files, only when no KiCad GUI runs (then nothing owns them);
    with `since`, only the files created after that time (this harness's own)."""
    if _other_kicad_running():
        return
    for f in glob.glob("/tmp/kicad/api*.lock") + glob.glob("/tmp/kicad/api*.sock"):
        if since is None or os.path.getmtime(f) >= since - 1:
            try:
                os.remove(f)
            except OSError:
                pass


class LiveKiCad:
    def __init__(self, board, workdir):
        self.board, self.workdir = os.path.abspath(board), os.path.abspath(workdir)
        self.proc, self.socket = None, None

    def __enter__(self):
        if not KIPY_PY or not os.path.isfile(KIPY_PY):
            raise RuntimeError("the live-board test needs kicad-python (kipy): set OPENPCBSI_KIPY_PYTHON to the Python "
                               f"of an environment that has it installed (now {KIPY_PY or 'unset'})")
        # a copy of the user's KiCad 10 settings (so no first-run wizard blocks the editor), with the API
        # server switched on in the copy only; the originals are only read
        import json
        cfg = os.path.join(self.workdir, "kicad_config", "10.0")
        if os.path.isdir(cfg):
            shutil.rmtree(cfg)
        src = os.path.expanduser("~/.config/kicad/10.0")
        if os.path.isdir(src):
            shutil.copytree(src, cfg)
        os.makedirs(cfg, exist_ok=True)
        cp = os.path.join(cfg, "kicad_common.json")
        common = json.load(open(cp)) if os.path.exists(cp) else {}
        common.setdefault("api", {})["enable_server"] = True
        json.dump(common, open(cp, "w"), indent=2)
        _clear_stale()
        before = {p: os.path.getmtime(p) for p in glob.glob("/tmp/kicad/*.sock")}
        t0 = self.t0 = time.time()
        self.proc = subprocess.Popen(["xvfb-run", "-a", "-s", "-screen 0 1600x1000x24", PCBNEW, self.board],
                                     env=dict(os.environ, KICAD_CONFIG_HOME=os.path.dirname(cfg)),
                                     stdout=open(os.path.join(self.workdir, "pcbnew.log"), "w"), stderr=subprocess.STDOUT,
                                     start_new_session=True)
        deadline = time.time() + 120
        while time.time() < deadline:
            time.sleep(2)
            socks = [p for p in glob.glob("/tmp/kicad/*.sock") if os.path.getmtime(p) >= t0 - 1 or p not in before]
            for s in socks:
                url = "ipc://" + s
                r = subprocess.run([KIPY_PY, "-c", f"from kipy import KiCad; print(KiCad(socket_path={url!r}).get_board().name)"],
                                   capture_output=True, text=True, timeout=30)
                if r.returncode == 0 and r.stdout.strip():
                    self.socket = url
                    return self
        self.__exit__(None, None, None)
        raise RuntimeError("the headless KiCad did not answer on its API within 120 s (see pcbnew.log)")

    def __exit__(self, *exc):
        if self.proc and self.proc.poll() is None:
            os.killpg(self.proc.pid, signal.SIGTERM)       # only the group this harness started
            try:
                self.proc.wait(timeout=20)
            except subprocess.TimeoutExpired:
                os.killpg(self.proc.pid, signal.SIGKILL)
        for _ in range(10):                              # let the group's processes exit, then tidy its files
            if not _other_kicad_running():
                break
            time.sleep(1)
        _clear_stale(since=getattr(self, "t0", None))
        return False
