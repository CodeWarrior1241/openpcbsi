"""G8: region solver run manager (openEMS).

  python3 -m openpcbsi region BOARD -c CONFIG REGION.toml

REGION.toml:
  [region]
  name = "em1_fcu_line"
  box = [x0, y0, x1, y1]           # mm, board coordinates
  bottom = "In1.Cu"                # the stack is cut at this reference plane (PEC below)
  ports = [ {name="p1", at=[x, y], layer="F.Cu"}, {name="p2", pad="U1.AB19"} ]
  fmax_ghz = 20
  res_mm = [0.05, 0.035]           # two grid resolutions; the result counts only if they agree
  tol_db = 0.5
  compare_z = {layer="F.Cu", w=0.12}   # optional: compare the extracted line impedance with G6

The job (geometry, materials, ports, solver settings) is written as JSON and hashed; the solve runs
in the openEMS venv (em_runner.py) and is cached under runs/em/<hash>/. Copper is modelled as PEC
sheets at each layer's height, vias as PEC cylinders, dielectrics as lossy slabs at f = fmax/2."""
import hashlib
import json
import math
import os
import subprocess
import sys

import numpy as np

VENV_PY = "/media/fpgadev/Dev_Tools/OpenEMS/install/venv/bin/python"


def _bbox_hit(pts, box, m=0.2):
    a = np.asarray(pts)
    return not (a[:, 0].max() < box[0] - m or a[:, 0].min() > box[2] + m or a[:, 1].max() < box[1] - m or a[:, 1].min() > box[3] + m)


def _clip(x1, y1, x2, y2, rect):
    """Liang-Barsky: the part of segment (x1,y1)-(x2,y2) inside rect = [x0, y0, x1, y1]."""
    dx, dy = x2 - x1, y2 - y1
    t0, t1 = 0.0, 1.0
    for p, q in ((-dx, x1 - rect[0]), (dx, rect[2] - x1), (-dy, y1 - rect[1]), (dy, rect[3] - y1)):
        if p == 0:
            if q < 0:
                return None
            continue
        r = q / p
        if p < 0:
            t0 = max(t0, r)
        else:
            t1 = min(t1, r)
        if t0 > t1:
            return None
    return x1 + t0 * dx, y1 + t0 * dy, x1 + t1 * dx, y1 + t1 * dy


def build_job(ctx, spec):
    b, st = ctx.b, ctx.b.stack
    box = spec["box"]
    bottom = spec.get("bottom")
    top_li = 0
    bot_li = b.L[bottom] if bottom else len(b.layers) - 1
    layers = b.layers[top_li:bot_li + 1]
    zmap = {l: -(st.cu(l)["z_top"] + st.cu(l)["z_bot"]) / 2 for l in layers}   # z up, 0 at the board top
    copper = []
    T = b.tracks
    # nets under a port end at the ports: their tracks are clipped to the ports' bounding box,
    # so no line runs on past a port into the absorbing boundary
    port_pts = []
    for p in spec["ports"]:
        if "at" in p:
            port_pts.append((p["at"][0], p["at"][1], b.L[p["layer"]]))
    port_nets = set()
    for (px, py, li) in port_pts:
        r = ctx.g.sample(b.raster(li), [px], [py])[0]
        if r:
            port_nets.add(int(r))
    if port_pts:
        xs_, ys_ = [q[0] for q in port_pts], [q[1] for q in port_pts]
        lng = (max(xs_) - min(xs_)) >= (max(ys_) - min(ys_))
        prect = [min(xs_), box[1], max(xs_), box[3]] if lng else [box[0], min(ys_), box[2], max(ys_)]
    for k in range(len(T["x1"])):
        ln = b.layers[T["layer"][k]]
        if ln not in zmap:
            continue
        x1, y1, x2, y2, w = T["x1"][k], T["y1"][k], T["x2"][k], T["y2"][k], T["w"][k]
        if not _bbox_hit([[x1, y1], [x2, y2]], box, w):
            continue
        if int(T["net"][k]) in port_nets and spec.get("clip_at_ports", True):
            c = _clip(x1, y1, x2, y2, prect)
            if c is None:
                continue
            x1, y1, x2, y2 = c
        L = math.hypot(x2 - x1, y2 - y1) or 1e-9
        nx_, ny_ = -(y2 - y1) / L * w / 2, (x2 - x1) / L * w / 2
        poly = [[x1 + nx_, y1 + ny_], [x2 + nx_, y2 + ny_], [x2 - nx_, y2 - ny_], [x1 - nx_, y1 - ny_]]
        copper.append({"z": zmap[ln], "poly": poly, "net": b.net_names[T["net"][k]]})
    for p in b.pads:
        for li in p["layer_ids"]:
            ln = b.layers[li]
            if ln in zmap and p["poly"]:
                for poly in p["poly"]:
                    if _bbox_hit(poly[0], box):
                        copper.append({"z": zmap[ln], "poly": poly[0], "net": p["net"]})
    for z in b.zones:
        if z["layer"] not in zmap:
            continue
        for poly in z["polys"]:
            if _bbox_hit(poly[0], box):
                # clip far-away vertices: keep the polygon (openEMS discretises only inside the mesh)
                copper.append({"z": zmap[z["layer"]], "poly": poly[0], "net": z["net"], "zone": True})
    vias = []
    V = b.vias
    for k in range(len(V["x"])):
        if box[0] - 0.3 <= V["x"][k] <= box[2] + 0.3 and box[1] - 0.3 <= V["y"][k] <= box[3] + 0.3:
            zt = zmap.get(b.layers[V["top"][k]], 0.0)
            zb = zmap.get(b.layers[min(V["bot"][k], bot_li)], zmap[layers[-1]])
            vias.append({"x": V["x"][k], "y": V["y"][k], "r": V["drill"][k] / 2 + 0.02, "z0": zb, "z1": zt})
    diel = []
    for d in st.dielectrics:
        if d["above"] in zmap and d.get("below") in zmap:
            diel.append({"z0": -d["z_bot"], "z1": -d["z_top"], "er": d["er"], "tand": d["tand"]})
    ports = []
    for p in spec["ports"]:
        if "pad" in p:
            ref, num = p["pad"].split(".")
            q = next(x for x in b.pads_by_ref[ref] if x["num"] == num)
            x, y, ln = q["x"], q["y"], b.layers[q["layer_ids"][0]]
        else:
            (x, y), ln = p["at"], p["layer"]
        refs = [b.layers[r] for r in ctx.ref_layers(b.L[ln])]
        rz = zmap[refs[0]] if refs and refs[0] in zmap else zmap[layers[-1]]
        ports.append({"name": p["name"], "x": x, "y": y, "z_sig": zmap[ln], "z_ref": rz, "R": p.get("R", 50.0)})
    mask = st.mask.get("F")
    job = {"box": box, "z_top_air": 1.0, "z_bottom": zmap[layers[-1]], "copper": copper, "vias": vias, "dielectrics": diel,
           "mask": {"t": mask["t"], "er": mask["er"]} if mask else None, "ports": ports,
           "fmax_ghz": spec.get("fmax_ghz", 20.0), "res_mm": spec.get("res_mm", [0.05, 0.035]),
           "tol_db": spec.get("tol_db", 0.5), "end_criteria": spec.get("end_criteria", 1e-4),
           "max_timesteps": spec.get("max_timesteps", 200000), "threads": spec.get("threads", 4),
           "board_sha": b.sha}
    return job


def run(ctx, spec, cache_dir, runs_dir, log=print, wait=True):
    job = build_job(ctx, spec)
    h = hashlib.sha256(json.dumps(job, sort_keys=True).encode()).hexdigest()[:16]
    od = os.path.join(runs_dir, "em", f"{spec.get('name', 'region')}_{h}")
    os.makedirs(od, exist_ok=True)
    res_path = os.path.join(od, "result.json")
    if os.path.exists(res_path):
        log(f"G8: cached result {res_path}")
        return json.load(open(res_path))
    with open(os.path.join(od, "job.json"), "w") as f:
        json.dump(job, f)
    log(f"G8: {len(job['copper'])} copper polygons, {len(job['vias'])} vias, {len(job['ports'])} ports -> {od}")
    if spec.get("compare_z"):
        from . import xsec
        cz = spec["compare_z"]
        li = ctx.b.L[cz["layer"]]
        refs = [ctx.b.layers[r] for r in ctx.ref_layers(li)]
        sec = xsec.build_section(ctx.b.stack, cz["layer"], refs, [cz["w"]], [], ctx.cfg["xsec"], ctx.cfg["xsec"]["mask_er"])
        r2 = xsec.solve_lines(sec, 5.0, (1e9,))["freqs"][0]
        job_cmp = {"Z0_2d": r2["Z0"], "eeff_2d": r2["eeff"]}
        with open(os.path.join(od, "compare_2d.json"), "w") as f:
            json.dump(job_cmp, f)
    runner = os.path.join(os.path.dirname(os.path.abspath(__file__)), "em_runner.py")
    env = dict(os.environ, OMP_NUM_THREADS=str(job["threads"]))
    cmd = [VENV_PY, "-I", runner, od]
    if not wait:
        return {"started": cmd, "dir": od}
    p = subprocess.run(cmd, env=env, capture_output=True, text=True)
    with open(os.path.join(od, "runner.log"), "w") as f:
        f.write(p.stdout + "\n" + p.stderr)
    if p.returncode != 0:
        log(p.stdout[-2000:] + p.stderr[-3000:])
        raise SystemExit(f"G8: openEMS run failed (log {od}/runner.log)")
    r = json.load(open(res_path))
    log(json.dumps({k: v for k, v in r.items() if k not in ("s",)}, indent=1)[:3000])
    return r
