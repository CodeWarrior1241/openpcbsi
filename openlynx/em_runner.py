"""G8 runner, executed with the openEMS venv's python:  python -I em_runner.py <job dir>

Builds the CSXCAD model from job.json, runs openEMS at each grid resolution, writes Touchstone
files and result.json with the convergence record. Kept free of openlynx imports (isolated mode)."""
import json
import math
import os
import sys
import time

import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import C0, EPS0


def build(job, res, sim_dir):
    fmax = job["fmax_ghz"] * 1e9
    FDTD = openEMS(NrTS=int(job["max_timesteps"]), EndCriteria=job["end_criteria"])
    FDTD.SetGaussExcite(fmax / 2, fmax / 2)
    FDTD.SetBoundaryCond(["MUR", "MUR", "MUR", "MUR", "PEC", "MUR"])     # PEC at the cut plane below
    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(1e-3)
    x0, y0, x1, y1 = job["box"]
    zb, zt = job["z_bottom"], job["z_top_air"]
    # dielectrics (lossy at fmax/2)
    w = 2 * math.pi * fmax / 2
    zl = [zb, 0.0, zt]
    for k, d in enumerate(job["dielectrics"]):
        m = CSX.AddMaterial(f"d{k}", epsilon=d["er"], kappa=w * EPS0 * d["er"] * d["tand"])
        m.AddBox([x0, y0, d["z0"]], [x1, y1, d["z1"]], priority=1)
        zl += [d["z0"], d["z1"]]
    if job.get("mask"):
        mk = CSX.AddMaterial("mask", epsilon=job["mask"]["er"])
        mk.AddBox([x0, y0, 0.0], [x1, y1, job["mask"]["t"]], priority=2)
        zl += [job["mask"]["t"]]
    cu = CSX.AddMetal("cu")
    xs, ys, zs = [x0, x1], [y0, y1], set()
    for c in job["copper"]:
        a = np.array(c["poly"], dtype=float)
        cu.AddPolygon([a[:, 0].tolist(), a[:, 1].tolist()], "z", c["z"], priority=10)
        zs.add(c["z"])
        if not c.get("zone"):
            inside = (a[:, 0] >= x0) & (a[:, 0] <= x1) & (a[:, 1] >= y0) & (a[:, 1] <= y1)
            xs += a[inside, 0].tolist()
            ys += a[inside, 1].tolist()
    for v in job["vias"]:
        cu.AddCylinder([v["x"], v["y"], v["z0"]], [v["x"], v["y"], v["z1"]], v["r"], priority=10)
    ports = []
    for k, p in enumerate(job["ports"]):
        d = res
        port = FDTD.AddLumpedPort(k + 1, p["R"], [p["x"] - d / 2, p["y"] - d / 2, p["z_ref"]],
                                  [p["x"] + d / 2, p["y"] + d / 2, p["z_sig"]], "z", excite=1.0 if k == 0 else 0, priority=20)
        ports.append(port)
        xs += [p["x"] - d / 2, p["x"] + d / 2]
        ys += [p["y"] - d / 2, p["y"] + d / 2]
    zl += sorted(zs)
    # mesh: key lines snapped and smoothed to the resolution
    def lines(vals, lo, hi):
        v = np.unique(np.round(np.clip(vals, lo, hi), 4))
        keep = [v[0]]
        for q in v[1:]:
            if q - keep[-1] >= res / 3:
                keep.append(q)
        return keep
    mesh.AddLine("x", lines(xs, x0, x1))
    mesh.AddLine("y", lines(ys, y0, y1))
    mesh.AddLine("z", lines(zl, zb, zt))
    mesh.SmoothMeshLines("x", res, 1.3)
    mesh.SmoothMeshLines("y", res, 1.3)
    mesh.SmoothMeshLines("z", min(res, 0.03), 1.3)
    ncell = len(mesh.GetLines("x")) * len(mesh.GetLines("y")) * len(mesh.GetLines("z"))
    return FDTD, CSX, ports, ncell


def main(od):
    job = json.load(open(os.path.join(od, "job.json")))
    fmax = job["fmax_ghz"] * 1e9
    f = np.linspace(fmax / 200, fmax, 200)
    record = {"runs": []}
    S = []
    for res in job["res_mm"]:
        sim = os.path.join(od, f"sim_{res:g}")
        os.makedirs(sim, exist_ok=True)
        FDTD, CSX, ports, ncell = build(job, res, sim)
        t = time.time()
        FDTD.Run(sim, cleanup=True, verbose=0, numThreads=int(job.get("threads", 4)))
        for p in ports:
            p.CalcPort(sim, f)
        a1 = ports[0].uf_inc
        s = np.array([p.uf_ref / a1 if k == 0 else p.uf_ref / a1 for k, p in enumerate(ports)])
        S.append(s)
        record["runs"].append({"res_mm": res, "cells": ncell, "seconds": round(time.time() - t, 1)})
        np.save(os.path.join(od, f"s_{res:g}.npy"), s)
        # Touchstone (first column of the S-matrix, port 1 excited)
        with open(os.path.join(od, f"port1_{res:g}.s{len(ports)}p"), "w") as fh:
            fh.write(f"! openlynx G8, port 1 excited only (column 1), res {res} mm\n# HZ S RI R 50\n")
            for i in range(len(f)):
                vals = []
                for k in range(len(ports)):
                    vals += [f"{s[k][i].real:.6e}", f"{s[k][i].imag:.6e}"]
                fh.write(f"{f[i]:.6e} " + " ".join(vals) + "\n")
    db = [20 * np.log10(np.abs(s) + 1e-12) for s in S]
    diff = float(np.max(np.abs(db[0] - db[-1]))) if len(S) > 1 else None
    record["max_db_difference"] = round(diff, 3) if diff is not None else None
    record["converged"] = bool(diff is not None and diff <= job["tol_db"])
    s = S[-1]
    out = {"freq_hz": f.tolist(), "S11_db": db[-1][0].tolist(), "convergence": record}
    if len(job["ports"]) == 2:
        out["S21_db"] = db[-1][1].tolist()
        # line impedance from the 2-port (assumes S12 = S21, S22 = S11 for a symmetric line)
        s11, s21 = s[0], s[1]
        z0 = 50.0
        num = (1 + s11) ** 2 - s21 ** 2
        den = (1 - s11) ** 2 - s21 ** 2
        zc = z0 * np.sqrt(num / den)
        band = (f > fmax / 20) & (f < fmax / 2)
        out["Z_line_ohm"] = round(float(np.median(np.abs(zc[band]))), 2)
        cmpf = os.path.join(od, "compare_2d.json")
        if os.path.exists(cmpf):
            c = json.load(open(cmpf))
            out["Z0_2d"] = round(c["Z0_2d"], 2)
            out["agreement_pct"] = round(100 * (out["Z_line_ohm"] - c["Z0_2d"]) / c["Z0_2d"], 2)
    with open(os.path.join(od, "result.json"), "w") as fh:
        json.dump(out, fh)
    print(json.dumps({k: v for k, v in out.items() if not isinstance(v, list)}))


if __name__ == "__main__":
    main(sys.argv[1])
