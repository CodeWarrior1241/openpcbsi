"""openEMS runner for one grid of a 3D region (#9-#16, #18). Runs in the openEMS python, not in
KiCad's: it reads the model written by region3d.py and writes the port waves.

    <openEMS python> -I em3d.py model.json out_dir threads

Units are um. Copper: tracks and pads as PEC prisms of their thickness (zero-thickness conducting
sheets with the copper's conductivity when the model is lossy); planes as sheets, holes and islands by
nested priority (no polygon library needed); vias as cylinders with pads. Ports: openEMS MSLPorts on
the synthetic lead-ins, fed beyond the PML, measured at the box edge.
"""
import json
import math
import os
import sys

import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS

EPS0 = 8.8541878128e-12
SIGMA_CU = 5.8e7


def main(mp, od, threads):
    od = os.path.abspath(od)                 # openEMS changes directory while it runs
    m = json.load(open(mp))
    fmax = m["fmax"]
    FDTD = openEMS(EndCriteria=m["end_criteria"], NrTS=m["max_timesteps"])
    FDTD.SetGaussExcite(fmax / 2, fmax / 2)
    FDTD.SetBoundaryCond([f"PML_{m['pml_cells']}"] * 4 + ["MUR", "MUR"])
    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    g = CSX.GetGrid()
    g.SetDeltaUnit(1e-6)
    for ax in ("x", "y", "z"):
        g.AddLine(ax, np.array(m["mesh"][ax]))
    x0, y0, x1, y1 = m["domain"]
    f0 = fmax / 2
    zs = sorted(set(m["z_sheets"].values()))
    ztop, zbot = max(zs), min(zs)
    # dielectrics and mask
    for i, d in enumerate(m["dielectrics"]):
        kw = {"epsilon": d["er"]}
        if m["lossy"]:
            kw["kappa"] = 2 * math.pi * f0 * EPS0 * d["er"] * d["tand"]
        CSX.AddMaterial(f"diel{i}", **kw).AddBox([x0, y0, d["z0"]], [x1, y1, d["z1"]], priority=1)
    for i, k in enumerate(m["mask"]):
        CSX.AddMaterial(f"mask{i}", epsilon=k["er"]).AddBox([x0, y0, k["z0"]], [x1, y1, k["z1"]], priority=2)

    def er_at(z):
        below = [d for d in m["dielectrics"] if abs(d["z1"] - z) < 1e-3]
        above = [d for d in m["dielectrics"] if abs(d["z0"] - z) < 1e-3]
        inside = [d for d in m["dielectrics"] if d["z0"] < z < d["z1"]]
        for c in (inside, below, above):
            if c:
                return c[0]["er"]
        return 1.0
    pec = CSX.AddMetal("pec")
    sheet_cu = CSX.AddConductingSheet("cu", conductivity=SIGMA_CU, thickness=35e-6) if m["lossy"] else None
    nh = 0
    for s in m["sheets"]:
        pts = np.array(s["pts"]).T
        z = s["z"]
        if not s["metal"]:
            nh += 1
            CSX.AddMaterial(f"hole{nh}", epsilon=er_at(z)).AddPolygon(pts, "z", z, priority=s["prio"])
            continue
        outer = z >= ztop - 1e-3 or z <= zbot + 1e-3
        if s["kind"] in ("track", "pad") and s["thick"] > 0 and not m["lossy"] and outer:   # inner copper: sheets
            t = s["thick"]
            z0 = z if z >= ztop - 1e-3 else (z - t if z <= zbot + 1e-3 else z - t / 2)
            pec.AddLinPoly(pts, "z", z0, t, priority=s["prio"])
            # conformal solder mask over the copper's top face, as the 2D solver has it
            mk = [k for k in m["mask"] if (k["z0"] >= z - 1e-3) == (z >= ztop - 1e-3)]
            if mk:
                # the copper's top face and sidewalls: the outline grown by the mask thickness, from the
                # copper's base to t_mask over its top (the copper, higher priority, keeps its own volume)
                tm = mk[0]["z1"] - mk[0]["z0"]
                gpts = np.array(s.get("mask_pts") or s["pts"]).T
                zc = z0 if z >= ztop - 1e-3 else z0 - tm
                nh += 1
                CSX.AddMaterial(f"maskcov{nh}", epsilon=mk[0]["er"]).AddLinPoly(gpts, "z", zc, t + tm, priority=3)
        elif sheet_cu is not None:
            sheet_cu.AddPolygon(pts, "z", z, priority=s["prio"])
        else:
            pec.AddPolygon(pts, "z", z, priority=s["prio"])
    for v in m["vias"]:
        pec.AddCylinder([v["x"], v["y"], v["z0"]], [v["x"], v["y"], v["z1"]], v["r"], priority=60)
        ring = np.array([[v["x"] + v["pad_r"] * math.cos(a), v["y"] + v["pad_r"] * math.sin(a)]
                         for a in np.linspace(0, 2 * math.pi, 24, endpoint=False)]).T
        for zp in v["pads"]:
            pec.AddPolygon(ring, "z", zp, priority=60)
    for i, e in enumerate(m["lumped"]):
        if e.get("body"):                        # a part modelled as its metal body (a large series C)
            bd = e["body"]
            pec.AddBox([bd["start"][0], bd["start"][1], bd["z0"]], [bd["stop"][0], bd["stop"][1], bd["z1"]], priority=60)
            continue
        kw = {k: e[k] for k in ("R", "C", "L") if e.get(k)}
        # R alone: openEMS's classic lumped resistor; with C or L: the series R-L-C element
        series = 1 if ("C" in kw or "L" in kw) else 0
        le = CSX.AddLumpedElement(f"le{i}_{e['ref']}", ny=0 if e["dir"] == "x" else 1, caps=True, LEtype=series, **kw)
        le.AddBox([e["start"][0], e["start"][1], e["z"]], [e["stop"][0], e["stop"][1], e["z"]], priority=70)
    # ports
    X, Y = np.array(m["mesh"]["x"]), np.array(m["mesh"]["y"])
    pml = {"xmin": X[m["pml_cells"]] - X[0], "xmax": X[-1] - X[-1 - m["pml_cells"]],
           "ymin": Y[m["pml_cells"]] - Y[0], "ymax": Y[-1] - Y[-1 - m["pml_cells"]]}
    P = []
    for p in sorted(m["ports"], key=lambda q: q["n"]):
        la = p["lead"]
        if p["side"] == "xmin":
            start, stop, d = [la[0], la[1], p["z"]], [la[2], la[3], p["z_ref"]], "x"
        elif p["side"] == "xmax":
            start, stop, d = [la[2], la[1], p["z"]], [la[0], la[3], p["z_ref"]], "x"
        elif p["side"] == "ymin":
            start, stop, d = [la[0], la[1], p["z"]], [la[2], la[3], p["z_ref"]], "y"
        else:
            start, stop, d = [la[0], la[3], p["z"]], [la[2], la[1], p["z_ref"]], "y"
        L = abs(stop[0] - start[0]) if d == "x" else abs(stop[1] - start[1])
        feed = pml[p["side"]] + 60.0
        # Feed_R = 50: every port terminated (the excitation's DC content drains instead of staying on
        # the floating strip, which kept the energy from ever decaying)
        # the measurement plane sits on the zero-thickness lead, 150 um short of the box edge: probes at
        # the edge would cut the thick track copper and measure only part of its current (as m3d.py)
        port = FDTD.AddMSLPort(p["n"], pec, start, stop, d, "z", excite=p["excite"], FeedShift=feed,
                               MeasPlaneShift=L - 150.0, Feed_R=50.0, priority=80)
        P.append((p["n"], port))
    os.makedirs(od, exist_ok=True)
    cells = len(X) * len(Y) * len(m["mesh"]["z"])
    print(f"em3d: {cells} cells, {len(P)} ports", flush=True)
    sim = os.path.join(od, "sim")
    FDTD.Run(sim, cleanup=True, verbose=0, numThreads=int(threads))
    f = np.linspace(0, fmax, 401)
    f[0] = 1e6
    out = {"freq_hz": f.tolist(), "cells": cells, "ports": {}}
    for n, port in P:
        port.CalcPort(sim, f)                    # measured line impedance of the lead
        zref = np.array(port.Z_ref)
        port.CalcPort(sim, f, ref_impedance=50)
        out["ports"][str(n)] = {"a_re": np.real(port.uf_inc).tolist(), "a_im": np.imag(port.uf_inc).tolist(),
                                "b_re": np.real(port.uf_ref).tolist(), "b_im": np.imag(port.uf_ref).tolist(),
                                "zref_re": np.real(zref).tolist()}
    json.dump(out, open(os.path.join(od, "ports.json"), "w"))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else 16)
