"""Two readings of one board compared item by item (plan phase 2): the file reader against KiCad's own
(the pcbnew golden snapshot of the reference corpus), or the cross-check parser against the primary.

    differences(a, b, fields=None) -> {"counts": {...}, "problems": [str], "compared": {...}}

`fields` limits the comparison to what the second reader can read (KiCadFiles reads a subset).
Tolerances: positions 1 um, sizes 2 um, polygon areas 1 % (KiCad's own polygons approximate arcs to
2 um), zone fill area 0.5 % per net and layer.
"""
import math
from collections import Counter, defaultdict

POS, SIZE, AREA_REL, ZONE_REL = 0.001, 0.002, 0.01, 0.005
TEARDROP_REL = 0.25
# KiCad's mask openings (pcbnew's TransformShapeToPolygon with ERROR_OUTSIDE) are grown by its 2 um arc
# tolerance on custom shapes: an opening may differ by that much per side between the two readings
MASK_GROW = 0.002


def _perim(polys):
    return sum(sum(math.hypot(r[i][0] - r[i - 1][0], r[i][1] - r[i - 1][1]) for i in range(len(r))) for poly in polys or [] for r in poly)


def _edge_tol(polys):
    return _perim(polys) * MASK_GROW * 1.1


def _area(polys):
    tot = 0.0
    for poly in polys or []:
        for k, r in enumerate(poly):
            if len(r) < 3:
                continue
            a = 0.5 * abs(sum(r[i][0] * r[i - 1][1] - r[i - 1][0] * r[i][1] for i in range(len(r))))
            tot += a if k == 0 else -a
    return tot


ALL = ("nets", "tracks", "vias", "zones", "footprints", "pads", "pad_polys", "mask", "paste", "outline", "stackup", "flash")


def differences(a, b, fields=None, limit=12):
    want = set(fields or ALL)
    probs, counts, compared = [], {}, defaultdict(int)

    def add(msg):
        if len(probs) < 400:
            probs.append(msg)
    if "nets" in want:
        na, nb = set(a["nets"]), set(b["nets"])
        counts["nets"] = (len(na), len(nb))
        for n in sorted(na ^ nb)[:limit]:
            add(f"net {n!r} only in {'the first' if n in na else 'the second'} reading")
    if "tracks" in want:
        # matched by distance (1 um), not by rounded keys: a coordinate like 176.37495 rounds to different
        # last digits through KiCad's integer nanometres and through the file's decimal text
        from scipy.spatial import cKDTree
        ta = [t for t in a["tracks"] if not t.get("arc")]
        tb = [t for t in b["tracks"] if not t.get("arc")]
        counts["track segments"] = (len(ta), len(tb))
        compared["tracks"] = len(ta)
        import numpy as np

        def ends(ts):
            return np.array([(t["x1"], t["y1"], t["x2"], t["y2"]) for t in ts]).reshape(-1, 4)
        A, B = ends(ta), ends(tb)
        used = np.zeros(len(tb), bool)
        miss_a = []
        if len(B):
            tree = cKDTree(np.r_[B, B[:, [2, 3, 0, 1]]])
            for i, t in enumerate(ta):
                hit = None
                for j in tree.query_ball_point(A[i], POS * 2):
                    jj = j % len(tb)
                    u = tb[jj]
                    if not used[jj] and u["layer"] == t["layer"] and u["net"] == t["net"] and abs(u["w"] - t["w"]) <= SIZE:
                        hit = jj
                        break
                if hit is None:
                    miss_a.append(t)
                else:
                    used[hit] = True
        else:
            miss_a = ta
        for t in miss_a[:limit]:
            add(f"track {(t['layer'], t['net'], round(t['x1'], 4), round(t['y1'], 4), round(t['x2'], 4), round(t['y2'], 4), t['w'])} only in the first reading")
        for j in np.flatnonzero(~used)[:limit]:
            t = tb[j]
            add(f"track {(t['layer'], t['net'], round(t['x1'], 4), round(t['y1'], 4), round(t['x2'], 4), round(t['y2'], 4), t['w'])} only in the second reading")
        la = sum(math.hypot(t["x2"] - t["x1"], t["y2"] - t["y1"]) for t in a["tracks"] if t.get("arc"))
        lb = sum(math.hypot(t["x2"] - t["x1"], t["y2"] - t["y1"]) for t in b["tracks"] if t.get("arc"))
        if abs(la - lb) > max(0.01, 0.001 * max(la, lb)):
            add(f"arc track length {la:.3f} mm against {lb:.3f} mm")
    if "vias" in want:
        key = lambda v: (round(v["x"] + 1e-7, 3), round(v["y"] + 1e-7, 3))     # vias sit on a coarser grid than tracks
        va, vb = {key(v): v for v in a["vias"]}, {key(v): v for v in b["vias"]}
        counts["vias"] = (len(va), len(vb))
        for k in set(va) | set(vb):
            if k not in va or k not in vb:
                add(f"via at {k} only in {'the first' if k in va else 'the second'} reading")
                continue
            x, y = va[k], vb[k]
            compared["vias"] += 1
            if abs(x["d"] - y["d"]) > SIZE or abs(x["drill"] - y["drill"]) > SIZE or x["net"] != y["net"]:
                add(f"via at {k}: {x['d']}/{x['drill']} {x['net']} against {y['d']}/{y['drill']} {y['net']}")
            if "flash" in want and set(x.get("flashed", [])) != set(y.get("flashed", [])):
                add(f"via at {k}: flashed on {sorted(x.get('flashed', []))} against {sorted(y.get('flashed', []))}")
    if "zones" in want:
        # teardrops are zones that KiCad regenerates when it refills (kicad-cli does, as the editor does);
        # where a net and layer have teardrops, their area may differ by up to TEARDROP_REL of it
        za, zb, td = defaultdict(float), defaultdict(float), defaultdict(float)
        for z in a["zones"]:
            if not z.get("rule_area"):
                za[(z["net"], z["layer"])] += _area(z["polys"])
                if z.get("teardrop"):
                    td[(z["net"], z["layer"])] = max(td[(z["net"], z["layer"])], 0.0) + _area(z["polys"])
        tdb = defaultdict(float)
        for z in b["zones"]:
            if not z.get("rule_area"):
                zb[(z["net"], z["layer"])] += _area(z["polys"])
                if z.get("teardrop"):
                    tdb[(z["net"], z["layer"])] += _area(z["polys"])
        counts["zone fill (net, layer)"] = (len(za), len(zb))
        for k in set(za) | set(zb):
            x, y = za.get(k, 0.0), zb.get(k, 0.0)
            compared["zone fills"] += 1
            if abs(x - y) > max(0.01, ZONE_REL * max(x, y), TEARDROP_REL * max(td.get(k, 0.0), tdb.get(k, 0.0))):
                add(f"zone fill {k}: {x:.3f} mm2 against {y:.3f} mm2")
    if "footprints" in want:
        fa, fb = {f["ref"]: f for f in a["footprints"]}, {f["ref"]: f for f in b["footprints"]}
        counts["footprints"] = (len(fa), len(fb))
        for r in set(fa) | set(fb):
            if r not in fa or r not in fb:
                add(f"footprint {r} only in {'the first' if r in fa else 'the second'} reading")
                continue
            x, y = fa[r], fb[r]
            compared["footprints"] += 1
            if math.hypot(x["x"] - y["x"], x["y"] - y["y"]) > POS or abs(((x["rot"] - y["rot"]) + 180) % 360 - 180) > 1e-3 or x["side"] != y["side"]:
                add(f"footprint {r}: at ({x['x']:.4f}, {x['y']:.4f}) {x['rot']} {x['side']} against ({y['x']:.4f}, {y['y']:.4f}) {y['rot']} {y['side']}")
            for k in ("dnp", "smd", "placed", "in_bom", "board_only"):
                if k in x and k in y and bool(x[k]) != bool(y[k]):
                    add(f"footprint {r}: {k} {x[k]} against {y[k]}")
            if x.get("value") != y.get("value") or x.get("fpid") != y.get("fpid"):
                add(f"footprint {r}: value/footprint {x.get('value')}/{x.get('fpid')} against {y.get('value')}/{y.get('fpid')}")
    if "pads" in want:
        key = lambda p: (p["ref"], str(p["num"]), round(p["x"], 4), round(p["y"], 4))
        pa, pb = defaultdict(list), defaultdict(list)
        for p in a["pads"]:
            pa[key(p)].append(p)
        for p in b["pads"]:
            pb[key(p)].append(p)
        counts["pads"] = (len(a["pads"]), len(b["pads"]))
        for k in set(pa) | set(pb):
            if len(pa.get(k, [])) != len(pb.get(k, [])):
                add(f"pad {k}: {len(pa.get(k, []))} against {len(pb.get(k, []))}")
                continue
            order = lambda p: (p.get("kind", ""), p.get("anchor", [0, 0])[0], p.get("anchor", [0, 0])[1],
                               tuple(sorted(p.get("layers") or [])) if "flash" in want else (),
                               tuple(sorted(p.get("mask") or {})) if "mask" in want else (), _area(p.get("poly")))
            for x, y in zip(sorted(pa[k], key=order), sorted(pb[k], key=order)):
                compared["pads"] += 1
                if x["kind"] != y["kind"] or x["net"] != y["net"]:
                    add(f"pad {k}: {x['kind']} {x['net']} against {y['kind']} {y['net']}")
                stol = SIZE + (2 * MASK_GROW if x.get("shape") == "custom" else 0.0)
                if abs(x["sx"] - y["sx"]) > stol or abs(x["sy"] - y["sy"]) > stol:
                    add(f"pad {k} ({x.get('shape')}): size {x['sx']:.4f} x {x['sy']:.4f} against {y['sx']:.4f} x {y['sy']:.4f}")
                if abs(x.get("drill", 0) - y.get("drill", 0)) > SIZE or abs((x.get("drill_y") or 0) - (y.get("drill_y") or 0)) > SIZE:
                    add(f"pad {k}: drill {x.get('drill')} x {x.get('drill_y')} against {y.get('drill')} x {y.get('drill_y')}")
                if "flash" in want and set(x.get("layers", [])) != set(y.get("layers", [])):
                    add(f"pad {k}: copper on {x['layers']} against {y['layers']}")
                if "pad_polys" in want:
                    ax, ay = _area(x["poly"]), _area(y["poly"])
                    # KiCad's own polygons lie up to 2 um inside an arc: about 1 % of a 0.45 mm round pad
                    if abs(ax - ay) > max(1e-4, AREA_REL * max(ax, ay), _perim(y["poly"]) * MASK_GROW):
                        add(f"pad {k} ({x.get('shape')}): copper area {ax:.4f} against {ay:.4f} mm2")
                if "mask" in want:
                    for m in set(x.get("mask", {})) | set(y.get("mask", {})):
                        ax, ay = _area(x.get("mask", {}).get(m, {}).get("polys")), _area(y.get("mask", {}).get(m, {}).get("polys"))
                        # pcbnew's mask polygons of custom pads with stroked primitives are grown further than KiCad
                        # plots them (checked against KiCad's own Gerber): up to 8 um per side is allowed there
                        grow_ = 0.008 if x.get("shape") == "custom" else MASK_GROW
                        tol_ = _perim(y.get("mask", {}).get(m, {}).get("polys")) * grow_ * 1.1
                        if abs(ax - ay) > max(1e-4, AREA_REL * max(ax, ay), tol_):
                            add(f"pad {k} ({x.get('shape')}): {m} opening {ax:.4f} against {ay:.4f} mm2")
                if "paste" in want:
                    for m in set(x.get("paste_poly", {})) | set(y.get("paste_poly", {})):
                        ax, ay = _area(x.get("paste_poly", {}).get(m)), _area(y.get("paste_poly", {}).get(m))
                        if abs(ax - ay) > max(1e-4, AREA_REL * max(ax, ay)):
                            add(f"pad {k} ({x.get('shape')}): {m} aperture {ax:.4f} against {ay:.4f} mm2")
    if "outline" in want:
        oa, ob = _area(a["outline"]), _area(b["outline"])
        counts["outline mm2"] = (round(oa, 3), round(ob, 3))
        if abs(oa - ob) > max(0.01, 0.0005 * max(oa, ob)):
            add(f"board outline {oa:.3f} mm2 against {ob:.3f} mm2")
    if "stackup" in want:
        ta = [(l["name"], l.get("thickness")) for l in a["stackup"]]
        tb = [(l["name"], l.get("thickness")) for l in b["stackup"]]
        if ta != tb:
            add(f"stackup differs: {ta} against {tb}")
        if abs((a.get("board_thickness") or 0) - (b.get("board_thickness") or 0)) > SIZE:
            add(f"board thickness {a.get('board_thickness')} against {b.get('board_thickness')}")
    return {"counts": counts, "problems": probs, "compared": dict(compared)}
