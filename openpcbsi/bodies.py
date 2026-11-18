"""Part bodies from the board's own 3D models.

`kicad-cli pcb export glb --no-board-body` writes every footprint's 3D model (STEP or VRML, as KiCad
resolves it) placed on the board, one glTF node per reference. From each part's triangles this module
takes its height above the board, the area it covers, and the volume it encloses (the signed volume of
its closed surfaces). The reflow ranking turns the volume into a heat capacity with one volumetric
heat capacity for moulded plastic and metal alike (RHOC, between epoxy's 1.9 and copper's 3.4 J/cm3K),
which the ranking's body variants then vary. A part whose model file is missing has no entry.
"""
import hashlib
import json
import os
import struct
import subprocess

import numpy as np

RHOC = 2.0e-3          # J/(mm3 K): volumetric heat capacity taken for a part's whole volume (see above)
_COMP = {5120: np.int8, 5121: np.uint8, 5122: np.int16, 5123: np.uint16, 5125: np.uint32, 5126: np.float32}
_NCOMP = {"SCALAR": 1, "VEC2": 2, "VEC3": 3, "VEC4": 4, "MAT4": 16}


def export(board, cache_dir):
    """The board's components as a GLB in the cache (by the board's hash), or None if kicad-cli fails."""
    h = hashlib.sha256(open(board, "rb").read()).hexdigest()[:16]
    os.makedirs(cache_dir, exist_ok=True)
    out = os.path.join(cache_dir, f"bodies_{h}.glb")
    if not os.path.exists(out):
        tmp = out + f".{os.getpid()}.tmp.glb"
        p = subprocess.run(["kicad-cli", "pcb", "export", "glb", "--no-board-body", "-f", "-o", tmp, board],
                           capture_output=True, text=True, timeout=1800)
        if p.returncode != 0 or not os.path.exists(tmp):
            return None
        os.replace(tmp, out)
    return out


def _read(path):
    b = open(path, "rb").read()
    _, _, _ = struct.unpack("<III", b[:12])
    off, j, bin_ = 12, None, None
    while off < len(b):
        ln, typ = struct.unpack("<II", b[off:off + 8])
        chunk = b[off + 8:off + 8 + ln]
        if typ == 0x4E4F534A:
            j = json.loads(chunk)
        elif typ == 0x004E4942:
            bin_ = chunk
        off += 8 + ln
    return j, bin_


def _accessor(j, bin_, k):
    a = j["accessors"][k]
    bv = j["bufferViews"][a["bufferView"]]
    n, comp = a["count"], _NCOMP[a["type"]]
    dt = np.dtype(_COMP[a["componentType"]])
    start = bv.get("byteOffset", 0) + a.get("byteOffset", 0)
    stride = bv.get("byteStride") or comp * dt.itemsize
    if stride == comp * dt.itemsize:
        arr = np.frombuffer(bin_, dtype=dt, count=n * comp, offset=start)
    else:
        arr = np.vstack([np.frombuffer(bin_, dtype=dt, count=comp, offset=start + i * stride) for i in range(n)])
    return arr.reshape(n, comp) if comp > 1 else arr


def _local(node):
    if "matrix" in node:
        return np.array(node["matrix"], float).reshape(4, 4).T
    M = np.eye(4)
    if "scale" in node:
        M = np.diag(list(node["scale"]) + [1.0]) @ M
    if "rotation" in node:
        x, y, z, w = node["rotation"]
        R = np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                      [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                      [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])
        RM = np.eye(4)
        RM[:3, :3] = R
        M = RM @ M
    if "translation" in node:
        T = np.eye(4)
        T[:3, 3] = node["translation"]
        M = T @ M
    return M


def read(path):
    """{reference: {height_mm, volume_mm3, area_mm2, bbox}} of the parts in a GLB written by kicad-cli.
    glTF is in metres with y up; KiCad's board x and y are its x and z (checked on placed parts)."""
    j, bin_ = _read(path)
    nodes = j.get("nodes", [])
    out = {}

    def walk(k, M, ref):
        node = nodes[k]
        M = M @ _local(node)
        name = node.get("name") or ""
        if ref is None and name and not name.startswith("=>"):
            ref = name
        if "mesh" in node and ref:
            e = out.setdefault(ref, {"pts": [], "vol": 0.0})
            for prim in j["meshes"][node["mesh"]]["primitives"]:
                P = _accessor(j, bin_, prim["attributes"]["POSITION"]).astype(float)
                P = (M @ np.c_[P, np.ones(len(P))].T).T[:, :3] * 1000.0          # mm
                e["pts"].append(P)
                if "indices" in prim and prim.get("mode", 4) == 4:
                    I = _accessor(j, bin_, prim["indices"]).astype(np.int64).reshape(-1, 3)
                    a, b_, c = P[I[:, 0]], P[I[:, 1]], P[I[:, 2]]
                    e["vol"] += float(np.einsum("ij,ij->i", a, np.cross(b_, c)).sum() / 6.0)
        for c in node.get("children", []):
            walk(c, M, ref)
    roots = j.get("scenes", [{}])[j.get("scene", 0)].get("nodes", range(len(nodes)))
    for r in roots:
        walk(r, np.eye(4), None)
    res = {}
    for ref, e in out.items():
        P = np.vstack(e["pts"]) if e["pts"] else np.zeros((0, 3))
        if not len(P):
            continue
        # glTF y is up: the part's height is its extent along y; its footprint is x by z
        h = float(P[:, 1].max() - P[:, 1].min())
        w, d = float(P[:, 0].max() - P[:, 0].min()), float(P[:, 2].max() - P[:, 2].min())
        vol = abs(e["vol"])
        res[ref] = {"height_mm": round(h, 3), "volume_mm3": round(vol, 3), "area_mm2": round(w * d, 3),
                    "bbox": [round(float(v), 3) for v in (P[:, 0].min(), P[:, 2].min(), P[:, 0].max(), P[:, 2].max())]}
    return res


def bodies(board, cache_dir):
    """The parts' bodies from their 3D models: read(), from the cached export; {} if there is none."""
    p = export(board, cache_dir)
    if not p:
        return {}
    j = p + ".json"
    if os.path.exists(j):
        try:
            return json.load(open(j))
        except ValueError:
            pass
    r = read(p)
    json.dump(r, open(j, "w"))
    return r
