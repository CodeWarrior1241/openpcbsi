"""G6 / G17: quasi-static 2D cross-section solver for single and coupled lines.

Finite volumes on a graded rectilinear grid (x across, z up). Solves Laplace's equation with
complex permittivity, so one solve per frequency gives C and G. Inductance comes from the vacuum
solve (L = mu0 eps0 C0^-1). Conductor loss follows Wheeler's incremental-inductance rule (all copper
surfaces receded by one grid step) with a Hammerstad-Groiss roughness factor. Dielectrics follow
a wideband Djordjevic-Sarkar model fitted to the stackup's er and tan(delta) at er_ref_ghz.

Models supported: trapezoidal conductors (etch factor), conformal solder mask, any number of
dielectric slabs, a ground plane below and/or above. Not modelled: plane holes, nearby coplanar
ground, frequency-dependent current crowding beyond Wheeler's rule.
"""
import cmath
import math

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.special import ellipk

EPS0 = 8.8541878128e-12
MU0 = 4e-7 * math.pi
C0 = 299792458.0
SIGMA_CU = 5.8e7


# ----------------------------------------------------------------------------------------------
# materials

def djordjevic_sarkar(er, tand, f_ref, m1=4.0, m2=12.0):
    """Return eps(f) (complex, eps' - j eps'') fitted so that eps(f_ref) = er (1 - j tand)."""
    w1, w2 = 2 * math.pi * 10 ** m1, 2 * math.pi * 10 ** m2

    def F(f):
        w = 2 * math.pi * f
        return cmath.log10((w2 + 1j * w) / (w1 + 1j * w))
    Fr = F(f_ref)
    # eps = einf + d/(m2-m1) * F ; Im(eps) = -er*tand
    d = (-er * tand) * (m2 - m1) / Fr.imag
    einf = er - (d / (m2 - m1) * Fr).real

    def eps(f):
        return einf + d / (m2 - m1) * F(max(f, 1.0))
    return eps


def skin_rs(f, sigma=SIGMA_CU):
    return math.sqrt(math.pi * f * MU0 / sigma)


def rough_factor(f, rms_um, sigma=SIGMA_CU):
    if rms_um <= 0:
        return 1.0
    delta = 1.0 / math.sqrt(math.pi * f * MU0 * sigma)
    return 1.0 + 2.0 / math.pi * math.atan(1.4 * (rms_um * 1e-6 / delta) ** 2)


# ----------------------------------------------------------------------------------------------
# geometry

class Section:
    """A cross-section. Units mm. z = 0 at the top of the lower reference plane."""

    def __init__(self):
        self.slabs = []        # (z0, z1, er, tand)
        self.planes = []       # ("below", z) | ("above", z)
        self.conductors = []   # dict(xc, wb, wt, z0, t)
        self.mask = None       # dict(t, er, tand, z_surface)
        self.air_above = None  # z of the top of the dielectric (microstrip)

    def add_conductor(self, xc, w, z0, t, etch=0.0):
        wt = max(w - 2 * etch * t, w * 0.2)
        self.conductors.append({"xc": xc, "wb": w, "wt": wt, "z0": z0, "t": t})


def _axis(keys, lo, hi, fine, coarse, fine_lo, fine_hi, grow=0.25):
    keys = sorted(set([lo, hi] + [k for k in keys if lo <= k <= hi]))
    pts = [lo]
    x = lo
    while x < hi:
        d = 0.0 if fine_lo <= x <= fine_hi else min(abs(x - fine_lo), abs(x - fine_hi))
        s = min(coarse, fine + grow * d)
        x = x + s
        pts.append(min(x, hi))
    allp = np.unique(np.round(np.r_[pts, keys], 9))
    keys = sorted(set(np.round(keys, 9).tolist()))
    # drop lines closer than fine/3 to a key line
    keep = []
    karr = np.array(keys)
    for p in allp:
        if np.min(np.abs(karr - p)) < 1e-12:
            keep.append(p)
        elif np.min(np.abs(karr - p)) > fine / 3:
            keep.append(p)
    return np.unique(np.round(np.array(keep), 9))


class Solver:
    def __init__(self, sec, res_um=5.0, recede=0.0):
        self.sec = sec
        self.res = res_um / 1000.0
        self.recede = recede
        self._build()

    def _build(self):
        s, h = self.sec, self.res
        rec = self.recede
        cx = [c["xc"] for c in s.conductors]
        span_lo = min(c["xc"] - c["wb"] / 2 for c in s.conductors)
        span_hi = max(c["xc"] + c["wb"] / 2 for c in s.conductors)
        zb = min([z for k, z in s.planes if k == "below"], default=None)
        za = max([z for k, z in s.planes if k == "above"], default=None)
        zc_lo = min(c["z0"] for c in s.conductors)
        zc_hi = max(c["z0"] + c["t"] for c in s.conductors)
        dists = ([zc_lo - zb] if zb is not None else []) + ([za - zc_hi] if za is not None else [])
        hmax = max(dists + [0.05])
        X = max(1.0, 12 * hmax)
        xlo, xhi = span_lo - X, span_hi + X
        zlo = zb - rec if zb is not None else zc_lo - 12 * hmax
        ztop = za + rec if za is not None else zc_hi + max(1.5, 15 * hmax)
        xkeys = []
        for c in s.conductors:
            for ww in (c["wb"], c["wt"]):
                for e in (c["xc"] - ww / 2, c["xc"] + ww / 2):
                    # sub-cell refinement at the edge (field singularity)
                    xkeys += [e + d * h for d in (-1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0)]
                xkeys += [c["xc"] - ww / 2 + rec, c["xc"] + ww / 2 - rec]
        zkeys = [zlo, ztop]
        for z0, z1, er, td in s.slabs:
            zkeys += [z0, z1]
        for c in s.conductors:
            zkeys += [c["z0"], c["z0"] + c["t"], c["z0"] + rec, c["z0"] + c["t"] - rec]
            for zz in (c["z0"], c["z0"] + c["t"]):
                zkeys += [zz + d * h for d in (-1.0, -0.5, -0.25, 0.25, 0.5, 1.0)]
        if s.mask:
            zkeys += [s.mask["z_surface"] + s.mask["t"], zc_hi + s.mask["t"]]
        fz = 3 * h
        self.xs = _axis(xkeys, xlo, xhi, h, 0.1, span_lo - 2 * hmax * 0 - 0.1, span_hi + 0.1)
        self.zs = _axis(zkeys, zlo, ztop, h, 0.05, max(zlo, zc_lo - 0.6 * hmax), min(ztop, zc_hi + 0.6 * hmax))
        NX, NZ = len(self.xs), len(self.zs)
        X, Z = np.meshgrid(self.xs, self.zs)        # nodes (NZ, NX)
        lab = -np.ones((NZ, NX), dtype=np.int64)    # -1 free, 0 ground, k conductor k
        if zb is not None:
            lab[Z <= zb - rec + 1e-12] = 0
        if za is not None:
            lab[Z >= za + rec - 1e-12] = 0
        for k, c in enumerate(s.conductors, start=1):
            z0, z1 = c["z0"] + rec, c["z0"] + c["t"] - rec
            if c["t"] <= 0:
                z0 = z1 = c["z0"]
            fr = np.clip((Z - c["z0"]) / (c["t"] or 1), 0, 1)
            half = (c["wb"] + (c["wt"] - c["wb"]) * fr) / 2 - rec
            inside = (Z >= z0 - 1e-12) & (Z <= z1 + 1e-12) & (np.abs(X - c["xc"]) <= half + 1e-12)
            lab[inside] = k
        self.lab = lab
        # cell centres
        self.xc = 0.5 * (self.xs[1:] + self.xs[:-1])
        self.zc = 0.5 * (self.zs[1:] + self.zs[:-1])
        self.n = len(s.conductors)

    def eps_cells(self, f):
        """Complex relative permittivity of every cell at frequency f (Hz); f=None: vacuum."""
        s = self.sec
        NZc, NXc = len(self.zc), len(self.xc)
        e = np.ones((NZc, NXc), dtype=complex)
        if f is None:
            return e
        Zc, Xc = np.meshgrid(self.zc, self.xc, indexing="ij")
        for z0, z1, model in self._slab_models():
            m = (Zc >= z0) & (Zc <= z1)
            e[m] = model(f)
        if s.mask:
            mk = s.mask
            ms = mk["model"](f)
            zs_ = mk["z_surface"]
            m = (Zc > zs_) & (Zc <= zs_ + mk["t"])
            for c in s.conductors:
                halfw = c["wb"] / 2 + mk["t"]
                m |= (Zc > c["z0"]) & (Zc <= c["z0"] + c["t"] + mk["t"]) & (np.abs(Xc - c["xc"]) <= halfw)
            e[m & (e == 1)] = ms
        return e

    def _slab_models(self):
        if not hasattr(self, "_sm"):
            self._sm = [(z0, z1, self.sec.models[k]) for k, (z0, z1, er, td) in enumerate(self.sec.slabs)]
        return self._sm

    def capacitance(self, f=None):
        """Maxwell capacitance matrix (F/m) of the n conductors, complex for lossy dielectrics."""
        eps = self.eps_cells(f) * EPS0
        xs, zs, lab = self.xs, self.zs, self.lab
        NZ, NX = lab.shape
        dx = np.diff(xs)
        dz = np.diff(zs)
        # horizontal edge (k, i)-(k, i+1): face height from cells (k-1, i) and (k, i)
        epad = np.zeros((NZ + 1, NX - 1), dtype=complex)
        epad[1:NZ, :] = eps
        dzpad = np.r_[0.0, dz, 0.0]
        gh = (epad[:-1] * dzpad[:-1, None] / 2 + epad[1:] * dzpad[1:, None] / 2) / dx[None, :]     # (NZ, NX-1)
        epad2 = np.zeros((NZ - 1, NX + 1), dtype=complex)
        epad2[:, 1:NX] = eps
        dxpad = np.r_[0.0, dx, 0.0]
        gv = (epad2[:, :-1] * dxpad[None, :-1] / 2 + epad2[:, 1:] * dxpad[None, 1:] / 2) / dz[:, None]  # (NZ-1, NX)
        idx = np.arange(NZ * NX).reshape(NZ, NX)
        a = np.r_[idx[:, :-1].ravel(), idx[:-1, :].ravel()]
        b = np.r_[idx[:, 1:].ravel(), idx[1:, :].ravel()]
        w = np.r_[gh.ravel(), gv.ravel()]
        N = NZ * NX
        A = sp.coo_matrix((np.r_[-w, -w], (np.r_[a, b], np.r_[b, a])), shape=(N, N)).tocsr()
        deg = -np.asarray(A.sum(axis=1)).ravel()
        A = A + sp.diags(deg)
        labf = lab.ravel()
        free = labf < 0
        fi = np.flatnonzero(free)
        Aff = A[fi][:, fi].tocsc()
        lu = spla.splu(Aff)
        Cm = np.zeros((self.n, self.n), dtype=complex)
        for k in range(1, self.n + 1):
            V = np.zeros(N, dtype=complex)
            V[labf == k] = 1.0
            rhs = -(A[fi] @ V)
            V[fi] = lu.solve(rhs)
            Q = A @ V
            for j in range(1, self.n + 1):
                Cm[j - 1, k - 1] = Q[labf == j].sum()
        return Cm


def build_section(stack, layer, ref_layers, widths, gaps, cfg_xsec, mask_er=3.8):
    """Cross-section of `widths` conductors (mm) separated by `gaps` on `layer` of a Stackup.

    Depths are KiCad's (0 at the board top, growing downwards). A layer referenced only from above
    (a bottom microstrip) is mirrored so that its plane is below in the solver frame."""
    cu = stack.cu(layer)
    refs = [stack.cu(r) for r in ref_layers]
    below = [r for r in refs if r["z_top"] > cu["z_top"]]
    above = [r for r in refs if r["z_top"] < cu["z_top"]]
    if not below and not above:
        raise ValueError("no reference plane")
    T = stack.thickness
    flip = not below

    def span(item):          # (shallow, deep) depth of an item in the solver's frame
        a, b = item["z_top"], item["z_bot"]
        return (T - b, T - a) if flip else (a, b)
    c_sh, c_dp = span(cu)
    if flip:
        below, above = above, []
    pb = min(below, key=lambda r: span(r)[0])
    zref = span(pb)[0]                         # top surface of the plane below
    pa = max(above, key=lambda r: span(r)[0]) if above else None
    lim = span(pa)[1] if pa else -1e9          # bottom surface of the plane above
    up = lambda depth: zref - depth
    f_ref = cfg_xsec["er_ref_ghz"] * 1e9
    sec = Section()
    sec.models = []
    for d in stack.dielectrics:
        sh, dp = span(d)
        if dp <= lim + 1e-12 or sh >= zref - 1e-12:
            continue
        sh, dp = max(sh, lim), min(dp, zref)
        sec.slabs.append((up(dp), up(sh), d["er"], d["tand"]))
        sec.models.append(djordjevic_sarkar(d["er"], d["tand"], f_ref))
    sec.planes.append(("below", 0.0))
    if pa:
        sec.planes.append(("above", up(lim)))
    t = cu["t"]
    zc0 = up(c_dp)
    centres, x = [], 0.0
    for k, w in enumerate(widths):
        centres.append(x + w / 2)
        x += w + (gaps[k] if k < len(widths) - 1 else 0)
    centres = np.array(centres) - x / 2
    for xc, w in zip(centres, widths):
        sec.add_conductor(float(xc), w, zc0, t, cfg_xsec.get("etch_factor", 0.0))
    if stack.is_outer(layer) and not pa:
        side = "F" if stack.index[layer] == 0 else "B"
        m = stack.mask.get(side)
        if m and m["t"] > 0:
            er = m.get("er") or mask_er
            sec.mask = {"t": m["t"], "er": er, "z_surface": zc0, "model": djordjevic_sarkar(er, 0.02, f_ref)}
    return sec


def solve_lines(sec, res_um=5.0, freqs_hz=(1e9,), rough_um=0.0, sigma=SIGMA_CU):
    """Line parameters for 1 or 2 conductors at each frequency. Results are cached on disk by a hash of
    the section and the settings (#30: an unchanged geometry is never solved twice)."""
    import hashlib
    import json
    import os
    key = json.dumps({"v": SOLVER_CACHE_VERSION, "slabs": sec.slabs, "planes": sec.planes, "cond": sec.conductors,
                      "mask": sec.mask, "air": sec.air_above, "res": res_um, "f": list(freqs_hz), "rough": rough_um,
                      "sigma": sigma}, sort_keys=True, default=str)
    h = hashlib.sha256(key.encode()).hexdigest()[:24]
    try:
        from .run import CACHE
        cdir = os.path.join(CACHE, "xsec")
    except Exception:
        cdir = None
    if cdir:
        cp = os.path.join(cdir, h + ".json")
        if os.path.exists(cp):
            try:
                return json.load(open(cp))
            except ValueError:
                pass
    out = _solve_lines(sec, res_um, freqs_hz, rough_um, sigma)
    if cdir:
        os.makedirs(cdir, exist_ok=True)
        tmp = f"{cp}.{os.getpid()}.tmp"
        with open(tmp, "w") as f:
            json.dump(out, f, default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o))
        os.replace(tmp, cp)
        out = json.load(open(cp))           # the same types whether cached or not
    return out


SOLVER_CACHE_VERSION = 1


def _solve_lines(sec, res_um=5.0, freqs_hz=(1e9,), rough_um=0.0, sigma=SIGMA_CU):
    """Line parameters for 1 or 2 conductors at each frequency."""
    solver = Solver(sec, res_um)
    C0m = solver.capacitance(None).real
    Lm = MU0 * EPS0 * np.linalg.inv(C0m)
    # Wheeler: recede all copper by one fine step
    dn = res_um / 1000.0
    try:
        sr = Solver(sec, res_um, recede=dn)
        L2 = MU0 * EPS0 * np.linalg.inv(sr.capacitance(None).real)
        dLdn = (L2 - Lm) / (dn * 1e-3)
    except Exception:
        dLdn = None
    n = len(sec.conductors)
    out = {"nodes": int(solver.lab.size), "freqs": []}
    t = min(c["t"] for c in sec.conductors)
    w = min(c["wb"] for c in sec.conductors)
    Rdc = 1.0 / (sigma * w * 1e-3 * max(t, 1e-6) * 1e-3)
    for f in freqs_hz:
        Cc = solver.capacitance(f)
        Cr = Cc.real
        G = -2 * math.pi * f * Cc.imag
        Rs = skin_rs(f, sigma) * rough_factor(f, rough_um, sigma)
        R = (Rs / MU0) * dLdn if dLdn is not None else np.zeros((n, n))
        R = np.sqrt(R ** 2 + np.diag(np.full(n, Rdc)) ** 2 * np.eye(n)) * np.sign(R + 1e-30) if n else R
        row = {"f_hz": f}
        if n == 1:
            L, C, Rv, Gv = Lm[0, 0], Cr[0, 0], R[0, 0], G[0, 0]
            row.update(_mode(L, C, Rv, Gv, f, C0m[0, 0]))
        else:
            # odd / even modes of a (near) symmetric pair
            Lo, Le = Lm[0, 0] - Lm[0, 1], Lm[0, 0] + Lm[0, 1]
            Co, Ce = Cr[0, 0] - Cr[0, 1], Cr[0, 0] + Cr[0, 1]
            Ro, Re = R[0, 0] - R[0, 1], R[0, 0] + R[0, 1]
            Go, Ge = G[0, 0] - G[0, 1], G[0, 0] + G[0, 1]
            C0o, C0e = C0m[0, 0] - C0m[0, 1], C0m[0, 0] + C0m[0, 1]
            odd = _mode(Lo, Co, Ro, Go, f, C0o)
            even = _mode(Le, Ce, Re, Ge, f, C0e)
            row.update({"Zodd": odd["Z0"], "Zeven": even["Z0"], "Zdiff": 2 * odd["Z0"], "Zcm": even["Z0"] / 2,
                        "eeff_odd": odd["eeff"], "eeff_even": even["eeff"],
                        "loss_diff_db_mm": odd["loss_db_mm"], "delay_diff_ps_mm": odd["delay_ps_mm"],
                        "loss_c_db_mm": odd["loss_c_db_mm"], "loss_d_db_mm": odd["loss_d_db_mm"],
                        # coupling (for G4): backward / forward coefficients of conductor 2 onto 1
                        "Kb": 0.25 * (Lm[0, 1] / Lm[0, 0] + (-Cr[0, 1]) / Cr[0, 0]),
                        "Kf_ps_mm": 0.5 * ((-Cr[0, 1]) / Cr[0, 0] - Lm[0, 1] / Lm[0, 0])
                        * math.sqrt(Lm[0, 0] * Cr[0, 0]) * 1e-3 * 1e12})
        out["freqs"].append(row)
    out["L"] = Lm.tolist()
    out["C0"] = C0m.tolist()
    return out


def _mode(L, C, R, G, f, C0_):
    w = 2 * math.pi * f
    Z = cmath.sqrt((R + 1j * w * L) / (G + 1j * w * C))
    gam = cmath.sqrt((R + 1j * w * L) * (G + 1j * w * C))
    Z0lossless = math.sqrt(L / C)
    alpha_c = R / (2 * Z0lossless)
    alpha_d = G * Z0lossless / 2
    db = 20 / math.log(10) * 1e-3
    return {"Z0": Z0lossless, "Zc_re": Z.real, "eeff": C / C0_, "delay_ps_mm": math.sqrt(L * C) * 1e-3 * 1e12,
            "loss_db_mm": gam.real * db, "loss_c_db_mm": alpha_c * db, "loss_d_db_mm": alpha_d * db}


# ----------------------------------------------------------------------------------------------
# closed forms (G17 references)

def hammerstad_jensen(w, h, er):
    """Zero-thickness microstrip: Z0 and eeff (Hammerstad & Jensen 1980)."""
    u = w / h
    a = 1 + math.log((u ** 4 + (u / 52) ** 2) / (u ** 4 + 0.432)) / 49 + math.log(1 + (u / 18.1) ** 3) / 18.7
    b = 0.564 * ((er - 0.9) / (er + 3)) ** 0.053
    ee = (er + 1) / 2 + (er - 1) / 2 * (1 + 10 / u) ** (-a * b)
    fu = 6 + (2 * math.pi - 6) * math.exp(-((30.666 / u) ** 0.7528))
    z01 = 60 * math.log(fu / u + math.sqrt(1 + (2 / u) ** 2))
    return z01 / math.sqrt(ee), ee


def _kratio(k):
    """K(k)/K(k') with modulus k."""
    return ellipk(k * k) / ellipk(1 - k * k)


def cohn_stripline(w, b, er):
    """Zero-thickness centred stripline, exact (Cohn 1954)."""
    k = 1 / math.cosh(math.pi * w / (2 * b))
    return 30 * math.pi / math.sqrt(er) * _kratio(k)


def cohn_coupled_stripline(w, s, b, er):
    """Zero-thickness edge-coupled centred stripline, exact (Cohn 1955): Z0e, Z0o."""
    ke = math.tanh(math.pi * w / (2 * b)) * math.tanh(math.pi * (w + s) / (2 * b))
    ko = math.tanh(math.pi * w / (2 * b)) / math.tanh(math.pi * (w + s) / (2 * b))
    kp = lambda k: math.sqrt(1 - k * k)
    z0e = 30 * math.pi / math.sqrt(er) * _kratio(kp(ke))
    z0o = 30 * math.pi / math.sqrt(er) * _kratio(kp(ko))
    return z0e, z0o


def stripline_alpha_c(w, b, t, er, f, sigma=SIGMA_CU):
    """Wheeler-based stripline conductor loss (Pozar 3.7), Np/m. For sqrt(er) Z0 < 120."""
    Rs = skin_rs(f, sigma)
    z0 = cohn_stripline(w, b, er)
    A = 1 + 2 * w / (b - t) + (b + t) / (math.pi * (b - t)) * math.log((2 * b - t) / t)
    return 2.7e-3 * Rs * er * z0 / (30 * math.pi * (b - t)) * A / 1e-3  # b in mm -> per m


def simple_section(kind, w, h, er, t=0.0, gap=None, b=None, tand=0.0):
    """Canonical sections for validation: 'microstrip' (plane below at h) or 'stripline' (b apart)."""
    sec = Section()
    model = lambda f, er=er, tand=tand: er * (1 - 1j * tand)
    if kind == "microstrip":
        sec.slabs = [(0.0, h, er, tand)]
        sec.models = [model]
        sec.planes = [("below", 0.0)]
        z0 = h
    else:
        sec.slabs = [(0.0, b, er, tand)]
        sec.models = [model]
        sec.planes = [("below", 0.0), ("above", b)]
        z0 = (b - t) / 2
    if gap is None:
        sec.add_conductor(0.0, w, z0, t)
    else:
        sec.add_conductor(-(w + gap) / 2, w, z0, t)
        sec.add_conductor((w + gap) / 2, w, z0, t)
    return sec
