"""G16: report files. report.html (one page per run), PNG per copper layer with the findings drawn
on the copper, Z(f) charts per rail, and a copy of the board with markers on User.9 (written with
pcbnew; the source board is never modified)."""
import html
import json
import os

import numpy as np

# reference palette (dataviz skill): categorical slots 1-2, status colours, neutral inks
SERIES = ["#2a78d6", "#eb6834"]
STATUS = {"error": "#d03b3b", "warning": "#fab219", "info": "#8a8a85"}
MARK = {"error": "X", "warning": "^", "info": "o"}
INK, INK2, GRID, SURF = "#22221f", "#6b6a63", "#e4e3de", "#fcfcfb"


def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                         "figure.facecolor": SURF, "axes.facecolor": SURF, "axes.spines.top": False,
                         "axes.spines.right": False, "legend.frameon": False})
    return plt


def all_findings(report):
    out = []
    for name, chk in report["checks"].items():
        items = chk if isinstance(chk, list) else [chk]
        for r in items:
            for f in r.get("findings", []):
                out.append(f)
    return out


def layer_images(report, board, ctx, out_dir, step=4):
    """One PNG per copper layer: copper in grey (reference nets darker), exclusion areas outlined,
    findings marked by severity (shape + colour + legend)."""
    plt = _mpl()
    from matplotlib.lines import Line2D
    g = board.grid()
    fs = all_findings(report)
    paths = {}
    areas = ctx.areas()[::step, ::step]
    for li, lname in enumerate(board.layers):
        r = board.raster(li)[::step, ::step]
        img = np.full(r.shape + (3,), 1.0)
        cu = r != 0
        ref = np.isin(r, ctx.ref_ids)
        img[cu & ~ref] = (0.80, 0.79, 0.76)
        img[ref] = (0.62, 0.61, 0.58)
        fig, ax = plt.subplots(figsize=(7.5, 7.5))
        ax.imshow(img, extent=g.extent(), interpolation="nearest")
        if areas.any():
            ax.contour(np.linspace(g.x0, g.x0 + g.nx * g.res, areas.shape[1]),
                       np.linspace(g.y0, g.y0 + g.ny * g.res, areas.shape[0]), areas > 0, levels=[0.5],
                       colors=[INK2], linewidths=0.8, linestyles="dashed")
        for ring in board.outline_rings():
            a = np.asarray(ring + ring[:1])
            ax.plot(a[:, 0], a[:, 1], color=INK, lw=1)
        counts = {}
        for sev in ("info", "warning", "error"):
            sel = [f for f in fs if f["severity"] == sev and f["x"] is not None and lname in (f.get("layer") or "")]
            if not sel:
                continue
            counts[sev] = len(sel)
            for f in sel:
                if f.get("x2") is not None:
                    ax.plot([f["x"], f["x2"]], [f["y"], f["y2"]], color=STATUS[sev], lw=2, solid_capstyle="round")
            ax.scatter([f["x"] for f in sel], [f["y"] for f in sel], marker=MARK[sev], s=36 if sev != "info" else 14,
                       c=STATUS[sev], edgecolors=SURF, linewidths=0.8, zorder=3)
        ax.set_xlim(board.bbox[0] - 1, board.bbox[2] + 1)
        ax.set_ylim(board.bbox[3] + 1, board.bbox[1] - 1)
        ax.set_aspect("equal")
        ax.set_title(f"{lname}: findings on the copper", color=INK, loc="left")
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        handles = [Line2D([], [], marker=MARK[s], ls="", color=STATUS[s], label=f"{s} ({counts[s]})")
                   for s in ("error", "warning", "info") if s in counts]
        handles.append(Line2D([], [], color=INK2, ls="dashed", label="exclusion area"))
        ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=8)
        fig.tight_layout()
        p = os.path.join(out_dir, f"layer_{lname.replace('.', '_')}.png")
        fig.savefig(p, dpi=110)
        plt.close(fig)
        paths[lname] = os.path.basename(p)
    return paths


def z_chart(rail, out_dir, other=None, labels=("A", "B")):
    """|Z(f)| at the load port with the target line (one axis, log-log)."""
    plt = _mpl()
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    f = np.array(rail["freqs_hz"])
    ax.loglog(f, np.array(rail["Z_ohm"]) * 1000, color=SERIES[0], lw=2, label=labels[0] if other else None)
    if other:
        ax.loglog(np.array(other["freqs_hz"]), np.array(other["Z_ohm"]) * 1000, color=SERIES[1], lw=2, label=labels[1])
        ax.legend(loc="upper left")
    zt = rail["target_ohm"] * 1000
    fm = rail.get("target_fmax_hz", f.max())
    ax.plot([f.min(), fm], [zt, zt], color=INK2, lw=1.2, ls="--")
    ax.text(f.min() * 1.15, zt * 1.12, f"target {zt:.1f} mOhm (to {fm/1e6:.0f} MHz)", color=INK2, fontsize=8)
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("|Z| at the load balls (mOhm)")
    ax.set_title(f"{rail['rail']}: PDN impedance", color=INK, loc="left")
    fig.tight_layout()
    from .si import inside, safe_name           # a rail name from a config is any text: never a path
    p = inside(out_dir, os.path.join(out_dir, f"pdn_{safe_name(rail['rail'])}.png"))
    fig.savefig(p, dpi=110)
    plt.close(fig)
    return os.path.basename(p)


def heatmap(T, mesh, board, path, title):
    plt = _mpl()
    g = mesh["g"]
    fig, ax = plt.subplots(figsize=(6.8, 6.0))
    ext = (g.x0, g.x0 + mesh["nx"] * mesh["cs"], g.y0 + mesh["ny"] * mesh["cs"], g.y0)
    im = ax.imshow(T, extent=ext, cmap="Oranges", interpolation="nearest")
    for ring in board.outline_rings():
        a = np.asarray(ring + ring[:1])
        ax.plot(a[:, 0], a[:, 1], color=INK, lw=1)
    ax.set_xlim(board.bbox[0] - 1, board.bbox[2] + 1)
    ax.set_ylim(board.bbox[3] + 1, board.bbox[1] - 1)
    ax.set_aspect("equal")
    ax.grid(False)
    ax.set_title(title, color=INK, loc="left")
    fig.colorbar(im, ax=ax, shrink=0.8, label="C")
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return os.path.basename(path)


def marker_board(report, board_path, out_dir, layer="User.9"):
    """Copy of the board with a circle and a text per finding (errors and warnings) on a user layer
    (`layer` if the board's layer table has it, else Cmts.User). Written as text into a copy of the board
    file, so nothing else in the file changes and no KiCad library is needed."""
    import re
    import uuid
    text = open(board_path, encoding="utf-8").read()
    head = text[:text.find("(setup")] if "(setup" in text else text[:20000]
    if f'"{layer}"' not in head:
        layer = "Cmts.User"
    esc = lambda t: t.replace("\\", "\\\\").replace('"', '\\"').replace("\n", " ")
    items, n = [], 0
    for f in all_findings(report):
        if f["severity"] == "info" or f["x"] is None:
            continue
        x, y = round(float(f["x"]), 4), round(float(f["y"]), 4)
        r = 0.6 if f["severity"] == "error" else 0.4
        items.append(f'\t(gr_circle (center {x} {y}) (end {round(x + r, 4)} {y}) (stroke (width 0.05) (type solid)) (fill no) '
                     f'(layer "{layer}") (uuid "{uuid.uuid4()}"))')
        msg = esc(f"{f['rule']}: {f['net']} {f['msg']}"[:90])
        items.append(f'\t(gr_text "{msg}" (at {round(x + 0.7, 4)} {y} 0) (layer "{layer}") (uuid "{uuid.uuid4()}") '
                     f'(effects (font (size 0.25 0.25) (thickness 0.04)) (justify left)))')
        n += 1
    end = text.rstrip().rfind(")")
    out_text = text[:end].rstrip() + "\n" + "\n".join(items) + "\n)\n"
    name = os.path.splitext(os.path.basename(board_path))[0] + "_openpcbsi.kicad_pcb"
    out = os.path.join(out_dir, name)
    open(out, "w", encoding="utf-8").write(out_text)
    return os.path.basename(out), n


# ----------------------------------------------------------------------------------------------

CSS = """
:root{--surface:#fcfcfb;--ink:#22221f;--ink2:#6b6a63;--rule:#e4e3de;--err:#d03b3b;--warn:#a87000;--ok:#0c7a0c;--chip:#f0efea}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--surface:#1a1a19;--ink:#f0efec;--ink2:#c3c2b7;--rule:#3a3a37;--err:#e66767;--warn:#fab219;--ok:#0ca30c;--chip:#262624}}
:root[data-theme="dark"]{--surface:#1a1a19;--ink:#f0efec;--ink2:#c3c2b7;--rule:#3a3a37;--err:#e66767;--warn:#fab219;--ok:#0ca30c;--chip:#262624}
body{background:var(--surface);color:var(--ink);font:14px/1.45 system-ui,sans-serif;margin:0 auto;max-width:1180px;padding:16px}
h1{font-size:22px;margin:8px 0}h2{font-size:17px;margin:28px 0 8px;border-bottom:1px solid var(--rule);padding-bottom:4px}
h3{font-size:14px;margin:16px 0 6px}.muted{color:var(--ink2)}
table{border-collapse:collapse;width:100%;margin:6px 0 12px;font-size:12.5px}
th,td{border-bottom:1px solid var(--rule);padding:3px 6px;text-align:left;vertical-align:top}
th{color:var(--ink2);font-weight:600}td.n{text-align:right;font-variant-numeric:tabular-nums}
.chip{display:inline-block;padding:0 6px;border-radius:9px;background:var(--chip);font-size:12px}
.error{color:var(--err)}.warning{color:var(--warn)}.ok{color:var(--ok)}
details{margin:4px 0}summary{cursor:pointer;color:var(--ink2)}
img{max-width:100%;height:auto;border:1px solid var(--rule)}
.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(340px,1fr));gap:10px}
.wrap{overflow-x:auto}
"""


def _t(rows, cols, num=()):
    h = ["<div class=wrap><table><tr>" + "".join(f"<th>{html.escape(c)}</th>" for c in cols) + "</tr>"]
    for r in rows:
        h.append("<tr>" + "".join(f"<td class={'n' if c in num else ''}>{html.escape(_fmt(r.get(c)))}</td>" for c in cols) + "</tr>")
    h.append("</table></div>")
    return "".join(h)


def _fmt(v):
    if v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.4g}"
    if isinstance(v, (dict, list)):
        return json.dumps(v)[:160]
    return str(v)


def sev_label(sev):
    icon = {"error": "✖", "warning": "▲", "info": "●"}.get(sev, "")
    return f"<span class={sev}>{icon} {sev}</span>"


def write_html(report, out_dir, images, zcharts):
    c = report["checks"]
    H = [f"<!doctype html><html><head><meta charset=utf-8><meta name=viewport content='width=device-width,initial-scale=1'>"
         f"<title>openpcbsi report</title><style>{CSS}</style></head><body>"]
    H.append(f"<h1>openpcbsi: {html.escape(os.path.basename(report['board']))}</h1>")
    H.append(f"<p class=muted>{html.escape(report['board'])}<br>sha256 {report['board_sha256'][:16]} · {report['date']} · "
             f"KiCad {html.escape(str(report.get('kicad')))} · config {html.escape(str(report.get('config')))} "
             f"({report['config_hash']}) · {report['seconds']} s</p>")
    v = report["validation"]
    H.append("<p>Solver validation (G17): " + " ".join(
        f"<span class=chip>{k}: <span class={'ok' if s == 'validated' else 'warning'}>{s}</span></span>" for k, s in v.items()) + "</p>")
    # summary
    H.append("<h2>Summary</h2>")
    rows = []
    for r in c.get("drc", []):
        s = r.get("summary", {})
        rows.append({"check": f"{r['spec']} {r['rule']}", "what": r["title"], "items": len(r["metrics"]),
                     "errors": s.get("error", 0), "warnings": s.get("warning", 0),
                     "in exclusion areas": sum(v for k, v in s.items() if k.endswith("_in_areas"))})
    H.append(_t(rows, ["check", "what", "items", "errors", "warnings", "in exclusion areas"], num=("items", "errors", "warnings", "in exclusion areas")))
    if images:
        H.append("<h2>Findings on the board</h2><div class=grid>" + "".join(
            f"<a href='{p}'><img src='{p}' alt='findings on {l}' loading=lazy></a>" for l, p in images.items()) + "</div>")
        if report.get("marker_board"):
            H.append(f"<p class=muted>Markers for KiCad: <code>{html.escape(report['marker_board'])}</code> (layer User.9).</p>")
    # DRC details
    if c.get("drc"):
        H.append("<h2>G1 electrical rules</h2>")
        for r in c["drc"]:
            H.append(f"<h3>{r['spec']} {html.escape(r['title'])} <span class=muted>({r['rule']})</span></h3>")
            fs = [f for f in r["findings"] if f["severity"] != "info"]
            fi = [f for f in r["findings"] if f["severity"] == "info"]
            if fs:
                H.append(_t([dict(f, severity_=f["severity"]) for f in fs[:200]],
                            ["severity", "net", "layer", "msg", "x", "y", "area"], num=("x", "y")))
            else:
                H.append("<p class=ok>No errors or warnings.</p>")
            if fi:
                H.append(f"<details><summary>{len(fi)} informational (inside exclusion areas, or listed by the rule)</summary>"
                         + _t(fi[:300], ["net", "layer", "msg", "x", "y", "area"], num=("x", "y")) + "</details>")
            mrows = [dict({"item": k}, **{kk: vv for kk, vv in v.items() if not isinstance(vv, (list, dict))})
                     for k, v in r["metrics"].items()]
            if mrows:
                cols = ["item"] + sorted({kk for m in mrows for kk in m if kk != "item"})
                H.append(f"<details><summary>metrics per item ({len(mrows)})</summary>" + _t(mrows, cols, num=cols[1:]) + "</details>")
    if c.get("impedance"):
        r = c["impedance"]
        H.append(f"<h2>G6 impedance planner <span class=chip>{r['validation']}</span></h2>")
        H.append(_t(r["rows"], ["layer", "kind", "w", "gap", "length_mm", "Z0", "Zdiff", "Zcm", "eeff", "delay_ps_mm",
                                "loss_db_mm", "target", "delta_target", "fab_z", "delta_fab"],
                    num=("w", "gap", "length_mm", "Z0", "Zdiff", "Zcm", "eeff", "delay_ps_mm", "target", "delta_target")))
    if c.get("dc"):
        r = c["dc"]
        H.append(f"<h2>G10 DC drop <span class=chip>engine: {r.get('engine', 'openpcbsi')}</span> <span class=chip>{r['validation']}</span></h2>")
        if r.get("engine") == "kipi":
            H.append("<p class=muted>KiPI (deps/, pinned). drop_mv: regulator target to the worst ball, supply minus ground: "
                     "shunts, shared upstream copper and the ground return included. drop_from_source_mv: from the rail's "
                     "shunt (or the sense point) to the worst ball, supply side only. Currents from the config (estimates unless stated).</p>")
            H.append(_t(r["rails"], ["rail", "net", "regulation", "current_a", "v_nom", "drop_mv", "drop_pct_of_window", "passed",
                                     "drop_from_source_mv", "v_ground_mv", "max_via_current_a", "max_J_A_mm2", "loss_w"],
                        num=("current_a", "v_nom", "drop_mv", "drop_pct_of_window", "drop_from_source_mv", "v_ground_mv",
                             "max_via_current_a", "max_J_A_mm2", "loss_w")))
            if r.get("cross_check"):
                H.append("<h3>Cross-check: KiPI against openpcbsi's own solver (same definition: from the shunt or sense point, supply side)</h3>")
                H.append(_t(r["cross_check"], ["rail", "kipi_mv", "openpcbsi_mv", "diff_pct", "kipi_via_a", "openpcbsi_via_a"],
                            num=("kipi_mv", "openpcbsi_mv", "diff_pct", "kipi_via_a", "openpcbsi_via_a")))
        else:
            H.append("<p class=muted>Supply net only (ground ideal). Currents from the config (estimates unless stated).</p>")
            H.append(_t(r["rails"], ["rail", "net", "regulation", "current_a", "v_nom", "drop_mv", "drop_pct_of_window",
                                     "worst_ball", "max_via_current_a", "max_track_J_A_mm2", "max_track_dT_C", "problems"],
                        num=("current_a", "v_nom", "drop_mv", "drop_pct_of_window", "max_via_current_a", "max_track_J_A_mm2", "max_track_dT_C")))
    if c.get("pdn"):
        r = c["pdn"]
        H.append(f"<h2>G11/G12 PDN impedance <span class=chip>engine: {r.get('engine', 'openpcbsi')}</span> <span class=chip>{r['validation']}</span></h2>")
        if r.get("engine") == "kipi":
            H.append("<p class=muted>|Z| from KiPI's field solve of rail and ground. The decoupling optimiser runs on openpcbsi's own "
                     "network. The charts overlay KiPI (A) and openpcbsi (B).</p>")
        rows = [dict(x, add=(x.get("optimiser") or {}).get("add"), met=(x.get("optimiser") or {}).get("met"),
                     removable=(x.get("optimiser") or {}).get("removable"),
                     resonances_mhz=[round(q["f_hz"] / 1e6, 1) for q in x.get("resonances", [])]) for x in r["rails"]]
        H.append(_t(rows, ["rail", "caps", "target_ohm", "max_Z_ohm", "f_at_max_hz", "max_ratio", "resonances_mhz", "add", "met", "removable"],
                    num=("caps", "target_ohm", "max_Z_ohm", "f_at_max_hz", "max_ratio")))
        H.append("<div class=grid>" + "".join(f"<img src='{html.escape(p, quote=True)}' alt='impedance of {html.escape(n, quote=True)}' "
                                              "loading=lazy>" for n, p in zcharts.items()) + "</div>")
    if c.get("crosstalk"):
        r = c["crosstalk"]
        H.append(f"<h2>G4 crosstalk <span class=chip>{r['validation']}</span></h2>")
        H.append(_t(r["top"], ["victim", "aggressor", "layer", "coupled_mm", "spacing_mm", "rise_ps", "Kb", "NEXT_pct", "FEXT_pct", "area"],
                    num=("coupled_mm", "spacing_mm", "rise_ps", "Kb", "NEXT_pct", "FEXT_pct")))
    if c.get("thermal"):
        r = c["thermal"]
        H.append(f"<h2>G13 thermal <span class=chip>{r['validation']}</span></h2>")
        H.append(f"<p>Total {r['total_w']} W ({r['joule_w']} W I²R of the rails). Board maximum {r['T_max_c']} C.</p>")
        H.append(_t(r["parts"], ["ref", "power_w", "T_board_c", "theta_jb", "Tj_est_c"], num=("power_w", "T_board_c", "theta_jb", "Tj_est_c")))
        if r.get("image"):
            H.append(f"<img src='{html.escape(str(r['image']), quote=True)}' alt='temperature map' loading=lazy style='max-width:560px'>")
    for key in ("vias", "si", "emc"):
        if c.get(key):
            r = c[key]
            H.append(f"<h2>{r.get('spec')} {html.escape(r.get('title', key))} <span class=chip>{r.get('validation', '')}</span></h2>")
            if r.get("rows"):
                cols = list(r["rows"][0].keys())[:14]
                H.append(_t(r["rows"][:80], cols, num=cols))
            if r.get("note"):
                H.append(f"<p class=muted>{html.escape(r['note'])}</p>")
    H.append("</body></html>")
    p = os.path.join(out_dir, "report.html")
    with open(p, "w") as f:
        f.write("".join(H))
    return p


def write_all(report, ctx_tuple, out_dir, markers=True, log=print):
    b, cfg, ctx = ctx_tuple
    images = layer_images(report, b, ctx, out_dir)
    zc = {}
    pchk = report["checks"].get("pdn", {})
    own = {r["rail"]: r for r in pchk.get("openpcbsi_rails", [])}
    for r in pchk.get("rails", []):
        if r.get("Z_ohm"):
            o = own.get(r["rail"]) if pchk.get("engine") == "kipi" else None
            zc[r["rail"]] = z_chart(r, out_dir, o if (o and o.get("Z_ohm")) else None, ("KiPI", "openpcbsi"))
    if markers:
        try:
            name, n = marker_board(report, b.path, out_dir)
            report["marker_board"] = name
            log(f"marker board: {os.path.join(out_dir, name)} ({n} markers on User.9)")
        except Exception as e:
            log(f"marker board not written: {e}")
    from .run import _jd
    from .paths import portable
    plain = portable(json.loads(json.dumps(report, default=_jd)), out_dir)
    p = write_html(plain, out_dir, images, zc)
    with open(os.path.join(out_dir, "report.json"), "w") as f:
        json.dump(plain, f, indent=1)
    log(f"html: {p}")
    return p
