"""Board generator: small KiCad 10 boards built through pcbnew with a known answer by construction (a
clean board and the same board with one planted defect). Used by the acceptance tests and by
`validate --3d`; runs in KiCad's system python3 (pcbnew).

    b = Board("x.kicad_pcb", layers=4, size=(40, 20))
    gnd = b.net("GND")
    b.plane("In1.Cu", gnd)
    b.pair("F.Cu", "LANE0_TX", x0=5, x1=35, y=10, w=0.12, gap=0.15)
    b.slot("In1.Cu", 19.75, 0, 20.25, 20)
    b.save()

The stackup is written into the file as text (KiCad 10's SWIG API does not expose it).
"""
import json
import os
import re

import pcbnew

MM = pcbnew.FromMM

# a symmetric stackup per copper count: (layer name, thickness mm, er, tan d) for dielectrics
STACKUPS = {
    2: [("dielectric 1", "core", 1.51, 4.5, 0.02)],
    4: [("dielectric 1", "prepreg", 0.12, 4.2, 0.02), ("dielectric 2", "core", 1.13, 4.5, 0.02),
        ("dielectric 3", "prepreg", 0.12, 4.2, 0.02)],
    6: [("dielectric 1", "prepreg", 0.12, 4.2, 0.02), ("dielectric 2", "core", 0.3, 4.4, 0.02),
        ("dielectric 3", "prepreg", 0.4, 4.0, 0.02), ("dielectric 4", "core", 0.3, 4.4, 0.02),
        ("dielectric 5", "prepreg", 0.12, 4.2, 0.02)],
}
OUTER_CU, INNER_CU, MASK = 0.035, 0.0175, (0.025, 3.8)

LAYER = {"F.Cu": pcbnew.F_Cu, "B.Cu": pcbnew.B_Cu, "In1.Cu": pcbnew.In1_Cu, "In2.Cu": pcbnew.In2_Cu,
         "In3.Cu": pcbnew.In3_Cu, "In4.Cu": pcbnew.In4_Cu, "In5.Cu": pcbnew.In5_Cu, "In6.Cu": pcbnew.In6_Cu,
         "Edge.Cuts": pcbnew.Edge_Cuts, "F.CrtYd": pcbnew.F_CrtYd, "B.CrtYd": pcbnew.B_CrtYd}


def V(x, y):
    return pcbnew.VECTOR2I(MM(x), MM(y))


class Board:
    def __init__(self, path, layers=4, size=(40.0, 20.0), stackup=None):
        self.path, self.n_layers, self.size = os.path.abspath(path), layers, size
        self.stackup = stackup or STACKUPS[layers]
        self.b = pcbnew.NewBoard(self.path)
        self.b.SetCopperLayerCount(layers)
        self.nets = {}
        self.netclass_patterns = []        # (netclass, pattern) for the .kicad_pro
        self.netclasses = {}               # name -> settings dict for the .kicad_pro
        self.dru = []                      # custom rule texts for the .kicad_dru
        w, h = size
        self.outline([(0, 0), (w, 0), (w, h), (0, h)])

    # ---- primitives ---------------------------------------------------------------------------
    def net(self, name):
        if name not in self.nets:
            n = pcbnew.NETINFO_ITEM(self.b, name)
            self.b.Add(n)
            self.nets[name] = n
        return self.nets[name]

    def outline(self, pts):
        for (x1, y1), (x2, y2) in zip(pts, pts[1:] + pts[:1]):
            s = pcbnew.PCB_SHAPE(self.b)
            s.SetShape(pcbnew.SHAPE_T_SEGMENT)
            s.SetLayer(pcbnew.Edge_Cuts)
            s.SetStart(V(x1, y1))
            s.SetEnd(V(x2, y2))
            s.SetWidth(MM(0.05))
            self.b.Add(s)

    def track(self, layer, pts, w, net):
        ni = self.net(net) if isinstance(net, str) else net
        for (x1, y1), (x2, y2) in zip(pts, pts[1:]):
            t = pcbnew.PCB_TRACK(self.b)
            t.SetLayer(LAYER[layer])
            t.SetWidth(MM(w))
            t.SetNet(ni)
            t.SetStart(V(x1, y1))
            t.SetEnd(V(x2, y2))
            self.b.Add(t)

    def via(self, x, y, net, drill=0.3, size=0.6, kind="through", top="F.Cu", bottom="B.Cu", backdrill_bottom=None):
        v = pcbnew.PCB_VIA(self.b)
        v.SetPosition(V(x, y))
        v.SetDrill(MM(drill))
        v.SetWidth(MM(size))
        v.SetViaType({"through": pcbnew.VIATYPE_THROUGH, "blind": pcbnew.VIATYPE_BLIND,
                      "buried": pcbnew.VIATYPE_BURIED, "micro": pcbnew.VIATYPE_MICROVIA}[kind])
        v.SetLayerPair(LAYER[top], LAYER[bottom])
        v.SetNet(self.net(net) if isinstance(net, str) else net)
        if backdrill_bottom:                     # drilled out from the bottom up to this layer
            v.SetBackdrillMode(pcbnew.BACKDRILL_MODE_BACKDRILL_BOTTOM)
            v.SetBottomBackdrillLayer(LAYER[backdrill_bottom])
            v.SetBottomBackdrillSize(MM(drill + 0.2))
        self.b.Add(v)
        return v

    def zone(self, layer, poly, net, priority=0, clearance=0.2):
        z = pcbnew.ZONE(self.b)
        z.SetLayer(LAYER[layer])
        z.SetNet(self.net(net) if isinstance(net, str) else net)
        z.SetAssignedPriority(priority)
        z.SetLocalClearance(MM(clearance))
        z.SetMinThickness(MM(0.1))
        z.SetPadConnection(pcbnew.ZONE_CONNECTION_FULL)
        o = z.Outline()
        o.NewOutline()
        for x, y in poly:
            o.Append(MM(x), MM(y))
        self.b.Add(z)
        return z

    def plane(self, layer, net, margin=0.5, priority=0):
        w, h = self.size
        return self.zone(layer, [(margin, margin), (w - margin, margin), (w - margin, h - margin), (margin, h - margin)],
                         net, priority)

    def slot(self, layer, x0, y0, x1, y1):
        """A rule area that keeps zone fill out: a gap in the plane of `layer`."""
        z = pcbnew.ZONE(self.b)
        z.SetLayer(LAYER[layer])
        z.SetIsRuleArea(True)
        z.SetDoNotAllowZoneFills(True)
        z.SetDoNotAllowTracks(False)
        z.SetDoNotAllowVias(False)
        z.SetDoNotAllowPads(False)
        z.SetDoNotAllowFootprints(False)
        o = z.Outline()
        o.NewOutline()
        for x, y in ((x0, y0), (x1, y0), (x1, y1), (x0, y1)):
            o.Append(MM(x), MM(y))
        self.b.Add(z)
        return z

    def pair(self, layer, name, x0, x1, y, w=0.12, gap=0.15, pts=None, ends=True):
        """Differential pair <name>_P / <name>_N running along x at centre y (or along `pts`)."""
        off = (w + gap) / 2
        if pts is None:
            pts = [(x0, y), (x1, y)]
        for sfx, s in (("_P", -1), ("_N", 1)):
            self.track(layer, [(x, yy + s * off) for x, yy in pts], w, name + sfx)
        if ends:                    # a driver and a receiver pad at each end (pads no wider than 1.5 x the track)
            pw = min(1.5 * w, 0.18)
            for ref, (x, yy) in ((f"U{name}_A", pts[0]), (f"U{name}_B", pts[-1])):
                self.footprint(ref.replace("_", ""), "IC", x, yy, [
                    {"num": 1, "net": name + "_P", "dx": 0, "dy": -off, "sx": pw},
                    {"num": 2, "net": name + "_N", "dx": 0, "dy": off, "sx": pw}],
                    side="F" if layer == "F.Cu" else "B")
        return name + "_P", name + "_N"

    def footprint(self, ref, value, x, y, pads, fields=None, courtyard=None, side="F", fpid=None):
        """pads: list of dicts {num, net, dx, dy, sx, sy, kind: smd|pth, drill, shape: rect|circle, pintype, func,
        mask}. `mask` (mm) opens the solder mask and the paste over an SMD pad, grown by that margin.
        `fpid` is 'Library:Name' (the manufacturing checks read the package from the name)."""
        fp = pcbnew.FOOTPRINT(self.b)
        fp.SetReference(ref)
        fp.SetValue(value)
        if fpid:
            lib, _, name = fpid.partition(":")
            fp.SetFPID(pcbnew.LIB_ID(lib, name))
            fp.SetAttributes(pcbnew.FP_SMD)
        fp.SetPosition(V(x, y))
        cu = pcbnew.F_Cu if side == "F" else pcbnew.B_Cu
        for d in pads:
            p = pcbnew.PAD(fp)
            p.SetNumber(str(d["num"]))
            kind = d.get("kind", "smd")
            p.SetAttribute(pcbnew.PAD_ATTRIB_SMD if kind == "smd" else pcbnew.PAD_ATTRIB_PTH)
            p.SetShape(pcbnew.PAD_SHAPE_CIRCLE if d.get("shape") == "circle" else pcbnew.PAD_SHAPE_RECTANGLE)
            p.SetSize(pcbnew.VECTOR2I(MM(d["sx"]), MM(d.get("sy", d["sx"]))))
            if kind == "smd":
                ls = pcbnew.LSET()
                ls.AddLayer(cu)
                if d.get("mask") is not None:
                    ls.AddLayer(pcbnew.F_Mask if side == "F" else pcbnew.B_Mask)
                    ls.AddLayer(pcbnew.F_Paste if side == "F" else pcbnew.B_Paste)
                    p.SetLocalSolderMaskMargin(MM(d["mask"]))
                p.SetLayerSet(ls)
            else:
                p.SetDrillSize(pcbnew.VECTOR2I(MM(d.get("drill", 0.3)), MM(d.get("drill", 0.3))))
                p.SetLayerSet(pcbnew.PAD.PTHMask())
            p.SetPosition(V(x + d["dx"], y + d["dy"]))
            if d.get("net"):
                p.SetNet(self.net(d["net"]))
            if d.get("func"):
                p.SetPinFunction(d["func"])
            if d.get("pintype"):
                p.SetPinType(d["pintype"])
            fp.Add(p)
        for k, v in (fields or {}).items():
            fp.SetField(k, str(v))
        if courtyard:
            x0, y0, x1, y1 = courtyard
            for (a, b_), (c, d_) in (((x0, y0), (x1, y0)), ((x1, y0), (x1, y1)), ((x1, y1), (x0, y1)), ((x0, y1), (x0, y0))):
                s = pcbnew.PCB_SHAPE(fp)
                s.SetShape(pcbnew.SHAPE_T_SEGMENT)
                s.SetLayer(pcbnew.F_CrtYd if side == "F" else pcbnew.B_CrtYd)
                s.SetStart(V(x + a, y + b_))
                s.SetEnd(V(x + c, y + d_))
                s.SetWidth(MM(0.05))
                fp.Add(s)
        self.b.Add(fp)
        return fp

    def two_pin(self, ref, value, x, y, net1, net2, pitch=1.0, pad=(0.6, 0.6), fields=None, fpid=None, mask=None):
        """An 0402/0603-like part along x with pads at x -/+ pitch/2."""
        return self.footprint(ref, value, x, y, [
            {"num": 1, "net": net1, "dx": -pitch / 2, "dy": 0, "sx": pad[0], "sy": pad[1], "mask": mask},
            {"num": 2, "net": net2, "dx": pitch / 2, "dy": 0, "sx": pad[0], "sy": pad[1], "mask": mask}],
            fields=fields, courtyard=(-pitch / 2 - 0.5, -0.5, pitch / 2 + 0.5, 0.5), fpid=fpid)

    # ---- project files ------------------------------------------------------------------------
    def netclass(self, name, patterns, **settings):
        """A KiCad netclass (written to the .kicad_pro) with net-name patterns."""
        self.netclasses[name] = settings
        for p in patterns:
            self.netclass_patterns.append((name, p))

    def rule(self, text):
        """A custom rule for the .kicad_dru, e.g. '(rule x (condition "...") (constraint diff_pair_uncoupled (max 2mm)))'."""
        self.dru.append(text)

    # ---- output -------------------------------------------------------------------------------
    def save(self, fill=True, version=None):
        if fill:
            pcbnew.ZONE_FILLER(self.b).Fill(self.b.Zones())
        pcbnew.SaveBoard(self.path, self.b)
        txt = open(self.path).read()
        txt = self._with_stackup(txt)
        if version is not None:
            txt = re.sub(r"\(version \d+\)", f"(version {version})", txt, count=1)
        open(self.path, "w").write(txt)
        self._write_project()
        return self.path

    def _with_stackup(self, txt):
        cu = ["F.Cu"] + [f"In{i}.Cu" for i in range(1, self.n_layers - 1)] + ["B.Cu"]
        L = ['\t\t(stackup',
             '\t\t\t(layer "F.SilkS" (type "Top Silk Screen"))',
             '\t\t\t(layer "F.Paste" (type "Top Solder Paste"))',
             f'\t\t\t(layer "F.Mask" (type "Top Solder Mask") (thickness {MASK[0]}) (epsilon_r {MASK[1]}))']
        for i, name in enumerate(cu):
            t = OUTER_CU if name in ("F.Cu", "B.Cu") else INNER_CU
            L.append(f'\t\t\t(layer "{name}" (type "copper") (thickness {t}))')
            if i < len(self.stackup):
                dn, typ, th, er, tand = self.stackup[i]
                L.append(f'\t\t\t(layer "{dn}" (type "{typ}") (thickness {th}) (material "FR4") '
                         f'(epsilon_r {er}) (loss_tangent {tand}))')
        L += [f'\t\t\t(layer "B.Mask" (type "Bottom Solder Mask") (thickness {MASK[0]}) (epsilon_r {MASK[1]}))',
              '\t\t\t(layer "B.Paste" (type "Bottom Solder Paste"))',
              '\t\t\t(layer "B.SilkS" (type "Bottom Silk Screen"))',
              '\t\t\t(copper_finish "None")', '\t\t\t(dielectric_constraints no)', '\t\t)']
        if "(stackup" in txt:
            return txt
        return txt.replace("\t(setup\n", "\t(setup\n" + "\n".join(L) + "\n", 1)

    def _write_project(self):
        stem = os.path.splitext(self.path)[0]
        classes = [{"name": "Default", "clearance": 0.2, "track_width": 0.2, "via_diameter": 0.6, "via_drill": 0.3,
                    "diff_pair_width": 0.2, "diff_pair_gap": 0.25, "schematic_color": "rgba(0, 0, 0, 0.000)",
                    "pcb_color": "rgba(0, 0, 0, 0.000)", "priority": 2147483647, "wire_width": 6, "bus_width": 12,
                    "line_style": 0, "microvia_diameter": 0.3, "microvia_drill": 0.1}]
        for i, (name, s) in enumerate(self.netclasses.items()):
            c = dict(classes[0], name=name, priority=i)
            c.update(s)
            classes.append(c)
        pro = {"meta": {"filename": os.path.basename(stem) + ".kicad_pro", "version": 3},
               "net_settings": {"classes": classes, "meta": {"version": 4},
                                "netclass_patterns": [{"netclass": n, "pattern": p} for n, p in self.netclass_patterns]},
               "board": {"design_settings": {}}}
        json.dump(pro, open(stem + ".kicad_pro", "w"), indent=2)
        if self.dru:
            open(stem + ".kicad_dru", "w").write("(version 1)\n" + "\n".join(self.dru) + "\n")
