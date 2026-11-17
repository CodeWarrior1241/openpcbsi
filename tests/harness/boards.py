"""Fixture boards for the acceptance tests. Each function builds one board (and its project files) and
returns its path. Coordinates in mm. Every board is small so that a check runs in seconds."""
import os

from gen import Board


def _p(out, name):
    d = os.path.join(out, name)
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, name + ".kicad_pcb")


def pair_over_plane(out, name="pair", slot=False, fill=True, layers=4, w=0.12, gap=0.15, net="LANE0_TX",
                    version=None, netclass=None, dru=None, length=30.0):
    """A pair on F.Cu over GND planes on In1/In2; `slot` cuts a 0.5 mm gap across In1 at x = 20."""
    b = Board(_p(out, name), layers=layers, size=(length + 10, 20))
    if layers >= 4:
        b.plane("In1.Cu", "GND")
        b.plane("In2.Cu", "GND")
    if slot:
        b.slot("In1.Cu", length / 2 + 4.75, 0, length / 2 + 5.25, 20)
    b.pair("F.Cu", net, 5, 5 + length, 10, w=w, gap=gap)
    if netclass:
        b.netclass(*netclass[:2], **(netclass[2] if len(netclass) > 2 else {}))
    for r in dru or []:
        b.rule(r)
    return b.save(fill=fill, version=version)


def line(out, name="line", w=0.22, layer="F.Cu", length=30.0, net="LINE0", layers=4):
    """One single-ended track over GND planes on In1/In2, with a pad at each end."""
    b = Board(_p(out, name), layers=layers, size=(length + 10, 20))
    b.plane("In1.Cu", "GND")
    b.plane("In2.Cu", "GND")
    b.track(layer, [(5, 10), (5 + length, 10)], w, net)
    pw = min(1.5 * w, 0.3)
    for ref, x in (("U1", 5), ("U2", 5 + length)):
        b.footprint(ref, "IC", x, 10, [{"num": 1, "net": net, "dx": 0, "dy": 0, "sx": pw}])
    return b.save()


def pair_uncoupled(out, name, net, split_len=3.0, layers=4, netclass=None, dru=None):
    """A pair whose two tracks separate (0.6 mm apart) for `split_len` mm in the middle."""
    b = Board(_p(out, name), layers=layers, size=(40, 20))
    b.plane("In1.Cu", "GND")
    b.plane("In2.Cu", "GND")
    w, gap, off = 0.12, 0.15, 0.135
    xa, xb = 20 - split_len / 2, 20 + split_len / 2
    for sfx, s in (("_P", -1), ("_N", 1)):
        y0, y1 = 10 + s * off, 10 + s * (off + 0.6)
        b.track("F.Cu", [(5, y0), (xa - 0.5, y0), (xa, y1), (xb, y1), (xb + 0.5, y0), (35, y0)], w, net + sfx)
    if netclass:
        b.netclass(*netclass[:2], **(netclass[2] if len(netclass) > 2 else {}))
    for r in dru or []:
        b.rule(r)
    return b.save()


def planeless(out, name="planeless"):
    """2-layer board, a pair and a single-ended net, no planes at all."""
    b = Board(_p(out, name), layers=2, size=(40, 20))
    b.pair("F.Cu", "LANE0_TX", 5, 35, 8)
    b.track("F.Cu", [(5, 14), (35, 14)], 0.2, "SIG1")
    b.track("B.Cu", [(5, 17), (35, 17)], 0.5, "GND")
    return b.save()


def rail_board(out, name="rail", current_field="1.5A", second_rail=True, unrouted=False):
    """+5V entry J1 -> regulator U1 (TPS62xxx) -> +1V8 -> load U2 (4 pins); optional +3V3 from U3 to U4.
    U2 carries a `Current` field when `current_field` is set."""
    b = Board(_p(out, name), layers=4, size=(40, 30))
    b.plane("In1.Cu", "GND")
    b.plane("In2.Cu", "GND")
    b.footprint("J1", "Conn_01x02", 3, 15, [
        {"num": 1, "net": "+5V", "dx": 0, "dy": -1.27, "sx": 1.7, "kind": "pth", "drill": 1.0},
        {"num": 2, "net": "GND", "dx": 0, "dy": 1.27, "sx": 1.7, "kind": "pth", "drill": 1.0}])
    def reg(ref, x, y, vin, vout):
        b.footprint(ref, "TPS62133", x, y, [
            {"num": 1, "net": vin, "dx": -1, "dy": -1, "sx": 0.6, "func": "VIN", "pintype": "power_in"},
            {"num": 2, "net": "GND", "dx": -1, "dy": 1, "sx": 0.6, "func": "GND", "pintype": "power_in"},
            {"num": 3, "net": vout, "dx": 1, "dy": -1, "sx": 0.6, "func": "VOUT", "pintype": "power_out"},
            {"num": 4, "net": "GND", "dx": 1, "dy": 1, "sx": 0.6, "func": "GND", "pintype": "power_in"}],
            courtyard=(-1.8, -1.8, 1.8, 1.8))
    def load(ref, x, y, net, field):
        pads = [{"num": i + 1, "net": net, "dx": -1.5 + i, "dy": -1, "sx": 0.5} for i in range(4)]
        pads += [{"num": i + 5, "net": "GND", "dx": -1.5 + i, "dy": 1, "sx": 0.5} for i in range(4)]
        b.footprint(ref, "FPGA_LOAD", x, y, pads, fields=({"Current": field} if field else None),
                    courtyard=(-2.5, -2, 2.5, 2))
    reg("U1", 12, 10, "+5V", "+1V8")
    load("U2", 30, 10, "+1V8", current_field)
    if not unrouted:
        b.track("F.Cu", [(3, 13.73), (11, 13.73), (11, 9)], 0.5, "+5V")
        b.track("F.Cu", [(13, 9), (18, 9)], 0.3, "+1V8")          # through two vias on B.Cu (via current)
        b.via(18, 9, "+1V8")
        b.track("B.Cu", [(18, 9), (24, 9)], 0.3, "+1V8")
        b.via(24, 9, "+1V8")
        b.track("F.Cu", [(24, 9), (28.5, 9)], 0.3, "+1V8")
        for i in range(4):
            b.track("F.Cu", [(28.5, 9), (28.5 + i, 9)], 0.3, "+1V8")
        for (x, y) in ((11, 11), (13, 11), (28.5, 11), (29.5, 11), (30.5, 11), (31.5, 11), (3, 16.27)):
            b.via(x, y + 0.8 if y == 11 else y, "GND")
            b.track("F.Cu", [(x, y), (x, y + 0.8 if y == 11 else y)], 0.3, "GND")
    if second_rail:
        reg("U3", 12, 22, "+5V", "+3V3")
        load("U4", 30, 22, "+3V3", None)
        if not unrouted:
            b.track("F.Cu", [(11, 13.73), (11, 21)], 0.5, "+5V")
            b.track("F.Cu", [(13, 21), (28.5, 21)], 0.3, "+3V3")
            for i in range(4):
                b.track("F.Cu", [(28.5, 21), (28.5 + i, 21)], 0.3, "+3V3")
    return b.save()


def thermal_part(out, name="thermal"):
    b = Board(_p(out, name), layers=4, size=(30, 30))
    b.plane("In1.Cu", "GND")
    b.plane("In2.Cu", "GND")
    b.footprint("U1", "HOT_IC", 15, 15, [{"num": i + 1, "net": "GND", "dx": -1.5 + i, "dy": 0, "sx": 0.6}
                                        for i in range(4)], fields={"Tj_max": "85"}, courtyard=(-2.5, -2, 2.5, 2))
    return b.save()


def via_grid(out, name="stitch", pitch=3.0, void=(15, 6, 27, 14)):
    """Two GND planes stitched by a via grid at `pitch`, with no vias inside `void` (x0, y0, x1, y1)."""
    b = Board(_p(out, name), layers=4, size=(42, 21))
    b.plane("In1.Cu", "GND")
    b.plane("In2.Cu", "GND")
    x = pitch / 2
    n = 0
    while x < 42:
        y = pitch / 2
        while y < 21:
            if not (void[0] <= x <= void[2] and void[1] <= y <= void[3]):
                b.via(x, y, "GND")
                n += 1
            y += pitch
        x += pitch
    return b.save()


def topology(out, name="topo"):
    """A T net (driver U1 -> two loads) and a fly-by net (driver -> three loads in a row)."""
    b = Board(_p(out, name), layers=4, size=(60, 30))
    b.plane("In1.Cu", "GND")
    b.plane("In2.Cu", "GND")
    one = lambda ref, x, y, net: b.footprint(ref, "R", x, y, [{"num": 1, "net": net, "dx": 0, "dy": 0, "sx": 0.6}])
    one("U1", 5, 10, "TNET")
    one("U2", 40, 5, "TNET")
    one("U3", 40, 15, "TNET")
    b.track("F.Cu", [(5, 10), (25, 10)], 0.2, "TNET")
    b.track("F.Cu", [(25, 10), (25, 5), (40, 5)], 0.2, "TNET")
    b.track("F.Cu", [(25, 10), (25, 15), (40, 15)], 0.2, "TNET")
    one("U4", 5, 24, "FLYBY")
    for i, x in enumerate((25, 37, 49)):
        one(f"U{5 + i}", x, 24, "FLYBY")
    b.track("F.Cu", [(5, 24), (55, 24)], 0.2, "FLYBY")
    return b.save()


def big_board(out, name="big", layers=6):
    """100 x 100 mm, 6 layers, planes on In1..In4, 40 pairs across: for timing budgets."""
    b = Board(_p(out, name), layers=layers, size=(100, 100))
    for l in ("In1.Cu", "In2.Cu", "In3.Cu", "In4.Cu"):
        b.plane(l, "GND")
    for i in range(40):
        b.pair("F.Cu" if i % 2 == 0 else "B.Cu", f"LANE{i}_TX", 5, 95, 3 + i * 2.35)
    return b.save()


# ---- manufacturing review (`openpcbsi dfm`, items 37-46) --------------------------------------
C0402 = "Capacitor_SMD:C_0402_1005Metric"
R0402 = "Resistor_SMD:R_0402_1005Metric"
# KiCad board-setup rules as strict as JLCPCB's figures for every measurand each rule covers
FAB_RULES = {"min_clearance": 0.10, "min_track_width": 0.09, "min_via_annular_width": 0.05, "min_via_diameter": 0.25,
             "min_through_hole_diameter": 0.15, "min_hole_to_hole": 0.45, "min_hole_clearance": 0.30,
             "min_copper_edge_clearance": 0.2, "solder_mask_min_width": 0.1}


def set_rules(board, rules, severities=None):
    """Write board-setup design rules into the board's .kicad_pro (the generator leaves them unset)."""
    import json
    pro = os.path.splitext(board)[0] + ".kicad_pro"
    d = json.load(open(pro))
    ds = d.setdefault("board", {}).setdefault("design_settings", {})
    ds["rules"] = dict(rules)
    ds["rule_severities"] = dict(severities or {})
    json.dump(d, open(pro, "w"), indent=2)
    return board


def acid_board(out, name="acid", planted=True):
    """Copper geometry with a planted defect each, or none:
    a 45 degree track exit from a pad; two tracks of one net 0.15 mm apart; two tracks of different
    nets 0.07 mm apart; a track 0.08 mm from a pad; a track 0.06 mm from a pad that has no net."""
    b = Board(_p(out, name), layers=4, size=(40, 20))
    b.plane("In1.Cu", "GND")
    b.plane("In2.Cu", "GND")
    b.footprint("U1", "PAD", 10, 10, [{"num": 1, "net": "SIG", "dx": 0, "dy": 0, "sx": 2.0}])
    if planted:
        b.track("F.Cu", [(10.5, 9.5), (14.5, 13.5)], 0.2, "SIG")
    else:
        b.track("F.Cu", [(10.5, 10), (15, 10)], 0.2, "SIG")
    pitch = 0.45 if planted else 0.65
    b.track("F.Cu", [(23, 5), (20, 5), (20, 5 + pitch), (23, 5 + pitch)], 0.3, "SIG2")
    b.track("F.Cu", [(26, 14), (34, 14)], 0.2, "SP_A")
    b.track("F.Cu", [(26, 14.27 if planted else 14.5), (34, 14.27 if planted else 14.5)], 0.2, "SP_B")
    b.footprint("U2", "PAD", 10, 17, [{"num": 1, "net": "PT_A", "dx": 0, "dy": 0, "sx": 1.0}])
    y = 16.5 - (0.08 if planted else 0.3) - 0.075
    b.track("F.Cu", [(8, y), (12, y)], 0.15, "PT_B")
    b.footprint("FID1", "Fiducial", 30, 5, [{"num": 1, "net": "", "dx": 0, "dy": 0, "sx": 1.0, "shape": "circle"}])
    x = 30.5 + (0.06 if planted else 0.3) + 0.075
    b.track("F.Cu", [(x, 3), (x, 7)], 0.15, "NN")
    return b.save()


def fab_board(out, name="fab", layers=4, w=0.08, copper=None, inner_track=None):
    """A track of width `w` on F.Cu; with `inner_track`, a track of that width on In1.Cu (its copper set by
    `copper`, e.g. {"In1.Cu": 0.07} for 2 oz)."""
    b = Board(_p(out, name), layers=layers, size=(40, 20), copper=copper)
    if layers >= 4:
        b.plane("In1.Cu", "GND")
        b.plane("In2.Cu", "GND")
    b.track("F.Cu", [(5, 10), (35, 10)], w, "SIG")
    b.track("F.Cu", [(5, 14), (35, 14)], 0.2, "SIG3")
    if inner_track:
        b.track("In1.Cu", [(5, 6), (35, 6)], inner_track, "SIG4")
    return b.save()


def mask_board(out, name="mask", planted=True):
    """Pad pairs with a planted solder-mask defect each (see the `kind` in the comments), or none."""
    b = Board(_p(out, name), layers=4, size=(60, 20))
    b.plane("In1.Cu", "GND")
    b.plane("In2.Cu", "GND")

    def pair(ref, x, pitch, margin):
        b.footprint(ref, "PADS", x, 10, [
            {"num": 1, "net": ref + "_A", "dx": -pitch / 2, "dy": 0, "sx": 0.5, "mask": margin},
            {"num": 2, "net": ref + "_B", "dx": pitch / 2, "dy": 0, "sx": 0.5, "mask": margin}])
    if planted:
        pair("U1", 8, 0.62, 0.0)        # pad_gap: copper 0.12 mm apart (< 0.15)
        pair("U2", 16, 0.75, 0.09)      # web: 0.25 mm apart, openings leave 0.07 mm (< 0.10)
        pair("U3", 24, 0.70, 0.12)      # no_dam: 0.20 mm apart, openings overlap
    else:
        pair("U1", 8, 0.80, 0.0)
        pair("U2", 16, 0.90, 0.05)
        pair("U3", 24, 1.00, 0.05)
    # one pad with a wide opening and a track of another net beside it
    b.footprint("U4", "PAD", 34, 10, [{"num": 1, "net": "U4_A", "dx": 0, "dy": 0, "sx": 1.0, "mask": 0.15 if planted else 0.0}])
    b.track("F.Cu", [(32, 10.65 if planted else 11.2), (36, 10.65 if planted else 11.2)], 0.1, "OTHER")   # 0.10 mm from the pad
    b.footprint("U5", "PAD", 42, 10, [{"num": 1, "net": "U5_A", "dx": 0, "dy": 0, "sx": 1.0, "mask": 0.05 if planted else 0.0}])
    b.track("F.Cu", [(40, 10.67 if planted else 11.2), (44, 10.67 if planted else 11.2)], 0.1, "OTHER2")  # 0.07 mm from the opening
    # a via in an SMD pad (planted) or beside it on a track
    b.footprint("U6", "PAD", 52, 10, [{"num": 1, "net": "GND", "dx": 0, "dy": 0, "sx": 1.0, "mask": 0.0}])
    if planted:
        b.via(52, 10, "GND", drill=0.3, size=0.5)
    else:
        b.track("F.Cu", [(52, 10), (52, 12)], 0.3, "GND")
        b.via(52, 12, "GND", drill=0.3, size=0.6)
    return b.save()


def drill_board(out, name="drill", planted=True):
    """Holes with a planted defect each, or none: a 0.1 mm via hole; a via ring of 0.03 mm; two holes
    0.12 mm apart wall to wall (0.42 mm centre to centre: noted, not failed); a via hole 0.12 mm from a
    track; a through-hole pad with a 0.1 mm ring; copper 0.1 mm from the board edge; two vias 0.18 mm
    apart centre to centre; two vias that overlap; a slot whose end is 0.25 mm from a track; two
    castellated holes too small (0.4 mm) and too close (0.48 mm edge to edge); copper 1 um beyond the edge clearance (marginal)."""
    b = Board(_p(out, name), layers=4, size=(40, 20))
    b.plane("In1.Cu", "GND")
    b.plane("In2.Cu", "GND")
    p = planted
    b.via(6, 6, "A", drill=0.1 if p else 0.3, size=0.4 if p else 0.6)
    b.via(10, 6, "B", drill=0.3, size=0.36 if p else 0.6)
    b.via(14, 6, "GND", drill=0.3, size=0.5)
    b.via(14.42 if p else 15.2, 6, "GND", drill=0.3, size=0.5)
    b.via(20, 6, "C", drill=0.3, size=0.45)
    x = 20.15 + (0.12 if p else 0.4) + 0.075
    b.track("F.Cu", [(x, 3), (x, 9)], 0.15, "D")
    b.footprint("J1", "PIN", 30, 10, [{"num": 1, "net": "E", "dx": 0, "dy": 0, "sx": 1.0 if p else 1.4, "kind": "pth", "drill": 0.8,
                                       "shape": "circle"}])
    y = 0.2 if p else 0.6
    b.track("F.Cu", [(25, y), (35, y)], 0.2, "F")
    # hole spacing judged centre to centre (the owner's reading of JLCPCB's 0.2 mm): 0.15 mm vias
    # 0.18 mm apart centre to centre fail; two 0.3 mm vias cutting into each other fail as such
    b.via(6, 12, "GND", drill=0.15, size=0.3)
    b.via(6.18 if p else 7.0, 12, "GND", drill=0.15, size=0.3)
    b.via(10, 12, "GND", drill=0.3, size=0.5)
    b.via(10.25 if p else 11.0, 12, "GND", drill=0.3, size=0.5)
    # a plated slot (0.6 x 1.6 mm, along y) with a track of another net 0.25 mm beyond the end of its
    # hole: a round hole of the slot's width at the pad centre would put the track 0.75 mm away
    b.footprint("J2", "SLOT", 15, 14, [{"num": 1, "net": "G", "dx": 0, "dy": 0, "sx": 1.0, "sy": 2.0, "kind": "pth", "drill": 0.6,
                                        "drill_y": 1.6}])
    yt = 14 + 0.8 + (0.25 if p else 0.6) + 0.075
    b.track("F.Cu", [(13.8, yt), (16.2, yt)], 0.15, "H")
    # castellated holes (plated holes the board edge cuts): 0.4 mm holes 0.4 mm apart edge to edge,
    # below JLCPCB's 0.5 mm for both
    for k in range(2):
        b.footprint(f"J{3 + k}", "CAST", 40, 8 + k * (0.88 if p else 1.4), [{"num": 1, "net": f"K{k}", "dx": 0, "dy": 0, "sx": 0.71 if p else 1.0,
                                                                             "kind": "pth", "drill": 0.4 if p else 0.6, "shape": "circle"}])
    # copper 1 um beyond the 0.2 mm edge clearance: inside the measuring slack, reported as marginal
    b.track("B.Cu", [(5, 0.301 if p else 0.6), (10, 0.301 if p else 0.6)], 0.2, "M")
    return b.save()


def chip_board(out, name="chips"):
    """Three 0402 capacitors on GND planes, across the belt (pads along y). C1: pad 1 on a thin track, pad 2
    on a short track to a via into the planes. C2: both pads on thin tracks (symmetric). C3: pad 1 in a
    pour, pad 2 on a thin track. A fiducial without a net, so that copper without a net is in the files."""
    b = Board(_p(out, name), layers=4, size=(60, 30))
    b.plane("In1.Cu", "GND")
    b.plane("In2.Cu", "GND")

    def chip(ref, x, n1, n2):
        b.footprint(ref, "100n", x, 15, [
            {"num": 1, "net": n1, "dx": 0, "dy": -0.5, "sx": 0.6, "mask": 0.0},
            {"num": 2, "net": n2, "dx": 0, "dy": 0.5, "sx": 0.6, "mask": 0.0}], courtyard=(-0.6, -1.1, 0.6, 1.1), fpid=C0402)
    chip("C1", 12, "SIG1", "GND")
    b.track("F.Cu", [(12, 14.5), (12, 9)], 0.15, "SIG1")
    b.track("F.Cu", [(12, 15.5), (12, 16.3)], 0.3, "GND")
    b.via(12, 16.3, "GND")
    chip("C2", 30, "SIGA", "SIGB")
    b.track("F.Cu", [(30, 14.5), (30, 9)], 0.15, "SIGA")
    b.track("F.Cu", [(30, 15.5), (30, 21)], 0.15, "SIGB")
    chip("C3", 48, "PWR", "SIGC")
    b.zone("F.Cu", [(44, 9), (52, 9), (52, 14.8), (44, 14.8)], "PWR")
    b.track("F.Cu", [(48, 15.5), (48, 21)], 0.15, "SIGC")
    b.footprint("FID1", "Fiducial", 5, 5, [{"num": 1, "net": "", "dx": 0, "dy": 0, "sx": 1.0, "shape": "circle", "mask": 0.5}])
    return b.save()


def reflow_board(out, name="reflow"):
    """A 3 g part (U1, 5 mm tall) on the left and 0402 resistors next to it (R1) and far from it (R2)."""
    b = Board(_p(out, name), layers=4, size=(80, 40))
    b.plane("In1.Cu", "GND")
    b.plane("In2.Cu", "GND")
    b.footprint("U1", "HEAVY", 14, 20, [{"num": i + 1, "net": "GND", "dx": -4 + 8 * (i % 2), "dy": -4 + 8 * (i // 2), "sx": 1.5, "mask": 0.0}
                                       for i in range(4)], fields={"Mass": "3", "Height": "5"}, courtyard=(-6, -6, 6, 6))
    b.two_pin("R1", "10k", 23, 20, "N1", "N2", fpid=R0402, mask=0.0)
    b.two_pin("R2", "10k", 70, 20, "N3", "N4", fpid=R0402, mask=0.0)
    return b.save()


def parts_board(out, name="parts", planted=True):
    """Parts with a planted defect each, or none: a paste aperture too narrow for the stencil; an 0402
    footprint whose pads are 1.4 mm apart; two 0402 parts 0.1 mm apart; an IC with a 0.3 mm pin pitch;
    an 0201 and a part on the back (for the economic service); no fiducial."""
    b = Board(_p(out, name), layers=4, size=(60, 30))
    b.plane("In1.Cu", "GND")
    b.plane("In2.Cu", "GND")
    p = planted
    b.footprint("U1", "IC", 8, 8, [{"num": i + 1, "net": f"P{i}", "dx": i * 1.5, "dy": 0, "sx": 0.14 if p else 0.3, "sy": 0.5, "mask": 0.0}
                                    for i in range(2)])
    b.two_pin("C1", "100n", 20, 8, "A1", "A2", pitch=1.9 if p else 1.0, pad=(0.5, 0.6), fpid=C0402, mask=0.0)
    b.two_pin("C2", "100n", 30, 8, "B1", "B2", fpid=C0402, mask=0.0)
    b.two_pin("C3", "100n", 31.7 if p else 32.5, 8, "B3", "B4", fpid=C0402, mask=0.0)
    pitch = 0.3 if p else 0.5
    b.footprint("U2", "IC8", 42, 8, [{"num": i + 1, "net": f"Q{i}", "dx": (i % 4) * pitch, "dy": 2.0 * (i // 4), "sx": 0.17, "sy": 0.8, "mask": 0.0}
                                      for i in range(8)], fpid="Package_SO:TSSOP-8_3x3mm")
    if p:
        b.two_pin("C4", "100n", 50, 20, "D1", "D2", pitch=0.5, pad=(0.28, 0.3), fpid="Capacitor_SMD:C_0201_0603Metric", mask=0.0)
        b.footprint("C5", "100n", 20, 20, [{"num": 1, "net": "E1", "dx": -0.5, "dy": 0, "sx": 0.6, "mask": 0.0},
                                           {"num": 2, "net": "E2", "dx": 0.5, "dy": 0, "sx": 0.6, "mask": 0.0}], side="B", fpid=C0402)
    else:
        b.footprint("FID1", "Fiducial", 5, 25, [{"num": 1, "net": "", "dx": 0, "dy": 0, "sx": 1.0, "shape": "circle", "mask": 0.5}])
    return b.save()


def custom_pad_board(out, name="custom"):
    """The pilot's U10 land pattern in miniature: a custom pad (0.2 mm anchor, 1.72 x 1.32 mm copper), a
    paste-only custom aperture (1.2 x 0.75 mm), a rectangular pad with a negative paste margin. Positive
    controls: a custom pad whose copper really is 0.2 x 0.2 mm (U2), and a paste-only custom aperture of
    0.2 x 0.2 mm (U3)."""
    b = Board(_p(out, name), layers=4, size=(40, 20))
    b.plane("In1.Cu", "GND")
    b.plane("In2.Cu", "GND")
    b.footprint("U1", "QFN", 10, 10, [
        {"num": 3, "net": "GND", "dx": 0, "dy": 0, "sx": 0.2, "sy": 0.2, "shape": "custom", "prims": [(-0.86, -0.66, 0.86, 0.66)], "mask": 0.0},
        {"num": "", "dx": 0, "dy": 2.5, "sx": 0.2, "sy": 0.2, "shape": "custom", "only": "paste", "prims": [(-0.6, -0.375, 0.6, 0.375)]},
        {"num": 1, "net": "A", "dx": -3, "dy": 0, "sx": 0.8, "sy": 0.4, "mask": 0.0, "paste_margin": -0.05}],
        fpid="Package_DFN_QFN:QFN-9_test", courtyard=(-4, -1.5, 1.5, 3.5))
    b.footprint("U2", "QFN", 25, 10, [
        {"num": 1, "net": "B", "dx": 0, "dy": 0, "sx": 0.2, "sy": 0.2, "shape": "custom", "prims": [(-0.1, -0.1, 0.1, 0.1)], "mask": 0.0}],
        fpid="Package_DFN_QFN:QFN-1_test", courtyard=(-1, -1, 1, 1))
    b.footprint("U3", "QFN", 32, 10, [
        {"num": 1, "net": "C", "dx": 0, "dy": 0, "sx": 1.0, "sy": 1.0, "mask": 0.0},
        {"num": "", "dx": 0, "dy": 0, "sx": 0.2, "sy": 0.2, "shape": "custom", "only": "paste", "prims": [(-0.1, -0.1, 0.1, 0.1)]}],
        fpid="Package_DFN_QFN:QFN-1b_test", courtyard=(-1, -1, 1, 1))
    return b.save()


def order_board(out, name="order", layers=6, ic_pitch=None, chips=("0402",), via_in_pad=True):
    """For the order options (plan phase 4): 0402-class chips (`chips` gives their EIA sizes), an IC with
    `ic_pitch` mm pin pitch if given, and a via in one chip pad when `via_in_pad`."""
    b = Board(_p(out, name), layers=layers, size=(40, 20))
    b.plane("In1.Cu", "GND")
    b.plane("In2.Cu", "GND")
    dims = {"01005": (0.4, 0.2, 0.2), "0201": (0.6, 0.3, 0.3), "0402": (1.0, 0.5, 0.5), "0603": (1.6, 0.8, 0.8)}
    lib = {"01005": "Capacitor_SMD:C_01005_0402Metric", "0201": "Capacitor_SMD:C_0201_0603Metric",
           "0402": C0402, "0603": "Capacitor_SMD:C_0603_1608Metric"}
    for k, code in enumerate(chips):
        L_, W_, _ = dims[code]
        b.two_pin(f"C{k + 1}", "100n", 6 + 5 * k, 6, f"N{k}", "GND", pitch=L_, pad=(L_ / 2, W_), fpid=lib[code], mask=0.0)
    if via_in_pad:
        L_ = dims[chips[0]][0]
        b.via(6 - L_ / 2, 6, "N0", drill=0.3, size=0.6)   # the pilot's vias in pads: 0.3 mm drill, 0.6 mm pad
    if ic_pitch:
        b.footprint("U1", "IC", 25, 12, [{"num": i + 1, "net": f"Q{i}", "dx": (i % 8) * ic_pitch, "dy": 2.0 * (i // 8),
                                           "sx": ic_pitch / 2, "sy": 0.8, "mask": 0.0} for i in range(16)],
                    fpid="Package_DFN_QFN:QFN-16_test", courtyard=(-0.6, -0.8, 8 * ic_pitch, 2.8))
    b.footprint("FID1", "Fiducial", 3, 17, [{"num": 1, "net": "", "dx": 0, "dy": 0, "sx": 1.0, "shape": "circle", "mask": 0.5}])
    return b.save()


def class_board(out, name="classes"):
    """Plan phase 5. C1: an 0402 whose pads are 1.05 mm apart (beyond the 1.0 mm body); C2: 0.95 mm (fine).
    U1: a part named 'XYZ_123' from KiCad's Package_DFN_QFN library; U2: 'MYPART' with 8 pads on two sides
    (classified by its pads); U3: a QFN-like part whose body (F.Fab, 6 x 6 mm) overhangs its pads by
    1 mm, with C3 0.6 mm from the body and 1.6 mm from the pads."""
    b = Board(_p(out, name), layers=4, size=(60, 30))
    b.plane("In1.Cu", "GND")
    b.plane("In2.Cu", "GND")
    b.two_pin("C1", "100n", 6, 6, "A1", "A2", pitch=1.05 + 0.5, pad=(0.5, 0.5), fpid=C0402, mask=0.0)
    b.two_pin("C2", "100n", 12, 6, "B1", "B2", pitch=0.95 + 0.5, pad=(0.5, 0.5), fpid=C0402, mask=0.0)
    b.footprint("U1", "IC", 22, 6, [{"num": i + 1, "net": f"P{i}", "dx": (i % 4) * 0.5, "dy": 2.0 * (i // 4), "sx": 0.25, "sy": 0.6, "mask": 0.0}
                                     for i in range(8)], fpid="Package_DFN_QFN:XYZ_123", courtyard=(-0.5, -0.6, 2.0, 2.6))
    b.footprint("U2", "IC", 32, 6, [{"num": i + 1, "net": f"Q{i}", "dx": (i % 4) * 0.65, "dy": 3.0 * (i // 4), "sx": 0.35, "sy": 1.0, "mask": 0.0}
                                     for i in range(8)], fpid="MyLib:MYPART", courtyard=(-0.5, -0.8, 2.5, 3.8))
    pads = [{"num": k + 1, "net": f"R{k}", "dx": -2 + (k % 5), "dy": -2 + 4 * (k // 5), "sx": 0.3, "sy": 0.6, "mask": 0.0} for k in range(10)]
    b.footprint("U3", "IC", 45, 15, pads, fpid="Package_DFN_QFN:QFN-10_6x6mm", fab=(-3, -3, 3, 3), courtyard=(-3.2, -3.2, 3.2, 3.2))
    b.two_pin("C3", "100n", 45 + 3 + 0.6 + 0.75, 15, "S1", "S2", pitch=1.0, pad=(0.5, 0.5), fpid=C0402, mask=0.0)
    return b.save()
