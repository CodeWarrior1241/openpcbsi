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
