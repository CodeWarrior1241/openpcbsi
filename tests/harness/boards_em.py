"""Fixture boards for the 3D-region (#9-#20) and pad-void tests in tests/acceptance_em.py. Each function
builds one board (and its project files) and returns its path. Coordinates in mm.

All are the same AC-coupled lane on a 4-layer board (GND on In1 and In2): a pair LANE0_RX_P/N on F.Cu that
fans out at 45 degrees into two series capacitors C1 / C2 at x = 20 (pads at y = 9 and 11), then
LANE0_RXC_P/N on to vias at x = 34 and on B.Cu. The planner
(`evaluate --em plan`) puts a series-cap box of about [16.75, 6.25, 23.25, 13.75] around C1 / C2 and a
layer-change box around the vias; the boards differ in what else is near those boxes."""
import os

from gen import Board


def _p(out, name):
    d = os.path.join(out, name)
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, name + ".kicad_pcb")


def _lane(b, prefix=""):
    """The AC-coupled lane; `prefix` is put before every lane net name (a hierarchical sheet path)."""
    for sfx, y, yc in (("_P", 9.865, 9.0), ("_N", 10.135, 11.0)):     # 45-degree jogs into and out of the pads
        b.track("F.Cu", [(5, y), (19.5 - abs(y - yc), y), (19.5, yc)], 0.12, prefix + "LANE0_RX" + sfx)
    b.two_pin("C1", "100n", 20, 9.0, prefix + "LANE0_RX_P", prefix + "LANE0_RXC_P", pitch=1.0, pad=(0.5, 0.5))
    b.two_pin("C2", "100n", 20, 11.0, prefix + "LANE0_RX_N", prefix + "LANE0_RXC_N", pitch=1.0, pad=(0.5, 0.5))
    for sfx, y, yc in (("_P", 9.865, 9.0), ("_N", 10.135, 11.0)):
        dy = -0.4 if sfx == "_P" else 0.4
        b.track("F.Cu", [(20.5, yc), (20.5 + abs(y - yc), y), (34 - abs(dy), y), (34, y + dy)], 0.12, prefix + "LANE0_RXC" + sfx)
        b.via(34, y + dy, prefix + "LANE0_RXC" + sfx)
        b.track("B.Cu", [(34, y + dy), (45, y)], 0.12, prefix + "LANE0_RXC" + sfx)


def _base(out, name):
    b = Board(_p(out, name), layers=4, size=(50, 20))
    b.plane("In1.Cu", "GND")
    b.plane("In2.Cu", "GND")
    return b


def hier_lane(out, name="hier", prefix="/serdes/"):
    """The lane with hierarchical net names ("/serdes/LANE0_RX_P"), or any other prefix (a traversal)."""
    b = _base(out, name)
    _lane(b, prefix)
    return b.save()


def diagonal_lane(out, name="diag"):
    """The lane with a 45-degree track of another net (SLOW_A, on B.Cu) crossing the series-cap box's
    xmin-ymax corner cleanly: it enters through the xmin side and leaves through the ymax side, with no
    end inside the box."""
    b = _base(out, name)
    _lane(b)
    b.track("B.Cu", [(14.0, 8.0), (22.0, 16.0)], 0.15, "SLOW_A")
    return b.save()


def stub_lane(out, name="stub", end_y=12.75):
    """The lane with a track of another net (SLOW_B, on F.Cu) coming down from the board's top edge
    region and ENDING at (18.5, end_y), inside the planned series-cap box near its ymax side: the planned
    box is refused; a box whose ymax side is moved out passes."""
    b = _base(out, name)
    _lane(b)
    b.track("F.Cu", [(18.5, 19.0), (18.5, end_y)], 0.15, "SLOW_B")
    return b.save()


def parallel_lane(out, name="parallel", y=12.6):
    """The lane with a track of another net (SLOW_C, on F.Cu) running along x over the whole board at `y`,
    parallel to the series-cap box's ymax side."""
    b = _base(out, name)
    _lane(b)
    b.track("F.Cu", [(1.0, y), (49.0, y)], 0.15, "SLOW_C")
    return b.save()


def fenced_lane(out, name="fenced", ys=(12.6, 14.1, 15.6, 17.1, 18.6)):
    """The lane with a bus of tracks of other nets running along x 1.5 mm apart above the capacitors
    (y = 12.6 ... 18.6): wherever the box's ymax side goes (out to +3 mm, or in towards the parts), one of
    them runs along it within the pre-flight's clearance, so no box passes."""
    b = _base(out, name)
    _lane(b)
    for i, y in enumerate(ys):
        b.track("F.Cu", [(1.0, y), (49.0, y)], 0.15, f"BUS{i}")
    return b.save()


def pad_void(out, name="padvoid", c2_margin=0.15, slot_x=10.0):
    """A lane whose tracks enter the capacitor pads: LANE0_RX_P along y = 9 into C1 (19.5, 9) and
    LANE0_RX_N along y = 11 into C2 (19.5, 11), 0.12 mm tracks, 0.5 mm pads, on F.Cu over GND on In1.
    In1 is voided under every capacitor pad by the margin pad_plane recommends (pad + 0.15 mm), except
    C2's pads, voided by `c2_margin` (larger: more void than asked for). A separate 0.5 mm slot in In1
    across the P track only at x = `slot_x` (a void that is not under any pad)."""
    b = _base(out, name)
    for sfx, y, ref in (("_P", 9.0, "C1"), ("_N", 11.0, "C2")):
        b.track("F.Cu", [(5, y), (19.5, y)], 0.12, "LANE0_RX" + sfx)
        b.two_pin(ref, "100n", 20, y, "LANE0_RX" + sfx, "LANE0_RXC" + sfx, pitch=1.0, pad=(0.5, 0.5))
        b.track("F.Cu", [(20.5, y), (35, y)], 0.12, "LANE0_RXC" + sfx)
        m = 0.15 if ref == "C1" else c2_margin
        for x in (19.5, 20.5):
            b.slot("In1.Cu", x - 0.25 - m, y - 0.25 - m, x + 0.25 + m, y + 0.25 + m)
    if slot_x:
        b.slot("In1.Cu", slot_x - 0.25, 8.4, slot_x + 0.25, 9.6)
    return b.save()


def broken_lane(out, name="broken", gap=(22.0, 22.3)):
    """The lane with LANE0_RXC_P's F.Cu track interrupted between x = gap[0] and gap[1] (a break the
    connection check must find)."""
    b = _base(out, name)
    for sfx, y, yc in (("_P", 9.865, 9.0), ("_N", 10.135, 11.0)):
        b.track("F.Cu", [(5, y), (19.5 - abs(y - yc), y), (19.5, yc)], 0.12, "LANE0_RX" + sfx)
    b.two_pin("C1", "100n", 20, 9.0, "LANE0_RX_P", "LANE0_RXC_P", pitch=1.0, pad=(0.5, 0.5))
    b.two_pin("C2", "100n", 20, 11.0, "LANE0_RX_N", "LANE0_RXC_N", pitch=1.0, pad=(0.5, 0.5))
    for sfx, y, yc in (("_P", 9.865, 9.0), ("_N", 10.135, 11.0)):
        if sfx == "_P":
            b.track("F.Cu", [(20.5, yc), (20.5 + abs(y - yc), y), (gap[0], y)], 0.12, "LANE0_RXC_P")
            b.track("F.Cu", [(gap[1], y), (34, y)], 0.12, "LANE0_RXC_P")
        else:
            b.track("F.Cu", [(20.5, yc), (20.5 + abs(y - yc), y), (34, y)], 0.12, "LANE0_RXC_N")
    return b.save()


def diagonal_pair(out, name="diagpair"):
    """A pair DIAG_P/N on F.Cu running at 45 degrees across the board, P along y = x - 8 and N 0.27 mm
    beside it (0.12 mm tracks, 0.15 mm gap)."""
    b = _base(out, name)
    o = 0.27 / 2 ** 0.5
    b.track("F.Cu", [(9.0, 1.0), (27.0, 19.0)], 0.12, "DIAG_P")
    b.track("F.Cu", [(9.0 + o, 1.0 - o), (27.0 + o, 19.0 - o)], 0.12, "DIAG_N")
    return b.save()


def pad_end(out, name="padend"):
    """A pair END_P/N on F.Cu from the board's left edge to the two pads of a connector J1 at x = 20 (the
    pair's own end: a one-port region in a board-only model)."""
    b = _base(out, name)
    for sfx, y in (("_P", 9.865), ("_N", 10.135)):
        b.track("F.Cu", [(3, y), (20, y)], 0.12, "END" + sfx)
    b.footprint("J1", "CONN", 20, 10, [{"num": 1, "net": "END_P", "dx": 0, "dy": -0.135, "sx": 0.2, "sy": 0.12},
                                       {"num": 2, "net": "END_N", "dx": 0, "dy": 0.135, "sx": 0.2, "sy": 0.12}])
    return b.save()
