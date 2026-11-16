"""Interface standards (#6): per-standard design limits, mapped to nets by name.

The values are typical layout guidance for each interface (impedance and tolerance, data rate and
edge, intra-pair skew, uncoupled length), not the text of the specifications, which several of these
are not free to read. Each entry names where its figures come from. A reviewed config and the board's
.kicad_dru override them; without either they replace openpcbsi's generic defaults.
"""
import re

# name: (net regex, fields). Order matters: the first match wins.
STANDARDS = [
    ("pcie_gen3", r"(?:^|[/_.])PCIE|(?:^|[/_.])PE[TR][PN]?\d|PCI_?E", dict(
        z_diff=85.0, z_tol=12.75, rise_ps=50.0, nyquist_ghz=4.0, uncoupled_max_mm=1.0, skew_max_mm=0.13,
        source="PCI Express base/CEM layout guidance: 85 ohm +/- 15 %, 8 GT/s, intra-pair skew <= 5 mil")),
    ("usb3", r"USB3|(?:^|[/_.])SS(?:TX|RX)|USB_?SS|U3[TR]X", dict(
        z_diff=90.0, z_tol=9.0, rise_ps=50.0, nyquist_ghz=2.5, uncoupled_max_mm=1.0, skew_max_mm=0.13,
        source="USB 3.x layout guidance: 90 ohm +/- 10 %, 5 Gb/s (Gen1)")),
    ("usb2", r"USB|(?:^|[/_.])D[_.]?[PN]$|(?:^|[/_.])DP$|(?:^|[/_.])DM$", dict(
        z_diff=90.0, z_tol=13.5, rise_ps=500.0, nyquist_ghz=0.24, uncoupled_max_mm=5.0, skew_max_mm=1.25,
        source="USB 2.0 high-speed layout guidance: 90 ohm +/- 15 %, 480 Mb/s")),
    ("hdmi_tmds", r"HDMI|TMDS", dict(
        z_diff=100.0, z_tol=15.0, rise_ps=75.0, nyquist_ghz=3.0, uncoupled_max_mm=1.5, skew_max_mm=0.25,
        source="HDMI 2.0 TMDS layout guidance: 100 ohm +/- 15 %, 6 Gb/s")),
    ("mipi_dphy", r"DPHY|CSI|DSI|MIPI", dict(
        z_diff=100.0, z_tol=10.0, rise_ps=100.0, nyquist_ghz=1.25, uncoupled_max_mm=2.0, skew_max_mm=0.5,
        source="MIPI D-PHY layout guidance: 100 ohm +/- 10 %, 2.5 Gb/s")),
    ("sata", r"SATA", dict(
        z_diff=100.0, z_tol=15.0, rise_ps=60.0, nyquist_ghz=3.0, uncoupled_max_mm=1.0, skew_max_mm=0.13,
        source="SATA 3 layout guidance: 85-115 ohm, 6 Gb/s")),
    ("ethernet_mdi", r"(?:^|[/_.])(?:MDI|TRD|TR?X?\d?_?[PN]?$)|ETH|ENET|RGMII", dict(
        z_diff=100.0, z_tol=10.0, rise_ps=2000.0, nyquist_ghz=0.0625, uncoupled_max_mm=10.0, skew_max_mm=2.5,
        source="1000BASE-T MDI layout guidance: 100 ohm +/- 10 %, 125 MBd")),
    ("serdes_mgt", r"(?:^|[/_.])(?:GT[HXYP]?|MGT)|REFCLK|GBTCLK|(?:^|[/_.])DP\d", dict(
        z_diff=100.0, z_tol=10.0, rise_ps=35.0, nyquist_ghz=6.25, uncoupled_max_mm=1.0, skew_max_mm=0.13,
        source="FPGA multi-gigabit transceiver layout guidance (e.g. AMD UG576/UG578): 100 ohm +/- 10 %, 12.5 Gb/s")),
    ("lvds", r"LVDS|(?:^|[/_.])LA\d|(?:^|[/_.])HA\d|(?:^|[/_.])HB\d", dict(
        z_diff=100.0, z_tol=10.0, rise_ps=300.0, nyquist_ghz=0.625, uncoupled_max_mm=5.0, skew_max_mm=2.5,
        source="TIA/EIA-644 LVDS layout guidance: 100 ohm +/- 10 %, ~1.25 Gb/s")),
]
_RX = [(n, re.compile(p, re.I), d) for n, p, d in STANDARDS]

# openpcbsi rule parameter <- standard field
RULE_FIELDS = {("pair_coupling", "limit_mm"): "uncoupled_max_mm", ("timing", "skew_mm"): "skew_max_mm"}


def of_net(net):
    """(name, fields) of the first standard whose pattern matches the net name, or None."""
    base = net.split("/")[-1]
    for n, rx, d in _RX:
        if rx.search(base) or rx.search(net):
            return n, d
    return None


def limit(rule, param, net):
    f = RULE_FIELDS.get((rule, param))
    s = of_net(net)
    if f and s and f in s[1]:
        return s[1][f], f"standard {s[0]} ({s[1]['source']})"
    return None
