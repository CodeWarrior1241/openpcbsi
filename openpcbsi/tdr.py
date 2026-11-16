"""Time-domain reflectometry from S11 (#17), as in Backplane's em/m3d.py.

The spectrum is filtered with a Gaussian step of 10-90 % rise time `rise` and zero-padded `pad` times
before the inverse FFT, so a via or capacitor dip is sampled many times rather than once.
"""
import numpy as np


def tdr(f, S11, Z0=100.0, rise=35e-12, pad=16):
    """Impedance seen from the port against time. f: uniform grid from (about) 0 Hz.
    The filter is delayed by 4 sigma and the time axis shifted back, so a reflection right at the port
    reads its full height (an undelayed Gaussian wraps half its impulse to negative time)."""
    f = np.asarray(f, float)
    df = f[1] - f[0]
    n = len(f)
    sig = rise / 2.563                              # Gaussian step 10-90 % = 2.563 sigma
    t0 = 4.0 * sig                                  # delay the filter so its impulse is causal (no wrap at t = 0)
    H = np.exp(-0.5 * (2 * np.pi * f * sig) ** 2) * np.exp(-2j * np.pi * f * t0)
    S = np.array(S11, complex)
    S[0] = np.real(S[1])
    N = (n - 1) * pad + 1
    Sp = np.zeros(N, complex)
    Sp[:n] = S * H
    h = np.fft.irfft(Sp, n=2 * (N - 1))             # impulse response
    step = np.cumsum(h)
    t = np.arange(len(step)) / (2 * (N - 1) * df) - t0
    rho = np.clip(step, -0.95, 0.95)
    return t, Z0 * (1 + rho) / (1 - rho)


def renorm(S11, Zl, Z0=100.0):
    """Reflection renormalised from Z0 to the line impedance Zl (one-port form)."""
    r = (Zl - Z0) / (Zl + Z0)
    S = np.asarray(S11)
    return (S - r) / (1 - r * S)


def tdr_line(f, S11, Zl, Z0=100.0, rise=35e-12):
    """TDR referenced to the model's own line: flat at Zl until the first discontinuity."""
    return tdr(f, renorm(S11, Zl, Z0), Zl, rise)


def excursions(t, z, z_line, t_max=None):
    """Largest dip and peak relative to the line impedance (%), up to t_max."""
    m = np.ones_like(t, bool) if t_max is None else t <= t_max
    zz = z[m]
    return {"min_pct": round(float((zz.min() - z_line) / z_line * 100), 2),
            "max_pct": round(float((zz.max() - z_line) / z_line * 100), 2)}
