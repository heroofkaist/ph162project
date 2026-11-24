"""Isolated hydrogen atom in a static electric field: 2s-2p mixing and the quench rate of 2s.

Units: energy in micro-eV, time in ns unless a name says otherwise.
Sources of the constants: CODATA 2018; Bezginov et al., Science 365, 1007 (2019) for the
Lamb shift; Wiese and Fuhr, J. Phys. Chem. Ref. Data 38, 565 (2009) for the 2p -> 1s rate.
"""
import numpy as np
from scipy.integrate import quad

HBAR_EV_S = 6.582119569e-16
H_EV_S = 4.135667696e-15
BOHR_M = 5.29177210903e-11
E_CHARGE_C = 1.602176634e-19
EPS0 = 8.8541878128e-12
HBAR = HBAR_EV_S * 1e6 * 1e9                 # micro-eV * ns

LAMB_SHIFT_MHZ = 1057.8298                   # 2s1/2 - 2p1/2 interval
DELTA = LAMB_SHIFT_MHZ * 1e6 * H_EV_S * 1e6  # 4.3748 micro-eV
GAMMA_2P_PER_S = 6.2649e8                    # 2p -> 1s decay rate (lifetime 1.596 ns)
GAMMA_2P = GAMMA_2P_PER_S * 1e-9             # 1/ns


def radial_integral_2s_2p():
    """Integral of R20 * R21 * r^3 dr in units of the Bohr radius (exact value -3*sqrt(3))."""
    r20 = lambda r: (1 / np.sqrt(2)) * (1 - r / 2) * np.exp(-r / 2)
    r21 = lambda r: (1 / np.sqrt(24)) * r * np.exp(-r / 2)
    return quad(lambda r: r20(r) * r21(r) * r ** 3, 0, np.inf)[0]


def angular_integral_s_p0():
    """Integral of Y00 * cos(theta) * Y10 over the sphere (exact value 1/sqrt(3))."""
    f = lambda th: (1 / np.sqrt(4 * np.pi)) * np.cos(th) * np.sqrt(3 / (4 * np.pi)) * np.cos(th) * np.sin(th)
    return 2 * np.pi * quad(f, 0, np.pi)[0]


def dipole_length_orbital():
    """<2s| z |2p, m=0> in Bohr radii: -3."""
    return radial_integral_2s_2p() * angular_integral_s_p0()


def dipole_length_fine_structure():
    """|<2s1/2| z |2p1/2>| in Bohr radii: sqrt(3).

    The state |2p1/2, mj=+1/2> contains |2p, m=0, spin up> with Clebsch-Gordan weight -sqrt(1/3).
    """
    return abs(dipole_length_orbital()) * np.sqrt(1 / 3)


def mixing_energy(field_V_per_m):
    """V = sqrt(3) e a0 E in micro-eV."""
    return dipole_length_fine_structure() * BOHR_M * field_V_per_m * 1e6


def quench_rate_weak_field(V, delta=DELTA, gamma=GAMMA_2P):
    """Gamma_eff = Gamma V^2 / (Delta^2 + (hbar Gamma / 2)^2), in 1/ns."""
    return gamma * V ** 2 / (delta ** 2 + (HBAR * gamma / 2) ** 2)


def quench_rate_exact(V, delta=DELTA, gamma=GAMMA_2P):
    """Decay rate of the 2s-like eigenstate of [[0, V], [V, -Delta - i hbar Gamma / 2]], in 1/ns.

    The 2s-like state is the eigenstate with the slower decay.  Its rate never exceeds Gamma / 2.
    """
    m = np.array([[0.0, V], [V, -delta - 0.5j * HBAR * gamma]])
    return float(-2 * np.linalg.eigvals(m).imag.max() / HBAR)


def survival_static(V, t, delta, gamma):
    """Probability that an atom held in a constant mixing V has not emitted by time t (starts in 2s)."""
    m = np.array([[0.0, V], [V, -delta - 0.5j * HBAR * gamma]])
    w, U = np.linalg.eig(m)
    a = U @ np.diag(np.exp(-1j * w * t / HBAR)) @ np.linalg.inv(U) @ np.array([1.0, 0.0])
    return float(np.vdot(a, a).real)


def coulomb_field(distance_m):
    """Field of a bare electron at a given distance, in V/m."""
    return E_CHARGE_C / (4 * np.pi * EPS0 * distance_m ** 2)
