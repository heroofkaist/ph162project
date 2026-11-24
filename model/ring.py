"""One electron on a three-site ring, coupled to three hydrogen-like atoms.

Each atom has three levels (2s, 2p, 1s).  The electron's position-dependent
electric field mixes 2s and 2p on the neighboring atoms, and 2p decays to 1s
by emitting a photon.  See docs/derivation.md for where every term comes from.

Units: energy in micro-eV, time in ns, rates in 1/ns.
Basis order: |electron site> x |atom 0> x |atom 1> x |atom 2>,
with the atomic levels ordered (2s, 2p, 1s).  The joint space has 3 * 3**3 = 81 states.
"""
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import csr_matrix, identity, kron
from scipy.sparse.linalg import expm_multiply

HBAR = 0.6582119569          # micro-eV * ns  (CODATA 2018)
N_SITES = 3
DIM = N_SITES * 3 ** N_SITES  # 81


@dataclass(frozen=True)
class Params:
    """Parameters of the ring model (illustrative values of the reference run)."""
    A: float = 1.0        # hopping energy (micro-eV)
    g: float = 1.6        # 2s-2p mixing energy when the electron is on a neighboring site (micro-eV)
    delta: float = 4.37   # 2s-2p gap, the Lamb shift rounded to three digits (micro-eV)
    gamma: float = 0.625  # 2p -> 1s decay rate, 1 / (1.6 ns)
    e1s: float = 0.0      # energy of 1s; it cancels in every result (see tests)


def _embed(electron_op, atom_ops):
    out = electron_op
    for op in atom_ops:
        out = np.kron(out, op)
    return out


def _on_atom(op, i):
    ops = [np.eye(3)] * N_SITES
    ops = list(ops)
    ops[i] = op
    return ops


def operators(p: Params = Params()):
    """Build the Hamiltonian, the three jump operators and the observables."""
    eye = np.eye(3)
    T = np.roll(eye, 1, axis=0)                       # T|j> = |j+1>  (indices mod 3)
    n = [np.diag(eye[j]) for j in range(N_SITES)]     # n_j = |j><j|
    X = np.zeros((3, 3)); X[0, 1] = X[1, 0] = 1.0     # X = |2s><2p| + |2p><2s|
    H_atom = np.diag([0.0, -p.delta, p.e1s])
    jump = np.zeros((3, 3)); jump[2, 1] = 1.0         # |1s><2p|

    H = _embed(-p.A * (T + T.T), [eye] * 3).astype(complex)          # hopping
    for i in range(N_SITES):
        H += _embed(eye, _on_atom(H_atom, i))                         # atomic energies
        field = n[(i + 1) % 3] - n[(i - 1) % 3]                       # sign of the field at atom i
        H += p.g * _embed(field, _on_atom(X, i))                      # field mixes 2s and 2p
    L = [np.sqrt(p.gamma) * _embed(eye, _on_atom(jump, i)) for i in range(N_SITES)]

    return {
        "H": H,
        "L": L,
        "current": _embed(1j * (T - T.T), [eye] * 3),                 # (hbar / A) * ring current
        "P1s": [_embed(eye, _on_atom(np.diag([0.0, 0.0, 1.0]), i)) for i in range(N_SITES)],
        "n": [_embed(nj, [eye] * 3) for nj in n],
    }


def initial_state(kind="circulating"):
    """Electron state times all atoms in 2s."""
    if kind == "circulating":
        electron = np.exp(2j * np.pi * np.arange(3) / 3) / np.sqrt(3)
    elif kind == "localized":
        electron = np.array([1.0, 0.0, 0.0], dtype=complex)
    else:
        raise ValueError(kind)
    atoms = np.zeros(27); atoms[0] = 1.0
    psi = np.kron(electron, atoms)
    return np.outer(psi, psi.conj())


def effective_hamiltonian(ops):
    """H_eff = H - (i hbar / 2) sum_i L_i^dagger L_i."""
    return ops["H"] - 0.5j * HBAR * sum(L.conj().T @ L for L in ops["L"])


def evolve(ops, rho0, times, rtol=1e-12, atol=1e-14):
    """Integrate the Lindblad equation with an explicit Runge-Kutta method (DOP853)."""
    H_eff, Ls = effective_hamiltonian(ops), ops["L"]

    def rhs(_, y):
        rho = y.reshape(DIM, DIM)
        d = (-1j / HBAR) * (H_eff @ rho - rho @ H_eff.conj().T)
        for L in Ls:
            d += L @ rho @ L.conj().T
        return d.ravel()

    sol = solve_ivp(rhs, (times[0], times[-1]), rho0.ravel().astype(complex), t_eval=times,
                    method="DOP853", rtol=rtol, atol=atol)
    return [sol.y[:, k].reshape(DIM, DIM) for k in range(len(times))]


def evolve_expm(ops, rho0, t):
    """Independent cross-check: matrix exponential of the 6561 x 6561 Liouvillian."""
    H_eff = csr_matrix(effective_hamiltonian(ops))
    eye = identity(DIM)
    liouv = (-1j / HBAR) * (kron(H_eff, eye) - kron(eye, H_eff.conj()))
    for L in ops["L"]:
        liouv = liouv + kron(csr_matrix(L), csr_matrix(L.conj()))
    return expm_multiply(liouv.tocsr() * t, rho0.ravel().astype(complex)).reshape(DIM, DIM)


# ----------------------------------------------------------------- observables
def reduced_electron(rho):
    """Trace out the atoms: 3 x 3 density matrix of the electron alone."""
    return np.einsum("iaja->ij", rho.reshape(3, 27, 3, 27))


def coherence(rho):
    """Sum of the three off-diagonal magnitudes of the electron matrix (1 at the start)."""
    r = reduced_electron(rho)
    return float(sum(abs(r[i, j]) for i in range(3) for j in range(i + 1, 3)))


def current(ops, rho):
    """Ring current in units of A / hbar (not normalized)."""
    return float(np.real(np.trace(ops["current"] @ rho)))


def ground_fraction(ops, rho):
    """Mean probability that an atom is in 1s."""
    return float(np.mean([np.real(np.trace(P @ rho)) for P in ops["P1s"]]))


def site_probabilities(rho):
    return np.real(np.diag(reduced_electron(rho)))


def record_statistics(rho):
    """Joint probability of (electron site, pattern of 1s atoms) and their mutual information.

    The pattern index is a bit mask: bit i is set if atom i is in 1s.
    Returns (mutual information in bits, joint[3, 8]).
    """
    d = np.real(np.diag(rho)).reshape(3, 3, 3, 3)
    joint = np.zeros((3, 8))
    for a in range(3):
        for b in range(3):
            for c in range(3):
                joint[:, (a == 2) + 2 * (b == 2) + 4 * (c == 2)] += d[:, a, b, c]
    ps, pm = joint.sum(1), joint.sum(0)
    info = sum(joint[j, m] * np.log2(joint[j, m] / (ps[j] * pm[m]))
               for j in range(3) for m in range(8) if joint[j, m] > 1e-15)
    return float(info), joint


def physicality(rho):
    """(|trace - 1|, largest anti-Hermitian element, smallest eigenvalue)."""
    ev = np.linalg.eigvalsh((rho + rho.conj().T) / 2)
    return abs(np.trace(rho) - 1), np.abs(rho - rho.conj().T).max(), ev.min()
