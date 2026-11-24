"""One electron on a long open chain of N sites with N hydrogen-like atoms (2s, 2p, 1s).

Same physics as model/ring.py, but the chain is long enough to ask whether the electron
stays near its launch site.  The joint space has N * 3**N states, which is too large for a
density matrix when N is larger than about 7.  This module therefore uses quantum
trajectories (Dalibard, Castin and Molmer 1992):

* between photon emissions the state vector evolves with H_eff = H - (i hbar / 2) sum L^dag L;
* its norm falls, and a photon is emitted when the squared norm reaches a random number;
* the emitting atom is then in 1s for the rest of the trajectory and drops out of the dynamics.

An atom in 1s is a classical record, so only the atoms that are still in 2s or 2p are kept in
the state vector.  Optionally the number of atoms that are in 2p at the same time is limited
to `max_2p` (the 2p level decays within 1.6 ns, so many simultaneous 2p atoms are very rare).
`exact_density_matrix` solves the full Lindblad equation for short chains and is used to
check the trajectories.

Units: energy in micro-eV, time in ns, rates in 1/ns.
"""
from dataclasses import dataclass
from itertools import combinations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import csr_matrix, identity

HBAR = 0.6582119569  # micro-eV * ns


@dataclass(frozen=True)
class ChainParams:
    N: int = 15           # number of sites and atoms
    A: float = 0.1        # hopping energy
    g: float = 1.6        # 2s-2p mixing energy for an electron on a neighboring site
    delta: float = 4.37   # 2s-2p gap
    gamma: float = 0.625  # 2p -> 1s decay rate
    max_2p: int = 3       # largest number of atoms that may be in 2p at the same time


def field_sign(site, atom):
    """Sign of the electron's field at an atom: +1 if the electron is on atom+1, -1 on atom-1."""
    if site == atom + 1:
        return 1.0
    if site == atom - 1:
        return -1.0
    return 0.0


class Sector:
    """Hamiltonian for a fixed set of atoms that have not emitted yet (the 'active' atoms)."""

    def __init__(self, p: ChainParams, active):
        self.p, self.active = p, tuple(active)
        n_act = len(self.active)
        kmax = min(p.max_2p, n_act)
        # subsets of active atoms that are in 2p, as bit masks over positions in `active`
        self.subsets = [sum(1 << b for b in c) for k in range(kmax + 1) for c in combinations(range(n_act), k)]
        self.index = {m: k for k, m in enumerate(self.subsets)}
        self.nsub = len(self.subsets)
        self.dim = p.N * self.nsub
        self.count = np.array([bin(m).count("1") for m in self.subsets])
        self.has = np.array([[(m >> b) & 1 for m in self.subsets] for b in range(n_act)], dtype=bool).reshape(n_act, self.nsub)
        pos = {a: b for b, a in enumerate(self.active)}

        rows, cols, vals = [], [], []
        for j in range(p.N):
            for k, m in enumerate(self.subsets):
                col = j * self.nsub + k
                rows.append(col); cols.append(col)
                vals.append(-p.delta * self.count[k] - 0.5j * HBAR * p.gamma * self.count[k])
                for j2 in (j - 1, j + 1):                       # hopping
                    if 0 <= j2 < p.N:
                        rows.append(j2 * self.nsub + k); cols.append(col); vals.append(-p.A)
                for atom in (j - 1, j + 1):                     # the field flips a neighboring atom
                    if atom in pos:
                        k2 = self.index.get(m ^ (1 << pos[atom]))
                        if k2 is not None:
                            rows.append(j * self.nsub + k2); cols.append(col); vals.append(p.g * field_sign(j, atom))
        H = csr_matrix((vals, (rows, cols)), shape=(self.dim, self.dim), dtype=complex)
        shift = -p.delta * kmax / 2                              # real shift: only a global phase
        self.M = ((-1j / HBAR) * (H - shift * identity(self.dim, dtype=complex, format="csr"))).tocsr()
        self.norm = abs(self.M).sum(axis=0).max()

    def propagate(self, psi, tau):
        """exp(M tau) psi by a scaled Taylor series (accurate to about 1e-13)."""
        if tau <= 0:
            return psi
        steps = max(1, int(np.ceil(self.norm * tau / 4.0)))
        h = tau / steps
        for _ in range(steps):
            term, out = psi, psi.copy()
            for n in range(1, 40):
                term = self.M @ term * (h / n)
                out += term
                if np.abs(term).max() < 1e-15:
                    break
            psi = out
        return psi

    def site_probabilities(self, psi):
        w = np.abs(psi.reshape(self.p.N, self.nsub)) ** 2
        return w.sum(axis=1)

    def weight_2p(self, psi):
        """Unnormalized probability that each active atom is in 2p."""
        w = (np.abs(psi.reshape(self.p.N, self.nsub)) ** 2).sum(axis=0)
        return self.has @ w


class SectorCache(dict):
    def __init__(self, p):
        super().__init__()
        self.p = p

    def __missing__(self, active):
        if len(self) > 400:
            self.clear()
        self[active] = Sector(self.p, active)
        return self[active]


def jump(old: Sector, new: Sector, psi, atom):
    """State after atom `atom` has emitted: keep the part where it was in 2p, then remove it."""
    b = old.active.index(atom)
    out = np.zeros((old.p.N, new.nsub), dtype=complex)
    src = psi.reshape(old.p.N, old.nsub)
    low = (1 << b) - 1
    for k, m in enumerate(old.subsets):
        if (m >> b) & 1:
            rest = m ^ (1 << b)
            m2 = (rest & low) | ((rest >> (b + 1)) << b)        # close the gap left by the removed atom
            out[:, new.index[m2]] = src[:, k]
    out = out.ravel()
    return out / np.linalg.norm(out)


def trajectory(p: ChainParams, start, times, rng, cache=None, return_rho=False):
    """One quantum trajectory.  Returns site probabilities [len(times), N] and the 1s pattern [len(times), N].

    With return_rho=True it also returns the density matrix of the electron alone, rho_el [len(times), N, N],
    obtained by tracing out the atoms.  Its average over trajectories is the electron matrix of the Lindblad equation.
    """
    cache = cache if cache is not None else SectorCache(p)
    sector = cache[tuple(range(p.N))]
    psi = np.zeros(sector.dim, dtype=complex)
    psi[start * sector.nsub] = 1.0
    threshold = rng.random()
    probs = np.zeros((len(times), p.N))
    dead = np.zeros((len(times), p.N), dtype=bool)
    rho = np.zeros((len(times), p.N, p.N), dtype=complex) if return_rho else None
    now = 0.0
    for it, t in enumerate(times):
        while now < t - 1e-12:
            trial = sector.propagate(psi, t - now)
            if np.vdot(trial, trial).real > threshold or p.gamma == 0.0:
                psi, now = trial, t
                break
            lo, hi = 0.0, t - now                               # find the emission time by bisection
            for _ in range(14):
                mid = 0.5 * (lo + hi)
                trial = sector.propagate(psi, mid)
                if np.vdot(trial, trial).real > threshold:
                    lo = mid
                else:
                    hi = mid
            psi = sector.propagate(psi, hi)
            now += hi
            w = sector.weight_2p(psi)
            atom = sector.active[rng.choice(len(w), p=w / w.sum())]
            new = cache[tuple(a for a in sector.active if a != atom)]
            psi, sector = jump(sector, new, psi, atom), new
            threshold = rng.random()
        pj = sector.site_probabilities(psi)
        probs[it] = pj / pj.sum()
        dead[it] = [a not in sector.active for a in range(p.N)]
        if return_rho:
            m = psi.reshape(p.N, sector.nsub)
            rho[it] = (m @ m.conj().T) / pj.sum()
    if return_rho:
        return probs, dead, rho
    return probs, dead


def no_emission(p: ChainParams, start, times):
    """Emission switched off (gamma = 0): one deterministic state vector."""
    q = ChainParams(**{**p.__dict__, "gamma": 0.0})
    return trajectory(q, start, times, np.random.default_rng(0))[0]


def free_hopping(N, A, start, times):
    """Mixing switched off: exact single-particle evolution on the open chain."""
    H = -A * (np.eye(N, k=1) + np.eye(N, k=-1))
    w, U = np.linalg.eigh(H)
    amp = U @ (np.exp(-1j * np.outer(w, times) / HBAR) * U[start][:, None])
    return (np.abs(amp) ** 2).T


# ------------------------------------------------------------------ exact check for short chains
def exact_density_matrix(p: ChainParams, start, times):
    """Full Lindblad equation with all N * 3**N states.  Only for N <= 4."""
    N = p.N
    eye = np.eye(3)
    X = np.zeros((3, 3)); X[0, 1] = X[1, 0] = 1.0
    Hat = np.diag([0.0, -p.delta, 0.0])
    Lj = np.zeros((3, 3)); Lj[2, 1] = 1.0

    def emb(el, ops):
        out = el
        for o in ops:
            out = np.kron(out, o)
        return out

    def on(op, i):
        ops = [eye] * N
        ops[i] = op
        return ops

    hop = -p.A * (np.eye(N, k=1) + np.eye(N, k=-1))
    H = emb(hop, [eye] * N).astype(complex)
    for i in range(N):
        H += emb(np.eye(N), on(Hat, i))
        fld = np.diag([field_sign(j, i) for j in range(N)])
        H += p.g * emb(fld, on(X, i))
    Ls = [np.sqrt(p.gamma) * emb(np.eye(N), on(Lj, i)) for i in range(N)]
    Heff = H - 0.5j * HBAR * sum(L.conj().T @ L for L in Ls)
    dim = N * 3 ** N
    psi = np.zeros(dim, dtype=complex); psi[start * 3 ** N] = 1.0
    rho0 = np.outer(psi, psi.conj())

    def rhs(_, y):
        r = y.reshape(dim, dim)
        d = (-1j / HBAR) * (Heff @ r - r @ Heff.conj().T)
        for L in Ls:
            d += L @ r @ L.conj().T
        return d.ravel()

    sol = solve_ivp(rhs, (0, times[-1]), rho0.ravel(), t_eval=times, method="DOP853", rtol=1e-9, atol=1e-11)
    out = np.zeros((len(times), N))
    for k in range(len(times)):
        r = sol.y[:, k].reshape(N, 3 ** N, N, 3 ** N)
        out[k] = np.real(np.einsum("iaia->i", r))
    return out
