"""Coherence of the electron on the long chain (Q1 of the project question).

The electron starts on the central site, so its coherence C = sum_{i<j} |rho_ij| starts at 0.
C grows when the electron spreads as a coherent wave.

Run:  python run_chain_coherence.py      (about 15 minutes on two cores)
Output: results/chain_coherence.csv and results/chain_coherence.json
"""
import csv
import json
import os
from multiprocessing import Pool

import numpy as np
from scipy.special import jv

from model.chain import HBAR, ChainParams, SectorCache, trajectory

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
BASE = dict(A=0.1, g=1.6, delta=4.37, gamma=0.625, max_2p=2)
TIMES = np.linspace(0, 100, 41)


def coherence(rho):
    iu = np.triu_indices(rho.shape[-1], 1)
    return np.abs(rho[..., iu[0], iu[1]]).sum(axis=-1)


def _batch(args):
    params, seeds = args
    p = ChainParams(**params)
    cache = SectorCache(p)
    s1 = np.zeros((len(TIMES), p.N, p.N), dtype=complex)
    s2 = np.zeros((len(TIMES), p.N, p.N))
    win = np.zeros(len(TIMES))
    c = p.N // 2
    for seed in seeds:
        pr, _, rho = trajectory(p, c, TIMES, np.random.default_rng(seed), cache, return_rho=True)
        s1 += rho
        s2 += np.abs(rho) ** 2
        win += pr[:, c - 2:c + 3].sum(1)
    return s1, s2, win, len(seeds)


def main():
    n = 400
    chunks = np.array_split(np.arange(n), 8)                 # the same seeds as the main run of run_chain.py
    with Pool(2) as pool:
        parts = pool.map(_batch, [(dict(BASE, N=31), list(ch)) for ch in chunks])
    s1 = sum(x[0] for x in parts) / n
    s2 = sum(x[1] for x in parts) / n
    win = sum(x[2] for x in parts) / n
    c_full = coherence(s1)
    # statistical floor: expected sum of |noise| of the off-diagonal elements
    var = np.maximum(s2 - np.abs(s1) ** 2, 0) / n
    iu = np.triu_indices(31, 1)
    floor = (np.sqrt(np.pi) / 2 * np.sqrt(var[:, iu[0], iu[1]])).sum(axis=1)

    # emission off: one deterministic state vector on a longer chain (no reflection before 100 ns)
    q = ChainParams(**dict(BASE, N=61, gamma=0.0))
    _, _, rho_off = trajectory(q, 30, TIMES, np.random.default_rng(0), return_rho=True)
    c_off = coherence(rho_off)

    # free hopping on an infinite chain: C_j = i^j J_j(2At/hbar), a pure state
    x = np.arange(-200, 201)
    amp = np.abs(jv(x[None, :], 2 * BASE["A"] * TIMES[:, None] / HBAR))
    c_free = (amp.sum(1) ** 2 - (amp ** 2).sum(1)) / 2

    with open(os.path.join(OUT, "chain_coherence.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time_ns", "coherence_full", "statistical_floor_full", "coherence_emission_off", "coherence_free", "window_full_check"])
        for row in zip(TIMES, c_full, floor, c_off, c_free, win):
            w.writerow([f"{v:.8g}" for v in row])
    k50 = int(np.argmin(abs(TIMES - 50)))
    res = {"trajectories": n,
           "coherence_100ns": {"full": float(c_full[-1]), "statistical_floor": float(floor[-1]), "emission_off": float(c_off[-1]), "free": float(c_free[-1])},
           "coherence_50ns": {"full": float(c_full[k50]), "emission_off": float(c_off[k50]), "free": float(c_free[k50])},
           "max_coherence_full": float(c_full.max()), "time_of_max_full": float(TIMES[int(np.argmax(c_full))]),
           "window_100ns_check": float(win[-1])}
    with open(os.path.join(OUT, "chain_coherence.json"), "w") as f:
        json.dump(res, f, indent=2)
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
