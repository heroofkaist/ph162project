"""Check the quantum trajectories against the exact Lindblad solution on a short chain.

Run:  python validate_chain.py      (about 10 minutes)
Output: results/chain_validation.json
"""
import json
import os

import numpy as np

from model.chain import ChainParams, SectorCache, exact_density_matrix, trajectory

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")


def main():
    p = ChainParams(N=4, A=0.3, max_2p=4)
    times = np.linspace(0, 40, 21)
    exact = exact_density_matrix(p, 1, times)
    n = 1500
    cache = SectorCache(p)
    runs = np.array([trajectory(p, 1, times, np.random.default_rng(100_000 + k), cache)[0] for k in range(n)])
    mean, err = runs.mean(axis=0), runs.std(axis=0) / np.sqrt(n)
    diff = np.abs(mean - exact)
    res = {"N": p.N, "A": p.A, "trajectories": n, "max_abs_difference": float(diff.max()),
           "max_difference_in_standard_errors": float((diff[1:] / err[1:]).max()),
           "rms_difference_in_standard_errors": float(np.sqrt(((diff[1:] / err[1:]) ** 2).mean()))}
    with open(os.path.join(OUT, "chain_validation.json"), "w") as f:
        json.dump(res, f, indent=2)
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
