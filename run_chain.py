"""Long-chain test: does the electron stay near its launch site?

Run:  python run_chain.py            (about 40 minutes on two cores)
      python run_chain.py --quick    (a few minutes, fewer trajectories)
Output: results/chain_*.csv and results/chain.json
"""
import csv
import json
import os
import sys
from multiprocessing import Pool

import numpy as np
from scipy.special import jv

from model.chain import HBAR, ChainParams, SectorCache, no_emission, trajectory

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
QUICK = "--quick" in sys.argv
WINDOW = 2                      # the fixed window is |site - launch site| <= 2


def write_csv(name, header, columns):
    with open(os.path.join(OUT, name), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for row in zip(*columns):
            w.writerow([f"{x:.8g}" for x in row])


def _batch(args):
    params, times, seeds = args
    p = ChainParams(**params)
    c, x = p.N // 2, np.arange(p.N) - p.N // 2
    cache = SectorCache(p)
    T = len(times)
    s = {k: np.zeros(T) for k in ("win", "win2", "x2", "x2sq", "between", "photons")}
    prob, dead_sum = np.zeros((T, p.N)), np.zeros((T, p.N))
    for seed in seeds:
        pr, dead = trajectory(p, c, times, np.random.default_rng(seed), cache)
        w = pr[:, c - WINDOW:c + WINDOW + 1].sum(1)
        x2 = (pr * x ** 2).sum(1)
        both = np.zeros_like(dead)
        both[:, 1:-1] = dead[:, :-2] & dead[:, 2:]          # both neighboring atoms have emitted
        s["win"] += w; s["win2"] += w ** 2; s["x2"] += x2; s["x2sq"] += x2 ** 2
        s["between"] += (pr * both).sum(1); s["photons"] += dead.sum(1)
        prob += pr; dead_sum += dead
    return s, prob, dead_sum, len(seeds)


def ensemble(params, times, n, seed0=0, procs=2):
    chunks = np.array_split(np.arange(seed0, seed0 + n), max(procs * 4, 1))
    with Pool(procs) as pool:
        parts = pool.map(_batch, [(params, times, list(ch)) for ch in chunks if len(ch)])
    tot = {k: sum(p[0][k] for p in parts) / n for k in parts[0][0]}
    prob = sum(p[1] for p in parts) / n
    dead = sum(p[2] for p in parts) / n
    return {"window": tot["win"], "window_err": np.sqrt(np.maximum(tot["win2"] - tot["win"] ** 2, 0) / n),
            "rms": np.sqrt(tot["x2"]), "x2": tot["x2"], "x2_err": np.sqrt(np.maximum(tot["x2sq"] - tot["x2"] ** 2, 0) / n),
            "between": tot["between"], "photons": tot["photons"], "prob": prob, "dead": dead}


def free_infinite(A, times, R=WINDOW):
    z = 2 * A * times / HBAR
    return sum(jv(n, z) ** 2 for n in range(-R, R + 1)), np.sqrt(2) * A * times / HBAR


def main():
    os.makedirs(OUT, exist_ok=True)
    res = {"window_half_width_sites": WINDOW, "max_2p": 2}
    base = dict(A=0.1, g=1.6, delta=4.37, gamma=0.625, max_2p=2)
    n_main, n_side = (40, 24) if QUICK else (400, 160)

    # ---- 1. main run: slow hopping, 31 sites -------------------------------------------------
    times = np.linspace(0, 100, 41)
    N = 31
    full = ensemble(dict(base, N=N), times, n_main)
    win_free, rms_free = free_infinite(base["A"], times)
    No = 61                                                   # longer chain: no reflection before 100 ns
    off = no_emission(ChainParams(**dict(base, N=No)), No // 2, times)
    xo = np.arange(No) - No // 2
    win_off = off[:, No // 2 - WINDOW:No // 2 + WINDOW + 1].sum(1)
    rms_off = np.sqrt((off * xo ** 2).sum(1))
    write_csv("chain_main.csv",
              ["time_ns", "window_full", "window_full_err", "rms_full", "between_two_1s_atoms", "photons_emitted",
               "window_emission_off", "rms_emission_off", "window_free", "rms_free"],
              [times, full["window"], full["window_err"], full["rms"], full["between"], full["photons"],
               win_off, rms_off, win_free, rms_free])
    x = np.arange(N) - N // 2
    sel = slice(N // 2 - 10, N // 2 + 11)
    k50 = int(np.argmin(abs(times - 50)))
    write_csv("chain_profile_100ns.csv", ["site_minus_launch", "full_model", "emission_off", "free_hopping", "atom_in_1s_full"],
              [x[sel], full["prob"][-1][sel], off[-1][No // 2 - 10:No // 2 + 11], jv(x[sel], 2 * base["A"] * 100 / HBAR) ** 2, full["dead"][-1][sel]])
    res["main"] = {"N": N, "trajectories": n_main, "A": base["A"],
                   "window_100ns": {"full": float(full["window"][-1]), "full_err": float(full["window_err"][-1]),
                                    "emission_off": float(win_off[-1]), "free": float(win_free[-1])},
                   "window_50ns": {"full": float(full["window"][k50]), "emission_off": float(win_off[k50]), "free": float(win_free[k50])},
                   "rms_100ns": {"full": float(full["rms"][-1]), "emission_off": float(rms_off[-1]), "free": float(rms_free[-1])},
                   "between_two_1s_atoms_100ns": float(full["between"][-1]),
                   "photons_emitted_100ns": float(full["photons"][-1])}
    print(json.dumps(res["main"], indent=1), flush=True)

    # ---- 2. the result does not depend on the chain length ----------------------------------
    res["chain_length"] = {}
    for n_sites in (11, 21):
        e = ensemble(dict(base, N=n_sites), times, n_main, seed0=10_000)
        res["chain_length"][str(n_sites)] = {"window_100ns": float(e["window"][-1]), "err": float(e["window_err"][-1]), "rms_100ns": float(e["rms"][-1])}
    res["chain_length"][str(N)] = {"window_100ns": float(full["window"][-1]), "err": float(full["window_err"][-1]), "rms_100ns": float(full["rms"][-1])}
    print(json.dumps(res["chain_length"], indent=1), flush=True)

    # ---- 3. the truncation to two simultaneous 2p atoms is converged ------------------------
    t_short = np.linspace(0, 60, 7)
    a = ensemble(dict(base, N=11, max_2p=2), t_short, 40, seed0=20_000)
    b = ensemble(dict(base, N=11, max_2p=4), t_short, 40, seed0=20_000)
    res["truncation_check_max_difference"] = float(np.abs(a["prob"] - b["prob"]).max())

    # ---- 4. scan of the hopping energy ----------------------------------------------------------
    from model import atom
    gq = atom.quench_rate_exact(base["g"], base["delta"], base["gamma"])
    scan = {"A": [], "ratio": [], "window": [], "window_err": [], "rms": [], "between": [], "photons": [], "window_free": []}
    for A in (0.05, 0.1, 0.2, 0.4):
        e = full if A == base["A"] else ensemble(dict(base, N=N, A=A), times, n_side, seed0=30_000)
        scan["A"].append(A); scan["ratio"].append(gq * HBAR / A)
        scan["window"].append(float(e["window"][-1])); scan["window_err"].append(float(e["window_err"][-1]))
        scan["rms"].append(float(e["rms"][-1])); scan["between"].append(float(e["between"][-1])); scan["photons"].append(float(e["photons"][-1]))
        scan["window_free"].append(float(free_infinite(A, times)[0][-1]))
        print("scan", A, scan["window"][-1], scan["rms"][-1], scan["between"][-1], flush=True)
    write_csv("chain_scan.csv", ["A_micro_eV", "quench_rate_over_hopping_rate", "window_100ns", "window_err", "rms_100ns", "between_two_1s_atoms", "photons_emitted", "window_free_100ns"],
              [scan[k] for k in ("A", "ratio", "window", "window_err", "rms", "between", "photons", "window_free")])
    res["scan"] = scan
    # energy of the electron next to one excited atom (level repulsion)
    res["energy_shift_per_excited_neighbor"] = float((-base["delta"] + np.sqrt(base["delta"] ** 2 + 4 * base["g"] ** 2)) / 2)
    with open(os.path.join(OUT, "chain.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("done", flush=True)


if __name__ == "__main__":
    main()
