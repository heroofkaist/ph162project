"""Reproduce every number and every curve of the project.

Run:  python run_all.py
Output: results/results.json and the CSV files in results/.
"""
import csv
import json
import os

import numpy as np

from model import atom
from model.ring import (HBAR, Params, coherence, current, evolve, evolve_expm, ground_fraction,
                        initial_state, operators, physicality, record_statistics, site_probabilities)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
os.makedirs(OUT, exist_ok=True)


def write_csv(name, header, columns):
    with open(os.path.join(OUT, name), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for row in zip(*columns):
            w.writerow([f"{x:.10g}" for x in row])


def time_series(params, kind, times):
    ops = operators(params)
    rhos = evolve(ops, initial_state(kind), times)
    j0 = current(ops, rhos[0]) if kind == "circulating" else 1.0
    return ops, rhos, j0


def main():
    res = {"units": "energy in micro-eV, time in ns, rates in 1/ns", "hbar": HBAR,
           "parameters": Params().__dict__}
    times = np.round(np.arange(0, 12.0001, 0.1), 10)
    settings = {"full": Params(), "emission_off": Params(gamma=0.0), "mixing_off": Params(g=0.0)}

    # ---- 1. circulating electron, three settings, 0 to 12 ns -------------------------------
    worst = [0.0, 0.0, 0.0]
    for name, p in settings.items():
        ops, rhos, j0 = time_series(p, "circulating", times)
        coh = [coherence(r) for r in rhos]
        cur = [current(ops, r) / j0 for r in rhos]
        gnd = [ground_fraction(ops, r) for r in rhos]
        write_csv(f"timeseries_{name}.csv", ["time_ns", "coherence", "current_over_initial", "fraction_1s"],
                  [times, coh, cur, gnd])
        late = times >= 6
        res[name] = {"coherence_12ns": coh[-1], "current_12ns": cur[-1], "fraction_1s_12ns": gnd[-1],
                     "site_probabilities_12ns": site_probabilities(rhos[-1]).tolist(),
                     "current_range_6_to_12ns": [float(np.min(np.array(cur)[late])), float(np.max(np.array(cur)[late]))]}
        for r in rhos:
            t, h, e = physicality(r)
            worst = [max(worst[0], float(t)), max(worst[1], float(h)), min(worst[2], float(e))]
        if name == "full":
            info, joint = record_statistics(rhos[-1])
            infos = [record_statistics(r)[0] for r in rhos]
            res["record"] = {"mutual_information_bits_12ns": info, "max_mutual_information_bits": max(infos),
                             "perfect_record_bits": float(np.log2(3)),
                             "prob_at_least_one_emission": float(1 - joint[:, 0].sum()),
                             "prob_only_atom0_decayed": float(joint[:, 1].sum()),
                             "electron_site_given_only_atom0": (joint[:, 1] / joint[:, 1].sum()).tolist()}
            write_csv("record_information.csv", ["time_ns", "mutual_information_bits"], [times, infos])
            # cross-checks
            rho_expm = evolve_expm(ops, rhos[0], 12.0)
            res["checks"] = {"two_solvers_max_abs_difference": float(np.abs(rho_expm - rhos[-1]).max())}
            ops2, rhos2, j02 = time_series(Params(e1s=-37.0), "circulating", np.array([0.0, 12.0]))
            res["checks"]["effect_of_1s_energy_offset"] = float(abs(current(ops2, rhos2[-1]) / j02 - cur[-1]))
            for nm, pp in (("full", Params(delta=atom.DELTA)), ("emission_off", Params(delta=atom.DELTA, gamma=0.0))):
                o3, r3, j3 = time_series(pp, "circulating", np.array([0.0, 12.0]))
                res["checks"][f"current_12ns_with_exact_gap_{nm}"] = current(o3, r3[-1]) / j3
    res["checks"].update({"max_trace_error": worst[0], "max_hermiticity_error": worst[1], "min_eigenvalue": worst[2]})

    # ---- 2. long time: the effect switches itself off --------------------------------------
    tl = np.round(np.arange(0, 100.0001, 0.5), 10)
    ops, rhos, j0 = time_series(Params(), "circulating", tl)
    cur = [current(ops, r) / j0 for r in rhos]
    gnd = [ground_fraction(ops, r) for r in rhos]
    write_csv("long_time.csv", ["time_ns", "current_over_initial", "fraction_1s"], [tl, cur, gnd])
    ops, rhos, j0 = time_series(Params(), "circulating", np.array([0.0, 60.0, 300.0]))
    res["long_time"] = {"fraction_1s_60ns": ground_fraction(ops, rhos[1]), "current_60ns": current(ops, rhos[1]) / j0,
                        "fraction_1s_300ns": ground_fraction(ops, rhos[2]), "current_300ns": current(ops, rhos[2]) / j0}

    # ---- 3. electron launched on site 0 ------------------------------------------------------
    loc = {}
    for name, p in (("full", Params()), ("free", Params(g=0.0))):
        ops = operators(p)
        loc[name] = [float(np.real(np.trace(ops["n"][0] @ r))) for r in evolve(ops, initial_state("localized"), times)]
    exact = (5 + 4 * np.cos(3 * Params().A * times / HBAR)) / 9
    write_csv("launch_site.csv", ["time_ns", "full_model", "free_hopping", "free_hopping_exact"], [times, loc["full"], loc["free"], exact])
    res["launch_site"] = {"t0.5ns": {"full": loc["full"][5], "free": loc["free"][5]},
                          "t1.0ns": {"full": loc["full"][10], "free": loc["free"][10]},
                          "free_vs_exact_max_difference": float(np.abs(np.array(loc["free"]) - exact).max())}

    # ---- 4. isolated hydrogen atom: quench rate against field -------------------------------
    fields = np.logspace(1, 6, 51)                                   # V/m
    exact_rate = [atom.quench_rate_exact(atom.mixing_energy(E)) * 1e9 for E in fields]
    weak_rate = [atom.quench_rate_weak_field(atom.mixing_energy(E)) * 1e9 for E in fields]
    write_csv("quench_rate_vs_field.csv", ["field_V_per_m", "rate_exact_per_s", "rate_weak_field_per_s"], [fields, exact_rate, weak_rate])
    V1 = atom.mixing_energy(1e3)
    rate1 = atom.quench_rate_weak_field(V1) * 1e9
    res["atom"] = {"lamb_gap_micro_eV": atom.DELTA, "dipole_length_orbital_a0": atom.dipole_length_orbital(),
                   "dipole_length_fine_structure_a0": atom.dipole_length_fine_structure(),
                   "V_at_1kV_per_m_micro_eV": V1, "hbar_gamma_over_2_micro_eV": HBAR * atom.GAMMA_2P / 2,
                   "quench_rate_1kV_per_m_per_s": rate1, "quench_time_1kV_per_m_micro_s": 1e6 / rate1,
                   "quench_rate_exact_1kV_per_m_per_s": atom.quench_rate_exact(V1) * 1e9,
                   "hop_time_ps_for_A_1meV": atom.HBAR_EV_S / 1e-3 * 1e12,
                   "rate_ratio_for_A_1meV": rate1 * atom.HBAR_EV_S / 1e-3,
                   "ceiling_rate_per_s": atom.GAMMA_2P_PER_S / 2, "field_where_V_equals_gap_kV_per_m": atom.DELTA / atom.mixing_energy(1.0) / 1e3,
                   "field_for_g_1.6_micro_eV_kV_per_m": 1.6 / V1,
                   "bare_electron_field_at_1nm_V_per_m": atom.coulomb_field(1e-9)}

    # ---- 5. ring model: emission within one hopping time against the ratio of the clocks ---
    p0 = Params()
    gq = atom.quench_rate_exact(p0.g, p0.delta, p0.gamma)
    ratios = 10 ** np.linspace(-2, 1, 13)
    full, frozen, As = [], [], []
    for r in ratios:
        A = HBAR * gq / r
        t = HBAR / A
        ops = operators(Params(A=A))
        rho = evolve(ops, initial_state("localized"), np.array([0.0, t]))[-1]
        full.append(1 - record_statistics(rho)[1][:, 0].sum())
        frozen.append(1 - atom.survival_static(p0.g, t, p0.delta, p0.gamma) ** 2)   # two neighboring atoms
        As.append(A)
    write_csv("emission_within_one_hop.csv", ["quench_rate_over_hopping_rate", "A_micro_eV", "full_model", "electron_held_on_site"],
              [ratios, As, full, frozen])
    res["clock_scan"] = {"quench_rate_g_1.6": gq, "hopping_rate_A_0.05": 0.05 / HBAR, "ratio_A_0.05": gq / (0.05 / HBAR),
                         "ratio_for_50_percent_emission": float(np.interp(0.5, full, ratios)),
                         "max_difference_full_vs_held": float(np.abs(np.array(full) - np.array(frozen)).max())}

    with open(os.path.join(OUT, "results.json"), "w") as f:
        json.dump(res, f, indent=2)
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
