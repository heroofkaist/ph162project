"""Draw the figures from the CSV files in results/.  Run after run_all.py."""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RES, FIG = os.path.join(HERE, "results"), os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)
TEAL, VIOLET, AMBER, GREY = "#0E8F9C", "#6C3FE0", "#D98E04", "#9AA3B5"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "#D5DAE5", "grid.linewidth": 0.6, "figure.dpi": 150})


def load(name):
    with open(os.path.join(RES, name)) as f:
        rows = list(csv.DictReader(f))
    return {k: [float(r[k]) for r in rows] for k in rows[0]}


def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, name))
    plt.close(fig)


full, off, free = load("timeseries_full.csv"), load("timeseries_emission_off.csv"), load("timeseries_mixing_off.csv")
fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharex=True)
for ax, key, title in zip(axes, ("coherence", "current_over_initial"), ("Coherence C of the electron", "Current / initial current")):
    ax.plot(free["time_ns"], free[key], color=GREY, ls=":", lw=1.5, label="Mixing off")
    ax.plot(off["time_ns"], off[key], color=VIOLET, ls="--", lw=1.1, label="Emission off")
    ax.plot(full["time_ns"], full[key], color=TEAL, lw=2.2, label="Full model")
    ax.set(title=title, xlabel="time (ns)", ylim=(0, 1.05), xlim=(0, 12))
axes[0].legend(frameon=False, loc="lower left")
save(fig, "coherence_and_current.png")

long = load("long_time.csv")
fig, ax = plt.subplots(figsize=(6.2, 3.8))
ax.plot(long["time_ns"], long["fraction_1s"], color=VIOLET, lw=2, label="Fraction of atoms in 1s")
ax.plot(long["time_ns"], long["current_over_initial"], color=TEAL, lw=2, label="Current / initial current")
ax.set(title="Full model over a longer time", xlabel="time (ns)", ylim=(0, 1.05), xlim=(0, 100))
ax.legend(frameon=False, loc="center right")
save(fig, "long_time.png")

q = load("quench_rate_vs_field.csv")
fig, ax = plt.subplots(figsize=(6.2, 3.8))
kv = [e / 1e3 for e in q["field_V_per_m"]]
ax.loglog(kv, q["rate_exact_per_s"], color=TEAL, lw=2.2, label="Quench rate of 2s (exact 2 x 2 model)")
ax.axhline(6.2649e8 / 2, color=VIOLET, ls="--", lw=1.5, label="Ceiling: half the 2p decay rate")
ax.axhline(1e-3 / 6.582119569e-16, color=AMBER, ls="--", lw=1.5, label="Hopping rate A / hbar for A = 1 meV")
ax.set(title="Hydrogen atom in a static field", xlabel="electric field (kV/m)", ylabel="rate (1/s)", ylim=(10, 1e13))
ax.legend(frameon=False, fontsize=9, loc="lower right")
save(fig, "quench_rate_vs_field.png")

s = load("emission_within_one_hop.csv")
fig, ax = plt.subplots(figsize=(6.2, 3.8))
ax.semilogx(s["quench_rate_over_hopping_rate"], s["electron_held_on_site"], color=VIOLET, ls="--", lw=1.5, label="Electron held on its site")
ax.semilogx(s["quench_rate_over_hopping_rate"], s["full_model"], "o-", color=TEAL, lw=2, ms=4, label="Full model")
ax.set(title="Probability of emission within one hopping time", xlabel="quench rate / hopping rate", ylim=(0, 1.05))
ax.legend(frameon=False, loc="upper left")
save(fig, "emission_within_one_hop.png")

l = load("launch_site.csv")
fig, ax = plt.subplots(figsize=(6.2, 3.8))
ax.plot(l["time_ns"], l["free_hopping"], color=GREY, ls=":", lw=1.5, label="Free hopping (mixing off)")
ax.plot(l["time_ns"], l["full_model"], color=TEAL, lw=2.2, label="Full model")
ax.set(title="Probability of finding the electron on its launch site", xlabel="time (ns)", ylim=(0, 1.05), xlim=(0, 12))
ax.legend(frameon=False, loc="upper right")
save(fig, "launch_site.png")
print("figures written to", FIG)
