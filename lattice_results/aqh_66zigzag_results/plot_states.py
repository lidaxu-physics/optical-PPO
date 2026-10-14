"""
Pendulum on the zigzag 6x6 lattice operated in its different nonlinear states: one panel per state (learning curve of
every seed, frozen evaluation as a diamond at the end), and a summary panel of the evaluations. Reads the Pendulum
result files <state folder>/mr_topo*seed*.json (pattern_free also holds LunarLander runs, which are skipped).

    python lattice_results/aqh_66zigzag_results/plot_states.py        (from the repository root)
"""
import glob, json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from microring import plot_style as ps

ps.apply()
HERE = os.path.dirname(os.path.abspath(__file__))
STATES = [("pattern_free", "below threshold", "linear_P500_D1.00", ps.NAVY),
          ("rolls", "Turing roll", "turing14_P2500_D0.10", ps.GREEN),
          ("static_soliton", "static soliton", "soliton1_P3500_D4.90", ps.ORANGE),
          ("nested_soliton", "nested soliton", "nested_soliton_P6500_D6.60", ps.RED),
          ("chaos", "chaotic comb", "chaos_P6500_D4.50", ps.BROWN)]
# reference bands: the ring-free policies at the SAME budget (baselines/, mean +- std over the seeds), and the best single
# ring of the main results (Pasha's budget, 512k steps) as a line
def _band(prefix):
    ev = [json.load(open(f))["eval"]["mean"] for f in glob.glob(os.path.join(HERE, "baselines", f"{prefix}_seed*.json"))]
    return (np.mean(ev) - np.std(ev), np.mean(ev) + np.std(ev)) if ev else None
REF = [(name, *b) for name, b in (("linear, no ring (3 seeds)", _band("linear")), ("poly2, no ring (3 seeds)", _band("poly2"))) if b]
SINGLE_RING = -218                                                       # mean of the main results' stationary ring, 3 seeds
YLIM = (-1400, 0)

runs = {}
for folder, _, _, _ in STATES:
    files = sorted(glob.glob(os.path.join(HERE, folder, "mr_topo*seed*.json")))
    runs[folder] = [r for r in (json.load(open(f)) for f in files)
                    if r["config"]["env"].startswith("Pendulum") and r["config"].get("init_state") and r["config"].get("mini_comb") != "all"]
ALL_LINES = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(HERE, "pattern_free", "mr_topo_pendulum_all_seed*.json")))]
if not any(runs.values()):
    sys.exit("no finished runs yet")

fig, axes = plt.subplots(2, 3, figsize=(12.5, 7.6), sharey=True)
axes = axes.ravel()
for ax, (folder, label, key, col) in zip(axes, STATES):
    for name, lo, hi in REF:
        ax.axhspan(lo, hi, color=ps.GRID, alpha=0.5, lw=0, zorder=0)
    ax.axhline(SINGLE_RING, color=ps.MUTED, lw=0.9, ls=(0, (4, 3)), zorder=0)
    evals = []
    for i, r in enumerate(runs[folder]):
        x = np.array([u["env_steps"] for u in r["updates"]]) / 1e3
        y = np.array([u["mean_return"] for u in r["updates"]], dtype=float)
        ok = ~np.isnan(y)
        ax.plot(x[ok], y[ok], color=col, lw=1.1, alpha=0.55 + 0.15 * (i == 0), label=f"seed {r['config']['seed']}")
        if "eval" in r:
            evals.append(r["eval"]["mean"])
            ax.plot([x[ok][-1] + 4], [r["eval"]["mean"]], marker="D", ms=5.5, color=col, ls="", markeredgecolor="white", zorder=5)
    health = runs[folder][0].get("ring_health") if runs[folder] else None
    lines = [f"{label}  ({key})"]
    if evals:
        lines.append(f"frozen eval {np.mean(evals):.0f} ± {np.std(evals):.0f}, {len(evals)} seeds"
                     + (f", resets {100 * health['reset_fraction']:.1f} % of calls" if health else ""))
    else:
        lines.append("running" if not runs[folder] else "training, no evaluation yet")
    ax.set_title("\n".join(lines), fontsize=9, loc="left", color=ps.INK if runs[folder] else ps.MUTED)
    ax.set_xlabel("environment steps (thousands)", fontsize=8.5)
    ax.set_xlim(0, 212); ax.set_ylim(*YLIM); ax.tick_params(labelsize=8)
    if not runs[folder]:
        ax.text(0.5, 0.5, "running", transform=ax.transAxes, ha="center", va="center", fontsize=12, color=ps.MUTED)
    if folder == "pattern_free":
        ax.legend(fontsize=7.5, loc="lower right", frameon=False)
        pass
for ax in (axes[0], axes[3]):
    ax.set_ylabel("mean return per update", fontsize=8.5)

ax = axes[5]                                                            # the summary of the frozen evaluations
N_CAT = len(STATES) + (1 if ALL_LINES else 0)                       # categories of the summary panel
for name, lo, hi in REF:
    ax.axhspan(lo, hi, color=ps.GRID, alpha=0.5, lw=0, zorder=0)
    above = name.startswith("poly2")
    ax.text(N_CAT - 0.55, hi + 10 if above else (hi + lo) / 2, name, fontsize=7, ha="right", va="bottom" if above else "center", color=ps.INK2)
ax.axhline(SINGLE_RING, color=ps.MUTED, lw=0.9, ls=(0, (4, 3)), zorder=0)
ax.text(N_CAT - 0.55, SINGLE_RING - 12, "best single ring, main results (512k steps)", fontsize=7, ha="right", va="top", color=ps.MUTED)
for k, (folder, label, key, col) in enumerate(STATES):
    ev = [r["eval"]["mean"] for r in runs[folder] if "eval" in r]
    if ev:
        ax.plot([k] * len(ev), ev, marker="o", ms=6, ls="", color=col, alpha=0.8)
        ax.plot([k - 0.25, k + 0.25], [np.mean(ev)] * 2, color=col, lw=1.8)
        ax.text(k, max(ev) + 40, f"{np.mean(ev):.0f}", ha="center", fontsize=8, color=col)
    else:
        ax.text(k, YLIM[0] + 60, "running", ha="center", fontsize=7, color=ps.MUTED, rotation=90, va="bottom")
names = [s[1] for s in STATES]
if ALL_LINES:                                                           # below threshold reading every line of the band
    k, ev = len(STATES), [r["eval"]["mean"] for r in ALL_LINES if "eval" in r]
    ax.plot([k] * len(ev), ev, marker="o", ms=6, ls="", color=ps.NAVY, alpha=0.5, mfc="white")
    ax.plot([k - 0.25, k + 0.25], [np.mean(ev)] * 2, color=ps.NAVY, lw=1.8, alpha=0.6)
    ax.text(k, max(ev) + 40, f"{np.mean(ev):.0f}", ha="center", fontsize=8, color=ps.NAVY)
    names.append("below threshold," + chr(10) + "all %d lines" % len(ALL_LINES[0]["ring"]["fine"]["lines"]))
ax.set_xticks(range(len(names))); ax.set_xticklabels(names, fontsize=7.5, rotation=25, ha="right")
ax.set_xlim(-0.6, len(names) - 0.4); ax.set_ylim(*YLIM)
ax.set_title("frozen evaluation (16 episodes), every seed\ndots: seeds, bar: mean", fontsize=9, loc="left")
ax.tick_params(labelsize=8)

fig.suptitle("Pendulum on the zigzag 6 × 6 lattice, one nonlinear state at a time (seeded from lattice_seeds; "
             "16 envs × 32 steps × 400 updates)", fontsize=10.5, x=0.01, ha="left", y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.965))
out = os.path.join(HERE, os.environ.get("STATES_OUT", "states_pendulum.png"))
fig.savefig(out, dpi=150)
print("saved", out, "| finished states:", [s[0] for s in STATES if runs[s[0]]])
