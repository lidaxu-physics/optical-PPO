"""
Does the lattice compute, or does detecting powers suffice? The same training with the optical power turned down.

    python lattice_results/characterization/07_power_control.py   (from the repository root)

Pump and tones are scaled together (F0^2 and eps^2 by 1/10 and 1/100), so the drive keeps its shape and every four-wave
mixing product weakens with it; the features are standardised before the readout, so nothing else changes. In the
limit the lattice is linear and each edge line is the square of its own tone, (1 + s_k)^2: powers without any mixing
between the inputs. Also shown: the same lattice ABOVE the comb threshold (topo_chaos, 64 modes per ring).
From the result files of finished runs: Pendulum (AQH 4 x 4) and, once finished, LunarLander (zigzag 6 x 6).
"""

import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
DIRS = {"Pendulum": "lattice_results/aqh_44_results", "LunarLander": "lattice_results/aqh_66zigzag_results/pattern_free"}     # lattice runs
from microring import plot_style as ps

ps.apply()
RUNS = [("mr_topo", "full power, below threshold", ps.NAVY), ("mr_topo_p0.1", "power / 10", ps.BLUE), ("mr_topo_p0.01", "power / 100", ps.GREEN),
        ("mr_topo_chaos", "power × 8: chaotic comb, 64 modes", ps.RED), ("linear", "no lattice: linear on the inputs", ps.AXIS)]


def load(env, name):
    path = f"{DIRS[env]}/{name}_seed0.json"
    return json.load(open(path)) if os.path.exists(path) else None


def panel(ax_c, ax_e, env, title):
    found = [(load(env, n), lab, col) for n, lab, col in RUNS]
    for r, lab, col in found:
        if r is None:
            continue
        x = np.array([u["env_steps"] for u in r["updates"]]) / 1e3
        y = np.array([u["mean_return"] for u in r["updates"]], dtype=float)
        ok = ~np.isnan(y)
        if ok.sum() > 120:                                      # many noisy updates: a running mean
            ax_c.plot(x[ok][7:-7], np.convolve(y[ok], np.ones(15) / 15, mode="valid"), color=col, lw=1.8, label=f"{lab}:  {r['eval']['mean']:.0f}")
        else:
            ax_c.plot(x[ok], y[ok], color=col, lw=1.8, marker="o", ms=2.5, label=f"{lab}:  {r['eval']['mean']:.0f}")
        ret = np.sort(r["eval"]["returns"])
        ax_e.plot(np.linspace(0, 100, len(ret)), ret, color=col, lw=1.8, marker="o", ms=3)
    ax_c.set_xlabel("env steps (thousands)", fontsize=9); ax_c.set_ylabel("mean return of the episodes finished in the update", fontsize=9)
    ax_c.legend(fontsize=8, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, title="run:  frozen evaluation", title_fontsize=8)
    ax_c.set_title(title, fontsize=10)
    ax_e.set_xlabel("greedy evaluation episodes, sorted (%)", fontsize=9); ax_e.set_ylabel("return", fontsize=9)
    ax_e.set_title("the evaluation episodes of each", fontsize=10)
    return [lab for r, lab, _ in found if r is not None]


envs = [("Pendulum", "Pendulum, AQH 4×4, 4 edge lines")]
if load("LunarLander", "mr_topo_p0.1") and load("LunarLander", "mr_topo_p0.01"):
    envs.append(("LunarLander", "LunarLander, zigzag 6×6, 9 edge lines"))
fig, ax = plt.subplots(len(envs), 2, figsize=(13, 4.6 * len(envs)), squeeze=False, gridspec_kw=dict(width_ratios=[1.25, 1]))
fig.suptitle("The same training with the optical power turned down, and turned up past the comb threshold", x=0.01, ha="left", fontsize=12)
for (env, title), row in zip(envs, ax):
    print(env, panel(row[0], row[1], env, title))
for x in ax.ravel():
    x.tick_params(labelsize=8)
fig.tight_layout(rect=(0, 0, 1, 0.95 if len(envs) == 2 else 0.91))
out = "lattice_results/characterization/07_power_control.png"
fig.savefig(out, dpi=130)
print("saved", out)
