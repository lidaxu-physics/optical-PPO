"""Waterfall of pump-mode power vs detuning: one panel per pump power P = F0^2,
stacked in a single column, LOWEST pump at the BOTTOM, shared detuning axis.

    python lattice_results/lattice_seeds/waterfall.py          # pump-mode power |b_{sigma_p,0}|^2
    python lattice_results/lattice_seeds/waterfall.py comb     # comb power: total supermode power - pump tooth
"""
import glob, os, sys
import numpy as np
sys.path.insert(0, ".")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from microring import plot_style as ps
ps.apply()

COMB = len(sys.argv) > 1 and sys.argv[1] == "comb"
import re
files = [f for f in glob.glob(os.path.join(os.path.dirname(__file__), "runs", "P*", "pumpmode_P*.npz"))
         if re.fullmatch(r"pumpmode_P[\d.]+\.npz", os.path.basename(f))]   # untagged (full-range) only
files.sort(key=lambda f: float(np.load(f)["F2"]))
if not files:
    sys.exit("no pumpmode npz found")
n = len(files)
fig, axes = plt.subplots(n, 1, figsize=(8.5, 0.95 * n + 1.4), sharex=True)
axes = np.atleast_1d(axes)
cmap = plt.cm.viridis
for k, f in enumerate(files):
    d = np.load(f)
    rows, F2 = d["rows"], float(d["F2"])
    ax = axes[n - 1 - k]                                  # lowest F at the bottom
    c = cmap(0.1 + 0.8 * k / max(n - 1, 1))
    y_mean = rows[:, 3] - rows[:, 1] if COMB else rows[:, 1]
    y_end = rows[:, 4] - rows[:, 2] if COMB else rows[:, 2]
    ax.plot(rows[:, 0], y_mean, color=c, lw=1.2)                          # hold mean
    ax.plot(rows[:, 0], y_end, color=c, lw=0.6, alpha=0.45)               # hold end (jitter = chaos)
    ax.text(1.005, 0.5, f"P = {F2:g}", transform=ax.transAxes, va="center", fontsize=8.5, color=c)
    ax.tick_params(labelsize=7)
    ax.margins(x=0)
axes[-1].set_xlabel(r"$\Delta_{eff}$")
ylab = (r"comb power $\sum|b|^2 - |b_{\sigma_p,\,\mu=0}|^2$ (per panel scale)" if COMB
        else r"pump-mode power $|b_{\sigma_p,\,\mu=0}|^2$ (per panel scale)")
fig.text(0.015, 0.5, ylab, va="center", rotation="vertical", fontsize=10)
fig.suptitle(("Comb" if COMB else "Pump-mode") +
             " power vs detuning -- cold blue->red sweeps, zigzag 6x6, J = 80",
             x=0.02, ha="left", fontweight="semibold")
fig.tight_layout(rect=(0.04, 0, 0.93, 0.96), h_pad=0.4)
out = os.path.join(os.path.dirname(__file__), "runs",
                   "combpower_waterfall.png" if COMB else "pumpmode_waterfall.png")
fig.savefig(out, dpi=150)
print(f"saved {out}", flush=True)
