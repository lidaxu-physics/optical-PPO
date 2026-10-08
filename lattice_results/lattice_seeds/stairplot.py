"""Staircase zoom of one pumpsweep npz: comb power vs detuning + drop-ring tooth heatmap.

    python lattice_results/lattice_seeds/stairplot.py lattice_results/lattice_seeds/runs/P6500/pumpmode_P6500_fine.npz
"""
import os, sys
import numpy as np
sys.path.insert(0, ".")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from microring import plot_style as ps
ps.apply()

d = np.load(sys.argv[1])
rows, F2 = d["rows"], float(d["F2"])
de = rows[:, 0]
fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True,
                         gridspec_kw={"height_ratios": [1, 1.3]})
axes[0].plot(de, rows[:, 3] - rows[:, 1], color=ps.BLUE, lw=1.3, label="comb power (hold mean)")
axes[0].plot(de, rows[:, 4] - rows[:, 2], color=ps.ORANGE, lw=0.8, alpha=0.7, label="hold end")
axes[0].set_ylabel(r"comb power $\sum|b|^2 - |b_{\sigma_p,0}|^2$")
axes[0].legend(fontsize=8)
if "teeth" in d:
    teeth = d["teeth"]                                   # (steps, N) drop-ring tooth powers, hold end
    N = teeth.shape[1]
    mu = np.fft.fftshift(np.fft.fftfreq(N, 1 / N)).astype(int)
    T = teeth[:, np.argsort(np.fft.fftfreq(N, 1 / N))]
    TdB = 10 * np.log10(T / max(T.max(), 1e-300) + 1e-12)
    dd = de[1] - de[0]
    de_grid = np.concatenate([de - dd / 2, [de[-1] + dd / 2]])
    mu_grid = np.linspace(mu[0] - .5, mu[-1] + .5, N + 1)
    im = axes[1].pcolormesh(de_grid, mu_grid, TdB.T, cmap="viridis", vmin=-80, vmax=0, shading="flat")
    axes[1].grid(False)
    axes[1].set_ylabel(r"tooth $\mu$ (drop ring)")
    fig.colorbar(im, ax=axes[1], pad=0.015, label="dB")
axes[1].set_xlabel(r"$\Delta_{eff}$")
fig.suptitle(f"Fine staircase -- $F_0^2$ = {F2:g}, comb power and drop comb vs detuning",
             x=0.02, ha="left", fontweight="semibold")
fig.tight_layout(rect=(0, 0, 1, 0.95))
out = sys.argv[1].replace(".npz", "_stairs.png")
fig.savefig(out, dpi=150)
print(f"saved {out}", flush=True)
