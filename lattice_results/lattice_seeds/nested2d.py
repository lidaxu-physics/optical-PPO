"""
2D dispersion heatmap of the nested comb: ALL teeth unfolded at once.

    python lattice_results/lattice_seeds/nested2d.py lattice_results/lattice_seeds/runs/sweep_F2400_J80.npz

One panel per dumped snapshot, the explorer's hmap convention: x = fine frequency nu in the
pump frame (slow-time FFT of a_{drop, mu}(t)), y = longitudinal mode mu, colour = power
(dB re the panel max), viridis.
A supermode sigma populated in tooth mu sits at nu = -(Delta + d2 mu^2 + lambda_sigma): the
thin dotted guides are these curves for the boundary supermodes -- the nested comb is the
lattice band structure drawn by the light itself, bent by the dispersion parabola d2 mu^2.
"""

import argparse, os, sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from microring import H_zigzag, boundary_sites
from microring import plot_style as ps

ps.apply()
p = argparse.ArgumentParser()
p.add_argument("npz")
p.add_argument("--floor", type=float, default=-70.0, help="dB floor of the colour scale")
args = p.parse_args()

d = np.load(args.npz)
nx, ny, J, phi, F2, N, dt, d2, kex = d["config"]
N = int(N)
lam_p, drop = float(d["lam_p"]), int(d["drop"])
H = H_zigzag(int(nx), int(ny), J=J, phi=phi)
lamH = np.linalg.eigvalsh(H)
w_edge = (np.abs(np.linalg.eigh(H)[1][boundary_sites(H)]) ** 2).sum(0)
edge_lam = lamH[w_edge > 0.6]
snaps = sorted((k for k in d.files if k.startswith("snap_")), key=lambda k: float(k.split("_")[1]))
sample_dt = 4 * dt
mu = np.fft.fftshift(np.fft.fftfreq(N, 1.0 / N)).astype(int)       # -N/2 .. N/2-1
mu_grid = np.linspace(mu[0] - 0.5, mu[-1] + 0.5, len(mu) + 1)

fig, axes = plt.subplots(2, 3, figsize=(15.5, 8), sharex=True, sharey=True)
for ax, key in zip(axes.ravel(), snaps):
    tr = d[key]                                                     # (n_t, R, N)
    sig = tr[:, drop, :] * np.hanning(tr.shape[0])[:, None]
    A = np.fft.fftshift(np.fft.fft(sig, axis=0), axes=0)            # (n_nu, N)
    nu = np.fft.fftshift(np.fft.fftfreq(tr.shape[0], d=sample_dt)) * 2 * np.pi
    P = np.abs(A) ** 2
    PdB = 10 * np.log10(P / P.max() + 10 ** (args.floor / 10))
    PdB = PdB[:, np.argsort(np.fft.fftfreq(N, 1.0 / N))]            # columns in mu order
    de = float(key.split("_")[1])
    Delta = de - lam_p
    dnu = nu[1] - nu[0]
    nu_grid = np.concatenate([nu - dnu / 2, [nu[-1] + dnu / 2]])
    im = ax.pcolormesh(nu_grid, mu_grid, PdB.T, cmap="viridis", vmin=args.floor, vmax=0, shading="flat")
    for l in edge_lam:                                              # guides: nu_sigma(mu)
        ax.plot(-(Delta + d2 * mu.astype(float) ** 2 + l), mu, color="w", lw=0.4, ls=":", alpha=0.5)
    ax.grid(False)
    ax.set_title(f"$\\Delta_{{eff}}$ = {de:+.2f}", fontsize=10)
    ax.set_xlim(nu.min(), nu.max())
for ax in axes[1]:
    ax.set_xlabel("fine frequency $\\nu$ ($\\kappa/2$)")
for ax in axes[:, 0]:
    ax.set_ylabel("longitudinal mode $\\mu$")
fig.suptitle(f"Nested comb, all teeth: drop-ring dispersion map -- $F_0^2$ = {F2:g}, J = {J:g} "
             f"(dotted: boundary-supermode curves $-(\\Delta + d_2\\mu^2 + \\lambda_\\sigma)$)",
             x=0.02, ha="left", fontweight="semibold")
cb = fig.colorbar(im, ax=axes, pad=0.015, fraction=0.03)
cb.set_label("dB re panel max"); cb.outline.set_visible(False)
out = args.npz.replace(".npz", "_nested2d.png")
fig.savefig(out, dpi=135)
print("saved", out)
