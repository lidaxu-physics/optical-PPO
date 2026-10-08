"""
The one-figure summary of a sweep: 2 x 3 panels.

    python lattice_results/lattice_seeds_256modes/summary.py lattice_results/lattice_seeds_256modes/runs/sweep_F2400_J80.npz [--mu 4] [--panels 1 -1]

Row 1 (vs Delta_eff of the pumped supermode): lattice power (hold mean + end), drop-ring
power (same), and the viridis spatiotemporal map |psi_drop(phi)|^2.
Row 2 (from the dumped traces): the fine spectrum of tooth mu at all snapshots (colour =
Delta_eff), and the 2D dispersion heatmap (x = fine frequency nu, y = longitudinal mode mu,
explorer hmap convention) at two chosen snapshots (--panels, indices into the snapshot list).
"""

import argparse, os, sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from microring import H_zigzag, boundary_sites
from microring import plot_style as ps

ps.apply()
p = argparse.ArgumentParser()
p.add_argument("npz")
p.add_argument("--mu", type=int, default=4)
p.add_argument("--panels", type=int, nargs=2, default=(1, -1), help="snapshot indices for the two 2D maps")
p.add_argument("--floor", type=float, default=-70.0)
p.add_argument("--rownorm", action="store_true", help="normalise each mu row of the dispersion maps separately")
args = p.parse_args()

d = np.load(args.npz)
nx, ny, J, phi, F2, N, dt, d2, kex = d["config"]
N = int(N)
lam_p, drop = float(d["lam_p"]), int(d["drop"])
rows, spatio, deltas = d["rows"], d["spatio"], d["deltas"]
H = H_zigzag(int(nx), int(ny), J=J, phi=phi)
lamH = np.linalg.eigvalsh(H)
w_edge = (np.abs(np.linalg.eigh(H)[1][boundary_sites(H)]) ** 2).sum(0)
edge_lam = lamH[w_edge > 0.6]
snaps = sorted((k for k in d.files if k.startswith("snap_")), key=lambda k: float(k.split("_")[1]))
sample_dt = 8 * dt
mu_sorted = np.fft.fftshift(np.fft.fftfreq(N, 1.0 / N)).astype(int)
mu_grid = np.linspace(mu_sorted[0] - 0.5, mu_sorted[-1] + 0.5, len(mu_sorted) + 1)
mu_idx = args.mu % N


def fine(key):
    """(nu, A(nu, mu-sorted) power) of the drop ring for one snapshot."""
    tr = d[key]
    sig = tr[:, drop, :] * np.hanning(tr.shape[0])[:, None]
    A = np.fft.fftshift(np.fft.fft(sig, axis=0), axes=0)
    nu = np.fft.fftshift(np.fft.fftfreq(tr.shape[0], d=sample_dt)) * 2 * np.pi
    return nu, np.abs(A) ** 2


fig, axes = plt.subplots(2, 3, figsize=(16, 8.6))

for ax, (cm, cf, lbl) in zip(axes[0, :2], ((1, 3, "whole lattice"), (2, 4, "drop ring"))):
    ax.plot(rows[:, 0], rows[:, cm], color=ps.BLUE, lw=1.6, label="hold mean")
    ax.plot(rows[:, 0], rows[:, cf], color=ps.ORANGE, lw=1.0, label="hold end")
    for k in snaps:
        ax.axvline(float(k.split("_")[1]), color=ps.GRID, lw=3, zorder=0)
    ax.set_xlabel("$\\Delta_{eff}$"); ax.set_ylabel(f"intracavity power, {lbl}")
    ax.legend(fontsize=8)

ax = axes[0, 2]
im = ax.imshow(spatio.T, origin="lower", aspect="auto", cmap="viridis",
               extent=[deltas[0], deltas[-1], 0, 2 * np.pi], interpolation="nearest")
ax.grid(False); ax.set_yticks([0, np.pi, 2 * np.pi]); ax.set_yticklabels(["0", "$\\pi$", "$2\\pi$"])
ax.set_xlabel("$\\Delta_{eff}$"); ax.set_ylabel("$\\varphi$ (drop ring)")
ax.set_title("$|\\psi_{drop}(\\varphi)|^2$ at hold end", fontsize=10)
fig.colorbar(im, ax=ax, pad=0.02).outline.set_visible(False)

ax = axes[1, 0]
cmap = plt.cm.plasma
for i, key in enumerate(snaps):
    nu, P = fine(key)
    de = float(key.split("_")[1])
    lam_axis = -nu - (de - lam_p) - d2 * args.mu ** 2
    o = np.argsort(lam_axis)
    PdB = 10 * np.log10(P[:, mu_idx] / max(P[:, mu_idx].max(), 1e-300) + 1e-12)
    ax.plot(lam_axis[o], PdB[o] + 15 * i, color=cmap(i / max(len(snaps) - 1, 1)), lw=0.9,
            label=f"$\\Delta_{{eff}}$ = {de:+.2f}")
for l in edge_lam:
    ax.axvline(l, color=ps.ORANGE, lw=0.6, alpha=0.5, zorder=0)
ax.set_xlim(lamH.min() - 8, lamH.max() + 8)
ax.set_xlabel("supermode eigenvalue $\\lambda$"); ax.set_ylabel("fine power of tooth $\\mu$ (dB, offset)")
ax.set_title(f"nested structure of tooth $\\mu$ = {args.mu} (orange: boundary ladder)", fontsize=10)
ax.legend(fontsize=6, loc="upper right")

for ax, pi in zip(axes[1, 1:], args.panels):
    key = snaps[pi]
    nu, P = fine(key)
    ref = P.max(axis=0, keepdims=True) if args.rownorm else P.max()
    PdB = 10 * np.log10(P / np.maximum(ref, 1e-300) + 10 ** (args.floor / 10))
    PdB = PdB[:, np.argsort(np.fft.fftfreq(N, 1.0 / N))]
    de = float(key.split("_")[1])
    dnu = nu[1] - nu[0]
    nu_grid = np.concatenate([nu - dnu / 2, [nu[-1] + dnu / 2]])
    im = ax.pcolormesh(nu_grid, mu_grid, PdB.T, cmap="viridis", vmin=args.floor, vmax=0, shading="flat")
    for l in edge_lam:
        ax.plot(-((de - lam_p) + d2 * mu_sorted.astype(float) ** 2 + l), mu_sorted,
                color="w", lw=0.4, ls=":", alpha=0.5)
    ax.grid(False); ax.set_xlim(nu.min(), nu.max()); ax.set_ylim(mu_grid[0], mu_grid[-1])
    ax.set_xlabel("fine frequency $\\nu$"); ax.set_ylabel("longitudinal mode $\\mu$")
    ax.set_title(f"dispersion map at $\\Delta_{{eff}}$ = {de:+.2f}", fontsize=10)
    fig.colorbar(im, ax=ax, pad=0.02).outline.set_visible(False)

fig.suptitle(f"Sweep summary -- zigzag {int(nx)}×{int(ny)}, J = {J:g}, $F_0^2$ = {F2:g} "
             f"(grey bands: snapshot detunings; dotted white: cold boundary-supermode curves)",
             x=0.02, ha="left", fontweight="semibold")
fig.tight_layout(rect=(0, 0, 1, 0.95))
out = args.npz.replace(".npz", "_summary.png")
fig.savefig(out, dpi=135)
print("saved", out)
