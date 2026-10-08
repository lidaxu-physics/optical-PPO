"""
Nested (fine) structure of one longitudinal tooth, from the snapshot traces of a sweep.

    python lattice_results/lattice_seeds/nested_mu.py lattice_results/lattice_seeds/runs/sweep_F2400_J80.npz --mu 4

For every dumped snapshot (6 detunings along the ramp), the complex drop-ring amplitude of
tooth mu, a_{drop, mu}(t), is Fourier-transformed over the trace. In the pump frame a
supermode sigma inside that tooth oscillates at nu = -(Delta + d2 mu^2 + lambda_sigma), so
the fine spectrum is drawn against lambda = -nu - Delta - d2 mu^2: peaks then fall ON the
supermode ladder of H (orange ticks: boundary supermodes), comparable across snapshots.
Rolls: a few locked fine lines; chaos: a filled band; nested soliton: the whole ladder lit
with a smooth (sech^2-like) envelope.
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
p.add_argument("--mu", type=int, default=4, help="longitudinal tooth to unfold")
args = p.parse_args()

d = np.load(args.npz)
nx, ny, J, phi, F2, N, dt, d2, kex = d["config"]
lam_p, drop = float(d["lam_p"]), int(d["drop"])
H = H_zigzag(int(nx), int(ny), J=J, phi=phi)
lamH = np.linalg.eigvalsh(H)
w_edge = (np.abs(np.linalg.eigh(H)[1][boundary_sites(H)]) ** 2).sum(0)
snaps = sorted((k for k in d.files if k.startswith("snap_")), key=lambda k: float(k.split("_")[1]))
mu_idx = args.mu % int(N)
sample_dt = 4 * dt                                        # evolve_trace(save_every=4)

fig, axes = plt.subplots(2, 3, figsize=(14.5, 7), sharex=True, sharey=True)
for ax, key in zip(axes.ravel(), snaps):
    tr = d[key]                                           # (n_t, R, N) complex
    sig = tr[:, drop, mu_idx]
    sig = sig * np.hanning(len(sig))
    A = np.fft.fftshift(np.fft.fft(sig))
    nu = np.fft.fftshift(np.fft.fftfreq(len(sig), d=sample_dt)) * 2 * np.pi
    de = float(key.split("_")[1])                         # Delta_eff of the pumped supermode
    Delta = de - lam_p                                    # laser detuning
    lam_axis = -nu - Delta - d2 * args.mu ** 2            # map fine frequency onto the supermode ladder
    PdB = 10 * np.log10(np.abs(A) ** 2 / max(np.abs(A).max() ** 2, 1e-300) + 1e-12)
    o = np.argsort(lam_axis)
    ax.plot(lam_axis[o], PdB[o], color=ps.BLUE, lw=0.9)
    for l, we in zip(lamH, w_edge):
        ax.vlines(l, 2, 7, color=ps.ORANGE if we > 0.6 else ps.MUTED, lw=1.2 if we > 0.6 else 0.7)
    ax.set_title(f"$\\Delta_{{eff}}$ = {de:+.2f}", fontsize=10)
    ax.set_ylim(-75, 9)
    ax.set_xlim(lamH.min() - 8, lamH.max() + 8)
for ax in axes[1]:
    ax.set_xlabel("supermode eigenvalue $\\lambda$ ($\\kappa/2$)")
for ax in axes[:, 0]:
    ax.set_ylabel("fine power (dB re max)")
fig.suptitle(f"Nested structure of tooth $\\mu$ = {args.mu} at the drop ring -- "
             f"$F_0^2$ = {F2:g}, J = {J:g} (ticks: supermode ladder, orange = boundary)",
             x=0.02, ha="left", fontweight="semibold")
fig.tight_layout(rect=(0, 0, 1, 0.94))
out = args.npz.replace(".npz", f"_nested_mu{args.mu}.png")
fig.savefig(out, dpi=140)
print("saved", out)
