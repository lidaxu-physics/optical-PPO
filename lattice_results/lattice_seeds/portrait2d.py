"""
Explorer-style 2D portrait + formation animation of the single edge soliton.

    python lattice_results/lattice_seeds/portrait2d.py

Re-performs the keep-1 pulse surgery on the verified 3-pulse state, then evolves while
recording frames: the carved state relaxing into the locked single soliton. Each frame is
the explorer's "2D Power" rendering on the zigzag geometry: every ring drawn as a circle,
arc colour = |psi_r(phi)|^2 (hot colormap, black background), pump ring outlined blue,
drop ring red. Outputs: runs/P3500/soliton_2d.png (final state) and runs/P3500/soliton_formation.gif/.mp4.
"""

import os, sys, time

import numpy as np
import torch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.collections import LineCollection
from microring import CoupledLLESolver, H_zigzag, boundary_sites, pump_supermode, zigzag_sites

RING_R = 0.62
nx = ny = 6
J, F2, de, dt = 80.0, 3500.0, 4.90, 0.0025
H = H_zigzag(nx, ny, J=J, phi=np.pi / 4)
R = len(H)
pump, drop = 0, R - (nx - 1)
bnd = boundary_sites(H)
lam_p, _, _ = pump_supermode(H, pump, edge=bnd, edge_min=0.6)

# ---- rebuild the carved (keep-1) initial condition from the verified 3-pulse state
a3 = np.load(os.path.join(os.path.dirname(__file__), "runs", "P3500", "verify_P3500_D+4.80_warm.npz"))["trace"][-1]
psi = np.fft.ifft(a3, axis=-1, norm="forward")
N = psi.shape[1]
prof = np.abs(psi[drop]) ** 2
thr = prof.min() + 0.5 * (prof.max() - prof.min())
peaks = [i for i in range(N) if prof[i] > thr and prof[i] >= prof[(i - 1) % N] and prof[i] >= prof[(i + 1) % N]]
merged = []
for i in peaks:
    if not merged or min((i - merged[-1]) % N, (merged[-1] - i) % N) > 3:
        merged.append(i)
kept = merged[:1]
w = np.zeros(N)
for c in kept:
    d_idx = np.minimum((np.arange(N) - c) % N, (c - np.arange(N)) % N)
    w = np.maximum(w, np.clip(1.5 - d_idx / (N / 16), 0, 1))
bg_region = np.ones(N, bool)
for c in merged:
    d_idx = np.minimum((np.arange(N) - c) % N, (c - np.arange(N)) % N)
    bg_region &= d_idx > N / 10
bg = psi[:, bg_region].mean(axis=1, keepdims=True)
psi0 = bg + (psi - bg) * np.clip(w, 0, 1)[None, :]
a = torch.as_tensor(np.fft.fft(psi0, axis=-1, norm="forward")[None], dtype=torch.complex128)

sol = CoupledLLESolver(H, N=N, dt=dt, Delta=de - lam_p, d2=0.0125, drive_modes=(1,),
                       pump_site=pump, kex_sites={pump: 1.0, drop: 1.0}, dtype=torch.complex128)
sol.set_drive(float(np.sqrt(F2)), torch.zeros(1, 1, dtype=torch.float64))

frames, times, t_now = [], [], 0.0
schedule = [(0.5, 40), (2.0, 20)]                       # (step in lifetimes, how many)
frames.append(np.abs(np.fft.ifft(a[0].numpy(), axis=-1, norm="forward")) ** 2); times.append(0.0)
t0 = time.time()
for step_lt, count in schedule:
    for _ in range(count):
        a, _ = sol.evolve(a, int(step_lt / dt))
        t_now += step_lt
        frames.append(np.abs(np.fft.ifft(a[0].numpy(), axis=-1, norm="forward")) ** 2)
        times.append(t_now)
print(f"{len(frames)} frames to t = {t_now:.0f} lifetimes [{time.time() - t0:.0f}s]", flush=True)

# ---- explorer-style rendering on the zigzag geometry
xy = zigzag_sites(nx, ny)
px = 2 * xy[:, 0] - 1 + (xy[:, 1] % 2)                  # the explorer's A_zigzag site_xy
py = xy[:, 1].astype(float)
phi_f = np.linspace(0, 2 * np.pi, N, endpoint=False)
vmax = max(f.max() for f in frames) / 1.3
norm2d = plt.Normalize(0, vmax)

fig, ax = plt.subplots(figsize=(9.5, 8.2))
fig.patch.set_facecolor("black")
ax.set_facecolor("black"); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
for sp in ax.spines.values():
    sp.set_edgecolor("#1e2230")
arc_lcs, out_lcs = [], []
for r in range(R):
    pts = np.array([px[r] + RING_R * np.cos(phi_f), py[r] + RING_R * np.sin(phi_f)]).T.reshape(-1, 1, 2)
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
    lc = LineCollection(segs, cmap=plt.cm.hot, norm=norm2d, linewidth=2.6, zorder=2)
    lc.set_array(frames[0][r, :-1]); ax.add_collection(lc); arc_lcs.append(lc)
    ot = np.linspace(0, 2 * np.pi, 160)
    op = np.array([px[r] + RING_R * np.cos(ot), py[r] + RING_R * np.sin(ot)]).T.reshape(-1, 1, 2)
    os_ = np.concatenate([op[:-1], op[1:]], axis=1)
    if r == pump:
        ol = LineCollection(os_, colors="#4a9eff", linewidths=2.2, alpha=0.75, zorder=3)
    elif r == drop:
        ol = LineCollection(os_, colors="#ff4a6e", linewidths=2.2, alpha=0.75, zorder=3)
    else:
        ol = LineCollection(os_, cmap=plt.cm.hot, norm=norm2d, linewidths=0.7, zorder=3)
        ol.set_array(np.interp(ot, phi_f, frames[0][r])[:-1])
    ax.add_collection(ol); out_lcs.append(ol)
ax.set_xlim(px.min() - 1.3, px.max() + 1.3); ax.set_ylim(py.min() - 1.6, py.max() + 1.3)
title = ax.set_title("", color="#c8d0e7", fontsize=11, pad=6)
ax.text(px[pump], py[pump] - 1.15, "IN", color="#4a9eff", ha="center", fontsize=9)
ax.text(px[drop], py[drop] + 0.95, "OUT", color="#ff4a6e", ha="center", fontsize=9)


def draw(i):
    f = frames[i]
    for r in range(R):
        arc_lcs[r].set_array(f[r, :-1])
        if r not in (pump, drop):
            out_lcs[r].set_array(np.interp(np.linspace(0, 2 * np.pi, 160), phi_f, f[r])[:-1])
    title.set_text(f"2D Power  |  single-soliton formation after pulse surgery   "
                   f"t = {times[i]:5.1f} lifetimes   ($F_0^2$={F2:g}, J={J:g}, $\\Delta_{{eff}}$={de})")
    return []


OUT = os.path.join(os.path.dirname(__file__), "runs", "P3500")
draw(len(frames) - 1)
fig.savefig(os.path.join(OUT, "soliton_2d.png"), dpi=150, facecolor="black")
anim = animation.FuncAnimation(fig, draw, frames=len(frames), blit=False)
anim.save(os.path.join(OUT, "soliton_formation.gif"), writer=animation.PillowWriter(fps=8), dpi=100,
          savefig_kwargs={"facecolor": "black"})
try:
    anim.save(os.path.join(OUT, "soliton_formation.mp4"),
              writer=animation.FFMpegWriter(fps=12, bitrate=2400), dpi=130,
              savefig_kwargs={"facecolor": "black"})
except Exception as e:
    print("mp4 skipped:", e)
print("saved runs/P3500/soliton_2d.png, soliton_formation.gif/.mp4", flush=True)
