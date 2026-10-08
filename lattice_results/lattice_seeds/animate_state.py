"""Explorer-style 2D portrait + circulation animation of any saved state.

    python lattice_results/lattice_seeds/animate_state.py --init "runs/P6500/verify_P6500_D+6.60_long.npz:trace" \
        --F2 6500 --de 6.6 --step_lt 1.0 --nfr 90 --out lattice_results/lattice_seeds/runs/P6500/nested_soliton \
        --label "single nested soliton circulating the topological edge"

Renders every ring as a circle on the zigzag geometry, arc colour = |psi_r(phi)|^2 (hot,
black background), pump ring blue, drop ring red. Saves <out>_2d.png (final frame),
<out>.gif and <out>.mp4.
"""
import argparse, os, sys, time
import numpy as np, torch
sys.path.insert(0, ".")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.collections import LineCollection
from microring import CoupledLLESolver, H_zigzag, boundary_sites, pump_supermode, zigzag_sites

p = argparse.ArgumentParser()
p.add_argument("--init", required=True, help="npz:key seed state")
p.add_argument("--F2", type=float, required=True)
p.add_argument("--de", type=float, required=True)
p.add_argument("--step_lt", type=float, default=1.0, help="lifetimes per frame")
p.add_argument("--nfr", type=int, default=90)
p.add_argument("--settle", type=float, default=5.0)
p.add_argument("--out", required=True, help="output prefix (no extension)")
p.add_argument("--label", default="state circulating the topological edge")
args = p.parse_args()

RING_R, dt = 0.62, 0.0025
H = H_zigzag(6, 6, J=80.0, phi=np.pi/4); R = len(H); pump, drop = 0, R-5
bnd = boundary_sites(H)
lam_p, _, _ = pump_supermode(H, pump, edge=bnd, edge_min=0.6)
path, key = args.init.split(":")
arr = np.load(path)[key]
a = torch.as_tensor((arr if arr.ndim == 2 else arr[-1])[None], dtype=torch.complex128)
N = a.shape[-1]
sol = CoupledLLESolver(H, N=N, dt=dt, Delta=args.de-lam_p, d2=0.0125, drive_modes=(1,),
                       pump_site=pump, kex_sites={pump:1.0, drop:1.0}, dtype=torch.complex128)
sol.set_drive(float(np.sqrt(args.F2)), torch.zeros(1,1,dtype=torch.float64))
a, _ = sol.evolve(a, int(args.settle/dt))
frames, times, t0 = [], [], time.time()
for i in range(args.nfr):
    frames.append(np.abs(np.fft.ifft(a[0].numpy(), axis=-1, norm="forward"))**2)
    times.append(i*args.step_lt)
    a, _ = sol.evolve(a, int(args.step_lt/dt))
print(f"{len(frames)} frames over {times[-1]:.1f} lifetimes [{time.time()-t0:.0f}s]", flush=True)

xy = zigzag_sites(6, 6)
px = 2*xy[:,0] - 1 + (xy[:,1] % 2); py = xy[:,1].astype(float)
phi_f = np.linspace(0, 2*np.pi, N, endpoint=False)
norm2d = plt.Normalize(0, max(f.max() for f in frames)/1.3)
fig, ax = plt.subplots(figsize=(9.5, 8.2)); fig.patch.set_facecolor("black")
ax.set_facecolor("black"); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
for sp in ax.spines.values(): sp.set_edgecolor("#1e2230")
arcs, outs = [], []
for r in range(R):
    pts = np.array([px[r]+RING_R*np.cos(phi_f), py[r]+RING_R*np.sin(phi_f)]).T.reshape(-1,1,2)
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
    lc = LineCollection(segs, cmap=plt.cm.hot, norm=norm2d, linewidth=2.6, zorder=2)
    lc.set_array(frames[0][r,:-1]); ax.add_collection(lc); arcs.append(lc)
    ot = np.linspace(0, 2*np.pi, 160)
    op = np.array([px[r]+RING_R*np.cos(ot), py[r]+RING_R*np.sin(ot)]).T.reshape(-1,1,2)
    os_ = np.concatenate([op[:-1], op[1:]], axis=1)
    if r == pump: ol = LineCollection(os_, colors="#4a9eff", linewidths=2.2, alpha=0.75, zorder=3)
    elif r == drop: ol = LineCollection(os_, colors="#ff4a6e", linewidths=2.2, alpha=0.75, zorder=3)
    else:
        ol = LineCollection(os_, cmap=plt.cm.hot, norm=norm2d, linewidths=0.7, zorder=3)
        ol.set_array(np.interp(ot, phi_f, frames[0][r])[:-1])
    ax.add_collection(ol); outs.append(ol)
ax.set_xlim(px.min()-1.3, px.max()+1.3); ax.set_ylim(py.min()-1.6, py.max()+1.3)
title = ax.set_title("", color="#c8d0e7", fontsize=11, pad=6)
ax.text(px[pump], py[pump]-1.15, "IN", color="#4a9eff", ha="center", fontsize=9)
ax.text(px[drop], py[drop]+0.95, "OUT", color="#ff4a6e", ha="center", fontsize=9)
def draw(i):
    f = frames[i]
    for r in range(R):
        arcs[r].set_array(f[r,:-1])
        if r not in (pump, drop):
            outs[r].set_array(np.interp(np.linspace(0,2*np.pi,160), phi_f, f[r])[:-1])
    title.set_text(f"2D Power  |  {args.label}   t = {times[i]:5.1f} lifetimes   "
                   f"($F_0^2$={args.F2:g}, $\\Delta_{{eff}}$={args.de:g})")
    return []
draw(len(frames)-1)
fig.savefig(args.out + "_2d.png", dpi=150, facecolor="black")
anim = animation.FuncAnimation(fig, draw, frames=len(frames), blit=False)
anim.save(args.out + ".gif", writer=animation.PillowWriter(fps=12), dpi=100,
          savefig_kwargs={"facecolor": "black"})
try:
    anim.save(args.out + ".mp4", writer=animation.FFMpegWriter(fps=15, bitrate=2400),
              dpi=130, savefig_kwargs={"facecolor": "black"})
except Exception as e:
    print("mp4 skipped:", e)
print(f"saved {args.out}_2d.png, .gif, .mp4", flush=True)
