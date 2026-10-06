"""Diagonal continuation in (pump P, detuning Δ): carry the P3500 single-pulse hero
up to P=8000 along a path that stays red of the chaotic zone, watching for the moment
the edge-uniform pulse localises along the boundary (envelope count -> 1 = the prize).

    python topo_comb_seed/diagwalk.py <init npz:key> <P0> <de0> <P1> <de1> <steps> [hold_lt] [tag]
"""
import sys, time
import numpy as np, torch
sys.path.insert(0, ".")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from microring import CoupledLLESolver, H_zigzag, boundary_sites, pump_supermode, zigzag_sites
from microring import plot_style as ps
ps.apply()

init = sys.argv[1]
P0, de0, P1, de1 = map(float, sys.argv[2:6])
steps = int(sys.argv[6])
hold = float(sys.argv[7]) if len(sys.argv) > 7 else 30.0
tag = sys.argv[8] if len(sys.argv) > 8 else ""

H = H_zigzag(6, 6, J=80.0, phi=np.pi/4); R = len(H); pump, drop = 0, R-5
bnd = boundary_sites(H)
lam_p, _, _ = pump_supermode(H, pump, edge=bnd, edge_min=0.6)
xy = zigzag_sites(6, 6).astype(float)
order = bnd[np.argsort(np.arctan2(xy[bnd,1]-xy[:,1].mean(), xy[bnd,0]-xy[:,0].mean()))]
path, key = init.split(":")
arr = np.load(path)[key]
a = torch.as_tensor((arr if arr.ndim == 2 else arr[-1])[None], dtype=torch.complex128)
N = a.shape[-1]

def count(a):
    psi2 = np.abs(np.fft.ifft(a[0].numpy(), axis=-1, norm="forward"))**2
    env = psi2[order].max(axis=1) - np.median(psi2[order], axis=1)
    c = env.max()/max(np.median(env),1e-12)
    if c < 2.0: return 0, c
    thr = env.min()+0.4*(env.max()-env.min()); ab = env > thr
    n = int(np.sum(np.diff(np.r_[ab, ab[0]].astype(int))==1)) if ab.any() and not ab.all() else 0
    return n, c

rows, states = [], {}
t0 = time.time()
for k in range(steps + 1):
    f = k/steps
    P, de = P0+(P1-P0)*f, de0+(de1-de0)*f
    sol = CoupledLLESolver(H, N=N, dt=0.0025, Delta=de-lam_p, d2=0.0125, drive_modes=(1,),
                           pump_site=pump, kex_sites={pump:1.0, drop:1.0}, dtype=torch.complex128)
    sol.set_drive(float(np.sqrt(P)), torch.zeros(1,1,dtype=torch.float64))
    a, _ = sol.evolve(a, int((hold-5.0)/0.0025))
    drops = []
    for _ in range(50):                      # 5 lt drop trace for pk-pk (travelling indicator)
        a, _ = sol.evolve(a, int(0.1/0.0025))
        psi = np.fft.ifft(a[0].numpy(), axis=-1, norm="forward")
        drops.append(float((np.abs(psi[drop])**2).mean()))
    psi2 = np.abs(np.fft.ifft(a[0].numpy(), axis=-1, norm="forward"))**2
    Ptot = float(psi2.mean(axis=-1).sum())
    # phi-pulse contrast on the drop ring (is the fast-time pulse still alive?)
    prof = psi2[drop]
    cphi = float(prof.max()/max(np.median(prof),1e-12))
    n, c = count(a)
    pk = 100*(max(drops)-min(drops))/max(np.mean(drops),1e-12)
    rows.append((P, de, Ptot, n, c, cphi, pk))
    states[f"state_{k:02d}"] = a[0].numpy()
    print(f"step {k:2d}: P={P:6.0f} de={de:5.2f} Ptot={Ptot:7.2f} bnd-env={n} (c={c:4.1f}) "
          f"phi-contrast={cphi:5.1f} drop pk-pk={pk:6.2f}% [{time.time()-t0:.0f}s]", flush=True)
    if cphi < 1.5 and c < 1.5:
        print("state melted to CW — stopping", flush=True); break

rows = np.array(rows)
out = f"topo_comb_seed/runs/P{P1:g}/diagwalk{tag}"
np.savez_compressed(out + ".npz", rows=rows, **states)
fig, axes = plt.subplots(1, 3, figsize=(14, 3.8))
axes[0].plot(rows[:,0], rows[:,2], color=ps.BLUE, lw=1.5)
axes[0].set_xlabel("pump P"); axes[0].set_ylabel("lattice power")
axes[1].step(rows[:,0], rows[:,3], where="mid", color=ps.ORANGE, lw=1.6)
axes[1].plot(rows[:,0], rows[:,5], color=ps.BLUE, lw=1.2, alpha=0.7, label=r"$\varphi$ contrast")
axes[1].set_xlabel("pump P"); axes[1].set_ylabel("boundary envelopes / contrast"); axes[1].legend()
axes[2].plot(rows[:,0], rows[:,6], color=ps.RED if hasattr(ps,'RED') else 'crimson', lw=1.5)
axes[2].set_xlabel("pump P"); axes[2].set_ylabel("drop pk-pk (%)")
fig.suptitle(f"Diagonal walk ({P0:g},{de0:g}) -> ({P1:g},{de1:g})", x=0.02, ha="left", fontweight="semibold")
fig.tight_layout(rect=(0,0,1,0.9)); fig.savefig(out + ".png", dpi=140)
print(f"saved {out}.png", flush=True)
