"""Continuation in PUMP POWER at fixed detuning: walk the pinned pair down in P."""
import sys, time
import numpy as np, torch
sys.path.insert(0, ".")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from microring import CoupledLLESolver, H_zigzag, boundary_sites, pump_supermode, zigzag_sites
from microring import plot_style as ps
ps.apply()
init, de, P0, P1, step = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
H = H_zigzag(6, 6, J=80.0, phi=np.pi/4); R = len(H); pump, drop = 0, R-5
bnd = boundary_sites(H)
lam_p, _, _ = pump_supermode(H, pump, edge=bnd, edge_min=0.6)
xy = zigzag_sites(6, 6).astype(float)
order = bnd[np.argsort(np.arctan2(xy[bnd,1]-xy[:,1].mean(), xy[bnd,0]-xy[:,0].mean()))]
path, key = init.split(":")
arr = np.load(path)[key]
a = torch.as_tensor((arr if arr.ndim == 2 else arr[-1])[None], dtype=torch.complex128)
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
Ps = np.arange(P0, P1 - 1e-9, -abs(step)) if P1 < P0 else np.arange(P0, P1 + 1e-9, abs(step))
sol0 = None
for P in Ps:
    sol = CoupledLLESolver(H, N=64, dt=0.0025, Delta=de-lam_p, d2=0.0125, drive_modes=(1,),
                           pump_site=pump, kex_sites={pump:1.0, drop:1.0}, dtype=torch.complex128)
    sol.set_drive(float(np.sqrt(P)), torch.zeros(1,1,dtype=torch.float64))
    a, mI = sol.evolve(a, int(25.0/0.0025), accumulate=True, sample_every=10)
    n, c = count(a)
    rows.append((P, float(mI.sum()), n, c))
    states[f"state_P{P:g}"] = a[0].numpy()
    if len(rows) % 8 == 1:
        print(f"P={P:g} Ptot={rows[-1][1]:7.2f} env={n} c={c:.1f} [{time.time()-t0:.0f}s]", flush=True)
    if c < 1.5: print(f"lost at P={P:g}", flush=True); break
rows = np.array(rows)
np.savez_compressed(f"topo_comb_seed/runs/P{P0:g}/powerwalk.npz", rows=rows, **states)
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
axes[0].plot(rows[:,0], rows[:,1], color=ps.BLUE, lw=1.5); axes[0].set_xlabel("pump P"); axes[0].set_ylabel("lattice power")
axes[1].step(rows[:,0], rows[:,2], where="mid", color=ps.ORANGE, lw=1.6); axes[1].set_xlabel("pump P"); axes[1].set_ylabel("envelopes")
fig.suptitle(f"Power-walk at $\\Delta_{{eff}}$ = {de}", x=0.02, ha="left", fontweight="semibold")
fig.tight_layout(rect=(0,0,1,0.92)); fig.savefig(f"topo_comb_seed/runs/P{P0:g}/powerwalk.png", dpi=140)
print(f"saved topo_comb_seed/runs/P{P0:g}/powerwalk.png", flush=True)
