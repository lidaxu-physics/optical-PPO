"""Adiabatic red-walk of the travelling envelope train: watch the envelope count."""
import sys, time
import numpy as np, torch
sys.path.insert(0, ".")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from microring import CoupledLLESolver, H_zigzag, boundary_sites, pump_supermode, zigzag_sites
from microring import plot_style as ps
ps.apply()
F2, de0, de1, step, dt = 8000.0, 7.0, 11.0, 0.05, 0.0025
H = H_zigzag(6, 6, J=80.0, phi=np.pi/4); R = len(H); pump, drop = 0, R-5
bnd = boundary_sites(H)
lam_p, _, _ = pump_supermode(H, pump, edge=bnd, edge_min=0.6)
xy = zigzag_sites(6, 6).astype(float)
order = bnd[np.argsort(np.arctan2(xy[bnd,1]-xy[:,1].mean(), xy[bnd,0]-xy[:,0].mean()))]
a = torch.as_tensor(np.load(sys.argv[1] if len(sys.argv)>1 else "topo_comb_seed/runs/P8000/verify_P8000_D+7.00_nested.npz")["trace"][-1][None], dtype=torch.complex128)
def env_count(a):
    psi2 = np.abs(np.fft.ifft(a[0].numpy(), axis=-1, norm="forward"))**2
    env = psi2[order].max(axis=1) - np.median(psi2[order], axis=1)
    c = env.max() / max(np.median(env), 1e-12)
    if c < 2.0: return 0, c                                    # no real envelopes: smooth state
    thr = env.min() + 0.4*(env.max()-env.min()); above = env > thr
    n = int(np.sum(np.diff(np.r_[above, above[0]].astype(int)) == 1)) if above.any() and not above.all() else 0
    return n, c
rows, states = [], {}
t0 = time.time()
for de in np.arange(de0, de1+1e-9, step):
    sol = CoupledLLESolver(H, N=64, dt=dt, Delta=de-lam_p, d2=0.0125, drive_modes=(1,),
                           pump_site=pump, kex_sites={pump:1.0, drop:1.0}, dtype=torch.complex128)
    sol.set_drive(float(np.sqrt(F2)), torch.zeros(1,1,dtype=torch.float64))
    a, mI = sol.evolve(a, int(25.0/dt), accumulate=True, sample_every=10)
    n, c = env_count(a)
    rows.append((de, float(mI.sum()), n, c))
    states[f"state_{de:+.2f}"] = a[0].numpy()
    if len(rows) % 8 == 1:
        print(f"de={de:+.2f} P={rows[-1][1]:7.2f} envelopes={n} contrast={c:.1f} [{time.time()-t0:.0f}s]", flush=True)
    if c < 1.5 and len(rows) > 3:
        print("train lost (uniform/collapsed), stopping", flush=True); break
rows = np.array(rows)
np.savez_compressed(f"topo_comb_seed/runs/P{F2:g}/trainwalk.npz", rows=rows, **states)
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
axes[0].plot(rows[:,0], rows[:,1], color=ps.BLUE, lw=1.5)
axes[0].set_xlabel("$\\Delta_{eff}$"); axes[0].set_ylabel("lattice power")
axes[1].step(rows[:,0], rows[:,2], where="mid", color=ps.ORANGE, lw=1.6)
axes[1].set_xlabel("$\\Delta_{eff}$"); axes[1].set_ylabel("envelopes on the boundary")
fig.suptitle("Adiabatic red-walk of the envelope train, $F_0^2$ = 8000", x=0.02, ha="left", fontweight="semibold")
fig.tight_layout(rect=(0,0,1,0.92))
fig.savefig(f"topo_comb_seed/runs/P{F2:g}/trainwalk.png", dpi=140)
print(f"saved topo_comb_seed/runs/P{F2:g}/trainwalk.png", flush=True)
