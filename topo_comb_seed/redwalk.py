"""Red-walk continuation of a verified locked state: count the phi-pulses vs detuning."""
import sys, time, os
import numpy as np, torch
sys.path.insert(0, ".")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from microring import CoupledLLESolver, H_zigzag, boundary_sites, pump_supermode
from microring import plot_style as ps
ps.apply()
init_npz, F2, de0, de1, step = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), 0.05
H = H_zigzag(6, 6, J=80.0, phi=np.pi/4); R = len(H); pump, drop = 0, R-5
lam_p, _, _ = pump_supermode(H, pump, edge=boundary_sites(H), edge_min=0.6)
F0, dt = float(np.sqrt(F2)), 0.0025
a = torch.as_tensor(np.load(init_npz)["trace"][-1][None], dtype=torch.complex128)
rows, states = [], {}
t0 = time.time()
des = np.arange(de0, de1 + 1e-9, step)
for de in des:
    sol = CoupledLLESolver(H, N=64, dt=dt, Delta=de - lam_p, d2=0.0125, drive_modes=(1,),
                           pump_site=pump, kex_sites={pump: 1.0, drop: 1.0}, dtype=torch.complex128)
    sol.set_drive(F0, torch.zeros(1, 1, dtype=torch.float64))
    a, mI = sol.evolve(a, int(25.0 / dt), accumulate=True, sample_every=10)
    prof = np.abs(np.fft.ifft(a[0, drop].numpy(), norm="forward"))**2
    thr = prof.min() + 0.35*(prof.max()-prof.min())
    above = prof > thr
    pulses = int(np.sum(np.diff(np.r_[above, above[0]].astype(int)) == 1)) if above.any() and not above.all() else 0
    rows.append((de, float(mI.sum()), pulses, float(prof.max()/max(prof.mean(),1e-12))))
    if abs(de*20 - round(de*20)) < 1e-6 and abs(de - round(de,1)) < 1e-6:  # every 0.1... keep every step cheap anyway
        pass
    states[f"state_{de:+.2f}"] = a[0].numpy()
    if len(rows) % 10 == 1:
        print(f"de={de:+.2f} P={rows[-1][1]:7.2f} pulses={pulses} contrast={rows[-1][3]:.1f} [{time.time()-t0:.0f}s]", flush=True)
    if rows[-1][1] < 13:   # collapsed to CW
        print("collapsed to CW, stopping", flush=True); break
rows = np.array(rows)
np.savez_compressed(f"topo_comb_seed/runs/P{F2:g}/redwalk_P{F2:g}.npz", rows=rows, lam_p=lam_p, **states)
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
axes[0].plot(rows[:,0], rows[:,1], color=ps.BLUE, lw=1.6)
axes[0].set_xlabel("$\\Delta_{eff}$"); axes[0].set_ylabel("lattice power")
axes[1].step(rows[:,0], rows[:,2], where="mid", color=ps.ORANGE, lw=1.6)
axes[1].set_xlabel("$\\Delta_{eff}$"); axes[1].set_ylabel("pulses per ring (drop)")
axes[1].set_ylim(-0.3, max(rows[:,2])+0.8)
fig.suptitle(f"Red walk from the locked state, $F_0^2$ = {F2:g}", x=0.02, ha="left", fontweight="semibold")
fig.tight_layout(rect=(0,0,1,0.93))
fig.savefig(f"topo_comb_seed/runs/P{F2:g}/redwalk_P{F2:g}.png", dpi=140)
print("saved topo_comb_seed/runs/P%g/redwalk_P%g.png" % (F2, F2), flush=True)
