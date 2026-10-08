"""Pulse surgery: carve N pulses out of the locked 3-pulse state, relax, hold, portrait."""
import sys, time
import numpy as np, torch
sys.path.insert(0, ".")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from microring import CoupledLLESolver, H_zigzag, boundary_sites, pump_supermode
from microring import plot_style as ps
ps.apply()
import argparse
_p = argparse.ArgumentParser()
_p.add_argument("keep", type=int, nargs="?", default=1)
_p.add_argument("--init", required=True, help="npz:key, or npz:frames@<Delta> for a sweep frame")
_p.add_argument("--de", type=float, default=4.90)
_p.add_argument("--F2", type=float, default=3500.0)
_p.add_argument("--T", type=float, default=500.0, help="lifetimes to relax/hold after surgery")
_a = _p.parse_args()
keep, de, F2, dt = _a.keep, _a.de, _a.F2, 0.0025
H = H_zigzag(6, 6, J=80.0, phi=np.pi/4); R = len(H); pump, drop = 0, R-5
bnd = boundary_sites(H)
lam_p, _, _ = pump_supermode(H, pump, edge=bnd, edge_min=0.6)
_path, _key = _a.init.split(":")
_z = np.load(_path)
if _key.startswith("frames@"):
    a0 = _z["frames"][int(np.argmin(np.abs(_z["deltas"] - float(_key.split("@")[1]))))]
else:
    a0 = _z[_key] if _z[_key].ndim == 2 else _z[_key][-1]
psi = np.fft.ifft(a0, axis=-1, norm="forward")                               # (R, N) fields
N = psi.shape[1]; phi = np.arange(N) * 2*np.pi/N
prof = np.abs(psi[drop])**2
# pulse positions: local maxima above mid-level
thr = prof.min() + 0.5*(prof.max()-prof.min())
peaks = [i for i in range(N) if prof[i] > thr and prof[i] >= prof[(i-1)%N] and prof[i] >= prof[(i+1)%N]]
# merge adjacent
merged = []
for i in peaks:
    if not merged or min((i - merged[-1]) % N, (merged[-1] - i) % N) > N // 20: merged.append(i)
print("pulse centres (phi idx):", merged, flush=True)
kept = merged[:keep]
# smooth keep-mask around kept pulses, width ~N/8
w = np.zeros(N)
for c in kept:
    d_idx = np.minimum((np.arange(N)-c) % N, (c-np.arange(N)) % N)
    w = np.maximum(w, np.clip(1.5 - d_idx/(N/16), 0, 1))
w = np.clip(w, 0, 1)
bg_region = np.ones(N, bool)
for c in merged:
    d_idx = np.minimum((np.arange(N)-c) % N, (c-np.arange(N)) % N)
    bg_region &= d_idx > N/10
bg = psi[:, bg_region].mean(axis=1, keepdims=True)                           # per-ring CW background
psi_new = bg + (psi - bg) * w[None, :]
a = torch.as_tensor(np.fft.fft(psi_new, axis=-1, norm="forward")[None], dtype=torch.complex128)
sol = CoupledLLESolver(H, N=N, dt=dt, Delta=de - lam_p, d2=0.0125, drive_modes=(1,),
                       pump_site=pump, kex_sites={pump: 1.0, drop: 1.0}, dtype=torch.complex128)
sol.set_drive(float(np.sqrt(F2)), torch.zeros(1, 1, dtype=torch.float64))
t0 = time.time(); powers = []
for c in range(int(_a.T / 5.0)):
    a, mI = sol.evolve(a, int(5.0/dt), accumulate=True, sample_every=10)
    powers.append(float(mI.sum()))
powers = np.array(powers)
psi2 = np.abs(np.fft.ifft(a[0].numpy(), axis=-1, norm="forward"))**2
from microring import zigzag_sites
xy = zigzag_sites(6, 6).astype(float)
order = bnd[np.argsort(np.arctan2(xy[bnd,1]-xy[:,1].mean(), xy[bnd,0]-xy[:,0].mean()))]
prof_f = psi2[drop]
thr = prof_f.min() + 0.35*(prof_f.max()-prof_f.min()); above = prof_f > thr
n_f = int(np.sum(np.diff(np.r_[above, above[0]].astype(int)) == 1)) if above.any() and not above.all() else 0
print(f"keep={keep}: final pulses={n_f}, P={powers[-1]:.2f}, drift={abs(powers[-1]-powers[len(powers)//5])/powers[-1]:.1e} [{time.time()-t0:.0f}s]", flush=True)
fig, axes = plt.subplots(1, 3, figsize=(14, 4))
axes[0].plot((np.arange(len(powers))+1)*5, powers, color=ps.BLUE, lw=1.4)
axes[0].set_xlabel("hold (lifetimes)"); axes[0].set_ylabel("lattice power"); axes[0].set_title(f"after surgery keep={keep}")
im = axes[1].imshow(psi2[order], origin="lower", aspect="auto", cmap="viridis",
                    extent=[0, 2*np.pi, 0, len(order)], interpolation="nearest")
axes[1].grid(False); axes[1].set_xlabel("$\\varphi$"); axes[1].set_ylabel("boundary ring")
axes[1].set_title(f"final: {n_f} pulse(s)"); fig.colorbar(im, ax=axes[1], pad=0.02).outline.set_visible(False)
axes[2].plot(phi, prof_f, color=ps.BLUE, lw=1.4)
axes[2].set_xlabel("$\\varphi$"); axes[2].set_ylabel("$|\\psi_{drop}(\\varphi)|^2$"); axes[2].set_title("drop-ring profile")
fig.suptitle(f"Pulse surgery at $\\Delta_{{eff}}$ = {de}, $F_0^2$ = {F2:g}", x=0.02, ha="left", fontweight="semibold")
fig.tight_layout(rect=(0,0,1,0.92))
fig.savefig(f"lattice_results/lattice_seeds_256modes/runs/P{F2:g}/surgery_keep{keep}.png", dpi=140)
np.savez_compressed(f"lattice_results/lattice_seeds_256modes/runs/P{F2:g}/surgery_keep{keep}.npz", a=a[0].numpy(), powers=powers, psi2=psi2)
print(f"saved lattice_results/lattice_seeds_256modes/runs/P{F2:g}/surgery_keep{keep}.png", flush=True)
