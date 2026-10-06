"""Boundary-direction surgery: carve the travelling envelope train (state #9) to ONE envelope."""
import sys, time
import numpy as np, torch
sys.path.insert(0, ".")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from microring import CoupledLLESolver, H_zigzag, boundary_sites, pump_supermode, zigzag_sites
from microring import plot_style as ps
ps.apply()
import argparse
_p = argparse.ArgumentParser(); _p.add_argument("--de", type=float, default=7.0)
_p.add_argument("--win", type=float, default=2.5, help="taper scale in rings")
_p.add_argument("--keep", type=int, default=1); _p.add_argument("--tag", type=str, default="")
_p.add_argument("--init", type=str, default="topo_comb_seed/runs/P8000/verify_P8000_D+7.00_nested.npz:trace")
_p.add_argument("--relax_F2", type=float, default=None, help="heal at this pump for 100 lt, then restore")
_p.add_argument("--kick", type=float, default=0.0, help="boundary phase-ramp winding m: dev *= exp(i 2pi m s/nb)")
_p.add_argument("--F2", type=float, default=8000.0)
_a = _p.parse_args()
F2, de, dt = _a.F2, _a.de, 0.0025
H = H_zigzag(6, 6, J=80.0, phi=np.pi/4); R = len(H); pump, drop = 0, R-5
bnd = boundary_sites(H)
lam_p, _, _ = pump_supermode(H, pump, edge=bnd, edge_min=0.6)
_path, _key = _a.init.split(":")
_arr = np.load(_path)[_key]
a0 = _arr if _arr.ndim == 2 else _arr[-1]
psi = np.fft.ifft(a0, axis=-1, norm="forward")
N = psi.shape[1]
# unrolled boundary order
xy = zigzag_sites(6, 6).astype(float)
order = bnd[np.argsort(np.arctan2(xy[bnd,1]-xy[:,1].mean(), xy[bnd,0]-xy[:,0].mean()))]
# per-ring background ~ phi-mean; envelope weight = peak deviation power
bg = psi.mean(axis=1, keepdims=True)
dev = psi - bg
w_r = np.abs(dev[order]).max(axis=1)**2
nb = len(order)
if _a.kick:
    ph = np.zeros(R)
    ph[order] = 2*np.pi*_a.kick*np.arange(nb)/nb
    dev = dev * np.exp(1j*ph)[:, None]
# find _a.keep strongest well-separated blobs along the boundary
cands = list(np.argsort(w_r)[::-1])
centres = []
for c in cands:
    if all(min((c-x) % nb, (x-c) % nb) > 2*_a.win for x in centres):
        centres.append(int(c))
    if len(centres) >= _a.keep: break
keep = np.zeros(nb)
for c in centres:
    dist = np.minimum((np.arange(nb)-c) % nb, (c-np.arange(nb)) % nb)
    keep = np.maximum(keep, np.clip(1.6 - dist/_a.win, 0, 1))
w_full = np.zeros(R); w_full[order] = keep                     # interior rings: suppress deviation too
psi_new = bg + dev * w_full[:, None]
a = torch.as_tensor(np.fft.fft(psi_new, axis=-1, norm="forward")[None], dtype=torch.complex128)
sol = CoupledLLESolver(H, N=N, dt=dt, Delta=de - lam_p, d2=0.0125, drive_modes=(1,),
                       pump_site=pump, kex_sites={pump: 1.0, drop: 1.0}, dtype=torch.complex128)
sol.set_drive(float(np.sqrt(F2)), torch.zeros(1, 1, dtype=torch.float64))
t0 = time.time(); powers = []
if _a.relax_F2:
    sol_lo = CoupledLLESolver(H, N=N, dt=dt, Delta=de - lam_p, d2=0.0125, drive_modes=(1,),
                              pump_site=pump, kex_sites={pump: 1.0, drop: 1.0}, dtype=torch.complex128)
    sol_lo.set_drive(float(np.sqrt(_a.relax_F2)), torch.zeros(1, 1, dtype=torch.float64))
    a, _ = sol_lo.evolve(a, int(100.0/dt))
    for Fmid in np.linspace(_a.relax_F2, F2, 10)[1:-1]:        # adiabatic restore
        sol_m = CoupledLLESolver(H, N=N, dt=dt, Delta=de - lam_p, d2=0.0125, drive_modes=(1,),
                                 pump_site=pump, kex_sites={pump: 1.0, drop: 1.0}, dtype=torch.complex128)
        sol_m.set_drive(float(np.sqrt(Fmid)), torch.zeros(1, 1, dtype=torch.float64))
        a, _ = sol_m.evolve(a, int(10.0/dt))
    print(f"healing at P={_a.relax_F2:g} + restore done [{time.time()-t0:.0f}s]", flush=True)
for c in range(60):                                            # 300 lifetimes
    a, mI = sol.evolve(a, int(5.0/dt), accumulate=True, sample_every=10)
    powers.append((float(mI.sum()), float(mI[0, drop].sum())))
powers = np.array(powers)
a_fin, tr = sol.evolve_trace(a, int(40.0/dt), save_every=4)
tr = np.asarray(tr)[:, 0]
psi2 = np.abs(np.fft.ifft(a_fin[0].numpy(), axis=-1, norm="forward"))**2
# blob count along boundary at final instant
env = (psi2[order] - np.median(psi2[order], axis=1, keepdims=True)).max(axis=1)
thr = env.min() + 0.35*(env.max()-env.min()); above = env > thr
blobs = int(np.sum(np.diff(np.r_[above, above[0]].astype(int)) == 1)) if above.any() and not above.all() else 0
print(f"final blobs along boundary: {blobs}, lattice P = {powers[-1][0]:.2f} [{time.time()-t0:.0f}s]", flush=True)
sig = tr[:, drop, :] * np.hanning(tr.shape[0])[:, None]
A = np.fft.fftshift(np.fft.fft(sig, axis=0), axes=0)
nu = np.fft.fftshift(np.fft.fftfreq(tr.shape[0], d=4*dt)) * 2*np.pi
P = np.abs(A)**2
fig, axes = plt.subplots(1, 4, figsize=(18, 4.2))
axes[0].plot((np.arange(len(powers))+1)*5, powers[:,0], color=ps.BLUE, lw=1.3, label="lattice")
axes[0].plot((np.arange(len(powers))+1)*5, powers[:,1]*powers[:,0].max()/max(powers[:,1].max(),1e-12),
             color=ps.ORANGE, lw=1.0, label="drop (scaled)")
axes[0].set_xlabel("hold (lifetimes)"); axes[0].set_ylabel("power"); axes[0].legend(fontsize=8)
axes[0].set_title(f"boundary surgery keep={_a.keep}, win={_a.win}, $\\Delta_{{eff}}$={de}", fontsize=10)
im = axes[1].imshow(psi2[order], origin="lower", aspect="auto", cmap="viridis",
                    extent=[0, 2*np.pi, 0, len(order)], interpolation="nearest")
axes[1].grid(False); axes[1].set_xlabel("$\\varphi$"); axes[1].set_ylabel("boundary ring (unrolled)")
axes[1].set_title(f"final: {blobs} envelope(s) on the edge", fontsize=10)
fig.colorbar(im, ax=axes[1], pad=0.02).outline.set_visible(False)
mu_idx = 4
PdB = 10*np.log10(P[:, mu_idx]/max(P[:, mu_idx].max(),1e-300)+1e-12)
lam_axis = -nu - (de - lam_p) - 0.0125*16
o = np.argsort(lam_axis)
axes[2].plot(lam_axis[o], PdB[o], color=ps.BLUE, lw=0.9)
axes[2].set_xlim(-200, 200); axes[2].set_ylim(-80, 5)
axes[2].set_xlabel("$\\lambda$"); axes[2].set_ylabel("dB"); axes[2].set_title("fine spectrum, tooth $\\mu$=4", fontsize=10)
mu_sorted = np.fft.fftshift(np.fft.fftfreq(64, 1/64)).astype(int)
mu_grid = np.linspace(mu_sorted[0]-.5, mu_sorted[-1]+.5, 65)
ref = P.max(axis=0, keepdims=True)
PdB2 = 10*np.log10(P/np.maximum(ref,1e-300)+1e-7)[:, np.argsort(np.fft.fftfreq(64,1/64))]
dnu = nu[1]-nu[0]; nu_grid = np.concatenate([nu-dnu/2, [nu[-1]+dnu/2]])
im = axes[3].pcolormesh(nu_grid, mu_grid, PdB2.T, cmap="viridis", vmin=-70, vmax=0, shading="flat")
axes[3].grid(False); axes[3].set_xlim(-120,120); axes[3].set_xlabel("$\\nu$"); axes[3].set_ylabel("$\\mu$")
axes[3].set_title("dispersion map (row-norm)", fontsize=10)
fig.colorbar(im, ax=axes[3], pad=0.02).outline.set_visible(False)
fig.suptitle(f"Boundary surgery -- $F_0^2$ = {F2:g}, $\\Delta_{{eff}}$ = {de}, hold 300 lt",
             x=0.02, ha="left", fontweight="semibold")
fig.tight_layout(rect=(0,0,1,0.92))
fig.savefig(f"topo_comb_seed/runs/P{F2:g}/bsurg{_a.tag}.png", dpi=140)
np.savez_compressed(f"topo_comb_seed/runs/P{F2:g}/bsurg{_a.tag}.npz", a=a_fin[0].numpy(), powers=powers, psi2=psi2, trace=tr)
print(f"saved topo_comb_seed/runs/P{F2:g}/bsurg{_a.tag}.png", flush=True)
