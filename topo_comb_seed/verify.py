"""
Stop-and-hold verification of a soliton candidate: ramp in from the blue, hold 500 lifetimes,
draw the portrait.

    python topo_comb_seed/verify.py --F2 3500 --target 4.8 [--ramp_from 3.6] [--T_hold 500]

Figure (runs/verify_F2<F2>_D<target>.png), four panels:
  1. power during the hold (lattice + drop; flat = stationary candidate);
  2. the nested-soliton portrait: |psi_r(phi)|^2 of the final state, boundary rings unrolled
     along the edge (y) x phi (x) -- ONE localized blob in both axes = single nested soliton,
     a stripe train = Turing/roll state, speckle = chaos;
  3. fine spectrum of tooth --mu over the last 40 lifetimes (locked lines vs noise);
  4. row-normalised dispersion map (nu x mu) of the same trace.
"""

import argparse, os, sys, time

import numpy as np
import torch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from microring import CoupledLLESolver, H_zigzag, boundary_sites, pump_supermode, zigzag_sites
from microring import plot_style as ps

ps.apply()
p = argparse.ArgumentParser()
p.add_argument("--F2", type=float, required=True)
p.add_argument("--target", type=float, required=True, help="Delta_eff to hold at")
p.add_argument("--ramp_from", type=float, default=None)
p.add_argument("--init", type=str, default=None, help="npz:snapKEY -- start from the last frame of that dumped trace instead of ramping")
p.add_argument("--T_hold", type=float, default=500.0)
p.add_argument("--mu", type=int, default=4)
p.add_argument("--J", type=float, default=80.0)
p.add_argument("--dt", type=float, default=0.0025)
p.add_argument("--nx", type=int, default=6)
p.add_argument("--ny", type=int, default=6)
p.add_argument("--tag", type=str, default="")
args = p.parse_args()

nx, ny, J, dt = args.nx, args.ny, args.J, args.dt
H = H_zigzag(nx, ny, J=J, phi=np.pi / 4)
R = len(H)
pump, drop = 0, R - (nx - 1)
bnd = boundary_sites(H)
lam_p, _, sp = pump_supermode(H, pump, edge=bnd, edge_min=0.6)
lamH = np.linalg.eigvalsh(H)
w_edge = (np.abs(np.linalg.eigh(H)[1][bnd]) ** 2).sum(0)
edge_lam = lamH[w_edge > 0.6]
F0 = float(np.sqrt(args.F2))
start = args.target - 1.2 if args.ramp_from is None else args.ramp_from
t0 = time.time()


def solver(de):
    s = CoupledLLESolver(H, N=64, dt=dt, Delta=de - lam_p, d2=0.0125, drive_modes=(1,),
                         pump_site=pump, kex_sites={pump: 1.0, drop: 1.0}, dtype=torch.complex128)
    s.set_drive(F0, torch.zeros(1, 1, dtype=torch.float64))
    return s


if args.init:
    path, key = args.init.split(":")
    tr0 = np.load(path)[key]
    frame = tr0 if tr0.ndim == 2 else tr0[-1]
    a = torch.as_tensor(frame[None, :, :], dtype=torch.complex128)
    print(f"initialised from {args.init} (last frame) [{time.time() - t0:.0f}s]", flush=True)
else:
    a = None
    for de in np.linspace(start, args.target, 40):               # approach ramp from the blue
        sol = solver(de)
        if a is None:
            a = sol.random_state(1, generator=torch.Generator().manual_seed(1))
        a, _ = sol.evolve(a, int(7.5 / dt))
    print(f"ramp {start:+.2f} -> {args.target:+.2f} done [{time.time() - t0:.0f}s]", flush=True)

sol = solver(args.target)
n_chunk, chunks = int(5.0 / dt), int(args.T_hold / 5.0)
powers = []
for _ in range(chunks):
    a, mI = sol.evolve(a, n_chunk, accumulate=True, sample_every=10)
    powers.append((float(mI.sum()), float(mI[0, drop].sum())))
powers = np.array(powers)
print(f"hold {args.T_hold:.0f} lifetimes done [{time.time() - t0:.0f}s]", flush=True)

a, tr = sol.evolve_trace(a, int(40.0 / dt), save_every=4)        # last 40 lifetimes
tr = np.asarray(tr)[:, 0]                                        # (n_t, R, N)
sample_dt = 4 * dt
psi2 = np.abs(np.fft.ifft(a[0].numpy(), axis=-1, norm="forward")) ** 2   # (R, N): |psi_r(phi)|^2

xy = zigzag_sites(nx, ny).astype(float)                          # boundary unrolled by angle
cx, cy = xy[:, 0].mean(), xy[:, 1].mean()
order = bnd[np.argsort(np.arctan2(xy[bnd, 1] - cy, xy[bnd, 0] - cx))]

sig = tr[:, drop, :] * np.hanning(tr.shape[0])[:, None]
A = np.fft.fftshift(np.fft.fft(sig, axis=0), axes=0)
nu = np.fft.fftshift(np.fft.fftfreq(tr.shape[0], d=sample_dt)) * 2 * np.pi
P = np.abs(A) ** 2

fig, axes = plt.subplots(1, 4, figsize=(18, 4.2))
ax = axes[0]
t_ax = (np.arange(len(powers)) + 1) * 5.0
ax.plot(t_ax, powers[:, 0], color=ps.BLUE, lw=1.4, label="lattice")
ax2 = powers[:, 1] * (powers[:, 0].max() / max(powers[:, 1].max(), 1e-12))
ax.plot(t_ax, ax2, color=ps.ORANGE, lw=1.1, label="drop (scaled)")
ax.set_xlabel("hold time (lifetimes)"); ax.set_ylabel("intracavity power")
ax.set_title(f"hold at $\\Delta_{{eff}}$ = {args.target:+.2f}", fontsize=10); ax.legend(fontsize=8)

ax = axes[1]
im = ax.imshow(psi2[order], origin="lower", aspect="auto", cmap="viridis",
               extent=[0, 2 * np.pi, 0, len(order)], interpolation="nearest")
ax.grid(False)
ax.set_xticks([0, np.pi, 2 * np.pi]); ax.set_xticklabels(["0", "$\\pi$", "$2\\pi$"])
ax.set_xlabel("$\\varphi$ in each ring"); ax.set_ylabel("boundary ring (unrolled)")
ax.set_title("final state $|\\psi_r(\\varphi)|^2$ on the edge", fontsize=10)
fig.colorbar(im, ax=ax, pad=0.02).outline.set_visible(False)

ax = axes[2]
mu_idx = args.mu % 64
lam_axis = -nu - (args.target - lam_p) - 0.0125 * args.mu ** 2
o = np.argsort(lam_axis)
PdB = 10 * np.log10(P[:, mu_idx] / max(P[:, mu_idx].max(), 1e-300) + 1e-12)
ax.plot(lam_axis[o], PdB[o], color=ps.BLUE, lw=0.9)
for l in edge_lam:
    ax.axvline(l, color=ps.ORANGE, lw=0.6, alpha=0.5, zorder=0)
ax.set_xlim(lamH.min() - 8, lamH.max() + 8); ax.set_ylim(-80, 5)
ax.set_xlabel("supermode eigenvalue $\\lambda$"); ax.set_ylabel("dB re max")
ax.set_title(f"fine spectrum of tooth $\\mu$ = {args.mu}", fontsize=10)

ax = axes[3]
mu_sorted = np.fft.fftshift(np.fft.fftfreq(64, 1.0 / 64)).astype(int)
mu_grid = np.linspace(mu_sorted[0] - 0.5, mu_sorted[-1] + 0.5, 65)
ref = P.max(axis=0, keepdims=True)
PdB2 = 10 * np.log10(P / np.maximum(ref, 1e-300) + 1e-7)
PdB2 = PdB2[:, np.argsort(np.fft.fftfreq(64, 1.0 / 64))]
dnu = nu[1] - nu[0]
nu_grid = np.concatenate([nu - dnu / 2, [nu[-1] + dnu / 2]])
im = ax.pcolormesh(nu_grid, mu_grid, PdB2.T, cmap="viridis", vmin=-70, vmax=0, shading="flat")
for l in edge_lam:
    ax.plot(-((args.target - lam_p) + 0.0125 * mu_sorted.astype(float) ** 2 + l), mu_sorted,
            color="w", lw=0.4, ls=":", alpha=0.5)
ax.grid(False); ax.set_xlim(-120, 120)
ax.set_xlabel("fine frequency $\\nu$"); ax.set_ylabel("$\\mu$")
ax.set_title("dispersion map (row-norm)", fontsize=10)
fig.colorbar(im, ax=ax, pad=0.02).outline.set_visible(False)

fig.suptitle(f"Stop-and-hold -- $F_0^2$ = {args.F2:g}, J = {J:g}, hold {args.T_hold:.0f} lifetimes",
             x=0.02, ha="left", fontweight="semibold")
fig.tight_layout(rect=(0, 0, 1, 0.92))
out = os.path.join(os.path.dirname(__file__), "runs", f"P{args.F2:g}",
                   f"verify_P{args.F2:g}_D{args.target:+.2f}" + (f"_{args.tag}" if args.tag else "") + ".png")
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=140)
np.savez_compressed(out.replace(".png", ".npz"), powers=powers, psi2=psi2, order=order,
                    trace=tr, nu=nu, target=args.target, F2=args.F2)
print("saved", out, flush=True)
