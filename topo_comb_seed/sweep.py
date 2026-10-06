"""
Detuning-ramp engine for the nested-soliton hunt (see README.md in this folder).

    python topo_comb_seed/sweep.py --F2 400                       # one coarse sweep, figure + dumps
    python topo_comb_seed/sweep.py --F2 400 --J 80 --dt 0.0025    # the resolved-ladder variant

Pump only, on the central edge supermode of the zigzag lattice; Delta_eff of that supermode
is ramped blue -> red in equal holds. Per hold: the mean and final total power (whole lattice
and drop ring); at `--snap` detunings, a short trace a(t) is dumped to npz for the nested
spectrum/phase diagnostics. Everything in the units of the repo (kappa/2 = 1).
"""

import argparse, os, sys, time

import numpy as np
import torch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from microring import CoupledLLESolver, H_zigzag, boundary_sites, pump_supermode
from microring import plot_style as ps

ps.apply()
p = argparse.ArgumentParser()
p.add_argument("--nx", type=int, default=6)
p.add_argument("--ny", type=int, default=6)
p.add_argument("--J", type=float, default=80.0)
p.add_argument("--phi", type=float, default=float(np.pi / 4))
p.add_argument("--F2", type=float, default=400.0, help="pump power F0^2")
p.add_argument("--N", type=int, default=64, help="longitudinal modes per ring")
p.add_argument("--dt", type=float, default=0.0025)
p.add_argument("--d2", type=float, default=0.0125)
p.add_argument("--kex", type=float, default=1.0)
p.add_argument("--start", type=float, default=-3.0, help="ramp start, Delta_eff of the pumped supermode")
p.add_argument("--stop", type=float, default=7.0, help="ramp end")
p.add_argument("--steps", type=int, default=120, help="number of holds")
p.add_argument("--T_hold", type=float, default=7.5, help="lifetimes per hold (3000 steps at dt = 0.0025)")
p.add_argument("--snap", type=float, nargs="*", default=None,
               help="Delta_eff values at which to dump a short trace (default: 6 evenly spaced)")
p.add_argument("--T_trace", type=float, default=40.0, help="length of each dumped trace")
p.add_argument("--tag", type=str, default="")
p.add_argument("--sigma", type=int, default=None, help="pump supermode index (default: auto edge pick)")
args = p.parse_args()

OUT = os.path.join(os.path.dirname(__file__), "runs", f"P{args.F2:g}")   # P = pump power F0^2
os.makedirs(OUT, exist_ok=True)
name = f"sweep_P{args.F2:g}_J{args.J:g}" + (f"_{args.tag}" if args.tag else "")

H = H_zigzag(args.nx, args.ny, J=args.J, phi=args.phi)
R = len(H)
bnd = boundary_sites(H)
pump, drop = 0, R - (args.nx - 1)                       # the mini-comb corners (make_ring convention)
lam_p, overlap, sigma_p = pump_supermode(H, pump, edge=bnd, edge_min=0.6)
if args.sigma is not None:
    import numpy as _np
    _lam, _v = _np.linalg.eigh(H)
    sigma_p, lam_p = args.sigma, float(_lam[args.sigma])
    overlap = float(_np.abs(_v[pump, args.sigma])**2)
F0 = float(np.sqrt(args.F2))
deltas = np.linspace(args.start, args.stop, args.steps)                 # Delta_eff of supermode sigma_p
snaps = list(args.snap) if args.snap else list(np.linspace(args.start + 1, args.stop - 0.5, 6))
print(f"zigzag {args.nx}x{args.ny} (R={R}), pump supermode sigma={sigma_p} (lambda={lam_p:+.2f}, "
      f"overlap {overlap:.3f}), F0^2={args.F2:g}, J={args.J:g}, dt={args.dt}", flush=True)

a, t0 = None, time.time()
trace_rows, phi_rows, dumps = [], [], {}
for i, de in enumerate(deltas):
    sol = CoupledLLESolver(H, N=args.N, dt=args.dt, Delta=de - lam_p, d2=args.d2,
                           drive_modes=(1,), pump_site=pump,                     # dummy tone mode, driven at 0
                           kex_sites={pump: args.kex, drop: args.kex}, dtype=torch.complex128)
    sol.set_drive(F0, torch.zeros(1, 1, dtype=torch.float64))
    if a is None:
        a = sol.random_state(1, generator=torch.Generator().manual_seed(1))
    n = int(round(args.T_hold / args.dt))
    a, mean_I = sol.evolve(a, n, accumulate=True, sample_every=10)
    P = a.abs().pow(2).sum(-1)[0]                                               # (R,) final ring powers
    trace_rows.append((de, float(mean_I.sum()), float(mean_I[0, drop].sum()),
                       float(P.sum()), float(P[drop])))
    phi_rows.append(sol.to_phi(a)[0, drop].abs().pow(2).numpy())               # (N,) |psi(phi)|^2 of the drop ring at hold end
    hit = [s for s in snaps if abs(s - de) <= (deltas[1] - deltas[0]) / 2]
    if hit:
        a_s, tr = sol.evolve_trace(a, int(args.T_trace / args.dt), save_every=4)
        dumps[f"snap_{hit[0]:+.2f}"] = tr[:, 0].numpy()                         # (n_t, R, N) complex
        a = a_s
    if i % 10 == 0:
        print(f"  {i:4d}/{args.steps}  Delta_eff={de:+.2f}  <P>={trace_rows[-1][1]:9.2f}  [{time.time()-t0:.0f}s]", flush=True)

rows = np.array(trace_rows)
spatio = np.array(phi_rows)                                                     # (steps, N): phi resolved, one row per hold
np.savez_compressed(os.path.join(OUT, name + ".npz"),
                    rows=rows, spatio=spatio, deltas=deltas, lam_p=lam_p, sigma_p=sigma_p, pump=pump, drop=drop,
                    config=np.array([args.nx, args.ny, args.J, args.phi, args.F2, args.N, args.dt, args.d2, args.kex]),
                    **dumps)

fig, axes = plt.subplots(1, 3, figsize=(15.5, 4))
for ax, (cm, cf, lbl) in zip(axes[:2], ((1, 3, "whole lattice"), (2, 4, "drop ring"))):
    ax.plot(rows[:, 0], rows[:, cm], color=ps.BLUE, lw=1.6, label="hold mean")
    ax.plot(rows[:, 0], rows[:, cf], color=ps.ORANGE, lw=1.1, label="hold end")
    for sn in snaps:
        ax.axvline(sn, color=ps.GRID, lw=4, zorder=0)
    ax.set_xlabel("$\\Delta_{eff}$ of the pumped edge supermode")
    ax.set_ylabel(f"intracavity power, {lbl}")
    ax.legend(fontsize=8)
ax = axes[2]
im = ax.imshow(spatio.T, origin="lower", aspect="auto", cmap="viridis",
               extent=[deltas[0], deltas[-1], 0, 2*np.pi], interpolation="nearest")
ax.grid(False)
ax.set_yticks([0, np.pi, 2*np.pi]); ax.set_yticklabels(["0", "$\\pi$", "$2\\pi$"])
ax.set_xlabel("$\\Delta_{eff}$"); ax.set_ylabel("$\\varphi$ (drop ring)")
ax.set_title("$|\\psi_{drop}(\\varphi)|^2$ at each hold end", fontsize=10)
fig.colorbar(im, ax=ax, pad=0.02).outline.set_visible(False)
fig.suptitle(f"Detuning ramp, zigzag {args.nx}×{args.ny}, $F_0^2$ = {args.F2:g}, J = {args.J:g} "
             f"(grey bands: dumped traces)", x=0.02, ha="left", fontweight="semibold")
fig.tight_layout(rect=(0, 0, 1, 0.93))
fig.savefig(os.path.join(OUT, name + ".png"))
print(f"done in {time.time()-t0:.0f}s -> topo_comb_seed/runs/{name}.png / .npz", flush=True)
