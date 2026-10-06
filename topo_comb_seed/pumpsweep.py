"""Cold detuning ramp recording the PUMP-MODE power: |b_{sigma_p, mu=0}|^2, the pumped
edge-supermode's pumped-tooth power (projection b = V^-1 a, threshold.py convention).

    python topo_comb_seed/pumpsweep.py --F2 3500 [--start -3 --stop 10 --steps 130 --T_hold 7.5]

Saves runs/P{F2}/pumpmode_P{F2}.npz with rows = (de, pump-mode mean, pump-mode end,
total supermode power mean, total end). Mean = 15 samples spread over the hold.
"""
import argparse, os, sys, time
import numpy as np, torch
sys.path.insert(0, ".")
from microring import CoupledLLESolver, H_zigzag, boundary_sites, pump_supermode

p = argparse.ArgumentParser()
p.add_argument("--F2", type=float, required=True)
p.add_argument("--start", type=float, default=-3.0)
p.add_argument("--stop", type=float, default=10.0)
p.add_argument("--steps", type=int, default=130)
p.add_argument("--T_hold", type=float, default=7.5)
p.add_argument("--dt", type=float, default=0.0025)
p.add_argument("--tag", type=str, default="")
args = p.parse_args()

H = H_zigzag(6, 6, J=80.0, phi=np.pi / 4); R = len(H)
pump, drop = 0, R - 5
bnd = boundary_sites(H)
lam_p, overlap, sp = pump_supermode(H, pump, edge=bnd, edge_min=0.6)
lam, V = np.linalg.eigh(H)
Vi = torch.as_tensor(np.linalg.inv(V), dtype=torch.complex128)
deltas = np.linspace(args.start, args.stop, args.steps)
F0 = float(np.sqrt(args.F2))
chunks = 15
n_chunk = int(round(args.T_hold / args.dt / chunks))
a, rows, teeth, t0 = None, [], [], time.time()
for i, de in enumerate(deltas):
    sol = CoupledLLESolver(H, N=64, dt=args.dt, Delta=de - lam_p, d2=0.0125, drive_modes=(1,),
                           pump_site=pump, kex_sites={pump: 1.0, drop: 1.0}, dtype=torch.complex128)
    sol.set_drive(F0, torch.zeros(1, 1, dtype=torch.float64))
    if a is None:
        a = sol.random_state(1, generator=torch.Generator().manual_seed(1))
    samples = []
    for _ in range(chunks):
        a, _ = sol.evolve(a, n_chunk)
        P = torch.einsum("sr,brn->bsn", Vi, a)[0].abs().pow(2)
        samples.append((float(P[sp, 0]), float(P.sum())))
    s = np.array(samples)
    rows.append((de, s[:, 0].mean(), s[-1, 0], s[:, 1].mean(), s[-1, 1]))
    teeth.append(np.abs(a[0, drop].numpy()) ** 2)        # drop-ring tooth powers at hold end
    if i % 20 == 0:
        print(f"  {i:4d}/{args.steps}  de={de:+.2f}  pumpmode={rows[-1][1]:9.2f}  [{time.time()-t0:.0f}s]", flush=True)
OUT = os.path.join(os.path.dirname(__file__), "runs", f"P{args.F2:g}")
os.makedirs(OUT, exist_ok=True)
name = f"pumpmode_P{args.F2:g}" + (f"_{args.tag}" if args.tag else "")
np.savez_compressed(os.path.join(OUT, name + ".npz"),
                    rows=np.array(rows), teeth=np.array(teeth), deltas=deltas,
                    lam_p=lam_p, sigma_p=sp, F2=args.F2)
print(f"saved {OUT}/{name}.npz [{time.time()-t0:.0f}s]", flush=True)
