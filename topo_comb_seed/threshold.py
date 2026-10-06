"""Comb-threshold finder: hold at fixed Delta_eff, watch the non-pumped-supermode fraction."""
import sys, time
import numpy as np, torch
sys.path.insert(0, ".")
from microring import CoupledLLESolver, H_zigzag, boundary_sites, pump_supermode
F2_LIST = [float(x) for x in (sys.argv[1:] or [2000, 5000, 12000, 30000])]
nx = ny = 6; J = 80.0; dt = 0.0025; de = 1.5
H = H_zigzag(nx, ny, J=J, phi=np.pi/4)
R = len(H); pump, drop = 0, R - (nx - 1)
lam_p, _, sp = pump_supermode(H, pump, edge=boundary_sites(H), edge_min=0.6)
lam, V = np.linalg.eigh(H)
Vi = torch.as_tensor(np.linalg.inv(V), dtype=torch.complex128)
t0 = time.time()
for F2 in F2_LIST:
    sol = CoupledLLESolver(H, N=64, dt=dt, Delta=de - lam_p, d2=0.0125, drive_modes=(1,),
                           pump_site=pump, kex_sites={pump: 1.0, drop: 1.0}, dtype=torch.complex128)
    sol.set_drive(float(np.sqrt(F2)), torch.zeros(1, 1, dtype=torch.float64))
    a = sol.random_state(1, generator=torch.Generator().manual_seed(1))
    a, _ = sol.evolve(a, int(150 / dt))
    b = torch.einsum("sr,brn->bsn", Vi, a)[0]            # supermode amplitudes (R, N)
    P = b.abs().pow(2)
    frac = float(1 - P[sp, 0] / P.sum())
    teeth = float(P[:, 1:].sum() / P.sum())              # power outside tooth mu=0
    print(f"F0^2={F2:7.0f}: non-pump-supermode fraction {frac:.4f}, power outside mu=0 {teeth:.4f}, "
          f"total {float(P.sum()):9.1f}  [{time.time()-t0:.0f}s]", flush=True)
