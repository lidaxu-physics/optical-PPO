"""Load any catalogued nonlinear state of the zigzag 6x6 topological lattice and keep computing.

Everything needed is in this folder: `states.npz` (the field of each state) and `states.json`
(the parameters that hold it alive). No run data, no re-sweeping — one clone of this repo on
any machine is enough.

    python topo_comb_seed/states/load.py --list
    python topo_comb_seed/states/load.py --check                       # verify the archive
    python topo_comb_seed/states/load.py nested_soliton_P6500_D6.60 --evolve 60 --plot out.png

From your own code:

    from topo_comb_seed.states.load import load_state
    a, sol, meta = load_state("nested_soliton_P6500_D6.60")      # torch (1, R, N), solver
    a, mean_I = sol.evolve(a, n_steps, accumulate=True)           # ... and carry on

`a` is the longitudinal-mode amplitude A[ring, mu]; the field inside ring r is
psi_r(phi) = ifft(a[r], norm="forward"). `sol` is a CoupledLLESolver already carrying the
right detuning and drive, so the state is stationary/periodic under `sol.evolve` by
construction.
"""
import argparse, json, os, sys

import numpy as np
import torch

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))                  # repo root
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
from microring import CoupledLLESolver, H_zigzag, boundary_sites, pump_supermode, zigzag_sites

with open(os.path.join(_HERE, "states.json")) as _f:
    DOC = json.load(_f)
LAT, CATALOGUE = DOC["lattice"], DOC["states"]


def lattice():
    """(H, pump, drop, boundary order counter-clockwise, lambda of the pumped supermode)."""
    H = H_zigzag(LAT["nx"], LAT["ny"], J=LAT["J"], phi=LAT["phi"])
    R = len(H)
    pump, drop = LAT["pump_site"], LAT["drop_site"]
    bnd = boundary_sites(H)
    lam_p, _, _ = pump_supermode(H, pump, edge=bnd, edge_min=0.6)
    xy = zigzag_sites(LAT["nx"], LAT["ny"]).astype(float)
    order = bnd[np.argsort(np.arctan2(xy[bnd, 1] - xy[:, 1].mean(), xy[bnd, 0] - xy[:, 0].mean()))]
    return H, pump, drop, order, float(lam_p)


def load_state(name, dt=None, F0_squared=None, Delta_eff=None):
    """Return (a, solver, meta). Override F0_squared/Delta_eff to walk the state elsewhere."""
    if name not in CATALOGUE:
        raise KeyError(f"unknown state {name!r}; options: {', '.join(CATALOGUE)}")
    meta = dict(CATALOGUE[name])
    a0 = np.load(os.path.join(_HERE, "states.npz"))[name]
    H, pump, drop, _order, lam_p = lattice()
    F2 = meta["F0_squared"] if F0_squared is None else F0_squared
    de = meta["Delta_eff"] if Delta_eff is None else Delta_eff
    sol = CoupledLLESolver(H, N=a0.shape[-1], dt=dt or LAT["dt"], Delta=de - lam_p,
                           d2=LAT["d2"], drive_modes=(1,), pump_site=pump,
                           kex_sites={pump: LAT["kex"], drop: LAT["kex"]},
                           dtype=torch.complex128)
    sol.set_drive(float(np.sqrt(F2)), torch.zeros(1, 1, dtype=torch.float64))
    return torch.as_tensor(a0[None], dtype=torch.complex128), sol, meta


def diagnostics(a):
    """Lattice power, drop-ring phi contrast, boundary envelope contrast — the fingerprint."""
    _H, _pump, drop, order, _lam = lattice()
    arr = a[0].numpy() if hasattr(a, "numpy") else np.asarray(a)[0] if np.ndim(a) == 3 else np.asarray(a)
    psi2 = np.abs(np.fft.ifft(arr, axis=-1, norm="forward")) ** 2
    prof, env = psi2[drop], psi2[order].max(axis=1) - np.median(psi2[order], axis=1)
    return {"lattice_power": float((np.abs(arr) ** 2).sum()),
            "drop_phi_contrast": float(prof.max() / max(np.median(prof), 1e-12)),
            "boundary_env_contrast": float(env.max() / max(np.median(env), 1e-12))}


def _print_catalogue():
    print(f"{'state':30s} {'F0^2':>6s} {'D_eff':>6s} {'env':>4s}  {'power':>7s}  description")
    print("-" * 118)
    for n, m in CATALOGUE.items():
        print(f"{n:30s} {m['F0_squared']:6g} {m['Delta_eff']:6.2f} {m['n_envelopes']:4d}  "
              f"{m['fingerprint']['lattice_power']:7.2f}  {m['description'].splitlines()[0][:52]}")


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("state", nargs="?", help="state name (see --list)")
    p.add_argument("--list", action="store_true")
    p.add_argument("--check", action="store_true", help="re-derive every fingerprint")
    p.add_argument("--evolve", type=float, default=0.0, help="lifetimes to evolve after loading")
    p.add_argument("--F2", type=float, default=None, help="override pump power F0^2")
    p.add_argument("--de", type=float, default=None, help="override Delta_eff")
    p.add_argument("--plot", type=str, default=None, help="save a boundary x phi portrait here")
    args = p.parse_args()

    if args.list or (not args.state and not args.check):
        _print_catalogue(); return
    if args.check:
        bad = 0
        for n, m in CATALOGUE.items():
            a, _sol, _ = load_state(n)
            got, want = diagnostics(a), m["fingerprint"]
            ok = all(abs(got[k] - want[k]) <= 1e-6 * max(1.0, abs(want[k])) for k in want)
            bad += not ok
            print(f"{'ok ' if ok else 'BAD'} {n:30s} power {got['lattice_power']:8.3f}")
        print(f"\n{len(CATALOGUE) - bad}/{len(CATALOGUE)} states verified")
        sys.exit(1 if bad else 0)

    a, sol, meta = load_state(args.state, F0_squared=args.F2, Delta_eff=args.de)
    print(f"{args.state}: F0^2={meta['F0_squared']:g}, Delta_eff={meta['Delta_eff']:g}, "
          f"{meta['n_envelopes']} envelope(s)\n{meta['description']}\n")
    print("loaded      ", {k: round(v, 4) for k, v in diagnostics(a).items()})
    if args.evolve:
        a, _ = sol.evolve(a, int(args.evolve / LAT["dt"]))
        print(f"after {args.evolve:g} lt", {k: round(v, 4) for k, v in diagnostics(a).items()})
    if args.plot:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        _H, _pump, _drop, order, _lam = lattice()
        psi2 = np.abs(np.fft.ifft(a[0].numpy(), axis=-1, norm="forward")) ** 2
        fig, ax = plt.subplots(figsize=(5.5, 4))
        im = ax.imshow(psi2[order], origin="lower", aspect="auto", cmap="viridis",
                       extent=[0, 2 * np.pi, 0, len(order)], interpolation="nearest")
        ax.set_xlabel(r"$\varphi$ in each ring"); ax.set_ylabel("boundary ring (unrolled)")
        ax.set_title(args.state, fontsize=10); fig.colorbar(im, ax=ax, pad=0.02)
        fig.tight_layout(); fig.savefig(args.plot, dpi=140)
        print("saved", args.plot)


if __name__ == "__main__":
    main()
