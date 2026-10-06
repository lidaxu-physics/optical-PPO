"""Extract the catalogued nonlinear states from the (huge) raw sweep/verify runs into one
small portable archive: states.npz (one (R, N) complex frame each) + states.json (parameters
and verification fingerprints).

Provenance only — it needs the ~10 GB `soliton/runs/` tree of the hunt worktree, which is NOT
in the repo. Run it from the hunt worktree root:

    python topo_comb_seed/states/build_states.py <out_dir>

Everyone else just uses `states.npz` via `load.py`; this file records exactly where each
state came from.
"""
import json, os, sys

import numpy as np

sys.path.insert(0, ".")
from microring import H_zigzag, boundary_sites, pump_supermode, zigzag_sites

OUT = sys.argv[1] if len(sys.argv) > 1 else "topo_comb_seed/states"

# --- the lattice these states live on (identical for every entry) -----------------------
NX = NY = 6
J, PHI, D2, KEX, N_MODES, DT = 80.0, float(np.pi / 4), 0.0125, 1.0, 64, 0.0025
H = H_zigzag(NX, NY, J=J, phi=PHI)
R = len(H)
PUMP, DROP = 0, R - (NX - 1)
BND = boundary_sites(H)
LAM_P, OVERLAP, SIGMA_P = pump_supermode(H, PUMP, edge=BND, edge_min=0.6)
xy = zigzag_sites(NX, NY).astype(float)
ORDER = BND[np.argsort(np.arctan2(xy[BND, 1] - xy[:, 1].mean(), xy[BND, 0] - xy[:, 0].mean()))]

# --- catalogue: name -> (source npz, key, F0^2, Delta_eff, n_envelopes, blurb) ----------
SPEC = [
    # P = 500 — below threshold, linear reference
    ("linear_P500_D1.00", "soliton/runs/P500/sweep_P500_J80_fine.npz", "snap_+1.00", 500, 1.00, 0,
     "Below threshold: linear tilted resonance, no comb. Reference state."),
    # P = 3500 — Turing-3 shelf and the stationary k=0 soliton family
    ("chaos_P3500_D3.10", "soliton/runs/P3500/sweep_P3500_J80.npz", "snap_+3.10", 3500, 3.10, 0,
     "Developed MI chaos: power wanders, fine spectra filled."),
    ("minicomb_P3500_D1.40", "soliton/runs/P3500/sweep_P3500_J80.npz", "snap_+1.40", 3500, 1.40, 0,
     "Locked mini-comb candidate inside the dynamic zone: discrete equidistant fine lines in "
     "tooth mu=4 under a smooth envelope. Nesting NOT established."),
    ("cw_P3500_D4.80", "soliton/runs/P3500/verify_P3500_D+4.80.npz", "trace", 3500, 4.80, 0,
     "CW lower branch at the same point as turing3_P3500_D4.80 (cold start) — multistability."),
    ("turing3_P3500_D4.80", "soliton/runs/P3500/verify_P3500_D+4.80_warm.npz", "trace", 3500, 4.80, 3,
     "Edge-uniform Turing-3: 3 equidistant phi-pulses, every boundary ring in phase (k=0)."),
    ("soliton2_P3500_D4.90", "soliton/runs/P3500/surgery_keep2.npz", "a", 3500, 4.90, 2,
     "Two stationary pulses per ring; power = CW + 2 x 3.6 (quantisation family)."),
    ("soliton1_P3500_D4.90", "soliton/runs/P3500/surgery_keep1.npz", "a", 3500, 4.90, 1,
     "HERO-1: single stationary soliton per boundary ring, k=0 along the edge (NOT nested). "
     "Prepared by pulse surgery on turing3; drift 1e-15 over 800 lifetimes."),
    # P = 6500 — the soliton staircase down to the single nested travelling soliton
    ("chaos_P6500_D1.00", "soliton/runs/P6500/verify_P6500_D+1.00_blue10.npz", "trace", 6500, 1.00, 0,
     "Chaos already owns the blue side at this pump (no Turing window here)."),
    ("chaos_P6500_D4.50", "soliton/runs/P6500/verify_P6500_D+4.50_chaos45.npz", "trace", 6500, 4.50, 0,
     "Developed chaos just before the collapse at 5.2."),
    ("gas_P6500_D5.50", "soliton/runs/P6500/verify_P6500_D+5.50_shelf55.npz", "trace", 6500, 5.50, 10,
     "Dense edge soliton gas, metastable: evaporated one envelope at t ~ 250 lifetimes."),
    ("gas_P6500_D5.90", "soliton/runs/P6500/verify_P6500_D+5.90_shelf59.npz", "trace", 6500, 5.90, 9,
     "Edge soliton gas: ~9 pulses localised in both ring and phi, disordered."),
    ("crystal_P6500_D6.30", "soliton/runs/P6500/verify_P6500_D+6.30_shelf63.npz", "trace", 6500, 6.30, 5,
     "Nested soliton crystal: ~5 pulses, clean periodic sawtooth in the power, rigid 2D comb."),
    ("nested_soliton_P6500_D6.60", "soliton/runs/P6500/verify_P6500_D+6.60_long.npz", "trace", 6500, 6.60, 1,
     "*** THE RESULT *** Single nested travelling edge soliton: ONE envelope localised in "
     "boundary-ring AND phi, circulating the 20-ring edge once every 0.56 lifetimes (35.7 rings/lt, chiral, one "
     "direction), rigid equidistant comb across all teeth. Verified 1100 lifetimes."),
    ("cw_P6500_D6.90", "soliton/runs/P6500/verify_P6500_D+6.90_shelf69.npz", "trace", 6500, 6.90, 0,
     "CW: one step red of the single soliton the shelf is over."),
]


def fingerprint(a):
    """Cheap, order-independent numbers a reload can be checked against."""
    psi2 = np.abs(np.fft.ifft(a, axis=-1, norm="forward")) ** 2
    prof = psi2[DROP]
    env = psi2[ORDER].max(axis=1) - np.median(psi2[ORDER], axis=1)
    return {
        "lattice_power": float((np.abs(a) ** 2).sum()),
        "drop_phi_contrast": float(prof.max() / max(np.median(prof), 1e-12)),
        "boundary_env_contrast": float(env.max() / max(np.median(env), 1e-12)),
    }


states, meta = {}, {}
for name, src, key, F2, de, n_env, blurb in SPEC:
    arr = np.load(src)[key]
    a = (arr if arr.ndim == 2 else arr[-1]).astype(np.complex128)
    states[name] = a
    meta[name] = {
        "F0_squared": F2, "Delta_eff": de, "Delta_raw": de - float(LAM_P),
        "n_envelopes": n_env, "description": blurb,
        "source": f"{os.path.basename(src)}:{key}" + ("[-1]" if arr.ndim == 3 else ""),
        "fingerprint": fingerprint(a),
    }
    print(f"{name:30s} P={meta[name]['fingerprint']['lattice_power']:8.2f}  "
          f"phi-contrast={meta[name]['fingerprint']['drop_phi_contrast']:6.1f}")

os.makedirs(OUT, exist_ok=True)
np.savez_compressed(os.path.join(OUT, "states.npz"), **states)
doc = {
    "lattice": {
        "model": "zigzag AQH (Haldane-type) coupled-ring lattice, LLE per ring",
        "nx": NX, "ny": NY, "R": int(R), "J": J, "phi": PHI, "spin": "+1",
        "d2": D2, "kex": KEX, "N_modes": N_MODES, "dt": DT,
        "pump_site": int(PUMP), "drop_site": int(DROP),
        "pump_supermode_sigma": int(SIGMA_P), "pump_supermode_lambda": float(LAM_P),
        "pump_supermode_overlap": float(OVERLAP),
        "boundary_sites_ccw": [int(i) for i in ORDER],
        "units": "kappa_in = kappa/2 = 1; time in 2/kappa (lifetimes); "
                 "Delta = (2/kappa)(omega_res - omega_pump), Delta > 0 is RED "
                 "(opposite sign to the Topological Photonics Explorer convention)",
        "equation": "dA_r/dt = [-(1 + i Delta) - i d2 mu^2] A_r,mu - i sum_r' H_rr' A_r',mu "
                    "- kappa_ex,r A_r,mu + i FT[|psi_r|^2 psi_r]_mu + F_r delta_mu,0",
        "drive": "single CW tone on ring 0, F0 = sqrt(F0_squared), placed so that the pumped "
                 "edge supermode sits at Delta_eff (solver Delta = Delta_eff - lambda_sigma)",
    },
    "states": meta,
}
with open(os.path.join(OUT, "states.json"), "w") as f:
    json.dump(doc, f, indent=2)
size = os.path.getsize(os.path.join(OUT, "states.npz")) / 1e6
print(f"\nwrote {len(states)} states to {OUT}/states.npz ({size:.1f} MB) and states.json")
