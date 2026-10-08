"""Extract the catalogued nonlinear states from the (huge) raw sweep/verify runs into one
small portable archive: states.npz (one (R, N) complex frame each) + states.json (parameters
and verification fingerprints).

Provenance only — it needs the raw runs in `lattice_results/lattice_seeds_256modes/runs/`
(several GB, NOT committed). Run it from the repo root:

    python lattice_results/lattice_seeds_256modes/states/build_states.py

Everyone else just uses `states.npz` via `load.py`; this file records exactly where each
state came from.
"""
import json, os, sys

import numpy as np

sys.path.insert(0, ".")
from microring import H_zigzag, boundary_sites, pump_supermode, zigzag_sites

OUT = sys.argv[1] if len(sys.argv) > 1 else "lattice_results/lattice_seeds_256modes/states"
R256 = "lattice_results/lattice_seeds_256modes/runs/"

# --- the lattice these states live on (identical for every entry) -----------------------
NX = NY = 6
J, PHI, D2, KEX, N_MODES, DT = 80.0, float(np.pi / 4), 0.0125, 1.0, 256, 0.0025
H = H_zigzag(NX, NY, J=J, phi=PHI)
R = len(H)
PUMP, DROP = 0, R - (NX - 1)
BND = boundary_sites(H)
LAM_P, OVERLAP, SIGMA_P = pump_supermode(H, PUMP, edge=BND, edge_min=0.6)
xy = zigzag_sites(NX, NY).astype(float)
ORDER = BND[np.argsort(np.arctan2(xy[BND, 1] - xy[:, 1].mean(), xy[BND, 0] - xy[:, 0].mean()))]

# --- catalogue: name -> (source npz, key, F0^2, Delta_eff, n_envelopes, blurb) ----------
SPEC = [
    # P = 2500 — static Turing roll, cold start at the instability onset
    ("turing15_P2500_D0.10", R256 + "P2500/verify_P2500_D+0.10_turing_cold.npz", "trace", 2500, 0.10, 15,
     "Static Turing roll: 15 equidistant phi-pulses in every boundary ring (k = 0), no drift "
     "(0.001 roll periods per lt), 100% of the comb on the k = 15 harmonics; a weak, non-growing "
     "beat with the supermode two ladder rungs away (omega ~ 21, drop-peak pk-pk 3.5%). "
     "Cold start, locks after ~120 lt."),
    # P = 3500 — chaotic comb, and the k = 0 soliton that is only metastable at 256 modes
    ("chaos_P3500_D3.00", R256 + "P3500/verify_P3500_D+3.00_chaos.npz", "trace", 3500, 3.00, 0,
     "Chaotic (MI) comb: lattice power wanders 53-55 without settling, speckled field, filled "
     "fine spectrum. Cold sweep."),
    ("metastable_k0_soliton_P3500_D4.10", R256 + "P3500/verify_P3500_D+4.10_upshift.npz", "trace", 3500, 4.10, 1,
     "METASTABLE, not an attractor: single phi-pulse per boundary ring, uniform along the edge "
     "(k = 0), peak 13.4. Upsampled from the 64-mode soliton1 and alive at t = 150 lt, but it "
     "decays to CW within the following 500 lt. No stable k = 0 soliton was found at 256 modes "
     "for P = 3500-5500, Delta_eff = 3.6-4.6 (lower pump/redder: decays to CW; higher pump or "
     "bluer: breaks into localised travelling envelopes)."),
    # P = 6500 — the nested states
    ("nested_pair_P6500_D6.20", R256 + "P6500/verify_P6500_D+6.20_shelf.npz", "trace", 6500, 6.20, 2,
     "Two-envelope nested travelling state, cold-sweep accessible (the shelf Delta_eff 5.2-6.22 "
     "ends 2 -> 0 at 6.23, with no single-envelope step at 256 modes)."),
    ("nested_soliton_P6500_D5.60", R256 + "P6500/verify_P6500_D+5.60_long.npz", "trace", 6500, 5.60, 1,
     "*** Single nested travelling edge soliton ***: ONE envelope localised along the boundary "
     "and in phi, circulating the 20-ring edge every 0.575 lt (34.8 rings/lt, one chirality), "
     "rigid equidistant comb across all teeth. One envelope in every sample over 650+ lt. NOT "
     "cold-sweep accessible: seeded from the 64-mode state upsampled to 256 modes, which "
     "survives only blue-shifted (the 64-mode point Delta_eff = 6.60 decays to CW at 256 modes)."),
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
        "nx": NX, "ny": NY, "R": int(R), "J": J, "phi": PHI, "spin": "-1",
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
