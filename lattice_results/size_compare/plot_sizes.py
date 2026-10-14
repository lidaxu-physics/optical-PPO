"""
Pendulum below threshold versus the size of the zigzag lattice: learning curves, frozen evaluation, and what the size
changes physically (mini FSR, intensity per boundary ring). Reads <protocol>/zz<n>/mr_topo_zz<n>_seed*.json for the
protocols fixed_pump and fixed_power, plus the 6 x 6 point of both from ../aqh_66zigzag_results/pattern_free/
mr_topo_pendulum_seed*.json. Sizes still running are left out. Boundary intensities are computed once (the pump alone,
warmed up) and cached in sizes_meta.json.

    python lattice_results/size_compare/plot_sizes.py        (from the repository root)
"""
import glob, json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", ".."))
from microring import plot_style as ps

ps.apply()
SIZES = [4, 6, 8, 10, 12]
PROTOCOLS = [("fixed_pump", "fixed pump, F0² = 500", "-", "o"), ("fixed_power", "fixed intensity per boundary ring", "--", "s")]
REF = [("linear policy, no ring", -1038, -718), ("best single ring", -252, -196)]
YLIM = (-1400, 0)


def load_runs(protocol, n):
    if n == 6:
        files = sorted(glob.glob(os.path.join(HERE, "..", "aqh_66zigzag_results", "pattern_free", "mr_topo_pendulum_seed*.json")))
    else:
        files = sorted(glob.glob(os.path.join(HERE, protocol, f"zz{n}", f"mr_topo_zz{n}_seed*.json")))
    return [json.load(open(f)) for f in files]


def boundary_intensity(n, F2):
    """Mean |a|^2 over the boundary rings of the pumped n x n zigzag lattice (pump alone, warmed up); as in run_sizes.py."""
    from microring import make_ring, TASKS
    from microring.lattice import boundary_sites
    task = TASKS["Pendulum-v1"]
    fm, _ = make_ring("topo", 1, task["obs_scale"], seed=0, squash=task["squash"], lattice="zigzag", nx=n, ny=n, J=80.0,
                      F0=float(np.sqrt(F2)), target_Delta=1.0, N=1, dt=0.0025, eps=0.0, T_relax=1.0, T_avg=1.2, T_warmup=20.0)
    I = (fm.a[0, :, 0].abs() ** 2).cpu().numpy()
    return float(I[boundary_sites(fm.solver.H)].mean())


runs = {p: {n: load_runs(p, n) for n in SIZES} for p, *_ in PROTOCOLS}
present = {p: [n for n in SIZES if runs[p][n]] for p in runs}
if not any(present.values()):
    sys.exit("no finished runs yet")

# the physics per size: mini FSR from the result files, boundary intensity computed and cached
meta_path = os.path.join(HERE, "sizes_meta.json")
meta = json.load(open(meta_path)) if os.path.exists(meta_path) else {}
for p in runs:
    for n in present[p]:
        r = runs[p][n][0]["ring"]; key = f"{p}_{n}"
        if key not in meta:
            meta[key] = dict(F2=float(r["F0"]) ** 2, delta=float(r["fine"]["delta"]), R=int(r["nx"]) * (int(r["ny"]) - 1) * 2,
                             I_boundary=boundary_intensity(n, float(r["F0"]) ** 2))
json.dump(meta, open(meta_path, "w"), indent=1)

cols = dict(zip(SIZES, [ps.NAVY, ps.BLUE, ps.AQUA, ps.ORANGE, ps.RED]))
fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4), gridspec_kw=dict(width_ratios=[1.35, 1, 0.9]))

ax = axes[0]                                                             # learning curves: mean over the seeds
for name, lo, hi in REF:
    ax.axhspan(lo, hi, color=ps.GRID, alpha=0.5, lw=0, zorder=0)
for p, label, ls, mk in PROTOCOLS:
    for n in present[p]:
        if n == 6 and p != "fixed_pump":
            continue                                                     # the 6 x 6 point is shared, drawn once
        curves = []
        for r in runs[p][n]:
            x = np.array([u["env_steps"] for u in r["updates"]]) / 1e3
            y = np.array([u["mean_return"] for u in r["updates"]], dtype=float)
            grid = np.arange(5, 206, 5)
            ok = ~np.isnan(y)
            curves.append(np.interp(grid, x[ok], y[ok]))
        m = np.mean(curves, 0)
        ax.plot(grid, m, color=cols[n], ls=ls, lw=1.4, label=f"{n}×{n}" + ("" if p == "fixed_pump" else " (fixed intensity)"))
ax.set_xlabel("environment steps (thousands)", fontsize=9); ax.set_ylabel("mean return per update, mean of 3 seeds", fontsize=9)
ax.set_ylim(*YLIM); ax.set_xlim(0, 206); ax.tick_params(labelsize=8)
ax.legend(fontsize=7.5, loc="lower right", frameon=False, ncol=2)
ax.set_title("learning curves", fontsize=10, loc="left")

ax = axes[1]                                                             # frozen evaluation versus size
for name, lo, hi in REF:
    ax.axhspan(lo, hi, color=ps.GRID, alpha=0.5, lw=0, zorder=0)
    ax.text(12.6, (hi + lo) / 2, name, fontsize=7, ha="right", va="center", color=ps.INK2)
for p, label, ls, mk in PROTOCOLS:
    xs, ms = [], []
    for n in present[p]:
        ev = [r["eval"]["mean"] for r in runs[p][n] if "eval" in r]
        if not ev:
            continue
        if n != 6 or p == "fixed_pump":                                   # the shared 6 x 6 point: drawn once
            ax.plot([n] * len(ev), ev, marker="o", ms=5, ls="", color=cols[n], alpha=0.7)
        xs.append(n); ms.append(np.mean(ev))
    if len(xs) > 1 or p == "fixed_pump":
        ax.plot(xs, ms, ls=ls, color=ps.INK2, lw=1.2, marker=mk, ms=7, mfc="white", label=label)
ax.set_xticks(SIZES); ax.set_xticklabels([f"{n}×{n}" for n in SIZES], fontsize=8)
ax.set_xlim(3, 13); ax.set_ylim(*YLIM); ax.tick_params(labelsize=8)
ax.set_ylabel("frozen evaluation (16 episodes)", fontsize=9)
ax.legend(fontsize=7.5, loc="lower left", frameon=False)
ax.set_title("frozen evaluation: seeds and their mean", fontsize=10, loc="left")

ax = axes[2]                                                             # what the size changes
ns = present["fixed_pump"]
ax.plot(ns, [meta[f"fixed_pump_{n}"]["I_boundary"] for n in ns], ls="-", marker="o", ms=6, color=ps.ORANGE, lw=1.2, label="fixed pump: |a|² per boundary ring")
if 6 in ns:
    ax.axhline(meta["fixed_pump_6"]["I_boundary"], ls="--", color=ps.ORANGE, lw=1.0, label="fixed intensity: held at the 6×6 value")
ax.set_ylabel("mean |a|² over the boundary rings", fontsize=9, color=ps.ORANGE); ax.set_ylim(0, None)
ax.set_xticks(SIZES); ax.set_xticklabels([f"{n}×{n}" for n in SIZES], fontsize=8); ax.set_xlim(3, 13); ax.tick_params(labelsize=8)
ax2 = ax.twinx()
ns = present["fixed_pump"]
ax2.plot(ns, [meta[f"fixed_pump_{n}"]["delta"] for n in ns], color=ps.NAVY, marker="D", ms=5, lw=1.2, label="mini FSR δ")
ax2.set_ylabel("mini FSR δ (units of κ/2)", fontsize=9, color=ps.NAVY); ax2.tick_params(labelsize=8); ax2.set_ylim(0, 22)
ax2.grid(False)
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, fontsize=7, loc="lower left", frameon=False)
ax.set_title("what the size changes", fontsize=10, loc="left")

done = {p: present[p] for p in present}
fig.suptitle("Pendulum below threshold (N = 1, J = 80, Δ_eff = 1) on zigzag lattices of different size; 16 envs × 32 steps × 400 updates, 3 seeds",
             fontsize=10.5, x=0.01, ha="left")
fig.tight_layout(rect=(0, 0, 1, 0.94))
out = os.path.join(HERE, os.environ.get("SIZES_OUT", "sizes_pendulum.png"))
fig.savefig(out, dpi=150)
print("saved", out, "| finished:", done)
