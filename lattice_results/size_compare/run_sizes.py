"""
Pendulum below threshold on zigzag AQH lattices of different size: launch the training runs.

Two protocols, both at J = 80, effective detuning 1.0, N = 1 (below threshold only mu = 0 carries light), dt = 0.0025,
eps = 1.3, T_relax = 12, T_avg = 1.2 -- the settings of aqh_66zigzag_results/pattern_free, whose 6 x 6 runs are the
6 x 6 point of both protocols (its archive state linear_P500_D1.00 is exactly this configuration):

    fixed_pump   the same pump F0^2 = 500 on every lattice. The light spreads over a longer boundary, so the
                 intensity per boundary ring, and with it the Kerr shift, falls with the size.
    fixed_power  F0^2 scaled so that the mean |a|^2 over the boundary rings equals the 6 x 6 one (0.664 at
                 F0^2 = 500): the same nonlinearity per ring, more rings.

    python lattice_results/size_compare/run_sizes.py --protocol fixed_pump --sizes 4 8 10 12 --seeds 0 1 2
    python lattice_results/size_compare/run_sizes.py --protocol fixed_power --sizes 4 8 10 12 --seeds 0 1 2

Runs go to lattice_results/size_compare/<protocol>/zz<n>/ (result mr_topo_zz<n>_seed<s>.json, checkpoints/), logs to
notes/size_<protocol>_zz<n>_s<s>.log. --dry prints the commands only. Run from the repository root.
"""
import argparse, os, subprocess, sys
import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
sys.path.insert(0, ROOT)
P0, REF_SIZE = 500.0, 6                   # the 6 x 6 point: F0^2 = 500, the archive's linear_P500_D1.00
COMMON = "--env Pendulum-v1 --regime topo --lattice zigzag --J 80 --target_Delta 1.0 --N 1 --dt 0.0025 --eps 1.3 --T_relax 12 --T_avg 1.2 --T_warmup 20"


def boundary_intensity(n, F2):
    """Mean |a|^2 over the boundary rings of the pumped n x n zigzag lattice (pump alone, warmed up)."""
    from microring import make_ring, TASKS
    from microring.lattice import boundary_sites
    task = TASKS["Pendulum-v1"]
    fm, _ = make_ring("topo", 1, task["obs_scale"], seed=0, squash=task["squash"], lattice="zigzag", nx=n, ny=n, J=80.0,
                      F0=float(np.sqrt(F2)), target_Delta=1.0, N=1, dt=0.0025, eps=0.0, T_relax=1.0, T_avg=1.2, T_warmup=20.0)
    I = (fm.a[0, :, 0].abs() ** 2).cpu().numpy()
    return float(I[boundary_sites(fm.solver.H)].mean())


p = argparse.ArgumentParser()
p.add_argument("--protocol", choices=["fixed_pump", "fixed_power"], required=True)
p.add_argument("--sizes", type=int, nargs="+", default=[4, 8, 10, 12])
p.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
p.add_argument("--n_envs", type=int, default=16)
p.add_argument("--n_steps", type=int, default=32)
p.add_argument("--n_updates", type=int, default=400)
p.add_argument("--threads", type=int, default=1)
p.add_argument("--resume", action="store_true")
p.add_argument("--dry", action="store_true")
args = p.parse_args()

env = dict(os.environ, PYTHONIOENCODING="utf-8", PYTHONUTF8="1", PYTHONDONTWRITEBYTECODE="1")
I_ref = boundary_intensity(REF_SIZE, P0) if args.protocol == "fixed_power" else None
for n in args.sizes:
    if args.protocol == "fixed_pump":
        F2 = P0
    else:                                   # below threshold |a|^2 is nearly linear in F0^2: one correction step is enough
        F2 = P0 * I_ref / boundary_intensity(n, P0)
        F2 *= I_ref / boundary_intensity(n, F2)
    out = f"lattice_results/size_compare/{args.protocol}/zz{n}"
    os.makedirs(os.path.join(ROOT, out, "checkpoints"), exist_ok=True)
    state = f"--nx {n} --ny {n} --F0 {np.sqrt(F2):.6f} --tag zz{n}"
    open(os.path.join(ROOT, out, "command.txt"), "w").write(f"{COMMON} {state}\n")
    print(f"{n}x{n}: F0^2 = {F2:.1f}" + (f" (boundary |a|^2 matched to {I_ref:.3f})" if I_ref else ""))
    for s in args.seeds:
        cmd = (f"{sys.executable} PPO_MR.py {COMMON} {state} --seed {s} --n_envs {args.n_envs} --n_steps {args.n_steps} "
               f"--n_updates {args.n_updates} --eval_episodes 16 --checkpoint_every 10 --threads {args.threads} --out_dir {out}"
               + (" --resume" if args.resume else ""))
        log = os.path.join(ROOT, "notes", f"size_{args.protocol}_zz{n}_s{s}")
        print("  ", cmd if args.dry else f"seed {s} -> {log}.log")
        if not args.dry:
            flags = subprocess.CREATE_NEW_PROCESS_GROUP | getattr(subprocess, "DETACHED_PROCESS", 0) if os.name == "nt" else 0
            subprocess.Popen(cmd.split(), cwd=ROOT, env=env, stdout=open(log + ".log", "w"), stderr=open(log + ".err", "w"),
                             creationflags=flags, start_new_session=os.name != "nt")
