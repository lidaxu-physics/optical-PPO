"""
Checks on a catalogued state of the zigzag 6x6 lattice before it is used as a feature map (run from the repository root):

    python lattice_results/aqh_66zigzag_results/state_checks.py stats     --state linear_P500_D1.00 --eps 1.3 [--periods 400]
    python lattice_results/aqh_66zigzag_results/state_checks.py tolerance --state soliton1_P3500_D4.90 --eps 0.02 0.05 0.1 0.2 0.5
    python lattice_results/aqh_66zigzag_results/state_checks.py relax     --state turing14_P2500_D0.10 --eps 0.5

stats     : pump + tones at a held input (0). After the relaxation, the fine-line coefficients c_n of the drop ring are
            computed period by period for `periods` periods: time series of |c_n| and arg c_n, histograms with a Gaussian
            fit (mean, sigma, CV, skewness), the trajectory in the complex plane, the autocorrelation time, and the gap
            between the time averages of two lanes started from different noise (ergodicity). A periodic response gives
            delta-function histograms; a chaotic one, their width and tau_c set T_avg and whether the mean field is usable.
            Figure: <out>/stats_<state>_eps<eps>.png
tolerance : for every eps, hold the inputs at their extremes (+1, -1, mixed) for T_hold lifetimes and report whether the
            state survived (its fingerprint: lattice power, phi contrast, boundary-envelope contrast) and the features'
            contrast between the extremes.
relax     : step response -- switch the input and record the per-period c_n until they settle; the decay rate sets T_relax.
"""
import argparse, os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import numpy as np, torch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from microring import make_ring, TASKS
from microring import plot_style as ps
from lattice_results.lattice_seeds.states.load import diagnostics

ps.apply()
p = argparse.ArgumentParser()
p.add_argument("what", choices=["stats", "tolerance", "relax"])
p.add_argument("--state", required=True)
p.add_argument("--eps", type=float, nargs="+", required=True)
p.add_argument("--periods", type=int, default=400)
p.add_argument("--T_hold", type=float, default=100.0)
p.add_argument("--T_relax", type=float, default=10.0)
p.add_argument("--T_warmup", type=float, default=20.0)
p.add_argument("--average", choices=["field", "power"], default="field")
p.add_argument("--tone_rate", type=float, default=None)
p.add_argument("--N", type=int, default=None)
p.add_argument("--threads", type=int, default=6)
p.add_argument("--out", default=os.path.dirname(os.path.abspath(__file__)))
args = p.parse_args()
torch.set_num_threads(args.threads)
np.set_printoptions(precision=4, suppress=True, linewidth=200)
task = TASKS["Pendulum-v1"]


def build(B, eps, seed=0):
    return make_ring("topo", B, task["obs_scale"], seed=seed, squash=task["squash"], init_state=args.state, eps=eps, N=args.N,
                     T_warmup=args.T_warmup, T_relax=args.T_relax, T_avg=1.2, average=args.average, tone_rate=args.tone_rate,
                     dtype=torch.complex128)


def per_period(fm, n_periods):
    """(n_periods, B, lines) complex: the Fourier coefficient of every fine line over each successive period."""
    sol, lines, out = fm.solver, fm.fine_lines, []
    for _ in range(n_periods):
        c = torch.zeros(fm.B, len(lines), dtype=torch.complex128)
        for _ in range(fm.n_period):
            fm.a, _ = sol.evolve(fm.a, fm.sample_every)
            c += fm.a[:, fm.readout_sites[0], 0, None].to(torch.complex128) * torch.as_tensor(np.exp(1j * fm.fine_delta * lines * sol.t))
        out.append((c / fm.n_period).numpy())
    return np.array(out)


def hold(fm, obs, T):
    fm.solver.set_drive(fm.F0, fm.encode(obs))
    fm.a, _ = fm.solver.evolve(fm.a, int(round(T / fm.solver.dt)))


if args.what == "stats":
    eps = args.eps[0]
    t0 = time.time()
    fm, cfg = build(2, eps)                                       # lane 1 starts from the same frame plus noise: a second realisation
    fm.a[1] += 1e-2 * fm.solver.random_state(1, generator=torch.Generator().manual_seed(1))[0].to(fm.a.dtype)
    hold(fm, np.zeros((2, 3)), args.T_relax)
    C = per_period(fm, args.periods)                              # (periods, 2, lines)
    lines = list(fm.fine_lines)
    print(f"{args.state} eps={eps}: {cfg['nx']}x{cfg['ny']} {cfg['lattice']} N={cfg['N']} pump sigma {cfg['pump_sigma']} tones {cfg['tone_sigma']} "
          f"delta {cfg['fine']['delta']:.4f} period {fm.period:.4f} lt; {args.periods} periods = {args.periods * fm.period:.0f} lt; {time.time() - t0:.0f} s")
    A, PH = np.abs(C[:, 0]), np.angle(C[:, 0])
    fig, ax = plt.subplots(3, len(lines), figsize=(3.6 * len(lines), 9.0), squeeze=False)
    t = np.arange(args.periods) * fm.period
    rows = []
    for k, n in enumerate(lines):
        a, ph = A[:, k], np.unwrap(PH[:, k])
        mu, sd = a.mean(), a.std()
        cv = sd / mu if mu > 0 else np.nan
        skew = ((a - mu) ** 3).mean() / sd ** 3 if sd > 0 else 0.0
        c_mean = C[:, 0, k].mean()
        ph_sd = np.abs(C[:, 0, k] - c_mean).std() / max(abs(c_mean), 1e-300)
        # autocorrelation time of |c_n|: integral of the normalised autocorrelation
        x = a - mu
        if sd > 0:
            ac = np.correlate(x, x, "full")[len(x) - 1:] / (x @ x)
            first_neg = np.argmax(ac < 0.05) if (ac < 0.05).any() else len(ac)
            tau = max(ac[:first_neg].sum() - 0.5, 0.5) * fm.period
        else:
            tau = 0.0
        # ergodicity: time averages of the two lanes vs the length of the average
        gap = [abs(np.abs(C[:m, 0, k]).mean() - np.abs(C[:m, 1, k]).mean()) / mu for m in (args.periods // 8, args.periods // 2, args.periods)]
        rows.append((n, mu, cv, skew, ph_sd, tau, gap))
        ax[0, k].plot(t, a, color=ps.BLUE, lw=0.7); ax[0, k].axhline(mu, color=ps.RED, lw=1, ls="--")
        ax[0, k].set_title(f"line n = {n:+d}:  |c|   mean {mu:.3g}, CV {cv:.2e}", fontsize=9); ax[0, k].set_xlabel("time (lifetimes)", fontsize=8)
        ax[1, k].plot(t, ph, color=ps.ORANGE, lw=0.7); ax[1, k].set_title(f"arg c (unwrapped), sd {np.std(ph):.2e} rad", fontsize=9); ax[1, k].set_xlabel("time (lifetimes)", fontsize=8)
        if sd > 0 and cv > 1e-6:
            ax[2, k].hist(a, bins=40, density=True, color=ps.BLUE, alpha=0.5)
            xs = np.linspace(a.min(), a.max(), 200); ax[2, k].plot(xs, np.exp(-0.5 * ((xs - mu) / sd) ** 2) / (sd * np.sqrt(2 * np.pi)), color=ps.RED, lw=1.2)
            ax[2, k].set_title(f"|c|: Gaussian fit, skew {skew:+.2f}, tau_c {tau:.2f} lt", fontsize=9)
        else:
            ax[2, k].text(0.5, 0.5, f"stationary\n(CV {cv:.1e})", ha="center", va="center", transform=ax[2, k].transAxes, fontsize=10)
            ax[2, k].set_xticks([]); ax[2, k].set_yticks([])
        for a_ in ax[:, k]:
            a_.tick_params(labelsize=7)
    fig.suptitle(f"{args.state}, eps = {eps}: statistics of the fine-line coefficients of the drop ring, period by period ({args.periods} periods)", fontsize=10, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = os.path.join(args.out, f"stats_{args.state}_eps{eps:g}.png")
    fig.savefig(out, dpi=130)
    print(f"{'line':>5s} {'mean |c|':>10s} {'CV':>9s} {'skew':>7s} {'|dc|/|c|':>9s} {'tau_c':>7s}  ergodicity gap at P/8, P/2, P")
    for n, mu, cv, skew, ph_sd, tau, gap in rows:
        print(f"{n:+5d} {mu:10.4g} {cv:9.2e} {skew:+7.2f} {ph_sd:9.2e} {tau:7.2f}  {np.array(gap)}")
    print("saved", out)

elif args.what == "tolerance":
    ext = np.array([[1.0, 1.0, 1.0], [-1.0, -1.0, -1.0], [1.0, -1.0, 1.0], [-1.0, 1.0, -1.0]])
    for eps in args.eps:
        t0 = time.time()
        fm, cfg = build(len(ext), eps)
        ref = diagnostics(fm.a[:1].cpu())
        hold(fm, fm.unsquashed(ext), args.T_hold)
        C = per_period(fm, 2)[-1]                                  # (B, lines) after the hold
        fp = [diagnostics(fm.a[i:i + 1].cpu()) for i in range(len(ext))]
        P = (fm.a.abs() ** 2).sum((1, 2)).numpy()
        alive = [abs(f["lattice_power"] / ref["lattice_power"] - 1) < 0.08 and abs(f["drop_phi_contrast"] / ref["drop_phi_contrast"] - 1) < 0.3 for f in fp]
        contrast = np.ptp(np.abs(C) ** 2, axis=0) / (np.abs(C) ** 2).mean(0)
        print(f"eps = {eps:<6g}: after {args.T_hold:g} lt at the extremes -- lattice power {P} (warm-up {ref['lattice_power']:.3f}), "
              f"phi contrast {[round(f['drop_phi_contrast'], 2) for f in fp]} (warm-up {ref['drop_phi_contrast']:.2f}), "
              f"state alive in {sum(alive)}/{len(ext)} lanes | feature contrast across the extremes per line {contrast} | {time.time() - t0:.0f} s")

else:
    eps = args.eps[0]
    fm, cfg = build(2, eps)
    hold(fm, np.zeros((2, 3)), 3 * args.T_relax)
    fm.solver.set_drive(fm.F0, fm.encode(fm.unsquashed(np.array([[0.7, -0.5, 0.6], [-0.7, 0.5, -0.6]]))))
    C = per_period(fm, int(round(40 / fm.period)))             # 40 lifetimes after the step
    A = np.abs(C) ** 2
    final = A[-5:].mean(0)
    dev = np.abs(A - final).max(axis=(1, 2)) / np.abs(final).max()
    t = (np.arange(len(dev)) + 1) * fm.period
    ok = (dev > 1e-7) & (dev < 0.5 * dev[0]) & (t > 2)
    rate = -np.polyfit(t[ok], np.log(dev[ok]), 1)[0] if ok.sum() > 4 else float("nan")
    print(f"{args.state} eps={eps}: residual of the previous input, relative to the largest line, after")
    for T in (1, 2, 3, 5, 10, 15, 20, 30):
        i = int(round(T / fm.period)) - 1
        if i < len(dev):
            print(f"   {T:3d} lt: {dev[i]:.2e}")
    print(f"   fitted relaxation rate {rate:.3f} per lifetime (T_relax for 1e-3: {np.log(1e3) / rate:.1f} lt)")
