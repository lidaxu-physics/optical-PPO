# size_compare: does the lattice size matter?

Pendulum below threshold (the pattern-free state) on zigzag AQH lattices of 4 × 4 to 12 × 12 rings, everything else as
in `../aqh_66zigzag_results/pattern_free/`: J = 80, effective detuning 1.0, N = 1 longitudinal mode (below threshold
only μ = 0 carries light), dt = 0.0025, ε = 1.3, T_relax = 12, T_avg = 1.2, the pump on the edge supermode the input
corner couples to best and one tone per observation on the three neighbouring edge supermodes, four fine lines read at
the drop ring. Budget 16 envs × 32 steps × 400 updates (205 k environment steps), 3 seeds, frozen evaluation over 16
greedy episodes. The 6 × 6 point is the pattern_free run itself (its archive state `linear_P500_D1.00` is exactly this
configuration: same pump σ = 29, tones [28, 30, 31], δ = 10.516).

Two protocols, because the same pump spread over a longer boundary gives each ring less light:

| protocol | pump | what is held fixed |
|---|---|---|
| `fixed_pump/` | F0² = 500 on every lattice | the drive; the intensity per boundary ring falls from 1.65 (4 × 4) to 0.11 (12 × 12) |
| `fixed_power/` | F0² = 188, 500, 834, 1291, 1832 for 4, 6, 8, 10, 12 | the mean |a|² over the boundary rings (0.655, the 6 × 6 value): the same Kerr shift per ring, more rings |

![sizes](sizes_pendulum.png)

## Result

Frozen evaluation, mean ± std over the three seeds (single ring −196 … −252, linear policy −718 … −1038):

| | 4 × 4 | 6 × 6 | 8 × 8 | 10 × 10 | 12 × 12 |
|---|---|---|---|---|---|
| fixed pump | −467 ± 85 | −347 ± 71 | −441 ± 115 | −678 ± 226 | −713 ± 230 |
| fixed intensity | −350 ± 49 | −347 ± 71 | −293 ± 101 | −370 ± 97 | −293 ± 53 |

- **At fixed intensity per ring the size does not matter**: all five sizes land at −290 … −370, within the seed
  scatter, and the learning curves of 8 × 8 to 12 × 12 lie on top of each other (they are even the fastest, because
  their pump is stronger and the tones with it).
- **At fixed pump the performance falls with the size** (−347 at 6 × 6 → −713 at 12 × 12): that is the intensity per
  ring falling, not the lattice getting bigger. 4 × 4 at F0² = 500 (intensity 1.65, 2.5 × the 6 × 6 value) is also worse
  than at the matched 188: too much Kerr shift is not better either.
- What the size does change (right panel): the mini FSR δ shrinks from 18.2 (4 × 4) to 4.7 (12 × 12) as more edge
  supermodes share the band, so the fine comb gets denser and one period 2π/δ longer (0.35 → 1.35 lifetimes); the
  read-out window was rounded to whole periods in every case.

So, for the pattern-free state, the policy quality is set by the nonlinearity per ring, and a larger lattice is
neither better nor worse once that is matched. Whether a larger lattice helps a task that needs more features (more
edge lines to read) is not tested here: all runs read the same four lines.

## Files

- `run_sizes.py` launches a protocol (`--protocol fixed_pump|fixed_power --sizes 4 8 10 12 --seeds 0 1 2`, `--dry` to print the commands, `--resume` to continue); it computes the fixed-intensity pumps itself.
- `plot_sizes.py` draws `sizes_pendulum.png` from whatever has finished; `sizes_meta.json` caches the boundary intensities and δ per size.
- `<protocol>/zz<n>/`: `mr_topo_zz<n>_seed<s>.json` (result, config, learning curve, evaluation), `command.txt` (the state arguments), `checkpoints/` (not tracked).
