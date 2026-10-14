# lattice_results

Everything of the coupled-ring (topological) lattice part of this repository.

| folder / file | what |
|---|---|
| `LATTICE_EXTENSION.md` | the write-up: model, mini-comb drive and read-out, validation, the first results, the five nonlinear states as Pendulum policies, the size scan, and where this leaves the lattice against the single ring |
| `lattice_seeds/` | the archive of nonlinear states of the zigzag 6 × 6 lattice at 64 longitudinal modes (`states/`: load any state and keep computing), their catalogue (`STATES.md`), the sweep / verification tools, and `summary/` — one movie per kind of state |
| `lattice_seeds_256modes/` | the same hunt at 256 modes: the Turing roll and the chaotic comb survive, the 64-mode static soliton is metastable and the nested soliton moves to Δ_eff = 5.60 (see its README) |
| `aqh_44_results/` | the first PPO runs on the AQH 4 × 4 mini-comb (Pendulum, one seed each): result files, animations |
| `aqh_66zigzag_results/` | the zigzag 6 × 6 lattice, one folder per state it is operated in, three seeds each: `pattern_free/` (below threshold; also the LunarLander runs and the all-lines Pendulum runs), `rolls/`, `static_soliton/`, `nested_soliton/`, `chaos/`; `baselines/` the budget-matched linear and poly2 policies without a ring; `plot_states.py` → `states_pendulum.png`; `state_checks.py` and `stats_*.png` the checks made before training; each state folder holds its `command.txt` and `topo_pendulum.mp4` |
| `size_compare/` | Pendulum below threshold on 4 × 4 … 12 × 12 zigzag lattices, fixed pump and fixed intensity per ring: `run_sizes.py`, `plot_sizes.py` → `sizes_pendulum.png`, its README, one animation per size |
| `characterization/` | `05` the mini-comb drive on the drop spectrum, `06` which lines to read and how many modes to simulate, `07` the optical-power controls |
| `animate_mini.py` | episode animation of a trained lattice policy (`--dir <run folder> --seed <s> [--tag <tag>] --out <file>.mp4`) |

Training writes into one of these folders with `python PPO_MR.py ... --out_dir lattice_results/<folder>`;
checkpoints (`<folder>/checkpoints/`) are not tracked; all scripts here are run from the repository root.
The logs of the runs live in `notes/` (not tracked).
