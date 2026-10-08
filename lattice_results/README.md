# lattice_results

Everything of the coupled-ring (topological) lattice part of this repository.

| folder / file | what |
|---|---|
| `LATTICE_EXTENSION.md` | the write-up: model, mini-comb drive and read-out, validation, results, the three questions |
| `lattice_seeds/` | the archive of nonlinear states of the zigzag 6 × 6 lattice (`states/`: load any state and keep computing), their catalogue (`STATES.md`), the sweep / verification tools, and `summary/` — one movie per kind of state (below threshold, static soliton, nested soliton, chaotic comb, soliton crystal, Turing roll) |
| `aqh_44_results/` | PPO runs on the AQH 4 × 4 mini-comb (Pendulum): result files, animations; checkpoints are not tracked |
| `aqh_66zigzag_results/` | PPO runs on the zigzag 6 × 6 lattice, one folder per kind of state the lattice is operated in: `pattern_free/` (below threshold; LunarLander so far). Planned: `chaos/`, `normal/`, `rolls/`, `static_soliton/`, `nested_soliton/`, seeded from `lattice_seeds/` |
| `size_compare/` | the same task on lattices of different size (planned) |
| `characterization/` | `05` the mini-comb drive on the drop spectrum, `06` which lines to read and how many modes to simulate, `07` the optical-power controls |
| `animate_mini.py` | episode animation of a trained lattice policy |

Training writes into one of these folders with `python PPO_MR.py ... --out_dir lattice_results/<folder>`;
all scripts here are run from the repository root.
