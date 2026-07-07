# vgsot-sim

A physics-based simulation toolkit for **VGSOT-MTJ (Voltage-Gated Spin-Orbit
Torque Magnetic Tunnel Junction) and SOT-MTJ** switching dynamics, calibrated to
the 80 nm β-W / CoFeB devices.

The toolkit ships **two engines for the same documented macrospin sLLG
physics**, cross-validated against each other trajectory by trajectory: a
**Python** package (`src/vgsot_sim/`) that drives the Monte-Carlo sweeps and the
chapter-figure pipeline, and a **Verilog-A** model (`va/llg/vgsot_llg.va`) that
integrates the magnetisation internally so the same switching event runs inside
an open-source SPICE transient. A result can be reproduced in either a Python
loop or in ngspice.

It couples three-terminal transport, magnetization dynamics, and thermal effects
into a single time-domain switching model, and reproduces the chapter §2.3.3
Sigmoid-shaped `P_sw(V)` curves at the 0.75 ns operating point, the Néel–Brown
`V_th(t_w)` trend, and the §2.3.5 process-variability budget.

**What it models:**

- Three-terminal transport (`V1, V2, V3` → `I_SOT`, `V_MTJ`) with bias-dependent `TMR(V)`
- LLG magnetization dynamics with STT and damping-/field-like SOT, plus VCMA
- Stochastic thermal-fluctuation field (FDT, Brown 1963) for Monte-Carlo switching statistics
- PMA and demagnetization (thin-disk or exact oblate-ellipsoid tensor)
- Coupled self-heating RC dynamics with `M_s(T)`, `K_i(T)` feedback
- Two LLG integrators: legacy spherical-Euler and a norm-preserving Cayley vector stepper

**How to run it:**

- **Python** — a command-line simulator (`vgsot-sim <case>`), a library
  (high-level case API + low-level kernels), and the chapter-figure pipeline
  under `scripts/`. See [Quick start — command line](#quick-start--command-line)
  and [Quick start — Python library](#quick-start--python-library).
- **Verilog-A** — a first-principles time-domain LLG solver for open-source
  SPICE (OpenVAF + ngspice ≥ 43), per-formula aligned to the Python kernels.
  See [Quick start — Verilog-A engine](#quick-start--verilog-a-engine-spice).

Full documentation lives in [`docs/`](docs/README.md).

## Installation

Clone the repository and install in editable mode:

```bash
pip install -e .
```

Requirements (installed automatically): numpy, scipy, matplotlib, tqdm. Tested on
Python 3.10+.

## Quick start — command line

All experiments are exposed through a unified CLI:

```
vgsot-sim <case_name>
```

By default, outputs go to `./result/` — figures as `Chapter02_local_*.png`, time
series and sweeps as `*.csv`. Change the output directory with
`--out_dir my_results`, or disable progress bars with `--no_progress`.

The four cases map to the chapter §2.3.3 Device A P→AP protocol (`t_w = 0.75 ns`
write pulse, `V_MTJ = 0`):

```
vgsot-sim terminal_voltage_control    # three-terminal voltage control (§2.1.1 T-circuit)
vgsot-sim sot_only_constant_current   # baseline: SOT-only, constant current pulse
vgsot-sim sot_switching_no_vcma       # no-VCMA SOT switching: sweep I_SOT, overlay mz(t)
vgsot-sim ser_sot_no_vcma_thermal     # Monte-Carlo P_sw vs |I_SOT| at 0.75 ns
```

See [docs/cases.md](docs/cases.md) for the physical meaning of each case and
[docs/parameters.md](docs/parameters.md) for the default parameters.

## Quick start — Python library

Import and run any case directly, overriding parameters as needed. The workflow
is: build a configuration dataclass → run a case → optionally save with
`result_io`.

Minimal example:

```python
from vgsot_sim import sot_only_constant_current

res = sot_only_constant_current()
print(res.mz[-1])
```

A deterministic SOT sweep at the §2.3.3 operating point:

```python
from vgsot_sim import sot_switching_no_vcma, SotSwitchingNoVcmaConfig

cfg = SotSwitchingNoVcmaConfig(
    # Defaults already use the §2.3.3 protocol; override I_SOT if you like.
    i_sot_list=[-1500e-6, -1200e-6, -1100e-6, -1000e-6, -800e-6],
)
result = sot_switching_no_vcma(cfg)

print(result.mz_curves.keys())
print(result.r_mtj_curves.keys())
```

The returned object:

| field             | type                    | description                |
| ----------------- | ----------------------- | -------------------------- |
| `time_s`          | `np.ndarray`            | simulation time axis       |
| `mz_curves`       | `dict[str, np.ndarray]` | magnetization trajectories |
| `r_mtj_curves`    | `dict[str, np.ndarray]` | MTJ resistance vs time     |
| `pulse_curves`    | `dict[str, np.ndarray]` | applied pulse waveform     |
| `switch_energy_j` | `dict[str, float]`      | switching energy           |
| `pulse_ylabel`    | `str`                   | label for the pulse plot   |

Each entry is keyed by the sweep parameter.

## Monte-Carlo SER

```python
from vgsot_sim import ser_sot_no_vcma_thermal, SerSotNoVcmaThermalConfig

cfg = SerSotNoVcmaThermalConfig(
    trials=500,
    i_sot_list=[-1500e-6, -1200e-6, -1100e-6, -1000e-6, -800e-6],
)
res = ser_sot_no_vcma_thermal(
    cfg,
    seed=2026,                # reproducible per-trial seeding
    enable_self_heating=True, # couple T(t) into M_s, K_i each step
    T_ambient_K=300.0,
)

print(res.x)      # I_SOT values (A)
print(res.ser)    # switching error rate
print(res.psw)    # switching success probability (1 − SER)
```

`SerResult` carries `x`, `ser`, the `psw` property (`1 − ser`, the chapter
convention), and `x_label`. For byte-reproducible runs and integrator selection
(`seed=`, `rng_mode=`, `integrator="cayley"` / `"euler_spherical"`) see
[docs/api.md](docs/api.md).

## Analysis and low-level API

The `vgsot_sim.analysis` sub-package exposes the §2.3 fits used by the chapter
scripts — `nb_fit` (Néel–Brown), `sigmoid_fit` (with Wilson CI), `variability`
(Brinkman D2D budget), and `sampling`. The §2.3.5 process-variability budget is
available in a single call as `vgsot_sim.ser_cases.variability_sweep`.

For full control over the excitation waveform, the low-level kernels
`run_piecewise_terminal_voltage` and `run_piecewise_direct_excitation` accept
per-stage voltages/currents, self-heating and integrator toggles, and an
explicit `rng`. Both are documented in [docs/api.md](docs/api.md).

## Quick start — Verilog-A engine (SPICE)

`va/llg/vgsot_llg.va` is a genuine time-domain LLG solver: it integrates the
magnetisation `m(t)` internally (Cartesian Landau–Lifshitz form, `ddt` state on
internal nodes) rather than replaying a precomputed compact curve. It is the
second engine of the same physics as the Python `switching_vector` stepper —
same effective field, same damping-/field-like SOT and STT torques, same
state-dependent MTJ read `R(m_z)` — so a trajectory can be reproduced in either.

**Toolchain.** An open-source flow of OpenVAF-Reloaded (compiles the `.va` to an
OSDI module) and ngspice ≥ 43 (loads the OSDI model and runs the transient). The
cross-validation script finds the executables from the `OPENVAF` / `NGSPICE`
environment variables, then from `PATH`.

Build the model and run the bundled deterministic switching testbench:

```sh
cd va/llg
openvaf vgsot_llg.va -o vgsot_llg.osdi     # compile to OSDI
ngspice_con -b tb_switch.spice             # .spiceinit loads the OSDI model
```

`tb_switch.spice` forces a super-threshold `I_SOT` (as a voltage across the
Ohmic SOT channel), holds `V_MTJ = 0`, starts `m` near `+z`, and writes `m(t)`
to `vgsot_llg_out.csv`. The magnetisation components `mx, my, mz` are exposed as
terminals, so the testbench sets the initial state with `.ic` and probes them
with `v(mx)` etc.; the transient runs with `uic`.

**Stochasticity (harness-owned RNG).** OpenVAF/OSDI does not reliably support
in-module random state, so the FDT thermal field is not drawn inside the module.
The input nodes `hx, hy, hz` instead carry the thermal field `H_th` [A/m]; the
Python harness drives them with the correctly scaled Brown-1963 white-noise
stream (same seed → identical noise in both engines). Grounding `hx/hy/hz` gives
the deterministic limit.

**Cross-validation.** `va/llg/cross_validate.py` drives the same deterministic
SOT setup through both engines and asserts their `m_z(t)` trajectories agree:

```sh
python va/llg/cross_validate.py
```

It compiles the model, runs the ngspice transient, compares against the Python
`switching_vector` reference, and writes `va/llg/cross_validate_summary.json`. If
OpenVAF or ngspice are absent the check skips gracefully (exit 0). The committed
run agrees to `max |Δm_z| ≈ 0.008` over 0–3 ns with `|m| = 1` conserved to
~1e-4 (exact at the current-on equilibrium); the small residual is fixed-step
Cayley versus ngspice's adaptive `ddt`.

See [`va/README.md`](va/README.md) for the model's physics, parameters, and
validation status.

## Notes

- Units: `I_SOT` in **Ampere**, `V_MTJ` in **Volt**, time in **seconds**.
- `pap=1` initialises in the anti-parallel state (`mz ≈ −1`); the resistance is
  flat at `R_AP` (not `R_P`) at the start of the run. The initialisation-angle
  and resistance conventions cancel within any single simulation.

## Documentation

- [Documentation index](docs/README.md) — start here
- [Simulation cases](docs/cases.md) · [Default parameters](docs/parameters.md) · [Parameter validation](docs/parameter_validation.md)
- [API reference](docs/api.md) · [Project structure](docs/structure.md) · [Technical details](docs/technical_details.md)
- [Reproducing the Chapter 2 figures and tables](docs/reproducing_figures.md) — regenerate every figure and table from its script
- [Verilog-A engine](va/README.md) — the second (SPICE) engine: build, run, and cross-validation against the Python solver
- Development: [implementation status](docs/maintenance/IMPLEMENTATION_STATUS.md) · [version notes](docs/maintenance/version_notes.md) · [claims backing matrix](CLAIMS.md)

## Citation

If you use this simulator in academic research, please cite:

```
Zhang Jincheng. (2026). VGSOT-SIM: A VGSOT switching simulation toolkit [Computer software]. GitHub.
```
