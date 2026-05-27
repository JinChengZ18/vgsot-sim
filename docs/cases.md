

---

## Simulation Cases and Their Physical Meaning

Below is a full description of each available case and what physical mechanism
it is designed to demonstrate. All cases share a common physics core (electronic
mapping → LLG step → TMR feedback) and differ only in which terminal
parameters they sweep, whether VCMA / thermal noise are enabled, and whether
they aggregate trials into a SER curve.

### Common toggles (apply to every case below)

The piecewise drivers `run_piecewise_direct_excitation` and
`run_piecewise_terminal_voltage` (which all cases ultimately call) accept the
following opt-in toggles. Defaults preserve the original chapter behaviour.

| Toggle | Default | Effect |
|---|---|---|
| `enable_self_heating` | `False` | Advance `T(t)` through the RC thermal network, recompute `M_s(T)` / `K_i(T)` each step, and feed them back into `anisotropy.field()`. Diagnostics populated in `SimResult.T_K`, `Ms_T`, `Ki_T`. |
| `T_ambient_K` | `300.0` | Initial / sink temperature for the thermal network. |
| `integrator` | `"euler_spherical"` | LLG step type. `"cayley"` selects a norm-preserving Cartesian step that takes an explicit `sigma_SH` 3-vector. |
| `sigma_SH` | `None` (= `[-1, 0, 0]`) | Spin-Hall polarisation direction; honoured only by the Cayley stepper. |
| `rng` | `None` | Pass an `np.random.Generator` to make a single run byte-reproducible end-to-end (init draw, thermal noise, every kernel). |
| `demag_mode` | `"ellipsoid"` | `"ellipsoid"` (exact oblate tensor with `D_elec`) or `"thin_disk"` (legacy approximation). |

The Monte-Carlo SER case (`ser_sot_no_vcma_thermal`) adds two more:

| Toggle | Default | Effect |
|---|---|---|
| `seed` | `None` | Master seed; per-trial seeds are derived as `_trial_seed(seed, i_sot, trial_idx)`. |
| `rng_mode` | `"legacy"` | `"legacy"` calls `np.random.seed(...)` once per trial (chapter-figure bit-stream); `"generator"` uses `np.random.default_rng(...)` per trial (multiprocessing-friendly, no global state). |

------

## `terminal_voltage_control` (device-level demo)

**Scenario:** Three-terminal voltage control (V1, V2, V3) — **CLI demo only**,
no counterpart in the §2.3.3 detailed-`P_sw` measurement (the actual experiment
drives `I_SOT` directly on the SOT channel).
**Physics focus:** Coupled electrical control of **V_MTJ** and **I_SOT**
through the electronic module.

### What happens

- Uses terminal voltages (V1, V2, V3) to compute (I_SOT, V_MTJ)
- Feeds (I_SOT, V_MTJ) into the switching dynamics
- Produces `mz(t)`, `R_MTJ(t)`, `V_MTJ(t)`, plus the V1/V2/V3 history

### Physical meaning

A **device-level control demo** showing how external terminal voltages map into:

- SOT current (torque-driven switching)
- MTJ bias (anisotropy modulation via VCMA, if enabled)

### Use this when

You want to test how a particular external pulse pattern (e.g. an actual write
circuit drive waveform) translates into the internal `(I_SOT, V_MTJ)` pair and
ultimately into a magnetisation trajectory.

### Reference figure

[`Chapter02_local_01.png`](../article/00_chapter_drafts/figs/Chapter02_local_01.png)
— three-terminal T-circuit schematic (chapter §2.1.1), produced by
[`demo/plot_t_circuit.py`](../demo/plot_t_circuit.py). This shows the static
mapping inside `electronic.electronic(V1, V2, V3, R_MTJ) → (I_SOT, V_MTJ)`.

------

## `sot_only_constant_current`

**Scenario:** Constant SOT current pulse only (`V_MTJ = 0`)
**Physics focus:** Pure **spin–orbit torque switching** baseline.

### Physical meaning

A **reference case** for:

- baseline switching trajectory without VCMA assistance
- comparisons against VCMA-assisted schemes

Pair with `enable_self_heating=True` to reproduce chapter §2.2.2.6 figures
showing how Joule heating shifts the SOT threshold as the pulse heats the
free layer.

### Reference figure

[`Chapter02_local_06.png`](../article/00_chapter_drafts/figs/Chapter02_local_06.png) —
`m_z(t)` + `T(t)` with vs without the coupled self-heating step, produced by
[`scripts/04_thermal_nonidealities/plot_self_heating_compare.py`](../scripts/04_thermal_nonidealities/plot_self_heating_compare.py).

------

## `sot_switching_no_vcma`

**Scenario:** SOT switching without VCMA (`VNV=0`), sweep `I_SOT`
**Physics focus:** switching driven by SOT only (no anisotropy modulation).

### Physical meaning

A **control experiment** answering:

> How much SOT current is needed for switching if VCMA is disabled?

### Default sweep window

With the calibrated `theta_SH=0.04` the deterministic threshold sits around
`|I_SOT| ≈ 1.1 mA`. The default `i_sot_list` therefore spans `−800 µA` to
`−200 µA` — both sides of threshold. Older configs targeting a literature
`theta_SH = 0.25` (threshold ~140 µA) will silently no-switch under the new
defaults; override `theta_SH` if you need that regime.

------

## `ser_sot_no_vcma_thermal`

**Scenario:** Monte-Carlo SER vs `I_SOT` without VCMA, with thermal noise
(`NON=1`)
**Physics focus:** thermal stochasticity and reliability.

### Physical meaning

Quantifies **switching error rate (SER)** under temperature-driven noise:

- switching becomes probabilistic
- useful for device-level reliability evaluation

### Result object

Returns `SerResult(x, ser, x_label)`. Since 2026-05 the result also exposes the
convenience property `psw = 1 − ser`, which is what chapter §2.3.3 reports
directly — figure scripts under `09_simulation_figures/` default to `psw`.

### Toggle highlights

- `enable_self_heating=True`: each MC trial advances `T(t)` independently and
  recomputes `M_s(T)`, `K_i(T)` per step; the resulting SER curve is
  temperature-corrected. Verified against the closed-form steady-state in
  `tests/test_toggles.py`.
- `seed=<int>` + `rng_mode={"legacy","generator"}`: deterministically derive
  the per-trial seed from `(seed, I_SOT, trial_idx)`; chapter figures use the
  `"legacy"` default; multiprocessing pipelines should pick `"generator"`.
- `integrator="cayley"`: forces every trial through the norm-preserving
  Cartesian stepper; useful for precision studies but reports a slightly higher
  threshold than the spherical-Euler default the calibration was performed on.

⚠️ Computationally heavy (Monte-Carlo trials). Defaults to `trials=200`; chapter
figures use `trials=80` per point with the sweep widened to capture the full
sigmoid.

### Reference figures

- [`Chapter02_local_08.png`](../article/00_chapter_drafts/figs/Chapter02_local_08.png)
  — four representative single trajectories across the threshold zone
  (`scripts/09_simulation_figures/plot_single_trajectory.py`).
- [`Chapter02_local_10.png`](../article/00_chapter_drafts/figs/Chapter02_local_10.png)
  — Monte-Carlo `P_sw(|I_SOT|)` at 0.75 ns, self-heating OFF vs ON, with
  Wilson 95 % CI + threshold-region inset
  (`scripts/09_simulation_figures/plot_ser_mc.py`).
- [`Chapter02_local_13.png`](../article/00_chapter_drafts/figs/Chapter02_local_13.png)
  — the experimental same-batch reference Sigmoid this case is calibrated against
  (`scripts/07_process_variability/sigmoid_fig.py`).

------

> **Removed in 2026-05.** The earlier `vcma_assisted_switching_{isot,vmtj}_sweep`,
> `optimized_vgsot_switching`, and `ser_optimized_vgsot` cases (plus the
> `run_two_pulse_optimized` low-level kernel that backed them) were dropped
> because they came from the upstream Verilog-A port and had no counterpart in
> the same-batch Device A measurements that the package is now calibrated
> against. The CLI surface is now the four SOT-only / device-demo cases above
> plus the Python-only `variability_sweep` below.

------

## `variability_sweep` (analysis case, no CLI)

**Scenario:** D2D Monte-Carlo over `CV(Δ)`, fit Sigmoid per wafer, return
broadening factor `F(CV)` and joint prediction `β^eff(CV) = η_c · F(CV) · β_NB^fit`.
**Physics focus:** chapter §2.3.5 process-variability budget.

### Physical meaning

A self-contained version of the analysis previously distributed across the
external `07_process_variability/variability_sim.py` driver. Wraps
`vgsot_sim.analysis.variability.transfer_function_F` so the whole §2.3.5
Brinkman/MC budget can be regenerated by one Python call:

```python
from vgsot_sim.ser_cases import variability_sweep

res = variability_sweep(
    Delta=5.15,       # NB fit Δ for the device
    Vc0=0.884,        # NB fit V_c0
    tw_ns=0.75,
    beta_meas=44.6,   # measured Sigmoid slope (V⁻¹) — Device A P→AP basis
)

# res.cv_sweep, res.beta_eff, res.vth_eff, res.F_func, res.eta_c, res.beta_combined
```

Returns a `VariabilitySweepResult` (see [`ser_cases.py`](../src/vgsot_sim/ser_cases.py)).
The lower-level math (`nb_fit`, `sigmoid_fit`, `variability`, `sampling`) is also
importable from `vgsot_sim.analysis`.

This case is not exposed through the CLI — call it directly from Python.

### Reference figures

- [`Chapter02_local_15.png`](../article/00_chapter_drafts/figs/Chapter02_local_15.png)
  — full §2.3.5 four-panel summary (Brinkman CV(Δ) budget, wafer-averaged
  Sigmoid family, transfer function F(CV_Δ), joint prediction vs experiment),
  produced by [`scripts/07_process_variability/variability_sim.py`](../scripts/07_process_variability/variability_sim.py).
- [`Chapter02_local_16.png`](../article/00_chapter_drafts/figs/Chapter02_local_16.png)
  — macrospin D2D Monte-Carlo sampling
  (`scripts/07_process_variability/macrospin_process_variability_mc.py`).
- [`Chapter02_local_17.png`](../article/00_chapter_drafts/figs/Chapter02_local_17.png)
  and [`Chapter02_local_18.png`](../article/00_chapter_drafts/figs/Chapter02_local_18.png)
  — Wilson coverage + MC sampling-size sensitivity (chapter §2.3.6,
  `scripts/08_sampling_effect/`).

------

## Chapter-figure pipelines (`scripts/`)

The cases above are the *physical experiments*. The thesis figures themselves
are produced by small driver scripts under `scripts/`, organised by chapter
section. The §2.3.3–§2.3.4 vgsot-sim figures live under
[`scripts/09_simulation_figures/`](../scripts/09_simulation_figures/) — see
that folder's [README](../scripts/09_simulation_figures/README.md) for the
full driver list (`plot_single_trajectory.py`, `plot_3d_trajectory.py`,
`plot_ser_mc.py`, `calibrate_to_experiment.py`) and the CLI toggles they expose.
