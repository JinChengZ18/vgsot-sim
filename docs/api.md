# API Reference

This document is the detailed Python API guide for `vgsot_sim`.
It is meant to complement the short examples in the root `README.md`.



## Import surface

The package re-exports the main configuration dataclasses, result dataclasses, high-level cases, low-level building blocks, and result IO helpers from `vgsot_sim.__init__`.

```python
from vgsot_sim import (
    # configs
    TerminalVoltageControlConfig,
    SotOnlyConstantCurrentConfig,
    SotSwitchingNoVcmaConfig,
    SerSotNoVcmaThermalConfig,

    # result types
    SimResult,
    SweepResult,
    SerResult,

    # high-level cases
    terminal_voltage_control,
    sot_only_constant_current,
    sot_switching_no_vcma,
    ser_sot_no_vcma_thermal,

    # low-level kernels
    run_piecewise_terminal_voltage,
    run_piecewise_direct_excitation,
)

# Process-variability analysis (Python-only, no CLI surface)
from vgsot_sim.ser_cases import variability_sweep, VariabilitySweepResult
```



## Mental model

The library has three layers.

1. **Config dataclasses** describe one experiment.
2. **Case functions** run one complete experiment and return structured results.
3. **Kernel functions** expose lower-level simulation entry points for custom waveforms.

A typical workflow is:

```python
from vgsot_sim import (
    SerSotNoVcmaThermalConfig,
    ser_sot_no_vcma_thermal,
)

# Defaults already match the chapter §2.3.3 protocol
# (0.75 ns write + 3.25 ns relax, |I_SOT| spans 0.7–1.5 mA).
cfg = SerSotNoVcmaThermalConfig(trials=80)
res = ser_sot_no_vcma_thermal(cfg, seed=2026, enable_self_heating=True)
print(res.psw)
```



## Configuration system

```
# View default parameters
cfg = SerSotNoVcmaThermalConfig()
print(cfg)

# Or convert to dict
from vgsot_sim import config_to_params
print(config_to_params(cfg))

# Override some parameters
cfg = SerSotNoVcmaThermalConfig(trials=200)
```



## Result objects

### `SimResult`

Returned by single-run simulations such as `terminal_voltage_control`, `sot_only_constant_current`, and `run_piecewise_direct_excitation`.

| Field | Type | Meaning |
|---|---|---|
| `time_s` | `np.ndarray` | Simulation time axis in seconds. |
| `mz` | `np.ndarray` | Magnetization z-component over time. |
| `r_mtj` | `np.ndarray` | MTJ resistance over time. |
| `v_mtj` | `np.ndarray` | MTJ voltage waveform actually used in simulation. |
| `i_sot` | `np.ndarray` | SOT current waveform actually used in simulation. |
| `switch_energy_j` | `float` | Estimated MTJ switching energy in joules. |
| `theta` | `np.ndarray \| None` | Polar angle history. |
| `phi` | `np.ndarray \| None` | Azimuth angle history. |
| `v1` | `np.ndarray \| None` | Terminal voltage V1 history, only for terminal-voltage mode. |
| `v2` | `np.ndarray \| None` | Terminal voltage V2 history, only for terminal-voltage mode. |
| `v3` | `np.ndarray \| None` | Terminal voltage V3 history, only for terminal-voltage mode. |
| `T_K` | `np.ndarray \| None` | Per-step temperature, populated when `enable_self_heating=True`. |
| `Ms_T` | `np.ndarray \| None` | Per-step `M_s(T)` fed into LLG, populated when `enable_self_heating=True`. |
| `Ki_T` | `np.ndarray \| None` | Per-step `K_i(T)` fed into LLG, populated when `enable_self_heating=True`. |

### `SweepResult`

Returned by sweep-style time-domain simulations.

| Field | Type | Meaning |
|---|---|---|
| `time_s` | `np.ndarray` | Shared time axis. |
| `mz_curves` | `dict[str, np.ndarray]` | `mz(t)` for each sweep condition. |
| `r_mtj_curves` | `dict[str, np.ndarray]` | `R_MTJ(t)` for each sweep condition. |
| `pulse_curves` | `dict[str, np.ndarray]` | Applied pulse waveform for each condition. Units depend on the case. |
| `switch_energy_j` | `dict[str, float]` | Switching energy per sweep condition. |
| `pulse_ylabel` | `str` | Suggested y-label for plotting `pulse_curves`. |
| `curves` | property | Backward-compatible alias for `mz_curves`. |

### `SerResult`

Returned by `ser_sot_no_vcma_thermal`.

| Field | Type | Meaning |
|---|---|---|
| `x` | `np.ndarray` | Sweep axis values. |
| `ser` | `np.ndarray` | Switching error rate for each point. |
| `psw` | `property → np.ndarray` | Switching success probability `P_sw = 1 − ser`. The thesis §2.3.3 reports `P_sw` directly; chapter figures default to this convention. |
| `x_label` | `str` | Suggested x-axis label. |

### `VariabilitySweepResult`

Returned by `variability_sweep` (chapter §2.3.5 D2D Monte-Carlo). See the
[Analysis sub-package](#analysis-sub-package-new-in-2026-05) section below
for fields and usage.

## Shared physical flags (PAP / NON / VNV)

These parameters appear in many APIs and are defined globally:

| Param | Meaning                                                 |
| ----- | ------------------------------------------------------- |
| `pap` | Initial magnetic state: 1 = **anti-parallel** start (m_z ≈ −1, R ≈ R_AP), 0 = parallel start (m_z ≈ +1, R = R_P). The name reads backwards; `tmr()` defines m_z = +1 as parallel. |
| `non` | Thermal noise toggle (1 = enable, 0 = disable)          |
| `vnv` | VCMA effect toggle (1 = enable, 0 = disable)            |



## High-level case APIs

These are the recommended entry points for most users.

### `terminal_voltage_control`

Three-terminal voltage driven simulation. The electrical submodel converts `(V1, V2, V3, R_MTJ)` into `(I_SOT, V_MTJ)` at each step.

```python
from vgsot_sim import TerminalVoltageControlConfig, terminal_voltage_control

cfg = TerminalVoltageControlConfig()
res = terminal_voltage_control(cfg)
```

**Config:** `TerminalVoltageControlConfig`

| Field | Default | Meaning |
|---|---:|---|
| `sim_start_step` | `1` | Simulation start index. |
| `sim_mid1_step` | `2000` | Boundary between stage 1 and stage 2. |
| `sim_mid2_step` | `None` | Boundary between stage 2 and stage 3. If `None`, stage 3 starts at `sim_mid1_step`. |
| `sim_end_step` | `5000` | Final simulation index. |
| `pap` | `1` | Initial magnetic state selector passed to `initialize.init`. |
| `v_stage1` | `(1.0, 0.0, 0.1)` | `(V1, V2, V3)` in stage 1. |
| `v_stage2` | `(-1.0, 0.0, 0.0)` | `(V1, V2, V3)` in stage 2. |
| `v_stage3` | `(0.0, 0.0, 0.0)` | `(V1, V2, V3)` in stage 3. |
| `estt_stage1/2/3` | varies | Enable STT contribution in each stage. |
| `esot_stage1/2/3` | varies | Enable SOT contribution in each stage. |
| `vnv` | `1` | Toggle VCMA-related term. Set to 0 to disable VCMA effect. |
| `non` | `1` | Toggle thermal-noise. Set to 0 to disable thermal effect. |
| `r_sot_fl_dl` | `0.83` | Damping-like ratio used by switching kernel. |
| `tick_spacing_s` | `5e-10` | Suggested plot tick spacing. |

**Returns:** `SimResult`

### `sot_only_constant_current`

Baseline SOT-only case using directly specified current and MTJ voltage values, typically one active pulse followed by relaxation.

```python
from vgsot_sim import SotOnlyConstantCurrentConfig, sot_only_constant_current

cfg = SotOnlyConstantCurrentConfig(i_sot_stage1=-95e-6)
res = sot_only_constant_current(cfg)
```

**Config:** `SotOnlyConstantCurrentConfig`

| Field | Default | Meaning |
|---|---:|---|
| `sim_start_step` | `1` | Simulation start index. |
| `sim_mid1_step` | `2000` | End of pulse stage. |
| `sim_end_step` | `5000` | Final simulation index. |
| `pap` | `1` | Initial state selector. |
| `i_sot_stage1` | `-400e-6` | SOT current during stage 1 (~3× sub-threshold for the calibrated `theta_SH=0.066`). |
| `i_sot_stage2` | `0.0` | SOT current during stage 2 and stage 3. |
| `v_mtj_stage1` | `0.0` | MTJ voltage during stage 1. |
| `v_mtj_stage2` | `0.0` | MTJ voltage during stage 2 and stage 3. |
| `vnv` | `1` | VCMA toggle. |
| `non` | `1` | Thermal-noise toggle. |
| `r_sot_fl_dl` | `0.83` | Damping-like ratio. |
| `tick_spacing_s` | `5e-10` | Suggested plot tick spacing. |

**Returns:** `SimResult`

### `sot_switching_no_vcma`

Sweeps `I_SOT` with `V_MTJ=0`, useful as the no-VCMA baseline.

```python
from vgsot_sim import SotSwitchingNoVcmaConfig, sot_switching_no_vcma

cfg = SotSwitchingNoVcmaConfig(i_sot_list=[-85e-6, -90e-6, -95e-6])
res = sot_switching_no_vcma(cfg)
print(res.mz_curves.keys())
```

**Config:** `SotSwitchingNoVcmaConfig`

| Field | Default | Meaning |
|---|---:|---|
| `i_sot_list` | `(-800e-6, -700e-6, -750e-6, -600e-6, -500e-6, -400e-6, -300e-6, -200e-6)` | SOT currents to sweep. Bracket the calibrated `theta_SH = 0.066` threshold near 1.1 mA. |
| `sim_start_step` | `1` | Simulation start index. |
| `sim_mid1_step` | `5000` | End of pulse stage (5 ns at `t_step = 1 ps`). |
| `sim_end_step` | `10000` | Final simulation index. |
| `pap` | `1` | Initial state selector. |
| `non` | `1` | Thermal-noise toggle. |
| `v_mtj` | `0.0` | Constant MTJ voltage for all stages. |
| `i_sot_relax` | `0.0` | Relaxation-stage current after the main pulse. |
| `vnv` | `0` | VCMA toggle. |
| `r_sot_fl_dl` | `0.83` | Damping-like ratio. |
| `tick_spacing_s` | `1e-9` | Suggested plot tick spacing. |

**Returns:** `SweepResult`

`pulse_curves` stores `I_SOT` in microamps for this case, and `pulse_ylabel` is already set accordingly.

> **Removed in 2026-05** — `vcma_assisted_switching_isot_sweep`,
> `vcma_assisted_switching_vmtj_sweep`, and `optimized_vgsot_switching` were
> dropped together with their configs. None of them had a counterpart in the
> same-batch Device A measurements that the package is now calibrated against.

### `ser_sot_no_vcma_thermal`

Monte-Carlo SER vs `I_SOT` under thermal noise.

```python
from vgsot_sim import SerSotNoVcmaThermalConfig, ser_sot_no_vcma_thermal

cfg = SerSotNoVcmaThermalConfig(trials=200, i_sot_list=[-100e-6, -96e-6, -92e-6])
res = ser_sot_no_vcma_thermal(cfg)
```

**Config:** `SerSotNoVcmaThermalConfig`

| Field | Default | Meaning |
|---|---:|---|
| `i_sot_list` | `(-800e-6, -700e-6, -750e-6, -600e-6, -500e-6, -400e-6, -300e-6, -200e-6)` | Currents to sweep — bracket the calibrated `theta_SH = 0.066` threshold near `\|I_SOT\| ≈ 1.1 mA`. |
| `trials` | `200` | Monte-Carlo trials per point. |
| `sim_start_step` | `1` | Simulation start index. |
| `sim_mid1_step` | `5000` | End of active pulse (5 ns at `t_step = 1 ps`). |
| `sim_end_step` | `10000` | Final simulation index. |
| `pap` | `1` | Initial state selector (1 ≡ AP). |
| `non` | `1` | Thermal-noise toggle. |
| `v_mtj` | `0.0` | MTJ bias. |
| `vnv` | `0` | VCMA toggle. |
| `r_sot_fl_dl` | `0.0` | Damping-like ratio. |
| `target_mz` | `1.0` | Desired final magnetization for success (P, since `pap=1` starts AP). |
| `failure_tol` | `1e-1` | Allowed deviation from `target_mz`. |

**Run-time keyword arguments** (passed to `ser_sot_no_vcma_thermal(cfg, ...)`):

| Argument | Default | Meaning |
|---|---|---|
| `show_progress` | `True` | Toggle the tqdm bar (CLI inherits `--no_progress`). |
| `enable_self_heating` | `False` | Per-trial advance of `T(t)` with `M_s(T)` / `K_i(T)` feedback. |
| `T_ambient_K` | `300.0` | Sink temperature for the thermal network. |
| `seed` | `None` | Master MC seed; per-trial seeds derived from `_trial_seed(seed, i_sot, idx)`. |
| `rng_mode` | `"legacy"` | `"legacy"` = `np.random.seed` (chapter-figure bit-stream); `"generator"` = `np.random.default_rng` (multiprocessing-friendly). Both share the seed stream. |
| `integrator` | `None` | `None` → stepper default (`"cayley"`); `"euler_spherical"` or `"cayley"` force the step type. |

**Returns:** `SerResult` (with `psw = 1 − ser` property).

> **Removed in 2026-05** — `ser_optimized_vgsot` and its `SerOptimizedVgsotConfig`
> / `SerOptimizedResult` types were dropped together with the
> `optimized_vgsot_switching` case. The two-pulse VGSOT protocol came from the
> upstream Verilog-A port and has no counterpart in the same-batch Device A
> experiments. If you need a similar SER vs pulse-shape sweep, drive
> `run_piecewise_direct_excitation` from a Python loop directly.

## Low-level kernels for custom waveforms

Use these when the built-in cases are not flexible enough.

### `run_piecewise_direct_excitation`

Most flexible low-level API for directly specifying up to three piecewise-constant stages of `V_MTJ` and `I_SOT`.

```python
from vgsot_sim import run_piecewise_direct_excitation

res = run_piecewise_direct_excitation(
    sim_start_step=1,
    sim_mid1_step=750,         # 0.75 ns write pulse at t_step = 1 ps
    sim_mid2_step=4000,
    sim_end_step=4000,
    pap=1,

    v_mtj_stage1=0.0, v_mtj_stage2=0.0, v_mtj_stage3=0.0,
    i_sot_stage1=-2000e-6, i_sot_stage2=0.0, i_sot_stage3=0.0,  # super-threshold

    estt_stage1=0, esot_stage1=1,
    estt_stage2=0, esot_stage2=1,
    estt_stage3=0, esot_stage3=1,

    vnv=0,
    non=1,
    r_sot_fl_dl=0.83,

    # Optional toggles (default: chapter behaviour)
    enable_self_heating=True,
    T_ambient_K=300.0,
    integrator="euler_spherical",   # or "cayley"
    sigma_SH=None,                  # Cayley only; default = [-1, 0, 0]
    rng=None,                       # pass np.random.default_rng(seed) for byte-repro
    demag_mode="ellipsoid",         # or "thin_disk"
)
```

**Key points**

- `sim_mid2_step` is part of the current API. Older examples that omit it are outdated.
- If `sim_mid2_step=None`, the code internally replaces it with `sim_end_step`.
- The stage layout is:

```text
stage1: [sim_start_step, sim_mid1_step)
stage2: [sim_mid1_step, sim_mid2_step)
stage3: [sim_mid2_step, sim_end_step)
```

- `integrator="cayley"` dispatches to `dynamic_switching_vector.switching_vector`,
  which is norm-preserving to machine precision and takes the spin-Hall
  polarisation `sigma_SH` as an explicit 3-vector. `"cayley"` is the default
  and the kernel the calibration was performed on (θ_SH = 0.066 was
  re-extracted against it after the FL-SOT sign fix); `"euler_spherical"`
  remains available for the §2.2.3.2 integrator comparison.
- With `enable_self_heating=True` the returned `SimResult` additionally
  populates `T_K`, `Ms_T`, `Ki_T` so the temperature trajectory and the
  T-corrected material parameters are available for post-hoc inspection.
- `rng` accepts an `np.random.Generator`; the same generator state reproduces
  the trajectory bit-for-bit across `field()` → `switching()`/`switching_vector()`
  → `stochastic()`.

This function is the recommended choice for user-defined protocols.

### `run_piecewise_terminal_voltage`

Symmetric to `run_piecewise_direct_excitation`, but the upstream drive signal
is the terminal-voltage triple `(V1, V2, V3)` (with the electronic module
computing `I_SOT`, `V_MTJ` at each step). Since 2026-05 it accepts the same
opt-in trio `integrator=`, `sigma_SH=`, `rng=` so terminal-voltage simulations
can use the Cayley path, arbitrary `sigma_SH`, and byte-reproducible noise:

```python
import numpy as np
from vgsot_sim import TerminalVoltageControlConfig, run_piecewise_terminal_voltage

rng = np.random.default_rng(seed=2026)
cfg = TerminalVoltageControlConfig(v_stage1=(1.0, 0.0, 0.1))
res = run_piecewise_terminal_voltage(
    cfg,
    integrator="cayley",
    sigma_SH=np.array([0.0, -1.0, 0.0]),   # rotate σ̂_SH off the default -x
    rng=rng,
    enable_self_heating=True, T_ambient_K=300.0,
)
```

> **Removed in 2026-05** — `run_two_pulse_optimized` was dropped together with
> the `optimized_vgsot_switching` / `ser_optimized_vgsot` cases it backed. If
> you need a two-stage protocol, call `run_piecewise_direct_excitation` with
> `sim_mid1_step` and `sim_mid2_step` placed at the pulse boundaries.

## Low-level physics modules

These are useful when you want to build your own solver loop or inspect submodels.

### `initialize.init(PAP)`

Initializes the simulation state.

**Returns:** `(R_MTJ0, theta0, mz0, phi0)`

### `electronic.electronic(V1, V2, V3, R_MTJ)`

Computes electrical quantities from terminal voltages and current MTJ resistance.

**Returns:** `(I_SOT, V_MTJ)`

### `dynamic_switching.switching(...)`

Performs one magnetization update step.

```python
mz_next, phi_next, theta_next = switching(
    V_MTJ, I_SOT, R_MTJ, theta, phi,
    ESTT, ESOT,
    VNV=1, NON=0,
    R_SOT_FL_DL=0.83,
)
```

### `anisotropy.field(theta, phi, V_MTJ, n, NON, ENE, VNV)`

Computes the effective field used by the switching model.

### `stochastic.stochastic(n)`

Generates thermal fluctuation terms.

### `tmr.tmr(V_MTJ, mz)`

Updates MTJ resistance from the current `mz` and MTJ voltage.

## Output helpers

The root README only shows the simplest save example. The current project also includes grouped CSV export helpers for sweep results.

### Single-run CSV

```python
from vgsot_sim import ensure_result_dir, save_timeseries_csv

out_dir = ensure_result_dir("result")
save_timeseries_csv(
    out_dir / "single_run.csv",
    res.time_s,
    {
        "mz": res.mz,
        "r_mtj": res.r_mtj,
        "v_mtj": res.v_mtj,
        "i_sot": res.i_sot,
    },
)
```

### Sweep CSV with readable grouped headers

```python
from vgsot_sim import ensure_result_dir, save_grouped_timeseries_csv

out_dir = ensure_result_dir("result")
save_grouped_timeseries_csv(
    out_dir / "sweep.csv",
    res.time_s,
    {
        "mz": res.mz_curves,
        "r_mtj": res.r_mtj_curves,
        "pulse": res.pulse_curves,
    },
)
```

This produces columns like:

```text
time_s,mz__I_SOT=-40.0uA,r_mtj__I_SOT=-40.0uA,pulse__I_SOT=-40.0uA
```

### Plot helpers

- `save_single_plot(...)`
- `save_two_panel_plot(...)`
- `save_three_panel_plot(...)`
- `build_stem(case, cfg)`
- `ensure_result_dir(path)`

## Customization cookbook

### Custom protocol 1: three-stage direct waveform

Use `run_piecewise_direct_excitation` when you want full control over pulse timing and amplitudes.

```python
from vgsot_sim import run_piecewise_direct_excitation

res = run_piecewise_direct_excitation(
    sim_start_step=1,
    sim_mid1_step=1500,
    sim_mid2_step=3000,
    sim_end_step=6000,
    pap=1,
    v_mtj_stage1=1.2,
    v_mtj_stage2=-0.6,
    v_mtj_stage3=0.0,
    i_sot_stage1=-40e-6,
    i_sot_stage2=-10e-6,
    i_sot_stage3=0.0,
    estt_stage1=0,
    esot_stage1=1,
    estt_stage2=0,
    esot_stage2=1,
    estt_stage3=0,
    esot_stage3=1,
    vnv=1,
    non=0,
    r_sot_fl_dl=0.83,
)
```

### Custom protocol 2: user-defined sweep around a low-level kernel

```python
from vgsot_sim import run_piecewise_direct_excitation

curves = {}
for i_sot in [-20e-6, -30e-6, -40e-6]:
    res = run_piecewise_direct_excitation(
        sim_start_step=1,
        sim_mid1_step=2000,
        sim_mid2_step=3500,
        sim_end_step=5000,
        pap=1,
        v_mtj_stage1=1.1,
        v_mtj_stage2=0.0,
        v_mtj_stage3=0.0,
        i_sot_stage1=i_sot,
        i_sot_stage2=0.0,
        i_sot_stage3=0.0,
        estt_stage1=0,
        esot_stage1=1,
        estt_stage2=0,
        esot_stage2=1,
        estt_stage3=0,
        esot_stage3=1,
        vnv=1,
        non=0,
        r_sot_fl_dl=0.83,
        show_progress=False,
    )
    curves[f"I_SOT={i_sot*1e6:.1f}uA"] = res.mz
```

### Custom protocol 3: two-stage SOT pulse (replaces the removed two-pulse VGSOT helper)

```python
# Drive the §2.3.3 protocol with a custom two-stage SOT current — the same
# building block that the old `run_two_pulse_optimized` helper used internally.
from vgsot_sim import run_piecewise_direct_excitation
from vgsot_sim.configs import PhysicalConstantsConfig

cc = PhysicalConstantsConfig()
for t1_s in [0.5e-9, 0.75e-9, 1.0e-9]:
    sim_end = int(round(4e-9 / cc.t_step))
    mid1    = int(round(t1_s / cc.t_step))
    res = run_piecewise_direct_excitation(
        sim_start_step=1,
        sim_mid1_step=mid1,
        sim_mid2_step=sim_end,
        sim_end_step=sim_end,
        pap=1,
        v_mtj_stage1=0.0, v_mtj_stage2=0.0, v_mtj_stage3=0.0,
        i_sot_stage1=-1500e-6, i_sot_stage2=0.0, i_sot_stage3=0.0,
        estt_stage1=0, esot_stage1=1, estt_stage2=0, esot_stage2=1,
        estt_stage3=0, esot_stage3=1,
        vnv=0, non=1, r_sot_fl_dl=0.83,
        show_progress=False, constants=cc,
    )
    print(t1_s, res.switch_energy_j, res.mz[-1])
```

## Analysis sub-package (new in 2026-05)

Math helpers previously inlined in chapter scripts 06/07/08 are now importable
from `vgsot_sim.analysis`:

```python
from vgsot_sim.analysis import nb_fit, sigmoid_fit, variability, sampling
```

| Module | Highlights |
|---|---|
| `nb_fit` | `extract_Vc(V, R)`, `loglin_fit(t_w, V_c)`, `nb_params(a, b, tau0)`, `psw_nb(V, t_ns, Δ, V_c0)` — Néel-Brown thermal-activation fitting from hysteresis loops |
| `sigmoid_fit` | `sigmoid4p`, `fit_sigmoid(V, P) → SigmoidFitResult(y0, L, V_th, k, β, R²)`, `wilson(p, n)` |
| `variability` | `cv_delta_budget(...)`, `transfer_function_F(cv_sweep, Δ, V_c0, ...)`, `wafer_average_psw(...)`, `fit_sigmoid_to_average(...)` |
| `sampling` | `exact_coverage(K, p, ε)`, `K_required_exact(p, ε)`, `clt_K_required(p, ε)`, MC sampling-size sensitivity helpers |

The `variability_sweep` case in `vgsot_sim.ser_cases` wraps
`variability.transfer_function_F` and returns a `VariabilitySweepResult`:

```python
from vgsot_sim.ser_cases import variability_sweep

res = variability_sweep(Delta=5.15, Vc0=0.884, tw_ns=0.75, beta_meas=44.6)
res.cv_sweep        # array of CV(Δ) scanned
res.beta_eff        # wafer-averaged Sigmoid slope at each CV
res.F_func          # β_eff(CV) / β_eff(CV=0)
res.eta_c           # β_meas / β_NB^fit (C2C calibration factor)
res.beta_combined   # η_c · F(CV) · β_NB^fit — joint D2D + C2C prediction
```

## Temperature-coupled cases (`temperature_cases.py`)

Produces the chapter §2.2.2 figures directly without the chapter-04 driver
scripts having to re-implement the math:

- `sweep_material_temperature(...)` — `M_s(T)`, `K_i(T)`, `η(T)`, `K_U^eff(T)`
  closed-form scalings.
- `tmr_voltage_sweep(...)` — `R_P(V)`, `R_AP(V)` under either Lorentzian or
  PDK TMR(V) selectable at the config level.
- `thermal_transient(...)` — pure RC thermal response without LLG coupling
  (the LLG-coupled version is `run_piecewise_direct_excitation(enable_self_heating=True)`).

## Compatibility notes

- Older examples may use `result.curves`. In current versions, prefer `result.mz_curves`. The old attribute is still available as a compatibility alias on `SweepResult`.
- Older examples may omit `sim_mid2_step` in `run_piecewise_direct_excitation`. Current code expects a three-stage interface, so include it explicitly.
- Sweep CSV export should use `save_grouped_timeseries_csv(...)` rather than flattening several dicts with overlapping labels.
- The `PhysicalConstantsConfig` defaults changed in 2026-05 (`theta_SH 0.25 → 0.04`, `TMR 1.19 → 1.0`, `RA 36e-12 → 16.6e-12`, plus new fields `D_elec`, `t_mtj`, `Cv`, `lambda_MgO`, `t_MgO`, `T_RT`, `T_C`, `cc_exponent`, `H_k_eff_RT`, `R_series`, `h_ex_{x,y,z}`, `tmr_model`, `a_tmr/b_tmr/c_tmr/k_tmr`). Old code that instantiated the full dataclass with every field set explicitly will keep working; code that relied on the old defaults producing a ~140 µA threshold needs to either set `theta_SH=0.25` explicitly or move its sweep range to the millisecond scale.
- `SerResult.psw = 1 − SER` is the chapter convention. New figure scripts default to plotting `psw`; the legacy `ser` field is still there for back-compat.
