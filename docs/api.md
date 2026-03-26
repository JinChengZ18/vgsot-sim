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
    VcmaAssistedSwitchingIsotSweepConfig,
    VcmaAssistedSwitchingVmtjSweepConfig,
    OptimizedVgsotSwitchingConfig,
    SerOptimizedVgsotConfig,

    # result types
    SimResult,
    SweepResult,
    SerResult,
    SerOptimizedResult,

    # high-level cases
    terminal_voltage_control,
    sot_only_constant_current,
    sot_switching_no_vcma,
    vcma_assisted_switching_isot_sweep,
    vcma_assisted_switching_vmtj_sweep,
    optimized_vgsot_switching,
    ser_sot_no_vcma_thermal,
    ser_optimized_vgsot,

    # low-level kernels
    run_piecewise_terminal_voltage,
    run_piecewise_direct_excitation,
    run_two_pulse_optimized,
)
```



## Mental model

The library has three layers.

1. **Config dataclasses** describe one experiment.
2. **Case functions** run one complete experiment and return structured results.
3. **Kernel functions** expose lower-level simulation entry points for custom waveforms.

A typical workflow is:

```python
from vgsot_sim import (
    VcmaAssistedSwitchingIsotSweepConfig,
    vcma_assisted_switching_isot_sweep,
)

cfg = VcmaAssistedSwitchingIsotSweepConfig(
    v_mtj=1.1,
    i_sot_list=[-40e-6, -30e-6, -20e-6],
)
res = vcma_assisted_switching_isot_sweep(cfg)
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

Returned by single-run simulations such as `terminal_voltage_control`, `sot_only_constant_current`, `run_piecewise_direct_excitation`, and `run_two_pulse_optimized`.

| Field | Type | Meaning |
|---|---|---|
| `time_s` | `np.ndarray` | Simulation time axis in seconds. |
| `mz` | `np.ndarray` | Magnetization z-component over time. |
| `r_mtj` | `np.ndarray` | MTJ resistance over time. |
| `v_mtj` | `np.ndarray` | MTJ voltage waveform actually used in simulation. |
| `i_sot` | `np.ndarray` | SOT current waveform actually used in simulation. |
| `switch_energy_j` | `float` | Estimated MTJ switching energy in joules. |
| `theta` | `np.ndarray | None` | Polar angle history. |
| `phi` | `np.ndarray | None` | Azimuth angle history. |
| `v1` | `np.ndarray | None` | Terminal voltage V1 history, only for terminal-voltage mode. |
| `v2` | `np.ndarray | None` | Terminal voltage V2 history, only for terminal-voltage mode. |
| `v3` | `np.ndarray | None` | Terminal voltage V3 history, only for terminal-voltage mode. |

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
| `x_label` | `str` | Suggested x-axis label. |

### `SerOptimizedResult`

Returned by `ser_optimized_vgsot`.

| Field | Type | Meaning |
|---|---|---|
| `t1_s` | `np.ndarray` | First-pulse duration values in seconds. |
| `ser` | `np.ndarray` | Switching error rate for each `t1`. |
| `mz_at_t1_avg` | `np.ndarray` | Average `mz` sampled at the end of the first pulse. |



## Shared physical flags (PAP / NON / VNV)

These parameters appear in many APIs and are defined globally:

| Param | Meaning                                                 |
| ----- | ------------------------------------------------------- |
| `pap` | Initial magnetic state: 1 = Parallel, 0 = Anti-parallel |
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
| `i_sot_stage1` | `-95e-6` | SOT current during stage 1. |
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
| `i_sot_list` | `(-85e-6, -90e-6, -95e-6, -100e-6)` | SOT currents to sweep. |
| `sim_start_step` | `1` | Simulation start index. |
| `sim_mid1_step` | `2000` | End of pulse stage. |
| `sim_end_step` | `5000` | Final simulation index. |
| `pap` | `1` | Initial state selector. |
| `non` | `1` | Thermal-noise toggle. |
| `v_mtj` | `0.0` | Constant MTJ voltage for all stages. |
| `i_sot_relax` | `0.0` | Relaxation-stage current after the main pulse. |
| `vnv` | `0` | VCMA toggle. |
| `r_sot_fl_dl` | `0.0` | Damping-like ratio. |
| `tick_spacing_s` | `5e-10` | Suggested plot tick spacing. |

**Returns:** `SweepResult`

`pulse_curves` stores `I_SOT` in microamps for this case, and `pulse_ylabel` is already set accordingly.

### `vcma_assisted_switching_isot_sweep`

Fixes `V_MTJ` and sweeps `I_SOT`.

```python
from vgsot_sim import (
    VcmaAssistedSwitchingIsotSweepConfig,
    vcma_assisted_switching_isot_sweep,
)

cfg = VcmaAssistedSwitchingIsotSweepConfig(
    v_mtj=1.1,
    i_sot_list=[-40e-6, -30e-6, -20e-6],
)
res = vcma_assisted_switching_isot_sweep(cfg)
```

**Config:** `VcmaAssistedSwitchingIsotSweepConfig`

| Field | Default | Meaning |
|---|---:|---|
| `v_mtj` | `1.2` | Constant MTJ bias during the full simulation. |
| `i_sot_list` | `(-90e-6, -30e-6, -18e-6, -16e-6)` | Currents to sweep. |
| `sim_start_step` | `1` | Simulation start index. |
| `sim_end_step` | `25000` | Final simulation index. |
| `pap` | `1` | Initial state selector. |
| `non` | `1` | Thermal-noise toggle. |
| `vnv` | `1` | VCMA toggle. |
| `r_sot_fl_dl` | `0.0` | Damping-like ratio. |
| `tick_spacing_s` | `5e-9` | Suggested plot tick spacing. |

**Returns:** `SweepResult`

### `vcma_assisted_switching_vmtj_sweep`

Fixes `I_SOT` and sweeps `V_MTJ`.

```python
from vgsot_sim import (
    VcmaAssistedSwitchingVmtjSweepConfig,
    vcma_assisted_switching_vmtj_sweep,
)

cfg = VcmaAssistedSwitchingVmtjSweepConfig(
    i_sot=-50e-6,
    v_mtj_list=[1.2, 1.3, 1.4],
)
res = vcma_assisted_switching_vmtj_sweep(cfg)
```

**Config:** `VcmaAssistedSwitchingVmtjSweepConfig`

| Field | Default | Meaning |
|---|---:|---|
| `i_sot` | `None` | SOT current. If `None`, the code computes an internal default from device constants. |
| `v_mtj_list` | `(1.3189, 1.3191, 1.333, 1.3489, 1.4937)` | MTJ voltages to sweep. |
| `sim_start_step` | `1` | Simulation start index. |
| `sim_end_step` | `25000` | Final simulation index. |
| `pap` | `1` | Initial state selector. |
| `non` | `1` | Thermal-noise toggle. |
| `vnv` | `1` | VCMA toggle. |
| `r_sot_fl_dl` | `0.0` | Damping-like ratio. |
| `tick_spacing_s` | `5e-9` | Suggested plot tick spacing. |

**Returns:** `SweepResult`

### `optimized_vgsot_switching`

Runs the proposed two-pulse VGSOT scheme for several `(t1, t2)` pairs.

```python
from vgsot_sim import OptimizedVgsotSwitchingConfig, optimized_vgsot_switching

cfg = OptimizedVgsotSwitchingConfig(
    t_pairs_s=[(1.4e-9, 1.6e-9), (1.8e-9, 1.2e-9)],
)
res = optimized_vgsot_switching(cfg)
```

**Config:** `OptimizedVgsotSwitchingConfig`

| Field | Default | Meaning |
|---|---:|---|
| `v_mtj_1` | `1.4937` | MTJ voltage during first pulse. |
| `v_mtj_2` | `-1.0` | MTJ voltage during second pulse. |
| `i_sot` | `None` | First-stage SOT current. If `None`, the code computes an internal default. |
| `t_pairs_s` | see code | Sequence of `(t1_s, t2_s)` pulse durations. |
| `sim_total_time_s` | `25e-9` | Total simulated time. |
| `pap` | `1` | Initial state selector. |
| `non` | `1` | Thermal-noise toggle. |
| `vnv` | `1` | VCMA toggle. |
| `r_sot_fl_dl` | `0.0` | Damping-like ratio. |
| `tick_spacing_s` | `5e-9` | Suggested plot tick spacing. |

**Returns:** `SweepResult`

`pulse_curves` stores the MTJ voltage waveform for this case.

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
| `i_sot_list` | `(-100e-6, -98e-6, -96e-6, -94e-6, -92e-6, -90e-6)` | Currents to sweep. |
| `trials` | `1000` | Monte-Carlo trials per point. |
| `sim_start_step` | `1` | Simulation start index. |
| `sim_mid1_step` | `2000` | End of active pulse. |
| `sim_end_step` | `5000` | Final simulation index. |
| `pap` | `1` | Initial state selector. |
| `non` | `1` | Thermal-noise toggle. |
| `v_mtj` | `0.0` | MTJ bias. |
| `vnv` | `0` | VCMA toggle. |
| `r_sot_fl_dl` | `0.0` | Damping-like ratio. |
| `target_mz` | `-1.0` | Desired final magnetization for success. |
| `failure_tol` | `1e-1` | Allowed deviation from `target_mz`. |

**Returns:** `SerResult`

### `ser_optimized_vgsot`

Monte-Carlo SER for the proposed two-pulse scheme.

```python
from vgsot_sim import SerOptimizedVgsotConfig, ser_optimized_vgsot

cfg = SerOptimizedVgsotConfig(
    iterations_num=100,
    t1_list_s=[1.4e-9, 1.5e-9, 1.6e-9],
    total_pulse_s=3e-9,
)
res = ser_optimized_vgsot(cfg)
```

**Config:** `SerOptimizedVgsotConfig`

| Field | Default | Meaning |
|---|---:|---|
| `iterations_num` | `300` | Monte-Carlo trials per `t1`. |
| `v_mtj_1` | `1.4937` | First-pulse MTJ voltage. |
| `v_mtj_2` | `-1.0` | Second-pulse MTJ voltage. |
| `i_sot` | `None` | First-stage current. If `None`, an internal default is computed. |
| `t1_list_s` | see code | Candidate first-pulse durations. |
| `total_pulse_s` | `3e-9` | Fixed `t1 + t2`. |
| `sim_total_time_s` | `25e-9` | Total simulated time. |
| `pap` | `1` | Initial state selector. |
| `non` | `1` | Thermal-noise toggle. |
| `vnv` | `1` | VCMA toggle. |
| `r_sot_fl_dl` | `0.0` | Damping-like ratio. |
| `target_final_mz` | `1.0` | Desired final magnetization. |
| `failure_tol` | `1e-1` | Allowed deviation from target. |

**Returns:** `SerOptimizedResult`

## Low-level kernels for custom waveforms

Use these when the built-in cases are not flexible enough.

### `run_piecewise_terminal_voltage`

Accepts a `TerminalVoltageControlConfig` and simulates up to three voltage stages.

```python
from vgsot_sim import TerminalVoltageControlConfig, run_piecewise_terminal_voltage

cfg = TerminalVoltageControlConfig(
    sim_mid1_step=2000,
    sim_mid2_step=3500,
    sim_end_step=5000,
    v_stage1=(1.0, 0.0, 0.1),
    v_stage2=(-1.0, 0.0, 0.0),
    v_stage3=(0.0, 0.0, 0.0),
)
res = run_piecewise_terminal_voltage(cfg)
```

This is the lowest-level entry point for cases driven by terminal voltages rather than directly specified `V_MTJ` and `I_SOT`.

### `run_piecewise_direct_excitation`

Most flexible low-level API for directly specifying up to three piecewise-constant stages of `V_MTJ` and `I_SOT`.

```python
from vgsot_sim import run_piecewise_direct_excitation

res = run_piecewise_direct_excitation(
    sim_start_step=1,
    sim_mid1_step=2000,
    sim_mid2_step=3500,
    sim_end_step=5000,
    pap=1,

    v_mtj_stage1=0.0,
    v_mtj_stage2=0.8,
    v_mtj_stage3=0.0,

    i_sot_stage1=-90e-6,
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

This function is the recommended choice for user-defined protocols.

### `run_two_pulse_optimized`

Convenience wrapper for the built-in two-pulse VGSOT protocol.

```python
from vgsot_sim import run_two_pulse_optimized

res = run_two_pulse_optimized(
    t1_s=1.5e-9,
    t2_s=1.5e-9,
    v_mtj_1=1.4937,
    v_mtj_2=-1.0,
    i_sot_1=-50e-6,
    i_sot_2=0.0,
    sim_total_time_s=25e-9,
    pap=1,
    non=0,
    vnv=1,
    r_sot_fl_dl=0.0,
)
```

Internally it converts `t1_s` and `t2_s` into step boundaries and then calls `run_piecewise_direct_excitation`.

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

### Custom protocol 3: two-pulse timing exploration

```python
from vgsot_sim import run_two_pulse_optimized

for t1_s in [1.4e-9, 1.6e-9, 1.8e-9]:
    t2_s = 3e-9 - t1_s
    res = run_two_pulse_optimized(
        t1_s=t1_s,
        t2_s=t2_s,
        v_mtj_1=1.4937,
        v_mtj_2=-1.0,
        i_sot_1=-50e-6,
        i_sot_2=0.0,
        sim_total_time_s=25e-9,
        pap=1,
        non=0,
        vnv=1,
        r_sot_fl_dl=0.0,
        show_progress=False,
    )
    print(t1_s, res.switch_energy_j, res.mz[-1])
```

## Compatibility notes

- Older examples may use `result.curves`. In current versions, prefer `result.mz_curves`. The old attribute is still available as a compatibility alias on `SweepResult`.
- Older examples may omit `sim_mid2_step` in `run_piecewise_direct_excitation`. Current code expects a three-stage interface, so include it explicitly.
- Sweep CSV export should use `save_grouped_timeseries_csv(...)` rather than flattening several dicts with overlapping labels.
