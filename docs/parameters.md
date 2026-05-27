

---

## Default parameters by case

> These are the **dataclass default values** in
> [`src/vgsot_sim/configs.py`](../src/vgsot_sim/configs.py). They determine
> what runs when you simply execute `vgsot-sim <case_name>` (the CLI just
> instantiates the matching `*Config()` with no overrides).
>
> **Sweep ranges are calibrated to `theta_SH = 0.04`** (chapter §2.3.4 default).
> The threshold sits around `|I_SOT| ≈ 1.1 mA`, so default sweeps span the
> millisecond range. Override `constants=PhysicalConstantsConfig(theta_SH=0.25)`
> if you need the literature β-W regime (~140 µA threshold).

### `TerminalVoltageControlConfig`

| Parameter               | Default                         | Meaning                                    |
| ----------------------- | ------------------------------- | ------------------------------------------ |
| `sim_start_step`        | `1`                             | Simulation start step index                |
| `sim_mid1_step`         | `2000`                          | Stage boundary (1 → 2) in steps            |
| `sim_mid2_step`         | `None` → equals `sim_end_step`  | Stage boundary (2 → 3) in steps            |
| `sim_end_step`          | `5000`                          | Simulation end step index                  |
| `pap`                   | `1`                             | Initial P/AP state selector (1 ≡ AP start) |
| `v_stage1`              | `(1.0, 0.0, 0.1)`               | (V1,V2,V3) in stage 1                      |
| `v_stage2`              | `(-1.0, 0.0, 0.0)`              | (V1,V2,V3) in stage 2                      |
| `v_stage3`              | `(0.0, 0.0, 0.0)`               | (V1,V2,V3) in stage 3                      |
| `estt_stage{1,2,3}`     | `0, 0, 1`                       | STT enable per stage                       |
| `esot_stage{1,2,3}`     | `1, 1, 1`                       | SOT enable per stage                       |
| `vnv`                   | `1`                             | VCMA enable flag passed to switching       |
| `non`                   | `1`                             | Thermal noise flag passed to switching     |
| `r_sot_fl_dl`           | `0.83`                          | Field-like / damping-like ratio            |
| `tick_spacing_s`        | `5e-10`                         | X-axis major tick spacing (s)              |
| `constants`             | `PhysicalConstantsConfig()`     | Physics defaults — see §below              |

------

### `SotOnlyConstantCurrentConfig`

Aligned to the chapter §2.3.3 detailed-`P_sw` protocol (0.75 ns write +
3.25 ns relaxation, super-threshold `I_SOT`, SOT-only).

| Parameter               | Default         | Meaning                                  |
| ----------------------- | --------------- | ---------------------------------------- |
| `sim_start_step`        | `1`             | Start step                               |
| `sim_mid1_step`         | `750`           | Pulse end step (0.75 ns @ `t_step=1 ps`) |
| `sim_end_step`          | `4000`          | End step (4 ns total → 3.25 ns relax)    |
| `pap`                   | `1`             | Initial state selector (1 ≡ AP start)    |
| `i_sot_stage1`          | `-1500e-6` A    | SOT current during pulse (~1.4× `I_th`)  |
| `i_sot_stage2`          | `0.0` A         | Relax stage SOT current                  |
| `v_mtj_stage1`          | `0.0` V         | MTJ bias during pulse (no VCMA)          |
| `v_mtj_stage2`          | `0.0` V         | MTJ bias during relax                    |
| `vnv`                   | `0`             | VCMA term disabled (matches experiment)  |
| `non`                   | `1`             | Thermal noise on                         |
| `r_sot_fl_dl`           | `0.83`          | FL/DL ratio                              |
| `tick_spacing_s`        | `5e-10`         | Tick spacing                             |

------

### `SotSwitchingNoVcmaConfig`

| Parameter               | Default                                                                | Meaning                                          |
| ----------------------- | ---------------------------------------------------------------------- | ------------------------------------------------ |
| `i_sot_list`            | `(-1500e-6, -1300e-6, -1200e-6, -1100e-6, -1000e-6, -900e-6, -700e-6)` A | Sweep brackets calibrated `I_th ≈ 1.1 mA`     |
| `sim_mid1_step`         | `750`                                                                  | Pulse end step (0.75 ns)                         |
| `sim_end_step`          | `4000`                                                                 | End step (4 ns)                                  |
| `v_mtj`                 | `0.0` V                                                                | MTJ bias (fixed)                                 |
| `i_sot_relax`           | `0.0` A                                                                | Relax stage current                              |
| `vnv`                   | `0`                                                                    | VCMA disabled                                    |
| `non`                   | `1`                                                                    | Thermal noise enabled                            |
| `r_sot_fl_dl`           | `0.83`                                                                 | FL/DL ratio                                      |
| `tick_spacing_s`        | `5e-10`                                                                | Tick spacing                                     |

------

### `SerSotNoVcmaThermalConfig`

| Parameter               | Default                                                                | Meaning                                              |
| ----------------------- | ---------------------------------------------------------------------- | ---------------------------------------------------- |
| `i_sot_list`            | `(-1500e-6, -1300e-6, -1200e-6, -1100e-6, -1000e-6, -900e-6, -700e-6)` A | Sweep brackets calibrated `I_th ≈ 1.1 mA`         |
| `trials`                | `100`                                                                  | Monte-Carlo trials per `I_SOT` (chapter fig uses 80) |
| `sim_mid1_step`         | `750`                                                                  | Pulse end step                                       |
| `sim_end_step`          | `4000`                                                                 | End step                                             |
| `v_mtj`                 | `0.0` V                                                                | MTJ bias (fixed)                                     |
| `vnv`                   | `0`                                                                    | VCMA disabled                                        |
| `non`                   | `1`                                                                    | Thermal noise enabled                                |
| `r_sot_fl_dl`           | `0.83`                                                                 | FL/DL ratio (unified with other SOT-only cases)      |
| `target_mz`             | `1.0`                                                                  | Successful final state (P, since `pap=1` starts AP)  |
| `failure_tol`           | `0.2`                                                                  | Failure threshold: `abs(mz_final − target_mz) > tol` |

Run-time toggles (passed to `ser_sot_no_vcma_thermal(cfg, ...)`):

| Argument | Default | Meaning |
|---|---|---|
| `enable_self_heating` | `False` | Couple T(t) → M_s/K_i feedback per step |
| `T_ambient_K` | `300.0` | Sink temperature for the thermal RC network |
| `seed` | `None` | Master MC seed (per-trial seeds derived from `(seed, I_SOT, idx)`) |
| `rng_mode` | `"legacy"` | `"legacy"` = `np.random.seed` (chapter figures); `"generator"` = `np.random.default_rng` |
| `integrator` | `None` | `None` → stepper default; `"euler_spherical"` or `"cayley"` to force |

------

> **Removed in 2026-05** — `VcmaAssistedSwitchingIsotSweepConfig`,
> `VcmaAssistedSwitchingVmtjSweepConfig`, `OptimizedVgsotSwitchingConfig`, and
> `SerOptimizedVgsotConfig` were dropped together with their case functions
> (and the `run_two_pulse_optimized` low-level kernel). They came from the
> upstream Verilog-A port and didn't map to the same-batch Device A
> measurements that the package is now calibrated against.

------

## `PhysicalConstantsConfig` — calibrated defaults

These are the device-physics parameters shared by every case via the
`constants` field on each config. They are calibrated to chapter §2.3.3 Device
A measurements (2026-05-16). The notable changes from older versions:

| Parameter | Default | Note |
|---|---|---|
| `theta_SH` | **`0.04`** | Calibrated to Device A P→AP V_th(0.75 ns) = 894 mV (was `0.25` from textbook β-W). For bare-material studies override explicitly. |
| `TMR` | **`1.0`** | ≈100 % to match the §2.3.3 hysteresis amplitude R_AP/R_P ≈ 2 (was `1.19`). |
| `RA` | **`16.6e-12`** Ω·m² | Calibrated to R_P ≈ 5 kΩ at D_elec = 65 nm (was `36e-12`). |
| `D_elec` | `65e-9` m | Electrical effective MTJ diameter (new; D_phys remains `80e-9`). |
| `t_mtj` | `20e-9` m | Total MTJ pillar thickness, used by the RC thermal time constant (new). |
| `T_RT`, `T_C`, `cc_exponent` | `300`, `1100`, `2.18` | Bloch / modified Callen-Callen scaling for `M_s(T)` / `K_i(T)` (new). |
| `H_k_eff_RT` | `5.0e4` A/m | Effective anisotropy field at RT, used by the T-scaling calibration formula (new). |
| `Cv`, `lambda_MgO`, `t_MgO` | `2.5e6`, `4.0`, `1.4e-9` | MgO thermal-network parameters (new). |
| `R_series` | `0.0` Ω | Optional parasitic contact / lead resistance (new opt-in). Set to non-zero to reproduce wafer-level `R(V)` measurements. |
| `h_ex_x, h_ex_y, h_ex_z` | `0, −50 Oe / (4π) along ŷ, 0` A/m | In-plane bias field; MUST be perpendicular to σ̂_SH = −x̂ for deterministic SOT switching (new). |
| `tmr_model` | `"pdk"` | Bias-dependent TMR(V) model selector (was implicit Lorentzian). Set to `"lorentzian"` for the legacy single-parameter form. |
| `k_tmr, a_tmr, b_tmr, c_tmr` | `1.2346, 0.1729, 0.1315, 0.4475` | Three-parameter PDK TMR(V) coefficients (new). |

See [`technical_details.md`](technical_details.md) §2.3 for the calibration story behind `theta_SH = 0.04` and the
diagnostic `09_simulation_figures/calibrate_to_experiment.py` script.

### Derived quantities

`PhysicalConstantsConfig` exposes these as `@property` (computed lazily so the
constants object stays single-source-of-truth):

| Property | Formula | Units |
|---|---|---|
| `gamma` | `2 μ_0 μ_B / ħ` | (A·s)/(kg·m²)·m/A — gyromagnetic ratio |
| `A2` | `d · w` (SOT channel cross-section) | m² |
| `R_W` | `ρ · l / (w · d)` (SOT channel resistance) | Ω |
| `A1` | `π D_elec² / 4` (electrical/magnetic effective area) | m² |
| `A1_phys` | `π D² / 4` (geometric area, thermal pillar) | m² |
| `v` | `t_f · A1` (magnetic active volume) | m³ |
| `Heff` | `2 K_i / (t_f M_s μ_0)` | A/m |
