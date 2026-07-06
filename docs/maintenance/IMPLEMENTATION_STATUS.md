# Implementation Status

Map of the **physical effects modeled in Chapter 2 of the thesis** versus what the current `vgsot-sim` code base actually implements. Items marked **Not implemented** are intentional gaps — the underlying physics is described in the thesis text but the corresponding numerical machinery is still pending.

Last updated: 2026-07-02 (RC deepening round: RTN node rigor fixes, `rtn/bridge.py` sLLG validation, `rtn/reservoir.py` reservoir layer + benchmarks, thesis §2.4 figures 2.19–2.21 / tables 2.14–2.15 via the `scripts/build_ppt_figs.py` panel flow — see `docs/maintenance/version_notes.md` 2026-07 and `scripts/10_rtn_reservoir/README.md`. Earlier 2026-05 round: opt-in toggles `R_series`, end-to-end `rng=…` plumbing, Psw↔SER display toggle, `rng_mode`/`integrator` selection; regression in `tests/test_toggles.py`.)

---

## Implemented

| Thesis section | Effect | Code location | Notes |
|---|---|---|---|
| 2.1.1 | Three-terminal T-circuit electronic mapping (V1, V2, V3) → (I_SOT, V_MTJ) | [`electronic.py`](../../src/vgsot_sim/electronic.py) | |
| 2.2.1.1 | LLG precession + Gilbert damping in spherical (θ, φ) coords | [`dynamic_switching.py`](../../src/vgsot_sim/dynamic_switching.py) | One-step explicit Euler |
| 2.2.1.1 | STT damping-like + field-like torque | `dynamic_switching.py:21-22, 30, 41` | Uses spin polarization $P$ |
| 2.2.1.1 | SOT damping-like + field-like torque | `dynamic_switching.py:23-24, 31-32, 42-43` | $\hat\sigma$ direction implicit — see "Known caveats" below |
| 2.2.1.2 | PMA effective field $H_{PMA} = 2K_i/(\mu_0 M_s t_f) m_z$ | [`anisotropy.py:35`](../../src/vgsot_sim/anisotropy.py) | Accepts optional `Ki_T`, `Ms_T` overrides for use inside a T-dependent loop |
| 2.2.1.2 | VCMA effective field $H_{VCMA} = -2\beta V_{MTJ} m_z / (\mu_0 M_s t_{ox} t_f)$ | `anisotropy.py:36` | β is the surface VCMA coefficient |
| 2.2.1.2 | Demagnetisation tensor — thin-disk **and** exact oblate-ellipsoid | [`demag.py`](../../src/vgsot_sim/demag.py), `anisotropy.field(... demag_mode=...)` | Default `"ellipsoid"`; uses `D_elec` for the in-plane diameter |
| 2.2.1.2 | Configurable in-plane bias field $H_{ex}$ | `anisotropy.py:30`, `configs.py:69-72` | `h_ex_x/y/z` — defaults to 200 Oe along +x |
| 2.2.1.2 | Thermal noise sampling — independent N(0,1) per axis, FDT amplitude | [`stochastic.py`](../../src/vgsot_sim/stochastic.py) + `anisotropy.py:35-37` | **Corrected from upstream** (see thesis §2.2.1.2 implementation note). `stochastic()` accepts optional `rng` for reproducibility. |
| 2.2.1.2 (P1) | Initial-state thermal sampling — Rayleigh `θ`, uniform `φ` | [`initialize.py:init`](../../src/vgsot_sim/initialize.py) | Each Monte-Carlo trial now starts from an independent thermal draw rather than a single deterministic offset. `init()` takes optional `rng`. |
| 2.2.3 (P1) | Polar-region numerical guard | `dynamic_switching.py:35-50` | `dphi/dt` frozen when `sin(θ) < 1e-8` to prevent `1/sin(θ)` roundoff blow-up |
| 2.2.2.3 (P1) | `compute_Rp()` separated from `init()` | `initialize.py:compute_Rp`, used by `tmr.py` | `tmr()` is no longer stochastic in $R_P$ (was previously coupled to the deterministic init angle) |
| — (P1) | CODATA-precision physical constants | `configs.py:24-30` | μ_0, e, ℏ, μ_B, k_B, m_e at 6–8 sig figs |
| 2.1.1, 2.2.2.3 | $R_P = RA/A_1$ (measured RA); genuine Simmons/BDR RA predictor where $\phi_{ox}$, $t_{ox}$ actually enter | [`initialize.py`](../../src/vgsot_sim/initialize.py) (`compute_Rp`, `resistance_area_bdr`) | RA is the calibration handle; the BDR predictor matches the measured RA at MgO $m^*\approx0.3\,m_e$ |
| 2.2.2.3 | TMR(V) – Lorentzian form | [`tmr.py`](../../src/vgsot_sim/tmr.py), `tmr_eff(model="lorentzian")` | |
| 2.2.2.3 | TMR(V) – PDK three-parameter form (default) | `tmr.py`, `tmr_eff(model="pdk")` | a/b/c/k_TMR exposed in config |
| 2.2.2.3 | Conductance-interpolation $R(m_z, V)$ | `tmr.py:tmr` | $m_z=+1 \to R_P$ |
| 2.1.3 | Joint $P_{sw}(t_w, I_{SOT}, V_{MTJ})$ via per-trial sLLG | `time_series_cases.py` + `ser_cases.py` | Monte-Carlo SER |
| 2.3.3 | Per-pulse switching probability via repeated trials | `ser_cases.py:ser_sot_no_vcma_thermal` etc. | |
| 2.2.2.4 | $D_{elec}$ vs $D_{phys}$ separation | `configs.py` (`D`, `D_elec`); `A1` uses $D_{elec}$, `A1_phys` for thermal | $R_P\approx5.0$ kΩ at the default $RA=16.6\times10^{-12}\,\Omega\cdot\mathrm{m}^2$, $D_{elec}=65$ nm |
| 2.2.2.2 | $M_s(T)$, $K_i(T)$, $\eta(T)$, $K_U^{\mathrm{eff}}(T)$ closed-form scalings | [`material_temperature.py`](../../src/vgsot_sim/material_temperature.py) | Bloch $T^{3/2}$ + modified Callen-Callen (exponent 2.18 default) |
| 2.2.2.1 | Self-heating RC dynamics — `q_mtj`, `q_sot`, `thermal_time_constant`, `steady_state_temperature`, `transient_temperature`, `self_heating_step` | [`thermal.py`](../../src/vgsot_sim/thermal.py) | Closed-form transient + per-step Euler update for use inside `dynamic_switching` loops |
| 2.2.2 (fig 2.13 / fig 2.12) | `sweep_material_temperature`, `tmr_voltage_sweep`, `thermal_transient` cases | [`temperature_cases.py`](../../src/vgsot_sim/temperature_cases.py) | Directly produces the data behind §2.2.2 figures; chapter-04 scripts are now thin wrappers around these |
| 2.2.2.6 (fig 2.14) | **Coupled self-heating in the LLG stepper** — `enable_self_heating=True` advances T(t), recomputes M_s(T), K_i(T) each step, threads them into `anisotropy.field()` | `time_series_cases.run_piecewise_direct_excitation` **and** `run_piecewise_terminal_voltage` + `dynamic_switching.switching(Ki_T=, Ms_T=)` | `SimResult.T_K / Ms_T / Ki_T` arrays record the diagnostics. Verified: closed-form T_eq agrees with stepper steady state to 0.01 K; M_s/K_i drifts agree with the analytic scaling laws. |
| 2.3.3 (fig 2.11) | **SER MC with self-heating + per-trial seeding** | `ser_cases.ser_sot_no_vcma_thermal(enable_self_heating=, T_ambient_K=, seed=)` | The `seed` argument deterministically derives per-trial seeds from `(seed, isot, trial_idx)`; the `enable_self_heating` flag threads through to the LLG stepper for temperature-corrected SER curves |
| 2.3.5 (fig 2.3-5) | **Built-in process-variability sweep** β_eff(CV_Δ) / F(CV_Δ) / η_c | [`ser_cases.variability_sweep(...)`](../../src/vgsot_sim/ser_cases.py) | Wraps `analysis.variability.transfer_function_F` so the whole §2.3.5 Brinkman/MC budget is regeneratable without the external `07_process_variability/variability_sim` script |
| 2.2.3.2 (fig 2.9 + technical_details §4.5) | **Cayley vector LLG integrator** with explicit `sigma_SH` 3-vector parameter | [`dynamic_switching_vector.py`](../../src/vgsot_sim/dynamic_switching_vector.py), `time_series_cases.run_piecewise_direct_excitation(integrator="cayley")` | Closed-form Cayley rotation: `m_{n+1} = m_n + (2s)/(1+s²|ω|²)·(ω×m_n + s ω×(ω×m_n))`. Norm-preserving to machine precision; converges to the spherical-Euler result as Δt → 0; reveals a small FL-SOT sign discrepancy in the Euler closed form near threshold |
| 2.2.4 (fig 2.9) | **Architecture-diagram generator** | [`demo/plot_architecture.py`](../../demo/plot_architecture.py) | Produces `Chapter02_local_07.png` for the chapter and the docs |
| 2.2.2.3 (toggle) | **Series resistance** `R_series` (contact / lead parasitic) | `PhysicalConstantsConfig.R_series` (default 0), `tmr(..., include_series=True)` | Off by default → byte-identical legacy behaviour; flip the config value to reproduce wafer-level R(V) that includes contact resistance |
| 2.2.3.3 (toggle) | **End-to-end `rng=…` plumbing** for byte-reproducible MC trials | `field()`, `switching()`, `run_piecewise_*` all accept an optional `rng=np.random.Generator` and forward it to `init()` + `stochastic()` | Legacy global `np.random.seed(...)` path unchanged when `rng=None` |
| 2.3 / 2.3.3 (display) | **Psw ↔ SER toggle** on `SerResult` | `SerResult.psw` property = `1 - ser`; figure scripts take `--metric=psw` (default) or `--metric=ser` | Chapter §2.3.3 reports P_sw directly, so the default switches to P_sw to align with experiment |
| 2.2.3.3 (toggle) | **`rng_mode` switch in `ser_sot_no_vcma_thermal`** | `rng_mode="legacy"` (default, `np.random.seed`) or `rng_mode="generator"` (`np.random.default_rng`, no global state) | Identical statistics, different bit-streams. Generator mode is multiprocessing-friendly. Per-trial seeds derived from a shared `_trial_seed(seed, i_sot, idx)` helper so the *stream of seeds* is the same in both modes. |
| 2.2.3.2 (toggle) | **Per-call `integrator` choice in `ser_sot_no_vcma_thermal`** | `integrator=None` (stepper default = **cayley**) / `"euler_spherical"` / `"cayley"` | **Default is now `"cayley"`** (the §2.2.3.2 method); θ_SH recalibrated to 0.066 against it after the FL-SOT fix. With the fix the two integrators share one RHS, so the prior "~10 %" euler-vs-Cayley note understated the true bug effect (~40 % on V_th). |
| 2.2 (regression) | **`tests/test_toggles.py`** covers R_series, rng plumbing, integrator switch, Psw alias, rng_mode validation, Cayley rng reproducibility, terminal-voltage integrator dispatch | [`tests/test_toggles.py`](../../tests/test_toggles.py) | 12 tests, all pass in 1.4 s |
| 2.2.3.2 (sym-symmetry) | **`run_piecewise_terminal_voltage` adds `integrator` / `sigma_SH` / `rng`** to match `run_piecewise_direct_excitation` | `time_series_cases.run_piecewise_terminal_voltage` | Cayley dispatch identical to direct_excitation; norm-preserving even under terminal-voltage piecewise driving |
| 2.2.3.2 (Cayley rng) | **rng= plumbed through Cayley path** | `switching_vector(rng=...)` → `field(rng=...)` → `stochastic(rng=...)` | Same `Generator(seed)` reproduces Cayley trajectory bit-for-bit; verified by `test_cayley_rng_byte_reproducible` |
| 2.2.2.7 (chapter doc) | **STT derivation §2.2.2.7** with Slonczewski form, LL conversion, STT↔SOT side-by-side table, and rationale for SOT-as-main-driver |  | Closes the §2.2.1.1 forward reference that earlier promised "STT 详细展开在 2.2.2 节给出" |
| 2.3.4 (calibration) | **θ_SH calibrated to V_th experiment** | `configs.PhysicalConstantsConfig.theta_SH=0.066` (effective; 0.04 pre-fix, 0.25–0.3 literature β-W); scan in `09_simulation_figures/calibrate_to_experiment.py` | Recalibrated (Cayley) against Device A P→AP `V_th(0.75 ns)=894 mV`: 50%-crossing `V_th^sim=0.895 V` (exact), jointly with `TMR=1.0`, `RA=16.6e-12 Ω·m²`. The corrected P_sw is broader (back-hopping plateau); near-threshold match preserved. `0.25–0.3` accepted for bare-material studies. |
| 2.3 (analysis tooling) | NB log-linear fit, Sigmoid fit + Wilson CI, Brinkman variability budget + D2D MC, binomial coverage + sampling sensitivity | [`analysis/`](../../src/vgsot_sim/analysis/) sub-package | Used by 06/07/08 chapter scripts in place of inline math |
| 2.4 (RC) | **Random-telegraph-noise (RTN) node primitive** — continuous-time 2-state Markov hopping with closed-form ⟨s⟩=tanh(ΔV/Vc0) and τ(V)=1/(r↑+r↓); rate floored at 1/τ0 (physical domain \|V\|<Vc0); exact propagator (V const within a step); `from_nb_fit` factory (ns units) | [`rtn/telegraph.py`](../../src/vgsot_sim/rtn/telegraph.py) | A single **candidate reservoir node** (not a reservoir: no input weights / coupling / readout); `tests/test_rtn_telegraph.py` (10). `V` is a phenomenological tilt, not `V_MTJ`. Complements (does not replace) the pulse-switching LLG. Bridge validation + reservoir layer: [`scripts/10_rtn_reservoir/`](../../scripts/10_rtn_reservoir/). |
| 2.2 (second engine) | **Full-dynamics LLG Verilog-A** — `m(t)` integrated internally (`ddt` state), per-formula aligned to `switching_vector`/`field`; PMA+VCMA+ellipsoid-demag+SOT(DL/FL)+STT; harness-driven thermal field (`hx,hy,hz`); compiles in OpenVAF, runs in ngspice | [`va/llg/vgsot_llg.va`](../../va/llg/vgsot_llg.va), [`va/README.md`](../../va/README.md) | Cross-engine: m_z(t) matches the Python engine to max ~0.006 over 0–3 ns (`va/llg/cross_validate_summary.json`); `tests/test_cross_engine_va.py` (skips without OpenVAF/ngspice) |

---

## Not implemented (modeled in thesis text, pending in code)

These items are described physically in Chapter 2 but do not yet have an executable counterpart in `vgsot-sim`. Calling the relevant cases today silently falls back to a constant value or an approximation.

### High priority (changes simulated output noticeably)

_(empty — all previously-tracked high-priority items moved to Implemented above.)_

### Medium priority

| Thesis section | Pending feature | Current behaviour | Needed |
|---|---|---|---|
| 2.2.1.1 | **Explicit spin polarisation vector** $\hat\sigma_{SH}$ (legacy spherical-Euler integrator only) | Available via the **Cayley vector integrator** (`integrator="cayley"`) which already takes `sigma_SH` at runtime. The spherical-Euler default still hard-codes $\hat\sigma = -\hat{x}$ (its FL-SOT `cos(theta)*cos(phi)` sign was corrected in 2026-06 and now matches Cayley for that fixed $\hat\sigma$; see `tests/test_integrator_consistency.py`). | Either fully retire the spherical expansion in favour of Cayley, or re-derive its closed form from Cartesian for arbitrary $\hat\sigma$ |

### Low priority / clean-up

| Thesis section | Pending | Notes |
|---|---|---|
| 2.3 (display) | Auto-detect `--metric` from saved SerResults | The toggle currently has to be passed on the command line; a future iteration could persist the metric inside `SerResult` and let downstream plot helpers pick it up automatically |

---

## Known caveats in implemented code

These are correct enough to use but worth being aware of when comparing against literature:

1. **`initialize.py` PAP–angle convention vs `tmr.py` resistance convention**. `initialize.py` sets `θ ≈ π` for `PAP=1` (labelled "Parallel"), whereas `tmr.py` returns $R_P$ at $m_z = +1$ (i.e. $\theta \approx 0$). The flat resistance at the start of a "PAP=1, parallel" simulation is therefore $R_{AP}$, not $R_P$. The two conventions cancel inside any single full simulation, but the labelling will surprise readers expecting one of the two conventions to dominate.

2. **SOT FL contribution to `dphi/dt` — sign corrected (2026-06)**. Hand-derivation against the standard Gilbert form $\tau_{DL} = -\gamma H_{DL} m\times(m\times\hat\sigma)$, $\tau_{FL} = +\gamma H_{FL} m\times\hat\sigma$ with $\hat\sigma = -\hat{x}$ (as commented in the upstream Verilog-A port, `m_s = (-1, 0, 0)`) matches the code for the H_eff, STT, and DL-SOT contributions to both `dtheta/dt` and `dphi/dt`, and for the FL-SOT contribution to `dtheta/dt`. The FL-SOT contribution to `dphi/dt` previously carried the wrong sign on the `cos(theta)*cos(phi)` term (equivalent to flipping $\hat\sigma_{FL}$ relative to $\hat\sigma_{DL}$), so the spherical-Euler azimuthal precession disagreed with the Cartesian/Cayley stepper by up to ~25%. It has been corrected to `+cos(theta)*cos(phi)`; both integrators now reduce to the SAME right-hand side, regression-tested in `tests/test_integrator_consistency.py` (dm/dt match to <1e-4 as dt→0, R_SOT_FL_DL=0.83). The remaining limitation of the spherical expansion is that it hard-codes $\hat\sigma = -\hat{x}$; for arbitrary $\hat\sigma_{SH}$ use the Cayley integrator, which takes `sigma_SH` as a runtime vector (see "Not implemented: Medium priority" — explicit spin polarisation vector).

3. **No magnetisation renormalisation**. Spherical coordinates analytically preserve $|\mathbf{m}|=1$, but the $1/\sin\theta$ in `dphi_dt` makes the integrator stiff near the poles. Long simulations or extreme parameter sweeps can accumulate error; see the Cayley integrator item in the high-priority gap list.

---

## How to use this document

When you find a discrepancy between thesis text and simulation output, first check whether the relevant effect is in the **Implemented** table. If it is in **Not implemented**, the thesis figure was likely produced by an analytic formula or by an external script — refer to the script paths under `07_process_variability/` and `08_sampling_effect/` rather than expecting the corresponding behaviour from `vgsot-sim`.

When you add a new physics feature to the code, please move its row from this file's **Not implemented** section to **Implemented** and update the corresponding thesis cross-reference.
