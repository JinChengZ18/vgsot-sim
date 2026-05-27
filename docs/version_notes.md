# Version Notes

This file records implementation-level changes that are useful for maintainers
and for reproducing older results. The README stays focused on installation and
day-to-day use.

## Recent corrections (2026-05)

This release fixes physics bugs inherited from the upstream Python port of the
Zhang et al. Verilog-A VGSOT-MTJ compact model, then adds temperature-coupled
stepping and the alternative Cayley integrator used by the chapter 2.2.3
discussion. If you migrate from a downstream fork, re-run prior SER and
Sigmoid-slope simulations because quantitative results will shift.

- **Thermal-field sampling now respects the FDT.** Previously `stochastic.py`
  returned `xi / |xi|`, a unit-length random direction, so `|H_TH|` was clamped
  to one value instead of being chi-square distributed. The corrected
  implementation draws three independent `N(0, 1)` samples. Therefore
  `E[|H_TH|^2]` matches the Brown 1963 / Garcia-Palacios theory and is about
  three times larger than before. Near-threshold `P_sw(V)` curves are
  measurably broader.
- **`H_ex` is now a configurable bias field.** The default is `-50 Oe` along
  `-y`, perpendicular to the implicit `sigma_SH = -x`. Direction matters:
  `H_ex` must be perpendicular to `sigma_SH`, or the SER curve collapses to a
  roughly `0.5` random-bit plateau. The chapter's "200 Oe along current
  direction" maps to `h_ex_y` in the simulator convention; see
  `docs/technical_details.md` section 2.7.
- **`TMR(V)` defaults to the PDK three-parameter form** from the Hikstor
  SOT-MRAM PDK extraction. Select with `constants.tmr_model = "pdk"` (default)
  or `"lorentzian"` for the previous behavior. PDK coefficients `a_tmr`,
  `b_tmr`, `c_tmr`, and `k_tmr` are exposed in the config.
- **Coupled self-heating LLG step.**
  `run_piecewise_direct_excitation(enable_self_heating=True)` and the `ser_*`
  cases now advance `T(t)` through the RC thermal network, recompute `M_s(T)`
  and `K_i(T)` each step, and feed the corrected material parameters back into
  `anisotropy.field()`. `SimResult.T_K`, `Ms_T`, and `Ki_T` carry diagnostics.
- **Cayley vector LLG integrator** in `dynamic_switching_vector.py` provides an
  exact-rotation Cartesian step with explicit `sigma_SH` 3-vector,
  machine-precision norm preservation, and no polar singularity. Opt in with
  `integrator="cayley"`.
- **Built-in process-variability case.** `ser_cases.variability_sweep` produces
  the chapter 2.3.5 Brinkman/MC budget in one call, replacing the external
  `07_process_variability` script.
- **`theta_SH`, `TMR_0`, and `RA` jointly calibrated to experiment**
  (2026-05-16). Defaults `theta_SH = 0.04`, `TMR = 1.0`, and
  `RA = 16.6 ohm um^2` give `R_P ~= 5 kohm`, `R_AP ~= 10 kohm`, and
  `V_th(0.75 ns) ~= 903 mV`, matching Device A P-to-AP detailed `P_sw` to
  within 1 percent. Literature beta-W `theta_SH = 0.25` is accepted as an
  override. See `docs/technical_details.md` section 2.3 and
  `scripts/09_simulation_figures/calibrate_to_experiment.py`.
- **New opt-in toggles** (2026-05-16): `R_series` parasitic resistance on
  `tmr()` (default `0`, byte-identical legacy behavior), full-pipeline `rng=...`
  plumbing through `field()`, `switching()`, and `run_piecewise_*` for
  byte-reproducible Monte Carlo, and a `Psw`/`SER` display flag
  (`SerResult.psw`, `--metric=psw|ser` on `plot_ser_mc.py`).
- **Case pruning + default realignment to §2.3.3 protocol** (2026-05-23). The
  CLI now exposes only the four cases that map to the same-batch Device A
  experiment: `terminal_voltage_control` (device-level T-circuit demo),
  `sot_only_constant_current`, `sot_switching_no_vcma`, and
  `ser_sot_no_vcma_thermal`. Removed: `vcma_assisted_switching_isot_sweep`,
  `vcma_assisted_switching_vmtj_sweep`, `optimized_vgsot_switching`,
  `ser_optimized_vgsot` (plus the `run_two_pulse_optimized` low-level kernel,
  `SerOptimizedResult` dataclass, and the matching configs). These came from
  the upstream Verilog-A port and had no counterpart in our measurements.
  Defaults of the remaining SOT-only cases are now aligned to the §2.3.3
  protocol: 0.75 ns write pulse (`sim_mid1_step=750`) + 3.25 ns relaxation
  (`sim_end_step=4000`), `i_sot_list` brackets the calibrated `I_th ≈ 1.1 mA`,
  `r_sot_fl_dl=0.83` unified across all three, `target_mz=1.0` /
  `failure_tol=0.2` for the SER case to match the chapter-figure protocol.
  Existing user code that explicitly instantiates the removed configs needs
  to migrate to `run_piecewise_direct_excitation` Python loops.
