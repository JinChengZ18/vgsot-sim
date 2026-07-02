# Version Notes

This file records implementation-level changes that are useful for maintainers
and for reproducing older results. The README stays focused on installation and
day-to-day use.

## 2026-06 — verification round, BDR/demag fixes, RTN primitive, Verilog-A engine

A correctness + over-claim audit drove this round. Quantitative LLG/SER results shift vs 2026-05 because of the FL-SOT fix; figures were re-run and `theta_SH` recalibrated 0.04 → 0.066 (below).

- **FL-SOT `dphi/dt` sign fixed** (`dynamic_switching.py`): the `cos θ·cos φ` term had the wrong sign, so the spherical-Euler integrator disagreed with the norm-preserving Cayley stepper by up to ~25% in azimuth. They now share one RHS (`tests/test_integrator_consistency.py`). **Consequence:** the SOT switching threshold rises (the bug had lowered V_th by ~40% — the 0.066/0.04 θ_SH ratio — not the ~10% the old euler-vs-Cayley note claimed). The effective `theta_SH` was recalibrated 0.04 → **0.066** against the Cayley integrator (now the default), restoring V_th(0.75 ns)=0.895 V (50% crossing exact). **Scientifically meaningful side effect:** the corrected P_sw(V) transition is somewhat *broader* with a more pronounced over-drive back-hopping plateau — the FL-SOT bug had been cosmetically *sharpening* the switching curve. The near-threshold operating-region match to the experimental Sigmoid is preserved; the deep-over-drive plateau (a genuine LLG feature) sits below the monotonic Sigmoid, as already documented.
- **`compute_Rp` is now the honest `R_P = RA/A1`** (`initialize.py`): the previous Brinkman-Dynes-Rowell-looking expression cancelled exactly to `RA/A1` and never depended on the barrier. A genuine barrier-dependent Simmons/BDR predictor `resistance_area_bdr` is provided separately (matches the measured RA at MgO `m* ≈ 0.3 m_e`).
- **`k_u_eff_of_T` uses the exact oblate-ellipsoid demag** (was thin-disk, inconsistent with the LLG field path).
- **FDT thermal-field amplitude tracks the self-heating temperature** `T(t)` via `field(T=...)`, and accepts an externally-injected field via `field(h_th_ext=...)` so the harness can drive the SAME Brown-1963 noise stream into both the Python engine and the Verilog-A engine.
- **New RTN node primitive** `vgsot_sim.rtn` — continuous-time two-state Markov hopping (⟨s⟩=tanh(ΔV/Vc0), τ(V) fading memory), a single *candidate reservoir node* (not a reservoir: no input weights / coupling / readout). Rate floored at the attempt frequency 1/τ0 so `|V|<Vc0` is the physical domain; `V` is a phenomenological tilt, not `V_MTJ`. Complements (does not replace) the pulse-switching LLG; deepening in `scripts/10_rtn_reservoir/`.
- **New Verilog-A engine** `va/llg/vgsot_llg.va` — a full-dynamics macrospin LLG that integrates `m(t)` internally on the open-source OpenVAF/ngspice toolchain, per-formula aligned to `switching_vector`/`field`, and cross-validated against the Python engine (`va/llg/cross_validate.py`; m_z(t) agrees to ~0.006). vgsot-sim is now a dual Python/Verilog-A codebase.
- **Backing matrix** `CLAIMS.md` maps every documented claim to its test or committed artifact. Test suite: 43 tests (`PYTHONPATH=src pytest tests/`).

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
  (2026-05-16; `theta_SH` updated 2026-06). Defaults `theta_SH = 0.066`
  (Cayley; was 0.04 pre-fix), `TMR = 1.0`, and `RA = 16.6 ohm um^2` give
  `R_P ~= 5 kohm`, `R_AP ~= 10 kohm`, and `V_th(0.75 ns) = 0.895 V` (50%
  crossing exact), matching Device A P-to-AP detailed `P_sw`. Literature
  beta-W `theta_SH = 0.25-0.3` is accepted as an override. See `docs/technical_details.md` section 2.3 and
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
