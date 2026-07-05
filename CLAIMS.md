# CLAIMS — backing matrix

Every substantive physics/engineering claim in the docs/code is mapped here to the test or committed artifact that backs it. Goal: **no unbacked assertion**. `backed` = a passing test or committed numeric artifact reproduces it; `gap` = asserted but not yet pinned by a dedicated check; `pending` = blocked on an open project item.

Run the test suite: `PYTHONPATH=src python -m pytest tests/ -q` (72 tests).

## Physics kernels (Python)

| Claim | Backing | Status |
|---|---|---|
| FDT thermal-field amplitude is Brown-1963 `sqrt(2 α kB T /(μ0 Ms γ V dt))`, scales as √T | `tests/test_fdt_temperature.py::test_fdt_amplitude_scales_with_sqrt_T` | backed |
| FDT noise temperature defaults to `constants.T`, tracks self-heating T(t) when supplied | `tests/test_fdt_temperature.py::test_fdt_default_T_is_constants_T` + threading in `time_series_cases.py` | backed |
| Spherical-Euler and Cayley integrators reduce to the SAME dm/dt (FL-SOT sign fixed) | `tests/test_integrator_consistency.py` (2) | backed |
| Cayley step preserves \|m\|=1 to ~1e-16 | `tests/test_toggles.py::test_integrator_toggle_runs`, `::test_terminal_voltage_integrator_runs` | backed |
| `k_u_eff_of_T` default demag is the exact oblate ellipsoid (not thin-disk) | `tests/test_kueff_demag.py` (2) | backed |
| `R_P = RA / A1` (honest; was a cancelling BDR expression) | `tests/test_bdr_resistance.py::test_compute_Rp_is_RA_over_A1`, `::test_compute_Rp_independent_of_barrier` | backed |
| Simmons/BDR predictor `resistance_area_bdr` genuinely depends on φ_ox, t_ox; matches measured RA at MgO m*≈0.3 | `tests/test_bdr_resistance.py::test_bdr_predictor_depends_on_barrier`, `::..._matches_measured_RA...` | backed |
| TMR(V) PDK 3-parameter (default) and Lorentzian forms; `R_series` additive | `tests/test_tmr_voltage.py` (4) + `tests/test_toggles.py::test_r_series_*` | backed |
| Self-heating RC couples M_s(T)/K_i(T) into the field each step | `IMPLEMENTATION_STATUS.md` states "verified to 0.01 K"; **no committed pytest** | gap (low) |
| Process-variability CV(Δ) budget partitions total variance; macrospin sample reproducible | `tests/test_variability.py` (2) | backed |
| End-to-end `rng=` byte-reproducibility (euler + cayley); `rng_mode` switch; Psw=1−SER alias | `tests/test_toggles.py` (rng/mode/psw tests) | backed |

## §2.2.3.2 numerical-methods audit (this work)

> Every strong §2.2.3.2 claim empirically pinned. The shipped `switching_vector` is the EXPLICIT-ω Cayley step (ω frozen at m_n), NOT the implicit midpoint the prose implied; the true midpoint is now selectable via `integrator="cayley_midpoint"` (`n_midpoint>0`). "artifact" = committed reproducible script; its numeric output lands in the gitignored `result/sec_2_2_3_2/`.

| Claim | Backing | Status |
|---|---|---|
| Norm \|m\|=1 preserved to ≤2e-16 over dt∈[1e-15,1e6]×\|ω\|∈[1,1e14]; the per-step renorm is decorative (raw 0.75 ns trajectory drift ≤1.9e-15) | `scripts/02_integrator/test_norm_preservation.py` | backed (artifact) |
| Published explicit-ω step is GLOBALLY FIRST-order (p=1.01); `cayley_midpoint` recovers SECOND order (p=2.00); a constant-ω control isolates the loss to ω(m(t)) time-variation, not damping/renorm | `tests/test_integrator_order.py`, `tests/test_integrator_midpoint.py` (10) | backed |
| Strong stochastic order = 0.97 [0.95,1.00] (GBM-gated harness, common-Brownian-path coupling, independent midpoint reference) | `scripts/11_strong_order/strong_order_brownian.py` | backed (artifact) |
| Weak stochastic order = 1.98–2.04 (self-Cayley + midpoint refs); survives a transverse symmetry-breaking field (pooled ≥1.8) | `scripts/10_weak_order/weak_order_audit.py`, `weak_order_symbreak.py` | backed (artifact) |
| Explicit-ω Cayley samples the correct Boltzmann law: T_eff/T = 0.979±0.014 (60 ns, Δt→0 extrap); true-midpoint = 0.979 (⇒ no explicit-ω / renorm stationary bias) | `scripts/11_boltzmann_teff/teff_dtsweep.py`, `merge_full_results.py` | backed (artifact) |
| Spherical-Euler equilibrates at EXACTLY T/2 (0.49, dt-independent) — missing Wong–Zakai drift D·cot θ; adding the term restores 0.94; RHS-consistency tests cannot catch this ensemble-limit defect | `scripts/11_boltzmann_teff/euler_drift_discrimination.py` | backed (artifact) |
| Figure provenance: fig 2.10 = cayley (driver default); fig 2.11 (SER) = euler by default but drive-dominated, so euler→cayley shifts I_50 only −0.5% (width ×1.16, euler sharper). runtime `theta_SH=0.066`: write-condition (self-heating ON) I_50≈1140 µA matches exp 1152 µA to ~1%; self-heating OFF ≈1300 µA (+13%). NB the ~12%/OFF figure alone would misstate the calibration | `scripts/09_simulation_figures/integrator_audit.py`, `audit_theta_sh.py`, `tests/test_figure_provenance.py` | backed |

## RTN node primitive (candidate reservoir node, `vgsot_sim.rtn`)

> Scope: a SINGLE free-running 2-state node, not a reservoir (no input weights, node
> coupling, or trained readout). `V` is a phenomenological tilt valid for `|V|<Vc0`,
> not literally `V_MTJ`. Building an actual reservoir + calibrating `V` to a device
> drive are in [`scripts/10_rtn_reservoir/`](scripts/10_rtn_reservoir/).

| Claim | Backing | Status |
|---|---|---|
| ⟨s⟩_inf = tanh(Δ V/Vc0) = (r↑−r↓)/(r↑+r↓) on \|V\|<Vc0 | `tests/test_rtn_telegraph.py::test_stationary_mean_equals_tanh_and_rate_ratio` | backed |
| τ(V)=1/(r↑+r↓) peaks at V=0 = τ0·exp(Δ)/2 | `::test_tau_peaks_at_zero_bias`, `::test_relaxation_time_is_inverse_total_rate` | backed |
| Exact 2-state propagator recovers the stationary mean; binary, seed-reproducible | `::test_exact_propagator_recovers_stationary_mean`, `::test_states_are_binary_and_seed_reproducible` | backed |
| `TelegraphParams.from_nb_fit` wires vgsot's own NB inversion (ns units) | `::test_from_nb_fit_copies_delta_vc0`, `::test_simulate_trace_shape_and_values` | backed |
| Escape rate floored at attempt freq 1/τ0; equals `nb_fit.psw_nb` implied rate; `step()` warns for \|V\|>Vc0 (physical domain `\|V\|<Vc0`) | `::test_rate_clipped_at_attempt_frequency`, `::test_rate_matches_psw_nb_instantaneous`, `::test_step_warns_outside_domain` | backed |

## RTN bridge + reservoir (stages 2–3, `vgsot_sim.rtn.bridge` / `reservoir`)

> Quantitative validation results in [`scripts/10_rtn_reservoir/README.md`](scripts/10_rtn_reservoir/) (+ committed `*_results.json`).

| Claim | Backing | Status |
|---|---|---|
| `ki_for_delta` keeps the well perpendicular (K_u^eff>0 incl. demag); PMA-only Ki goes in-plane | `tests/test_rtn_bridge.py::test_ki_delta_roundtrip`, `::test_low_barrier_is_perpendicular_not_in_plane` | backed |
| `bias→V` tilt slope = μ0·Ms·v/(kB·T); free-run / dwell / PSD helpers run | `::test_tilt_per_field_value`, `::test_free_run_shape_and_range`, `::test_dwell_times_*`, `::test_psd_lorentzian_runs_on_synthetic` | backed |
| Low-Δ free-running sLLG → exponential dwell (CV≈0.99), τ0≈15–27 ns (≠1 ns), tanh-shape ⟨m_z⟩ amplitude-compressed A≈0.7 | numbers quoted in `scripts/10_rtn_reservoir/README.md`; regenerate deterministically (fixed seeds) via `validate_bridge.py` (`*_results.json` is gitignored as reproducible output) | backed (script + README) |
| Heterogeneous W_in reservoir MC≈8 ≫ broadcast-identical baseline MC≈0.6 (ridge readout recovers a linear map) | `tests/test_rtn_reservoir.py::test_heterogeneous_beats_broadcast_baseline`, `::test_ridge_recovers_linear_map`, `::test_memory_capacity_decays_with_delay` | backed |
| NARMA-10 NRMSE≈0.55; stochastic single-device MC≈0.36 (needs replica-averaging) | numbers quoted in `scripts/10_rtn_reservoir/README.md`; regenerate deterministically (fixed seeds) via `benchmark_reservoir.py` (`*_results.json` is gitignored as reproducible output) | backed (script + README) |
| Ring delay-line (`ring_reservoir`) beats the filter-bank MC ceiling (MC 19→37 @ n=25→200, grows with n); random coupling lowers MC at every tested (radius, scale) | `tests/test_rtn_reservoir.py::test_ring_reservoir_beats_filter_bank_mc` + committed fig `article/figs/Chapter02_local_21.png` (panel script recomputes with fixed seeds) | backed |
| Ring trade-off: IPC deg-2 = 0, parity at chance; Mackey-Glass ranking reverses the MC ranking (qualitative only — h=84 magnitudes numerically fragile, Gram cond 1e14–1e20) | `scripts/10_rtn_reservoir/README.md` (D1/D5 sections incl. conditioning caveat); `::test_mackey_glass_series_and_task` (smoke) | backed (documented) |
| Kernel/generalization rank + ESP convergence; device-budget optimum is few-nodes × deep-averaging at every tested B (64/256/1024 → n=8) | `::test_kernel_quality_and_esp`; `scripts/10_rtn_reservoir/budget_allocation.json` (committed run) | backed |

## Verilog-A engine (`va/llg/vgsot_llg.va`)

| Claim | Backing | Status |
|---|---|---|
| Full-LLG `.va` compiles (OpenVAF) and integrates m(t) in ngspice with \|m\|≈1 | `va/llg/tb_switch.spice` run (\|m\| drift ~1e-4) | backed (manual) |
| `.va` m_z(t) matches the Python `switching_vector` engine | demonstrated: max \|Δm_z\|≈0.006 over 0–3 ns, exact at equilibrium | partial (`pending`: committed regression + pytest = project #15) |
| Stochastic MC via harness-driven `h_th` (same seed → Python==VA) | injection path present; seeded-equivalence harness | pending (#15) |

## Consuming project (`smtj_pbnn_sim/eda/`) — cross-references

| Claim | Backing | Status |
|---|---|---|
| `smtj_sot.va` OSDI↔numpy regression R²=1.0, max\|err\|=3.5e-4 (tool self-consistency) | `eda/testbenches/regression_summary.json` (actually executed) | backed |
| Behavioural sigmoid vs measured Device-A P→AP R²=0.992 | `eda/testbenches/golden_summary.json` | backed |
| `.va` params ↔ device formulas ↔ committed golden stay consistent | `tests/test_golden_pins_device_source.py` (2) | backed |
| LLG ↔ behavioural threshold cross-validation (0.25 mV) | `eda/testbenches/llg_validate_summary.json` | **pending**: STALE after the FL-SOT fix; will not reproduce until θ_SH recalibration (#12) |
| `θ_SH=0.04` gives V_th(0.75 ns)≈894 mV (~1%) | was reproduced by the **buggy** euler; **invalidated** by the FL-SOT fix → needs recal (#12) | pending |

## Clarifications (residual doc reconciles)

- **γ folds μ0.** Code `gamma = 2 u0 uB / ħ ≈ 2.21e5 m/(A·s)` (= μ0·g·μB/ħ, g=2). Doc formulas for `H_SOT`/FDT that omit an explicit μ0 use fields in A/m where γ absorbs μ0 — dimensionally self-consistent, not a discrepancy.
- **H_ex default.** `h_ex_y = −50·1000/(4π) ≈ −3979 A/m` = −50 Oe along −y, **perpendicular** to σ_SH=−x̂ (required for deterministic SOT switching). The chapter "200 Oe along the current" maps to `h_ex_y`; the 50-vs-200 Oe magnitude is a documented simulator-tuned value (`configs.py` §5). Sign/units are internally consistent.
- **Two resistance forms.** `initialize.py` sets the starting resistance via a fixed-TMR R(θ); `tmr.tmr()` uses the bias-dependent TMR_eff(V). Both coexist (documented caveat #1) and cancel within a single simulation; unifying on `tmr.tmr()` is a low-priority cleanup.
- **`compute_Rp` BDR.** Corrected: it is now the honest `RA/A1`; the genuine barrier-physics predictor is `resistance_area_bdr` (see initialize.py). The old "~10% Cayley-vs-euler threshold" doc note is **wrong** (the bug-vs-correct gap is ~43%) and will be removed when #12 recalibrates.
