# Reproducing the Chapter 2 figures and tables

This guide maps every figure and table in Chapter 2 (the sMTJ device-modelling,
integrator, process-variability, sampling, and RTN-reservoir study built on
`vgsot-sim`) to the script that generates it, the exact command, its input data,
and its rough runtime. It lets a reader regenerate the research conclusions
behind each result.

## Prerequisites

Install the package once in editable mode from the repository root, then every
command below is a plain `python <script>`:

```bash
pip install -e .
```

On Windows, prefix commands that print Chinese/Unicode with
`PYTHONIOENCODING=utf-8` (PowerShell: `$env:PYTHONIOENCODING="utf-8"`). If you
prefer not to install, replace `python` with `PYTHONPATH=src python` and run from
the repository root.

How figures are assembled: each script renders one or more clean panels (no
`(a)/(b)` letters, no figure numbers — see the figure conventions). The final
composite in [`article/figs/`](../article/figs/) is assembled from those panels
in a PowerPoint step that adds the panel letters. Running a script therefore
reproduces the **data and panels**; the published composite adds only cosmetic
labels on top.

## Figure map

| Figure | File | What it shows | Generator | Runtime |
|---|---|---|---|---|
| 图2.1 | `Chapter02_local_01.png` | 3-terminal SOT-sMTJ structure + T-equivalent circuit | `demo/plot_t_circuit.py` | <5 s |
| 图2.2 | `Chapter02_local_02.png` | Layered behavioural model for probabilistic computing | `demo/plot_behavioral_model_layers.py` | <5 s |
| 图2.3 | `Chapter02_local_03.png` | Thermal / non-ideality modelling overview | `demo/plot_thermal_nonideal_overview.py` | <5 s |
| 图2.4 | `Chapter02_local_04.png` | Self-heating temperature transients T(t) | `scripts/04_thermal_nonidealities/plot_thermal_transients.py` | <5 s |
| 图2.5 | `Chapter02_local_05.png` | M_s(T), K_i(T), η(T) and TMR(V) | `scripts/04_thermal_nonidealities/plot_material_params.py` | <5 s |
| 图2.6 | `Chapter02_local_06.png` | Self-heating ON/OFF switching trajectory | `scripts/04_thermal_nonidealities/plot_self_heating_compare.py` | ~seconds |
| 图2.7 | `Chapter02_local_07.png` | Norm-preserving Cayley integrator verification | `scripts/09_simulation_figures/fig_integrator_verification_panels.py` | assembly <10 s (needs audit data first) |
| 图2.8 | `Chapter02_local_08.png` | vgsot-sim three-layer architecture | `demo/plot_architecture.py` | <5 s |
| 图2.9 | `Chapter02_local_09.png` | Single m_z(t) / R_MTJ(t) switching event | `scripts/09_simulation_figures/plot_single_trajectory.py` | <10 s |
| 图2.10 | `Chapter02_local_10.png` | 3-D magnetisation trajectory on the unit sphere | `scripts/09_simulation_figures/plot_3d_trajectory.py` | <10 s |
| 图2.11 | `Chapter02_local_11.png` | Monte-Carlo P_sw(\|I_SOT\|) sweep at 0.75 ns | `scripts/09_simulation_figures/plot_ser_mc.py` | ~minutes |
| 图2.12 | `Chapter02_local_12.png` | Measurement platform + device photos | *(experimental — not code-generated)* | — |
| 图2.13 | `Chapter02_local_13.png` | Write characteristics + Néel–Brown joint P_sw | `scripts/06_psw_t_fitting/fit_from_loops.py` (figure 1) | <5 s |
| 图2.14 | `Chapter02_local_14.png` | Sigmoid measurement vs Néel–Brown extrapolation | `scripts/07_process_variability/sigmoid_fig.py` | <5 s |
| 图2.15 | `Chapter02_local_15.png` | Device-to-device Néel–Brown consistency | `scripts/06_psw_t_fitting/fit_from_loops.py` (figure 2) | <5 s |
| 图2.16 | `Chapter02_local_16.png` | Process-variability composite (CV(Δ), F(CV), β) | `scripts/07_process_variability/variability_sim.py` | ~5–15 s |
| 图2.17 | `Chapter02_local_17.png` | Macrospin process-variability MC P_sw(V) | `scripts/07_process_variability/macrospin_process_variability_mc.py` | ~16 h at the committed defaults |
| 图2.18 | `Chapter02_local_18.png` | MC sampling-size sensitivity of F̂(N, CV_Δ) | `scripts/08_sampling_effect/sampling_sensitivity_sim.py` | ~2–5 min |
| 图2.19 | `Chapter02_local_19.png` | Hardware Bernoulli-sampling reliability | `scripts/08_sampling_effect/hw_sampling_reliability.py` | ~2–4 min |
| 图2.20 | `Chapter02_local_20.png` | Low-barrier RTN node characteristics | `scripts/10_rtn_reservoir/plot_node_figs.py` → `scripts/build_ppt_figs.py` | ~minutes |
| 图2.21 | `Chapter02_local_21.png` | sLLG free-running validation of the RTN abstraction | `scripts/10_rtn_reservoir/plot_bridge_figs.py` → `scripts/build_ppt_figs.py` | ~minutes (sLLG dwell) |
| 图2.22 | `Chapter02_local_22.png` | RTN-node reservoir minimal validation | `scripts/10_rtn_reservoir/plot_reservoir_figs.py` → `scripts/build_ppt_figs.py` | ~seconds (mean-field) |

## §2.1–2.2.4 — device model, thermal effects, architecture (图2.1–2.11)

The schematic and self-contained-simulation figures need no external data:

```bash
# Schematics (demo/ scripts import no packages; run from anywhere)
python demo/plot_t_circuit.py                 # 图2.1
python demo/plot_behavioral_model_layers.py   # 图2.2
python demo/plot_thermal_nonideal_overview.py # 图2.3
python demo/plot_architecture.py              # 图2.8

# Self-contained vgsot_sim runs (thermal + trajectories)
python scripts/04_thermal_nonidealities/plot_thermal_transients.py     # 图2.4
python scripts/04_thermal_nonidealities/plot_material_params.py        # 图2.5
python scripts/04_thermal_nonidealities/plot_self_heating_compare.py   # 图2.6
python scripts/09_simulation_figures/plot_single_trajectory.py         # 图2.9
python scripts/09_simulation_figures/plot_3d_trajectory.py             # 图2.10
python scripts/09_simulation_figures/plot_ser_mc.py                    # 图2.11
```

`plot_ser_mc.py` accepts `--metric=psw|ser`, `--rng-mode=legacy|generator`, and
`--integrator=cayley|euler_spherical` (defaults reproduce the chapter figure).
图2.6 and 图2.9–2.11 run the LLG stepper with self-heating on and the calibrated
`theta_SH = 0.066` / Cayley defaults.

## §2.2.3.2 — Cayley integrator verification (图2.7)

图2.7 is assembled from precomputed audit results under `result/sec_2_2_3_2/`,
which is **gitignored** — regenerate it first. These runs are the integrator
audit ("显式-ω 几何格式" convergence, Boltzmann effective temperature, norm
preservation, pole/Euler behaviour, strong/weak stochastic order, θ_SH
provenance) that backs the section's claims:

```bash
python scripts/02_integrator/order_of_accuracy.py                       # A  deterministic order (fast)
python scripts/02_integrator/test_norm_preservation.py                  # C  |m|=1 preservation (fast)
python scripts/09_simulation_figures/c12_norm_stability_sweep.py        # D1 Euler norm blow-up (fast)
python scripts/11_pole_singularity/probe_pole_singularity.py            # D2 spherical pole guard (fast)
python scripts/09_simulation_figures/integrator_audit.py                # F  figure-provenance audit (fast)
python scripts/09_simulation_figures/audit_theta_sh.py --full --n-seeds 200  # F  θ_SH ensemble (minutes)
python scripts/11_strong_order/strong_order_brownian.py                 # E  strong order (minutes)
python scripts/10_weak_order/weak_order_audit.py                        # E  weak order (minutes)
python scripts/11_boltzmann_teff/teff_dtsweep.py --full                 # B  effective temperature (overnight)

# Then assemble the figure from the A/B audit data:
python scripts/09_simulation_figures/fig_integrator_verification_panels.py   # 图2.7
```

The panel-B effective-temperature sweep (`teff_dtsweep.py --full`) is the only
overnight run; a fast pilot (drop `--full`) suffices to see the T_eff/T ≈ 1
(Cayley) vs 0.5 (spherical-Euler) split, though the tight confidence intervals in
the published panel need the full 60 ns × 48-trajectory ensemble. A byte-level
regression pins the fast subset: `python -m pytest tests/test_integrator_order.py
tests/test_figure_provenance.py -q`.

## §2.3.3–2.3.5 — experimental fitting and process variability (图2.13–2.17)

图2.13 and 图2.15 are produced together by one run of `fit_from_loops.py`, which
reads seven committed R–V hysteresis-loop files:

```bash
python scripts/06_psw_t_fitting/fit_from_loops.py       # 图2.13 (fig 1) + 图2.15 (fig 2)
python scripts/07_process_variability/sigmoid_fig.py            # 图2.14
python scripts/07_process_variability/variability_sim.py       # 图2.16 (analytic MC, seed=42)
python scripts/07_process_variability/macrospin_process_variability_mc.py   # 图2.17 (macrospin LLG MC; defaults 16 devices x 32 trials, Cayley, ~16 h)
```

Experimental data provenance: the raw same-batch wafer measurements are committed
under [`scripts/05_experimental_raw_data/raw/`](../scripts/05_experimental_raw_data/raw/)
(`device2_paperB/`, `device4_paperA/`; the space-named originals). `fit_from_loops.py`
reads ASCII-renamed copies of the seven full-loop files that live in
`scripts/07_process_variability/` (`device{2,4}pulse_width_*_ns_200_Oe.txt`); the
two sets are byte-identical loop data. 图2.14 (`sigmoid_fig.py`) hard-codes its
same-batch 0.75 ns P_sw points inline. 图2.17 runs a fresh macrospin ensemble
(`--devices`, `--thermal-trials`, `--seed` tune cost; self-heating on by default,
which is the state SER/P_sw thresholds must be reported in). Its defaults are the
committed figure's settings — 16 devices × 32 trials (512 pooled, Wilson 95%
half-width ±4.3% at p = 0.5 for the nominal curve) on the Cayley kernel; the
committed run took 15.6 h in one process, and cost scales as 2 × devices ×
trials — the nominal baseline is run at the same pooled trial count as the whole
mismatch ensemble, so `--devices` sets the nominal curve's statistics as well. It also writes `macrospin_variability_summary.json`
(crossings, fit-free 0.25→0.50 spans, logistic slopes under two baseline
conventions, per-device crossings) next to the figure. It writes only into its
own folder, so copy `Chapter02_local_17.png` to `article/figs/` after a
regeneration. The stdout of the committed run is kept as
`macrospin_variability_run_16x32.log` (whitelisted in `.gitignore`) because the
chapter quotes numbers derived from its tables.

### §2.3.3–2.3.4 verification runs (footnotes + 图2.13(d) panel)

Two reviewer-round verification experiments back the quantitative claims in
the §2.3.3 `note-tauret-freerun` and §2.3.4 `note-eta-sllg` footnotes and the
pulse-width-transfer panel offered for 图2.13(d):

```bash
# E1 — pulse-width transfer of the single-point θ_SH calibration + zero-drive
# retention bound (per-width shards parallelize; 3/5 ns need the re-centered
# _hi supplements; ~30-60 min per shard)
python scripts/12_pulse_width_transfer/run_e1.py --widths 0.75   # ... 1, 1.5, 2, 3, 5
python scripts/12_pulse_width_transfer/run_e1.py --widths 3 --center-ua 973 --suffix _hi
python scripts/12_pulse_width_transfer/run_e1.py --widths 5 --center-ua 877 --suffix _hi
python scripts/12_pulse_width_transfer/run_e1.py --freerun-only
python scripts/12_pulse_width_transfer/run_e1.py --analyze       # e1_results.json + panels

# E2 — measured-Sigmoid overlay and η_c slope-ladder decomposition
python scripts/13_sigmoid_overlay/run_e2.py --dense-shard 0      # ... shards 1-7
python scripts/13_sigmoid_overlay/run_e2.py --overdrive
python scripts/13_sigmoid_overlay/run_e2.py --relaxwin 3.25      # ... 10, 50
python scripts/13_sigmoid_overlay/run_e2.py --analyze            # e2_results.json + panels
```

Both use the fig 2.11 calibration workpoint (Cayley integrator, self-heating
ON, deterministic per-trial seeds); the e-series result JSONs are gitignored
and regenerate from these commands. See the two scripts' READMEs for the
acceptance criteria and panel descriptions.

## §2.3.6 — Monte-Carlo sampling budget (图2.18–2.19)

Both are self-contained Monte-Carlo studies with fixed seeds:

```bash
python scripts/08_sampling_effect/sampling_sensitivity_sim.py   # 图2.18 (~2–5 min)
python scripts/08_sampling_effect/hw_sampling_reliability.py    # 图2.19 (~2–4 min)
```

## §2.4 — RTN-node reservoir (图2.20–2.22)

These figures are assembled by `scripts/build_ppt_figs.py` from panel stems
`ch02_19/20/21`, which the plotting scripts emit (note the +1 offset: stem
`ch02_19` → 图2.20, from the 图2.7 insertion). The plotting scripts consume
committed result JSON produced by the compute scripts; see
[`scripts/10_rtn_reservoir/README.md`](../scripts/10_rtn_reservoir/README.md) for
the full sequence and expected runtimes.

```bash
# (Re)compute the reservoir/bridge results, then render panels:
python scripts/10_rtn_reservoir/validate_bridge.py        # -> bridge_results.json  (~minutes: free-running sLLG dwell)
python scripts/10_rtn_reservoir/benchmark_reservoir.py    # -> reservoir_results.json  (mean-field, fast)
python scripts/10_rtn_reservoir/plot_node_figs.py         # ch02_19_{a,b,c}  -> 图2.20
python scripts/10_rtn_reservoir/plot_bridge_figs.py       # ch02_20_{a,b,c}  -> 图2.21
python scripts/10_rtn_reservoir/plot_reservoir_figs.py    # ch02_21_{a,b,c}  -> 图2.22
python scripts/build_ppt_figs.py                          # assemble _20/_21/_22
```

## Tables

Most tables carry numbers emitted to **stdout** by the same script that renders
the neighbouring figure; a few are input-parameter or measured-data tables that
are sourced rather than regenerated. Prerequisites are the same as above
(`pip install -e .`, then run from the repository root).

### Computed tables

| Table | Content | Reproduce |
|---|---|---|
| 表2.6 | Four Sigmoid fits at t_w = 0.75 ns | `scripts/07_process_variability/sigmoid_fig.py` → stdout (same run as 图2.14) |
| 表2.7 | Device A model-independent observables | `scripts/06_psw_t_fitting/fit_from_loops.py` → stdout (same run as 图2.13 / 图2.15) |
| 表2.8 | Inverted Néel–Brown parameters (Δ, V_c0; τ₀ = 1 ns) | `scripts/06_psw_t_fitting/fit_from_loops.py` → stdout |
| 表2.9 | Néel–Brown extrapolation vs Sigmoid measurement | `scripts/06_psw_t_fitting/fit_from_loops.py` + `scripts/07_process_variability/sigmoid_fig.py` → stdout |
| 表2.10 | CV(Δ) variance-budget decomposition | `scripts/07_process_variability/variability_sim.py` → stdout |
| 表2.11 | Process margin: array β^eff retention vs CV_Δ | `scripts/07_process_variability/variability_sim.py` → stdout |
| 表2.12 | Estimator performance over (CV_Δ, N), R = 40, N_ref = 50000 | `scripts/08_sampling_effect/sampling_sensitivity_sim.py` → stdout (same run as 图2.18) |
| 表2.13 | K_req: exact vs Monte-Carlo vs CLT | `scripts/08_sampling_effect/hw_sampling_reliability.py` → stdout (same run as 图2.19) |
| 表2.14 | Low-barrier sLLG dwell statistics + attempt time τ₀ | `scripts/10_rtn_reservoir/validate_bridge.py` → `scripts/10_rtn_reservoir/bridge_results.json` |
| 表2.15 | Reservoir benchmark capacity (three configurations) | `scripts/10_rtn_reservoir/benchmark_reservoir.py` → `scripts/10_rtn_reservoir/reservoir_results.json` |

The 表2.14 / 表2.15 numbers are read from the committed result JSON
(`bridge_results.json`, `reservoir_results.json`); re-run the compute script to
regenerate them (`validate_bridge.py` runs the free-running sLLG dwell, ~minutes;
the reservoir mean-field benchmark is fast).

### Input and measured tables

These are not regenerated by a run — they are the sourced inputs and measurements
the model consumes:

| Table | Content | Source |
|---|---|---|
| 表2.1 | STT vs SOT driving-mechanism comparison | Qualitative comparison stated in the chapter text — no data source |
| 表2.2 | SOT-channel transport, thermal, and VCMA parameters | Literature values; the ones the simulator uses are the `PhysicalConstantsConfig` defaults in `src/vgsot_sim/configs.py` (VCMA β) plus the thermal parameters in `src/vgsot_sim/thermal.py` |
| 表2.3 | 80 nm SOT-MTJ core magnetic and dimensional parameters | `PhysicalConstantsConfig` defaults in `src/vgsot_sim/configs.py` (M_s, K_i, D_elec, t_f, …) |
| 表2.4 | TMR bias-decay PDK coefficients + port-level calibration | PDK coefficients (`k/a/b/c_tmr`) and the calibrated `TMR`, `RA`, `theta_SH` in `configs.py`; the θ_SH / RA calibration is reproduced by `scripts/09_simulation_figures/calibrate_to_experiment.py` |
| 表2.5 | Device-array test statistics (batch means) | Experimental same-batch wafer measurements; raw data under `scripts/05_experimental_raw_data/` |

## Caveats

- **图2.12 is not code-generated.** It is the experimental test-platform block
  diagram and device photographs; there is no script for it.
- **图2.7 needs `result/sec_2_2_3_2/` first.** That directory is gitignored, so a
  fresh clone must run the §2.2.3.2 audit commands above before assembling the
  figure.
- **Long runs:** 图2.7 panel B and 图2.17 (overnight; 图2.17 is ~16 h at the
  committed defaults); 图2.18–2.19, the §2.2.3.2 audit runs,
  and the §2.4 free-running sLLG dwell (`validate_bridge.py`) take minutes. All
  other results render in seconds.
- Scripts render panels into their own folder; the published composites in
  `article/figs/` add panel `(a)/(b)` letters via the PowerPoint assembly step.
