# Reproducing the Chapter 2 figures

This guide maps every figure in Chapter 2 (the sMTJ device-modelling, integrator,
process-variability, sampling, and RTN-reservoir study built on `vgsot-sim`) to
the script that generates it, the exact command, its input data, and its rough
runtime. It lets a reader regenerate the research conclusions behind each figure.

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
| 图2.17 | `Chapter02_local_17.png` | Macrospin process-variability MC P_sw(V) | `scripts/07_process_variability/macrospin_process_variability_mc.py` | minutes |
| 图2.18 | `Chapter02_local_18.png` | MC sampling-size sensitivity of F̂(N, CV_Δ) | `scripts/08_sampling_effect/sampling_sensitivity_sim.py` | ~2–5 min |
| 图2.19 | `Chapter02_local_19.png` | Hardware Bernoulli-sampling reliability | `scripts/08_sampling_effect/hw_sampling_reliability.py` | ~2–4 min |
| 图2.20 | `Chapter02_local_20.png` | Low-barrier RTN node characteristics | `scripts/10_rtn_reservoir/plot_node_figs.py` → `scripts/build_ppt_figs.py` | ~minutes |
| 图2.21 | `Chapter02_local_21.png` | sLLG free-running validation of the RTN abstraction | `scripts/10_rtn_reservoir/plot_bridge_figs.py` → `scripts/build_ppt_figs.py` | long (LLG dwell) |
| 图2.22 | `Chapter02_local_22.png` | RTN-node reservoir minimal validation | `scripts/10_rtn_reservoir/plot_reservoir_figs.py` → `scripts/build_ppt_figs.py` | long (benchmarks) |

## §2.1–2.2.4 — device model, thermal effects, architecture (图2.1–2.11)

The schematic and self-contained-simulation figures need no external data:

```bash
# Schematics (demo/ scripts import no packages; run from anywhere)
python demo/plot_t_circuit.py                 # 图2.1
python demo/plot_behavioral_model_layers.py   # 图2.2
python demo/plot_thermal_nonideal_overview.py # 图2.3
python demo/plot_architecture.py              # 图2.8  (see caveat below)

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
python scripts/07_process_variability/macrospin_process_variability_mc.py   # 图2.17 (macrospin LLG MC)
```

Experimental data provenance: the raw same-batch wafer measurements are committed
under [`scripts/05_experimental_raw_data/raw/`](../scripts/05_experimental_raw_data/raw/)
(`device2_paperB/`, `device4_paperA/`; the space-named originals). `fit_from_loops.py`
reads ASCII-renamed copies of the seven full-loop files that live in
`scripts/07_process_variability/` (`device{2,4}pulse_width_*_ns_200_Oe.txt`); the
two sets are byte-identical loop data. 图2.14 (`sigmoid_fig.py`) hard-codes its
same-batch 0.75 ns P_sw points inline. 图2.17 runs a fresh macrospin ensemble
(`--devices`, `--thermal-trials`, `--seed` tune cost; self-heating on by default,
which is the state SER/P_sw thresholds must be reported in).

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
python scripts/10_rtn_reservoir/validate_bridge.py        # -> bridge_results.json
python scripts/10_rtn_reservoir/benchmark_reservoir.py    # -> reservoir_results.json  (slow)
python scripts/10_rtn_reservoir/plot_node_figs.py         # ch02_19_{a,b,c}  -> 图2.20
python scripts/10_rtn_reservoir/plot_bridge_figs.py       # ch02_20_{a,b,c}  -> 图2.21
python scripts/10_rtn_reservoir/plot_reservoir_figs.py    # ch02_21_{a,b,c}  -> 图2.22
python scripts/build_ppt_figs.py                          # assemble _20/_21/_22
```

## Caveats

- **图2.12 is not code-generated.** It is the experimental test-platform block
  diagram and device photographs; there is no script for it.
- **`plot_architecture.py` writes `Chapter02_local_07.png`, not `_08`.** Its
  output number is stale from before the 图2.7 integrator figure was inserted
  (which shifted architecture to 图2.8). The canonical `article/figs/Chapter02_local_08.png`
  is produced in the compositing step; regenerate the architecture panel from the
  script but expect the `_07` filename locally.
- **图2.7 needs `result/sec_2_2_3_2/` first.** That directory is gitignored, so a
  fresh clone must run the §2.2.3.2 audit commands above before assembling the
  figure.
- **Long runs:** 图2.7 panel B (overnight), 图2.17–2.19 (minutes), and the §2.4
  reservoir benchmarks (slow). All other figures render in seconds.
- Scripts render panels into their own folder; the published composites in
  `article/figs/` add panel `(a)/(b)` letters via the PowerPoint assembly step.
