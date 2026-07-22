# 13 — Measured-Sigmoid overlay and η_c decomposition (E2)

Fig 2.11 aligns simulation and experiment at the 50% threshold only. This experiment overlays the full measured Device-A AP→P Sigmoid (fig 2.14b data, mirrored from `../07_process_variability/sigmoid_fig.py`) on a densely sampled simulated transition (13 points × 1000 pooled trials), extracts β_s^sLLG with a 95% CI, and places the macrospin engine on the slope ladder β^NB (7.9 V⁻¹) → β_s^sLLG → β^meas (44.6 V⁻¹) — i.e. how much of the C2C narrowing factor η_c ≈ 5.6 the macrospin dynamics explains.

It also confronts the simulated over-drive back-hopping plateau (fig 2.11, ~0.8) with the measured clean saturation to 1.0 by ~1020 mV, and separates judgement convention from device physics by scanning the post-pulse relaxation window (3.25 / 10 / 50 ns at 1400 and 2000 µA): if the plateau height moves with the window, the "plateau" is a readout-timing artifact of the fixed-window m_z criterion rather than device physics.

## Run

```bash
# Transition scan sharded by seed block (pool = 8 x 125 = 1000 trials/point)
python run_e2.py --dense-shard 0   # ... repeat for shards 1-7

python run_e2.py --overdrive          # 1250-1600 uA, 200 trials/point
python run_e2.py --relaxwin 3.25      # window sensitivity at 1400/2000 uA
python run_e2.py --relaxwin 10
python run_e2.py --relaxwin 50

python run_e2.py --analyze            # pool -> fits + e2_overlay.png + e2_relaxwin.png + e2_results.json
```

Settings mirror the fig 2.11 calibration workpoint: Cayley integrator, self-heating ON, AP→P (`pap=1`, `target_mz=1.0`), 0.75 ns pulse, deterministic per-trial seeds.

## Outputs

- `e2_dense_shard<k>.json`, `e2_overdrive.json`, `e2_relaxwin_<w>.json` — raw scans.
- `e2_results.json` — pooled fit, slope ladder (β^NB / β^sLLG / β^meas, η_c shares), overdrive + window data.
- `e2_overlay.png` — measured points/fit vs simulated points/fit vs uncorrected NB curve; clean panel, composition in PPT.
- `e2_relaxwin.png` — plateau height vs judgement window.
