# 12 — Pulse-width transferability of the single-point θ_SH calibration (E1)

The chapter calibrates the effective spin-Hall angle θ_SH = 0.066 against a single experimental anchor: Device A, P→AP, t_w = 0.75 ns, V_th = 894 mV (Cayley integrator, self-heating ON). This experiment tests whether that one-point calibration *transfers* across the pulse widths the experiment actually measured, and whether the deep-barrier macrospin model (intrinsic Δ₀ ≈ 48.5 at the calibrated parameter set) spontaneously reproduces the *shallow* effective Néel-Brown parameters (Δ_pulse ≈ 4.9, b ≈ 175 mV, Table 2.8) when its simulated V_th(ln t_w) is inverted with the same `analysis.nb_fit` machinery used on the measured data.

A zero-drive free-evolution run at the same calibrated parameter set (|H_ex| = 50 Oe model default and 200 Oe experimental bias, 5 seeds × 2 µs each) bounds the actual retention of the simulated device: zero flips over 10 µs per field, where τ_ret = 135 ns would predict ~74, shows the Table-2.7 "retention time" is a linear-barrier extrapolation artifact of the nanosecond write regime, not a physical dwell time.

## Run

```bash
# One process per pulse width (parallel-friendly; ~20-45 min each at 200 trials x 9 points)
python run_e1.py --widths 0.75
python run_e1.py --widths 1
python run_e1.py --widths 1.5
python run_e1.py --widths 2
python run_e1.py --widths 3
python run_e1.py --widths 5

# Zero-drive free evolution (2 fields x 5 seeds x 2 us, ~15 min)
python run_e1.py --freerun-only

# Merge everything -> NB inversion, figures, e1_results.json + console table
python run_e1.py --analyze

# Tiny end-to-end sanity pass
python run_e1.py --smoke
```

All Monte-Carlo settings mirror fig 2.11's calibration workpoint: `integrator="cayley"`, self-heating ON, `pap=1`/`target_mz=1.0` (P→AP), 3.25 ns relaxation tail, deterministic per-trial seeds (`rng_mode="generator"`, base seed 20260719 + 100·t_w).

## Outputs

- `e1_w<width>.json` — per-width P_sw(V) scan + 4-parameter logistic V_th (with 1σ) + interpolated 50% crossing.
- `e1_freerun.json` — flip counts / m_z extrema per (field, seed).
- `e1_results.json` — merged summary: simulated (a, b, Δ_pulse, V_c0, τ_ret) vs experimental Table-2.7/2.8 references, per-width deviations.
- `e1_vth_vs_tw.png` — V_th vs t_w overlay (sim points vs experimental log-linear law); clean panel, no baked figure number, composition happens in PPT.
- `e1_psw_family.png` — the underlying P_sw(V) transition family.

Experimental reference constants in the script header use Table-2.7 precision (a = 0.795 V, b = 175 mV; the prose rounds to 0.79/0.18).
