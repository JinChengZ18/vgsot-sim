# 14 — VCMA-branch activation (E3)

The platform is named vgsot-sim, §2.1.3 builds the SOT-VCMA joint drive model, and §2.2.4 promises a dual-pulse scenario "used to quantify the energy-optimization benefit" — yet no chapter figure ever exercised the VCMA path (`vnv=1`). This experiment runs it in two configurations at the fig 2.11 calibration workpoint (0.75 ns, P→AP, Cayley, self-heating ON):

1. **Simultaneous assist** — V_MTJ held during the SOT write stage only (stage layout: idle/precharge → write → relax). Output: I_th(V_MTJ) for V ∈ {−0.8, −0.4, 0, +0.4, +0.8} V against the analytic static scaling I_c0(V)/I_c0(0) = 1 − ΔK_v/K_U^eff with ΔK_v = β_VCMA·V/(t_ox·t_f).
2. **Sequential dual pulse** — the §2.2.4 preset's literal description: a 1 ns VCMA-only pre-pulse (+0.8 V, I_SOT = 0) followed by an SOT-only write. Linear VCMA has no persistence, so this configuration is an honest null test of the preset's premise.

Write energy at the p = 0.5 workpoint is compared in `--analyze` from the extracted thresholds: E = I_th²·R_W·t_w + V_MTJ²/R_P·t_V (R_P = RA/A_elec ≈ 5 kΩ).

Scan windows are wide and sign-dependent, deliberately **not** centered on the analytic scaling — whether the dynamic 0.75 ns threshold follows the static K_U^eff scaling is what the experiment measures. At +0.8 V the transition collapses into a stochastic band (P_sw ≈ 0.6–0.9 across 0.25–1.23 mA with no clean 50% crossing): strong bias removes so much barrier that the deterministic write window disappears, so VCMA assist has an optimum bias range rather than a monotone benefit.

## Run

```bash
# One detached process per bias point (~25-35 min each at 150 trials x 9-10 points)
python run_e3.py --vmtj -0.8
python run_e3.py --vmtj -0.4
python run_e3.py --vmtj 0
python run_e3.py --vmtj 0.4
python run_e3.py --vmtj 0.8      # + low-current supplement, see below

python run_e3.py --dualpulse     # sequential 1 ns pre-pulse variant

python run_e3.py --analyze       # merge -> e3_results.json + two panels
```

The +0.8 V low-current supplement (`e3_v0.8_lo.json`, 250–530 µA) is produced by calling `mc_psw` directly with the same seed base (see git history); `--analyze` merges all files sharing a `v_mtj`.

## Outputs

- `e3_v<V>.json` (+ `_lo` supplements) — P_sw(I_SOT) scans per bias.
- `e3_dualpulse.json` — sequential dual-pulse scan.
- `e3_results.json` — thresholds, sim-vs-analytic scaling, energy table, dual-pulse benefit.
- `e3_ith_vs_vmtj.png`, `e3_psw_family.png` — clean panels (Arial, English, no baked figure numbers; composition happens in PPT).

Deterministic per-trial seeds (`rng_mode`-equivalent generator seeding, base 20260719+3000×); sign convention verified against the chapter (+V assists, −V hinders: P_sw at 1160 µA goes 0.33 → 0.67 → 0.00 for 0/+0.8/−0.8 V).
