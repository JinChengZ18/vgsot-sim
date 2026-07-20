# 15 — Symmetry-breaking bias-field calibre check (E4)

The experiment measured at H_x = 200 Oe (§2.3.1/2.3.2) while the simulator's default exchange bias is −50 Oe along −ŷ (§2.2.1.2 footnote) — a 4× gap the chapter previously did not disclose. Maintenance docs recorded a ~0.4 SER tail at 200 Oe attributed to unmodeled channels (interfacial DMI, roughness). This experiment maps both knobs at the fig 2.11 workpoint (0.75 ns, P→AP, Cayley, self-heating ON, vnv=0):

1. **Field magnitude** — P_sw(I_SOT) at |H_ex| ∈ {0, 25, 50, 100, 150, 200, 300} Oe along −ŷ (⊥σ̂_SH), 11 points × 150 trials each. Findings: H = 0 keeps a switching threshold but the over-drive P_sw saturates at ≈0.5 instead of 1 (the unbiased-coin limit); 25 Oe reaches the highest peak (0.87, still rising at the grid end) while 50 Oe gives a flat plateau (0.71–0.74); degradation sets in from 100 Oe, the 150–200 Oe over-drive P_sw drops *below* 0.5 (200 Oe band 0.35–0.53) and 300 Oe degenerates to non-monotone fluctuation about 0.5; the nominal first-crossing at 200 Oe sits ~12% above the 50 Oe threshold (no stable above-0.5 plateau exists there). This reproduces (and sharpens) the documented 200 Oe degradation with committed data, and quantifies what the effective θ_SH calibration absorbed.
2. **Alignment tolerance** — rotate the 50 Oe bias from ⊥σ̂_SH toward σ̂_SH (0–45°), probing P_sw at 1160/1400 µA × 200 trials. The collapse is gradual, not abrupt: at 1400 µA, P_sw = 0.72/0.59/0.58/0.46/0.41 for 0/10/20/30/45° — a ±10–15° usable window, random beyond ~30°. The §2.2.1.2 footnote's "any deviation collapses SER to ~0.5" claim was refined accordingly.

## Run

```bash
# One detached process per field magnitude (~25 min each)
python run_e4.py --field 0     # ... 25, 50, 100, 150, 200, 300
python run_e4.py --angles      # angle probe at 50 Oe (~40 min)
python run_e4.py --analyze     # -> e4_results.json + panels
python run_e4.py --smoke
```

## Outputs

- `e4_h<H>.json` — P_sw(I) scans per field magnitude; `e4_angles.json` — angle probe records.
- `e4_results.json` — per-field 50% crossings, top-band statistics, 200/50 Oe threshold ratio, angle records.
- `e4_field_family.png`, `e4_angle_tolerance.png` — clean panels (Arial, English, no baked figure numbers).

Deterministic per-trial seeds (base 20260719+4000×). E5 (reference-layer stray-field asymmetry via `h_ex_z`) is designed to share this directory.
