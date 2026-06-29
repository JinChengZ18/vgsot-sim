# Parameter validation against fabrication literature

The vgsot-sim default parameters (`configs.PhysicalConstantsConfig`) are cross-checked here against two Hikstor β-W/CoFeB SOT-MTJ fabrication papers, to confirm they sit in physically defensible ranges and to flag the deliberately *effective* values.

> **Scope (important).** These papers are used ONLY as **material- and process-node reference** for the β-W/CoFeB stack (spin-Hall angle, M_s, anisotropy, RA, channel geometry *ranges*); they describe **different devices** and are **not** a validation target. The model's calibration target and reported results — V_th, R_P, R_AP, and the P_sw transition — are aligned to the **same-batch Device-A experimental measurements** of Chapter §2.3 (V_th(0.75 ns)=894 mV, R_P=4.9 kΩ, R_AP=10 kΩ, Sigmoid slope β_s=44.6 V⁻¹), not to the papers. Where the modelled device differs from a paper's device (e.g. RA 16.6 vs 36 Ω·µm²), the model follows the **Device-A experiment**.

Sources:
- **[P1]** E. Liu, W. Yang, …, S. He, *"A novel Channel-less SOT-MRAM with 115% TMR, 2 ns Switching, and High Bit Yield (>99.9%)"*, IEDM 2024 (Hikstor). Fig. 2(a) micromagnetic parameters + measured statistics.
- **[P2]** W. Yang et al., *"Achieving High Yield of Perpendicular SOT-MTJ Manufactured on 300 mm Wafers"*, IEEE TED 2024 (Hikstor).

| Quantity | vgsot-sim default | Literature (P1 / P2) | Assessment |
|---|---|---|---|
| Spin Hall angle θ_SH | **effective**, 0.04 → **0.066** (recalibrated, Cayley) | **0.3** (β-W, both P1 & P2: P1 Fig.2a; P2 "θ_SH = −0.3") | The model uses a small *effective* θ_SH because the single-macrospin + lumped-thermal model omits SOT-channel dissipation (Néel–Edelstein, spin-memory loss, parasitic series R, current-direction spread). The FL-SOT integrator fix raises the required effective θ_SH from 0.04 to 0.066 (V_th(0.75 ns)=0.895 V), moving it **toward** the literature 0.3 though still ~4.5× below it — i.e. the lumped model still absorbs a large spin-torque-efficiency deficit. |
| Saturation magnetisation Ms | 0.625 MA/m | 1.0–1.1 MA/m (P1: 1100 kA/m) | On the low side for CoFeB; acceptable for a thin pinned FL but flagged. A higher Ms would raise demag/H_k. |
| Anisotropy | K_i = 0.32 mJ/m² (interfacial); K_U,eff ≈ K_i/t_f − ½μ₀Ms²(N_z−N_x) | K_u = 0.845 MJ/m³ (P1) | Comparable order (K_i/t_f ≈ 0.29 MJ/m³ before demag); consistent. |
| Gilbert damping α | 0.05 | 0.1 (P1 micromagnetic); 0.005–0.01 intrinsic CoFeB | Between intrinsic and the micromagnetic-effective value; reasonable. |
| TMR | 100% (TMR=1.0, R_AP/R_P=2) | 115% (P1) / 119–130% (P2) | Slightly conservative; matches the chapter Device-A hysteresis amplitude R_AP/R_P≈2. |
| Resistance-area RA | 16.6 Ω·µm² (→ R_P≈5 kΩ) | **36 Ω·µm²** (P2) → R_P≈10.85 kΩ | Both physical: 16.6 matches the chapter Device-A low-R state (~4.9 kΩ); 36 Ω·µm² is the 300 mm-wafer device (R_P≈10.85 kΩ). The genuine Simmons/BDR predictor `resistance_area_bdr` reproduces RA from the barrier at MgO m*≈0.3 m_e. |
| Parallel resistance R_P | 5.0 kΩ | 5–9 kΩ (P1 Fig.8/19/21) | ✓ in range. |
| SOT channel R_SOT | 776 Ω | 300–800 Ω (P1 Fig.11; CHL ~500–600) | ✓ in range. |
| MTJ diameter | D_phys=80 nm, D_elec=65 nm | CD 80–110 nm (P1/P2 ~80–100) | ✓; D_elec accounts for ~7.5 nm edge-etch damage band. |
| SOT channel thickness | d = 4.3 nm | ~5 nm (P2); β-W 4–5 nm window (P1 Fig.3) | ✓ in range. |
| Critical current I_c | ~1.15 mA @ 0.75 ns | 660 µA @ 5 ns (P1), 680 µA @ 2 ns (P2) | Higher, consistent with the shorter 0.75 ns pulse (V_c rises monotonically as pulse shortens, P1 Fig.13) and the small effective θ_SH. |
| Néel–Brown switching Δ | 4.91 (smtj_pbnn) / 5.15 (vgsot NB fit) | retention Δ = 50–55 (P1, switching-field method) | **Different quantities**: the model's Δ is the Néel-Brown *switching-law* exponent governing P_sw(V,t_w) near threshold; the literature 50–55 is the *retention* barrier E_b/k_BT. They are not directly comparable; both are reported for completeness. |

**Conclusion.** All electrical/geometric parameters (R_P, R_SOT, RA, CD, channel thickness) sit squarely within the Hikstor β-W/CoFeB process ranges, confirming the inputs are physically reasonable. The notable *effective* parameter is θ_SH: the model deliberately uses a port-level effective value (well below the literature 0.3) to absorb dissipation channels the lumped model does not resolve; Ms is on the low side. These are documented modelling choices, not errors. To restate the scope: this literature comparison establishes only that the **material/process inputs** are credible — the model's **outputs** (R_P=5 kΩ↔Device-A 4.9 kΩ, V_th=0.895 V↔894 mV, R_AP/R_P≈2↔hysteresis amplitude) are validated against the §2.3 **Device-A measurements**, which remain the sole quantitative alignment target.
