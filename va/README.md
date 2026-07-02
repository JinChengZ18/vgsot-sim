# `va/` — vgsot-sim Verilog-A engine

This is the **second engine** of vgsot-sim. The Python package (`src/vgsot_sim/`) and this Verilog-A model implement the **same documented macrospin sLLG physics**; they are cross-validated against each other so a result can be reproduced in either a Python loop or an open-source SPICE transient.

## `llg/vgsot_llg.va` — full-dynamics macrospin LLG

A genuine time-domain LLG solver in Verilog-A: it integrates the magnetisation `m(t)` **internally** (Cartesian Landau–Lifshitz form, `ddt` state on internal nodes), not a precomputed compact curve. It is the Verilog-A counterpart of `src/vgsot_sim/{anisotropy,dynamic_switching_vector,tmr,electronic}.py`, per-formula aligned:

- effective field `H_eff = H_PMA + H_VCMA + H_demag(exact oblate ellipsoid) + H_ex + H_th`
- `dm/dt = -g (m×H) - αg (m×(m×H)) + DL/FL spin-orbit torque + DL/FL spin-transfer torque`, `g = γ/(1+α²)`, `σ_SH = -x̂`, `σ_STT = +ẑ`
- electrical: Ohmic SOT write channel `R_sot`; state-dependent MTJ read `R(m_z) = R_P(1+TMR)/(1+TMR(1+m_z)/2)`

### Open-source stochasticity (harness-owned RNG)

OpenVAF/OSDI does not reliably support in-module `$rdist`/persistent random state, so the FDT thermal field is **not drawn inside the module**. Instead the three input nodes `hx, hy, hz` carry the thermal field `H_th` [A/m]; the Python harness drives them with the correctly scaled Brown-1963 white-noise stream (same seed → identical noise in both engines). Grounding `hx/hy/hz` gives the deterministic limit. This is the macrospin generalisation of the `eda/models/smtj_sot.va` harness-RNG design.

### Build & run (open-source: OpenVAF-Reloaded + ngspice ≥ 43)

```sh
cd va/llg
openvaf vgsot_llg.va -o vgsot_llg.osdi          # compile to OSDI
ngspice_con -b tb_switch.spice                   # .spiceinit loads the OSDI model
```

`tb_switch.spice` forces a super-threshold `I_SOT` (via a voltage across the SOT channel), holds `V_MTJ = 0`, starts `m` near `+z`, and writes `m(t)` to `vgsot_llg_out.csv`. State nodes `mx, my, mz` are terminals so the testbench sets `.ic` and probes `v(mx)` etc.; the transient uses `uic`.

### Status

- **Deterministic dynamics: validated.** `|m|=1` conserved to ~1e-4 with no manual renormalisation; m_z(t) matches the Python `switching_vector` trajectory to **max ~0.006 over 0–3 ns** (exact agreement at the current-on equilibrium; the small transient difference is fixed-step Cayley vs ngspice adaptive `ddt`).
- **Stochastic via harness:** the `h_th` injection path is in place; the seeded Python↔VA Monte-Carlo equivalence harness and its committed regression are pending (see project `#15`).
- **`theta_SH` default** is the recalibrated `0.066` (Cayley, post-FL-SOT-fix; `0.04` was the pre-fix value). The `.va` exposes `theta_SH` as a parameter, so the calibrated value can be set without editing the model.

The compact behavioural Verilog-A (operating-point sigmoid + telegraph observables) lives separately in the consuming project at `eda/models/smtj_sot.va`; it is a fast circuit surrogate, not a replacement for this first-principles LLG engine.
