

---

## 📊 Simulation Cases and Their Physical Meaning

Below is a full description of each available case and what physical mechanism it is designed to demonstrate.

------

## 🔹 `terminal_voltage_control`

**Scenario:** Three-terminal voltage control (V1, V2, V3)
**Physics focus:** Coupled electrical control of **V_MTJ** and **I_SOT** through the electronic module.

### What happens

- Uses terminal voltages (V1, V2, V3) to compute (I_SOT, V_MTJ)
- Feeds (I_SOT, V_MTJ) into the switching dynamics
- Produces mz(t), R_MTJ(t), V_MTJ(t)

### Physical meaning

A **device-level control demo** showing how external terminal voltages map into:

- SOT current (torque-driven switching)
- MTJ bias (anisotropy modulation via VCMA, if enabled)

------

## 🔹 `sot_only_constant_current`

**Scenario:** Constant SOT current pulse only (V_MTJ = 0)
**Physics focus:** Pure **spin–orbit torque switching** baseline.

### Physical meaning

A **reference case** for:

- baseline switching trajectory without VCMA assistance
- comparisons against VCMA-assisted schemes

------

## 🔹 `sot_switching_no_vcma`

**Scenario:** SOT switching without VCMA (VNV=0), sweep I_SOT
**Physics focus:** switching driven by SOT only (no anisotropy modulation).

### Physical meaning

A **control experiment** answering:

> How much SOT current is needed for switching if VCMA is disabled?

------

## 🔹 `ser_sot_no_vcma_thermal`

**Scenario:** Monte-Carlo SER vs I_SOT without VCMA, with thermal noise (NON=1)
**Physics focus:** thermal stochasticity and reliability.

### Physical meaning

Quantifies **switching error rate (SER)** under temperature-driven noise:

- switching becomes probabilistic
- useful for device-level reliability evaluation

⚠️ Computationally heavy (Monte-Carlo trials).

------

## 🔹 `vcma_assisted_switching_isot_sweep`

**Scenario:** VCMA enabled (VNV=1), fix V_MTJ and sweep I_SOT
**Physics focus:** interplay of VCMA barrier modulation and SOT strength.

### Physical meaning

Demonstrates that **VCMA reduces required I_SOT** (lower energy / faster switching).

------

## 🔹 `vcma_assisted_switching_vmtj_sweep`

**Scenario:** VCMA enabled (VNV=1), fix I_SOT and sweep V_MTJ
**Physics focus:** voltage control of anisotropy barrier via VCMA.

### Physical meaning

Shows how **stronger V_MTJ (VCMA)** changes switching trajectory and speed by modulating the anisotropy energy landscape.

------

## 🔹 `optimized_vgsot_switching`

**Scenario:** Optimized two-pulse VCMA+SOT scheme, sweep (t1,t2)
**Physics focus:** pulse engineering for energy-efficient switching.

### Physical meaning

A performance showcase where:

- VCMA temporarily lowers the barrier
- SOT drives reversal more efficiently
- pulse timing (t1,t2) trades off speed vs margin

------

## 🔹 `ser_optimized_vgsot`

**Scenario:** Monte-Carlo SER vs t1 for optimized scheme (t2 = total_pulse - t1)
**Physics focus:** reliability of the optimized scheme under thermal noise.

### Physical meaning

Answers:

> Does the optimized VCMA+SOT pulse remain reliable with thermal fluctuations?

Outputs SER(t1) and the average mz at the end of the first pulse.

⚠️ Computationally heavy (Monte-Carlo trials).

