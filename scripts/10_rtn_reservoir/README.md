# RTN reservoir computing — bridge validation (Stage 2) + reservoir benchmarks (Stage 3)

## Stage 2 — sLLG ↔ RTN bridge validation

Does the repo's macrospin **stochastic LLG** engine, run free-running at a low
barrier, reproduce the two-state RTN statistics that
[`vgsot_sim.rtn.telegraph`](../../src/vgsot_sim/rtn/telegraph.py) asserts? This
directory drives that test via [`vgsot_sim.rtn.bridge`](../../src/vgsot_sim/rtn/bridge.py).

Run: `PYTHONPATH=src python scripts/10_rtn_reservoir/validate_bridge.py`
(≈20 min; writes `bridge_results.json`). Fast unit tests of the helpers:
`pytest tests/test_rtn_bridge.py`.

## Setup (and one correction to the earlier audit)

The perpendicular barrier is `Δ = K_u^eff·v/(kB·T)` with
`K_u^eff = Ki/tf − ½·μ0·Ms²·(Nz−Nx)`. The shape-anisotropy term is **large** here
(`Nz−Nx ≈ 0.96`, ~208 kB·T), so a "PMA-only" low-`Ki` (e.g. the `Ki=6.24e-6` used
in the first audit pass) drives `K_u^eff` **negative** → the easy axis goes
in-plane and there is **no** perpendicular telegraph. `bridge.ki_for_delta`
inverts the correct relation: Δ=2 needs `Ki≈2.62e-4` (default `Ki` → Δ=48.5).

The symmetry-breaking bias is a longitudinal field `h_ex_z` (a Zeeman tilt). In
the two-state limit `⟨m_z⟩ = tanh(μ0·Ms·h_ex_z·v/(kB·T))`, so the RTN's abstract
`V` maps linearly to `h_ex_z`. (`V_MTJ` enters the LLG only via *symmetric* VCMA
and cannot tilt the wells; the SOT current drives via a damping torque, not an
energy tilt — so there is no one-to-one `V_MTJ → V` identity.)

Integration uses `t_step = 4 ps`, verified to give the same escape statistics as
1 ps (τ_dwell 88.9 vs 91.1 ns) for a 4× speedup; 10 ps drifts (177 ns) and is too
coarse.

## Results (2026-06-30, `bridge_results.json`)

**(b) Dwell-time distribution → exponential: CONFIRMED.** Δ=2, 8 µs, 68 flips:
τ_dwell = 112.8 ns, **CV = 0.99** (ideal Poisson → 1.0). The CV estimate is noisy
at ~30–50 flips (0.85–1.27) and converges to ~1 with more flips. Equivalent to a
Lorentzian PSD.

**Attempt time τ0 ≈ 15–27 ns — NOT the model's 1 ns.** From `τ0 = τ_dwell/exp(Δ)`:

| Δ | flips | τ_dwell (ns) | CV | τ0 (ns) |
|---|---|---|---|---|
| 1.5 | 40 | 121.8 | 0.89 | 27.2 |
| 2.0 | 31 | 146.4 | 1.05 | 19.8 |
| 2.5 | 26 | 174.8 | 0.85 | 14.4 |
| 2.0 (8 µs) | 68 | 112.8 | 0.99 | 15.3 |

The device's own low-barrier macrospin (α=0.05) hops on a **tens-of-ns** attempt
time, ~15–27× the `τ0 = 1 ns` assumed by the RTN model and `nb_fit`. The reduced
model's *timescales do not transfer* from the LLG without recalibration — the
central calibration finding of Stage 2.

**(a) ⟨m_z⟩(bias) → tanh: shape CONFIRMED, amplitude compressed.** Δ=2:

| h_ex_z (A/m) | tilt | ⟨m_z⟩ | tanh(tilt) |
|---|---|---|---|
| −3000 | −2.08 | −0.69 | −0.97 |
| −1500 | −1.04 | −0.51 | −0.78 |
| 0 | 0 | −0.07 | 0 |
| +1500 | +1.04 | +0.43 | +0.78 |
| +3000 | +2.08 | +0.57 | +0.97 |

⟨m_z⟩ is a monotonic sigmoid through the origin, but saturates at ~0.6–0.7, not
tanh's ~0.97, because the macrospin spends ~30 % of the time in the equatorial
cone (`frac|m_z|>0.5 ≈ 0.71`): even fully polarised, the in-well `m_z ≈ ±0.85`,
not ±1. So the LLG realises `⟨m_z⟩ ≈ A·tanh(g·tilt)` with `A ≈ 0.7` — a
**calibratable gain/amplitude**, not the literal `tanh` identity. This is the
physically-understood version of the first audit's "7× off at V=0.3" (which was an
artifact of its in-plane-collapsed setup, not a real device result).

**(c) Lorentzian PSD.** Established indirectly by the exponential dwell (CV≈1).
The direct PSD-corner fit (`bridge.psd_lorentzian`) is only good to ~2× because the
spectrum hits a white-noise floor a few × f_c above the corner; the authoritative
corner is `τ_c = τ_dwell/2 ≈ 56 ns`.

## Takeaways for the thesis (§2.4 / calibration)

- The reduced RTN tanh/τ structure **does** emerge from the full sLLG at low Δ —
  exponential dwell and a tanh-shaped transfer — so the device-physics grounding is
  real, not just asserted.
- But two honest corrections must be carried into §2.4: (i) `τ0` is **tens of ns**
  for this device, not 1 ns; report timescales in units of the measured `τ0` or
  recalibrate. (ii) the transfer is `A·tanh(g·V)` with `A<1` (thermal-cone
  compression), so a reservoir readout sees a reduced output swing.
- These are good "试错-修正" material: the bridge revealed that the clean reduced
  model needs recalibration before it can be claimed as "the same device".

## Stage 3 — reservoir benchmarks

[`vgsot_sim.rtn.reservoir`](../../src/vgsot_sim/rtn/reservoir.py) builds the layer
the audit found missing: random input weights `W_in` project the scalar input to
each node (`V_j = a_in·w_j·u + b_j`), nodes are heterogeneous in `(Delta_j, Vc0_j)`
(a spread of fading-memory timescales), and a trained ridge readout maps the
node-state design matrix to the target. Run:
`PYTHONPATH=src python scripts/10_rtn_reservoir/benchmark_reservoir.py`
(`reservoir_results.json`). Mean-field state path is vectorised → benchmarks run in
~2 s (they use the analytic telegraph propagator, not the slow LLG).

### Results (2026-06-30, `reservoir_results.json`)

**The W_in reservoir works — and decisively beats the audit's collapse.** Linear
short-term memory capacity (Jaeger; `MC = Σ_k corr²(û[t−k], u[t−k])`):

| reservoir | MC |
|---|---|
| **A2 baseline** — identical nodes, one broadcast input | **0.58** |
| heterogeneous W_in, n=25 | 6.26 |
| heterogeneous W_in, n=50 | 7.56 |
| heterogeneous W_in, n=100 | 7.89 |
| heterogeneous W_in, n=200 | 8.11 |

So per-node input projection + heterogeneity lifts MC ~14× over the broadcast
baseline → **"supports reservoir computing" is now backed by a real number**, not
asserted. `MC_1 ≈ 0.9` (clean immediate recall), decaying with delay (proper fading
memory).

**Filter-bank ceiling (honest limit).** MC saturates at ~8 even at n=200: the nodes
are *independent* leaky integrators (no recurrence), so they form a filter bank
whose linear memory is rank-limited. Pushing MC toward `n` needs inter-node
coupling or the single-node time-delay (Appeltant) architecture — Stage 3.2.

**NARMA-10** (n=200, needs nonlinearity + 10-step memory): NRMSE ≈ 0.55, R² ≈ 0.70
— a working but modest reservoir, consistent with the MC~8 memory depth.

**Device stochasticity is the dominant cost.** With actual binary device states
(`mode="stochastic"`), MC = 0.36 (1 device/node) → 0.67 (4) → 0.94 (16) → 1.23 (64
replicas) — still far below the mean-field 7.9. A hardware sMTJ reservoir therefore
needs heavy device-averaging (or a noise-tolerant readout); single-shot device
states alone are too noisy. This is a concrete hardware-design takeaway.

### Nonlinear capacity, parity, and the delay reservoir

Beyond linear memory, the reservoir has real **nonlinear** computing power. Approximate
information-processing capacity (Dambre 2012, degree-1+2 Legendre basis) and delayed
N-bit parity (a task a linear readout *cannot* solve without nonlinear mixing):

| reservoir (n=100) | linear MC | IPC total | IPC deg-1 | IPC deg-2 | parity-2/3/4 |
|---|---|---|---|---|---|
| filter bank (independent W_in) | 7.6 | 6.60 | 4.97 | 1.62 | 1.00 / 1.00 / 0.91 |
| delay reservoir (Appeltant) | 5.58 | 6.82 | 3.71 | 3.10 | 1.00 / 1.00 / 0.92 |

Both solve parity-2/3/4 well, confirming genuine nonlinear memory (not just linear
recall). The single-node **delay reservoir** (`DelayReservoir`: one RTN node
time-multiplexed into virtual nodes with masked input + delayed feedback) was added
to test whether recurrence beats the filter-bank's linear-MC ceiling. It does **not**:
its linear MC is *lower* (5.6 vs 7.6) and does not scale with `n_virtual`. Instead it
sits at a different operating point — *more* nonlinear capacity (IPC deg-2 3.1 vs 1.6)
for *less* linear memory, comparable total IPC (~6.7). Genuinely raising the linear-MC
ceiling needs a larger/better-tuned coupled network — an honest negative result for the
"recurrence beats the filter bank" hypothesis as implemented.

### Device grounding (energy / latency / window)

Translating the RC operating point to device numbers (`device_grounding.py`,
`device_grounding.json`; `R_MTJ ≈ 5.0 kΩ`, `R_SOT ≈ 776 Ω`, §2.3 write energy
0.78 pJ). The RTN `V` is a phenomenological tilt with an unresolved electrical
mapping (Stage 2), so per-node energy is an order-of-magnitude estimate between two
bounding bias paths, not a calibrated figure.

| operating point | latency / throughput | memory horizon | energy / node-step |
|---|---|---|---|
| **faithful** (τ0≈20 ns, V≈0.2 V) | dt≈20 ns → **~50 MOPS/node** | τ_mem(Δ=2) ≈ 67 ns | ≈ 0.16 pJ (MTJ path) … 1.0 pJ (SOT path) |
| optimistic (model τ0=1 ns) | dt≈1 ns → ~1 GOPS/node | τ_mem ≈ 3 ns | ≈ 8 fJ … 0.05 pJ |

- **Latency is the binding device constraint**: the Stage-2 LLG attempt time
  (~20 ns, not 1 ns) caps a faithful reservoir at ~50 MOPS/node, with a tens-of-ns
  memory horizon — the model's 1 ns / GOPS figure is optimistic by ~20×.
- **Energy** per node-step (~0.04–0.36 pJ on the MTJ path) is at or below a single
  0.78 pJ write — promising, but gated on resolving the bias→V mapping.
- **The real cost driver is device-averaging.** Single-shot device MC = 0.36, so a
  usable node needs ~16–64 physical devices averaged (Stage 3), multiplying area and
  energy by that factor. This dominates any per-step figure.
- **Usable window**: `Δ ≈ 1–3` (τ_max 1.4–10 ns, tanh graded over `|V| ≲ 0.6–1.8 V`),
  `|V| < Vc0`. Higher Δ buys longer memory but a narrower bias window.

### Coupled networks: topology decides (D1)

Adding inter-node coupling to the mean-field reservoir (`ReservoirConfig
coupling_*`, `input_mode`) and scanning:

- **Random sparse Gaussian coupling HURTS.** Every (spectral radius 0.3–1.1 ×
  bias scale 0.3–0.8) tested lowers MC from 7.9 to 5.5–6.7 — recurrence through
  the saturating tanh trades linear memory away, same lesson as the Appeltant
  delay loop.
- **A simple-cycle "ring delay-line" WINS**: homogeneous low-Δ nodes
  (`Δ=1.0`) on a unidirectional ring, input injected into node 0 only
  (`ring_reservoir()` helper; cf. Rodan & Tiño's simple cycle reservoir). With
  the per-hop small-signal gain `g = hop_gain·Δ/Vc0` just below 1 the chain
  propagates the input almost losslessly: **MC = 19/25/33/37 at n=25/50/100/200**
  — ~4× past the filter-bank ceiling, and still growing with n. Instability
  (MC collapse) sets in as `g` crosses ~1, e.g. Δ=1.5 dies between hop products
  0.55 and 0.62 — the tuning rule is physical and sharp.
- **The price is total loss of nonlinearity**: the tuned ring has IPC deg-2 =
  0.00 and parity-2/3/4 at chance (0.53/0.49/0.48) — it maximises linear memory
  precisely by never engaging the tanh. The three constructions populate a
  memory↔nonlinearity frontier (filter bank 7.9/1.6, delay loop 5.6/3.1, ring
  32.6/0.0 as MC/IPC-deg2), mirroring at network level the single-node bias
  trade-off — total capacity is conserved (Dambre), the topology chooses the
  allocation.

### Reservoir-quality metrics (D2) and the device-budget allocation (D3)

**Kernel / generalization rank + ESP** (`kernel_quality`, `esp_convergence`;
n=100, fixed seeds): filter bank — kernel rank 18, generalization rank 14, ESP
distance 2.6e-16 (leaky nodes forget the initial condition to machine
precision); ring delay-line — kernel rank 36 (richer separation, consistent
with its MC), generalization rank 31 (it deliberately remembers the far past),
ESP distance 1.1e-2 after 400 steps (convergence is slow because the memory is
~n hops deep — the flip side of its long memory).

**Budget allocation** (`budget_allocation.py` + json): with a total budget of
`B = n x R` physical devices (n logical nodes, R binary devices averaged per
node, stochastic mode), the optimum at EVERY tested budget is few-nodes ×
deep-averaging — n=8 wins at B=64/256/1024 with best MC = 0.47/0.76/1.21.
Noise reduction (~1/sqrt(R)) beats added dimensionality throughout, and the
best MC grows only ~B^(1/3) over the tested range: approaching the mean-field
limit (7.9) by brute replication would take tens of thousands of devices per
reservoir. Design implication: single-shot binary readout is the wrong regime —
use longer per-step time-averaging (R can equivalently be realised in time) or
readouts robust to binary states.

### Next steps (for thesis write-up, recommended with author review)

- **Thesis write-back** — a §2.4.5 bridge subsection (τ0 ≈ tens of ns; the
  `A·tanh(g·V)` transfer) and an RC chapter, positioned against the spintronic-RC
  literature (Jaeger & Haas 2004, Appeltant 2011, Torrejon 2017, Grollier 2020,
  Welbourne 2021, Dambre 2012; Camsari p-bits are *not* RC — verify DOIs).
- **Optional extensions** — coupled-network reservoir (to actually raise the MC
  ceiling), kernel/generalisation rank, spoken-digit benchmark, ESP test.

