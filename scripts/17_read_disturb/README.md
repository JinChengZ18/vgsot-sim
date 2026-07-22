# 17 — Read-disturb bound for the Bernoulli sampling interface (E8b)

Every hardware sample in §2.3.6 is a write-then-read cycle, and the sampling budget there implicitly assumes the read leaves the sampled bit alone. The chapter never bounded that, so this experiment measures it: the read bias drives both a tunnelling current (spin transfer on the free layer) and a VCMA barrier modulation, and either could in principle flip the state.

Setup: `I_SOT = 0`, read bias held across the MTJ for 10 ns (far longer than a real read window), STT and VCMA both enabled, self-heating on, Cayley kernel, 300 trials per point; a trial counts as disturbed when m_z crosses the equator into the other well.

Findings:

States are named by m_z rather than P/AP on purpose: the chapter's §2.2.2.3 resistance mapping and the kernel's `tmr()` disagree about which sign of m_z is the parallel state (see the backlog item on that contradiction), and the read-disturb result must not depend on which is right.

- **m_z ≈ −1 start**: 0/300 at 10 mV and at 0.2 V; 13/300 at 0.4 V; 299/300 at 0.6 V; 300/300 at 0.8 and 1.0 V.
- **m_z ≈ +1 start**: 0/300 at every bias from 10 mV to 1.0 V. The asymmetry follows the STT polarisation (+z): a positive read bias drives the magnetisation toward +z, so only the m_z ≈ −1 well is at risk *at this polarity*.
- At the chapter's 10 mV read bias the zero count gives a 1.3% 95% upper bound per read; that bound is set by the trial count. The onset is bracketed in (0.2, 0.4] V, so the data-supported drive margin is ~21× in tunnelling current (1.0 µA at the read point vs 21 µA at the highest bias with no observed disturb), rising to ~46× against the first bias where disturb appears.
- **Channel isolation** (0.8 V, m_z ≈ −1): VCMA only (`--channel vcma`, estt = 0) gives 0/120; STT only (`--channel stt`, vnv = 0) gives 119/120, matching the both-channels 300/300. So the disturb is carried by the tunnelling current's spin-transfer torque and VCMA alone cannot flip the device. The control sits at a saturated bias, so it does **not** resolve whether VCMA barrier lowering shifts the disturb threshold — that needs the same comparison at the 0.4 V onset.

## Run

```bash
# Main sweep: one detached process per (bias, state), ~25 min each
python run_e8b.py --bias 0.01 --state p     # ... 0.2, 0.4, 0.6, 0.8, 1.0 x {p, ap}

# Channel-isolation controls at the saturated bias
python run_e8b.py --bias 0.8 --state p --channel vcma --trials 120
python run_e8b.py --bias 0.8 --state p --channel stt  --trials 120

python run_e8b.py --analyze                 # -> e8b_results.json + panel
python run_e8b.py --smoke
```

`--analyze` keeps the channel-isolation runs out of the main sweep (they carry `estt`/`vnv` gates in their JSON) and reports them separately.

## Outputs

- `e8b_<state>_v<bias>.json`, `e8b_p_v0.8_{vcma,stt}.json` — raw counts.
- `e8b_results.json` — per-point probabilities with Wilson 95% upper bounds, the onset bias, and the controls.
- `e8b_read_disturb.png` — disturb probability vs read bias, both states, log axis; clean panel, composition in PPT.
