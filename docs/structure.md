

## Project structure

The package is organized in three layers: **physics kernels**, **simulation cases/configurations**, and **CLI / result export**. The public API is then re-exported through `__init__.py`. The current code-level structure is:

```python
vgsot_sim/
    __init__.py             # package public API: re-export configs, cases, kernels, IO helpers
    __main__.py             # `python -m vgsot_sim` entry, forwards to `cli.main()`

    constants.py            # physical constants, device geometry, material parameters, time step
    stochastic.py           # normalized Gaussian thermal-noise direction sampler
    initialize.py           # initial MTJ state and initial resistance
    electronic.py           # terminal-voltage -> (I_SOT, V_MTJ) electrical conversion
    anisotropy.py           # effective-field construction: PMA, VCMA, demag, exchange, thermal field
    dynamic_switching.py    # one-step LLG analytical update for theta / phi / mz
    tmr.py                  # MTJ resistance feedback from `V_MTJ` and `mz`

    configs.py              # dataclass configs for all experiment / sweep / SER cases

    time_series_cases.py    # core time-domain simulation kernels and sweep-style cases
    ser_cases.py            # Monte-Carlo SER cases built on top of time-series kernels
    cases.py                # unified case registry and aggregated exports

    result_io.py            # result directory creation, filename stem building, CSV export, plotting
    cli.py                  # command-line entry for selecting cases and saving outputs
```

### Module responsibilities

- **Low-level physics modules**: `constants.py`, `stochastic.py`, `initialize.py`, `electronic.py`, `anisotropy.py`, `dynamic_switching.py`, `tmr.py`
- **Experiment definition layer**: `configs.py`, `time_series_cases.py`, `ser_cases.py`, `cases.py`
- **User interface / output layer**: `cli.py`, `result_io.py`, `__main__.py`, `__init__.py`

### Runtime dependency

- `cli.py` parses the case name, instantiates a config dataclass, runs a case through `cases.py`, and writes CSV / PNG outputs through `result_io.py`.
- `cases.py` combines two families of cases:
  - `time_series_cases.py`: deterministic waveform-driven simulations returning `SimResult` / `SweepResult`
  - `ser_cases.py`: Monte-Carlo SER simulations returning `SerResult` / `SerOptimizedResult`
- `time_series_cases.py` is the main orchestration layer. It:
  - reads constants from `constants.py`
  - initializes the device state via `initialize.py`
  - optionally converts terminal voltages to `I_SOT` and `V_MTJ` via `electronic.py`
  - advances magnetization by calling `dynamic_switching.py`
  - updates resistance through `tmr.py`
- `dynamic_switching.py` calls `anisotropy.py` to build the effective magnetic field.
- `anisotropy.py` calls `stochastic.py` when thermal noise is enabled.
- `ser_cases.py` repeatedly invokes `run_piecewise_direct_excitation()` or `run_two_pulse_optimized()` from `time_series_cases.py` to estimate switching error rate.