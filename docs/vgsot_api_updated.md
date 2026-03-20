# API Reference

This document is the detailed Python API guide for `vgsot_sim`.
It is meant to complement the short examples in the root `README.md`.

---

## Quick start (minimal example)

```python
from vgsot_sim import sot_only_constant_current

res = sot_only_constant_current()
print(res.mz[-1])
```

---

## Import surface

```python
from vgsot_sim import (
    # configs
    TerminalVoltageControlConfig,
    SotOnlyConstantCurrentConfig,
    SotSwitchingNoVcmaConfig,
    SerSotNoVcmaThermalConfig,
    VcmaAssistedSwitchingIsotSweepConfig,
    VcmaAssistedSwitchingVmtjSweepConfig,
    OptimizedVgsotSwitchingConfig,
    SerOptimizedVgsotConfig,

    # result types
    SimResult,
    SweepResult,
    SerResult,
    SerOptimizedResult,

    # high-level cases
    terminal_voltage_control,
    sot_only_constant_current,
    sot_switching_no_vcma,
    vcma_assisted_switching_isot_sweep,
    vcma_assisted_switching_vmtj_sweep,
    optimized_vgsot_switching,
    ser_sot_no_vcma_thermal,
    ser_optimized_vgsot,

    # low-level kernels
    run_piecewise_terminal_voltage,
    run_piecewise_direct_excitation,
    run_two_pulse_optimized,
)
```

---

## Mental model

The library has three layers.

1. **Config dataclasses** → define one experiment  
2. **Case functions** → run one experiment  
3. **Kernel functions** → full control over waveform  

👉 Recommendation:
- Most users → use high-level APIs  
- Advanced users → use low-level kernels  

---

## Configuration system

### Default usage

```python
res = ser_sot_no_vcma_thermal()
```

Equivalent to:

```python
cfg = SerSotNoVcmaThermalConfig()
res = ser_sot_no_vcma_thermal(cfg)
```

---

### Inspect defaults

```python
cfg = SerSotNoVcmaThermalConfig()
print(cfg)
```

---

### Convert to dict

```python
from vgsot_sim import config_to_params
print(config_to_params(cfg))
```

---

### Override parameters

```python
cfg = SerSotNoVcmaThermalConfig(trials=200)
```

---

## Shared physical flags (PAP / NON / VNV)

These parameters appear in many APIs and are defined globally:

| Param | Meaning |
|------|--------|
| `pap` | Initial magnetic state: 1 = Parallel, 0 = Anti-parallel |
| `non` | Thermal noise toggle (1 = enable, 0 = disable) |
| `vnv` | VCMA effect toggle (1 = enable, 0 = disable) |

---

## Units

- Time: seconds (s)
- Current: amperes (A)
- Voltage: volts (V)

---

## Result objects

### SimResult

Single-run result container.

Fields:
- time_s
- mz
- r_mtj
- v_mtj
- i_sot
- switch_energy_j

---

## Notes on defaults

Some parameters (e.g., `i_sot=None`) are computed internally from device physics constants.

---

## Low-level API highlight

### run_piecewise_direct_excitation

This is the MOST flexible and recommended low-level API.

Use it when:
- You need custom waveforms
- You want full control over timing

---

