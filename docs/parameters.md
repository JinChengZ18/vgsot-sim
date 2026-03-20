

## Default parameters by case

> These are the **function default values** in `src/vgsot_sim/cli.py`. They determine what runs when you simply execute `vgsot-sim <case_name>`.

### `terminal_voltage_control`

| Parameter               | Default                         | Meaning                                    |
| ----------------------- | ------------------------------- | ------------------------------------------ |
| `out_dir`               | `"result"`                      | Output folder for PNG/CSV                  |
| `sim_start_step`        | `1`                             | Simulation start step index                |
| `sim_mid1_step`         | `2000`                          | Stage boundary (1 → 2) in steps            |
| `sim_mid2_step`         | `None` → equals `sim_mid1_step` | Stage boundary (2 → 3) in steps            |
| `sim_end_step`          | `5000`                          | Simulation end step index                  |
| `pap`                   | `1`                             | Initial P/AP state selector                |
| `v_stage1`              | `(1.0, 0.0, 0.1)`               | (V1,V2,V3) in stage 1                      |
| `v_stage2`              | `(-1.0, 0.0, 0.0)`              | (V1,V2,V3) in stage 2 (default duration 0) |
| `v_stage3`              | `(0.0, 0.0, 0.0)`               | (V1,V2,V3) in stage 3                      |
| `vnv`                   | `1`                             | VCMA enable flag passed to switching       |
| `non`                   | `0`                             | Thermal noise flag passed to switching     |
| `r_sot_fl_dl`           | `0.83`                          | Field-like / damping-like ratio            |
| `tick_spacing_s`        | `5e-10`                         | X-axis major tick spacing (s)              |
| `save_csv` / `save_png` | `True` / `True`                 | Save CSV and PNG                           |

------

### `sot_only_constant_current`

| Parameter               | Default         | Meaning                                  |
| ----------------------- | --------------- | ---------------------------------------- |
| `sim_start_step`        | `1`             | Start step                               |
| `sim_mid1_step`         | `2000`          | Pulse end step                           |
| `sim_end_step`          | `5000`          | End step                                 |
| `pap`                   | `1`             | Initial P/AP                             |
| `i_sot_stage1`          | `-95e-6` A      | SOT current during pulse                 |
| `i_sot_stage2`          | `0.0` A         | Relax stage SOT current                  |
| `v_mtj_stage1`          | `0.0` V         | MTJ bias during pulse                    |
| `v_mtj_stage2`          | `0.0` V         | MTJ bias during relax                    |
| `vnv` / `non`           | `1` / `0`       | VCMA on, noise off (passed to switching) |
| `r_sot_fl_dl`           | `0.83`          | FL/DL ratio                              |
| `tick_spacing_s`        | `5e-10`         | Tick spacing                             |
| `save_csv` / `save_png` | `True` / `True` | Save CSV and PNG                         |

------

### `sot_switching_no_vcma`

| Parameter               | Default                               | Meaning                                  |
| ----------------------- | ------------------------------------- | ---------------------------------------- |
| `i_sot_list`            | `(-85e-6, -90e-6, -95e-6, -100e-6)` A | Sweep list for I_SOT                     |
| `sim_mid1_step`         | `2000`                                | Pulse end step                           |
| `sim_end_step`          | `5000`                                | End step                                 |
| `v_mtj`                 | `0.0` V                               | MTJ bias (fixed)                         |
| `i_sot_relax`           | `0.0` A                               | Relax stage current                      |
| `vnv`                   | `0`                                   | VCMA disabled                            |
| `non`                   | `0`                                   | Noise disabled (set `non=1` for thermal) |
| `r_sot_fl_dl`           | `0.0`                                 | FL/DL ratio (kept 0 here)                |
| `tick_spacing_s`        | `5e-10`                               | Tick spacing                             |
| `save_csv` / `save_png` | `True` / `True`                       | Save CSV and PNG                         |

------

### `ser_sot_no_vcma_thermal`

| Parameter               | Default                                               | Meaning                                              |
| ----------------------- | ----------------------------------------------------- | ---------------------------------------------------- |
| `i_sot_list`            | `(-100e-6, -98e-6, -96e-6, -94e-6, -92e-6, -90e-6)` A | Sweep list for SER curve                             |
| `trials`                | `1000`                                                | Monte-Carlo trials per I_SOT                         |
| `sim_mid1_step`         | `2000`                                                | Pulse end step                                       |
| `sim_end_step`          | `5000`                                                | End step                                             |
| `v_mtj`                 | `0.0` V                                               | MTJ bias (fixed)                                     |
| `vnv`                   | `0`                                                   | VCMA disabled                                        |
| `non`                   | `1`                                                   | Thermal noise enabled                                |
| `target_mz`             | `-1.0`                                                | Successful final target state                        |
| `failure_tol`           | `1e-1`                                                | Failure threshold: `abs(mz_final - target_mz) > tol` |
| `save_csv` / `save_png` | `True` / `True`                                       | Save CSV and PNG                                     |

------

### `vcma_assisted_switching_isot_sweep`

| Parameter               | Default                              | Meaning                       |
| ----------------------- | ------------------------------------ | ----------------------------- |
| `v_mtj`                 | `1.2` V                              | Fixed MTJ bias (VCMA control) |
| `i_sot_list`            | `(-90e-6, -30e-6, -18e-6, -16e-6)` A | Sweep list for I_SOT          |
| `sim_end_step`          | `25000`                              | End step (longer window)      |
| `vnv` / `non`           | `1` / `0`                            | VCMA on, noise off            |
| `r_sot_fl_dl`           | `0.0`                                | FL/DL ratio                   |
| `tick_spacing_s`        | `5e-9`                               | Tick spacing                  |
| `save_csv` / `save_png` | `True` / `True`                      | Save CSV and PNG              |

------

### `vcma_assisted_switching_vmtj_sweep`

| Parameter               | Default                                     | Meaning                                       |
| ----------------------- | ------------------------------------------- | --------------------------------------------- |
| `i_sot`                 | `None` → computed by formula                | Fixed I_SOT; if None, computed from constants |
| `v_mtj_list`            | `(1.3189, 1.3191, 1.333, 1.3489, 1.4937)` V | Sweep list for V_MTJ                          |
| `sim_end_step`          | `25000`                                     | End step                                      |
| `vnv` / `non`           | `1` / `0`                                   | VCMA on, noise off                            |
| `r_sot_fl_dl`           | `0.0`                                       | FL/DL ratio                                   |
| `tick_spacing_s`        | `5e-9`                                      | Tick spacing                                  |
| `save_csv` / `save_png` | `True` / `True`                             | Save CSV and PNG                              |

------

### `optimized_vgsot_switching`

| Parameter               | Default                                                      | Meaning                     |
| ----------------------- | ------------------------------------------------------------ | --------------------------- |
| `v_mtj_1` / `v_mtj_2`   | `1.4937` V / `-1.0` V                                        | Two-pulse MTJ voltages      |
| `i_sot`                 | `None` → computed by formula                                 | Pulse SOT current           |
| `t_pairs_s`             | `((25e-9,0), (1.4e-9,1.6e-9), (1.8e-9,1.2e-9), (2.2e-9,0.8e-9))` | Sweep (t1,t2) pairs         |
| `sim_total_time_s`      | `25e-9` s                                                    | Total simulated time window |
| `vnv` / `non`           | `1` / `0`                                                    | VCMA on, noise off          |
| `tick_spacing_s`        | `5e-9`                                                       | Tick spacing                |
| `save_csv` / `save_png` | `True` / `True`                                              | Save CSV and PNG            |

------

### `ser_optimized_vgsot`

| Parameter               | Default                                                    | Meaning                              |
| ----------------------- | ---------------------------------------------------------- | ------------------------------------ |
| `iterations_num`        | `300`                                                      | Monte-Carlo trials per t1            |
| `v_mtj_1` / `v_mtj_2`   | `1.4937` V / `-1.0` V                                      | Two-pulse MTJ voltages               |
| `i_sot`                 | `None` → computed by formula                               | Pulse SOT current                    |
| `t1_list_s`             | `(1.3e-9, 1.4e-9, 1.5e-9, 1.6e-9, 1.7e-9, 1.8e-9, 1.9e-9)` | Sweep list for t1                    |
| `total_pulse_s`         | `3e-9` s                                                   | Constraint `t2 = total_pulse_s - t1` |
| `sim_total_time_s`      | `25e-9` s                                                  | Total simulated time window          |
| `vnv` / `non`           | `1` / `1`                                                  | VCMA on, thermal noise on            |
| `target_final_mz`       | `1.0`                                                      | Successful final target state        |
| `failure_tol`           | `1e-1`                                                     | Failure threshold                    |
| `save_csv` / `save_png` | `True` / `True`                                            | Save CSV and PNG                     |