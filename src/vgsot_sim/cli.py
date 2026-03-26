from __future__ import annotations

import argparse

from . import cases
from .configs import (
    OptimizedVgsotSwitchingConfig,
    SerOptimizedVgsotConfig,
    SerSotNoVcmaThermalConfig,
    SotOnlyConstantCurrentConfig,
    SotSwitchingNoVcmaConfig,
    TerminalVoltageControlConfig,
    VcmaAssistedSwitchingIsotSweepConfig,
    VcmaAssistedSwitchingVmtjSweepConfig,
)
from .result_io import (
    build_stem,
    ensure_result_dir,
    save_grouped_timeseries_csv,
    save_single_plot,
    save_three_panel_plot,
    save_timeseries_csv,
    save_two_panel_plot,
    save_xy_csv,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="vgsot-sim",
        description="VGSOT MTJ switching simulation demos.",
    )
    parser.add_argument(
        "case",
        choices=list(cases.ALL_CASES),
        help="Which simulation case to run.",
    )
    parser.add_argument(
        "--out_dir",
        default="result",
        help="Output directory for CSV/PNG files.",
    )
    parser.add_argument(
        "--no_progress",
        action="store_true",
        help="Disable progress bars.",
    )
    args = parser.parse_args()

    out_dir = ensure_result_dir(args.out_dir)
    show_progress = not args.no_progress

    if args.case == "terminal_voltage_control":
        cfg = TerminalVoltageControlConfig()
        res = cases.terminal_voltage_control(cfg, show_progress=show_progress)
        stem = build_stem(args.case, cfg)

        save_timeseries_csv(
            out_dir / f"{stem}.csv",
            res.time_s,
            {
                "mz": res.mz,
                "r_mtj": res.r_mtj,
                "v_mtj": res.v_mtj,
                "i_sot": res.i_sot,
                "v1": res.v1,
                "v2": res.v2,
                "v3": res.v3,
                "theta": res.theta,
                "phi": res.phi,
            },
        )
        save_three_panel_plot(
            out_dir / f"{stem}.png",
            res.time_s,
            {"mz": res.mz},
            {r"$R_{\mathrm{MTJ}}$": res.r_mtj},
            {"V1": res.v1, "V2": res.v2, "V3": res.v3},
            ylabel_top="mz",
            ylabel_mid=r"$R_{\mathrm{MTJ}}$",
            ylabel_bot="Terminal Voltage (V)",
            tick_spacing_s=cfg.tick_spacing_s,
            x_is_time=True,
            switch_energy_j=res.switch_energy_j,
        )

    elif args.case == "sot_only_constant_current":
        cfg = SotOnlyConstantCurrentConfig()
        res = cases.sot_only_constant_current(cfg, show_progress=show_progress)
        stem = build_stem(args.case, cfg)

        save_timeseries_csv(
            out_dir / f"{stem}.csv",
            res.time_s,
            {
                "mz": res.mz,
                "r_mtj": res.r_mtj,
                "v_mtj": res.v_mtj,
                "i_sot": res.i_sot,
                "theta": res.theta,
                "phi": res.phi,
            },
        )
        save_three_panel_plot(
            out_dir / f"{stem}.png",
            res.time_s,
            {"mz": res.mz},
            {r"$R_{\mathrm{MTJ}}$": res.r_mtj},
            {r"$I_{\mathrm{SOT}}$": res.i_sot * 1e6},
            ylabel_top="mz",
            ylabel_mid=r"$R_{\mathrm{MTJ}}$",
            ylabel_bot=r"$I_{\mathrm{SOT}}$ ($\mu$A)",
            tick_spacing_s=cfg.tick_spacing_s,
            x_is_time=True,
            switch_energy_j=res.switch_energy_j,
        )

    elif args.case == "sot_switching_no_vcma":
        cfg = SotSwitchingNoVcmaConfig()
        res = cases.sot_switching_no_vcma(cfg, show_progress=show_progress)
        stem = build_stem(args.case, cfg)

        save_grouped_timeseries_csv(
            out_dir / f"{stem}.csv",
            res.time_s,
            {
                "mz": res.mz_curves,
                "r_mtj": res.r_mtj_curves,
                "pulse": res.pulse_curves,
            },
        )
        save_three_panel_plot(
            out_dir / f"{stem}.png",
            res.time_s,
            res.mz_curves,
            res.r_mtj_curves,
            res.pulse_curves,
            ylabel_top="mz",
            ylabel_mid=r"$R_{\mathrm{MTJ}}$",
            ylabel_bot=res.pulse_ylabel,
            tick_spacing_s=cfg.tick_spacing_s,
            legend_title=r"$I_{\mathrm{SOT}}$ sweep",
            x_is_time=True,
            switch_energy_j=res.switch_energy_j,
        )

    elif args.case == "ser_sot_no_vcma_thermal":
        cfg = SerSotNoVcmaThermalConfig()
        res = cases.ser_sot_no_vcma_thermal(cfg, show_progress=show_progress)
        stem = build_stem(args.case, cfg)

        save_xy_csv(
            out_dir / f"{stem}.csv",
            ["i_sot_A", "ser"],
            res.x,
            res.ser,
        )
        save_single_plot(
            out_dir / f"{stem}.png",
            res.x,
            {"SER": res.ser},
            xlabel=res.x_label,
            ylabel="SER",
            x_is_time=False,
        )

    elif args.case == "vcma_assisted_switching_isot_sweep":
        cfg = VcmaAssistedSwitchingIsotSweepConfig()
        res = cases.vcma_assisted_switching_isot_sweep(cfg, show_progress=show_progress)
        stem = build_stem(args.case, cfg)

        save_grouped_timeseries_csv(
            out_dir / f"{stem}.csv",
            res.time_s,
            {
                "mz": res.mz_curves,
                "r_mtj": res.r_mtj_curves,
                "pulse": res.pulse_curves,
            },
        )
        save_three_panel_plot(
            out_dir / f"{stem}.png",
            res.time_s,
            res.mz_curves,
            res.r_mtj_curves,
            res.pulse_curves,
            ylabel_top="mz",
            ylabel_mid=r"$R_{\mathrm{MTJ}}$",
            ylabel_bot=res.pulse_ylabel,
            tick_spacing_s=cfg.tick_spacing_s,
            legend_title=r"$I_{\mathrm{SOT}}$ sweep",
            x_is_time=True,
            switch_energy_j=res.switch_energy_j,
        )

    elif args.case == "vcma_assisted_switching_vmtj_sweep":
        cfg = VcmaAssistedSwitchingVmtjSweepConfig()
        res = cases.vcma_assisted_switching_vmtj_sweep(cfg, show_progress=show_progress)
        stem = build_stem(args.case, cfg)

        save_grouped_timeseries_csv(
            out_dir / f"{stem}.csv",
            res.time_s,
            {
                "mz": res.mz_curves,
                "r_mtj": res.r_mtj_curves,
                "pulse": res.pulse_curves,
            },
        )
        save_three_panel_plot(
            out_dir / f"{stem}.png",
            res.time_s,
            res.mz_curves,
            res.r_mtj_curves,
            res.pulse_curves,
            ylabel_top="mz",
            ylabel_mid=r"$R_{\mathrm{MTJ}}$",
            ylabel_bot=res.pulse_ylabel,
            tick_spacing_s=cfg.tick_spacing_s,
            legend_title=r"$V_{\mathrm{MTJ}}$ sweep",
            x_is_time=True,
            switch_energy_j=res.switch_energy_j,
        )

    elif args.case == "optimized_vgsot_switching":
        cfg = OptimizedVgsotSwitchingConfig()
        res = cases.optimized_vgsot_switching(cfg, show_progress=show_progress)
        stem = build_stem(args.case, cfg)

        save_grouped_timeseries_csv(
            out_dir / f"{stem}.csv",
            res.time_s,
            {
                "mz": res.mz_curves,
                "r_mtj": res.r_mtj_curves,
                "pulse": res.pulse_curves,
            },
        )
        save_three_panel_plot(
            out_dir / f"{stem}.png",
            res.time_s,
            res.mz_curves,
            res.r_mtj_curves,
            res.pulse_curves,
            ylabel_top="mz",
            ylabel_mid=r"$R_{\mathrm{MTJ}}$",
            ylabel_bot=res.pulse_ylabel,
            tick_spacing_s=cfg.tick_spacing_s,
            legend_title="pulse timing sweep",
            x_is_time=True,
            switch_energy_j=res.switch_energy_j,
        )

    elif args.case == "ser_optimized_vgsot":
        cfg = SerOptimizedVgsotConfig()
        res = cases.ser_optimized_vgsot(cfg, show_progress=show_progress)
        stem = build_stem(args.case, cfg)

        save_xy_csv(
            out_dir / f"{stem}.csv",
            ["t1_s", "ser", "mz_at_t1_avg"],
            res.t1_s,
            res.ser,
            res.mz_at_t1_avg,
        )
        save_two_panel_plot(
            out_dir / f"{stem}.png",
            res.t1_s,
            res.ser,
            res.mz_at_t1_avg,
            xlabel1="t1",
            ylabel1="SER",
            xlabel2="t1",
            ylabel2="mz_at_t1_avg",
            x_is_time=True,
        )


if __name__ == "__main__":
    main()