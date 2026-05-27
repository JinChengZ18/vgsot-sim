from __future__ import annotations

import argparse

from . import cases
from .configs import (
    SerSotNoVcmaThermalConfig,
    SotOnlyConstantCurrentConfig,
    SotSwitchingNoVcmaConfig,
    TerminalVoltageControlConfig,
)
from .result_io import (
    build_stem,
    ensure_result_dir,
    save_grouped_timeseries_csv,
    save_single_plot,
    save_three_panel_plot,
    save_timeseries_csv,
    save_xy_csv,
)

DEFAULT_FIGURE_FILENAMES = {
    "terminal_voltage_control": "Chapter02_local_24.png",
    "sot_only_constant_current": "Chapter02_local_25.png",
    "sot_switching_no_vcma": "Chapter02_local_26.png",
    "ser_sot_no_vcma_thermal": "Chapter02_local_27.png",
}


def _default_figure_path(out_dir, case: str):
    return out_dir / DEFAULT_FIGURE_FILENAMES[case]


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="vgsot-sim",
        description="VGSOT-MTJ switching simulation cases (chapter §2.3.3 protocol).",
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
            _default_figure_path(out_dir, args.case),
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
            _default_figure_path(out_dir, args.case),
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
            _default_figure_path(out_dir, args.case),
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
            ["i_sot_A", "ser", "psw"],
            res.x,
            res.ser,
            res.psw,
        )
        save_single_plot(
            _default_figure_path(out_dir, args.case),
            -res.x * 1e6,                         # plot vs |I_SOT| in µA, chapter convention
            {r"$P_{\mathrm{sw}}$": res.psw},
            xlabel=r"$|I_{\mathrm{SOT}}|$ ($\mu$A)",
            ylabel=r"$P_{\mathrm{sw}}$",
            x_is_time=False,
        )

    elif args.case == "variability_sweep":
        # Specialised analysis case — not part of the CLI surface because it
        # requires NB-fit inputs (Δ, V_c0, β_meas) that aren't on the same
        # footing as the device-physics cases above. Call it from Python:
        #
        #     from vgsot_sim.ser_cases import variability_sweep
        #     res = variability_sweep(Delta=5.15, Vc0=0.884, ...)
        #
        # The case name is still in `cases.ALL_CASES` so the registry stays
        # truthful, but exposing it through argparse would require a much
        # richer CLI than the device-physics cases warrant.
        raise SystemExit(
            "variability_sweep is a Python-only API — see "
            "`docs/cases.md` / `docs/api.md` for usage."
        )


if __name__ == "__main__":
    main()
