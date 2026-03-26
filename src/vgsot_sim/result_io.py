from __future__ import annotations

import csv
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np


IGNORE_STEM_KEYS = {
    "tick_spacing_s",
    "sim_start_step",
    "estt_stage1",
    "estt_stage2",
    "estt_stage3",
    "esot_stage1",
    "esot_stage2",
    "esot_stage3",
}

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 12
plt.rcParams["axes.linewidth"] = 1.2
plt.rcParams["lines.linewidth"] = 2.0
plt.rcParams["legend.frameon"] = False
plt.rcParams["figure.dpi"] = 200
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["xtick.major.width"] = 1.0
plt.rcParams["ytick.major.width"] = 1.0
plt.rcParams["xtick.major.size"] = 4
plt.rcParams["ytick.major.size"] = 4


def ensure_result_dir(path: str | Path = "result") -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def _sanitize_token(s: str) -> str:
    return str(s).strip().replace(" ", "_").replace("-", "_")


def _fmt_eng(value: float, unit: str) -> str:
    v = float(value)
    u = unit.lower()
    if u == "s":
        if abs(v) >= 1e-3:
            return f"{v:.3g}s"
        if abs(v) >= 1e-6:
            return f"{v*1e3:.3g}ms"
        if abs(v) >= 1e-9:
            return f"{v*1e6:.3g}us"
        if abs(v) >= 1e-12:
            return f"{v*1e9:.3g}ns"
        return f"{v*1e12:.3g}ps"
    if u == "v":
        return f"{v:.4g}V"
    if u == "a":
        if abs(v) >= 1e-3:
            return f"{v*1e3:.4g}mA"
        return f"{v*1e6:.4g}uA"
    return f"{v:.4g}{unit}"


def _summarize_seq(xs: Sequence[float], unit: str) -> str:
    if len(xs) == 0:
        return "empty"
    mn = float(np.min(xs))
    mx = float(np.max(xs))
    n = len(xs)
    if mn == mx:
        return f"{_fmt_eng(mn, unit)}_n{n}"
    return f"{_fmt_eng(mn, unit)}_{_fmt_eng(mx, unit)}_n{n}"


def _format_energy_fj(value_j: float) -> str:
    return f"{float(value_j) * 1e15:.2f} fJ"


def _build_energy_text(
    switch_energy_j: float | Mapping[str, float] | None = None,
    energy_text: str | None = None,
) -> str | None:
    if energy_text is not None:
        return energy_text
    if switch_energy_j is None:
        return None
    if isinstance(switch_energy_j, Mapping):
        lines = [r"$E_{\mathrm{MTJ}}$"]
        for label, value in switch_energy_j.items():
            lines.append(f"{label}: {_format_energy_fj(value)}")
        return "\n".join(lines)
    return rf"$E_{{\mathrm{{MTJ}}}}$ = {_format_energy_fj(switch_energy_j)}"


def _add_energy_box(ax: plt.Axes, text: str | None) -> None:
    if not text:
        return
    ax.text(
        0.02,
        0.98,
        text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox=dict(
            boxstyle="round,pad=0.3",
            facecolor="white",
            edgecolor="gray",
            alpha=0.85,
        ),
    )


def config_to_params(cfg: Any) -> dict[str, Any]:
    if is_dataclass(cfg):
        return asdict(cfg)
    if isinstance(cfg, Mapping):
        return dict(cfg)
    raise TypeError("cfg must be a dataclass instance or mapping.")


def build_stem(case: str, cfg: Any, max_len: int = 120) -> str:
    params = config_to_params(cfg)
    parts = [_sanitize_token(case)]
    for k in sorted(params.keys()):
        if k in IGNORE_STEM_KEYS:
            continue
        v = params[k]
        if v is None:
            continue
        key = _sanitize_token(k)
        if isinstance(v, (list, tuple, np.ndarray)):
            if "isot" in key or key.endswith("_a"):
                val = _summarize_seq(v, "a")
            elif "vmtj" in key or key.endswith("_v"):
                val = _summarize_seq(v, "v")
            elif key.startswith("t") or key.endswith("_s"):
                val = _summarize_seq(v, "s")
            else:
                val = f"n{len(v)}"
        elif isinstance(v, bool):
            val = str(v)
        elif isinstance(v, int):
            val = str(v)
        elif isinstance(v, float):
            if "isot" in key or key.endswith("_a"):
                val = _fmt_eng(v, "a")
            elif "vmtj" in key or key.endswith("_v"):
                val = _fmt_eng(v, "v")
            elif key.startswith("t") or key.endswith("_s"):
                val = _fmt_eng(v, "s")
            else:
                val = f"{v:.4g}"
        else:
            val = str(v)
        token = f"{key}={_sanitize_token(val)}"
        candidate = "__".join(parts + [token])
        if len(candidate) > max_len:
            break
        parts.append(token)
    return "__".join(parts)


def _format_time_xlabel(xlabel: str) -> str:
    lower = xlabel.lower()
    if "ns" in lower:
        return xlabel
    if "time" in lower:
        return "Time (ns)"
    return f"{xlabel} (ns)"


def _as_1d_array(name: str, values: np.ndarray | Sequence[float]) -> np.ndarray:
    arr = np.asarray(values).reshape(-1)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional after flattening.")
    return arr


def _write_csv(path: Path, header: Sequence[str], columns: Sequence[np.ndarray]) -> None:
    if not columns:
        raise ValueError("At least one data column is required.")
    normalized = [_as_1d_array(name, col) for name, col in zip(header, columns)]
    n_rows = len(normalized[0])
    if len(header) != len(normalized):
        raise ValueError("Header length must match the number of data columns.")
    for name, col in zip(header, normalized):
        if len(col) != n_rows:
            raise ValueError(f"Column '{name}' has length {len(col)}, expected {n_rows}.")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for row in zip(*normalized):
            writer.writerow(row)


def save_timeseries_csv(path: str | Path, time_s: np.ndarray, series: Mapping[str, np.ndarray]) -> None:
    path = Path(path)
    header = ["time_s"] + list(series.keys())
    columns = [time_s] + [series[k] for k in series.keys()]
    _write_csv(path, header, columns)


def save_grouped_timeseries_csv(path: str | Path, time_s: np.ndarray, groups: Mapping[str, Mapping[str, np.ndarray]]) -> None:
    path = Path(path)
    header = ["time_s"]
    columns: list[np.ndarray] = [time_s]
    for group_name, group_series in groups.items():
        for label, values in group_series.items():
            header.append(f"{group_name}__{label}")
            columns.append(values)
    _write_csv(path, header, columns)


def save_xy_csv(path: str | Path, header: Sequence[str], *arrays: np.ndarray) -> None:
    path = Path(path)
    _write_csv(path, list(header), list(arrays))


def _finalize_figure(fig: plt.Figure, reserve_right: bool = False) -> None:
    if reserve_right:
        fig.tight_layout(rect=(0.0, 0.0, 0.80, 1.0))
    else:
        fig.tight_layout()


def save_single_plot(
    path: str | Path,
    x: np.ndarray,
    y_dict: Mapping[str, np.ndarray],
    xlabel: str,
    ylabel: str,
    tick_spacing_s: float | None = None,
    legend_title: str | None = None,
    x_is_time: bool = False,
    switch_energy_j: float | Mapping[str, float] | None = None,
    energy_text: str | None = None,
) -> None:
    x_plot = x * 1e9 if x_is_time else x
    xlabel_plot = _format_time_xlabel(xlabel) if x_is_time else xlabel
    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    for label, y in y_dict.items():
        ax.plot(x_plot, y, label=label)
    ax.set_xlabel(xlabel_plot)
    ax.set_ylabel(ylabel)
    if tick_spacing_s is not None:
        spacing = tick_spacing_s * 1e9 if x_is_time else tick_spacing_s
        ax.xaxis.set_major_locator(ticker.MultipleLocator(spacing))
    ax.grid(alpha=0.3, linestyle="--")
    ax.tick_params(direction="in", length=4, width=1.0)
    use_outer_legend = len(y_dict) > 1
    if use_outer_legend:
        ax.legend(title=legend_title, loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=9)
    _add_energy_box(ax, _build_energy_text(switch_energy_j, energy_text))
    _finalize_figure(fig, reserve_right=use_outer_legend)
    fig.savefig(str(path), dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_two_panel_plot(
    path: str | Path,
    x: np.ndarray,
    y1: np.ndarray,
    y2: np.ndarray,
    xlabel1: str,
    ylabel1: str,
    xlabel2: str,
    ylabel2: str,
    tick_spacing_s: float | None = None,
    x_is_time: bool = False,
    switch_energy_j: float | Mapping[str, float] | None = None,
    energy_text: str | None = None,
) -> None:
    x_plot = x * 1e9 if x_is_time else x
    xlabel1_plot = _format_time_xlabel(xlabel1) if x_is_time else xlabel1
    xlabel2_plot = _format_time_xlabel(xlabel2) if x_is_time else xlabel2
    fig, axs = plt.subplots(2, figsize=(7.6, 6.2), sharex=x_is_time)
    axs[0].plot(x_plot, y1)
    axs[0].set(xlabel=xlabel1_plot, ylabel=ylabel1)
    axs[0].grid(alpha=0.3, linestyle="--")
    axs[1].plot(x_plot, y2)
    axs[1].set(xlabel=xlabel2_plot, ylabel=ylabel2)
    axs[1].grid(alpha=0.3, linestyle="--")
    if tick_spacing_s is not None:
        spacing = tick_spacing_s * 1e9 if x_is_time else tick_spacing_s
        axs[0].xaxis.set_major_locator(ticker.MultipleLocator(spacing))
        axs[1].xaxis.set_major_locator(ticker.MultipleLocator(spacing))
    for ax in axs:
        ax.tick_params(direction="in", length=4, width=1.0)
    _add_energy_box(axs[0], _build_energy_text(switch_energy_j, energy_text))
    _finalize_figure(fig)
    fig.savefig(str(path), dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_three_panel_plot(
    path: str | Path,
    x: np.ndarray,
    y_top: Mapping[str, np.ndarray],
    y_mid: Mapping[str, np.ndarray],
    y_bot: Mapping[str, np.ndarray],
    ylabel_top: str,
    ylabel_mid: str,
    ylabel_bot: str,
    tick_spacing_s: float | None = None,
    legend_title: str | None = None,
    x_is_time: bool = True,
    switch_energy_j: float | Mapping[str, float] | None = None,
    energy_text: str | None = None,
) -> None:
    x_plot = x * 1e9 if x_is_time else x
    xlabel_plot = "Time (ns)" if x_is_time else "X"
    fig, axs = plt.subplots(3, 1, figsize=(9.2, 8.2), sharex=True)
    for label, y in y_top.items():
        axs[0].plot(x_plot, y, label=label)
    axs[0].set_ylabel(ylabel_top)
    axs[0].grid(alpha=0.3, linestyle="--")
    for label, y in y_mid.items():
        axs[1].plot(x_plot, y, label=label)
    axs[1].set_ylabel(ylabel_mid)
    axs[1].grid(alpha=0.3, linestyle="--")
    for label, y in y_bot.items():
        axs[2].plot(x_plot, y, label=label)
    axs[2].set_xlabel(xlabel_plot)
    axs[2].set_ylabel(ylabel_bot)
    axs[2].grid(alpha=0.3, linestyle="--")
    if tick_spacing_s is not None:
        spacing = tick_spacing_s * 1e9 if x_is_time else tick_spacing_s
        locator = ticker.MultipleLocator(spacing)
        for ax in axs:
            ax.xaxis.set_major_locator(locator)
    for ax in axs:
        ax.tick_params(direction="in", length=4, width=1.0)
    use_outer_legend = len(y_top) > 1
    if use_outer_legend:
        axs[0].legend(
            title=legend_title,
            fontsize=9,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            borderaxespad=0.0,
        )
    _add_energy_box(axs[0], _build_energy_text(switch_energy_j, energy_text))
    _finalize_figure(fig, reserve_right=use_outer_legend)
    fig.savefig(str(path), dpi=300, bbox_inches="tight")
    plt.close(fig)
