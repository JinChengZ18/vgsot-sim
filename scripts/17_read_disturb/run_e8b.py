"""E8b — read-disturb bound for the Bernoulli sampling interface.

Every hardware sample in section 2.3.6 is a write-then-read cycle, and the
sampling budget there implicitly assumes the read leaves the state alone.
The chapter never bounds that: the read applies a bias across the MTJ,
which drives both a tunnelling current (STT on the free layer) and a VCMA
barrier modulation. This experiment measures the disturb probability of a
read directly.

Both mechanisms are enabled (estt = 1, vnv = 1) and I_SOT = 0, so the only
drive is the read bias itself. Both initial states are covered but only the
positive bias polarity is swept: the STT polarisation is +z, so a positive
bias pushes the magnetisation toward +z and only the m_z ~ -1 well is at
risk. That is the worst case for this polarity; reversing the read polarity
moves the vulnerability to the other well by symmetry. The measurement point
is the chapter's 10 mV read bias; the sweep continues far above it to locate
the onset and therefore the margin.

The channel-isolation controls (--channel) run at a saturated bias, so they
establish that VCMA alone cannot flip the device while STT alone reproduces
the full result. They do NOT resolve whether VCMA barrier lowering shifts the
disturb threshold — that would need the comparison at the onset bias.

Modes:
  python run_e8b.py --bias 0.4 --state p     # -> e8b_p_v0.4.json
  python run_e8b.py --analyze                # -> e8b_results.json + panel
  python run_e8b.py --smoke
"""
from __future__ import annotations

import argparse
import glob
import json
import time
from datetime import datetime
from pathlib import Path

import numpy as np

from vgsot_sim.configs import PhysicalConstantsConfig
from vgsot_sim.time_series_cases import run_piecewise_direct_excitation

HERE = Path(__file__).resolve().parent
BASE_SEED = 20260719

READ_NS = 10.0          # read window; the chapter's read is far shorter
BIASES_V = (0.01, 0.2, 0.4, 0.6, 0.8, 1.0)
TRIALS = 300


def _trial_seed(seed, v, trial):
    return (int(abs(seed)) * 0x9E3779B1 + int(round(v * 1e6)) * 0x85EBCA6B
            + int(trial) * 0xC2B2AE35) & 0x7FFFFFFF


def run_point(v_mtj: float, state: str, trials: int = TRIALS, read_ns: float = READ_NS,
              estt: int = 1, vnv: int = 1):
    """Hold `v_mtj` across the MTJ for `read_ns` with no SOT drive and count
    trials whose final m_z has left the starting well.

    `estt` / `vnv` gate the two channels a read bias can act through (spin
    transfer from the tunnelling current, and VCMA barrier modulation), so
    turning one off isolates the other."""
    cc = PhysicalConstantsConfig()
    steps = int(round(read_ns * 1e-9 / cc.t_step))
    pap = 1 if state == "p" else 0          # pap=1 starts at m_z ~ -1
    start = -1.0 if pap == 1 else +1.0
    disturbed = 0
    t0 = time.time()
    for tr in range(trials):
        rng = np.random.default_rng(_trial_seed(BASE_SEED + 8000 + 100 * estt + 10 * vnv,
                                               v_mtj, tr))
        res = run_piecewise_direct_excitation(
            sim_start_step=1, sim_mid1_step=steps, sim_mid2_step=steps,
            sim_end_step=steps, pap=pap,
            v_mtj_stage1=v_mtj, v_mtj_stage2=0.0, v_mtj_stage3=0.0,
            i_sot_stage1=0.0, i_sot_stage2=0.0, i_sot_stage3=0.0,
            estt_stage1=estt, esot_stage1=0, estt_stage2=estt, esot_stage2=0,
            estt_stage3=estt, esot_stage3=0,
            vnv=vnv, non=1, r_sot_fl_dl=0.83, show_progress=False,
            constants=cc, enable_self_heating=True, T_ambient_K=300.0,
            integrator="cayley", rng=rng,
        )
        if float(res.mz[steps]) * start < 0:      # crossed to the other well
            disturbed += 1
    return dict(v_mtj=v_mtj, state=state, read_ns=read_ns, trials=trials,
                estt=estt, vnv=vnv,
                disturbed=disturbed, p_disturb=disturbed / trials,
                wall_s=time.time() - t0,
                timestamp=datetime.now().isoformat(timespec="seconds"))


def wilson_upper(k, n, z=1.96):
    """Upper end of the Wilson interval - the useful quantity when k = 0."""
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return float(centre + half)


def analyze():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    allrecs = []
    for f in sorted(glob.glob(str(HERE / "e8b_*_v*.json"))):
        allrecs.append(json.loads(Path(f).read_text(encoding="utf-8")))
    if not allrecs:
        raise SystemExit("No e8b shard JSONs found.")
    # Channel-isolation controls carry estt/vnv gates; keep them out of the
    # main sweep and report them separately.
    recs = [r for r in allrecs if r.get("estt", 1) == 1 and r.get("vnv", 1) == 1]
    controls = [r for r in allrecs if r not in recs]
    recs.sort(key=lambda d: (d["state"], d["v_mtj"]))
    controls.sort(key=lambda d: (d["v_mtj"], d.get("estt", 1)))

    for r in recs:
        r["p_upper95"] = wilson_upper(r["disturbed"], r["trials"])
    at_read = [r for r in recs if abs(r["v_mtj"] - 0.01) < 1e-9]
    worst_read = max((r["p_upper95"] for r in at_read), default=float("nan"))

    for c in controls:
        c["p_upper95"] = wilson_upper(c["disturbed"], c["trials"])
    summary = dict(
        read_ns=READ_NS, trials=TRIALS, records=recs, channel_controls=controls,
        read_bias_V=0.01, read_disturb_upper95=worst_read,
        onset_V=next((r["v_mtj"] for r in recs if r["disturbed"] > 0), None),
        timestamp=datetime.now().isoformat(timespec="seconds"),
    )
    (HERE / "e8b_results.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    CHARCOAL, NAVY, CRIMSON = "#2B2B2B", "#1F5FA8", "#A82038"
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans"],
        "font.size": 12, "axes.labelsize": 13, "axes.titlesize": 13.5,
        "mathtext.fontset": "stix", "axes.linewidth": 0.9,
        "axes.edgecolor": CHARCOAL, "axes.spines.top": False, "axes.spines.right": False,
        "xtick.direction": "in", "ytick.direction": "in",
        "figure.dpi": 150, "savefig.dpi": 300,
    })
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    for state, colour, lbl in (("p", NAVY, "start P"), ("ap", CRIMSON, "start AP")):
        rs = [r for r in recs if r["state"] == state]
        if not rs:
            continue
        v = [r["v_mtj"] for r in rs]
        p = [max(r["p_disturb"], 1e-4) for r in rs]
        up = [r["p_upper95"] for r in rs]
        ax.plot(v, p, "o-", ms=6, lw=1.4, color=colour, mfc="white", mew=1.5, label=lbl)
        ax.plot(v, up, ":", lw=1.2, color=colour, alpha=0.8)
    ax.set_yscale("log")
    ax.set_xlabel(r"read bias $|V_{\mathrm{MTJ}}|$ (V)")
    ax.set_ylabel("disturb probability per read")
    ax.set_title(f"Read disturb over a {READ_NS:g} ns hold, no SOT drive")
    ax.legend(fontsize=9.5)
    fig.tight_layout()
    fig.savefig(HERE / "e8b_read_disturb.png", bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print("=" * 68)
    print(f"E8b read disturb ({READ_NS:g} ns hold, I_SOT = 0, STT + VCMA on)")
    print("=" * 68)
    print(f"{'state':>6} {'V_MTJ (V)':>10} {'disturbed':>10} {'p':>8} {'95% upper':>10}")
    for r in recs:
        print(f"{r['state']:>6} {r['v_mtj']:>10.2f} {r['disturbed']:>4}/{r['trials']:<5} "
              f"{r['p_disturb']:>8.4f} {r['p_upper95']:>10.4f}")
    if controls:
        print("\n  channel isolation (same bias, one channel disabled):")
        for c in controls:
            chan = "VCMA only (estt=0)" if c.get("estt", 1) == 0 else "STT only (vnv=0)"
            print(f"    {c['state']:>3} @ {c['v_mtj']:.2f} V, {chan:>18}: "
                  f"{c['disturbed']:>4}/{c['trials']:<5} (p = {c['p_disturb']:.4f})")
    print(f"\nAt the 10 mV read bias: disturb probability < {worst_read:.4f} "
          f"per read (95% upper bound, {TRIALS} trials)")
    print(f"First bias with any observed disturb: {summary['onset_V']}")
    print("Wrote e8b_results.json, e8b_read_disturb.png")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bias", type=float, default=None)
    ap.add_argument("--state", choices=("p", "ap"), default="p")
    ap.add_argument("--trials", type=int, default=TRIALS)
    ap.add_argument("--channel", choices=("both", "vcma", "stt"), default="both",
                    help="which read-bias channel is active: both (default), "
                         "vcma (estt=0) or stt (vnv=0)")
    ap.add_argument("--analyze", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()

    if args.smoke:
        for st in ("p", "ap"):
            out = run_point(1.0, st, trials=6, read_ns=2.0)
            print(f"[smoke] {st} @1.0 V, 2 ns: disturbed {out['disturbed']}/6 "
                  f"({out['wall_s']:.0f} s)")
        return
    if args.bias is not None:
        estt, vnv = {"both": (1, 1), "vcma": (0, 1), "stt": (1, 0)}[args.channel]
        print(f"[E8b] {args.state} @ {args.bias} V x {args.trials} trials "
              f"(channel={args.channel})", flush=True)
        out = run_point(args.bias, args.state, trials=args.trials, estt=estt, vnv=vnv)
        suffix = "" if args.channel == "both" else f"_{args.channel}"
        p = HERE / f"e8b_{args.state}_v{args.bias:g}{suffix}.json"
        p.write_text(json.dumps(out, indent=2), encoding="utf-8")
        print(f"[E8b] wrote {p.name} (disturbed {out['disturbed']}/{out['trials']}, "
              f"wall {out['wall_s']/60:.1f} min)", flush=True)
        return
    if args.analyze:
        analyze()
        return
    ap.error("choose one of --bias/--state, --analyze, or --smoke")


if __name__ == "__main__":
    main()
