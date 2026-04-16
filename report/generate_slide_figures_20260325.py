#!/usr/bin/env python3
"""Generate the unified 2026-03-25 slide figures."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter, MaxNLocator

ROOT = Path(__file__).resolve().parent
DATADIR = ROOT / "data"
FIGDIR = ROOT / "figures"
FIGDIR.mkdir(exist_ok=True)

QMC_J2_0_L8 = -0.673487
QMC_J2_0_L16 = -0.669976
YUBIN_GO_J2_0_D5 = -0.6733
YUBIN_GO_J2_0_L16 = -0.6696

MEASURE = {
    (0.0, 5, "CU"): (-0.671811, 0.000033),
    (0.0, 5, "SR"): (-0.673108, 0.000023),
    (0.5, 5, "CU"): (-0.494904, 0.000056),
    (0.5, 5, "SR"): (-0.497786, 0.000047),
    (0.0, 8, "CU"): (-0.671875, 0.000256),
    (0.5, 8, "CU"): (-0.496481, 0.000278),
    (0.0, 8, "SR"): (-0.673316, 0.000020),
}

TRACE_ENDPOINTS = {
    "j2_0_stepselector": -0.673489,
    "j2_0_minsr": -0.673402,
    "j2_05_sr_cu": -0.49814,
    "j2_05_sr_su": -0.49752,
    "j2_05_stepselector_560212": -0.49850,
    "j2_05_minsr": -0.49855,
    "yubin_go_d8": -0.4978,
    "yubin_go_d9": -0.4980,
}

COLORS = {
    "d5": "#2878B5",
    "d8": "#D65F5F",
    "stepsel": "#2C6E49",
    "minsr": "#C98C10",
    "su": "#8A63D2",
    "qmc": "#222222",
    "yubin_d8": "#7C6FD6",
    "yubin_d9": "#C05A8A",
    "yubin_l16": "#6E59A5",
    "muted": "#6C757D",
}

LR_LABELS = {
    "sr_start": "start lr = 0.1",
    "stepsel_j2_0": "adaptive lr: 0.3 -> 0.15 -> 0.075",
    "stepsel_j2_05": "green trace = original StepSel + job 560212",
    "minsr": "lr = 0.1",
}


def set_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 220,
            "savefig.dpi": 220,
            "figure.facecolor": "white",
            "axes.facecolor": "#FCFDFE",
            "font.size": 11,
            "axes.titlesize": 14,
            "axes.labelsize": 12,
            "axes.edgecolor": "#C7CDD4",
            "axes.linewidth": 1.0,
            "axes.titleweight": "bold",
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "xtick.color": "#38424C",
            "ytick.color": "#38424C",
            "text.color": "#18212A",
            "legend.frameon": False,
            "legend.fontsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def load_csv_energy(path: Path, denom: float = 64.0) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    iters, energies, errors, grads = [], [], [], []
    with path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            iters.append(int(row["iteration"]))
            energies.append(float(row["energy"]) / denom)
            errors.append(float(row["energy_error"]) / denom)
            grads.append(float(row["gradient_norm"]))
    return np.array(iters), np.array(energies), np.array(errors), np.array(grads)


def load_simple_trace(path: Path) -> tuple[np.ndarray, np.ndarray]:
    iters, energies = [], []
    with path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            iters.append(int(row["iteration"]))
            energies.append(float(row["energy_per_site"]))
    return np.array(iters), np.array(energies)


def load_jsonl_trace(path: Path, denom: float, *, offset: int = 0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    iters, energies, errors = [], [], []
    with path.open() as f:
        for line in f:
            row = json.loads(line)
            iters.append(int(row["iter"]) + offset)
            energies.append(float(row["energy"]) / denom)
            errors.append(float(row["energy_error"]) / denom)
    return np.array(iters), np.array(energies), np.array(errors)


def load_jsonl_tail(path: Path, denom: float, nrows: int, *, offset: int = 0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = []
    with path.open() as f:
        for line in f:
            row = json.loads(line)
            rows.append(
                (
                    int(row["iter"]) + offset,
                    float(row["energy"]) / denom,
                    float(row["energy_error"]) / denom,
                )
            )
    rows = rows[-nrows:]
    iters = np.array([row[0] for row in rows])
    energies = np.array([row[1] for row in rows])
    errors = np.array([row[2] for row in rows])
    return iters, energies, errors


def align_digitized_endpoint(energies: np.ndarray, exact_endpoint: float) -> tuple[np.ndarray, float]:
    shift = exact_endpoint - float(energies[-1])
    return energies + shift, shift


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        0.015,
        0.985,
        label,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        color="white",
        va="top",
        bbox={"boxstyle": "round,pad=0.28", "facecolor": "#233445", "edgecolor": "none"},
    )


def soften_grid(ax: plt.Axes) -> None:
    ax.grid(axis="y", color="#DDE3EA", linewidth=0.8, alpha=0.9)
    ax.set_axisbelow(True)


def format_energy_axis(ax: plt.Axes, decimals: int = 4, nbins: int = 5) -> None:
    ax.xaxis.set_major_locator(MaxNLocator(nbins=nbins))
    ax.xaxis.set_major_formatter(FormatStrFormatter(f"%.{decimals}f"))


def add_legend(ax: plt.Axes, **kwargs) -> None:
    legend = ax.legend(loc="upper right", frameon=True, fancybox=True, framealpha=0.94, **kwargs)
    legend.get_frame().set_facecolor("white")
    legend.get_frame().set_edgecolor("#DDE3EA")
    legend.get_frame().set_linewidth(0.8)


def plot_with_band(
    ax: plt.Axes,
    iters: np.ndarray,
    energies: np.ndarray,
    errors: np.ndarray | None,
    *,
    color: str,
    label: str,
    lw: float = 2.0,
    ls: str = "-",
    marker: str | None = None,
    alpha_band: float = 0.12,
) -> None:
    if errors is not None:
        ax.fill_between(iters, energies - errors, energies + errors, color=color, alpha=alpha_band, linewidth=0)
    ax.plot(iters, energies, color=color, lw=lw, ls=ls, marker=marker, markersize=4 if marker else 0, label=label)


def add_footer(fig: plt.Figure, text: str) -> None:
    fig.text(
        0.5,
        0.05,
        text,
        ha="center",
        fontsize=10,
        color=COLORS["muted"],
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "#F3F6F9", "edgecolor": "#DDE3EA"},
    )


def fig_j2_0_overview() -> None:
    d5_it, d5_e, d5_err, _ = load_csv_energy(DATADIR / "8x8_J2=0_D5_vmc_trajectory.csv")
    d8_it, d8_e, d8_err, _ = load_csv_energy(DATADIR / "8x8_J2=0_D8_vmc_trajectory.csv")
    ss_it, ss_e, ss_err, _ = load_csv_energy(DATADIR / "8x8_J2=0_D8_stepselector_vmc.csv")
    ms_it, ms_e, ms_err, _ = load_csv_energy(DATADIR / "8x8_J2=0_D8_minsr_lr01.csv")

    fig = plt.figure(figsize=(12.8, 5.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.65, 1.0], left=0.07, right=0.98, top=0.90, bottom=0.16, wspace=0.22)
    ax = fig.add_subplot(gs[0, 0])
    axr = fig.add_subplot(gs[0, 1])

    plot_with_band(ax, d5_it, d5_e, d5_err, color=COLORS["d5"], lw=2.0, label=f"D=5 SR ({LR_LABELS['sr_start']})")
    plot_with_band(ax, d8_it, d8_e, d8_err, color=COLORS["d8"], lw=2.1, label=f"D=8 SR ({LR_LABELS['sr_start']})")
    plot_with_band(
        ax,
        ss_it,
        ss_e,
        ss_err,
        color=COLORS["stepsel"],
        lw=2.3,
        label=f"D=8 StepSelector ({LR_LABELS['stepsel_j2_0']})",
    )
    plot_with_band(ax, ms_it, ms_e, ms_err, color=COLORS["minsr"], lw=2.0, label=f"D=8 MinSR ({LR_LABELS['minsr']})")
    ax.axhline(QMC_J2_0_L8, color=COLORS["qmc"], lw=1.6, ls="--", label="QMC")
    ax.axhline(YUBIN_GO_J2_0_D5, color=COLORS["yubin_d8"], lw=1.2, ls="-.", alpha=0.85, label="Yu-Bin GO D=5")
    ax.axhline(MEASURE[(0.0, 5, "CU")][0], color=COLORS["d5"], lw=1.2, ls=":", alpha=0.8)
    ax.axhline(MEASURE[(0.0, 8, "CU")][0], color=COLORS["d8"], lw=1.2, ls=":", alpha=0.8)
    ax.text(
        0.015,
        0.03,
        "shaded band = +/- 1sigma MC error; dotted = CU baselines",
        transform=ax.transAxes,
        fontsize=8.5,
        color=COLORS["muted"],
    )
    ax.set_xlim(-1, 60)
    ax.set_ylim(-0.67378, -0.67155)
    ax.set_xlabel("Optimization step")
    ax.set_ylabel("Energy per site")
    ax.set_title("Optimization traces", loc="left", pad=10)
    soften_grid(ax)
    add_panel_label(ax, "A")
    add_legend(ax, fontsize=8.2)

    axr.axvline(QMC_J2_0_L8, color=COLORS["qmc"], lw=1.5, ls="--")
    axr.axvline(YUBIN_GO_J2_0_D5, color=COLORS["yubin_d8"], lw=1.2, ls="-.")
    rows = [1, 0]
    labels = ["D=8", "D=5"]
    baselines = [MEASURE[(0.0, 8, "CU")][0], MEASURE[(0.0, 5, "CU")][0]]
    finals = [MEASURE[(0.0, 8, "SR")][0], MEASURE[(0.0, 5, "SR")][0]]
    row_colors = [COLORS["d8"], COLORS["d5"]]
    for y, label, base, final, color in zip(rows, labels, baselines, finals, row_colors):
        axr.plot([base, final], [y, y], color=color, lw=3.2, solid_capstyle="round")
        axr.scatter([base], [y], color="white", edgecolor=color, s=70, zorder=3, linewidth=2)
        axr.scatter([final], [y], color=color, s=70, zorder=3)
        axr.text(final, y + 0.12, f"{final:.6f}", fontsize=9, color=color, ha="center")
    axr.scatter([TRACE_ENDPOINTS["j2_0_minsr"], TRACE_ENDPOINTS["j2_0_stepselector"]], [1.22, 1.42], color=[COLORS["minsr"], COLORS["stepsel"]], s=60, zorder=4)
    axr.text(TRACE_ENDPOINTS["j2_0_minsr"], 1.30, "MinSR", color=COLORS["minsr"], fontsize=9, ha="center")
    axr.text(TRACE_ENDPOINTS["j2_0_stepselector"], 1.50, "StepSel", color=COLORS["stepsel"], fontsize=9, ha="center")
    axr.text(YUBIN_GO_J2_0_D5, 1.62, "Yu-Bin GO D=5", color=COLORS["yubin_d8"], fontsize=8.5, ha="center")
    axr.text(QMC_J2_0_L8, -0.38, "QMC", color=COLORS["qmc"], fontsize=9, ha="center")
    axr.set_yticks(rows)
    axr.set_yticklabels(labels)
    axr.set_ylim(-0.1, 1.72)
    axr.set_xlim(-0.67378, -0.67155)
    axr.set_xlabel("Energy per site")
    axr.set_title("Final energies", loc="left", pad=10)
    format_energy_axis(axr, decimals=4, nbins=5)
    soften_grid(axr)
    add_panel_label(axr, "B")
    add_footer(fig, "Raising D alone barely moves the CU baseline. VMC is what closes the gap to QMC.")
    fig.savefig(FIGDIR / "slide_fig_j2_0_overview.png")
    plt.close(fig)


def fig_j2_05_initial_state() -> None:
    d5_it, d5_e, d5_err, _ = load_csv_energy(DATADIR / "8x8_J2=0.5_D5_vmc_trajectory.csv")
    cu_it, cu_e, cu_err, _ = load_csv_energy(DATADIR / "8x8_J2=0.5_D8_yubin_vmc_partial.csv")
    su_it, su_e, su_err, _ = load_csv_energy(DATADIR / "8x8_J2=0.5_D8_su_vmc_partial.csv")

    fig = plt.figure(figsize=(12.8, 5.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.65, 1.0], left=0.07, right=0.98, top=0.90, bottom=0.16, wspace=0.22)
    ax = fig.add_subplot(gs[0, 0])
    axr = fig.add_subplot(gs[0, 1])

    plot_with_band(ax, d5_it, d5_e, d5_err, color=COLORS["d5"], lw=2.1, label=f"D=5 SR ({LR_LABELS['sr_start']})")
    plot_with_band(ax, cu_it, cu_e, cu_err, color=COLORS["d8"], lw=2.2, label=f"D=8 from Yu-Bin CU ({LR_LABELS['sr_start']})")
    plot_with_band(ax, su_it, su_e, su_err, color=COLORS["su"], lw=2.0, label=f"D=8 from SU ({LR_LABELS['sr_start']})")
    ax.axhline(TRACE_ENDPOINTS["yubin_go_d8"], color=COLORS["yubin_d8"], lw=1.2, ls="-.", alpha=0.85, label="Yu-Bin GO D=8")
    ax.axhline(MEASURE[(0.5, 5, "CU")][0], color=COLORS["d5"], lw=1.1, ls=":", alpha=0.8)
    ax.axhline(MEASURE[(0.5, 8, "CU")][0], color=COLORS["d8"], lw=1.1, ls=":", alpha=0.8)
    ax.text(
        0.015,
        0.03,
        "shaded band = +/- 1sigma MC error; dotted = measured CU baselines",
        transform=ax.transAxes,
        fontsize=8.5,
        color=COLORS["muted"],
    )
    ax.set_xlim(-1, 30)
    ax.set_ylim(-0.49895, -0.4944)
    ax.set_xlabel("Optimization step")
    ax.set_ylabel("Energy per site")
    ax.set_title("How the runs move", loc="left", pad=10)
    soften_grid(ax)
    add_panel_label(ax, "A")
    add_legend(ax, fontsize=8.3)

    summary_labels = ["D=5 CU", "D=5 SR", "D=8 CU", "D=8 SR from CU", "D=8 SR from SU"]
    summary_vals = [
        MEASURE[(0.5, 5, "CU")][0],
        MEASURE[(0.5, 5, "SR")][0],
        MEASURE[(0.5, 8, "CU")][0],
        TRACE_ENDPOINTS["j2_05_sr_cu"],
        TRACE_ENDPOINTS["j2_05_sr_su"],
    ]
    summary_colors = [COLORS["d5"], COLORS["d5"], COLORS["d8"], COLORS["d8"], COLORS["su"]]
    ypos = np.arange(len(summary_labels))[::-1]
    axr.axvline(TRACE_ENDPOINTS["yubin_go_d8"], color=COLORS["yubin_d8"], lw=1.2, ls="-.")
    for y, label, val, color in zip(ypos, summary_labels, summary_vals, summary_colors):
        axr.scatter(val, y, color=color, s=62, zorder=3)
        axr.text(val + 0.00005, y, f"{val:.5f}", va="center", fontsize=10, color=color)
    axr.plot([MEASURE[(0.5, 5, "CU")][0], MEASURE[(0.5, 5, "SR")][0]], [ypos[0], ypos[1]], color=COLORS["d5"], lw=2.5, alpha=0.7)
    axr.plot([MEASURE[(0.5, 8, "CU")][0], TRACE_ENDPOINTS["j2_05_sr_cu"]], [ypos[2], ypos[3]], color=COLORS["d8"], lw=2.5, alpha=0.7)
    axr.set_yticks(ypos)
    axr.set_yticklabels(summary_labels)
    axr.set_xlim(-0.49895, -0.4944)
    axr.set_xlabel("Energy per site")
    axr.set_title("Endpoint comparison", loc="left", pad=10)
    format_energy_axis(axr, decimals=4, nbins=5)
    soften_grid(axr)
    add_panel_label(axr, "B")
    axr.text(TRACE_ENDPOINTS["yubin_go_d8"], 4.25, "Yu-Bin GO D=8", color=COLORS["yubin_d8"], fontsize=8.5, ha="center")
    add_footer(fig, "CU starts better than SU, and our D = 8 SR trajectory already moves below the Yu-Bin GO reference.")
    fig.savefig(FIGDIR / "slide_fig_j2_05_initial_state.png")
    plt.close(fig)


def fig_j2_05_optimizer_merged() -> None:
    y8_it, y8_e_raw = load_simple_trace(DATADIR / "yubin_j05_D8_trajectory.csv")
    y9_it, y9_e_raw = load_simple_trace(DATADIR / "yubin_j05_D9_trajectory.csv")
    y8_e, d8_shift = align_digitized_endpoint(y8_e_raw, TRACE_ENDPOINTS["yubin_go_d8"])
    y9_e, d9_shift = align_digitized_endpoint(y9_e_raw, TRACE_ENDPOINTS["yubin_go_d9"])

    sr_it, sr_e, sr_err, _ = load_csv_energy(DATADIR / "8x8_J2=0.5_D8_yubin_vmc_partial.csv")
    ss_old_it, ss_old_e, ss_old_err, _ = load_csv_energy(DATADIR / "8x8_J2=0.5_D8_stepselector_vmc.csv")
    old_mask = ss_old_it < 30
    ss_old_it = ss_old_it[old_mask]
    ss_old_e = ss_old_e[old_mask]
    ss_old_err = ss_old_err[old_mask]
    ss_new_it, ss_new_e, ss_new_err = load_jsonl_tail(
        ROOT / "../run/8x8J2=0.5D8/20260307_yubin_vmc_stepselector/vmc/energy/optimization_log.jsonl",
        64.0,
        30,
        offset=30,
    )
    ms_it, ms_e, ms_err = load_jsonl_trace(
        ROOT / "../run/8x8J2=0.5D8/20260319_icpx_minsr_60iter/vmc/energy/optimization_log.jsonl",
        64.0,
    )

    sr_best_idx = int(np.argmin(sr_e))
    ss_combined_it = np.concatenate([ss_old_it, ss_new_it])
    ss_combined_e = np.concatenate([ss_old_e, ss_new_e])
    ss_combined_err = np.concatenate([ss_old_err, ss_new_err])
    ss_best_idx = int(np.argmin(ss_combined_e))
    ms_best_idx = int(np.argmin(ms_e))

    fig = plt.figure(figsize=(12.8, 5.25))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.65, 1.0], left=0.07, right=0.98, top=0.90, bottom=0.16, wspace=0.22)
    ax = fig.add_subplot(gs[0, 0])
    axr = fig.add_subplot(gs[0, 1])

    ax.plot(y8_it, y8_e, color=COLORS["yubin_d8"], lw=2.0, ls="-.", label="Yu-Bin GO D=8")
    ax.plot(y9_it, y9_e, color=COLORS["yubin_d9"], lw=1.8, ls="--", label="Yu-Bin GO D=9")
    plot_with_band(ax, sr_it, sr_e, sr_err, color=COLORS["d8"], lw=2.1, label=f"Our SR ({LR_LABELS['sr_start']})")
    plot_with_band(
        ax,
        ss_old_it,
        ss_old_e,
        ss_old_err,
        color=COLORS["stepsel"],
        lw=2.2,
        marker="D",
        label="Our StepSel (green trace; extended by job 560212)",
    )
    plot_with_band(
        ax,
        ss_new_it,
        ss_new_e,
        ss_new_err,
        color=COLORS["stepsel"],
        lw=2.2,
        marker="D",
        label="_nolegend_",
    )
    plot_with_band(ax, ms_it, ms_e, ms_err, color=COLORS["minsr"], lw=2.0, label=f"Our MinSR ({LR_LABELS['minsr']})")
    ax.axvline(29.5, color="#BFC6CF", lw=1.2, ls="--", alpha=0.9)
    ax.text(
        31.0,
        -0.49618,
        "job 560212 starts\n12824 samples, lr = 0.1",
        fontsize=8.3,
        color=COLORS["muted"],
        va="top",
        bbox={"boxstyle": "round,pad=0.2", "facecolor": "white", "edgecolor": "none", "alpha": 0.8},
    )
    ax.axhline(TRACE_ENDPOINTS["yubin_go_d8"], color=COLORS["yubin_d8"], lw=1.0, ls=":", alpha=0.65)
    ax.axhline(TRACE_ENDPOINTS["j2_05_sr_cu"], color=COLORS["d8"], lw=1.0, ls=":", alpha=0.65)
    ax.text(
        0.015,
        0.03,
        f"shaded band = +/- 1sigma MC error; digitized GO curves aligned by {d8_shift:+.1e} (D=8) and {d9_shift:+.1e} (D=9)",
        transform=ax.transAxes,
        fontsize=8.2,
        color=COLORS["muted"],
    )
    ax.set_xlim(-1, 60)
    ax.set_ylim(-0.49905, -0.4953)
    ax.set_xlabel("Optimization step")
    ax.set_ylabel("Energy per site")
    ax.set_title("Canonical D = 8 benchmark plot", loc="left", pad=10)
    soften_grid(ax)
    add_panel_label(ax, "A")
    add_legend(ax, fontsize=7.8)

    rank_labels = ["GO D=8", "GO D=9", "baseline SR", "StepSel 560212", "MinSR"]
    rank_vals = [
        TRACE_ENDPOINTS["yubin_go_d8"],
        TRACE_ENDPOINTS["yubin_go_d9"],
        float(sr_e[sr_best_idx]),
        float(ss_combined_e[ss_best_idx]),
        float(ms_e[ms_best_idx]),
    ]
    rank_errs = [0.0, 0.0, float(sr_err[sr_best_idx]), float(ss_combined_err[ss_best_idx]), float(ms_err[ms_best_idx])]
    rank_colors = [COLORS["yubin_d8"], COLORS["yubin_d9"], COLORS["d8"], COLORS["stepsel"], COLORS["minsr"]]
    ypos = np.arange(len(rank_labels))[::-1]
    for y, label, val, err, color in zip(ypos, rank_labels, rank_vals, rank_errs, rank_colors):
        axr.hlines(y, xmin=-0.49875, xmax=val, color="#E6EAF0", lw=2)
        if err > 0:
            axr.errorbar(val, y, xerr=err, fmt="o", color=color, ecolor=color, elinewidth=1.4, capsize=3, ms=6.5)
        else:
            axr.scatter(val, y, color=color, s=66, zorder=3)
        axr.text(val + 0.000045, y, f"{val:.5f}", va="center", fontsize=9.6, color=color)
    axr.text(float(ss_combined_e[ss_best_idx]), 1.18, "canonical green trace", fontsize=8.8, color=COLORS["stepsel"], ha="center")
    axr.text(-0.49872, 4.35, "560212 replaces the old dark-blue continuation", fontsize=8.8, color=COLORS["stepsel"], ha="left")
    axr.set_yticks(ypos)
    axr.set_yticklabels(rank_labels)
    axr.set_xlim(-0.49875, -0.49745)
    axr.set_ylim(-0.2, 4.6)
    axr.set_xlabel("Best energy reached")
    axr.set_title("Best values used in the merged deck", loc="left", pad=10)
    format_energy_axis(axr, decimals=4, nbins=5)
    soften_grid(axr)
    add_panel_label(axr, "B")

    add_footer(
        fig,
        "This is the canonical J2 = 0.5 D = 8 comparison: the old blue low-sample continuation is removed, and job 560212 extends the green StepSelector trace.",
    )
    fig.savefig(FIGDIR / "slide_fig_j2_05_optimizer_merged.png")
    plt.close(fig)


def fig_16x16_j2_0_benchmark() -> None:
    yb_it, yb_e_raw = load_simple_trace(DATADIR / "yubin_heis_grad_L16_auto.csv")
    yb_e, yb_shift = align_digitized_endpoint(yb_e_raw, YUBIN_GO_J2_0_L16)

    sr_it, sr_e, sr_err = load_jsonl_trace(ROOT / "../run/16x16J2=0D8/20260319_icpx_sr_fresh/vmc/energy/optimization_log.jsonl", 256.0)
    ms_it, ms_e, ms_err = load_jsonl_trace(ROOT / "../run/16x16J2=0D8/20260319_icpx_minsr_fresh/vmc/energy/optimization_log.jsonl", 256.0)
    sr_best_idx = int(np.argmin(sr_e))
    ms_best_idx = int(np.argmin(ms_e))

    fig = plt.figure(figsize=(12.8, 5.3))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.65, 1.0], left=0.07, right=0.98, top=0.90, bottom=0.16, wspace=0.22)
    ax = fig.add_subplot(gs[0, 0])
    axr = fig.add_subplot(gs[0, 1])

    ax.plot(yb_it, yb_e, color=COLORS["yubin_l16"], lw=2.0, ls="-.", label="Yu-Bin GO D=5")
    plot_with_band(ax, sr_it, sr_e, sr_err, color=COLORS["d8"], label="Our SR, lr = 0.1", lw=2.2)
    plot_with_band(ax, ms_it, ms_e, ms_err, color=COLORS["minsr"], label="Our MinSR, lr = 0.1", lw=2.1)
    ax.axhline(QMC_J2_0_L16, color=COLORS["qmc"], lw=1.6, ls="--", label="QMC")
    ax.scatter(sr_it[sr_best_idx], sr_e[sr_best_idx], color=COLORS["d8"], s=34, zorder=5)
    ax.scatter(ms_it[ms_best_idx], ms_e[ms_best_idx], color=COLORS["minsr"], s=34, zorder=5)
    ax.text(
        0.015,
        0.03,
        f"shaded band = +/- 1sigma MC error; Yu-Bin curve shifted by {yb_shift:+.2e} to match the exact draft endpoint",
        transform=ax.transAxes,
        fontsize=8.4,
        color=COLORS["muted"],
    )
    ax.set_xlim(-1, 60)
    ax.set_ylim(-0.67003, -0.66902)
    ax.set_xlabel("Optimization step")
    ax.set_ylabel("Energy per site")
    ax.set_title("Trajectory benchmark", loc="left", pad=10)
    soften_grid(ax)
    add_panel_label(ax, "A")
    add_legend(ax, fontsize=8.0)

    labels = ["Yu-Bin GO D=5", "Our SR D=8", "Our MinSR D=8", "QMC"]
    vals = [YUBIN_GO_J2_0_L16, float(sr_e[sr_best_idx]), float(ms_e[ms_best_idx]), QMC_J2_0_L16]
    errs = [0.0, float(sr_err[sr_best_idx]), float(ms_err[ms_best_idx]), 0.0]
    row_colors = [COLORS["yubin_l16"], COLORS["d8"], COLORS["minsr"], COLORS["qmc"]]
    ypos = np.arange(len(labels))[::-1]
    for y, label, val, err, color in zip(ypos, labels, vals, errs, row_colors):
        axr.hlines(y, xmin=-0.67003, xmax=val, color="#E6EAF0", lw=2.0)
        if err > 0:
            axr.errorbar(val, y, xerr=err, fmt="o", color=color, ecolor=color, elinewidth=1.5, capsize=3, ms=6.5)
        else:
            marker = "D" if label == "QMC" else "o"
            axr.scatter(val, y, color=color, s=66, marker=marker, zorder=3)
        axr.text(val + 0.000015, y, f"{val:.6f}", va="center", fontsize=9.4, color=color)
    axr.text(float(sr_e[sr_best_idx]), 1.25, "iter 45", fontsize=8.6, color=COLORS["d8"], ha="center")
    axr.text(float(ms_e[ms_best_idx]), 0.25, "iter 34", fontsize=8.6, color=COLORS["minsr"], ha="center")
    axr.set_yticks(ypos)
    axr.set_yticklabels(labels)
    axr.set_xlim(-0.67003, -0.66952)
    axr.set_xlabel("Best energy reached")
    axr.set_title("Endpoint benchmark", loc="left", pad=10)
    format_energy_axis(axr, decimals=4, nbins=5)
    soften_grid(axr)
    add_panel_label(axr, "B")

    add_footer(fig, "16x16 is now a completed benchmark: SR and MinSR both beat Yu-Bin D=5 and land 4.6e-4 above QMC.")
    fig.savefig(FIGDIR / "slide_fig_16x16_j2_0_benchmark.png")
    plt.close(fig)


def main() -> None:
    set_style()
    fig_j2_0_overview()
    fig_j2_05_initial_state()
    fig_j2_05_optimizer_merged()
    fig_16x16_j2_0_benchmark()
    print("Saved unified 0325 slide figures to", FIGDIR)


if __name__ == "__main__":
    main()
