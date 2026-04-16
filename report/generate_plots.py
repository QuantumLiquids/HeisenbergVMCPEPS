#!/usr/bin/env python3
"""Generate plots for PEPS+VMC research report (updated 2026-03-09)."""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import csv
import json
import os

matplotlib.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'legend.fontsize': 10,
    'figure.dpi': 150,
    'savefig.dpi': 150,
})

DATADIR = os.path.join(os.path.dirname(__file__), 'data')
FIGDIR = os.path.join(os.path.dirname(__file__), 'figures')
os.makedirs(FIGDIR, exist_ok=True)

# QMC benchmarks (J2=0, PBC)
QMC_J2_0 = {8: -0.673487, 12: -0.670685, 16: -0.669976}

# Measure results {(J2, D, source): (E/site, err/site)}
MEASURE = {
    (0.0, 5, 'CU'):       (-0.671811, 0.000033),
    (0.0, 5, 'CU+VMC'):   (-0.673108, 0.000023),
    (0.5, 5, 'CU'):       (-0.494904, 0.000056),
    (0.5, 5, 'CU+VMC'):   (-0.497786, 0.000047),
    (0.0, 8, 'CU'):       (-0.671875, 0.000256),
    (0.5, 8, 'CU'):       (-0.496481, 0.000278),
    (0.0, 8, 'CU+VMC'):   (-0.673316, 0.000020),
}


def load_csv(fname):
    path = os.path.join(DATADIR, fname)
    iters, energies, errors, grads = [], [], [], []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            iters.append(int(row['iteration']))
            energies.append(float(row['energy']))
            errors.append(float(row['energy_error']))
            grads.append(float(row['gradient_norm']))
    return np.array(iters), np.array(energies), np.array(errors), np.array(grads)


def load_jsonl_trace(path, denom, tail=None, offset=0):
    rows = []
    with open(path) as f:
        for line in f:
            row = json.loads(line)
            rows.append(
                (
                    int(row['iter']) + offset,
                    float(row['energy']) / denom,
                    float(row['energy_error']) / denom,
                )
            )
    if tail is not None:
        rows = rows[-tail:]
    iters = np.array([row[0] for row in rows])
    energies = np.array([row[1] for row in rows])
    errors = np.array([row[2] for row in rows])
    return iters, energies, errors


# ──────────────────────────────────────────────────────────────────────
# Figure 1: Energy convergence for 8x8 J2=0 (D=5 and D=8, 30-iter runs)
# ──────────────────────────────────────────────────────────────────────
def fig1_j2_0_convergence():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    for D, color, marker in [(5, '#2196F3', 'o'), (8, '#E91E63', 's')]:
        fname = f'8x8_J2=0_D{D}_vmc_trajectory.csv'
        it, E, err, _ = load_csv(fname)
        Esite = E / 64
        errsite = err / 64
        ax1.errorbar(it, Esite, yerr=errsite, fmt=f'{marker}-', color=color,
                     markersize=4, linewidth=1.5, capsize=2,
                     label=f'D={D} VMC trajectory')
        cu_e = MEASURE[(0.0, D, 'CU')][0]
        ax1.axhline(cu_e, color=color, linestyle=':', alpha=0.5,
                    label=f'D={D} CU baseline ({cu_e:.5f})')

    ax1.axhline(QMC_J2_0[8], color='black', linestyle='--', linewidth=1.5,
                label=f'QMC ({QMC_J2_0[8]:.6f})')
    ax1.set_xlabel('VMC Iteration')
    ax1.set_ylabel('Energy per site')
    ax1.set_title('8×8 J2=0 PBC: VMC Convergence')
    ax1.legend(loc='upper right', fontsize=9)
    ax1.set_xlim(-0.5, 30)

    for D, color, marker in [(5, '#2196F3', 'o'), (8, '#E91E63', 's')]:
        fname = f'8x8_J2=0_D{D}_vmc_trajectory.csv'
        it, _, _, grad = load_csv(fname)
        ax2.semilogy(it, grad, f'{marker}-', color=color,
                     markersize=4, linewidth=1.5, label=f'D={D}')

    ax2.set_xlabel('VMC Iteration')
    ax2.set_ylabel('Gradient norm')
    ax2.set_title('8×8 J2=0: Gradient Norm')
    ax2.legend()
    ax2.set_xlim(-0.5, 30)

    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, 'fig1_8x8_J2=0_convergence.png'))
    plt.close(fig)
    print('  fig1 done')


# ──────────────────────────────────────────────────────────────────────
# Figure 2: Energy convergence for 8x8 J2=0.5 (D=5 complete + D=8 partial)
# ──────────────────────────────────────────────────────────────────────
def fig2_j2_05_convergence():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    it5, E5, err5, grad5 = load_csv('8x8_J2=0.5_D5_vmc_trajectory.csv')
    ax1.errorbar(it5, E5/64, yerr=err5/64, fmt='o-', color='#2196F3',
                 markersize=4, linewidth=1.5, capsize=2,
                 label='D=5 CU→VMC (30 iters, done)')

    it8y, E8y, err8y, grad8y = load_csv('8x8_J2=0.5_D8_yubin_vmc_partial.csv')
    ax1.errorbar(it8y, E8y/64, yerr=err8y/64, fmt='s-', color='#E91E63',
                 markersize=5, linewidth=1.5, capsize=2,
                 label=f'D=8 CU→VMC ({len(it8y)-1}/30 iters, running)')

    it8s, E8s, err8s, grad8s = load_csv('8x8_J2=0.5_D8_su_vmc_partial.csv')
    ax1.errorbar(it8s, E8s/64, yerr=err8s/64, fmt='^-', color='#FF9800',
                 markersize=5, linewidth=1.5, capsize=2,
                 label=f'D=8 SU→VMC ({len(it8s)-1}/30 iters, running)')

    ax1.axhline(MEASURE[(0.5, 5, 'CU')][0], color='#2196F3', linestyle=':', alpha=0.5,
                label=f'D=5 CU baseline ({MEASURE[(0.5, 5, "CU")][0]:.5f})')
    ax1.axhline(MEASURE[(0.5, 8, 'CU')][0], color='#E91E63', linestyle=':', alpha=0.5,
                label=f'D=8 CU baseline ({MEASURE[(0.5, 8, "CU")][0]:.5f})')

    ax1.set_xlabel('VMC Iteration')
    ax1.set_ylabel('Energy per site')
    ax1.set_title('8×8 J2=0.5 PBC: VMC Convergence')
    ax1.legend(loc='upper right', fontsize=8)

    ax2.semilogy(it5, grad5, 'o-', color='#2196F3', markersize=4, linewidth=1.5, label='D=5')
    ax2.semilogy(it8y, grad8y, 's-', color='#E91E63', markersize=5, linewidth=1.5, label='D=8 CU→VMC')
    ax2.semilogy(it8s, grad8s, '^-', color='#FF9800', markersize=5, linewidth=1.5, label='D=8 SU→VMC')
    ax2.set_xlabel('VMC Iteration')
    ax2.set_ylabel('Gradient norm')
    ax2.set_title('8×8 J2=0.5: Gradient Norm')
    ax2.legend()

    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, 'fig2_8x8_J2=0.5_convergence.png'))
    plt.close(fig)
    print('  fig2 done')


# ──────────────────────────────────────────────────────────────────────
# Figure 3: VMC improvement bar chart (measured energies)
# ──────────────────────────────────────────────────────────────────────
def fig3_energy_comparison():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

    # J2 = 0
    labels_j0 = ['D=5\nCU', 'D=5\nCU+VMC', 'D=8\nCU', 'D=8\nCU+VMC']
    vals_j0 = [MEASURE[(0.0, 5, 'CU')], MEASURE[(0.0, 5, 'CU+VMC')],
               MEASURE[(0.0, 8, 'CU')], MEASURE[(0.0, 8, 'CU+VMC')]]
    colors_j0 = ['#90CAF9', '#2196F3', '#F48FB1', '#E91E63']
    y_j0 = [v[0] for v in vals_j0]
    yerr_j0 = [v[1] for v in vals_j0]

    bars = ax1.bar(labels_j0, y_j0, yerr=yerr_j0, capsize=5, color=colors_j0,
                   edgecolor='gray', linewidth=0.5)
    ax1.axhline(QMC_J2_0[8], color='black', linestyle='--', linewidth=1.5,
                label=f'QMC = {QMC_J2_0[8]:.6f}')
    ax1.set_ylabel('Energy per site')
    ax1.set_title('8×8 J2=0 PBC')
    ax1.legend(loc='upper left')
    ax1.set_ylim(-0.6740, -0.6710)
    for bar, val in zip(bars, y_j0):
        ax1.text(bar.get_x() + bar.get_width()/2, val + 0.00008,
                 f'{val:.5f}', ha='center', va='bottom', fontsize=9, fontweight='bold')

    # J2 = 0.5 — D=8 now at iter 26: E/site = -31.880758/64 = -0.498137
    e_d8_j05_current = -31.880758 / 64
    labels_j05 = ['D=5\nCU', 'D=5\nCU+VMC', 'D=8\nCU', f'D=8\nCU+VMC\n(iter 26)']
    y_j05 = [MEASURE[(0.5, 5, 'CU')][0], MEASURE[(0.5, 5, 'CU+VMC')][0],
             MEASURE[(0.5, 8, 'CU')][0], e_d8_j05_current]
    yerr_j05 = [MEASURE[(0.5, 5, 'CU')][1], MEASURE[(0.5, 5, 'CU+VMC')][1],
                MEASURE[(0.5, 8, 'CU')][1], 0.00005]
    colors_j05 = ['#90CAF9', '#2196F3', '#F48FB1', '#E91E63']

    bars2 = ax2.bar(labels_j05, y_j05, yerr=yerr_j05, capsize=5, color=colors_j05,
                    edgecolor='gray', linewidth=0.5)
    bars2[3].set_hatch('//')
    bars2[3].set_alpha(0.7)
    ax2.set_ylabel('Energy per site')
    ax2.set_title('8×8 J2=0.5 PBC')
    ax2.set_ylim(-0.4990, -0.4935)
    for bar, val in zip(bars2, y_j05):
        ax2.text(bar.get_x() + bar.get_width()/2, val + 0.0001,
                 f'{val:.5f}', ha='center', va='bottom', fontsize=9, fontweight='bold')
    ax2.text(bars2[3].get_x() + bars2[3].get_width()/2, y_j05[3] - 0.0003,
             '(running)', ha='center', va='top', fontsize=8, fontstyle='italic', color='gray')

    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, 'fig3_energy_comparison.png'))
    plt.close(fig)
    print('  fig3 done')


# ──────────────────────────────────────────────────────────────────────
# Figure 4: Large systems — 16x16 trajectories + finite-size comparison
# ──────────────────────────────────────────────────────────────────────
def fig4_large_system():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Left: 16x16 J2=0 and J2=0.5 early trajectories
    it16j0, E16j0, err16j0, _ = load_csv('16x16_J2=0_D8_vmc_partial.csv')
    ax1.errorbar(it16j0, E16j0/256, yerr=err16j0/256, fmt='D-', color='#4CAF50',
                 markersize=8, linewidth=1.5, capsize=3,
                 label=f'16×16 J2=0 D=8 (iter 0: {E16j0[0]/256:.5f})')
    ax1.axhline(QMC_J2_0[16], color='#4CAF50', linestyle='--', linewidth=1.2,
                label=f'QMC J2=0 = {QMC_J2_0[16]:.6f}')

    it16j05, E16j05, err16j05, _ = load_csv('16x16_J2=0.5_D8_vmc_partial.csv')
    ax1.errorbar(it16j05, E16j05/256, yerr=err16j05/256, fmt='s-', color='#9C27B0',
                 markersize=8, linewidth=1.5, capsize=3,
                 label=f'16×16 J2=0.5 D=8 (iter 0: {E16j05[0]/256:.5f})')

    ax1.set_xlabel('VMC Iteration')
    ax1.set_ylabel('Energy per site')
    ax1.set_title('16×16 D=8 PBC: Early Convergence\n(2-node, 112 ranks)')
    ax1.legend(fontsize=9)
    ax1.set_xlim(-0.5, 3)

    # Right: J2=0 finite-size comparison
    ax2.plot([8, 12, 16], [QMC_J2_0[L] for L in [8, 12, 16]],
             'k--D', markersize=8, linewidth=1.5, label='QMC')
    ax2.plot([8], [MEASURE[(0.0, 8, 'CU+VMC')][0]], 's', color='#E91E63', markersize=10,
             label=f'D=8 CU+VMC 30-iter ({MEASURE[(0.0, 8, "CU+VMC")][0]:.5f})')
    ax2.plot([8], [MEASURE[(0.0, 8, 'CU')][0]], 'o', color='#F48FB1', markersize=8,
             label=f'D=8 CU baseline')
    ax2.plot([16], [E16j0[0]/256], '^', color='#4CAF50', markersize=10,
             label=f'16×16 iter 0 ({E16j0[0]/256:.5f})')
    ax2.plot([8], [MEASURE[(0.0, 5, 'CU+VMC')][0]], 'v', color='#2196F3', markersize=8,
             label=f'D=5 CU+VMC ({MEASURE[(0.0, 5, "CU+VMC")][0]:.5f})')

    ax2.set_xlabel('L (L×L lattice)')
    ax2.set_ylabel('Energy per site')
    ax2.set_title('J2=0 PBC: Finite-Size Scaling')
    ax2.legend(fontsize=8)
    ax2.set_xticks([8, 12, 16])

    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, 'fig4_large_system.png'))
    plt.close(fig)
    print('  fig4 done')


# ──────────────────────────────────────────────────────────────────────
# Figure 5: D=5 vs D=8 comparison at 8x8 J2=0
# ──────────────────────────────────────────────────────────────────────
def fig5_D_comparison():
    fig, ax = plt.subplots(figsize=(8, 5))

    it5, E5, err5, _ = load_csv('8x8_J2=0_D5_vmc_trajectory.csv')
    it8, E8, err8, _ = load_csv('8x8_J2=0_D8_vmc_trajectory.csv')

    ax.errorbar(it5, E5/64, yerr=err5/64, fmt='o-', color='#2196F3',
                markersize=3, linewidth=1.2, capsize=1.5, alpha=0.8, label='J2=0 D=5')
    ax.errorbar(it8, E8/64, yerr=err8/64, fmt='s-', color='#E91E63',
                markersize=3, linewidth=1.2, capsize=1.5, alpha=0.8, label='J2=0 D=8')

    ax.axhline(QMC_J2_0[8], color='black', linestyle='--', linewidth=1.5,
               label=f'QMC = {QMC_J2_0[8]:.6f}')
    ax.axhline(MEASURE[(0.0, 5, 'CU+VMC')][0], color='#2196F3', linestyle=':', alpha=0.4)
    ax.axhline(MEASURE[(0.0, 8, 'CU+VMC')][0], color='#E91E63', linestyle=':', alpha=0.4)

    ax.set_xlabel('VMC Iteration')
    ax.set_ylabel('Energy per site')
    ax.set_title('8×8 J2=0 PBC: D=5 vs D=8 VMC Convergence')
    ax.legend(loc='upper right')
    ax.annotate(f'D=5 meas: {MEASURE[(0.0, 5, "CU+VMC")][0]:.5f}',
                xy=(29, MEASURE[(0.0, 5, 'CU+VMC')][0]), fontsize=9,
                xytext=(20, -0.6726), arrowprops=dict(arrowstyle='->', color='#2196F3'),
                color='#2196F3')
    ax.annotate(f'D=8 meas: {MEASURE[(0.0, 8, "CU+VMC")][0]:.5f}',
                xy=(29, MEASURE[(0.0, 8, 'CU+VMC')][0]), fontsize=9,
                xytext=(20, -0.6738), arrowprops=dict(arrowstyle='->', color='#E91E63'),
                color='#E91E63')

    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, 'fig5_D_comparison_J2=0.png'))
    plt.close(fig)
    print('  fig5 done')


# ──────────────────────────────────────────────────────────────────────
# Figure 6: Summary energy table
# ──────────────────────────────────────────────────────────────────────
def fig6_summary_table():
    fig, ax = plt.subplots(figsize=(13, 5.5))
    ax.axis('off')

    columns = ['System', 'J2', 'D', 'Method', 'E/site', 'Error', 'vs QMC', 'Status']
    data = [
        ['8×8', '0', '5', 'CU baseline', '-0.67181', '3e-5', '-0.25%', 'Done'],
        ['8×8', '0', '5', 'CU + 30 VMC', '-0.67311', '2e-5', '-0.06%', 'Done'],
        ['8×8', '0', '8', 'CU baseline', '-0.67188', '3e-4', '-0.24%', 'Done'],
        ['8×8', '0', '8', 'CU + 30 VMC', '-0.67332', '2e-5', '-0.02%', 'Done'],
        ['8×8', '0', '8', 'CU + StepSel (60/60)', '-0.67349', '8e-4', '<0.01%', 'Done'],
        ['8×8', '0.5', '5', 'CU baseline', '-0.49490', '6e-5', '--', 'Done'],
        ['8×8', '0.5', '5', 'CU + 30 VMC', '-0.49779', '5e-5', '--', 'Done'],
        ['8×8', '0.5', '8', 'CU baseline', '-0.49648', '3e-4', '--', 'Done'],
        ['8×8', '0.5', '8', 'CU + 30 VMC', '-0.49814', '~5e-5', '--', 'Done'],
        ['8×8', '0.5', '8', 'SU + 30 VMC', '-0.49752', '~6e-5', '--', 'Done'],
        ['8×8', '0.5', '8', 'CU + StepSel (60/60)', '-0.49856', '1e-2', '--', 'Done'],
        ['8×8', '0', '8', 'MinSR lr=0.1 (60/60)', '-0.67340', '4e-3', '-0.01%', 'Done'],
        ['8×8', '0.5', '8', 'MinSR lr=0.1 (60/60)', '-0.49855', '1e-2', '--', 'Done'],
        ['12×12', '0', '8', 'CU + SR (12/60)', '-0.67007', '6e-3', '-0.09%', 'Running'],
        ['16×16', '0', '8', 'CU + SR (60/60)', '-0.66967', '5e-3', '+0.05%', 'Done'],
        ['16×16', '0', '8', 'CU + MinSR (60/60)', '-0.66967', '5e-3', '+0.05%', 'Done'],
        ['16×16', '0.5', '8', 'CU + SR (60/60)', '-0.49588', '5e-3', '--', 'Done'],
        ['16×16', '0.5', '8', 'CU + MinSR (60/60)', '-0.49578', '5e-3', '--', 'Done'],
        ['16×16', '0.5', '8', 'CU + VMC (iter 0)', '-0.49505', '~3e-5', '--', 'Crashed'],
    ]

    colors_cell = []
    for row in data:
        row_colors = ['white'] * len(row)
        if row[-1] == 'Running':
            row_colors = ['#FFF9C4'] * len(row)
        elif row[-1] == 'Failed':
            row_colors = ['#FFCDD2'] * len(row)
        elif row[-1] in ('Resubmit', 'Submitted'):
            row_colors = ['#FFE0B2'] * len(row)
        elif row[-1] == 'Crashed':
            row_colors = ['#FFCDD2'] * len(row)
        if 'vs QMC' in columns:
            idx = columns.index('vs QMC')
            if row[idx] not in ['--', '']:
                try:
                    pct = float(row[idx].replace('%', ''))
                    if abs(pct) < 0.1:
                        row_colors[idx] = '#C8E6C9'
                except ValueError:
                    pass
        colors_cell.append(row_colors)

    table = ax.table(cellText=data, colLabels=columns, loc='center',
                     cellColours=colors_cell,
                     colColours=['#E3F2FD']*len(columns))
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.5)

    ax.set_title('PEPS+VMC Energy Results Summary (PBC, Square Heisenberg) — 2026-03-16',
                 fontsize=13, fontweight='bold', pad=20)

    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, 'fig6_summary_table.png'))
    plt.close(fig)
    print('  fig6 done')


# ──────────────────────────────────────────────────────────────────────
# Figure 7: StepSelector comparison for 8x8 J2=0 D=8
# ──────────────────────────────────────────────────────────────────────
def fig7_stepselector_comparison():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # Old 30-iter run (no initial step selector)
    it_old, E_old, err_old, grad_old = load_csv('8x8_J2=0_D8_vmc_trajectory.csv')
    ax1.errorbar(it_old, E_old/64, yerr=err_old/64, fmt='s-', color='#E91E63',
                 markersize=4, linewidth=1.5, capsize=2, alpha=0.8,
                 label='30-iter, no StepSelector (lr=0.1)')

    # New 60-iter run with step selector (16 iters so far)
    it_ss, E_ss, err_ss, grad_ss = load_csv('8x8_J2=0_D8_stepselector_vmc.csv')

    # Load lr column for coloring
    lrs = []
    with open(os.path.join(DATADIR, '8x8_J2=0_D8_stepselector_vmc.csv')) as f:
        reader = csv.DictReader(f)
        for row in reader:
            lrs.append(float(row['lr']))
    lrs = np.array(lrs)

    # Plot lr phases separately
    mask03 = lrs == 0.3
    mask015 = lrs == 0.15
    mask0075 = lrs == 0.075
    ax1.errorbar(it_ss[mask03], E_ss[mask03]/64, yerr=err_ss[mask03]/64, fmt='D-', color='#4CAF50',
                 markersize=5, linewidth=1.5, capsize=2,
                 label=f'StepSelector: lr=0.3 phase')
    ax1.errorbar(it_ss[mask015], E_ss[mask015]/64, yerr=err_ss[mask015]/64, fmt='D-', color='#1565C0',
                 markersize=5, linewidth=1.5, capsize=2,
                 label=f'StepSelector: lr=0.15 phase')
    if mask0075.any():
        ax1.errorbar(it_ss[mask0075], E_ss[mask0075]/64, yerr=err_ss[mask0075]/64, fmt='D-', color='#0D47A1',
                     markersize=5, linewidth=1.5, capsize=2,
                     label=f'StepSelector: lr=0.075 phase')

    ax1.axhline(QMC_J2_0[8], color='black', linestyle='--', linewidth=1.5,
                label=f'QMC = {QMC_J2_0[8]:.6f}')
    ax1.axhline(MEASURE[(0.0, 8, 'CU')][0], color='gray', linestyle=':', alpha=0.6,
                label=f'CU baseline ({MEASURE[(0.0, 8, "CU")][0]:.5f})')

    # Annotate LR halving at iter 10
    ax1.axvline(10, color='#1565C0', linestyle=':', alpha=0.5)
    ax1.annotate('AutoSelector:\nlr 0.3→0.15',
                 xy=(10, E_ss[10]/64), xytext=(12, -0.6725),
                 arrowprops=dict(arrowstyle='->', color='#1565C0'),
                 fontsize=9, color='#1565C0')

    ax1.set_xlabel('VMC Iteration')
    ax1.set_ylabel('Energy per site')
    ax1.set_title(f'8×8 J2=0 D=8: StepSelector ({len(it_ss)}/60 iters)')
    ax1.legend(fontsize=8, loc='upper right')
    ax1.set_xlim(-0.5, max(len(it_old), len(it_ss)) + 1)

    # Right: gradient norm comparison
    ax2.semilogy(it_old, grad_old, 's-', color='#E91E63',
                 markersize=4, linewidth=1.5, alpha=0.8,
                 label='No StepSelector (lr=0.1)')
    ax2.semilogy(it_ss[mask03], grad_ss[mask03], 'D-', color='#4CAF50',
                 markersize=5, linewidth=1.5,
                 label='StepSelector lr=0.3')
    ax2.semilogy(it_ss[mask015], grad_ss[mask015], 'D-', color='#1565C0',
                 markersize=5, linewidth=1.5,
                 label='StepSelector lr=0.15')
    if mask0075.any():
        ax2.semilogy(it_ss[mask0075], grad_ss[mask0075], 'D-', color='#0D47A1',
                     markersize=5, linewidth=1.5,
                     label='StepSelector lr=0.075')
    ax2.axvline(10, color='#1565C0', linestyle=':', alpha=0.5)

    ax2.set_xlabel('VMC Iteration')
    ax2.set_ylabel('Gradient norm')
    ax2.set_title('Gradient Norm Comparison')
    ax2.legend(fontsize=9)
    ax2.set_xlim(-0.5, max(len(it_old), len(it_ss)) + 1)

    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, 'fig7_stepselector_comparison.png'))
    plt.close(fig)
    print('  fig7 done')


# ──────────────────────────────────────────────────────────────────────
# Figure 8: Comparison with Yubin's draft results
# ──────────────────────────────────────────────────────────────────────

# Yubin's data (extracted from yubin_draft/main_prb.tex and figures)
YUBIN_J2_0 = {  # L: (E/site GO D=5 Dc=20)
    6: -0.6787, 8: -0.6733, 12: -0.6704, 16: -0.6696,
}
YUBIN_J2_0_CU = {8: -0.6720}  # CU baseline from text
YUBIN_J2_0_TDL = -0.6689  # from heis_fss.png

YUBIN_J2_05_D8 = -0.4978   # from j05l8D.png, L=8 D=8 converged
YUBIN_J2_05_D9 = -0.4980   # from j05l8D.png, L=8 D=9 converged
YUBIN_J2_05_FSS = {  # from E-Dc.png at Dc=20 (approximate)
    4: -0.529, 6: -0.503, 8: -0.499, 12: -0.497, 16: -0.496,
}
YUBIN_J2_05_TDL_PBC = -0.4951  # from fss4-24.png
YUBIN_J2_05_TDL_OBC = -0.4962

# QMC benchmarks
QMC_J2_0_TDL = -0.669437  # from heis_fss.png


def load_simple_csv(fname):
    """Load CSV with iteration and energy_per_site columns."""
    path = os.path.join(DATADIR, fname)
    iters, energies = [], []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            iters.append(int(row['iteration']))
            energies.append(float(row['energy_per_site']))
    return np.array(iters), np.array(energies)


def fig8_yubin_comparison():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

    # ── Left: J2=0 finite-size comparison ──
    Ls = [6, 8, 12, 16]
    qmc_vals = [QMC_J2_0.get(L, None) for L in Ls]

    # Yubin GO D=5
    yubin_vals = [YUBIN_J2_0[L] for L in Ls]
    ax1.plot(Ls, yubin_vals, 'v-', color='#9C27B0', markersize=9, linewidth=1.5,
             label='Yubin GO D=5 Dc=20')

    # QMC
    ax1.plot(Ls, qmc_vals, 'k--D', markersize=8, linewidth=1.5, label='QMC')

    # Our results
    ax1.plot([8], [MEASURE[(0.0, 5, 'CU+VMC')][0]], 'o', color='#2196F3', markersize=10,
             label=f'Ours: D=5 SR 30-iter ({MEASURE[(0.0, 5, "CU+VMC")][0]:.5f})')
    ax1.plot([8], [MEASURE[(0.0, 8, 'CU+VMC')][0]], 's', color='#E91E63', markersize=10,
             label=f'Ours: D=8 SR 30-iter ({MEASURE[(0.0, 8, "CU+VMC")][0]:.5f})')
    it_ss, E_ss, _, _ = load_csv('8x8_J2=0_D8_stepselector_vmc.csv')
    best_ss = min(E_ss/64)
    ax1.plot([8], [best_ss], '^', color='#4CAF50', markersize=10,
             label=f'Ours: D=8 StepSel best ({best_ss:.5f})')

    sr16_it, sr16_e, _ = load_jsonl_trace(
        os.path.join(os.path.dirname(__file__), '..', 'run', '16x16J2=0D8', '20260319_icpx_sr_fresh', 'vmc', 'energy', 'optimization_log.jsonl'),
        256.0,
    )
    ms16_it, ms16_e, _ = load_jsonl_trace(
        os.path.join(os.path.dirname(__file__), '..', 'run', '16x16J2=0D8', '20260319_icpx_minsr_fresh', 'vmc', 'energy', 'optimization_log.jsonl'),
        256.0,
    )
    best16 = min(np.min(sr16_e), np.min(ms16_e))
    ax1.plot([16], [best16], 'P', color='#FF9800', markersize=10,
             label=f'Ours: 16×16 SR/MinSR best ({best16:.5f})')

    ax1.set_xlabel('L (L×L lattice)')
    ax1.set_ylabel('Energy per site')
    ax1.set_title('$J_2=0$ PBC: Ours vs Yubin')
    ax1.legend(fontsize=8, loc='lower right')
    ax1.set_xticks([6, 8, 12, 16])
    ax1.set_xlim(5, 17)

    # ── Right: J2=0.5 optimization trajectory comparison ──
    # Yubin D=8 GO
    yit8, ye8 = load_simple_csv('yubin_j05_D8_trajectory.csv')
    ax2.plot(yit8, ye8, 'o-', color='#CE93D8', markersize=4, linewidth=1.5,
             label='Yubin GO D=8 (15 steps)')

    # Yubin D=9 GO
    yit9, ye9 = load_simple_csv('yubin_j05_D9_trajectory.csv')
    ax2.plot(yit9, ye9, 'v-', color='#9C27B0', markersize=3, linewidth=1.2, alpha=0.7,
             label='Yubin GO D=9 (40 steps)')

    # Ours: CU→SR D=8
    it8y, E8y, err8y, _ = load_csv('8x8_J2=0.5_D8_yubin_vmc_partial.csv')
    ax2.errorbar(it8y, E8y/64, yerr=err8y/64, fmt='s-', color='#E91E63',
                 markersize=4, linewidth=1.5, capsize=2,
                 label=f'Ours: D=8 CU→SR (30 iters)')

    # Ours: SU→SR D=8
    it8s, E8s, err8s, _ = load_csv('8x8_J2=0.5_D8_su_vmc_partial.csv')
    ax2.errorbar(it8s, E8s/64, yerr=err8s/64, fmt='^-', color='#FF9800',
                 markersize=4, linewidth=1.5, capsize=2,
                 label=f'Ours: D=8 SU→SR (30 iters)')

    # Ours: StepSelector D=8
    it_ss05, E_ss05, err_ss05, _ = load_csv('8x8_J2=0.5_D8_stepselector_vmc.csv')
    mask_03 = it_ss05 < 30
    ax2.errorbar(it_ss05[mask_03], E_ss05[mask_03]/64, yerr=err_ss05[mask_03]/64,
                 fmt='D-', color='#4CAF50', markersize=6, linewidth=1.8, capsize=2,
                 label='Ours: D=8 StepSel + job 560212')
    it_560212, E_560212, err_560212 = load_jsonl_trace(
        os.path.join(os.path.dirname(__file__), '..', 'run', '8x8J2=0.5D8', '20260307_yubin_vmc_stepselector', 'vmc', 'energy', 'optimization_log.jsonl'),
        64.0,
        tail=30,
        offset=30,
    )
    ax2.errorbar(it_560212, E_560212, yerr=err_560212,
                 fmt='D-', color='#4CAF50', markersize=5, linewidth=1.8, capsize=2,
                 label='_nolegend_')
    ax2.axvline(29.5, color='gray', linestyle='--', alpha=0.4, linewidth=1)
    ax2.text(30.5, -0.49555, 'job 560212 starts\n12824 samples, lr=0.1', fontsize=8, color='gray', va='top')

    # Reference lines
    ax2.axhline(YUBIN_J2_05_D8, color='#CE93D8', linestyle=':', alpha=0.5)
    ax2.axhline(-0.49814, color='#E91E63', linestyle=':', alpha=0.5)

    # Best energy shown in the canonical merged figure
    ss_all_it = np.concatenate([it_ss05[mask_03], it_560212])
    ss_all_e = np.concatenate([E_ss05[mask_03]/64, E_560212])
    best_idx = np.argmin(ss_all_e)
    best_e = ss_all_e[best_idx]
    best_it = ss_all_it[best_idx]
    ax2.annotate(f'Yubin D=8: {YUBIN_J2_05_D8:.4f}',
                 xy=(14, YUBIN_J2_05_D8), xytext=(18, -0.4965),
                 arrowprops=dict(arrowstyle='->', color='#CE93D8'),
                 fontsize=9, color='#9C27B0')
    ax2.annotate(f'Best shown: {best_e:.5f} (job 560212)',
                 xy=(best_it, best_e), xytext=(40, -0.4995),
                 arrowprops=dict(arrowstyle='->', color='#2E7D32'),
                 fontsize=9, color='#2E7D32', fontweight='bold')

    ax2.set_xlabel('Optimization Step')
    ax2.set_ylabel('Energy per site')
    ax2.set_title('$J_2=0.5$ 8×8 D=8 PBC: Trajectory Comparison')
    ax2.legend(fontsize=8, loc='upper right')
    ax2.set_xlim(-1, 62)
    ax2.set_ylim(-0.5000, -0.4953)

    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, 'fig8_yubin_comparison.png'))
    plt.close(fig)
    print('  fig8 done')


# ──────────────────────────────────────────────────────────────────────
# Figure 10: MinSR optimizer comparison (8x8 J2=0 D=8)
# ──────────────────────────────────────────────────────────────────────
def load_minsr_csv(fname):
    path = os.path.join(DATADIR, fname)
    iters, energies, errors, grads, sr_ngrads = [], [], [], [], []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            iters.append(int(row['iteration']))
            energies.append(float(row['energy']))
            errors.append(float(row['energy_error']))
            grads.append(float(row['gradient_norm']))
            sr_ngrads.append(float(row['sr_ngrad_norm']))
    return (np.array(iters), np.array(energies), np.array(errors),
            np.array(grads), np.array(sr_ngrads))


def fig10_minsr_comparison():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

    # MinSR lr=0.3 (diverges)
    it03, E03, err03, _, ngrad03 = load_minsr_csv('8x8_J2=0_D8_minsr_lr03.csv')
    # MinSR lr=0.1 (converges)
    it01, E01, err01, _, ngrad01 = load_minsr_csv('8x8_J2=0_D8_minsr_lr01.csv')

    # SR baseline (30-iter, lr=0.1)
    it_sr, E_sr, err_sr, _ = load_csv('8x8_J2=0_D8_vmc_trajectory.csv')

    # StepSelector (best result)
    it_ss, E_ss, err_ss, _ = load_csv('8x8_J2=0_D8_stepselector_vmc.csv')

    # Left: Energy trajectory
    ax1.errorbar(it03, E03/64, yerr=err03/64, fmt='x-', color='#F44336',
                 markersize=4, linewidth=1.2, capsize=1.5, alpha=0.8,
                 label=f'MinSR lr=0.3 (diverges)')
    ax1.errorbar(it01, E01/64, yerr=err01/64, fmt='D-', color='#4CAF50',
                 markersize=3, linewidth=1.5, capsize=1.5,
                 label=f'MinSR lr=0.1 (final: {E01[-1]/64:.5f})')
    ax1.errorbar(it_sr, E_sr/64, yerr=err_sr/64, fmt='s-', color='#E91E63',
                 markersize=3, linewidth=1.2, capsize=1.5, alpha=0.7,
                 label=f'SR lr=0.1 (final: {E_sr[-1]/64:.5f})')
    ax1.plot(it_ss, E_ss/64, '-', color='#1565C0', linewidth=1.0, alpha=0.5,
             label=f'StepSel best: {min(E_ss/64):.5f}')

    ax1.axhline(QMC_J2_0[8], color='black', linestyle='--', linewidth=1.5,
                label=f'QMC = {QMC_J2_0[8]:.6f}')

    ax1.set_xlabel('Optimization Step')
    ax1.set_ylabel('Energy per site')
    ax1.set_title('8×8 $J_2=0$ D=8: MinSR vs SR Optimizer')
    ax1.legend(fontsize=8, loc='upper right')
    ax1.set_xlim(-1, 61)
    ax1.set_ylim(-0.6750, -0.6500)

    # Right: Zoomed convergence region (MinSR lr=0.1 vs SR)
    ax2.errorbar(it01, E01/64, yerr=err01/64, fmt='D-', color='#4CAF50',
                 markersize=4, linewidth=1.5, capsize=2,
                 label=f'MinSR lr=0.1')
    ax2.errorbar(it_sr, E_sr/64, yerr=err_sr/64, fmt='s-', color='#E91E63',
                 markersize=4, linewidth=1.5, capsize=2, alpha=0.8,
                 label=f'SR lr=0.1')
    ax2.plot(it_ss, E_ss/64, '-', color='#1565C0', linewidth=1.0, alpha=0.5,
             label=f'StepSelector')

    ax2.axhline(QMC_J2_0[8], color='black', linestyle='--', linewidth=1.5,
                label=f'QMC = {QMC_J2_0[8]:.6f}')

    # Annotate final values
    minsr_final = np.mean(E01[-10:]/64)
    sr_final = E_sr[-1]/64
    ax2.annotate(f'MinSR avg: {minsr_final:.5f}',
                 xy=(55, minsr_final), xytext=(35, -0.6728),
                 fontsize=9, color='#4CAF50',
                 arrowprops=dict(arrowstyle='->', color='#4CAF50'))
    ax2.annotate(f'SR final: {sr_final:.5f}',
                 xy=(29, sr_final), xytext=(10, -0.6726),
                 fontsize=9, color='#E91E63',
                 arrowprops=dict(arrowstyle='->', color='#E91E63'))

    ax2.set_xlabel('Optimization Step')
    ax2.set_ylabel('Energy per site')
    ax2.set_title('Zoomed: Convergence Comparison')
    ax2.legend(fontsize=8, loc='upper right')
    ax2.set_xlim(-1, 61)
    ax2.set_ylim(-0.6742, -0.6710)

    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, 'fig10_minsr_comparison.png'))
    plt.close(fig)
    print('  fig10 done')


# ──────────────────────────────────────────────────────────────────────
# Figure 11: 16×16 J2=0 trajectory comparison (Yubin + SR/MinSR lr=0.1 + SR lr=0.2)
# ──────────────────────────────────────────────────────────────────────
def fig11_16x16_j2_0_trajectories():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

    N = 256.0  # 16x16

    # Yubin GO D=5 Dc=20
    yit, ye = load_simple_csv('yubin_heis_grad_L16_auto.csv')
    ax1.plot(yit, ye, 'v-', color='#9C27B0', markersize=5, linewidth=1.5,
             label=f'Yubin GO D=5 Dc=20 ({ye[-1]:.5f})')

    # Our SR lr=0.1 (60 iters)
    sr_it, sr_e, sr_err = load_jsonl_trace(
        os.path.join(os.path.dirname(__file__), '..', 'run', '16x16J2=0D8',
                     '20260319_icpx_sr_fresh', 'vmc', 'energy', 'optimization_log.jsonl'),
        N)
    ax1.errorbar(sr_it, sr_e, yerr=sr_err, fmt='s-', color='#E91E63',
                 markersize=3, linewidth=1.5, capsize=1.5,
                 label=f'SR lr=0.1, 1280 samp (best: {np.min(sr_e):.5f})')

    # Our MinSR lr=0.1 (60 iters)
    ms_it, ms_e, ms_err = load_jsonl_trace(
        os.path.join(os.path.dirname(__file__), '..', 'run', '16x16J2=0D8',
                     '20260319_icpx_minsr_fresh', 'vmc', 'energy', 'optimization_log.jsonl'),
        N)
    ax1.errorbar(ms_it, ms_e, yerr=ms_err, fmt='D-', color='#4CAF50',
                 markersize=3, linewidth=1.5, capsize=1.5,
                 label=f'MinSR lr=0.1, 1280 samp (best: {np.min(ms_e):.5f})')

    # Our SR lr=0.2, 2560 samples (60 iters)
    sr2_path = os.path.join(DATADIR, '16x16_J2=0_D8_sr_lr02_2560_vmc.csv')
    sr2_it, sr2_e, sr2_err, _ = [], [], [], []
    with open(sr2_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            sr2_it.append(int(row['iter']))
            sr2_e.append(float(row['energy']) / N)
            sr2_err.append(float(row['energy_error']) / N)
    sr2_it, sr2_e, sr2_err = np.array(sr2_it), np.array(sr2_e), np.array(sr2_err)
    ax1.errorbar(sr2_it, sr2_e, yerr=sr2_err, fmt='^-', color='#FF9800',
                 markersize=3, linewidth=1.5, capsize=1.5,
                 label=f'SR lr=0.2, 2560 samp (best: {np.min(sr2_e):.5f})')

    # QMC benchmark
    ax1.axhline(QMC_J2_0[16], color='black', linestyle='--', linewidth=1.5,
                label=f'QMC = {QMC_J2_0[16]:.6f}')

    ax1.set_xlabel('Optimization Step')
    ax1.set_ylabel('Energy per site')
    ax1.set_title('16×16 $J_2=0$ D=8 PBC: Trajectory Comparison')
    ax1.legend(fontsize=8, loc='upper right')
    ax1.set_xlim(-1, 62)

    # ── Right: zoomed view of convergence region ──
    ax2.errorbar(sr_it, sr_e, yerr=sr_err, fmt='s-', color='#E91E63',
                 markersize=4, linewidth=1.5, capsize=2,
                 label='SR lr=0.1, 1280 samp')
    ax2.errorbar(ms_it, ms_e, yerr=ms_err, fmt='D-', color='#4CAF50',
                 markersize=4, linewidth=1.5, capsize=2,
                 label='MinSR lr=0.1, 1280 samp')
    ax2.errorbar(sr2_it, sr2_e, yerr=sr2_err, fmt='^-', color='#FF9800',
                 markersize=4, linewidth=1.5, capsize=2,
                 label='SR lr=0.2, 2560 samp')

    ax2.axhline(QMC_J2_0[16], color='black', linestyle='--', linewidth=1.5,
                label=f'QMC = {QMC_J2_0[16]:.6f}')

    # Annotations for final averages (last 10 iters)
    sr_avg = np.mean(sr_e[-10:])
    ms_avg = np.mean(ms_e[-10:])
    sr2_avg = np.mean(sr2_e[-10:])
    ax2.annotate(f'SR lr=0.1 avg: {sr_avg:.5f}',
                 xy=(55, sr_avg), xytext=(30, sr_avg + 0.00015),
                 fontsize=9, color='#E91E63',
                 arrowprops=dict(arrowstyle='->', color='#E91E63'))
    ax2.annotate(f'SR lr=0.2 avg: {sr2_avg:.5f}',
                 xy=(55, sr2_avg), xytext=(30, sr2_avg - 0.00020),
                 fontsize=9, color='#FF9800',
                 arrowprops=dict(arrowstyle='->', color='#FF9800'))

    ax2.set_xlabel('Optimization Step')
    ax2.set_ylabel('Energy per site')
    ax2.set_title('Zoomed: Convergence Comparison')
    ax2.legend(fontsize=8, loc='upper right')
    ax2.set_xlim(-1, 62)
    # Auto-zoom to converged region
    all_e = np.concatenate([sr_e[5:], ms_e[5:], sr2_e[5:]])
    ymin, ymax = np.min(all_e) - 0.0003, np.max(all_e) + 0.0003
    ax2.set_ylim(ymin, ymax)

    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, 'fig11_16x16_j2_0_trajectories.png'))
    plt.close(fig)
    print('  fig11 done')


# ──────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    print('Generating plots...')
    fig1_j2_0_convergence()
    fig2_j2_05_convergence()
    fig3_energy_comparison()
    fig4_large_system()
    fig5_D_comparison()
    fig6_summary_table()
    fig7_stepselector_comparison()
    fig8_yubin_comparison()
    fig10_minsr_comparison()
    fig11_16x16_j2_0_trajectories()
    print(f'All figures saved to {FIGDIR}/')
