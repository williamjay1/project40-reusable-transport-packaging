"""Redraw the five v39 SI figures from frozen, auditable evidence.

No optimizer, resampling, smoothing or fitted probability is introduced here.
This file provides functions to redraw_supplement.py; it is not a standalone CLI
and does not inspect or compile the manuscript.
"""
from pathlib import Path
from decimal import Decimal, ROUND_HALF_UP
import argparse
import csv
import hashlib
import json
import math
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import findfont
from matplotlib.patches import FancyArrowPatch, Rectangle
from matplotlib.text import Text
import numpy as np
from PIL import Image

# Only the read-only audit is imported; never call the main save/main helpers.
from plot_epse_nature_v39 import layout_audit

ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT
OUT = ROOT / 'figures'
MANIFEST = ROOT / 'results/epse_v39/supplement_figure_provenance.json'
STYLE = ROOT / 'code/vendor/scienceplots_v39'
TEX = ROOT / 'figure_data/figure_semantics.tex'
STRUCTURAL = FROZEN / 'results/lifecycle_v38/extension/structural_results.json'
OWNERSHIP = ROOT / 'benchmarks/historical/calendar_thresholds_v23/strong_baseline_comparisons.csv'
COLLECTION = ROOT / 'benchmarks/historical/calendar_collection_v25/collection_energy_fate_frontier.csv'
BLUE, ORANGE, GREY = '#0072B2', '#D55E00', '#6B7075'
STYLE_NAMES = ['science.mplstyle', 'nature.mplstyle', 'no-latex.mplstyle']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def identity(path):
    result = dict(path=str(path), bytes=path.stat().st_size, sha256=sha(path))
    if path.is_relative_to(ROOT):
        result['project_relative_path'] = path.relative_to(ROOT).as_posix()
    return result


def read_csv(path):
    with path.open(encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream))


def panel(ax, letter, title):
    ax.text(-.02, 1.14, letter, transform=ax.transAxes, fontsize=8, fontweight='bold')
    ax.text(.045, 1.14, title, transform=ax.transAxes, fontsize=7)
    ax.tick_params(direction='out', length=2.6, width=.5, pad=3)


def key(fig, handles, x=.52):
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(x, .985),
               ncol=2, frameon=False, handlelength=1.8, columnspacing=2.1)


def ownership_key(fig):
    key(fig, [plt.Line2D([], [], color=ORANGE, marker='s', linewidth=1, label='Dedicated'),
              plt.Line2D([], [], color=BLUE, marker='o', linewidth=1, label='Shared')])


def si1(structural):
    # Exactly the v38 lifecycle_foreground_v38.py figure recipe, cross-checked
    # at all matching stored grid rows (L=2,4,8); all 486 stored checks pass.
    lives = list(range(1, 10))
    dedicated = [3 * math.ceil(3 / min(5, 1 + (L - 1) // 3)) for L in lives]
    pooled = [math.ceil(9 / min(5, L)) for L in lives]
    assert all(r['exact_births'] == r['formula'] for r in structural['formula_rows'])
    matching = [r for r in structural['formula_rows']
                if (r['N'], r['owners'], r['K'], r['spacing'], r['delay']) == (9, 3, 5, 1, 1)]
    assert len(matching) == 6
    for row in matching:
        values = pooled if row['ownership'] == 'pooled' else dedicated
        assert values[row['L'] - 1] == row['exact_births']
    boxes = {r['ownership']: r for r in structural['state_boxes']}
    for box in boxes.values():
        assert len(box['states']) == 162
        assert [min(r['new_births'] for r in box['states']),
                max(r['new_births'] for r in box['states'])] == [box['lower'], box['upper']]
    fig = plt.figure(figsize=(7, 2.85))
    a = fig.add_axes([.075, .23, .43, .49])
    b = fig.add_axes([.68, .23, .28, .49])
    for values, color, marker in [(dedicated, ORANGE, 's'), (pooled, BLUE, 'o')]:
        a.step(lives, values, where='post', color=color, linewidth=1)
        a.plot(lives, values, linestyle='none', marker=marker, color=color, markersize=3)
    a.set(xlabel='Ready-end asset life (intervals)', ylabel='New asset entries (count)',
          xlim=(.7, 9.3), ylim=(.5, 9.5), xticks=[1, 3, 5, 7, 9], yticks=[1, 3, 5, 7, 9])
    for mode, y, color, marker in [('dedicated', 0, ORANGE, 's'), ('pooled', 1, BLUE, 'o')]:
        box = boxes[mode]
        b.plot([box['lower'], box['upper']], [y, y], color=color, linewidth=1.6,
               marker=marker, markersize=3)
    b.set(xlabel='Additional new entries (count)', xlim=(-.2, 3.2), ylim=(-.45, 1.45),
          xticks=[0, 1, 2, 3], yticks=[0, 1], yticklabels=['Dedicated', 'Shared'])
    panel(a, 'a', 'Equal-spacing construction')
    panel(b, 'b', 'Admitted old-state ranges')
    ownership_key(fig)
    return fig, dict(lifetimes=lives, dedicated=dedicated, shared=pooled,
                     state_ranges={m: [boxes[m]['lower'], boxes[m]['upper']] for m in boxes},
                     conditions=dict(jobs=9, owners=3, K=5, spacing=1, delay=1, synthetic=True),
                     comparison=dict(original_recipe_18_values_identical=True,
                                     matching_stored_rows=len(matching), structural_grid_checks=486,
                                     state_settings=324, maximum_numeric_difference=0),
                     interpretation='Old-state extrema, not confidence intervals or observed histories')


def si2(rows):
    order = ['flat_tau1', 'staggered_tau1', 'staggered_tau3',
             'staggered_tau1_H24', 'staggered_tau1_H48']
    labels = ['Flat\n12 weeks, delay 1', 'Staggered\n12 weeks, delay 1',
              'Staggered\n12 weeks, delay 3', 'Staggered\n24 weeks, delay 1',
              'Staggered\n48 weeks, delay 1']
    by_name = {row['scenario']: row for row in rows}
    assert set(by_name) == set(order)
    fig = plt.figure(figsize=(7, 3.6))
    axes = [fig.add_axes([left, .19, .20, .55]) for left in [.22, .485, .755]]
    fields = ['initial_assets', 'cost_per_service', 'maximum_upper_impact_ratio']
    data = dict(scenarios=order, panels=[], comparison=dict(fields_checked=30, maximum_numeric_difference=0))
    for idx, (ax, field) in enumerate(zip(axes, fields)):
        record = dict(field=field, dedicated=[], shared=[], display_divisor=1000 if idx == 0 else 1)
        for mode, display_mode, offset, color in [('dedicated', 'dedicated', -.17, ORANGE),
                                                 ('pooled', 'shared', .17, BLUE)]:
            original = [float(by_name[n][f'{mode}_{field}']) for n in order]
            values = [v / record['display_divisor'] for v in original]
            ax.barh(np.arange(5) + offset, values, height=.31, color=color, linewidth=0)
            record[display_mode] = original
        ax.set(ylim=(4.6, -.6), yticks=np.arange(5),
               yticklabels=labels if idx == 0 else [''] * 5)
        ax.tick_params(axis='y', length=0)
        if idx == 0:
            ax.set(xlabel='Opening stock\n(thousand equivalents)', xlim=(0, 220), xticks=[0, 100, 200])
        elif idx == 1:
            ax.set(xlabel='Partial cost\n(EUR per service)', xlim=(0, 2.7), xticks=[0, 1, 2])
        else:
            ax.set(xlabel='Maximum upper\nimpact ratio', xlim=(0, 1.4), xticks=[0, .5, 1])
            ax.axvline(1, color=GREY, linestyle='--', linewidth=.7)
        panel(ax, 'abc'[idx], ['Stock', 'Partial cost', 'Impact'][idx])
        data['panels'].append(record)
    ownership_key(fig)
    return fig, data


def si3(rows):
    assert len(rows) == 124
    fig = plt.figure(figsize=(7, 2.9))
    axes = [fig.add_axes([.085, .24, .38, .48]), fig.add_axes([.585, .24, .37, .48])]
    data = dict(curves=[], comparison=dict(csv_rows=124, plotted_points=124,
                                          full_grid_retained=True, maximum_numeric_difference=0))
    for ax, delay in zip(axes, [1, 3]):
        for fate, color, marker in [(0, BLUE, 'o'), (1, GREY, 's')]:
            selected = sorted([r for r in rows if int(r['return_delay_weeks']) == delay and
                               float(r['missing_fate_virgin_production_equivalent_stress']) == fate],
                              key=lambda r: float(r['collection_fraction']))
            x = [float(r['collection_fraction']) for r in selected]
            y = [float(r['design_prior_95pct_maximum_ratio']) for r in selected]
            assert len(x) == 31 and x == [round(.70 + j / 100, 2) for j in range(31)]
            ax.plot(x, y, color=color, marker=marker, markersize=2, linewidth=1)
            first = min(xx for xx, yy in zip(x, y) if yy < 1)
            ax.axvline(first, color=color, linestyle=':', linewidth=.7)
            data['curves'].append(dict(delay_weeks=delay, fate_stress=fate,
                                       collection_fraction=x, design_prior_q95_maximum_ratio=y,
                                       first_grid_point_below_one=first))
        ax.axhline(1, color='black', linestyle='--', linewidth=.65)
        ax.set(xlabel='Declared collection fraction q', xlim=(.695, 1.005),
               ylim=(0, 4.4), xticks=[.7, .8, .9, 1], yticks=[0, 1, 2, 3, 4])
        panel(ax, 'a' if delay == 1 else 'b', f'Return delay: {delay} week' + ('s' if delay > 1 else ''))
    axes[0].set_ylabel('95th-percentile maximum\nupper impact ratio')
    key(fig, [plt.Line2D([], [], color=BLUE, marker='o', linewidth=1, label='No additional fate charge'),
              plt.Line2D([], [], color=GREY, marker='s', linewidth=1, label='Production-equivalent fate stress')])
    data['interpretation'] = 'Original declared-prior q95 and both fate panels; not fitted collection probabilities'
    return fig, data


def si4(source):
    # Validate each semantic card against its original frozen TikZ text.
    labels = [('Service 1', 'product delivered'), ('Return + pass', 'ready for reuse'),
              ('Service 2', 'product delivered'), ('Return + pass', 'ready for reuse'),
              ('Service 3', 'threshold B = 3')]
    for title, subtitle in labels:
        assert title in source
        assert subtitle in source or subtitle == 'threshold B = 3' and 'threshold $B=3$' in source
    fig = plt.figure(figsize=(7, 1.55))
    ax = fig.add_axes([.01, .05, .98, .90])
    ax.set(xlim=(0, 5), ylim=(0, 1))
    ax.set_axis_off()
    cards = []
    for i, (title, subtitle) in enumerate(labels):
        box = Rectangle((i + .035, .51), .88, .37, linewidth=.6,
                        edgecolor=GREY, facecolor='#F5F6F7', gid='workflow-card')
        ax.add_patch(box)
        texts = [ax.text(i + .475, .752, title, ha='center', va='center', fontsize=7),
                 ax.text(i + .475, .622, subtitle, ha='center', va='center', fontsize=6.5)]
        cards.append((box, texts))
        if i < 4:
            ax.add_patch(FancyArrowPatch((i + .925, .695), (i + 1.025, .695),
                                        arrowstyle='->', mutation_scale=7, linewidth=.65, color=GREY))
    ax.text(2.5, .31, 'L = B - 1 = 2 completed return-to-ready events before service 3',
            ha='center', va='center', fontsize=7)
    ax.text(2.5, .125, 'A post-service-3 return is counted only if the terminal boundary includes it (delta = 1, L = 3).',
            ha='center', va='center', fontsize=6.5)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for box, texts in cards:
        boundary = box.get_window_extent(renderer)
        for text in texts:
            bb = text.get_window_extent(renderer)
            assert boundary.contains(bb.x0 - 2, bb.y0 - 2) and boundary.contains(bb.x1 + 2, bb.y1 + 2)
    return fig, dict(cards=[list(x) for x in labels], connections=[[i, i + 1] for i in range(4)],
                     notes=['L=B-1=2 before service 3', 'delta=1 gives L=3 only with terminal return'],
                     comparison=dict(original_five_events_identical=True, no_probabilities_added=True),
                     node_audit=dict(text_contained_in_each_card=True, arrows_between_disjoint_card_boundaries=True))


def si5(source):
    table = re.search(r'\\label\{tab:scenarios\}(.*?)\\end\{table\}', source, re.S).group(1)
    matches = re.findall(r'^(IFEU[^\n]*|R3PACK[^\n]*) & (0\.\d+) & (0\.\d+) \\\\', table, re.M)
    assert len(matches) == 4
    q = ['0.90', '0.96', '0.97', '0.646']
    labels = ['IFEU individual', 'IFEU open pool', 'IFEU managed pool', 'R3PACK scenario']
    records = []
    for (raw_label, second, third), value, name in zip(matches, q, labels):
        assert value in raw_label
        for exponent, number in [(2, second), (3, third)]:
            assert Decimal(value).__pow__(exponent).quantize(Decimal('.001'), rounding=ROUND_HALF_UP) == Decimal(number)
        records.append(dict(scenario=name, q=value, theta2=float(second), theta3=float(third),
                            printed_theta2=second, printed_theta3=third))
    fig = plt.figure(figsize=(7, 2.5))
    ax = fig.add_axes([.27, .235, .68, .49])
    ax.plot([r['theta2'] for r in records], np.arange(4) - .10, linestyle='none',
            marker='s', markersize=3.5, color=BLUE)
    ax.plot([r['theta3'] for r in records], np.arange(4) + .10, linestyle='none',
            marker='o', markersize=3.5, color=GREY)
    ax.set(ylim=(3.5, -.5), yticks=np.arange(4), yticklabels=labels,
           xlim=(0, 1.02), xticks=[0, .2, .4, .6, .8, 1],
           xlabel='Conditional cumulative attainment (scenario)')
    ax.tick_params(direction='out', length=2.6, width=.5, pad=3)
    key(fig, [plt.Line2D([], [], color=BLUE, marker='s', linestyle='none', label='L = 2 (q²)'),
              plt.Line2D([], [], color=GREY, marker='o', linestyle='none', label='L = 3 (q³)')], x=.60)
    return fig, dict(records=records, comparison=dict(original_table_values_checked=8, maximum_numeric_difference=0),
                     interpretation='Scenario hazards and chain-rule calculation, not measured probabilities',
                     historical_graphic_coordinate_correction=dict(
                         point='R3PACK L=3', old_TikZ_coordinate=2.545,
                         old_axis_implied_value=(2.545 - 1.35) / 4.5, retained_table_value=.270,
                         reason='Replace inaccurate manual graphic coordinate with unchanged frozen table value'))


def setup():
    plt.style.use([str(STYLE / n) for n in STYLE_NAMES])
    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
                         'font.size': 7, 'axes.labelsize': 7, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
                         'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
                         'axes.grid': False, 'axes.spines.top': False, 'axes.spines.right': False,
                         'xtick.top': False, 'ytick.right': False, 'text.usetex': False,
                         'savefig.bbox': None, 'figure.dpi': 120, 'mathtext.fontset': 'dejavusans'})
    from matplotlib.font_manager import FontProperties
    return Path(findfont(FontProperties(family=plt.rcParams['font.sans-serif'])))
