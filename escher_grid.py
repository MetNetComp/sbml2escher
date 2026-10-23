#!/usr/bin/env python3
"""Sampled or structural SBML maps with inorganic species hidden and grid placement."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np

from flux_sampling2 import sample_fluxes
from annotation_context import hidden_metabolite
from grid_layout2 import optimize_grid
from major_backbones2 import RED, backbone_plan, pathway_definitions, validate_backbones
from model_input import load_model
from reaction_network import build_network, validate_network
from stoichiometric_network import stoichiometric_matrix
from exchange_routes import ORANGE, shared_coenzymes, connector_plan, connector_length

BASE = Path(__file__).resolve().parent
# Explicit biochemical display policy; acyl-CoA metabolites remain shared grid nodes.
COENZYMES = {'atp', 'adp', 'amp', 'gtp', 'gdp', 'gmp', 'utp', 'udp', 'ump',
             'ctp', 'cdp', 'cmp', 'nad', 'nadh', 'nadp', 'nadph', 'fad', 'fadh2',
             'fmn', 'fmnh2', 'coa', 'q8', 'q8h2', 'q10', 'q10h2', 'mqn8', 'mql8',
             'thf', 'dhf', '10fthf', '5mthf', 'mlthf', 'pydx5p', 'btn', 'thmpp'}
INORGANIC_IDS = {'h', 'h2o', 'o2', 'co2', 'co', 'hco3', 'co3', 'nh4', 'nh3',
                 'pi', 'ppi', 'so4', 'so3', 'h2s', 'no2', 'no3', 'h2', 'n2',
                 'na1', 'k', 'cl', 'ca2', 'mg2', 'fe2', 'fe3', 'cu2', 'zn2'}
INORGANIC_CARBON_FORMULAS = {'CO', 'CO2', 'CO3', 'CHO3', 'HCO3', 'CH2O3', 'H2CO3'}


def dump(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n')


def model_hash(model):
    return hashlib.sha256(json.dumps(model, sort_keys=True).encode()).hexdigest()


def classify_metabolites(model, override=None):
    overrides = override or {}
    ids = {m['id'] for m in model['metabolites']}
    if set(overrides) - ids or set(overrides.values()) - {'inorganic', 'coenzyme', 'organic'}:
        raise ValueError('Classification overrides must map known metabolite IDs to inorganic/coenzyme/organic')
    result = {}
    for m in model['metabolites']:
        mid = m['id']
        base = re.sub(r'_[a-z][0-9]?$', '', mid).lower()
        formula = m.get('chemicalFormula', m.get('formula', '')) or ''
        elements = set(re.findall(r'[A-Z][a-z]?', formula))
        if mid in overrides:
            kind, reason = overrides[mid], 'explicit override'
        elif base in INORGANIC_IDS or formula in INORGANIC_CARBON_FORMULAS:
            kind, reason = 'inorganic', 'known inorganic ID or carbon oxide/carbonate formula'
        elif elements and 'C' not in elements:
            kind, reason = 'inorganic', 'chemical formula contains no carbon'
        elif base in COENZYMES:
            kind, reason = 'coenzyme', 'explicit coenzyme/nucleotide/electron-carrier ID list'
        else:
            kind, reason = 'organic', 'retained; no inorganic/coenzyme rule matched'
        result[mid] = dict(kind=kind, reason=reason, name=m.get('name', ''), formula=formula,
                           formula_missing=not bool(formula))
    return result


def select_reactions(model, stats, classification, coverage, zero_tolerance=1e-7, major_reactions=()):
    if not math.isfinite(coverage) or not 0 <= coverage <= 100:
        raise ValueError('coverage must be between 0 and 100')
    if not math.isfinite(zero_tolerance) or zero_tolerance <= 0:
        raise ValueError('zero-tolerance must be positive and finite')
    major_reactions = set(major_reactions)
    lookup = {r['id']: r for r in model['reactions']}
    records, eligible = [], []
    if {s['reaction_id'] for s in stats} != set(lookup) or len(stats) != len(lookup):
        raise ValueError('Statistics must contain each source reaction exactly once')
    for stat in stats:
        r = lookup[stat['reaction_id']]
        median = float(stat['median'])
        if not math.isfinite(median):
            raise ValueError('Nonfinite median')
        visible = [mid for mid in r['metabolites'] if not hidden_metabolite(classification[mid])]
        record = dict(stat, priority=abs(median), selected=False, rank=None, median_rank=None,
                      is_major_pathway=r['id'] in major_reactions,
                      visible_metabolites=visible, omitted_inorganic=[m for m in r['metabolites'] if m not in visible])
        record['exclusion_reason'] = ('zero_median' if abs(median) <= zero_tolerance else
                                      'inorganic_only' if not visible else 'below_coverage_cutoff')
        omitted_placeholders = [m for m in r['metabolites'] if classification[m].get('boundary_placeholder')]
        if omitted_placeholders:
            record['omitted_boundary_placeholders'] = omitted_placeholders
            record['omitted_inorganic'] = [m for m in record['omitted_inorganic'] if m not in omitted_placeholders]
        records.append(record)
        if abs(median) > zero_tolerance and visible:
            eligible.append(record)
    eligible.sort(key=lambda r: (-r['priority'], r['reaction_id']))
    for rank, record in enumerate(eligible, 1):
        record['median_rank'] = rank
    eligible.sort(key=lambda r: (not r['is_major_pathway'], -r['priority'], r['reaction_id']))
    count = math.ceil(coverage * len(eligible) / 100)
    for rank, record in enumerate(eligible, 1):
        record['rank'] = rank
        if rank <= count:
            record.update(selected=True, exclusion_reason='')
    selected = {r['reaction_id'] for r in eligible[:count]}
    return selected, records, dict(coverage_percent=coverage, eligible_reactions=len(eligible),
                                  selected_reactions=len(selected), zero_tolerance=zero_tolerance,
                                  ranking='Major pathways first, then descending abs(signed median flux) within each group; reaction ID breaks ties',
                                  major_eligible_reactions=sum(r['is_major_pathway'] for r in eligible),
                                  major_selected_reactions=sum(r['is_major_pathway'] for r in eligible[:count]),
                                  count_rule='ceil(coverage / 100 * eligible_count); 0 gives empty map')


def display_graph(model, selected, classification, required_metabolites=(), shared_metabolites=()):
    reactions = [dict(r, metabolites={m: c for m, c in r['metabolites'].items()
                                     if not hidden_metabolite(classification[m])})
                 for r in model['reactions'] if r['id'] in selected]
    used = {m for r in reactions for m in r['metabolites']}
    used.update(required_metabolites)
    projected = dict(model, metabolites=[m for m in model['metabolites'] if m['id'] in used], reactions=reactions)
    coenzymes = {m for m in used if classification[m]['kind'] == 'coenzyme'} - set(shared_metabolites)
    nodes, edges = build_network(projected, coenzymes)
    return projected, nodes, edges


def select_structural_reactions(model, classification, coverage, major_reactions=()):
    """Select by pathway membership and ID, without inventing flux statistics."""
    if not math.isfinite(coverage) or not 0 <= coverage <= 100:
        raise ValueError('coverage must be between 0 and 100')
    major_reactions = set(major_reactions)
    records, eligible = [], []
    for reaction in model['reactions']:
        visible = [mid for mid in reaction['metabolites'] if not hidden_metabolite(classification[mid])]
        record = dict(reaction_id=reaction['id'], selected=False, rank=None,
                      is_major_pathway=reaction['id'] in major_reactions,
                      visible_metabolites=visible,
                      omitted_inorganic=[mid for mid in reaction['metabolites'] if mid not in visible],
                      lower_bound=reaction['lower_bound'], upper_bound=reaction['upper_bound'],
                      exclusion_reason='below_coverage_cutoff' if visible else 'inorganic_only')
        omitted_placeholders = [m for m in reaction['metabolites'] if classification[m].get('boundary_placeholder')]
        if omitted_placeholders:
            record['omitted_boundary_placeholders'] = omitted_placeholders
            record['omitted_inorganic'] = [m for m in record['omitted_inorganic'] if m not in omitted_placeholders]
        records.append(record)
        if visible:
            eligible.append(record)
    eligible.sort(key=lambda record: (not record['is_major_pathway'], record['reaction_id']))
    count = math.ceil(coverage * len(eligible) / 100)
    for rank, record in enumerate(eligible, 1):
        record['rank'] = rank
        if rank <= count:
            record.update(selected=True, exclusion_reason='')
    selected = {record['reaction_id'] for record in eligible[:count]}
    return selected, records, dict(mode='structure', coverage_percent=coverage,
                                  eligible_reactions=len(eligible), selected_reactions=len(selected),
                                  major_eligible_reactions=sum(r['is_major_pathway'] for r in eligible),
                                  major_selected_reactions=sum(r['is_major_pathway'] for r in eligible[:count]),
                                  ranking='Major pathways first, then reaction ID; no flux or activity ranking',
                                  count_rule='ceil(coverage / 100 * eligible_count); 0 gives empty map')


def attach_coenzymes(nodes, xy, radius=.45, grid_shape=None):
    groups = {}
    for key, node in nodes.items():
        if node['kind'] == 'currency':
            groups.setdefault(node['parent_id'], []).append(key)
    for parent, keys in sorted(groups.items()):
        if grid_shape is None:
            angles = 2 * math.pi * (np.arange(len(keys)) + .25) / len(keys)
        else:
            angles = 2 * math.pi * (np.arange(4096) + .25) / 4096
        vectors = np.column_stack([np.cos(angles), np.sin(angles)])
        vectors *= radius / np.abs(vectors).sum(axis=1, keepdims=True)
        candidates = np.array(xy[parent]) + vectors
        if grid_shape is not None:
            cols, rows = grid_shape
            valid = candidates[(candidates[:, 0] >= 0) & (candidates[:, 0] <= cols - 1) &
                               (candidates[:, 1] >= 0) & (candidates[:, 1] <= rows - 1)]
            if len(valid) < len(keys):
                raise ValueError('Grid cannot contain off-grid coenzyme copies; use at least 2 rows and columns')
            candidates = valid[np.arange(len(keys)) * len(valid) // len(keys)]
        for key, position in zip(sorted(keys), candidates):
            xy[key] = position.tolist()


def render(output, nodes, edges, xy, records, selection, model_id, backbones=None, program_name='sbml2escher', sampling_report=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.cm import ScalarMappable
    stats = {r['reaction_id']: r for r in records}
    structural = selection.get('mode') == 'structure'
    major_mids = {mid for route in (backbones or {}).get('routes', []) for mid in route['metabolite_ids']}
    connectors = (backbones or {}).get('exchange_connectors')
    orange_edges = set(connectors['edge_indices']) if connectors else set()
    orange_nodes = {key for i in orange_edges for key in (edges[i]['source'], edges[i]['target'])}
    selected = [r for r in records if r['selected']]
    priorities = [] if structural else [r['priority'] for r in selected if not r.get('is_major_pathway', False)]
    vmin, vmax = (min(priorities), max(priorities)) if priorities else (1, 10)
    vmax = max(vmax, vmin * 1.01)
    norm = LogNorm(vmin=max(vmin, 1e-12), vmax=vmax)
    cmap = plt.get_cmap('viridis')
    span = np.ptp(np.array(list(xy.values())), axis=0) if xy else np.array([8, 8])
    size = (max(10, min(22, .8 * span[0] + 3)), max(10, min(30, .8 * span[1] + 3)))
    fig, ax = plt.subplots(figsize=size)
    ax.set_aspect('equal')
    if xy:
        arr = np.array(list(xy.values()))
        ax.set_xticks(np.arange(math.floor(arr[:, 0].min()), math.ceil(arr[:, 0].max()) + 1))
        ax.set_yticks(np.arange(math.floor(arr[:, 1].min()), math.ceil(arr[:, 1].max()) + 1))
    ax.grid(color='#d7dee5', alpha=.55, linewidth=.5, zorder=0)
    # Defer autoscaling until all artists exist; scatter otherwise rescans
    # previous collections for every added node on genome-scale maps.
    ax.set_autoscale_on(False)
    for j, edge in enumerate(edges):
        rid = edge['reaction_id']
        a, b = np.array(xy[edge['source']]), np.array(xy[edge['target']])
        reverse = edge['direction'] == 'reverse' if structural else stats[rid]['median'] < 0
        if reverse:
            a, b = b, a
        major = stats[rid].get('is_major_pathway', False)
        orange = j in orange_edges
        color = RED if major else ORANGE if orange else '#cbd3d8' if connectors or structural else cmap(norm(stats[rid]['priority']))
        coenzyme = nodes[edge['source']]['kind'] == 'currency' or nodes[edge['target']]['kind'] == 'currency'
        elbow = np.array([b[0], a[1]])
        path, = ax.plot([a[0], elbow[0], b[0]], [a[1], elbow[1], b[1]],
                        color=color, lw=.7 if coenzyme else 2.3 if major or orange else 1.0,
                        alpha=.95 if major or orange else .35 if connectors else .65 if coenzyme else .55,
                        zorder=2 if major or orange else 1)
        path.set_gid(f'edge-{j}')
        direction = b - elbow if np.linalg.norm(b - elbow) > 1e-10 else b - a
        length = np.linalg.norm(direction)
        if length and (not structural or edge['direction'] != 'blocked'):
            start = b - direction / length * min(.25, length * .6)
            arrow = ax.annotate('', xy=b, xytext=start, arrowprops=dict(arrowstyle='->', color=color, lw=.8), zorder=2)
            arrow.set_gid(f'arrow-{j}')
            if structural:
                arrow.arrow_patch.set_gid(f'arrow-{j}')
        if structural and edge['direction'] == 'reversible':
            initial = elbow - a if np.linalg.norm(elbow - a) > 1e-10 else b - a
            initial_length = np.linalg.norm(initial)
            if initial_length:
                start = a + initial / initial_length * min(.25, initial_length * .6)
                arrow = ax.annotate('', xy=a, xytext=start, arrowprops=dict(arrowstyle='->', color=color, lw=.8), zorder=2)
                arrow.set_gid(f'reverse-arrow-{j}')
                arrow.arrow_patch.set_gid(f'reverse-arrow-{j}')
    for j, (key, node) in enumerate(nodes.items()):
        p = xy[key]; rid = node['original_id'] if node['kind'] == 'reaction' else node['parent_id'].removeprefix('rxn:')
        is_reaction = node['kind'] == 'reaction'; coenzyme = node['kind'] == 'currency'
        major = stats[rid].get('is_major_pathway', False) if is_reaction else node['original_id'] in major_mids
        color = RED if major and not coenzyme else ORANGE if key in orange_nodes else '#83949d' if connectors or structural else cmap(norm(stats[rid]['priority'])) if is_reaction else '#83949d' if coenzyme else '#173c53'
        artist = ax.scatter(*p, marker='s' if is_reaction else 'o', s=25 if is_reaction else 9 if coenzyme else 32,
                            c=[color], zorder=3, edgecolor='white', linewidth=.4)
        artist.set_gid(f'node-{j}')
        label = node['original_id']
        if is_reaction and not structural:
            label += '\n' + f"{stats[rid]['median']:.3g}"
        text = ax.annotate(label, p, xytext=(4, 5) if not coenzyme else (3, -8), textcoords='offset points',
                           fontsize=5 if coenzyme else 6 if is_reaction else 7,
                           color=color, fontweight='normal' if coenzyme else 'bold', zorder=4)
        text.set_gid(f'label-{j}')
    provisional = bool(not structural and sampling_report and not sampling_report['mixing_checks_passed'])
    title_prefix = 'PROVISIONAL: sampling convergence criteria not met\n' if provisional else ''
    ax.set_title(title_prefix + f"{model_id} | {program_name} | top {selection['coverage_percent']:g}%: "
                 f"{len(selected)}/{selection['eligible_reactions']} reactions\n"
                 + ('No flux sampling; arrows show allowed bounds' if structural else
                  'Inorganic species hidden; integer-grid nodes; arrows follow signed median flux'), fontsize=13)
    ax.set_xlabel('Grid x'); ax.set_ylabel('Grid y')
    if priorities and not connectors:
        fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=.6, label='Other reactions: |median flux| (model flux units)')
    if not nodes:
        ax.text(.5, .5, 'No reactions selected', transform=ax.transAxes, ha='center')
    ax.set_autoscale_on(True)
    ax.margins(.08)
    oxidation_label='pyruvate oxidation' if program_name=='sbml2escher' else 'PDH'
    connector_label = (f'Orange: shared undirected EX connectors; unique orange L1 length = {connectors["total_length"]:g}.\n' if connectors else '')
    fig.text(.02, .018, f'Red: glucose exchange / glycolysis / {oxidation_label} / TCA. Gray nodes: local coenzymes.\n' + connector_label +
             ('Orthogonal branches: L1 length. Connectivity and bounds do not establish steady-state activity.' if structural else
              'Orthogonal branches: L1 length. Median describes the sampled feasible space, not a measured flux.'), fontsize=9)
    fig.savefig(output / 'layout.png', dpi=180, bbox_inches='tight')
    fig.savefig(output / 'layout.svg', bbox_inches='tight')
    if connectors:
        visible = set()
        for j, edge in enumerate(edges):
            if j in orange_edges or stats[edge['reaction_id']].get('is_major_pathway',False):
                visible.update((f'edge-{j}', f'arrow-{j}'))
                visible.add(f'reverse-arrow-{j}')
        for j,(key,node) in enumerate(nodes.items()):
            if key in orange_nodes or (node['kind']=='reaction' and stats[node['original_id']].get('is_major_pathway',False)) or (node['kind']=='metabolite' and node['original_id'] in major_mids):
                visible.update((f'node-{j}',f'label-{j}'))
        for artist in ax.get_children():
            if artist.get_gid() and artist.get_gid() not in visible:
                artist.set_visible(False)
        fig.savefig(output/'orange_routes.png',dpi=180,bbox_inches='tight')
        fig.savefig(output/'orange_routes.svg',bbox_inches='tight')
    plt.close(fig)
    # Preserve SVG IDs for filtering without changing exported scientific coordinates.
    svg = (output / 'layout.svg').read_text()
    svg = svg[svg.index('<svg'):]
    registry = []
    for j, edge in enumerate(edges):
        rank = stats[edge['reaction_id']]['rank']
        is_co = nodes[edge['source']]['kind'] == 'currency' or nodes[edge['target']]['kind'] == 'currency'
        highlight = j in orange_edges or stats[edge['reaction_id']].get('is_major_pathway',False)
        registry.extend([dict(id=f'edge-{j}', rank=rank, coenzyme=is_co, highlight=highlight), dict(id=f'arrow-{j}', rank=rank, coenzyme=is_co, highlight=highlight)])
        if structural:
            registry.append(dict(id=f'reverse-arrow-{j}', rank=rank, coenzyme=is_co, highlight=highlight))
    for j, (key, node) in enumerate(nodes.items()):
        ranks = [stats[e['reaction_id']]['rank'] for e in edges if key in (e['source'], e['target'])]
        if node['kind'] == 'reaction':
            ranks = [stats[node['original_id']]['rank']]
        rank = min(ranks) if ranks else 1  # Always-visible glucose reference, except at 0%.
        highlight = key in orange_nodes or (node['kind']=='reaction' and stats[node['original_id']].get('is_major_pathway',False)) or (node['kind']=='metabolite' and node['original_id'] in major_mids)
        registry.extend([dict(id=f'node-{j}', rank=rank, coenzyme=node['kind'] == 'currency',highlight=highlight),
                         dict(id=f'label-{j}', rank=rank, coenzyme=node['kind'] == 'currency',highlight=highlight)])
    page = '''<!doctype html><html lang="ja"><meta charset="utf-8"><title>sbml2escher</title>
<style>body{font:15px system-ui;color:#203744;margin:20px}#viewer{height:78vh;overflow:auto;border:1px solid #ccd;background:white}#viewer svg{width:1600px;height:auto}header{margin-bottom:12px}a{margin-right:14px}</style>
<header><h1>sbml2escher — __MODEL__</h1><p>無機物質を除外。赤：最重要経路（グルコース交換・解糖系・PDH・TCA）。最重要経路を優先し、その中とその他の反応は |flux median| 順。正方形：反応、丸：化合物、灰色小円：補酵素。矢印は中央値の符号方向。</p>
<label>上位 <input id="coverage" type="range" min="0" max="__MAX__" step="1" value="__MAX__"><span id="value"></span>%</label>
<label><input id="coenzymes" type="checkbox" checked>補酵素</label>
<label>拡大 <input type="range" min="600" max="4200" value="1600" step="100" oninput="document.querySelector('#viewer svg').style.width=this.value+'px'"></label>
<p>スライダーは保存済み座標の表示を絞ります。再最適化にはCLIの --coverage を使用してください。</p>
<a href="layout.png">PNG</a><a href="layout.svg">SVG</a><a href="escher_preview.html">Escher</a><a href="report.json">検証結果</a><a href="reaction_ranking.csv">優先順位</a></header>
<div id="viewer">__SVG__</div>
<script>const registry=__REGISTRY__,eligible=__ELIGIBLE__;
function update(){const v=Number(document.getElementById('coverage').value),limit=Math.ceil(v*eligible/100),co=document.getElementById('coenzymes').checked;document.getElementById('value').textContent=v;for(const item of registry){const el=document.getElementById(item.id);if(el)el.style.display=(item.rank<=limit&&(co||!item.coenzyme))?'':'none';}}
document.getElementById('coverage').oninput=update;document.getElementById('coenzymes').onchange=update;update();</script></html>'''
    import html
    if structural:
        page = page.replace('最重要経路を優先し、その中とその他の反応は |flux median| 順。',
                            'フラックスサンプリングなし。最重要経路を優先し、各群内は反応ID順（活性の順位ではありません）。')
        page = page.replace('矢印は中央値の符号方向。', '矢印はモデルの上下限による許可方向（可逆は双方向、上下限0は矢印なし）。')
    if provisional:
        notice=(f'<p class="sampling-warning" style="background:#fff0cc;padding:12px">暫定結果：samplingの収束基準が未達です。優先順位は暫定です。最大split R-hat={sampling_report["max_split_rhat"]:.3f}、最小ESS={sampling_report["min_effective_sample_size"]:.1f}。</p>')
        page=page.replace('</h1>', '</h1>'+notice,1)
    if program_name=='sbml2escher':
        page=page.replace('・PDH・','・ピルビン酸酸化・')
    if connectors:
        page=page.replace('</h1>', '</h1><p>オレンジ：主要経路と各EX反応を結ぶ無向の接続経路。共有枝を1回だけ数えた長さの総和を最小化します。灰色：その他の枝。外周格子点はEX反応専用です。表示を絞ると接続経路の一部が隠れます。</p>',1)
        page=page.replace('<label>拡大','<label><input id="focus" type="checkbox" checked>主要経路・EX接続のみ</label><label>拡大',1)
        page=page.replace('(item.rank<=limit&&(co||!item.coenzyme))', '(item.rank<=limit&&(co||!item.coenzyme)&&(!document.getElementById("focus").checked||item.highlight))')
        page=page.replace('update();</script>', 'document.getElementById("focus").onchange=update;update();</script>')
        if structural:
            page=page.replace('id="focus" type="checkbox" checked', 'id="focus" type="checkbox"')
            missing = len(connectors.get('disconnected_exchange_node_ids', []))
            page=page.replace('</h1>', f'</h1><p>主要経路から接続できないEX反応：{missing}件。これらは灰色で保持し、接続を補いません。</p>',1)
    page = page.replace('sbml2escher', html.escape(program_name)).replace('__MODEL__', html.escape(model_id)).replace('__MAX__', str(selection['coverage_percent'])).replace('__SVG__', svg).replace('__REGISTRY__', json.dumps(registry)).replace('__ELIGIBLE__', str(selection['eligible_reactions']))
    (output / 'index.html').write_text(page)


SAMPLING_FILES = ('flux_samples.npz', 'flux_statistics.json', 'sampling_report.json', 'escher_flux_medians.json')


def validate_sampling_options(args):
    """Protect sampling outputs before any structural-mode files are written."""
    if not getattr(args, 'no_sampling', False):
        return
    if args.reuse_samples:
        raise ValueError('--no-sampling cannot be combined with --reuse-samples')
    if any((args.output_dir / name).exists() for name in SAMPLING_FILES):
        raise ValueError('Use a separate output directory for --no-sampling; sampling files already exist')


def sampling_data(args, model, digest):
    if args.reuse_samples:
        source = args.reuse_samples
        sampling_report = json.loads((source / 'sampling_report.json').read_text())
        if sampling_report.get('model_sha256') != digest:
            raise ValueError('Saved sampling model hash differs from the current model')
        sample_path = source / 'flux_samples.npz'
        if sampling_report.get('samples_sha256') != hashlib.sha256(sample_path.read_bytes()).hexdigest():
            raise ValueError('Saved sample checksum mismatch')
        with np.load(sample_path, allow_pickle=False) as saved:
            samples = saved['fluxes']
            if saved['reaction_ids'].tolist() != [r['id'] for r in model['reactions']]:
                raise ValueError('Saved sample reaction order mismatch')
        if samples.ndim != 3 or samples.shape[2] != len(model['reactions']) or not np.isfinite(samples).all():
            raise ValueError('Invalid saved sample dimensions/values')
        S, _, _ = stoichiometric_matrix(model)
        lo = np.array([r['lower_bound'] for r in model['reactions']]); hi = np.array([r['upper_bound'] for r in model['reactions']])
        if np.abs(samples @ np.asarray(S).T).max() > 1e-7 or max(0, (lo - samples).max(), (samples - hi).max()) > 1e-7:
            raise ValueError('Saved samples violate model constraints')
        stats = json.loads((source / 'flux_statistics.json').read_text())
        medians = np.median(samples.reshape(-1, samples.shape[-1]), axis=0)
        if [s['reaction_id'] for s in stats] != [r['id'] for r in model['reactions']] or not np.allclose([s['median'] for s in stats], medians, rtol=0, atol=1e-12):
            raise ValueError('Saved statistics differ from samples')
    else:
        samples, stats, sampling_report = sample_fluxes(model, args.samples, args.chains, args.thinning,
                                                       args.burn_in, args.seed)
        sampling_report['model_sha256'] = digest
    return samples, stats, sampling_report


def run(args):
    from export_escher import create_map, validate_map, write_preview
    validate_sampling_options(args)
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    model = load_model(args.model)
    strategy = getattr(args, '_strategy', None)
    program_name = getattr(args, 'program_name', 'sbml2escher')
    digest = model_hash(model)
    no_sampling = getattr(args, 'no_sampling', False)
    if no_sampling:
        sampling_report = dict(performed=False, mode='disabled', model_sha256=digest)
    else:
        try:
            reason = getattr(args, '_sampling_unavailable_reason', None)
            if reason:
                raise ValueError(reason)
            samples, stats, sampling_report = sampling_data(args, model, digest)
        except (ValueError, RuntimeError, OSError, ImportError, MemoryError, ArithmeticError, KeyError, TypeError, np.linalg.LinAlgError) as exc:
            # Only sampling acquisition is recoverable; layout/export failures propagate.
            reason = f'{type(exc).__name__}: {exc}'
            print(f'Flux sampling unavailable ({reason}); switching to no sample mode.', flush=True)
            no_sampling = True
            # Preserve existing sampled maps and statistics.
            if any((output / name).exists() for name in SAMPLING_FILES):
                output = output / 'no_sampling'
            structural_args = argparse.Namespace(**vars(args))
            structural_args.output_dir = output
            structural_args.no_sampling = True
            structural_args.reuse_samples = None
            validate_sampling_options(structural_args)
            output.mkdir(parents=True, exist_ok=True)
            sampling_report = dict(performed=False, mode='fallback', reason=reason, model_sha256=digest)
    if not no_sampling:
        if args.reuse_samples is None or args.reuse_samples.resolve() != output.resolve():
            np.savez_compressed(output / 'flux_samples.npz', fluxes=samples,
                                reaction_ids=[r['id'] for r in model['reactions']])
        sampling_report['samples_sha256'] = hashlib.sha256((output / 'flux_samples.npz').read_bytes()).hexdigest()
        dump(output / 'sampling_report.json', sampling_report); dump(output / 'flux_statistics.json', stats)
        if not no_sampling and not sampling_report['mixing_checks_passed']:
            print('Sampling diagnostic thresholds not met; see sampling_report.json. Increase samples/thinning/burn-in.', flush=True)
    override = json.loads(args.classification.read_text()) if args.classification else None
    classification = strategy.classify(model, override) if strategy else classify_metabolites(model, override)
    routes, skipped_routes = strategy.pathways(model) if strategy else pathway_definitions(model)
    major_reactions = {rid for route in routes for rid in route['reaction_ids']}
    if no_sampling:
        selected, records, selection = select_structural_reactions(model, classification, args.coverage, major_reactions)
    else:
        selected, records, selection = select_reactions(model, stats, classification, args.coverage, args.zero_tolerance, major_reactions)
    metabolite_ids = {m['id'] for m in model['metabolites']}
    glucose_id = strategy.glucose_id if strategy else args.glucose_id or next((m for m in ('glc__D_e', 'glc__D_c') if m in metabolite_ids), None)
    if glucose_id is not None and glucose_id not in metabolite_ids:
        raise ValueError('Unknown --glucose-id: ' + glucose_id)
    if glucose_id is not None and classification[glucose_id]['kind'] != 'organic':
        raise ValueError('Glucose must be classified as organic for its fixed grid position')
    references = {glucose_id} if selected and glucose_id else set()
    shared = set()
    roots = {'met:'+mid for route in routes for mid in route['metabolite_ids']}
    roots |= {'rxn:'+rid for rid in major_reactions if rid in selected}
    if strategy and selected:
        visible_exchanges = {'rxn:'+e['reaction_id'] for e in strategy.extraction['exchanges'] if e['reaction_id'] in selected}
        shared = shared_coenzymes(model, selected, classification, roots, visible_exchanges, allow_disconnected=no_sampling)
    projected, nodes, edges = display_graph(model, selected, classification, references, shared)
    grid_ids = sorted(key for key, node in nodes.items() if node['kind'] != 'currency')
    grid_pairs = [(e['source'], e['target']) for e in edges
                  if nodes[e['source']]['kind'] != 'currency' and nodes[e['target']]['kind'] != 'currency']
    glucose_node = 'met:' + glucose_id if references else None
    exchange_nodes = strategy.exchanges(grid_ids) if strategy else {key for key in grid_ids if nodes[key]['kind'] == 'reaction' and nodes[key]['original_id'].startswith('EX_')}
    fixed_offsets, backbones = strategy.plan(routes, grid_ids) if strategy else backbone_plan(routes, grid_ids)
    backbones['skipped_routes'] = skipped_routes
    grid_options = strategy.grid_options(nodes, backbones) if strategy else {}
    connectors = connector_plan(nodes, edges, roots, exchange_nodes, major_reactions, allow_disconnected=no_sampling) if strategy and selected else None
    history = []
    best = None
    previous_xy = None
    for iteration in range(2 if connectors else 1):
        objective_pairs = [(edges[i]['source'], edges[i]['target']) for i in connectors['edge_indices']] if connectors else grid_pairs
        candidate_xy, candidate_placement = optimize_grid(grid_ids, objective_pairs, args.grid_cols, args.grid_rows, args.time_limit, args.seed + iteration,
                                      top_left_node=glucose_node, boundary_nodes=exchange_nodes, fixed_offsets=fixed_offsets,
                                      initial_positions=previous_xy, **grid_options)
        if not connectors:
            xy, placement = candidate_xy, candidate_placement
            break
        candidate_plan = connector_plan(nodes, edges, roots, exchange_nodes, major_reactions, candidate_xy, allow_disconnected=no_sampling)
        before_reroute = connector_length(connectors, edges, candidate_xy)
        if connector_length(candidate_plan, edges, candidate_xy) > before_reroute:
            candidate_plan = connectors
        value = connector_length(candidate_plan, edges, candidate_xy)
        history.append(dict(iteration=iteration+1, before_reroute=before_reroute, orange_length=value,
                            optimized_fixed_network_length=candidate_placement['objective']))
        if best is None or value < best[0]:
            best = (value, candidate_xy, candidate_placement, candidate_plan)
        connectors = candidate_plan
        previous_xy = candidate_xy
    if connectors:
        value, xy, placement, connectors = best
        connectors.update(total_length=value, shared_coenzyme_ids=sorted(shared), optimization_history=history)
        backbones['exchange_connectors'] = connectors
        placement['fixed_network_solver_objective'] = placement['objective']
        placement['objective'] = value
        placement['objective_metric'] = connectors['objective']
        placement['optimality_proven'] = False
        # Fixed-forest MILP bounds cannot certify the joint routing problem.
        placement['fixed_network_lower_bound'] = placement.pop('global_lower_bound', None)
        placement.pop('relative_gap', None)
    cols, rows = placement['grid_cols'], placement['grid_rows']
    attach_coenzymes(nodes, xy, grid_shape=(cols, rows))
    backbone_validation = validate_backbones(backbones, xy, rows)
    backbones['validation'] = backbone_validation
    if glucose_node and xy[glucose_node] != placement['position_constraints']['top_left_position']:
        raise RuntimeError('Glucose is not at its required top-left position')
    if any(xy[key][0] not in (0, cols - 1) and xy[key][1] not in (0, rows - 1) for key in exchange_nodes):
        raise RuntimeError('An EX reaction is not on the grid perimeter')
    if strategy and any(xy[key][0] in (0, cols - 1) or xy[key][1] in (0, rows - 1) for key in set(grid_ids)-exchange_nodes):
        raise RuntimeError('A non-EX grid node occupies the perimeter')
    if any(not (0 <= p[0] <= cols - 1 and 0 <= p[1] <= rows - 1) for p in xy.values()):
        raise RuntimeError('A displayed node lies outside the grid rectangle')
    network = validate_network(projected, nodes, edges, xy)
    for key in grid_ids:
        if not np.allclose(xy[key], np.rint(xy[key]), atol=1e-9, rtol=0):
            raise RuntimeError('A non-coenzyme node is off-grid')
    if len({tuple(xy[k]) for k in grid_ids}) != len(grid_ids):
        raise RuntimeError('Grid nodes overlap')
    if any(hidden_metabolite(classification[n['original_id']]) for n in nodes.values() if n['kind'] != 'reaction'):
        raise RuntimeError('Inorganic species in display')
    length = float(sum(np.abs(np.array(xy[e['source']]) - xy[e['target']]).sum() for e in edges))
    coenzyme_count = sum(n['kind'] == 'currency' for n in nodes.values())
    grid_length = sum(sum(abs(a-b) for a,b in zip(xy[a_key],xy[b_key])) for a_key,b_key in grid_pairs)
    if not math.isclose(length, grid_length + .45 * coenzyme_count, abs_tol=1e-7):
        raise RuntimeError('Edge-length objective differs from drawn orthogonal geometry')
    if connectors and not math.isclose(connector_length(connectors,edges,xy),placement['objective'],abs_tol=1e-7):
        raise RuntimeError('Orange objective differs from drawn orthogonal geometry')
    original = {r['id']: r for r in model['reactions']}
    for reaction in projected['reactions']:
        expected = {m: c for m, c in original[reaction['id']]['metabolites'].items() if not hidden_metabolite(classification[m])}
        if reaction['metabolites'] != expected:
            raise RuntimeError('Visible stoichiometry differs from source')
    metadata = dict(source_model_sha256=digest, coverage=selection, sampling=sampling_report, major_backbones=backbones,
                    original_reactions={r['id']: r for r in model['reactions']}, classification=classification,
                    note='Display projection only: inorganic participants omitted; source_model.json retains the original source model.')
    if nodes:
        data = create_map(projected, nodes, edges, xy, scale=180, padding=80,
                          layout_method=program_name + (' structural grid map' if no_sampling else ' flux-median-ranked grid map'))
        data[0]['metadata'][program_name] = metadata
        data[0]['map_description'] = ('Structural display projection. No sampling or flux statistics; directions from source bounds.' if no_sampling else
                                      'Flux median ranked display projection. Inorganic participants hidden. Sampling uses the full original model.')
        escher_check = validate_map(data, projected, nodes, xy)
        reaction_data = None if no_sampling else {r['reaction_id']: r['median'] * data[0]['metadata']['reactions'][r['reaction_id']]['coefficient_multiplier']
                         for r in records if r['selected']}
        write_preview(output / 'escher_preview.html', data, reaction_data)
        preview = output / 'escher_preview.html'
        preview_html = preview.read_text()
        css = []
        orange_pairs = {frozenset((edges[i]['source'],edges[i]['target'])) for i in connectors['edge_indices']} if connectors else set()
        source_ids = data[0]['metadata']['source_node_ids']
        for rid, reaction in data[1]['reactions'].items():
            if no_sampling and data[0]['metadata']['reactions'][reaction['bigg_id']]['direction'] == 'blocked':
                css.append(f'#r{rid} .arrowheads {{display:none !important;}}')
            if connectors:
                css.append(f'#r{rid} .segment {{stroke:#cbd3d8 !important;}}')
                css.append(f'#r{rid} .arrowheads path {{stroke:#cbd3d8 !important;fill:#cbd3d8 !important;}}')
                for sid, segment in reaction['segments'].items():
                    pair = frozenset((source_ids[segment['from_node_id']],source_ids[segment['to_node_id']]))
                    if pair in orange_pairs:
                        css.extend([f'#r{rid} #s{sid} .segment {{stroke:{ORANGE} !important;}}',
                                    f'#r{rid} #s{sid} .arrowheads path {{stroke:{ORANGE} !important;fill:{ORANGE} !important;}}'])
            if reaction['bigg_id'] not in major_reactions:
                continue
            css.extend([f'#r{rid} .segment {{stroke:{RED} !important;}}',
                        f'#r{rid} .arrowheads path {{stroke:{RED} !important;fill:{RED} !important;}}',
                        f'#r{rid} text {{fill:{RED} !important;}}'])
        style = '<style id="major-pathway-colors">' + '\n'.join(css) + '</style>'
        preview_html = preview_html.replace('<header>', style + '<header>', 1)
        if not no_sampling and not sampling_report['mixing_checks_passed']:
            notice=(f'<p class="sampling-warning" style="background:#fff0cc;padding:12px">暫定結果：samplingの収束基準が未達です。最大split R-hat={sampling_report["max_split_rhat"]:.3f}、最小ESS={sampling_report["min_effective_sample_size"]:.1f}。</p>')
            preview_html=preview_html.replace('<header>','<header>'+notice,1)
        preview_html = preview_html.replace('：FBAフラックスのEscher表示', '：flux中央値のEscher表示')
        if no_sampling:
            preview_html = preview_html.replace('<header>', '<header><p>フラックスサンプリングなし。赤：最重要経路。矢印はモデルの上下限による許可方向です。上下限0の反応は停止表示です。</p>', 1)
        preview_html = preview_html.replace('<p>線の太さ・色は', '<p>赤色は最重要経路（グルコース交換・解糖系・PDH・TCA）です。その他の線の太さ・色は')
        if program_name=='sbml2escher':
            preview_html=preview_html.replace('・PDH・','・ピルビン酸酸化・')
        if connectors:
            preview_html=preview_html.replace('<header>','<header><p>オレンジ：主要経路からEXへの無向接続。灰色：その他。長さの最小化はPNG/SVGの直交枝が対象です。</p>',1)
        preview.write_text(preview_html)
    else:
        data = [dict(schema='https://escher.github.io/escher/jsonschema/1-0-0#', homepage='https://escher.github.io',
                     map_id=digest[:24], map_name='sbml2escher empty', map_description='No reactions selected'),
                dict(nodes={}, reactions={}, text_labels={}, canvas=dict(x=0, y=0, width=500, height=500))]
        from jsonschema import Draft4Validator
        from export_escher import SCHEMA
        Draft4Validator(json.loads(SCHEMA.read_text())).validate(data)
        escher_check = dict(schema_valid=True, reactions=0)
        reaction_data = {}
        (output / 'escher_preview.html').write_text('<!doctype html><meta charset="utf-8"><p>No reactions selected</p>')
    dump(output / 'escher_map.json', data)
    if not no_sampling:
        dump(output / 'escher_flux_medians.json', reaction_data)
    dump(output / 'source_model.json', model); dump(output / 'display_model.json', projected)
    dump(output / 'metabolite_classification.json', classification); dump(output / 'reaction_ranking.json', records)
    dump(output / 'escher_map_metadata.json', metadata)
    dump(output / 'major_backbones.json', backbones)
    with (output / 'reaction_ranking.csv').open('w', newline='') as stream:
        fields = (['reaction_id', 'rank', 'is_major_pathway', 'selected', 'exclusion_reason', 'lower_bound', 'upper_bound'] if no_sampling else
                  ['reaction_id', 'rank', 'is_major_pathway', 'median_rank', 'median', 'priority', 'selected', 'exclusion_reason', 'split_rhat', 'ess', 'q05', 'q95'])
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction='ignore'); writer.writeheader(); writer.writerows(records)
    with (output / 'coordinates.csv').open('w', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(['id', 'kind', 'original_id', 'x', 'y', 'on_grid'])
        for key, node in nodes.items():
            writer.writerow([key, node['kind'], node['original_id'], *xy[key], node['kind'] != 'currency'])
    with (output / 'edges.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['source', 'target', 'reaction_id', 'metabolite_id', 'coefficient', 'direction'])
        writer.writeheader(); writer.writerows(edges)
    report = dict(source_model=str(args.model.resolve()), source_model_sha256=digest,
                  selection=selection, sampling=sampling_report, placement=placement, major_backbones=backbones,
                  source_reactions=len(model['reactions']), drawn_reactions=len(selected),
                  grid_nodes=len(grid_ids), coenzyme_copies=coenzyme_count, edges=len(edges),
                  glucose_metabolite_id=glucose_id, glucose_anchor_active=bool(glucose_node),
                  glucose_reference_only=bool(glucose_node and not any(e['metabolite_id'] == glucose_id for e in edges)),
                  exchange_reaction_ids=sorted(nodes[k]['original_id'] for k in exchange_nodes),
                  total_orthogonal_edge_length=length, coenzyme_edge_length=.45 * coenzyme_count,
                  inorganic_metabolites=[m for m, c in classification.items() if c['kind'] == 'inorganic'],
                  validation=dict(all_noncoenzyme_nodes_on_distinct_grid_points=True, no_inorganic_nodes=True,
                                  glucose_top_left_verified=True, exchange_reactions_on_perimeter=True,
                                  all_nodes_within_grid_rectangle=True, major_backbones=backbone_validation,
                                  no_zero_median_reactions=True, visible_stoichiometry_preserved=True,
                                  full_model_used_for_sampling=True, network=network, escher=escher_check),
                  limitations=['Global length optimality is reported, not assumed.',
                               'Crossings, coincident edge paths, and label overlap are not minimized.',
                               'Display model omits inorganic participants and is not suitable for FBA.',
                               'Escher draws direct segments; orthogonal length objective describes layout.svg and layout.png.',
                               'Sampling uses original bounds without a growth floor or loopless constraints; internal cycles may dominate.'])
    if no_sampling:
        report['validation'].pop('no_zero_median_reactions')
        report['validation'].pop('full_model_used_for_sampling')
        report['validation']['no_sampling_or_flux_statistics_used'] = True
        report['validation']['arrows_follow_source_bounds'] = True
        report['limitations'][-1] = 'Structural connectivity includes inactive and zero-bound reactions; steady-state feasibility is not evaluated.'
    placeholders = [mid for mid, c in classification.items() if c.get('boundary_placeholder')]
    if placeholders:
        report['boundary_placeholder_metabolites'] = placeholders
        report['validation']['no_boundary_placeholder_nodes'] = True
        report['limitations'].append('CellDesigner source/sink placeholders are omitted only in the display; original stoichiometry and sampling constraints are unchanged.')
    report['program_name'] = program_name
    if strategy:
        report['extraction_status'] = strategy.extraction['status']
        report['validation']['glucose_is_top_left_compound'] = True
        report['validation']['perimeter_exclusive_to_exchanges'] = True
        if connectors:
            report['exchange_connectors'] = connectors
            report['orange_orthogonal_edge_length'] = connectors['total_length']
            report['validation']['all_exchanges_connected_to_major_pathways'] = connectors['all_exchanges_connected']
            dump(output/'exchange_connectors.json',connectors)
        strategy.save(output)
    dump(output / 'report.json', report)
    render(output, nodes, edges, xy, records, selection, model['id'], backbones, program_name, sampling_report)
    print(json.dumps(dict(output=str(output.resolve()), selection=selection, sampling_passed=sampling_report.get('mixing_checks_passed'),
                          placement=placement, total_orthogonal_edge_length=length), indent=2), flush=True)
    return report


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model', type=Path, default=BASE / 'e_coli_core.xml')
    p.add_argument('--output-dir', type=Path, default=BASE / 'e_coli_core_sbml2escher_output')
    p.add_argument('--coverage', type=float, default=100, help='Top percentage of eligible reactions (0..100); --no-sampling uses pathway membership and IDs')
    p.add_argument('--zero-tolerance', type=float, default=1e-7)
    p.add_argument('--samples', type=int, default=2000, help='Retained samples per chain')
    p.add_argument('--chains', type=int, default=4)
    p.add_argument('--thinning', type=int, default=100)
    p.add_argument('--burn-in', type=int, default=10000)
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--time-limit', type=float, default=60, help='Grid MILP time limit in seconds')
    p.add_argument('--glucose-id', help='Glucose species to anchor; default glc__D_e, then glc__D_c if present')
    p.add_argument('--grid-cols', type=int)
    p.add_argument('--grid-rows', type=int)
    p.add_argument('--classification', type=Path, help='JSON metabolite ID -> inorganic/coenzyme/organic')
    sampling = p.add_mutually_exclusive_group()
    sampling.add_argument('--reuse-samples', type=Path, help='Reuse a validated sampling output directory; reoptimize layout')
    sampling.add_argument('--no-sampling', action='store_true', help='Draw source structure with no flux sampling/statistics; rank major pathways first, then IDs')
    return p
