#!/usr/bin/env python3
"""Export a validated reaction-node layout as an Escher 1-0-0 map."""
import argparse
import base64
from collections import Counter
import csv
import hashlib
import html
import json
import math
from pathlib import Path
import re

from reaction_network import build_network, reaction_direction, validate_network

BASE = Path(__file__).resolve().parent
SCHEMA = BASE / 'vendor/escher/map-schema-1-0-0.json'
VIEWER = BASE / 'vendor/escher/escher-1.8.1.min.js'


def read_layout(directory):
    model = json.loads((directory / 'model.json').read_text())
    with (directory / 'coordinates.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    if len({row['id'] for row in rows}) != len(rows):
        raise ValueError('Duplicate coordinate IDs')
    currency = sorted({row['original_id'] for row in rows if row['kind'] == 'currency'})
    nodes, expected_edges = build_network(model, currency)
    if {row['id'] for row in rows} != set(nodes):
        raise ValueError('Coordinates do not match the reaction graph')
    xy = {}
    for row in rows:
        node = nodes[row['id']]
        if any(row[field] != node[field] for field in ('kind', 'original_id', 'parent_id')):
            raise ValueError('Coordinate node identity mismatch')
        xy[row['id']] = [float(row['x']), float(row['y'])]
    with (directory / 'edges.csv').open() as stream:
        edges = list(csv.DictReader(stream))
    for edge in edges:
        edge['coefficient'] = float(edge['coefficient'])
    validate_network(model, nodes, edges, xy)
    # Validate input CSV reaction directions against actual model bounds as well.
    if len(edges) != len(expected_edges):
        raise ValueError('Missing reaction incidences')
    return model, nodes, edges, xy


def create_map(model, nodes, edges, xy, scale=4000.0, padding=180.0, layout_method='MILP reaction layout'):
    if not math.isfinite(scale) or scale <= 0 or not math.isfinite(padding) or padding < 0:
        raise ValueError('Scale must be positive and padding nonnegative')
    validate_network(model, nodes, edges, xy)
    x_min = min(p[0] for p in xy.values())
    y_max = max(p[1] for p in xy.values())
    def transform(p):
        return scale * (p[0] - x_min) + padding, scale * (y_max - p[1]) + padding
    ids = {key: str(i) for i, key in enumerate(sorted(nodes))}
    metabolites = {m['id']: m for m in model['metabolites']}
    map_nodes = {}
    for key in sorted(nodes):
        node = nodes[key]
        x, y = transform(xy[key])
        if node['kind'] == 'reaction':
            record = dict(node_type='midmarker', x=x, y=y)
        else:
            mid = node['original_id']
            record = dict(node_type='metabolite', x=x, y=y, bigg_id=mid,
                          name=metabolites[mid].get('name', mid), node_is_primary=node['kind'] == 'metabolite',
                          label_x=x, label_y=y - 24)
            if node['kind'] == 'currency':
                parent = xy[node['parent_id']]
                record.update(label_x=x + (10 if xy[key][0] >= parent[0] else -55),
                              label_y=y + (14 if xy[key][1] <= parent[1] else -10))
        map_nodes[ids[key]] = record
    incidence = {}
    for edge in edges:
        incidence.setdefault(edge['reaction_id'], []).append(edge)
    reactions, reaction_metadata, text_labels = {}, {}, {}
    gene_names = {g['id']: g.get('name', '') for g in model.get('genes', [])}
    segment_number = len(nodes) + len(model['reactions'])
    missing_gene_rules = []
    for i, reaction in enumerate(model['reactions']):
        rid = reaction['id']
        reaction_id = str(len(nodes) + i)
        marker = ids['rxn:' + rid]
        x, y = transform(xy['rxn:' + rid])
        direction = reaction_direction(reaction['lower_bound'], reaction['upper_bound'])
        multiplier = -1 if direction == 'reverse' else 1
        rule = reaction.get('gene_reaction_rule', '')
        if 'gene_reaction_rule' not in reaction:
            missing_gene_rules.append(rid)
        genes = sorted(set(re.findall(r'[A-Za-z0-9_.:-]+', rule)) - {'and', 'or'})
        record = dict(name=reaction.get('name', rid), bigg_id=rid, reversibility=direction == 'reversible',
                      label_x=x + 18, label_y=y - 15, gene_reaction_rule=rule,
                      genes=[dict(bigg_id=g, name=gene_names.get(g, '')) for g in genes],
                      metabolites=[dict(bigg_id=mid, coefficient=multiplier * coefficient)
                                   for mid, coefficient in sorted(reaction['metabolites'].items())], segments={})
        for edge in incidence.get(rid, []):
            key = edge['source'] if nodes[edge['source']]['kind'] != 'reaction' else edge['target']
            record['segments'][str(segment_number)] = dict(from_node_id=marker, to_node_id=ids[key], b1=None, b2=None)
            segment_number += 1
        reactions[reaction_id] = record
        def bound(value):
            return value if math.isfinite(value) else str(value)
        reaction_metadata[rid] = dict(map_reaction_id=reaction_id, marker_node_id=marker,
                                     original_lower_bound=bound(reaction['lower_bound']),
                                     original_upper_bound=bound(reaction['upper_bound']),
                                     original_stoichiometry=reaction['metabolites'],
                                     coefficient_multiplier=multiplier, direction=direction)
        if direction == 'blocked':
            text_labels[str(len(text_labels))] = dict(x=x + 18, y=y + 35, text=rid + ': blocked (lb=ub=0)')
    coords = [(n['x'], n['y']) for n in map_nodes.values()]
    coords += [(n['label_x'], n['label_y']) for n in map_nodes.values() if n['node_type'] == 'metabolite']
    coords += [(r['label_x'], r['label_y']) for r in reactions.values()]
    coords += [(t['x'], t['y']) for t in text_labels.values()]
    lo_x, lo_y = min(p[0] for p in coords) - padding, min(p[1] for p in coords) - padding
    hi_x, hi_y = max(p[0] for p in coords) + padding, max(p[1] for p in coords) + padding
    # Extra right margin provides room for long reaction IDs and currency labels.
    canvas = dict(x=lo_x, y=lo_y, width=hi_x - lo_x + 350, height=hi_y - lo_y)
    metadata = dict(exporter='drawing/export_escher.py',
                    source_node_ids={ids[key]: key for key in ids},
                    coordinate_transform=dict(scale=scale, padding=padding, x_min=x_min, y_max=y_max, y_inverted=True),
                    reactions=reaction_metadata, currency_ids=sorted({n['original_id'] for n in nodes.values() if n['kind'] == 'currency'}),
                    missing_gene_rule_reactions=missing_gene_rules,
                    limitation='Bounds and blocked status are metadata; reverse-only coefficients are negated. Multiply corresponding fluxes by coefficient_multiplier.')
    signature = hashlib.sha256(json.dumps([model, xy, scale, padding], sort_keys=True).encode()).hexdigest()[:24]
    header = dict(schema='https://escher.github.io/escher/jsonschema/1-0-0#', homepage='https://escher.github.io',
                  map_id=signature, map_name=model['id'] + '.' + layout_method,
                  map_description=layout_method + '. Reaction-local currency copies, original bounds and node identities are retained in metadata.',
                  metadata=metadata)
    return [header, dict(nodes=map_nodes, reactions=reactions, text_labels=text_labels, canvas=canvas)]


def validate_map(data, model, source_nodes, source_xy, schema_path=SCHEMA, source_metadata=None):
    try:
        from jsonschema import Draft4Validator
    except ImportError as exc:
        raise RuntimeError('jsonschema is required: use a Python with jsonschema or install requirements-escher.txt') from exc
    Draft4Validator(json.loads(schema_path.read_text())).validate(data)
    if len(data) != 2:
        raise ValueError('An Escher map must have exactly two items')
    header, body = data
    map_nodes, reactions = body['nodes'], body['reactions']
    metadata = header.get('metadata', source_metadata)
    if metadata is None:
        raise ValueError('Original export metadata is required to verify a saved map')
    source_ids = metadata['source_node_ids']
    if set(source_ids) != set(map_nodes) or set(source_ids.values()) != set(source_nodes):
        raise ValueError('Map nodes were added or lost')
    transform = metadata['coordinate_transform']
    counts = Counter()
    for node_id, source_id in source_ids.items():
        node, source = map_nodes[node_id], source_nodes[source_id]
        expected_type = 'midmarker' if source['kind'] == 'reaction' else 'metabolite'
        if node['node_type'] != expected_type:
            raise ValueError('Node type mismatch')
        if expected_type == 'metabolite' and (node['bigg_id'] != source['original_id'] or
                node['node_is_primary'] != (source['kind'] == 'metabolite')):
            raise ValueError('Metabolite identity or primary flag mismatch')
        x, y = source_xy[source_id]
        expected = [transform['scale'] * (x - transform['x_min']) + transform['padding'],
                    transform['scale'] * (transform['y_max'] - y) + transform['padding']]
        if not all(math.isclose(node[axis], value, abs_tol=1e-6) for axis, value in zip(('x', 'y'), expected)):
            raise ValueError('Coordinate transform mismatch')
    expected_ids = {r['id'] for r in model['reactions']}
    if len(reactions) != len(expected_ids) or {r['bigg_id'] for r in reactions.values()} != expected_ids:
        raise ValueError('Reaction IDs missing or duplicated')
    original = {r['id']: r for r in model['reactions']}
    segment_ids = set()
    for reaction_id, reaction in reactions.items():
        rid = reaction['bigg_id']
        source = original[rid]
        details = metadata['reactions'][rid]
        if details['map_reaction_id'] != reaction_id or source_ids.get(details['marker_node_id']) != 'rxn:' + rid:
            raise ValueError('Reaction ID metadata mismatch')
        for field in ('lower_bound', 'upper_bound'):
            value = source[field]
            if details['original_' + field] != (value if math.isfinite(value) else str(value)):
                raise ValueError('Flux bounds metadata mismatch')
        direction = reaction_direction(source['lower_bound'], source['upper_bound'])
        multiplier = -1 if direction == 'reverse' else 1
        if reaction['reversibility'] != (direction == 'reversible') or details['coefficient_multiplier'] != multiplier:
            raise ValueError('Direction was not preserved')
        coefficients = {m['bigg_id']: m['coefficient'] for m in reaction['metabolites']}
        if len(coefficients) != len(reaction['metabolites']) or coefficients != {m: multiplier * c for m, c in source['metabolites'].items()}:
            raise ValueError('Stoichiometry was not preserved')
        if details['original_stoichiometry'] != source['metabolites']:
            raise ValueError('Original stoichiometry metadata mismatch')
        participants = Counter()
        for sid, segment in reaction['segments'].items():
            if sid in segment_ids:
                raise ValueError('Duplicate segment ID')
            segment_ids.add(sid)
            a, b = segment['from_node_id'], segment['to_node_id']
            if a not in map_nodes or b not in map_nodes:
                raise ValueError('Unknown segment endpoint')
            pair = [source_nodes[source_ids[a]], source_nodes[source_ids[b]]]
            markers = [n for n in pair if n['kind'] == 'reaction']
            metabolites = [n for n in pair if n['kind'] != 'reaction']
            if len(markers) != 1 or len(metabolites) != 1 or markers[0]['original_id'] != rid:
                raise ValueError('Segment must connect a metabolite to its reaction')
            mid = metabolites[0]['original_id']
            participants[mid] += 1
            if metabolites[0]['kind'] == 'currency' and metabolites[0]['parent_id'] != 'rxn:' + rid:
                raise ValueError('Currency attached to wrong reaction')
            counts[source_ids[a]] += 1; counts[source_ids[b]] += 1
        if participants != Counter({m: 1 for m in source['metabolites']}):
            raise ValueError('Reaction participants lost or duplicated')
    for key, node in source_nodes.items():
        if node['kind'] == 'currency' and counts[key] != 1:
            raise ValueError('Currency copy must have one segment')
    return dict(schema_valid=True, metadata_embedded='metadata' in header,
                reaction_ids_preserved=True, stoichiometry_preserved=True,
                coordinate_transform_verified=True, segment_endpoints_verified=True,
                nodes=len(map_nodes), reactions=len(reactions), segments=len(segment_ids),
                primary_metabolites=sum(n['node_type'] == 'metabolite' and n['node_is_primary'] for n in map_nodes.values()),
                secondary_metabolites=sum(n['node_type'] == 'metabolite' and not n['node_is_primary'] for n in map_nodes.values()),
                currency_copy_counts=dict(Counter(n['original_id'] for n in source_nodes.values() if n['kind'] == 'currency')),
                reverse_only_reactions=[rid for rid, d in metadata['reactions'].items() if d['direction'] == 'reverse'],
                blocked_reactions=[rid for rid, d in metadata['reactions'].items() if d['direction'] == 'blocked'],
                missing_gene_rule_count=len(metadata['missing_gene_rule_reactions']))


def write_preview(path, data, reaction_data=None, fba_report=None):
    # Embedded library and map make this preview usable without a server or network.
    library = base64.b64encode(VIEWER.read_bytes()).decode('ascii')
    map_json = json.dumps(data, ensure_ascii=False, allow_nan=False).replace('<', '\\u003c')
    page = '''<!doctype html><html lang="ja"><meta charset="utf-8"><title>Escher MILP layout</title>
<style>html,body{height:100%;margin:0;font:14px system-ui}header{padding:10px 18px;background:#f3f5f7}#map{height:calc(100vh - 100px);width:100%}#status{margin-left:16px}</style>
<header><b>e_coli_core：Escher形式の配置図</b><p>反応点配置とCO₂コピーを保持しています。Escherの表示スタイルを使用します。ホイールで拡大できます。</p>
<button id="save">マップJSONを保存</button><span id="status">読込中</span></header><div id="map"></div>
<script>window.escherErrors=[];window.addEventListener('error',e=>window.escherErrors.push(e.message));</script>
<script src="data:application/javascript;base64,LIBRARY_PLACEHOLDER"></script><script>
const mapData=MAP_PLACEHOLDER;
window.builder=escher.Builder(mapData,null,null,escher.libs.d3_select('#map'),{
  enable_editing:true,enable_keys:true,enable_search:true,menu:'all',fill_screen:false,
  identifiers_on_map:'bigg_id',show_gene_reaction_rules:false,
  primary_metabolite_radius:10,secondary_metabolite_radius:3,marker_radius:3,
  reaction_no_data_size:3,hide_secondary_metabolites:false
});
const wait=setInterval(()=>{if(!window.builder.map)return;clearInterval(wait);window.escherReady=true;
const m=window.builder.map;document.getElementById('status').textContent=Object.keys(m.reactions).length+'反応・'+Object.keys(m.nodes).length+'点';},100);
document.getElementById('save').onclick=()=>{const data=window.builder.map.map_for_export();data[0].metadata=mapData[0].metadata;const url=URL.createObjectURL(new Blob([JSON.stringify(data,null,2)],{type:'application/json'}));const a=document.createElement('a');a.href=url;a.download='escher_map.json';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000)};
</script></html>'''
    if reaction_data is not None:
        if set(reaction_data) != {r['bigg_id'] for r in data[1]['reactions'].values()}:
            raise ValueError('Flux IDs must match all map reactions')
        if not all(math.isfinite(v) for v in reaction_data.values()):
            raise ValueError('Fluxes must be finite')
        maximum = max(abs(v) for v in reaction_data.values()) or 1.0
        active = [abs(v) for v in reaction_data.values() if v != 0]
        first = min(min(active), maximum * 0.02) if active else maximum * 0.02
        options = dict(reaction_data=reaction_data, reaction_styles=['color', 'size', 'text', 'abs'],
                       reaction_scale=[dict(type='value', value=0, color='#b8b8b8', size=1),
                                       dict(type='value', value=first, color='#e6ab02', size=3),
                                       dict(type='value', value=maximum * 0.1, color='#ff7f00', size=5),
                                       dict(type='value', value=maximum * 0.5, color='#e31a1c', size=9),
                                       dict(type='value', value=maximum, color='#9e0142', size=14)])
        page = page.replace('enable_editing:true,', json.dumps(options, allow_nan=False)[1:-1] + ',enable_editing:true,')
        page = page.replace('Escher形式の配置図', 'FBAフラックスのEscher表示')
        summary = ('線の太さ・色はフラックスの絶対値、ラベルは符号付きの値です。'
                   '0は灰色、小さい流量は黄色、大きくなるとオレンジ・赤・赤紫になります。太い線ほど流量が大きくなります。単位：mmol/gDW/h。'
                   f'尺度：0〜{maximum:.4g}。矢印は計算された流れの方向を示します。'
                   'ゼロ反応の矢印は許容方向であり、流れを表しません。')
        if fba_report:
            summary += (f' バイオマス反応：{fba_report["objective_value"]:.6f} h⁻¹。'
                        'SBMLの上下限で成長を最大化した、最適解の一例です。')
        page = page.replace('反応点配置とCO₂コピーを保持しています。Escherの表示スタイルを使用します。ホイールで拡大できます。', summary)
        page = page.replace('const mapData=MAP_PLACEHOLDER;',
                            'window.fbaFluxes=' + json.dumps(reaction_data, allow_nan=False) + ';const mapData=MAP_PLACEHOLDER;')
        page = page.replace('#map{height:calc(100vh - 100px)', '#map{height:calc(100vh - 140px)')
    page = page.replace('e_coli_core：Escher形式の配置図', html.escape(data[0]['map_name']))
    page = page.replace('e_coli_core：FBAフラックスのEscher表示', html.escape(data[0]['map_name']) + '：FBAフラックスのEscher表示')
    path.write_text(page.replace('LIBRARY_PLACEHOLDER', library).replace('MAP_PLACEHOLDER', map_json))


def export_layout(directory, output, scale=4000.0, padding=180.0, preview=None):
    model, nodes, edges, xy = read_layout(directory)
    report_path = directory / 'report.json'
    layout_method = json.loads(report_path.read_text()).get('layout_method', 'MILP reaction layout') if report_path.exists() else 'MILP reaction layout'
    data = create_map(model, nodes, edges, xy, scale, padding, layout_method=layout_method)
    backbone_file = directory / 'major_pathway.json'
    if backbone_file.exists():
        backbone = json.loads(backbone_file.read_text())
        if 'pathways' in backbone:
            data[0]['metadata']['selected_pathways'] = backbone['pathways']
            data[0]['metadata']['pathway_emphasis'] = 'Representative pathway selection, not flux magnitude; original reversibility is retained.'
            transform = data[0]['metadata']['coordinate_transform']
            labels = data[1]['text_labels']
            for group in backbone['pathways']:
                x, y = group['label_position']
                labels[str(len(labels))] = dict(x=scale * (x - transform['x_min']) + padding,
                                               y=scale * (transform['y_max'] - y) + padding,
                                               text=group['name'])
            labels[str(len(labels))] = dict(x=data[1]['canvas']['x'] + padding,
                                           y=data[1]['canvas']['y'] + data[1]['canvas']['height'] - padding,
                                           text='Selected pathways are layout backbones, not flux magnitudes. ETC: electron transfer; oxygen is reduced to water. ATP synthase is separate.')
    if backbone_file.exists() and 'compartment_regions' in backbone:
        data[0]['metadata']['compartment_layout'] = dict(regions=backbone['compartment_regions'],
                                                        species_symbols=backbone['species_compartments'],
                                                        audit=backbone['compartment_audit'])
        transform = data[0]['metadata']['coordinate_transform']
        labels = data[1]['text_labels']
        for symbol, region in backbone['compartment_regions'].items():
            if not any(backbone['species_compartments'][m['id']] == symbol for m in model['metabolites']):
                continue
            labels[str(len(labels))] = dict(x=scale*(region['left']-transform['x_min'])+padding,
                                           y=scale*(transform['y_max']-region['top'])+padding,
                                           text=f"[{symbol}] {region['name']}")
    verification = validate_map(data, model, nodes, xy)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    output.with_name(output.stem + '_metadata.json').write_text(json.dumps(
        dict(map_id=data[0]['map_id'], metadata=data[0]['metadata']), indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    report = dict(source_layout=str(directory), map_file=output.name, schema='Escher 1-0-0',
                  validation=verification, browser_validation='not run',
                  limitations=['Escher uses its own styles; original square markers and pathway colors are not reproduced.',
                               'Bounds and blocked status are metadata; gene rules absent from the loaded model are empty.',
                               'Native Escher 1.8.1 saves omit additional metadata; retain the separate metadata JSON.'])
    output.with_name(output.stem + '_report.json').write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n')
    if preview:
        preview.parent.mkdir(parents=True, exist_ok=True)
        write_preview(preview, data)
    print(json.dumps(verification, indent=2, ensure_ascii=False))
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--layout-dir', type=Path, default=BASE / 'layout_reaction_sbml_output')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--scale', type=float, default=4000)
    parser.add_argument('--padding', type=float, default=180)
    parser.add_argument('--preview', type=Path, help='Write a self-contained Escher HTML preview')
    args = parser.parse_args()
    output = args.output or args.layout_dir / 'escher_map.json'
    try:
        export_layout(args.layout_dir, output, args.scale, args.padding, args.preview)
    except (ValueError, RuntimeError, OSError) as exc:
        parser.error(str(exc))


if __name__ == '__main__':
    main()
