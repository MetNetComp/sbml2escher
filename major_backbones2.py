"""Validated E. coli core backbones on integer lines and a square perimeter.

Offsets are measured from the top-left corner: x rightward, y upward (<= 0).
No reaction is added here: the sampling/coverage selection still owns membership.
"""
import json
import math
from pathlib import Path

BASE = Path(__file__).resolve().parent
RED = '#d62728'


def pathway_definitions(model):
    """Enable known routes only when their complete IDs exist in the source model."""
    sources = [('glycolysis', BASE / 'major_pathway.json'),
               ('pdh', BASE / 'pathway_output/pyruvate_to_acetyl_coa.json'),
               ('tca', BASE / 'pathway_output/tca_cycle.json')]
    reactions = {r['id']: r for r in model['reactions']}
    species = {m['id'] for m in model['metabolites']}
    routes, skipped = [], []
    for name, filename in sources:
        data = json.loads(filename.read_text())
        mids = data['metabolite_ids']
        rids = data.get('reaction_ids') or [r['reaction_id'] for r in data['reactions']]
        if not set(mids) <= species or not set(rids) <= set(reactions):
            skipped.append(dict(id=name, reason='Complete canonical pathway IDs are absent from the source model'))
            continue
        pairs = list(zip(mids, mids[1:] + mids[:1])) if name == 'tca' else list(zip(mids, mids[1:]))
        if len(pairs) != len(rids):
            raise ValueError('One reaction is required per backbone metabolite pair: ' + name)
        for rid, (a, b) in zip(rids, pairs):
            reaction = reactions[rid]
            ca = reaction['metabolites'].get(a, 0); cb = reaction['metabolites'].get(b, 0)
            if ca * cb >= 0:
                raise ValueError(f'Invalid canonical backbone incidence: {name}/{rid}/{a}/{b}')
            if reaction['lower_bound'] == reaction['upper_bound'] == 0:
                continue  # A blocked primary reaction is excluded by the median-zero rule.
            if (ca < 0 and reaction['upper_bound'] <= 0) or (ca > 0 and reaction['lower_bound'] >= 0):
                raise ValueError(f'Canonical backbone direction forbidden by model bounds: {name}/{rid}')
        routes.append(dict(id=name, name=data['name'], metabolite_ids=mids, reaction_ids=rids,
                           shape='square' if name == 'tca' else 'vertical' if name == 'glycolysis' else 'horizontal'))
    exchange_id, glucose_id = 'EX_glc__D_e', 'glc__D_e'
    if exchange_id in reactions and glucose_id in species:
        reaction = reactions[exchange_id]
        if set(reaction['metabolites']) != {glucose_id} or reaction['metabolites'][glucose_id] == 0:
            raise ValueError('Glucose exchange must contain only glc__D_e: ' + exchange_id)
        routes.append(dict(id='glucose_exchange', name='Glucose exchange',
                           metabolite_ids=[glucose_id], reaction_ids=[exchange_id], shape='exchange'))
    return routes, skipped


def backbone_plan(routes, displayed_nodes):
    """Fix shortest integer line slots and distribute the TCA nodes around a square."""
    all_offsets, route_offsets = {}, {}
    def store(key, point, local):
        point = list(point)
        if key in all_offsets and all_offsets[key] != point:
            raise ValueError('Conflicting backbone position: ' + key)
        all_offsets[key] = point; local[key] = point
    for route in routes:
        if route['id'] == 'tca':
            continue
        mids, rids = route['metabolite_ids'], route['reaction_ids']
        if route['id'] == 'glucose_exchange':
            local = {}
            store('met:' + mids[0], [0, 0], local)
            store('rxn:' + rids[0], [1, 0], local)
            route_offsets[route['id']] = local
            continue
        start = all_offsets.get('met:' + mids[0], [0, 0] if route['id'] == 'glycolysis' else [0, -2])
        direction = [0, -1] if route['shape'] == 'vertical' else [1, 0]
        chain = []
        for mid, rid in zip(mids, rids):
            chain.extend(['met:' + mid, 'rxn:' + rid])
        chain.append('met:' + mids[-1])
        local = {}
        for step, key in enumerate(chain):
            store(key, [start[0] + step * direction[0], start[1] + step * direction[1]], local)
        route_offsets[route['id']] = local
    square = None
    for route in routes:
        if route['id'] != 'tca':
            continue
        mids, rids = route['metabolite_ids'], route['reaction_ids']
        chain = [key for mid, rid in zip(mids, rids) for key in ('met:' + mid, 'rxn:' + rid)]
        side = math.ceil(len(chain) / 4)
        entry = all_offsets.get('met:accoa_c', [1, -3])
        left, top = entry[0] + 1, entry[1] + 1
        perimeter = ([(left, top - t) for t in range(side)] +
                     [(left + t, top - side) for t in range(side)] +
                     [(left + side, top - side + t) for t in range(side)] +
                     [(left + side - t, top) for t in range(side)])
        local = {}
        slots = [i * len(perimeter) // len(chain) for i in range(len(chain))]
        for key, slot in zip(chain, slots):
            store(key, perimeter[slot], local)
        route_offsets['tca'] = local
        square = dict(left=left, top_offset=top, side_length=side, perimeter_slots=slots,
                      total_perimeter_length=4 * side, corners_occupied=all(c in slots for c in [0, side, 2 * side, 3 * side]))
    fixed = {key: point for key, point in all_offsets.items() if key in displayed_nodes}
    report_routes = [dict(route, displayed_reaction_ids=[r for r in route['reaction_ids'] if 'rxn:' + r in displayed_nodes],
                         displayed_node_ids=[key for key in route_offsets[route['id']] if key in displayed_nodes],
                         node_offsets=route_offsets[route['id']]) for route in routes]
    return fixed, dict(routes=report_routes, square=square, color=RED,
                       fixed_node_offsets=fixed, geometry_scope='Displayed parts only; excluded reactions are not drawn')


def validate_backbones(plan, xy, rows):
    for key, (x, y) in plan['fixed_node_offsets'].items():
        actual = xy[key]
        if abs(actual[0] - x) > 1e-7 or abs(actual[1] - (rows - 1 + y)) > 1e-7:
            raise RuntimeError('Backbone position was not preserved: ' + key)
    checks = []
    for route in plan['routes']:
        points = [xy[key] for key in route['displayed_node_ids']]
        if points and route['shape'] in ('vertical', 'horizontal'):
            axis = 0 if route['shape'] == 'vertical' else 1
            if max(p[axis] for p in points) - min(p[axis] for p in points) > 1e-7:
                raise RuntimeError('Major pathway is not straight: ' + route['id'])
        if points and route['shape'] == 'square':
            square = plan['square']; left = square['left']; top = rows - 1 + square['top_offset']; side = square['side_length']
            if any(not (left <= x <= left + side and top - side <= y <= top and
                        (x in (left, left + side) or y in (top - side, top))) for x, y in points):
                raise RuntimeError('TCA node is not on its square perimeter')
        checks.append(dict(id=route['id'], shape=route['shape'], checked_nodes=len(points), verified=True))
    return dict(all_displayed_backbone_positions_verified=True, routes=checks)
