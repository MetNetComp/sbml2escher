"""Reaction incidence graph retaining stoichiometry and per-reaction currency copies."""
import math


def reaction_direction(lower, upper):
    if lower > upper or math.isnan(lower) or math.isnan(upper):
        raise ValueError('Invalid reaction bounds')
    if lower < 0 < upper:
        return 'reversible'
    if upper > 0:
        return 'forward'
    if lower < 0:
        return 'reverse'
    return 'blocked'


def build_network(model, currency):
    currency = set(currency)
    nodes = {}
    edges = []
    for metabolite in model['metabolites']:
        mid = metabolite['id']
        if mid not in currency:
            nodes['met:' + mid] = dict(kind='metabolite', original_id=mid,
                                       name=metabolite.get('name', ''), parent_id='')
    for reaction in model['reactions']:
        rid = reaction['id']
        key = 'rxn:' + rid
        direction = reaction_direction(reaction['lower_bound'], reaction['upper_bound'])
        nodes[key] = dict(kind='reaction', original_id=rid, name=reaction.get('name', ''),
                          parent_id='', direction=direction,
                          lower_bound=reaction['lower_bound'], upper_bound=reaction['upper_bound'])
        for mid, coefficient in sorted(reaction['metabolites'].items()):
            if not math.isfinite(coefficient) or coefficient == 0:
                raise ValueError(f'Invalid stoichiometry: {rid}/{mid}')
            mkey = 'met:' + mid
            if mid in currency:
                mkey = 'copy:' + rid + ':' + mid
                nodes[mkey] = dict(kind='currency', original_id=mid, name=mid, parent_id=key)
            if mkey not in nodes:
                raise ValueError(f'Unknown metabolite: {mid}')
            source, target = (mkey, key) if coefficient < 0 else (key, mkey)
            edges.append(dict(source=source, target=target, reaction_id=rid,
                              metabolite_id=mid, coefficient=coefficient, direction=direction))
    return nodes, edges


def validate_network(model, nodes, edges, xy=None):
    expected = {(r['id'], mid): coefficient for r in model['reactions']
                for mid, coefficient in r['metabolites'].items()}
    actual = {}
    degrees = {}
    for edge in edges:
        pair = (edge['reaction_id'], edge['metabolite_id'])
        if pair in actual:
            raise ValueError('Duplicate reaction incidence')
        actual[pair] = edge['coefficient']
        a, b = nodes[edge['source']], nodes[edge['target']]
        if (a['kind'] == 'reaction') == (b['kind'] == 'reaction'):
            raise ValueError('Graph is not bipartite')
        reaction, metabolite = (a, b) if a['kind'] == 'reaction' else (b, a)
        if reaction['original_id'] != pair[0] or metabolite['original_id'] != pair[1]:
            raise ValueError('Incidence endpoint identity mismatch')
        if (a['kind'] == 'reaction') != (edge['coefficient'] > 0):
            raise ValueError('Incidence sign mismatch')
        if reaction['direction'] != edge['direction']:
            raise ValueError('Direction mismatch')
        for key in (edge['source'], edge['target']):
            degrees[key] = degrees.get(key, 0) + 1
    if expected != actual:
        raise ValueError('Reaction stoichiometry was not preserved')
    expected_ids = {r['id'] for r in model['reactions']}
    if {n['original_id'] for n in nodes.values() if n['kind'] == 'reaction'} != expected_ids:
        raise ValueError('Reaction nodes missing')
    for reaction in model['reactions']:
        node = nodes['rxn:' + reaction['id']]
        if node['direction'] != reaction_direction(reaction['lower_bound'], reaction['upper_bound']):
            raise ValueError('Flux bounds direction was not preserved')
    if any(degrees.get(key) != 1 for key, node in nodes.items() if node['kind'] == 'currency'):
        raise ValueError('Currency copy is not degree one')
    if xy is not None and (set(xy) != set(nodes) or
                           any(not all(math.isfinite(v) for v in p) for p in xy.values())):
        raise ValueError('Every node must have finite coordinates')
    return dict(all_reactions_preserved=True, stoichiometry_preserved=True,
                directions_verified=True, bipartite=True, currency_copies_degree_one=True,
                reactions=len(expected_ids), incidence_edges=len(edges))
