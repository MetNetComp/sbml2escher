"""Shared undirected connectors from central pathways to visible exchanges.

Only real, displayed organic incidences are traversed. This is graph connectivity,
not a claim about feasible directed flux or atom transfer.
"""
import heapq
from collections import defaultdict
from annotation_context import hidden_metabolite

ORANGE = '#f28e2b'


def shortest_paths(adjacency, roots, weight):
    distance = {key: 0 for key in roots}
    parent = {}
    queue = [(0, key) for key in sorted(roots)]
    heapq.heapify(queue)
    while queue:
        cost, key = heapq.heappop(queue)
        if cost != distance[key]:
            continue
        for other in sorted(adjacency[key]):
            proposed = cost + weight(key, other)
            if proposed < distance.get(other, float('inf')):
                distance[other] = proposed
                parent[other] = key
                heapq.heappush(queue, (proposed, other))
    return distance, parent


def shared_coenzymes(model, selected, classification, roots, exchanges, *, allow_disconnected=False):
    """Promote only cofactors required by otherwise disconnected EX endpoints."""
    adjacency = defaultdict(set)
    for reaction in model['reactions']:
        if reaction['id'] not in selected:
            continue
        key = 'rxn:' + reaction['id']
        for mid in reaction['metabolites']:
            if not hidden_metabolite(classification[mid]):
                other = 'met:' + mid
                adjacency[key].add(other); adjacency[other].add(key)
    available = set(adjacency)
    roots = set(roots) & available
    # A coenzyme costs more than any simple organic path. Prefer existing
    # shared organic connections before introducing additional shared cofactors.
    penalty = len(available) + 1
    def weight(a, b):
        mid = a[4:] if a.startswith('met:') else b[4:]
        return 1 + penalty * (classification[mid]['kind'] == 'coenzyme')
    distance, parent = shortest_paths(adjacency, roots, weight)
    promoted = set()
    for terminal in sorted(exchanges):
        if terminal not in distance:
            if allow_disconnected:
                continue
            raise ValueError('EX is disconnected from all displayed major pathways: ' + terminal)
        key = terminal
        while key not in roots:
            if key.startswith('met:') and classification[key[4:]]['kind'] == 'coenzyme':
                promoted.add(key[4:])
            key = parent[key]
    return promoted


def connector_plan(nodes, edges, roots, exchanges, major_reactions, xy=None, *, allow_disconnected=False):
    """Greedily grow a shared shortest-path forest; count each orange edge once."""
    adjacency = defaultdict(set)
    edge_lookup = {}
    for i, edge in enumerate(edges):
        a, b = edge['source'], edge['target']
        if nodes[a]['kind'] == 'currency' or nodes[b]['kind'] == 'currency':
            continue
        adjacency[a].add(b); adjacency[b].add(a)
        edge_lookup[frozenset((a, b))] = i
    roots = set(roots) & set(nodes)
    connected = set(roots)
    pending = set(exchanges) - connected
    chosen = set()
    disconnected = set()
    def weight(a, b):
        edge = edges[edge_lookup[frozenset((a, b))]]
        if edge['reaction_id'] in major_reactions:
            return 0
        return sum(abs(x-y) for x, y in zip(xy[a], xy[b])) if xy else 1
    while pending:
        distance, parent = shortest_paths(adjacency, connected, weight)
        unreachable = pending - set(distance)
        if unreachable:
            if not allow_disconnected:
                raise ValueError('No visible organic path from major pathways to EX: ' + ', '.join(sorted(unreachable)))
            disconnected.update(unreachable)
            pending -= unreachable
            if not pending:
                break
        terminal = min(pending, key=lambda key: (distance[key], key))
        key = terminal
        while key not in connected:
            other = parent[key]
            index = edge_lookup[frozenset((key, other))]
            if edges[index]['reaction_id'] not in major_reactions:
                chosen.add(index)
            connected.add(key)
            key = other
        pending -= connected
    # Independently reconstruct a certificate using only orange and red edges.
    colored = defaultdict(set)
    for i, edge in enumerate(edges):
        if i in chosen or (edge['reaction_id'] in major_reactions and
                           nodes[edge['source']]['kind'] != 'currency' and nodes[edge['target']]['kind'] != 'currency'):
            colored[edge['source']].add(edge['target']); colored[edge['target']].add(edge['source'])
    distance, parent = shortest_paths(colored, roots, lambda a, b: 1)
    paths = {}
    for terminal in sorted(exchanges):
        if terminal in disconnected:
            continue
        path = [terminal]
        while path[-1] not in roots:
            path.append(parent[path[-1]])
        paths[terminal] = path[::-1]
    return dict(color=ORANGE, edge_indices=sorted(chosen), root_node_ids=sorted(roots),
                exchange_paths=paths, all_exchanges_connected=not disconnected,
                disconnected_exchange_node_ids=sorted(disconnected),
                connectivity='Undirected displayed stoichiometric incidence; coenzyme copies excluded; red backbone allowed',
                objective='Sum of Manhattan lengths of unique orange incidence edges; shared edges counted once',
                method='Greedy shared shortest-path forest alternating with constrained grid placement',
                global_optimality_proven=False)


def connector_length(plan, edges, xy):
    return sum(sum(abs(a-b) for a, b in zip(xy[edges[i]['source']], xy[edges[i]['target']]))
               for i in plan['edge_indices'])
