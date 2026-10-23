"""Distinct integer-grid placement minimizing unweighted incidence L1 length."""
import math
import numpy as np
from scipy.optimize import Bounds, LinearConstraint, linear_sum_assignment, milp
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee
from scipy.spatial import cKDTree


def edge_length(points, edges):
    return int(np.abs(points[edges[:, 0]] - points[edges[:, 1]]).sum()) if len(edges) else 0


def large_grid_initialization(n, edges, cells, allowed):
    """Seed large grids without dense eigendecomposition or assignment costs.

    Grid constraints form disjoint boundary/interior sets with nested margins
    and singleton anchors. Assign the most restricted sets first, then choose
    nearby unoccupied cells using a spatial index. This is a feasible heuristic,
    not a minimum-cost assignment.
    """
    if len(edges):
        a, b = edges.T
        adjacency = coo_matrix((np.ones(2 * len(edges)),
                                (np.r_[a, b], np.r_[b, a])), shape=(n, n)).tocsr()
        order = reverse_cuthill_mckee(adjacency, symmetric_mode=True)
    else:
        order = np.arange(n)
    # A serpentine grid traversal keeps successive graph-order nodes nearby.
    traversal = np.lexsort((np.where(cells[:, 1] % 2, -cells[:, 0], cells[:, 0]), cells[:, 1]))
    targets = np.empty((n, 2))
    targets[order] = cells[traversal[np.linspace(0, len(cells) - 1, n).astype(int)]]
    groups = {}
    for i, mask in enumerate(allowed):
        key = np.packbits(mask).tobytes()
        if key not in groups:
            groups[key] = [np.flatnonzero(mask), []]
        groups[key][1].append(i)
    occupied = np.zeros(len(cells), dtype=bool)
    points = np.empty((n, 2), dtype=int)
    for candidates, indices in sorted(groups.values(), key=lambda group: len(group[0])):
        available = candidates[~occupied[candidates]]
        if len(available) < len(indices):
            raise ValueError('Grid cannot satisfy position constraints without overlapping nodes')
        tree = cKDTree(cells[available])
        nearest = min(8, len(available))
        _, proposals = tree.query(targets[indices], k=nearest)
        proposals = np.asarray(proposals).reshape(len(indices), nearest)
        for i, proposed in zip(indices, proposals):
            choices = available[proposed]
            free = choices[~occupied[choices]]
            k = nearest
            while not len(free):
                k = min(len(available), k * 2)
                _, proposed = tree.query(targets[i], k=k)
                choices = available[np.atleast_1d(proposed)]
                free = choices[~occupied[choices]]
                if k == len(available) and not len(free):
                    raise ValueError('Grid cannot satisfy position constraints without overlapping nodes')
            cell = free[0]
            occupied[cell] = True
            points[i] = cells[cell]
    return points


def heuristic_layout(n, edges, cols, rows, seed=0, allowed=None, initial_points=None):
    rng = np.random.default_rng(seed)
    cells = np.array([(x, y) for y in range(rows) for x in range(cols)])
    if n == 0:
        return np.empty((0, 2), dtype=int), dict(initial_length=0, final_length=0)
    if allowed is None:
        allowed = np.ones((n, len(cells)), dtype=bool)
    if initial_points is None and n > 5000:
        points = large_grid_initialization(n, edges, cells, allowed)
    elif initial_points is None:
        adjacency = np.zeros((n, n))
        for a, b in edges:
            adjacency[a, b] += 1; adjacency[b, a] += 1
        laplacian = np.diag(adjacency.sum(axis=1)) - adjacency
        _, vectors = np.linalg.eigh(laplacian)
        embedding = np.zeros((n, 2))
        for axis in range(min(2, max(0, n - 1))):
            v = vectors[:, axis + 1]
            embedding[:, axis] = (v - v.min()) / max(np.ptp(v), 1e-12) * ([cols, rows][axis] - 1)
        costs = np.abs(embedding[:, None] - cells[None]).sum(axis=2)
        try:
            _, assignment = linear_sum_assignment(np.where(allowed, costs, np.inf))
        except ValueError as exc:
            raise ValueError('Grid cannot satisfy top-left/boundary constraints without overlapping nodes') from exc
        points = cells[assignment].copy()
    else:
        points = np.array(initial_points, dtype=int).copy()
    initial = edge_length(points, edges)
    neighbors = [[] for _ in range(n)]
    for a, b in edges:
        neighbors[a].append(b); neighbors[b].append(a)
    incident = [np.flatnonzero((edges == i).any(axis=1)) for i in range(n)] if len(edges) else [np.array([], dtype=int) for _ in range(n)]
    owner = {tuple(p): i for i, p in enumerate(points)}
    best, best_value = points.copy(), initial
    current_value = initial
    for iteration in range(max(20000, n * 600)):
        i = int(rng.integers(n))
        p = cells[int(rng.integers(len(cells)))]
        j = owner.get(tuple(p))
        if j == i or not allowed[i, int(p[1]) * cols + int(p[0])]:
            continue
        if j is not None and not allowed[j, int(points[i, 1]) * cols + int(points[i, 0])]:
            continue
        affected = incident[i] if j is None else np.union1d(incident[i], incident[j])
        before = edge_length(points, edges[affected])
        old = points[i].copy()
        points[i] = p
        if j is not None:
            points[j] = old
        delta = edge_length(points, edges[affected]) - before
        temperature = 1.8 * (1 - iteration / max(20000, n * 600)) ** 2 + .02
        if delta <= 0 or rng.random() < math.exp(-delta / temperature):
            owner.pop(tuple(old)); owner[tuple(p)] = i
            if j is not None:
                owner[tuple(old)] = j
            current_value += delta
            if current_value < best_value:
                best_value, best = current_value, points.copy()
        else:
            points[i] = old
            if j is not None:
                points[j] = p
    points = best
    # Exhaustive swaps and moves to empty cells until a local minimum.
    while True:
        owner = {tuple(p): i for i, p in enumerate(points)}
        improvement = 0
        for i in rng.permutation(n):
            if n > 300:
                # Genome-scale maps: polish near the L1 median of incident nodes,
                # their current point and connected nodes rather than every cell.
                adjacent = neighbors[i]
                target = np.median(points[adjacent], axis=0) if adjacent else points[i]
                centers = [np.rint(target).astype(int), points[i]] + [points[j] for j in adjacent]
                candidates = set()
                for center in centers:
                    for dx in range(-2, 3):
                        for dy in range(-2, 3):
                            x, y = int(center[0]) + dx, int(center[1]) + dy
                            if 0 <= x < cols and 0 <= y < rows:
                                candidates.add(y * cols + x)
                candidate_points = cells[sorted(candidates)]
            else:
                candidate_points = cells
            for p in candidate_points:
                j = owner.get(tuple(p))
                if j == i or not allowed[i, int(p[1]) * cols + int(p[0])]:
                    continue
                if j is not None and not allowed[j, int(points[i, 1]) * cols + int(points[i, 0])]:
                    continue
                affected = incident[i] if j is None else np.union1d(incident[i], incident[j])
                before = edge_length(points, edges[affected])
                old = points[i].copy(); points[i] = p
                if j is not None:
                    points[j] = old
                delta = edge_length(points, edges[affected]) - before
                if delta < 0:
                    owner.pop(tuple(old)); owner[tuple(p)] = i
                    if j is not None:
                        owner[tuple(old)] = j
                    improvement += -delta
                else:
                    points[i] = old
                    if j is not None:
                        points[j] = p
        if not improvement:
            break
    return points, dict(method=(('Graph-order spatial assignment, seeded annealing, sparse-neighborhood swaps/empty moves' if n > 5000 else
                                 'Spectral assignment, seeded annealing, sparse-neighborhood swaps/empty moves') if n > 300 else
                                'Spectral assignment, seeded annealing, exhaustive swaps/empty moves'),
                        polishing_scope='Local neighborhoods of medians/current/adjacent points' if n > 300 else 'All grid cells',
                        initial_length=initial, final_length=edge_length(points, edges))


def optimize_grid(node_ids, edge_pairs, cols=None, rows=None, time_limit=60, seed=0, progress=print,
                  *, top_left_node=None, boundary_nodes=(), fixed_offsets=None, top_left_margin=0, node_top_margins=None, max_milp_nodes=None,
                  exclusive_boundary=False, initial_positions=None):
    auto_cols, auto_rows = cols is None, rows is None
    n = len(node_ids)
    if len(set(node_ids)) != n:
        raise ValueError('Duplicate grid node IDs')
    boundary_nodes = set(boundary_nodes)
    node_top_margins = dict(node_top_margins or {})
    if top_left_margin < 0 or int(top_left_margin) != top_left_margin:
        raise ValueError('Top-left margin must be a nonnegative integer')
    if set(node_top_margins) - set(node_ids) or any(v < 0 or int(v) != v for v in node_top_margins.values()):
        raise ValueError('Node top margins must map known nodes to nonnegative integers')
    fixed_offsets = dict(fixed_offsets or {})
    if set(fixed_offsets) - set(node_ids):
        raise ValueError('Unknown fixed backbone node')
    for key, point in fixed_offsets.items():
        if len(point) != 2 or any(not math.isfinite(v) or int(v) != v for v in point) or point[0] < 0 or point[1] > 0:
            raise ValueError('Backbone offsets must be integer [x>=0, y<=0] values')
    if len(set(map(tuple, fixed_offsets.values()))) != len(fixed_offsets):
        raise ValueError('Backbone nodes share a fixed position')
    anchor_x = int(exclusive_boundary)
    if top_left_node in fixed_offsets and fixed_offsets[top_left_node] != [anchor_x, -top_left_margin]:
        raise ValueError('Backbone position conflicts with the top-left anchor')
    if boundary_nodes - set(node_ids) or (top_left_node is not None and top_left_node not in node_ids):
        raise ValueError('Unknown top-left or boundary node')
    index = {key: i for i, key in enumerate(node_ids)}
    edges = np.array([(index[a], index[b]) for a, b in edge_pairs], dtype=int).reshape(-1, 2)
    if cols is None and rows is None:
        cols = max(2, math.ceil(math.sqrt(max(n, 1) * 1.15)))
        rows = max(2, math.ceil(max(n, 1) * 1.15 / cols))
    elif cols is None:
        cols = max(2, math.ceil(max(n, 1) / rows))
    elif rows is None:
        rows = max(2, math.ceil(max(n, 1) / cols))
    if fixed_offsets:
        # Leave a one-cell margin below and to the right of the backbone.
        min_cols = max(x for x, y in fixed_offsets.values()) + 2
        min_rows = -min(y for x, y in fixed_offsets.values()) + 2
        if auto_cols:
            cols = max(cols, min_cols)
        elif cols < min_cols - 1:
            raise ValueError('Grid has too few columns for the straight/square backbones')
        if auto_rows:
            rows = max(rows, min_rows)
        elif rows < min_rows - 1:
            raise ValueError('Grid has too few rows for the straight/square backbones')
    if rows <= max([top_left_margin] + list(node_top_margins.values())):
        raise ValueError('Grid has too few rows for the exchange/compound top margin')
    if cols < 1 or rows < 1 or cols * rows < n:
        raise ValueError('Grid must contain at least one distinct point per non-coenzyme node')
    if n == 0:
        return {}, dict(grid_cols=cols, grid_rows=rows, objective=0, optimality_proven=True, method='Empty graph')
    def occupied_boundary_slots():
        fixed = {key: [x, rows - 1 + y] for key, (x, y) in fixed_offsets.items()}
        if top_left_node is not None:
            fixed[top_left_node] = [anchor_x, rows - 1 - top_left_margin]
        return sum(key not in boundary_nodes and (x in (0, cols - 1) or y in (0, rows - 1))
                   for key, (x, y) in fixed.items())
    def perimeter():
        return cols * rows - max(cols - 2, 0) * max(rows - 2, 0)
    while perimeter() < len(boundary_nodes) + occupied_boundary_slots():
        if not auto_cols and not auto_rows:
            raise ValueError('Grid perimeter has too few points for EX reactions and the top-left anchor')
        if auto_cols and (not auto_rows or cols <= rows):
            cols += 1
        else:
            rows += 1
    if exclusive_boundary:
        if n > len(boundary_nodes) and ((not auto_cols and cols < 3) or (not auto_rows and rows < 3)):
            raise ValueError('Grid interior requires at least 3 columns and rows')
        while max(cols - 2, 0) * max(rows - 2, 0) < n - len(boundary_nodes):
            if not auto_cols and not auto_rows:
                raise ValueError('Grid interior has too few points; perimeter is EX-only')
            if auto_cols and (not auto_rows or cols <= rows):
                cols += 1
            else:
                rows += 1
    cells = np.array([(x, y) for y in range(rows) for x in range(cols)])
    boundary_mask = ((cells[:, 0] == 0) | (cells[:, 0] == cols - 1) |
                     (cells[:, 1] == 0) | (cells[:, 1] == rows - 1))
    allowed = np.ones((n, len(cells)), dtype=bool)
    for key in boundary_nodes:
        allowed[index[key]] = boundary_mask
    if exclusive_boundary:
        for key in set(node_ids) - boundary_nodes:
            allowed[index[key]] = ~boundary_mask
    if top_left_node is not None:
        allowed[index[top_left_node]] &= (cells[:, 0] == anchor_x) & (cells[:, 1] == rows - 1 - top_left_margin)
    for key, margin in node_top_margins.items():
        allowed[index[key]] &= cells[:, 1] <= rows - 1 - margin
    for key, (x, y) in fixed_offsets.items():
        fixed_mask = (cells[:, 0] == x) & (cells[:, 1] == rows - 1 + y)
        allowed[index[key]] &= fixed_mask
    def validate_placement(candidate):
        if (len(set(map(tuple, candidate))) != n or np.any(candidate < 0) or
                np.any(candidate[:, 0] >= cols) or np.any(candidate[:, 1] >= rows)):
            raise RuntimeError('Invalid grid placement')
        if not all(allowed[i, int(y) * cols + int(x)] for i, (x, y) in enumerate(candidate)):
            raise RuntimeError('Grid placement violates top-left or EX boundary constraints')
    def constraint_report():
        return dict(top_left_node=top_left_node,
                    top_left_position=[anchor_x, rows - 1 - top_left_margin] if top_left_node is not None else None,
                    top_left_margin=top_left_margin, node_top_margins=node_top_margins,
                    boundary_node_ids=sorted(boundary_nodes),
                    exclusive_boundary=exclusive_boundary,
                    top_left_verified=True, boundary_nodes_verified=True,
                    fixed_backbone_positions={k: [x, rows - 1 + y] for k, (x, y) in fixed_offsets.items()},
                    fixed_backbone_positions_verified=True)
    initial_points = None
    if initial_positions is not None:
        initial_points = np.array([initial_positions[key] for key in node_ids], dtype=float).reshape(-1,2)
        if not np.isfinite(initial_points).all() or not np.equal(initial_points,np.rint(initial_points)).all():
            raise ValueError('Initial grid coordinates must be finite integers')
        validate_placement(initial_points)
    progress(f'Grid: initializing {n} nodes on {cols} x {rows} cells', flush=True)
    points, heuristic = heuristic_layout(n, edges, cols, rows, seed, allowed, initial_points)
    heuristic['initialization'] = ('Validated previous placement' if initial_points is not None else
                                   'Graph-order spatial assignment' if n > 5000 else 'Spectral assignment')
    validate_placement(points)
    incumbent = edge_length(points, edges)
    progress(f'Grid: {n} nodes, {len(edges)} edges, {cols} x {rows}, initial length {incumbent}', flush=True)
    if not len(edges):
        return dict(zip(node_ids, points.tolist())), dict(grid_cols=cols, grid_rows=rows, objective=0,
                                                        optimality_proven=True, method='Edgeless graph',
                                                        position_constraints=constraint_report())
    if max_milp_nodes is not None and n > max_milp_nodes:
        bound=float(len(edges))
        return dict(zip(node_ids, points.tolist())), dict(grid_cols=cols, grid_rows=rows,
                    objective=incumbent, objective_metric='Unweighted Manhattan incidence length in grid units',
                    method='Constraint-preserving heuristic; global MILP exceeds configured node budget',
                    selected_solution='validated heuristic incumbent', heuristic=heuristic,
                    optimality_proven=incumbent==bound, global_lower_bound=bound,
                    relative_gap=(incumbent-bound)/max(incumbent,1), max_milp_nodes=max_milp_nodes,
                    position_constraints=constraint_report(), distinct_grid_points_verified=True)
    # x_i, y_i integer; d_xe, d_ye continuous; four directional binaries per pair.
    e = len(edges)
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    boundary_offset = 2 * n + 2 * e + 4 * len(pairs)
    size = boundary_offset + 4 * len(boundary_nodes)
    cost = np.zeros(size); cost[2 * n:2 * n + 2 * e] = 1
    lb, ub = np.zeros(size), np.ones(size)
    ub[:n] = cols - 1; ub[n:2 * n] = rows - 1
    ub[2 * n:2 * n + e] = cols - 1; ub[2 * n + e:2 * n + 2 * e] = rows - 1
    if top_left_node is not None:
        anchor = index[top_left_node]
        lb[anchor] = ub[anchor] = anchor_x
        lb[n + anchor] = ub[n + anchor] = rows - 1 - top_left_margin
    for key, (x, y) in fixed_offsets.items():
        anchor = index[key]
        lb[anchor] = ub[anchor] = x
        lb[n + anchor] = ub[n + anchor] = rows - 1 + y
    for key, margin in node_top_margins.items():
        ub[n + index[key]] = min(ub[n + index[key]], rows - 1 - margin)
    if exclusive_boundary:
        for key in set(node_ids) - boundary_nodes:
            i = index[key]
            lb[i] = max(lb[i], 1); ub[i] = min(ub[i], cols - 2)
            lb[n + i] = max(lb[n + i], 1); ub[n + i] = min(ub[n + i], rows - 2)
    integrality = np.ones(size); integrality[2 * n:2 * n + 2 * e] = 0
    ii, jj, vv, lows, highs = [], [], [], [], []
    def add(terms, low=-np.inf, high=np.inf):
        k = len(lows)
        for j, v in terms:
            ii.append(k); jj.append(j); vv.append(v)
        lows.append(low); highs.append(high)
    for k, (a, b) in enumerate(edges):
        for offset, d in [(0, 2 * n + k), (n, 2 * n + e + k)]:
            add([(a + offset, 1), (b + offset, -1), (d, -1)], high=0)
            add([(a + offset, -1), (b + offset, 1), (d, -1)], high=0)
    for k, (a, b) in enumerate(pairs):
        z = 2 * n + 2 * e + 4 * k
        for axis, big_m in [(0, cols), (1, rows)]:
            offset = axis * n
            add([(a + offset, 1), (b + offset, -1), (z + axis * 2, big_m)], high=big_m - 1)
            add([(b + offset, 1), (a + offset, -1), (z + axis * 2 + 1, big_m)], high=big_m - 1)
        add([(z + j, 1) for j in range(4)], low=1)
    add([(j, 1) for j in range(2 * n, 2 * n + 2 * e)], high=incumbent)
    for k, key in enumerate(sorted(boundary_nodes)):
        i = index[key]; z = boundary_offset + 4 * k
        add([(i, 1), (z, cols)], high=cols)
        add([(i, -1), (z + 1, cols)], high=1)
        add([(n + i, 1), (z + 2, rows)], high=rows)
        add([(n + i, -1), (z + 3, rows)], high=1)
        add([(z + j, 1) for j in range(4)], low=1)
    # A fixed corner removes reflection symmetry; reflecting would move the anchor.
    if n > 1 and top_left_node is None and not fixed_offsets:
        add([(0, 1), (1, -1)], high=0)
        add([(n, 1), (n + 1, -1)], high=0)
    matrix = coo_matrix((vv, (ii, jj)), shape=(len(lows), size)).tocsc()
    result = milp(cost, integrality=integrality, bounds=Bounds(lb, ub),
                  constraints=LinearConstraint(matrix, lows, highs),
                  options={'time_limit': time_limit, 'mip_rel_gap': 0, 'presolve': True})
    used = 'validated heuristic incumbent'
    if result.x is not None:
        proposal = np.rint(np.column_stack([result.x[:n], result.x[n:2 * n]])).astype(int)
        validate_placement(proposal)
        if edge_length(proposal, edges) <= incumbent:
            points = proposal; used = 'HiGHS MILP incumbent'
    validate_placement(points)
    objective = edge_length(points, edges)
    bound = getattr(result, 'mip_dual_bound', None)
    bound = max(float(len(edges)), float(bound)) if bound is not None and np.isfinite(bound) else float(len(edges))
    proven = bool(result.status == 0 or objective <= bound + 1e-6)
    report = dict(method='Integer-coordinate MILP / SciPy HiGHS', solver_status=int(result.status),
                  solver_message=result.message, time_limit_seconds=time_limit,
                  grid_cols=cols, grid_rows=rows, variables=size, constraints=len(lows),
                  objective=objective, objective_metric='Unweighted Manhattan incidence length in grid units',
                  global_lower_bound=bound, relative_gap=(objective - bound) / max(objective, 1),
                  optimality_proven=proven, selected_solution=used, heuristic=heuristic,
                  distinct_grid_points_verified=len(set(map(tuple, points))) == n,
                  position_constraints=constraint_report())
    return dict(zip(node_ids, points.tolist())), report
