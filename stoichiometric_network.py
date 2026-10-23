#!/usr/bin/env python3
"""BiGG JSON / SBML → 化学量論行列 S → ネットワークの隣接行列（標準ライブラリのみ）。

S の行は代謝物、列は反応。負値は基質、正値は生成物。
B = (S != 0) とすると、無向の代謝物隣接行列は (B B^T > 0)、
反応隣接行列は (B^T B > 0)。どちらも対角成分を 0 にする。
有向の場合は反応の許容方向ごとに基質→生成物を結ぶ。
反応ネットワークでは、ある反応の生成物が別の反応の基質なら結ぶ。
有向行列の A[i][j] = 1 は i → j を表す。

これは構造上の接続であり、定常状態でのフラックス実現可能性や
原子の移動を表すものではない。ATP、水なども通常の代謝物として含む。
"""

import argparse
import csv
from pathlib import Path
from model_input import load_model


def stoichiometric_matrix(model):
    """BiGG JSON の辞書から (S, metabolite_ids, reaction_ids) を返す。"""
    metabolites = [m['id'] for m in model['metabolites']]
    reactions = [r['id'] for r in model['reactions']]
    index = {name: i for i, name in enumerate(metabolites)}
    S = [[0.0] * len(reactions) for _ in metabolites]
    for j, reaction in enumerate(model['reactions']):
        for name, coefficient in reaction['metabolites'].items():
            S[index[name]][j] = float(coefficient)
    return S, metabolites, reactions


def adjacency_from_stoichiometry(S, network='metabolite', directed=False,
                                 bounds=None):
    """S から 0/1 の隣接行列を返す（自己ループなし）。

    network: 'metabolite' または 'reaction'。
    bounds: 反応順の (lower_bound, upper_bound) リスト。
      有向の場合のみ使用。省略時は全反応を順方向と仮定する。
      上限 > 0 なら順方向、下限 < 0 なら逆方向を含む。
      上下限がともに 0 の反応は有向グラフで孤立扱いになる。
      無向の場合は上下限に関係なく S の非ゼロ構造を使用する。
    """
    if network not in ('metabolite', 'reaction'):
        raise ValueError('network must be metabolite or reaction')
    if not S or not S[0] or any(len(row) != len(S[0]) for row in S):
        raise ValueError('S must be a nonempty rectangular matrix')
    n_met, n_rxn = len(S), len(S[0])
    size = n_met if network == 'metabolite' else n_rxn
    A = [[0] * size for _ in range(size)]

    if not directed:
        # 非ゼロの集合を使い、符号付き S の積による相殺を避ける。
        groups = (
            [{i for i in range(n_met) if S[i][j] != 0}
             for j in range(n_rxn)] if network == 'metabolite' else
            [{j for j in range(n_rxn) if S[i][j] != 0}
             for i in range(n_met)]
        )
        for group in groups:
            for source in group:
                for target in group:
                    if source != target:
                        A[source][target] = 1
        return A

    if bounds is None:
        bounds = [(0, 1)] * n_rxn
    if len(bounds) != n_rxn or any(lo > hi for lo, hi in bounds):
        raise ValueError('bounds must contain one valid (lower, upper) per reaction')

    consumers = [set() for _ in range(n_met)]
    producers = [set() for _ in range(n_met)]
    for j, (lower, upper) in enumerate(bounds):
        substrates = {i for i in range(n_met) if S[i][j] < 0}
        products = {i for i in range(n_met) if S[i][j] > 0}
        directions = []
        if upper > 0:
            directions.append((substrates, products))
        if lower < 0:
            directions.append((products, substrates))
        for sources, targets in directions:
            if network == 'metabolite':
                for source in sources:
                    for target in targets:
                        if source != target:
                            A[source][target] = 1
            else:
                for i in sources:
                    consumers[i].add(j)
                for i in targets:
                    producers[i].add(j)

    if network == 'reaction':
        for producing, consuming in zip(producers, consumers):
            for source in producing:
                for target in consuming:
                    if source != target:
                        A[source][target] = 1
    return A


def write_matrix(path, matrix, row_ids, column_ids):
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['id', *column_ids])
        for name, row in zip(row_ids, matrix):
            writer.writerow([name, *row])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', type=Path,
                        default=Path(__file__).with_name('e_coli_core.json'))
    parser.add_argument('--network', choices=['metabolite', 'reaction'],
                        default='metabolite')
    parser.add_argument('--directed', action='store_true',
                        help='反応の上下限を反映した有向隣接行列を作成')
    parser.add_argument('--output-dir', type=Path, default=Path('network_output'))
    args = parser.parse_args()
    model = load_model(args.model)
    S, metabolites, reactions = stoichiometric_matrix(model)
    bounds = [(r['lower_bound'], r['upper_bound']) for r in model['reactions']]
    A = adjacency_from_stoichiometry(S, args.network, args.directed, bounds)
    ids = metabolites if args.network == 'metabolite' else reactions
    args.output_dir.mkdir(parents=True, exist_ok=True)
    prefix = args.model.stem
    s_path = args.output_dir / f'{prefix}_stoichiometric.csv'
    mode = 'directed' if args.directed else 'undirected'
    a_path = args.output_dir / f'{prefix}_{args.network}_{mode}_adjacency.csv'
    write_matrix(s_path, S, metabolites, reactions)
    write_matrix(a_path, A, ids, ids)
    edges = sum(map(sum, A)) // (1 if args.directed else 2)
    print(f'S: {len(metabolites)} x {len(reactions)} -> {s_path}')
    print(f'A: {len(ids)} x {len(ids)}, edges: {edges} -> {a_path}')


if __name__ == '__main__':
    main()
