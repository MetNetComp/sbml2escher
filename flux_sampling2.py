"""Bounded-polytope hit-and-run sampling using NumPy, SciPy and HiGHS.

FVA identifies the affine hull, including reactions fixed by mass balance.
Warm-up covariance is adapted only before production. Frozen, symmetric random
line directions and uniform chord draws have a uniform stationary distribution.
Finite chains are correlated; split R-hat and ESS diagnose incomplete mixing.
"""
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import linprog
from stoichiometric_network import stoichiometric_matrix


def diagnostics(chains):
    """Classical split R-hat and initial-positive/monotone-sequence ESS."""
    m, n, p = chains.shape
    h = n // 2
    split = np.concatenate([chains[:, :h], chains[:, -h:]], axis=0)
    m, n, p = split.shape
    within = split.var(axis=1, ddof=1).mean(axis=0)
    between = n * split.mean(axis=1).var(axis=0, ddof=1)
    variance = (n - 1) / n * within + between / n
    active = variance > 1e-20
    rhat = np.ones(p)
    rhat[active] = np.sqrt(variance[active] / np.maximum(within[active], 1e-30))
    centered = split - split.mean(axis=1, keepdims=True)
    fft = np.fft.rfft(centered, n=2 * n, axis=1)
    acov = np.fft.irfft(fft * np.conj(fft), n=2 * n, axis=1)[:, :n].mean(axis=0) / n
    rho = 1 - (within[None] - acov) / np.maximum(variance[None], 1e-30)
    rho[0] = 1
    ess = np.full(p, float(m * n))
    for j in np.flatnonzero(active):
        pairs = rho[:n // 2 * 2, j].reshape(-1, 2).sum(axis=1)
        positive = np.flatnonzero(pairs <= 0)
        pairs = pairs[:positive[0]] if len(positive) else pairs
        tau = max(1.0, -1.0 + 2 * np.minimum.accumulate(pairs).sum()) if len(pairs) else 1.0
        ess[j] = min(m * n, m * n / tau)
    return rhat, ess


def sample_fluxes(model, samples=2000, chains=4, thinning=100, burn_in=10000,
                  seed=0, tolerance=1e-7, progress=print):
    if samples < 20 or chains < 2 or thinning < 1 or burn_in < 1:
        raise ValueError('Require samples >= 20, chains >= 2, thinning >= 1, burn-in >= 1')
    S, mids, rids = stoichiometric_matrix(model)
    S = np.asarray(S, dtype=float)
    lower = np.array([r['lower_bound'] for r in model['reactions']], dtype=float)
    upper = np.array([r['upper_bound'] for r in model['reactions']], dtype=float)
    if not np.all(np.isfinite([lower, upper])):
        raise ValueError('Uniform sampling requires finite flux bounds')
    n = len(rids)
    def lp(c):
        result = linprog(c, A_eq=S, b_eq=np.zeros(len(mids)), bounds=list(zip(lower, upper)),
                         method='highs', options={'primal_feasibility_tolerance': 1e-9,
                                                  'dual_feasibility_tolerance': 1e-9})
        if not result.success:
            raise ValueError('Flux-space LP failed: ' + result.message)
        return result.x
    warmup, fmin, fmax = [], [], []
    for j in range(n):
        unit = np.eye(1, n, j)[0]
        a, b = lp(unit), lp(-unit)
        warmup.extend([a, b]); fmin.append(a[j]); fmax.append(b[j])
    fmin, fmax = np.array(fmin), np.array(fmax)
    fixed = fmax - fmin <= tolerance
    affine = np.vstack([S, np.eye(n)[fixed]])
    N = null_space(affine)
    origin = np.mean(warmup, axis=0)
    progress(f'Sampling: affine dimension {N.shape[1]}, {int(fixed.sum())} fixed reactions', flush=True)
    rng = np.random.default_rng(seed)
    if N.shape[1] == 0:
        result = np.broadcast_to(origin, (chains, samples, n)).copy()
    else:
        warm = (np.array(warmup) - origin) @ N
        def transform(z):
            covariance = np.atleast_2d(np.cov(z, rowvar=False))
            values, vectors = np.linalg.eigh(covariance)
            floor = max(float(values.max()) * 1e-10, 1e-14)
            return N @ (vectors * np.sqrt(np.maximum(values, floor)))
        B = transform(warm)
        states = np.broadcast_to(origin, (chains, n)).copy()
        def step():
            direction = rng.normal(size=(chains, B.shape[1])) @ B.T
            norm = np.linalg.norm(direction, axis=1)
            if np.any(norm < 1e-15):
                raise RuntimeError('Degenerate sampling direction')
            direction /= norm[:, None]
            movable = np.abs(direction) > 1e-14
            a = np.divide(lower - states, direction, out=np.full_like(states, -np.inf), where=movable)
            b = np.divide(upper - states, direction, out=np.full_like(states, np.inf), where=movable)
            lo = np.max(np.where(movable, np.minimum(a, b), -np.inf), axis=1)
            hi = np.min(np.where(movable, np.maximum(a, b), np.inf), axis=1)
            if np.any(~np.isfinite(lo + hi)) or np.any(hi < lo - tolerance):
                raise RuntimeError('Invalid hit-and-run chord; increase numerical precision')
            alpha = lo + rng.random(chains) * np.maximum(hi - lo, 0)
            states[:] += alpha[:, None] * direction
        # Covariance adaptation uses discarded draws only. Production freezes B.
        for phase in range(4):
            history = []
            for i in range(max(1000, burn_in // 2)):
                step()
                if i % 10 == 0:
                    history.extend((states - origin) @ N)
            B = transform(np.array(history))
        for _ in range(burn_in):
            step()
        result = np.empty((chains, samples, n))
        for i in range(samples):
            for _ in range(thinning):
                step()
            result[:, i] = states
            if (i + 1) % max(1, samples // 4) == 0:
                progress(f'Sampling: {i + 1}/{samples} samples per chain', flush=True)
    residual = float(np.max(np.abs(result @ S.T)))
    violation = float(max(0, np.max(lower - result), np.max(result - upper)))
    if residual > tolerance or violation > tolerance or not np.all(np.isfinite(result)):
        raise RuntimeError(f'Sample validation failed: mass balance {residual}, bounds {violation}')
    rhat, ess = diagnostics(result)
    medians = np.median(result.reshape(-1, n), axis=0)
    report = dict(method='Random-direction hit-and-run with discarded covariance adaptation',
                  stationary_target='Uniform volume in the affine feasible flux polytope',
                  bounds_modified=False, objective_conditioning=False, seed=seed,
                  chains=chains, samples_per_chain=samples, total_samples=chains * samples,
                  thinning=thinning, burn_in=burn_in, adaptation_phases=4,
                  affine_dimension=N.shape[1], fixed_reaction_ids=[rids[i] for i in np.flatnonzero(fixed)],
                  mass_balance_max_abs=residual, bound_violation_max=violation,
                  max_split_rhat=float(rhat.max()), min_effective_sample_size=float(ess.min()),
                  diagnostic_criteria=dict(max_split_rhat=1.05, min_effective_sample_size=100),
                  mixing_checks_passed=bool(rhat.max() <= 1.05 and ess.min() >= 100),
                  limitation='Finite correlated MCMC samples, not proof of convergence; bounds and internal cycles affect medians.')
    stats = [dict(reaction_id=rid, median=float(medians[j]), priority=abs(float(medians[j])),
                  q05=float(np.quantile(result[:, :, j], .05)), q95=float(np.quantile(result[:, :, j], .95)),
                  fva_min=float(fmin[j]), fva_max=float(fmax[j]), split_rhat=float(rhat[j]),
                  ess=float(ess[j]), chain_medians=np.median(result[:, :, j], axis=1).tolist())
             for j, rid in enumerate(rids)]
    return result, stats, report
