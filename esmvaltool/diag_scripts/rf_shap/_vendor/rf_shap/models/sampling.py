"""Training-sample designs for grid-cell x time records.

``sample_train`` is the single entry point. Existing methods keep their
previous behaviour. The designs below follow the cited sampling rules and
select rows that already exist in the pool; they do not invent new covariate
combinations.

Literature
----------
area-weighted
    Hansen & Hurwitz (1943), Ann. Math. Statist. 14, 333-362.
    Cochran (1977), Sampling Techniques, 3rd ed. Wiley.
    Equal-angle cell area is proportional to cos(latitude).
process-stratified
    Cochran (1977), stratified random sampling.
    Metzger et al. (2005), Global Ecol. Biogeogr. 14, 549-563.
    de Gruijter, Brus, Bierkens & Knotters (2006), Sampling for Natural
    Resource Monitoring. Springer.
cLHS
    Minasny & McBratney (2006), Comput. Geosci. 32, 1378-1388.
    doi:10.1016/j.cageo.2005.12.009
    Objective O = w1*O1 + w2*O2 + w3*O3, minimised by swap annealing.
    O1 is marginal-stratum occupancy, O2 class proportions, O3 the
    correlation matrix.
feature-space coverage
    Brus, de Gruijter & van Groenigen (2007), k-means coverage.
    Brus (2019), Geoderma 338, 464-480. doi:10.1016/j.geoderma.2018.07.036
medoid
    Kaufman & Rousseeuw (1990), Finding Groups in Data. Wiley.
maximin / farthest point
    Kennard & Stone (1969), Technometrics 11, 137-148.
    doi:10.1080/00401706.1969.10490666
    Johnson, Moore & Ylvisaker (1990), J. Stat. Plan. Inference 26, 131-148.
    doi:10.1016/0378-3758(90)90122-B
extreme / threshold
    Coles (2001), An Introduction to Statistical Modeling of Extreme Values.
    Springer. Peaks-over-threshold and tail enrichment. Inverse-probability
    weights put the sample back on the population.
pool-based active learning
    Lewis & Gale (1994), SIGIR, uncertainty sampling.
    Seung, Opper & Sompolinsky (1992), Query by committee.
    Abe & Mamitsuka (1998), query-by-bagging.
    Settles (2009), Active Learning Literature Survey, UW-Madison TR 1648.
acLHS
    Le & Vargas (2024), Comput. Geosci. 184, 105539.
    doi:10.1016/j.cageo.2024.105539
    The R package vargaslab/acLHS scores a marginal histogram, Pearson /
    Spearman / Kendall correlation, and a classical variogram. Here the
    marginal and correlation terms are the multivariate cLHS terms on X, and
    the variogram term is the classical Matheron estimator of y along time
    or space. Proposals are unique pool rows (swap annealing), not the
    continuous differential-evolution search used by that R package.
"""

from __future__ import annotations

import numpy as np
from sklearn.cluster import MiniBatchKMeans
from sklearn.preprocessing import StandardScaler


def _rng(seed: int) -> np.random.Generator:
    return np.random.default_rng(seed)


def _cap(n_samples: int | None, n: int) -> int:
    if n_samples is None or n_samples >= n:
        return n
    if n_samples <= 0:
        raise ValueError("n_samples must be positive")
    return int(n_samples)


def _as_float(values) -> np.ndarray:
    return np.asarray(values, dtype=np.float64)


def random_sample(n: int, n_samples: int | None, seed: int) -> np.ndarray:
    k = _cap(n_samples, n)
    return _rng(seed).choice(n, size=k, replace=False)


def stratified_sample(y: np.ndarray, n_samples: int | None, n_bins: int = 10, seed: int = 42) -> np.ndarray:
    n = len(y)
    k = _cap(n_samples, n)
    bins = np.digitize(y, np.nanquantile(y, np.linspace(0, 1, n_bins + 1)[1:-1]))
    rng = _rng(seed)
    chosen = []
    remaining = k
    groups = [np.where(bins == b)[0] for b in np.unique(bins)]
    groups = [g for g in groups if len(g)]
    for i, group in enumerate(groups):
        take = remaining if i == len(groups) - 1 else max(1, int(round(k * len(group) / n)))
        take = min(take, len(group), remaining)
        chosen.append(rng.choice(group, size=take, replace=False))
        remaining -= take
        if remaining <= 0:
            break
    idx = np.concatenate(chosen) if chosen else random_sample(n, k, seed)
    if len(idx) < k:
        extra = np.setdiff1d(np.arange(n), idx)
        idx = np.concatenate([idx, rng.choice(extra, size=min(k - len(idx), len(extra)), replace=False)])
    return idx[:k]


def cluster_sample(
    x: np.ndarray,
    n_samples: int | None,
    seed: int = 42,
    *,
    draw: str = "centroid",
    n_clusters: int | None = None,
) -> np.ndarray:
    """Cluster in standardized feature space and keep real pool rows.

    ``draw="centroid"`` is feature-space coverage sampling: one row nearest
    each k-means centre (Brus 2019). ``draw="medoid"`` keeps the member with
    the smallest total distance to the other members (Kaufman & Rousseeuw
    1990); clusters larger than 800 fall back to the nearest-to-centre row.
    ``draw="within"`` draws real members from each cluster instead of only
    the centre. Empty clusters are skipped.
    """
    n = len(x)
    k = _cap(n_samples, n)
    if k == n:
        return np.arange(n)
    draw = draw.lower()
    if draw not in {"centroid", "medoid", "within"}:
        raise ValueError("cluster draw must be centroid, medoid, or within")
    n_clust = k if draw != "within" else min(int(n_clusters or max(2, k // 5)), k, n)
    scaled = StandardScaler().fit_transform(x)
    km = MiniBatchKMeans(n_clusters=n_clust, random_state=seed, batch_size=min(2048, n), n_init=3)
    labels = km.fit_predict(scaled)
    centers = km.cluster_centers_
    rng = _rng(seed)
    groups = []
    for c in range(n_clust):
        members = np.where(labels == c)[0]
        if len(members):
            groups.append((c, members))
    if not groups:
        return random_sample(n, k, seed)
    if draw == "within":
        return _draw_from_groups(n, k, [members for _c, members in groups], rng, balance="equal")
    chosen = []
    for c, members in groups:
        if draw == "medoid" and len(members) <= 800:
            block = scaled[members]
            total = np.sum(np.linalg.norm(block[:, None, :] - block[None, :, :], axis=2), axis=1)
            chosen.append(int(members[int(np.argmin(total))]))
        else:
            dist = np.linalg.norm(scaled[members] - centers[c], axis=1)
            chosen.append(int(members[int(np.argmin(dist))]))
    return np.unique(chosen)


def fs_coverage_sample(x: np.ndarray, n_samples: int | None, seed: int = 42) -> np.ndarray:
    """Feature-space coverage sampling: the real row nearest each centre.

    Brus (2019) minimises the mean squared shortest distance in standardized
    covariate space by k-means, then keeps one existing record per cluster.
    """
    return cluster_sample(x, n_samples, seed, draw="centroid")


def kennard_stone(x: np.ndarray, n_samples: int | None, seed: int = 42) -> np.ndarray:
    n = len(x)
    k = _cap(n_samples, n)
    if k == n:
        return np.arange(n)
    scaled = StandardScaler().fit_transform(x)
    rng = _rng(seed)
    start = int(rng.integers(0, n))
    d0 = np.linalg.norm(scaled - scaled[start], axis=1)
    selected = [start, int(np.argmax(d0))]
    remaining = set(range(n)) - set(selected)
    min_dist = np.minimum(
        np.linalg.norm(scaled - scaled[selected[0]], axis=1),
        np.linalg.norm(scaled - scaled[selected[1]], axis=1),
    )
    while len(selected) < k and remaining:
        nxt = max(remaining, key=lambda i: min_dist[i])
        selected.append(nxt)
        remaining.remove(nxt)
        min_dist = np.minimum(min_dist, np.linalg.norm(scaled - scaled[nxt], axis=1))
    return np.asarray(selected)


def maximin_sample(x: np.ndarray, n_samples: int | None, seed: int = 42) -> np.ndarray:
    """Greedy maximin / farthest-point sample of existing rows.

    Kennard & Stone (1969) start from the farthest pair, then repeatedly add
    the candidate whose nearest already selected neighbour is as far as
    possible. That is a sequential construction of a maximin distance design
    (Johnson, Moore & Ylvisaker 1990). Pools larger than 2500 rows use two
    farthest-point steps from a random seed instead of the exact diameter,
    which is quadratic.
    """
    n = len(x)
    k = _cap(n_samples, n)
    if k == n:
        return np.arange(n)
    if k == 1:
        return np.array([int(_rng(seed).integers(0, n))])
    scaled = StandardScaler().fit_transform(x)
    first, second = _farthest_pair(scaled, seed)
    selected = [first]
    if second != first:
        selected.append(second)
    min_dist = np.linalg.norm(scaled - scaled[selected[0]], axis=1)
    for row in selected[1:]:
        min_dist = np.minimum(min_dist, np.linalg.norm(scaled - scaled[row], axis=1))
    min_dist[selected] = -np.inf
    while len(selected) < k:
        nxt = int(np.argmax(min_dist))
        selected.append(nxt)
        min_dist = np.minimum(min_dist, np.linalg.norm(scaled - scaled[nxt], axis=1))
        min_dist[nxt] = -np.inf
    return np.asarray(selected, dtype=int)


def _farthest_pair(scaled: np.ndarray, seed: int) -> tuple[int, int]:
    n = len(scaled)
    if n <= 2500:
        best = (-1.0, 0, min(1, n - 1))
        block = 400
        for i0 in range(0, n, block):
            i1 = min(n, i0 + block)
            dist = np.linalg.norm(scaled[i0:i1, None, :] - scaled[None, :, :], axis=2)
            flat = int(np.argmax(dist))
            row, col = divmod(flat, n)
            if dist[row, col] > best[0]:
                best = (float(dist[row, col]), i0 + row, col)
        return int(best[1]), int(best[2])
    start = int(_rng(seed).integers(0, n))
    second = int(np.argmax(np.linalg.norm(scaled - scaled[start], axis=1)))
    first = int(np.argmax(np.linalg.norm(scaled - scaled[second], axis=1)))
    return first, second


def _marginal_strata(x: np.ndarray, k: int) -> np.ndarray:
    from scipy.stats import rankdata

    n, p = x.shape
    strata = np.empty((n, p), dtype=int)
    for j in range(p):
        ranks = rankdata(x[:, j], method="average") - 1.0
        strata[:, j] = np.minimum((ranks * k / max(n, 1)).astype(int), k - 1)
    return strata


def _correlation(x: np.ndarray) -> np.ndarray:
    if x.shape[1] == 1 or len(x) < 3:
        return np.eye(x.shape[1])
    corr = np.corrcoef(x, rowvar=False)
    if np.ndim(corr) == 0:
        return np.eye(1)
    return np.asarray(corr, dtype=float)


def _correlation_cost(sample: np.ndarray, population: np.ndarray) -> float:
    corr = _correlation(sample)
    if corr.shape != population.shape:
        return 0.0
    diff = np.abs(corr - population)
    if diff.shape[0] < 2:
        return 0.0
    upper = diff[np.triu_indices(diff.shape[0], k=1)]
    upper = upper[np.isfinite(upper)]
    return float(upper.sum()) if len(upper) else 0.0


def _occupancy_cost(strata: np.ndarray, choice: np.ndarray, k: int) -> float:
    cost = 0.0
    for j in range(strata.shape[1]):
        counts = np.bincount(strata[choice, j], minlength=k)
        cost += float(np.abs(counts[:k] - 1).sum())
    return cost


def _class_cost(labels: np.ndarray, choice: np.ndarray, proportions: np.ndarray) -> float:
    counts = np.bincount(labels[choice], minlength=len(proportions)).astype(float)
    return float(np.abs(counts / max(len(choice), 1) - proportions).sum())


def _anneal(score_fn, n: int, k: int, rng: np.random.Generator, n_iter: int, cooling: float = 0.95) -> np.ndarray:
    """Swap one selected pool row at a time. Worse swaps are accepted with exp(-delta/T)."""
    selected = rng.choice(n, size=k, replace=False)
    pool = np.setdiff1d(np.arange(n), selected)
    current = float(score_fn(selected))
    best = selected.copy()
    best_score = current
    temperature = max(current, 1e-6)
    cycle = max(1, int(n_iter) // 20)
    for step in range(int(n_iter)):
        if len(pool) == 0:
            break
        drop_i = int(rng.integers(0, k))
        add_i = int(rng.integers(0, len(pool)))
        candidate = selected.copy()
        candidate[drop_i] = pool[add_i]
        value = float(score_fn(candidate))
        delta = value - current
        if delta <= 0.0 or rng.random() < np.exp(-delta / temperature):
            pool[add_i] = selected[drop_i]
            selected = candidate
            current = value
            if current < best_score:
                best = selected.copy()
                best_score = current
                if best_score <= 1e-12:
                    break
        if (step + 1) % cycle == 0:
            temperature = max(temperature * cooling, 1e-8)
    return np.asarray(best, dtype=int)


def clhs_sample(
    x: np.ndarray,
    n_samples: int | None,
    seed: int = 42,
    n_iter: int = 2000,
    *,
    categorical: np.ndarray | None = None,
    objective_weights: tuple[float, float, float] | None = None,
) -> np.ndarray:
    """Conditioned Latin hypercube on existing rows (Minasny & McBratney 2006).

    Each continuous column is split into ``n_samples`` equal-probability
    strata. Annealing swaps pool rows to make every stratum appear once (O1),
    to match class proportions when ``categorical`` labels are given (O2),
    and to match the population correlation matrix (O3).
    """
    x = _as_float(x)
    n, _p = x.shape
    k = _cap(n_samples, n)
    if k == n:
        return np.arange(n)
    weights = (1.0, 1.0, 1.0) if objective_weights is None else tuple(objective_weights)
    strata = _marginal_strata(x, k)
    population = _correlation(x)
    labels = None if categorical is None else _as_labels(categorical)
    proportions = None if labels is None else np.bincount(labels).astype(float) / n

    def score(choice: np.ndarray) -> float:
        o1 = _occupancy_cost(strata, choice, k)
        o2 = 0.0 if labels is None else _class_cost(labels, choice, proportions)
        o3 = _correlation_cost(x[choice], population)
        return weights[0] * o1 + weights[1] * o2 + weights[2] * o3

    return _anneal(score, n, k, _rng(seed), n_iter)


def _as_labels(categorical) -> np.ndarray:
    values = np.asarray(categorical)
    _uniq, inverse = np.unique(values, return_inverse=True)
    return inverse.astype(int)


def spatial_sample(lat: np.ndarray, lon: np.ndarray, n_samples: int | None, seed: int = 42) -> np.ndarray:
    coords = np.column_stack([lat, lon])
    return cluster_sample(coords, n_samples, seed)


def systematic_sample(n: int, n_samples: int | None, seed: int = 42) -> np.ndarray:
    k = _cap(n_samples, n)
    step = n / k
    start = float(_rng(seed).uniform(0, step))
    idx = (start + step * np.arange(k)).astype(int)
    return np.unique(np.clip(idx, 0, n - 1))


def grid_area_weights(
    lat: np.ndarray,
    time_weights: np.ndarray | None = None,
) -> np.ndarray:
    """Equal-angle cell area, proportional to cos(latitude), times optional row weights."""
    weights = np.cos(np.deg2rad(_as_float(lat)))
    weights = np.where(np.isfinite(weights), np.clip(weights, 0.0, None), 0.0)
    if time_weights is not None:
        extra = _as_float(time_weights)
        if len(extra) != len(weights):
            raise ValueError("time_weights must have one value per row")
        weights = weights * np.where(np.isfinite(extra), extra, 0.0)
    if float(weights.sum()) <= 0.0:
        weights = np.ones(len(lat), dtype=float)
    return weights


def pps_weights(area_weights: np.ndarray, indices: np.ndarray) -> np.ndarray:
    """Approximate inverse inclusion weights for a PPS sample.

    Uses the with-replacement inclusion formula 1-(1-p)^k. This reweights a
    sample toward the population; it is not an exact Horvitz-Thompson weight
    for sampling without replacement.
    """
    weights = _as_float(area_weights)
    weights = np.where(np.isfinite(weights) & (weights > 0), weights, 0.0)
    total = float(weights.sum())
    if total <= 0.0:
        raise ValueError("area weights must have a positive sum")
    probability = weights / total
    k = len(indices)
    inclusion = 1.0 - (1.0 - probability[np.asarray(indices, dtype=int)]) ** k
    return 1.0 / np.maximum(inclusion, 1e-12)


def area_weighted_sample(
    n: int,
    n_samples: int | None,
    seed: int,
    *,
    lat: np.ndarray | None = None,
    area_weights: np.ndarray | None = None,
    time_weights: np.ndarray | None = None,
) -> np.ndarray:
    """Draw existing rows with probability proportional to area x time weight.

    Hansen & Hurwitz (1943) and Cochran (1977). If ``area_weights`` is omitted,
    equal-angle grid area is cos(latitude). A zero weight is never drawn.
    """
    k = _cap(n_samples, n)
    if area_weights is None:
        if lat is None:
            raise ValueError("area-weighted sampling needs lat or area_weights")
        area_weights = grid_area_weights(lat, time_weights)
    elif time_weights is not None:
        area_weights = _as_float(area_weights) * _as_float(time_weights)
    weights = _as_float(area_weights)
    if len(weights) != n:
        raise ValueError("area_weights must have one value per row")
    weights = np.where(np.isfinite(weights) & (weights > 0), weights, 0.0)
    if int(np.count_nonzero(weights)) < k:
        raise ValueError("not enough rows with a positive area weight")
    return _rng(seed).choice(n, size=k, replace=False, p=weights / weights.sum())


def _season_codes(time: np.ndarray) -> np.ndarray:
    values = np.asarray(time)
    if np.issubdtype(values.dtype, np.datetime64):
        months = values.astype("datetime64[M]").astype(int) % 12 + 1
        codes = np.full(len(values), 3, dtype=int)
        codes[np.isin(months, (12, 1, 2))] = 0
        codes[np.isin(months, (3, 4, 5))] = 1
        codes[np.isin(months, (6, 7, 8))] = 2
        return codes
    numeric = _as_float(values)
    edges = np.unique(np.nanquantile(numeric, (0.25, 0.5, 0.75)))
    return np.digitize(numeric, edges)


def _combine_codes(parts: list[np.ndarray]) -> np.ndarray:
    ident = np.zeros(len(parts[0]), dtype=np.int64)
    for code in parts:
        code = np.asarray(code, dtype=np.int64)
        ident = ident * (int(code.max()) + 1 if len(code) else 1) + code
    return ident


def stratum_codes(
    x: np.ndarray | None = None,
    *,
    strata: np.ndarray | None = None,
    stratum_columns: list[int] | None = None,
    categorical_columns: list[int] | None = None,
    n_bins: int = 4,
    time: np.ndarray | None = None,
) -> np.ndarray:
    """Build one stratum id per row.

    Explicit ``strata`` are used as given. Otherwise the requested covariate
    columns are quantile-binned, except columns listed as categorical.
    With no columns, datetime ``time`` becomes meteorological seasons and
    numeric time becomes four quantile bins. Empty combinations never become
    groups, so a high-dimensional cross-classification does not invent rows.
    """
    if strata is not None:
        return _as_labels(strata)
    parts: list[np.ndarray] = []
    categorical = set(categorical_columns or [])
    if stratum_columns:
        if x is None:
            raise ValueError("stratum_columns needs the covariate matrix")
        for column in stratum_columns:
            values = np.asarray(x[:, int(column)])
            if int(column) in categorical:
                parts.append(_as_labels(values))
            else:
                edges = np.unique(np.nanquantile(_as_float(values), np.linspace(0, 1, n_bins + 1)[1:-1]))
                parts.append(np.digitize(_as_float(values), edges))
    if time is not None and not stratum_columns and strata is None:
        parts.append(_season_codes(time))
    if not parts:
        raise ValueError("process stratification needs strata, stratum_columns, or time")
    return _combine_codes(parts)


def _draw_from_groups(n: int, k: int, groups: list[np.ndarray], rng: np.random.Generator, balance: str) -> np.ndarray:
    balance = balance.lower()
    if balance not in {"equal", "proportional"}:
        raise ValueError("balance must be equal or proportional")
    chosen: list[np.ndarray] = []
    if balance == "proportional":
        remaining = k
        for i, group in enumerate(groups):
            take = remaining if i == len(groups) - 1 else max(1, int(round(k * len(group) / n)))
            take = min(take, len(group), remaining)
            if take > 0:
                chosen.append(rng.choice(group, size=take, replace=False))
                remaining -= take
            if remaining <= 0:
                break
    else:
        if k < len(groups):
            pick = rng.choice(len(groups), size=k, replace=False)
            return np.asarray([int(rng.choice(groups[int(i)])) for i in pick], dtype=int)
        base, extra = divmod(k, len(groups))
        order = np.argsort([-len(group) for group in groups])
        remaining = k
        for rank, group_i in enumerate(order):
            take = base + (1 if rank < extra else 0)
            take = min(take, len(groups[int(group_i)]), remaining)
            if take > 0:
                chosen.append(rng.choice(groups[int(group_i)], size=take, replace=False))
                remaining -= take
    idx = np.concatenate(chosen) if chosen else np.array([], dtype=int)
    if len(idx) < k:
        leftover = np.setdiff1d(np.arange(n), idx)
        if len(leftover):
            idx = np.concatenate([idx, rng.choice(leftover, size=min(k - len(idx), len(leftover)), replace=False)])
    return idx[:k]


def process_stratified_sample(
    x: np.ndarray,
    n_samples: int | None,
    seed: int = 42,
    n_bins: int = 4,
    *,
    strata: np.ndarray | None = None,
    stratum_columns: list[int] | None = None,
    categorical_columns: list[int] | None = None,
    time: np.ndarray | None = None,
    balance: str = "equal",
) -> np.ndarray:
    """Stratify by process or environment, then sample inside occupied strata.

    ``balance="equal"`` gives each occupied stratum the same count, which is
    the balance across known states. ``balance="proportional"`` follows the
    stratum's share of the pool. Strata with no rows are skipped.
    """
    n = len(x) if x is not None else len(strata if strata is not None else time)
    k = _cap(n_samples, n)
    codes = stratum_codes(
        x,
        strata=strata,
        stratum_columns=stratum_columns,
        categorical_columns=categorical_columns,
        n_bins=n_bins,
        time=time,
    )
    groups = [np.where(codes == code)[0] for code in np.unique(codes)]
    return _draw_from_groups(n, k, groups, _rng(seed), balance)


def _extreme_mask(
    score: np.ndarray,
    *,
    quantile: float,
    side: str,
    threshold: float | None,
    threshold_band: float | None,
) -> np.ndarray:
    score = _as_float(score)
    side = side.lower()
    if side not in {"both", "high", "low"}:
        raise ValueError("extreme side must be both, high, or low")
    if threshold is None:
        high = np.nanquantile(score, quantile)
        low = np.nanquantile(score, 1.0 - quantile)
        mask = np.zeros(len(score), dtype=bool)
        if side in {"both", "high"}:
            mask |= score >= high
        if side in {"both", "low"}:
            mask |= score <= low
        return mask
    if threshold_band is None:
        if side == "low":
            return score <= float(threshold)
        return score >= float(threshold)
    return np.abs(score - float(threshold)) <= float(threshold_band)


def extreme_sample(
    x: np.ndarray,
    y: np.ndarray,
    n_samples: int | None,
    seed: int = 42,
    *,
    fraction: float = 0.5,
    quantile: float = 0.9,
    side: str = "both",
    threshold: float | None = None,
    threshold_band: float | None = None,
    column: int | None = None,
    lat: np.ndarray | None = None,
    area_weights: np.ndarray | None = None,
    time_weights: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Oversample tails or a threshold neighbourhood, then fill from the rest.

    Coles (2001) models exceedances over a high threshold. ``fraction`` of the
    design is drawn from that pool and the rest from its complement, so the
    training distribution is intentionally not the population. The returned
    weights are inverse draws from the two pools (equal probability inside a
    pool, or PPS when area weights are supplied) and can be passed to
    ``model.fit(..., sample_weight=)`` when the target is the population.
    """
    y = _as_float(y)
    n = len(y)
    k = _cap(n_samples, n)
    score = y if column is None else _as_float(x)[:, int(column)]
    mask = _extreme_mask(score, quantile=quantile, side=side, threshold=threshold, threshold_band=threshold_band)
    extreme_rows = np.where(mask)[0]
    bulk_rows = np.where(~mask)[0]
    if len(extreme_rows) == 0 or len(bulk_rows) == 0:
        idx = random_sample(n, k, seed)
        return idx, np.ones(len(idx), dtype=float)
    n_extreme = min(len(extreme_rows), max(1, int(round(k * float(fraction)))))
    n_extreme = min(n_extreme, k)
    n_bulk = min(len(bulk_rows), k - n_extreme)
    if n_extreme + n_bulk < k:
        n_extreme = min(len(extreme_rows), k - n_bulk)
    weights = None
    if area_weights is not None or lat is not None:
        weights = grid_area_weights(lat, time_weights) if area_weights is None else _as_float(area_weights)
        if time_weights is not None and area_weights is not None:
            weights = weights * _as_float(time_weights)
    rng = _rng(seed)
    picked_ext = _weighted_choice(extreme_rows, n_extreme, weights, rng)
    picked_bulk = _weighted_choice(bulk_rows, n_bulk, weights, rng)
    idx = np.concatenate([picked_ext, picked_bulk])
    inclusion = np.empty(len(idx), dtype=float)
    inclusion[: len(picked_ext)] = _pool_weight(extreme_rows, picked_ext, weights)
    inclusion[len(picked_ext) :] = _pool_weight(bulk_rows, picked_bulk, weights)
    return idx, inclusion


def _weighted_choice(rows: np.ndarray, take: int, weights: np.ndarray | None, rng: np.random.Generator) -> np.ndarray:
    if take <= 0:
        return np.array([], dtype=int)
    if weights is None:
        return rng.choice(rows, size=take, replace=False)
    local = np.where(np.isfinite(weights[rows]) & (weights[rows] > 0), weights[rows], 0.0)
    if float(local.sum()) <= 0.0 or int(np.count_nonzero(local)) < take:
        return rng.choice(rows, size=take, replace=False)
    return rng.choice(rows, size=take, replace=False, p=local / local.sum())


def _pool_weight(pool: np.ndarray, picked: np.ndarray, weights: np.ndarray | None) -> np.ndarray:
    if len(picked) == 0:
        return np.array([], dtype=float)
    if weights is None:
        probability = len(picked) / max(len(pool), 1)
        return np.full(len(picked), 1.0 / max(probability, 1e-12))
    return _pps_inside(pool, picked, weights)


def _pps_inside(pool: np.ndarray, picked: np.ndarray, weights: np.ndarray) -> np.ndarray:
    local = np.where(np.isfinite(weights[pool]) & (weights[pool] > 0), weights[pool], 0.0)
    total = float(local.sum())
    if total <= 0.0:
        probability = len(picked) / max(len(pool), 1)
        return np.full(len(picked), 1.0 / max(probability, 1e-12))
    position = {int(row): i for i, row in enumerate(pool)}
    probability = local[np.array([position[int(row)] for row in picked])] / total
    inclusion = 1.0 - (1.0 - probability) ** len(picked)
    return 1.0 / np.maximum(inclusion, 1e-12)


def active_learning_sample(
    x: np.ndarray,
    y: np.ndarray,
    n_samples: int | None,
    seed: int = 42,
    *,
    rounds: int = 3,
    batch_size: int | None = None,
) -> np.ndarray:
    """Pool-based batch active learning on records whose target is already known.

    A small extra-trees committee is fit on the rows chosen so far. The next
    batch maximises committee disagreement (the standard deviation across
    trees; Seung et al. 1992, Abe & Mamitsuka 1998) times distance to rows
    already chosen, so near-duplicates are not added together (Settles 2009).
    The committee is refit after every batch. High uncertainty is not by
    itself evidence that the new row improves the emulator.
    """
    from sklearn.ensemble import ExtraTreesRegressor

    x = _as_float(x)
    y = _as_float(y)
    n = len(y)
    k = _cap(n_samples, n)
    if k == n:
        return np.arange(n)
    rng = _rng(seed)
    rounds = max(1, int(rounds))
    initial = min(max(k // (rounds + 1), 8), max(1, k // 2))
    initial = min(initial, k)
    selected = list(rng.choice(n, size=initial, replace=False))
    scaled = StandardScaler().fit_transform(x)
    while len(selected) < k:
        model = ExtraTreesRegressor(
            n_estimators=30,
            max_depth=8,
            min_samples_leaf=1,
            n_jobs=1,
            random_state=seed,
        )
        model.fit(x[selected], y[selected])
        remaining = np.setdiff1d(np.arange(n), np.asarray(selected, dtype=int))
        take = min(int(batch_size or max(1, int(np.ceil((k - initial) / rounds)))), k - len(selected), len(remaining))
        uncertainty = _committee_std(model, x[remaining])
        if len(remaining) > 4000:
            keep = np.argpartition(uncertainty, -4000)[-4000:]
            uncertainty = uncertainty[keep]
            candidate_rows = remaining[keep]
        else:
            candidate_rows = remaining
        picked = _diverse_argmax(scaled, np.asarray(selected, dtype=int), candidate_rows, uncertainty, take)
        selected.extend(int(row) for row in picked)
    return np.asarray(selected[:k], dtype=int)


def _committee_std(model, rows: np.ndarray) -> np.ndarray:
    if not hasattr(model, "estimators_"):
        raise TypeError("active learning needs an ensemble with estimators_")
    predictions = np.stack([estimator.predict(rows) for estimator in model.estimators_], axis=0)
    return predictions.std(axis=0)


def _diverse_argmax(
    scaled: np.ndarray,
    selected: np.ndarray,
    candidates: np.ndarray,
    uncertainty: np.ndarray,
    take: int,
) -> np.ndarray:
    chosen: list[int] = []
    scores = uncertainty.astype(float).copy()
    min_dist = _min_distance(scaled[candidates], scaled[selected])
    for _ in range(take):
        alive = np.isfinite(scores)
        if not np.any(alive):
            break
        reach = min_dist / (float(np.max(min_dist[alive])) + 1e-12)
        peak = float(np.max(scores[alive]))
        utility = np.full(len(scores), -np.inf)
        utility[alive] = (scores[alive] / (abs(peak) + 1e-12)) * reach[alive]
        pick = int(np.argmax(utility))
        chosen.append(int(candidates[pick]))
        scores[pick] = np.nan
        step = np.linalg.norm(scaled[candidates] - scaled[candidates[pick]], axis=1)
        min_dist = np.minimum(min_dist, step)
        min_dist[pick] = 0.0
    return np.asarray(chosen, dtype=int)


def _min_distance(candidates: np.ndarray, selected: np.ndarray) -> np.ndarray:
    if len(selected) == 0:
        return np.ones(len(candidates))
    nearest = np.full(len(candidates), np.inf)
    block = 512
    for i0 in range(0, len(selected), block):
        dist = np.linalg.norm(candidates[:, None, :] - selected[i0 : i0 + block][None, :, :], axis=2)
        nearest = np.minimum(nearest, dist.min(axis=1))
    return nearest


def _numeric_time(time: np.ndarray) -> np.ndarray:
    values = np.asarray(time)
    if np.issubdtype(values.dtype, np.datetime64):
        return values.astype("datetime64[s]").astype(np.float64)
    return _as_float(values)


def _dependence_coordinates(
    n: int,
    time: np.ndarray | None,
    lat: np.ndarray | None,
    lon: np.ndarray | None,
    dependence: str,
) -> np.ndarray:
    mode = dependence.lower()
    have_space = lat is not None and lon is not None
    have_time = time is not None
    if mode == "auto":
        if have_space:
            mode = "space"
        elif have_time:
            mode = "time"
        else:
            raise ValueError("acLHS needs time or lat/lon")
    if mode not in {"time", "space", "both"}:
        raise ValueError("dependence must be auto, time, space, or both")
    if mode in {"time", "both"} and not have_time:
        raise ValueError("acLHS time dependence needs a time coordinate")
    if mode in {"space", "both"} and not have_space:
        raise ValueError("acLHS space dependence needs lat and lon")
    columns = []
    if mode in {"time", "both"}:
        columns.append(_numeric_time(time))
    if mode in {"space", "both"}:
        lat_v = _as_float(lat)
        lon_v = _as_float(lon)
        columns.append(lon_v * np.cos(np.deg2rad(float(np.nanmean(lat_v)))))
        columns.append(lat_v)
    coords = np.column_stack(columns)
    if mode == "both":
        centre = coords.mean(axis=0)
        scale = coords.std(axis=0)
        scale[scale < 1e-8] = 1.0
        coords = (coords - centre) / scale
    if coords.shape[1] == 1:
        coords = np.column_stack([coords[:, 0], np.zeros(n)])
    return coords


def _variogram_edges(coords: np.ndarray, n_lags: int) -> np.ndarray:
    probe = coords if len(coords) <= 600 else coords[np.linspace(0, len(coords) - 1, 600).astype(int)]
    dist = np.linalg.norm(probe[:, None, :] - probe[None, :, :], axis=2)
    positive = dist[dist > 0]
    min_val = float(np.min(positive)) if len(positive) else 1.0
    edges = np.arange(min_val / 2.0, min_val * (n_lags + 1) + min_val * 0.5, min_val)
    if len(edges) < 2:
        edges = np.array([0.0, max(min_val, 1.0)])
    return edges


def _classical_variogram(coords: np.ndarray, values: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Matheron estimator used by geoR's classical variogram: mean half squared difference."""
    n_bins = len(edges) - 1
    gamma = np.full(n_bins, np.nan)
    n = len(values)
    if n < 2:
        return gamma
    if n <= 500:
        dist = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=2)
        half = 0.5 * (values[:, None] - values[None, :]) ** 2
        upper = np.triu_indices(n, k=1)
        bins = np.digitize(dist[upper], edges) - 1
        half = half[upper]
    else:
        # Bound the cost of a large selected set. Pairs are a fixed random
        # subset of the upper triangle, seeded from the coordinates so a
        # given subset has a stable variogram during annealing.
        left = np.repeat(np.arange(n - 1), np.arange(n - 1, 0, -1))
        right = np.concatenate([np.arange(i + 1, n) for i in range(n - 1)])
        take = min(len(left), 20000)
        pick = _rng(n * 17 + int(np.round(np.nansum(values[: min(8, n)]) * 1000))).choice(len(left), size=take, replace=False)
        left, right = left[pick], right[pick]
        dist_v = np.linalg.norm(coords[left] - coords[right], axis=1)
        half = 0.5 * (values[left] - values[right]) ** 2
        bins = np.digitize(dist_v, edges) - 1
    for bin_id in range(n_bins):
        mask = bins == bin_id
        if np.any(mask):
            gamma[bin_id] = float(np.mean(half[mask]))
    return gamma


def _variogram_cost(coords: np.ndarray, values: np.ndarray, edges: np.ndarray, reference: np.ndarray) -> float:
    gamma = _classical_variogram(coords, values, edges)
    both = np.isfinite(gamma) & np.isfinite(reference)
    if not np.any(both):
        return 0.0
    return float(np.abs(gamma[both] - reference[both]).sum())


def aclhs_sample(
    x: np.ndarray,
    y: np.ndarray,
    n_samples: int | None,
    seed: int = 42,
    n_iter: int = 2000,
    *,
    time: np.ndarray | None = None,
    lat: np.ndarray | None = None,
    lon: np.ndarray | None = None,
    dependence: str = "auto",
    n_lags: int = 8,
    objective_weights: tuple[float, float, float] | None = None,
) -> np.ndarray:
    """cLHS plus a classical variogram penalty (Le & Vargas 2024).

    O1 and the correlation term are the Minasny & McBratney criteria on the
    covariate matrix. The third term is the absolute difference between the
    classical variogram of ``y`` in the sample and in a reference subset of
    the pool. ``dependence`` is ``time``, ``space``, or ``both``; ``auto``
    uses lat/lon when they are present and otherwise time. Each term is
    divided by its value on the initial sample, so the default weights
    ``(1, 1, 1)`` give them equal influence. This is not the unnormalised
    ``(10, 1000, 0.001)`` example from the R package, and the search swaps
    real rows instead of running differential evolution.
    """
    x = _as_float(x)
    y = _as_float(y)
    n = len(y)
    k = _cap(n_samples, n)
    if k == n:
        return np.arange(n)
    coords = _dependence_coordinates(n, time, lat, lon, dependence)
    weights = (1.0, 1.0, 1.0) if objective_weights is None else tuple(float(v) for v in objective_weights)
    strata = _marginal_strata(x, k)
    population = _correlation(x)
    rng = _rng(seed)
    reference_n = min(n, 800)
    reference = rng.choice(n, size=reference_n, replace=False) if reference_n < n else np.arange(n)
    edges = _variogram_edges(coords[reference], n_lags)
    gamma_ref = _classical_variogram(coords[reference], y[reference], edges)

    def components(choice: np.ndarray) -> tuple[float, float, float]:
        o1 = _occupancy_cost(strata, choice, k)
        ocorr = _correlation_cost(x[choice], population)
        ovar = _variogram_cost(coords[choice], y[choice], edges, gamma_ref)
        return o1, ocorr, ovar

    initial = rng.choice(n, size=k, replace=False)
    base = components(initial)
    scale = tuple(value if value > 1e-8 else 1.0 for value in base)

    def score(choice: np.ndarray) -> float:
        parts = components(choice)
        return sum(weights[i] * parts[i] / scale[i] for i in range(3))

    # Continue from the same initial draw used to scale the objective.
    return _anneal_from(score, initial, n, rng, n_iter)


def _anneal_from(score_fn, selected: np.ndarray, n: int, rng: np.random.Generator, n_iter: int, cooling: float = 0.95) -> np.ndarray:
    selected = np.asarray(selected, dtype=int).copy()
    k = len(selected)
    pool = np.setdiff1d(np.arange(n), selected)
    current = float(score_fn(selected))
    best = selected.copy()
    best_score = current
    temperature = max(current, 1e-6)
    cycle = max(1, int(n_iter) // 20)
    for step in range(int(n_iter)):
        if len(pool) == 0:
            break
        drop_i = int(rng.integers(0, k))
        add_i = int(rng.integers(0, len(pool)))
        candidate = selected.copy()
        candidate[drop_i] = pool[add_i]
        value = float(score_fn(candidate))
        delta = value - current
        if delta <= 0.0 or rng.random() < np.exp(-delta / temperature):
            pool[add_i] = selected[drop_i]
            selected = candidate
            current = value
            if current < best_score:
                best = selected.copy()
                best_score = current
                if best_score <= 1e-12:
                    break
        if (step + 1) % cycle == 0:
            temperature = max(temperature * cooling, 1e-8)
    return np.asarray(best, dtype=int)


def sample_train(
    x: np.ndarray,
    y: np.ndarray,
    method: str = "random",
    n_samples: int | None = None,
    n_bins: int = 10,
    seed: int = 42,
    lat: np.ndarray | None = None,
    lon: np.ndarray | None = None,
    *,
    time: np.ndarray | None = None,
    area_weights: np.ndarray | None = None,
    time_weights: np.ndarray | None = None,
    strata: np.ndarray | None = None,
    stratum_columns: list[int] | None = None,
    categorical_columns: list[int] | None = None,
    categorical: np.ndarray | None = None,
    balance: str = "equal",
    n_clusters: int | None = None,
    extreme_fraction: float = 0.5,
    extreme_quantile: float = 0.9,
    extreme_side: str = "both",
    threshold: float | None = None,
    threshold_band: float | None = None,
    threshold_column: int | None = None,
    n_iter: int = 2000,
    active_rounds: int = 3,
    active_batch: int | None = None,
    dependence: str = "auto",
    objective_weights: tuple[float, float, float] | None = None,
    return_weights: bool = False,
):
    """Select training-row indices. ``return_weights=True`` also returns design weights.

    Weights are None for equal-probability designs. Area-weighted and extreme
    designs return inverse-probability weights so a later ``sample_weight``
    can put summaries back on the population. The fitted model is unchanged
    unless the caller passes those weights on.
    """
    method = method.lower()
    n = len(y)
    weights = None
    if method in {"none", "all"}:
        idx = np.arange(n)
    elif method == "random":
        idx = random_sample(n, n_samples, seed)
    elif method == "stratified":
        idx = stratified_sample(y, n_samples, n_bins, seed)
    elif method in {"fscs", "cluster", "coverage"}:
        idx = fs_coverage_sample(x, n_samples, seed)
    elif method in {"cluster_medoid", "medoid"}:
        idx = cluster_sample(x, n_samples, seed, draw="medoid")
    elif method in {"cluster_within"}:
        idx = cluster_sample(x, n_samples, seed, draw="within", n_clusters=n_clusters)
    elif method in {"clhs", "lhs"}:
        idx = clhs_sample(
            x,
            n_samples,
            seed,
            n_iter,
            categorical=categorical,
            objective_weights=objective_weights,
        )
    elif method in {"kennard_stone", "ks"}:
        idx = kennard_stone(x, n_samples, seed)
    elif method in {"maximin", "farthest", "farthest_point"}:
        idx = maximin_sample(x, n_samples, seed)
    elif method == "spatial":
        if lat is None or lon is None:
            raise ValueError("spatial sampling needs lat/lon")
        idx = spatial_sample(lat, lon, n_samples, seed)
    elif method == "systematic":
        idx = systematic_sample(n, n_samples, seed)
    elif method in {"area", "area_weighted", "pps"}:
        idx = area_weighted_sample(
            n,
            n_samples,
            seed,
            lat=lat,
            area_weights=area_weights,
            time_weights=time_weights,
        )
        used = grid_area_weights(lat, time_weights) if area_weights is None else _as_float(area_weights)
        if time_weights is not None and area_weights is not None:
            used = used * _as_float(time_weights)
        weights = pps_weights(used, idx)
    elif method in {"process", "process_stratified", "env_stratified"}:
        idx = process_stratified_sample(
            x,
            n_samples,
            seed,
            n_bins=n_bins,
            strata=strata,
            stratum_columns=stratum_columns,
            categorical_columns=categorical_columns,
            time=time,
            balance=balance,
        )
    elif method in {"extreme", "threshold", "tail"}:
        idx, weights = extreme_sample(
            x,
            y,
            n_samples,
            seed,
            fraction=extreme_fraction,
            quantile=extreme_quantile,
            side=extreme_side,
            threshold=threshold,
            threshold_band=threshold_band,
            column=threshold_column,
            lat=lat,
            area_weights=area_weights,
            time_weights=time_weights,
        )
    elif method in {"active", "active_learning", "pool"}:
        idx = active_learning_sample(
            x,
            y,
            n_samples,
            seed,
            rounds=active_rounds,
            batch_size=active_batch,
        )
    elif method in {"aclhs"}:
        idx = aclhs_sample(
            x,
            y,
            n_samples,
            seed,
            n_iter,
            time=time,
            lat=lat,
            lon=lon,
            dependence=dependence,
            objective_weights=objective_weights,
        )
    else:
        raise ValueError(f"Unknown sampling method: {method}")
    if return_weights:
        return idx, weights
    return idx
