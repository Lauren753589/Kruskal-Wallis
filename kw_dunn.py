"""
Kruskal-Wallis test with Dunn's post-hoc test.

Shared by kruskal_dunn_app.py and kruskal_dunn.py. Needs numpy, pandas and scipy only.

Groups are passed as a dict {group name: array of values}; the dict order is the group order.
"""
import itertools

import numpy as np
import pandas as pd
from scipy import stats


# =============================================================================
# Multiple-comparison corrections
# =============================================================================
def holm(pvals):
    """Holm step-down adjustment."""
    p = np.asarray(pvals, float)
    m = len(p)
    adj, running = np.empty(m), 0.0
    for rank, i in enumerate(np.argsort(p, kind="stable")):
        running = max(running, min(1.0, (m - rank) * p[i]))
        adj[i] = running
    return adj


def bonferroni(pvals):
    p = np.asarray(pvals, float)
    return np.minimum(1.0, p * len(p))


def benjamini_hochberg(pvals):
    """Benjamini-Hochberg step-up adjustment (controls the false discovery rate)."""
    p = np.asarray(pvals, float)
    m = len(p)
    adj, running = np.empty(m), 1.0
    for k, i in enumerate(np.argsort(p, kind="stable")[::-1]):  # largest p first
        running = min(running, p[i] * m / (m - k))
        adj[i] = running
    return adj


CORRECTIONS = {
    "Holm": holm,
    "Bonferroni": bonferroni,
    "Benjamini-Hochberg": benjamini_hochberg,
    "None": lambda p: np.asarray(p, float),
}


# =============================================================================
# Getting groups out of a table
# =============================================================================
def guess_columns(df):
    """Best guess at (group column, value column, specimen column or None)."""
    numeric = [c for c in df.columns if pd.to_numeric(df[c], errors="coerce").notna().mean() > 0.9]
    text = [c for c in df.columns if c not in numeric]
    value = numeric[-1] if numeric else None
    group = text[0] if text else next((c for c in df.columns if c != value), None)
    specimen = next((c for c in df.columns if c not in (group, value)
                     and any(w in str(c).lower() for w in ("specimen", "sample", "replicate", "subject"))), None)
    return group, value, specimen


def groups_from_table(df, group_col, value_col, specimen_col=None, order=None):
    """Return {group: values}.

    If `specimen_col` is given, values that share a group and specimen are averaged first,
    so each specimen counts once. `order` picks and orders the groups (default: file order).
    """
    d = df[[c for c in (group_col, specimen_col, value_col) if c]].copy()
    d[value_col] = pd.to_numeric(d[value_col], errors="coerce")
    d = d.dropna(subset=[group_col, value_col])
    d[group_col] = d[group_col].astype(str).str.strip()
    if specimen_col:
        d = d.groupby([group_col, specimen_col], sort=False)[value_col].mean().reset_index()
    names = order if order is not None else list(dict.fromkeys(d[group_col]))
    out = {g: d.loc[d[group_col] == g, value_col].to_numpy(float) for g in names}
    return {g: v for g, v in out.items() if len(v)}


def comparisons(names, control=None):
    """All pairs, or each group against `control`."""
    if control is None:
        return list(itertools.combinations(names, 2))
    return [(control, g) for g in names if g != control]


# =============================================================================
# The tests
# =============================================================================
def _pooled_ranks(groups):
    """All values, and each group's ranks within the pooled data."""
    values = np.concatenate([np.asarray(v, float) for v in groups.values()])
    ranks = stats.rankdata(values)
    out, start = {}, 0
    for g, v in groups.items():
        out[g] = ranks[start:start + len(v)]
        start += len(v)
    return values, out


def describe(groups):
    """n, median, quartiles, range and mean rank per group."""
    _, ranks = _pooled_ranks(groups)
    rows = []
    for g, v in groups.items():
        v = np.asarray(v, float)
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        rows.append({"Group": g, "n": len(v), "Median": med, "Q1": q1, "Q3": q3,
                     "Min": v.min(), "Max": v.max(), "Mean rank": ranks[g].mean()})
    return pd.DataFrame(rows).set_index("Group")


def kruskal_wallis(groups):
    """H, degrees of freedom, p, N and epsilon squared (H / (N - 1))."""
    H, p = stats.kruskal(*groups.values())
    N = sum(len(v) for v in groups.values())
    return {"H": H, "df": len(groups) - 1, "p": p, "N": N, "epsilon2": H / (N - 1)}


def dunn(groups, pairs, correction="Holm"):
    """Dunn's test on ranks pooled across all groups, with tie correction."""
    values, ranks = _pooled_ranks(groups)
    N = len(values)
    _, counts = np.unique(values, return_counts=True)
    tie = np.sum(counts.astype(float) ** 3 - counts)

    rows = []
    for a, b in pairs:
        se = np.sqrt((N * (N + 1) / 12 - tie / (12 * (N - 1))) * (1 / len(ranks[a]) + 1 / len(ranks[b])))
        diff = ranks[a].mean() - ranks[b].mean()
        z = diff / se
        rows.append({"Group A": a, "Group B": b, "Mean rank difference": diff,
                     "z": z, "p": 2 * stats.norm.sf(abs(z))})
    res = pd.DataFrame(rows, columns=["Group A", "Group B", "Mean rank difference", "z", "p"])
    res["p (adjusted)"] = CORRECTIONS[correction](res["p"]) if len(res) else []
    return res


def stars(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


# =============================================================================
# Plain-text report
# =============================================================================
def report(groups, comparison_sets, correction="Holm", title=None):
    """Text summary. `comparison_sets` is a list of (heading, pairs)."""
    w = max(len(g) for g in groups)
    lines = [f"===== {title} =====" if title else ""]
    for g, r in describe(groups).iterrows():
        lines.append(f"{g:{w}s}  n = {int(r['n']):2d}  median = {r['Median']:.4g}  "
                     f"IQR = {r['Q1']:.4g} to {r['Q3']:.4g}")

    kw = kruskal_wallis(groups)
    lines.append(f"Kruskal-Wallis H({kw['df']}) = {kw['H']:.2f}, p = {kw['p']:.3g}, "
                 f"epsilon^2 = {kw['epsilon2']:.2f}")

    for heading, pairs in comparison_sets:
        res = dunn(groups, pairs, correction)
        lines.append(f" Dunn's post-hoc, {correction} correction: {heading} ({len(res)} comparisons)")
        for _, r in res.iterrows():
            lines.append(f"   {r['Group A']:{w}s} vs {r['Group B']:{w}s}  z = {r['z']:6.2f}  "
                         f"p = {r['p']:.4f}  p(adj) = {r['p (adjusted)']:.4f}")
    _, ranks = _pooled_ranks(groups)
    lines.append(" Mean ranks: " + ", ".join(f"{g} {r.mean():.1f}" for g, r in ranks.items()))
    return "\n".join(lines).strip("\n")
