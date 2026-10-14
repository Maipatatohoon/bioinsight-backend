# backend/app/analysis/enrichment.py
"""
Over-representation analysis (hypergeometric test) with Benjamini-Hochberg FDR.
Pure-python + scipy so it is unit-testable and has no network dependency.
"""
from typing import Any, Dict, Iterable, List, Optional, Sequence

from scipy.stats import hypergeom


def benjamini_hochberg(pvalues: Sequence[float]) -> List[float]:
    """Return BH-adjusted p-values (monotonic in the raw p order)."""
    m = len(pvalues)
    if m == 0:
        return []
    order = sorted(range(m), key=lambda i: pvalues[i])
    adjusted = [0.0] * m
    prev = 1.0
    for rank in range(m, 0, -1):
        idx = order[rank - 1]
        value = pvalues[idx] * m / rank
        prev = min(prev, value)
        adjusted[idx] = min(prev, 1.0)
    return adjusted


def _term_genes(value: Any) -> set:
    if isinstance(value, dict):
        return set(value.get("genes") or [])
    return set(value or [])


def _term_name(term: str, value: Any) -> str:
    if isinstance(value, dict):
        return value.get("name") or term
    return term


def enrich_genes(
    genes: Iterable[str],
    annotations: Dict[str, Any],
    universe: Optional[Iterable[str]] = None,
    min_term_size: int = 3,
    max_term_size: int = 500,
) -> List[Dict[str, Any]]:
    """Hypergeometric over-representation of `genes` within `annotations`.

    `annotations` maps a term id -> either {"name", "genes"} or a list of genes.
    `universe` is the background gene set; if omitted it is the union of all
    annotated genes.
    """
    query = {g for g in genes if g}
    term_sets = {t: _term_genes(v) for t, v in annotations.items()}

    if universe is None:
        background = set()
        for s in term_sets.values():
            background |= s
    else:
        background = {g for g in universe if g}
    if not background:
        return []

    query = query & background
    if not query:
        return []

    n_query = len(query)
    n_background = len(background)

    results: List[Dict[str, Any]] = []
    for term, geneset in term_sets.items():
        geneset = geneset & background
        size = len(geneset)
        if size < min_term_size or size > max_term_size:
            continue
        overlap = query & geneset
        k = len(overlap)
        if k == 0:
            continue
        # P(X >= k) under the hypergeometric null
        pvalue = float(hypergeom.sf(k - 1, n_background, size, n_query))
        results.append({
            "term": term,
            "name": _term_name(term, annotations[term]),
            "pvalue": pvalue,
            "overlap": k,
            "termSize": size,
            "querySize": n_query,
            "backgroundSize": n_background,
            "genes": sorted(overlap),
        })

    if not results:
        return []

    for r, adj in zip(results, benjamini_hochberg([r["pvalue"] for r in results])):
        r["padj"] = adj
    results.sort(key=lambda r: r["pvalue"])
    return results
