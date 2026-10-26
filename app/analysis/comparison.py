"""Saved DESeq2-vs-edgeR comparison.

The comparison is deliberately computed from the two result tables already
stored on the job. It never re-runs either R analysis, and it does not invent
a combined p-value: a gene is a consensus call only when both methods call it
significant (adjusted p-value and fold-change thresholds) and agree on the
direction of change.
"""
from __future__ import annotations

import math
import time
from typing import Any, Dict, Iterable, List, Optional, Tuple


def _number(value: Any) -> Optional[float]:
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return number


def _gene(row: Dict[str, Any]) -> Optional[str]:
    for key in ("Gene", "gene", "Gene_ID", "gene_id", "id", "ID"):
        value = row.get(key)
        if value is not None and str(value).strip():
            return str(value).strip()
    return None


def _index(rows: Iterable[Dict[str, Any]], label: str) -> Tuple[Dict[str, Dict[str, Any]], int]:
    indexed: Dict[str, Dict[str, Any]] = {}
    duplicates = 0
    for row in rows or []:
        gene = _gene(row)
        if not gene:
            continue
        if gene in indexed:
            duplicates += 1
            continue
        indexed[gene] = row
    return indexed, duplicates


def _call(row: Optional[Dict[str, Any]], padj_threshold: float, fc_threshold: float) -> Dict[str, Any]:
    if not row:
        return {
            "log2fc": None,
            "pvalue": None,
            "padj": None,
            "significant": False,
            "direction": None,
        }
    log2fc = _number(row.get("log2FoldChange", row.get("logFC", row.get("log2fc"))))
    pvalue = _number(row.get("pvalue", row.get("PValue", row.get("P.Value"))))
    padj = _number(row.get("padj", row.get("FDR", row.get("adj.P.Val"))))
    # Significance requires the adjusted p-value. Falling back to a raw
    # p-value would silently evade multiple-testing correction.
    metric = padj
    significant = (
        log2fc is not None
        and metric is not None
        and metric <= padj_threshold
        and abs(log2fc) >= fc_threshold
        and log2fc != 0
    )
    direction = None
    if significant and log2fc > 0:
        direction = "up"
    elif significant and log2fc < 0:
        direction = "down"
    return {
        "log2fc": log2fc,
        "pvalue": pvalue,
        "padj": padj,
        "significant": bool(significant),
        "direction": direction,
    }


def compute_comparison(
    deseq2_rows: List[Dict[str, Any]],
    edger_rows: List[Dict[str, Any]],
    padj_threshold: float = 0.05,
    log2fc_threshold: float = 1.0,
) -> Dict[str, Any]:
    if not math.isfinite(padj_threshold) or not 0 < padj_threshold <= 1:
        raise ValueError("padj_threshold must be between 0 and 1.")
    if not math.isfinite(log2fc_threshold) or log2fc_threshold < 0:
        raise ValueError("log2fc_threshold must be non-negative.")

    deseq2, deseq2_duplicates = _index(deseq2_rows, "DESeq2")
    edger, edger_duplicates = _index(edger_rows, "edgeR")
    genes = sorted(set(deseq2) | set(edger))

    rows: List[Dict[str, Any]] = []
    summary = {
        "genes_compared": len(genes),
        "deseq2_significant": 0,
        "edger_significant": 0,
        "shared_significant": 0,
        "deseq2_only": 0,
        "edger_only": 0,
        "direction_agree": 0,
        "direction_disagree": 0,
        "consensus": 0,
        "consensus_up": 0,
        "consensus_down": 0,
        "not_significant": 0,
        "deseq2_missing_in_edger": 0,
        "edger_missing_in_deseq2": 0,
    }

    for gene in genes:
        d = _call(deseq2.get(gene), padj_threshold, log2fc_threshold)
        e = _call(edger.get(gene), padj_threshold, log2fc_threshold)
        d_sig = d["significant"]
        e_sig = e["significant"]
        shared = d_sig and e_sig
        agree = shared and d["direction"] == e["direction"]
        disagree = shared and d["direction"] != e["direction"]
        consensus = bool(agree and d["direction"] in ("up", "down"))

        if d_sig:
            summary["deseq2_significant"] += 1
        if e_sig:
            summary["edger_significant"] += 1
        if shared:
            summary["shared_significant"] += 1
        if d_sig and not e_sig:
            summary["deseq2_only"] += 1
        if e_sig and not d_sig:
            summary["edger_only"] += 1
        if agree:
            summary["direction_agree"] += 1
        if disagree:
            summary["direction_disagree"] += 1
        if consensus:
            summary["consensus"] += 1
            if d["direction"] == "up":
                summary["consensus_up"] += 1
            else:
                summary["consensus_down"] += 1
        if not d_sig and not e_sig:
            summary["not_significant"] += 1
        if gene not in edger:
            summary["deseq2_missing_in_edger"] += 1
        if gene not in deseq2:
            summary["edger_missing_in_deseq2"] += 1

        rows.append({
            "gene": gene,
            "deseq2_log2fc": d["log2fc"],
            "deseq2_pvalue": d["pvalue"],
            "deseq2_padj": d["padj"],
            "deseq2_significant": d_sig,
            "deseq2_direction": d["direction"],
            "edger_log2fc": e["log2fc"],
            "edger_pvalue": e["pvalue"],
            "edger_padj": e["padj"],
            "edger_significant": e_sig,
            "edger_direction": e["direction"],
            "shared_significant": shared,
            "consensus": consensus,
            "direction_agreement": (
                "agree" if agree else "disagree" if disagree
                else "one_method" if (d_sig or e_sig) else "not_significant"
            ),
        })

    return {
        "created_at": time.time(),
        "parameters": {
            "padj_threshold": padj_threshold,
            "log2fc_threshold": log2fc_threshold,
            "significance_metric": "adjusted p-value (BH FDR); raw p-values are reported but never used for significance when padj is unavailable",
            "consensus_rule": (
                "Consensus = significant in both DESeq2 and edgeR at the stated "
                "thresholds and the same log2 fold-change direction. Genes "
                "significant in only one method remain method-specific."
            ),
        },
        "summary": summary,
        "duplicates": {
            "deseq2": deseq2_duplicates,
            "edger": edger_duplicates,
        },
        "rows": rows,
    }
