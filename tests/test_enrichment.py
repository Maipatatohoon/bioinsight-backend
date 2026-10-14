# backend/tests/test_enrichment.py
"""
Over-representation (hypergeometric) enrichment + Benjamini-Hochberg, and the
/api/enrichment endpoint. Uses synthetic annotations (no network).
"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app.analysis.enrichment import enrich_genes, benjamini_hochberg


def _annotations():
    # universe = 100 genes; TERM_A is enriched in the query; TERM_B is not.
    universe = [f"G{i}" for i in range(100)]
    term_a = {"name": "Enriched term", "genes": [f"G{i}" for i in range(10)]}        # query hits 8/10
    term_b = {"name": "Background term", "genes": [f"G{i}" for i in range(50, 60)]}  # query hits 0
    term_small = {"name": "Tiny term", "genes": ["G0", "G1"]}                        # filtered by min size
    return {"GO:A": term_a, "GO:B": term_b, "GO:SMALL": term_small}, universe


class TestEnrichment(unittest.TestCase):
    def test_returns_only_enriched_terms_with_correct_shape(self):
        annotations, universe = _annotations()
        query = [f"G{i}" for i in range(8)] + [f"G{i}" for i in range(20, 30)]
        results = enrich_genes(query, annotations, universe=universe, min_term_size=3)
        terms = {r["term"] for r in results}
        self.assertIn("GO:A", terms)
        self.assertNotIn("GO:B", terms)       # no overlap
        self.assertNotIn("GO:SMALL", terms)   # below min size
        top = results[0]
        self.assertEqual(top["term"], "GO:A")
        self.assertEqual(top["overlap"], 8)
        self.assertLess(top["pvalue"], 0.05)
        self.assertIn("padj", top)
        self.assertEqual(sorted(top["genes"]), [f"G{i}" for i in range(8)])

    def test_empty_query_returns_empty(self):
        annotations, universe = _annotations()
        self.assertEqual(enrich_genes([], annotations, universe=universe), [])

    def test_benjamini_hochberg_is_monotonic_and_bounded(self):
        p = [0.001, 0.01, 0.02, 0.5, 0.9]
        adj = benjamini_hochberg(p)
        self.assertEqual(len(adj), len(p))
        self.assertTrue(all(0 <= a <= 1 for a in adj))
        # adjusted p-values are non-decreasing when sorted by raw p
        ordered = [a for _, a in sorted(zip(p, adj))]
        self.assertEqual(ordered, sorted(ordered))
        self.assertGreaterEqual(adj[0], p[0])


class TestEnrichmentEndpoint(unittest.TestCase):
    def setUp(self):
        from app.main import app
        from starlette.testclient import TestClient
        self.client = TestClient(app)

    def test_endpoint_returns_results(self):
        annotations, universe = _annotations()
        resp = self.client.post("/api/enrichment", json={
            "genes": [f"G{i}" for i in range(8)],
            "annotations": annotations,
            "universe": universe,
        })
        self.assertEqual(resp.status_code, 200, resp.text)
        body = resp.json()
        self.assertIn("results", body)
        self.assertTrue(any(r["term"] == "GO:A" for r in body["results"]))


if __name__ == "__main__":
    unittest.main()
