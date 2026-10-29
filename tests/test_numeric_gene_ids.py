import io
import unittest

import pandas as pd

from app.analysis.validator import validate_analysis_input


class NumericGeneIdentifierTest(unittest.TestCase):
    """Regression test for GSE60450: Entrez IDs must not become sample names.

    The GEO count matrix uses numeric Entrez Gene IDs in a ``Gene`` column.
    pandas infers that column as int64, so a dtype-based check previously left
    it in the sample matrix and dropped real samples.
    """

    def setUp(self):
        self.counts = pd.read_csv(io.StringIO(
            "Gene,S1,S2,S3,S4\n"
            "497097,10,12,14,11\n"
            "100503874,3,4,2,5\n"
            "100038431,0,1,0,0\n"
            "19888,40,38,44,41\n"
        ))
        self.meta = pd.read_csv(io.StringIO(
            "Sample,condition\n"
            "S1,Control\nS2,Control\nS3,Treated\nS4,Treated\n"
        ))

    def test_gene_column_is_index_not_sample(self):
        counts, meta, warnings = validate_analysis_input(
            self.counts,
            self.meta,
            {"condition_col": "condition", "control_group": "Control", "treated_group": "Treated"},
        )
        self.assertEqual(list(counts.columns), ["S1", "S2", "S3", "S4"])
        self.assertEqual(counts.index.tolist(), ["497097", "100503874", "100038431", "19888"])
        self.assertEqual(meta.index.tolist(), ["S1", "S2", "S3", "S4"])
        self.assertFalse(any("no metadata row" in w for w in warnings))

    def test_integer_counts_survive(self):
        counts, _, _ = validate_analysis_input(
            self.counts,
            self.meta,
            {"condition_col": "condition", "control_group": "Control", "treated_group": "Treated"},
        )
        self.assertEqual(int(counts.loc["497097", "S1"]), 10)
        self.assertTrue((counts.dtypes == int).all())


if __name__ == "__main__":
    unittest.main()
