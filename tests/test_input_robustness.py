import unittest

import pandas as pd

from app.analysis.validator import DataValidationError, validate_analysis_input


class InputRobustnessTest(unittest.TestCase):
    def _counts(self, samples, genes=("G1", "G2", "G3", "G4")):
        return pd.DataFrame(
            {s: [10 + i, 20 + i, 5 + i, 30 + i] for i, s in enumerate(samples)},
            index=pd.Index(genes, name="Gene"),
        )

    def test_featurecounts_length_column_is_annotation_not_sample(self):
        counts = self._counts(["S1", "S2", "S3", "S4"])
        counts.insert(0, "Length", [100, 200, 300, 400])
        counts.index.name = "Geneid"
        meta = pd.DataFrame({"Sample": ["S1", "S2", "S3", "S4"],
                             "Group": ["Control", "Control", "Treated", "Treated"]})
        clean, clean_meta, warnings = validate_analysis_input(
            counts, meta,
            {"condition_col": "Group", "control_group": "Control", "treated_group": "Treated"},
        )
        self.assertNotIn("Length", clean.columns)
        self.assertEqual(list(clean.columns), ["S1", "S2", "S3", "S4"])
        self.assertTrue(any("annotation" in w.lower() for w in warnings))
        self.assertFalse(any("no metadata row" in w for w in warnings))

    def test_third_group_is_excluded_from_selected_contrast(self):
        counts = self._counts(["S1", "S2", "S3", "S4", "S5", "S6"])
        meta = pd.DataFrame({
            "Sample": ["S1", "S2", "S3", "S4", "S5", "S6"],
            "Group": ["Control", "Control", "Treated", "Treated", "Other", "Other"],
        })
        clean, clean_meta, warnings = validate_analysis_input(
            counts, meta,
            {"condition_col": "Group", "control_group": "Control", "treated_group": "Treated"},
        )
        self.assertEqual(list(clean.columns), ["S1", "S2", "S3", "S4"])
        self.assertEqual(clean_meta.index.tolist(), ["S1", "S2", "S3", "S4"])
        self.assertEqual(set(clean_meta["condition"]), {"Control", "Treated"})
        self.assertTrue(any("outside the selected" in w for w in warnings))

    def test_unselected_group_with_one_sample_does_not_block(self):
        counts = self._counts(["S1", "S2", "S3", "S4", "S5"])
        meta = pd.DataFrame({
            "Sample": ["S1", "S2", "S3", "S4", "S5"],
            "Group": ["Control", "Control", "Treated", "Treated", "Solo"],
        })
        clean, clean_meta, warnings = validate_analysis_input(
            counts, meta,
            {"condition_col": "Group", "control_group": "Control", "treated_group": "Treated"},
        )
        self.assertEqual(clean_meta.index.tolist(), ["S1", "S2", "S3", "S4"])
        self.assertEqual(clean.shape[1], 4)


    def test_single_replicate_reports_real_reason(self):
        counts = self._counts(["S1", "S2", "S3"])
        meta = pd.DataFrame({"Sample": ["S1", "S2", "S3"],
                             "Group": ["Control", "Control", "Treated"]})
        with self.assertRaises(DataValidationError) as ctx:
            validate_analysis_input(
                counts, meta,
                {"condition_col": "Group", "control_group": "Control", "treated_group": "Treated"},
            )
        message = str(ctx.exception)
        self.assertIn("treatment group 'treated' has 1 sample", message.lower())
        self.assertNotIn("Fewer than 2 matching samples", message)


if __name__ == "__main__":
    unittest.main()
