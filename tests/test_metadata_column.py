# backend/tests/test_metadata_column.py
"""Regression: metadata whose grouping column is not literally `condition`.

The R engines default to a `~ condition` design, so the validator must rename
whatever grouping column it resolves. Without this, the app's own demo dataset
(`group` column) failed every engine with "object 'condition' not found".
"""
import os
import sys
import unittest

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app.analysis.validator import validate_analysis_input


def _counts():
    return pd.DataFrame(
        {
            "S1": [10, 20, 30, 40],
            "S2": [11, 21, 31, 41],
            "S3": [90, 80, 70, 60],
            "S4": [91, 81, 71, 61],
        },
        index=pd.Index(["G1", "G2", "G3", "G4"], name="gene"),
    )


def _meta(column):
    return pd.DataFrame(
        {
            "sample": ["S1", "S2", "S3", "S4"],
            column: ["Control", "Control", "Treated", "Treated"],
        }
    ).set_index("sample")


class MetadataColumnTest(unittest.TestCase):
    def test_group_column_is_renamed_to_condition(self):
        _, meta, warnings = validate_analysis_input(_counts(), _meta("group"), {})
        self.assertIn("condition", meta.columns)
        self.assertEqual(list(meta["condition"]), ["Control", "Control", "Treated", "Treated"])
        self.assertTrue(any("Using 'group'" in w for w in warnings))

    def test_treatment_column_is_renamed_to_condition(self):
        _, meta, _ = validate_analysis_input(_counts(), _meta("treatment"), {})
        self.assertIn("condition", meta.columns)

    def test_condition_column_is_left_alone(self):
        _, meta, warnings = validate_analysis_input(_counts(), _meta("condition"), {})
        self.assertEqual(list(meta.columns), ["condition"])
        self.assertFalse(any("Using" in w for w in warnings))

    def test_explicit_condition_col_is_honoured(self):
        _, meta, _ = validate_analysis_input(_counts(), _meta("group"), {"condition_col": "group"})
        self.assertIn("condition", meta.columns)

    def test_single_replicate_group_still_rejected(self):
        meta = _meta("group")
        meta.loc["S3", "group"] = "Control"
        with self.assertRaises(Exception) as ctx:
            validate_analysis_input(_counts(), meta, {})
        self.assertIn("only 1 sample", str(ctx.exception))


if __name__ == "__main__":
    unittest.main()
