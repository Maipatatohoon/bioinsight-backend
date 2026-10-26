# backend/app/analysis/validator.py
"""
Data Validation Engine for BioInsight 360 Backend.
Enforces scientific correctness, sample consistency, and input integrity.
"""
import pandas as pd
import numpy as np
from typing import Dict, Any, Tuple, List

class DataValidationError(Exception):
    """Custom exception for scientific data validation errors."""
    def __init__(self, message: str, details: Dict[str, Any] = None):
        super().__init__(message)
        self.message = message
        self.details = details or {}

def validate_analysis_input(counts_df: pd.DataFrame, meta_df: pd.DataFrame, params: Dict[str, Any]) -> Tuple[pd.DataFrame, pd.DataFrame, List[str]]:
    """
    Validates count matrix and metadata for DE analysis.
    Returns cleaned (counts_df, meta_df, warnings).
    """
    warnings = []
    
    # 1. Basic Structure Validation
    if counts_df.empty:
        raise DataValidationError("Count matrix is empty.")
    if meta_df.empty:
        raise DataValidationError("Metadata matrix is empty.")

    # Normalize gene identifier column if counts_df came with gene column instead of index
    if counts_df.columns[0].lower() in ['gene', 'gene_id', 'gene_name', 'symbol', 'id', 'ensembl_id'] or counts_df.index.dtype == 'int64':
        first_col = counts_df.columns[0]
        # Check if first column contains non-numeric identifiers while other columns are numeric/counts
        if not pd.api.types.is_numeric_dtype(counts_df[first_col]):
            counts_df = counts_df.set_index(first_col)

    # Normalize sample identifier column if meta_df came with sample column instead of index
    if meta_df.columns[0].lower() in ['sample', 'sample_id', 'id', 'run', 'gsm', 'srr'] or meta_df.index.dtype == 'int64':
        first_col = meta_df.columns[0]
        if not pd.api.types.is_numeric_dtype(meta_df[first_col]):
            meta_df = meta_df.set_index(first_col)

    # Clean string indices / column headers
    counts_df.index = counts_df.index.astype(str).str.strip()
    counts_df.columns = counts_df.columns.astype(str).str.strip()
    meta_df.index = meta_df.index.astype(str).str.strip()

    # 2. Check Duplicate Gene IDs / Sample IDs
    if counts_df.index.duplicated().any():
        dup_count = counts_df.index.duplicated().sum()
        warnings.append(f"Found {dup_count} duplicate gene IDs. Aggregating duplicates by summing counts.")
        counts_df = counts_df.groupby(counts_df.index).sum()
    if counts_df.columns.duplicated().any():
        raise DataValidationError("Count matrix contains duplicate sample columns. Rename them so every sample is unique.")
    if meta_df.index.duplicated().any():
        duplicated = meta_df.index[meta_df.index.duplicated()].unique().tolist()[:5]
        raise DataValidationError(
            "Metadata contains duplicate sample identifiers: " + ", ".join(map(str, duplicated))
        )

    # 3. Sample Overlap Check
    common_samples = counts_df.columns.intersection(meta_df.index)
    counts_only = [s for s in counts_df.columns if s not in common_samples]
    meta_only = [s for s in meta_df.index if s not in common_samples]
    if counts_only:
        warnings.append(
            f"{len(counts_only)} count-matrix sample(s) had no metadata row and were dropped: "
            + ", ".join(map(str, counts_only[:8]))
            + (" ..." if len(counts_only) > 8 else "")
        )
    if meta_only:
        warnings.append(
            f"{len(meta_only)} metadata sample(s) had no count column and were dropped: "
            + ", ".join(map(str, meta_only[:8]))
            + (" ..." if len(meta_only) > 8 else "")
        )
    if len(common_samples) < 4:
        raise DataValidationError(
            f"Fewer than 4 matching samples between count matrix and metadata (found {len(common_samples)}: {list(common_samples)}). "
            "RNA-seq differential expression requires matching sample names across files."
        )

    counts_df = counts_df[common_samples]
    meta_df = meta_df.loc[common_samples]

    # 4. Check Non-Numeric or Missing Values in Counts
    if counts_df.isnull().any().any():
        affected = counts_df.isnull().sum()
        affected = affected[affected > 0].index.astype(str).tolist()[:8]
        raise DataValidationError(
            "Count matrix contains missing values (NaN/blank). A missing RNA-seq count is not "
            "zero and is not imputed by BioInsight 360. Remove or fix the affected sample "
            "column(s): " + ", ".join(affected)
        )

    try:
        counts_numeric = counts_df.apply(pd.to_numeric)
    except Exception as e:
        raise DataValidationError("Count matrix contains non-numeric values that cannot be converted to counts.")

    # 5. Check Negative / Non-finite Values
    if not np.isfinite(counts_numeric.to_numpy(dtype=float)).all():
        raise DataValidationError(
            "Count matrix contains infinite or non-finite values. Fix those cells before analysis."
        )
    if (counts_numeric < 0).any().any():
        raise DataValidationError(
            "Count matrix contains negative values. Differential expression engines (DESeq2, edgeR) require non-negative raw integer count data."
        )

    # 6. Reject normalized expression values (TPM/FPKM) instead of silently rounding
    rounded = counts_numeric.round()
    if ((counts_numeric - rounded).abs() > 1e-6).any().any():
        raise DataValidationError(
            "Non-integer expression values detected. DESeq2 and edgeR require raw, un-normalized integer counts; "
            "TPM/FPKM or other normalized values are not valid inputs. Upload the gene-level raw count matrix."
        )
    counts_numeric = counts_numeric.round().astype(int)

    # 7. Biological Replicates Check
    condition_col = params.get("condition_col", "condition")
    if condition_col not in meta_df.columns:
        # Try to find first non-sample-id column
        possible_cols = [c for c in meta_df.columns if c.lower() not in ['sample', 'sample_id', 'id']]
        if possible_cols:
            condition_col = possible_cols[0]
            warnings.append(f"Condition column '{params.get('condition_col')}' not found. Using '{condition_col}' as design factor.")
        else:
            raise DataValidationError(f"Condition column '{condition_col}' not found in metadata columns: {list(meta_df.columns)}.")

    # Explicitly excluded samples must not become a factor level in the model.
    meta_df[condition_col] = meta_df[condition_col].astype(str).str.strip()
    group_values = meta_df[condition_col]
    excluded_mask = group_values.str.lower().isin(
        ["exclude", "excluded", "unassigned", "na", "nan", "none", ""]
    )
    if excluded_mask.any():
        excluded_samples = meta_df.index[excluded_mask].tolist()
        warnings.append(
            f"Excluded {len(excluded_samples)} sample(s) marked Exclude/Unassigned: "
            + ", ".join(map(str, excluded_samples[:8]))
            + (" ..." if len(excluded_samples) > 8 else "")
        )
        meta_df = meta_df.loc[~excluded_mask]
        counts_numeric = counts_numeric[meta_df.index]

    if meta_df.empty:
        raise DataValidationError("No samples remain after removing excluded/Unassigned samples.")

    # The R engines build their design formula around a column literally named
    # `condition` (the --design default in app/*/*.R). Normalise whichever
    # grouping column we resolved so metadata using `group`, `treatment`, etc.
    # does not fail with "object 'condition' not found".
    if condition_col != "condition":
        meta_df = meta_df.rename(columns={condition_col: "condition"})
        condition_col = "condition"

    counts_per_group = meta_df[condition_col].astype(str).value_counts()
    for group_name, count in counts_per_group.items():
        if count < 2:
            raise DataValidationError(
                f"Experimental group '{group_name}' has only {count} sample(s). "
                "Differential expression tools (DESeq2, edgeR) require at least 2 biological replicates per group to estimate dispersion."
            )

    control_group = params.get("control_group")
    treated_group = params.get("treated_group")
    if control_group is not None and treated_group is not None:
        control_group = str(control_group)
        treated_group = str(treated_group)
        if control_group == treated_group:
            raise DataValidationError("Control and treatment groups must be different.")
        for label, level in (("control", control_group), ("treatment", treated_group)):
            count = int((meta_df[condition_col].astype(str) == level).sum())
            if count < 2:
                available = ", ".join(map(str, counts_per_group.index))
                raise DataValidationError(
                    f"The {label} group '{level}' has {count} sample(s); at least 2 are required. "
                    f"Available groups: {available}."
                )

    return counts_numeric.astype(int), meta_df, warnings
