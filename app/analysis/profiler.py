import pandas as pd
import numpy as np
from typing import Dict, List, Tuple
from .models import AnalysisUnit, ProfilingReport
import uuid

# Common biological column name variations
COLUMN_SYNONYMS = {
    "organism": ["organism", "species", "species_name", "taxon", "organism_id"],
    "condition": ["condition", "group", "treatment", "phenotype"],
    "time": ["time", "timepoint", "hours"],
    "tissue": ["tissue", "organ", "sample_type"]
}

def infer_column_mappings(df: pd.DataFrame) -> Dict[str, str]:
    """
    Deterministically map DataFrame columns to standard biological metadata concepts.
    Does not guess silently if ambiguous.
    """
    mapping = {}
    lower_cols = {col.lower(): col for col in df.columns}
    
    for concept, synonyms in COLUMN_SYNONYMS.items():
        found = [lower_cols[s] for s in synonyms if s in lower_cols]
        if len(found) == 1:
            mapping[concept] = found[0]
        # If multiple found, we leave it out of automatic mapping to force user confirmation
    return mapping

def detect_input_type(df: pd.DataFrame) -> Tuple[bool, List[str]]:
    """
    Detects if the input is raw counts or existing DE results.
    Returns (is_raw_counts, detected_engines)
    """
    # Check if there are float values in numeric columns which shouldn't be in raw counts
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    
    # If all numeric columns are integers (or floats that are whole numbers), it's likely raw counts
    # But wait, DESeq2 results have log2FoldChange (floats), pvalue (floats).
    
    engines = []
    lower_cols = [c.lower() for c in df.columns]
    
    # Look for signatures of existing DE results
    if any("log2foldchange" in c for c in lower_cols) or any("padj" in c for c in lower_cols):
        engines.append("DESeq2")
    
    if any("logfc" in c for c in lower_cols) and any("fdr" in c for c in lower_cols):
        engines.append("edgeR")
        
    if any("aveexpr" in c for c in lower_cols) and any("adj.p.val" in c for c in lower_cols):
        engines.append("limma")
        
    if len(engines) > 0:
        return False, engines
        
    # Check if raw counts
    is_raw = True
    for col in numeric_cols:
        if not np.all(np.mod(df[col].dropna(), 1) == 0):
            is_raw = False
            break
            
    return is_raw, engines

def profile_dataset(df: pd.DataFrame) -> ProfilingReport:
    """
    Profiles the dataset and partitions it into Analysis Units.
    """
    inferred_mapping = infer_column_mappings(df)
    is_raw_counts, detected_engines = detect_input_type(df)
    
    analysis_units = []
    warnings = []
    
    # Identify sample columns (numeric columns not identified as metadata)
    metadata_cols = list(inferred_mapping.values())
    if "gene" in [c.lower() for c in df.columns]:
        gene_col_idx = [c.lower() for c in df.columns].index("gene")
        metadata_cols.append(df.columns[gene_col_idx])
        
    sample_cols = [c for c in df.columns if c not in metadata_cols and np.issubdtype(df[c].dtype, np.number)]
    
    # Partitioning logic (Simplified for initial implementation)
    # If we have organism information, we MUST partition by organism.
    if "organism" in inferred_mapping:
        org_col = inferred_mapping["organism"]
        unique_orgs = df[org_col].dropna().unique()
        
        for org in unique_orgs:
            unit = AnalysisUnit(
                unit_id=str(uuid.uuid4()),
                organism=str(org),
                contrast="Detected automatically",
                samples=sample_cols,
                is_raw_counts=is_raw_counts,
                detected_engines=detected_engines
            )
            analysis_units.append(unit)
    else:
        # Default unit if no partitioning metadata found
        warnings.append("No organism column detected. Assuming a single organism dataset.")
        unit = AnalysisUnit(
            unit_id=str(uuid.uuid4()),
            organism="Unknown (Please specify)",
            contrast="Default Contrast",
            samples=sample_cols,
            is_raw_counts=is_raw_counts,
            detected_engines=detected_engines
        )
        analysis_units.append(unit)
        
    if not is_raw_counts and not detected_engines:
        warnings.append("Input appears to contain normalized expression values rather than raw integer counts. DESeq2 cannot be safely run on this input.")
        
    return ProfilingReport(
        total_genes=len(df),
        total_samples=len(sample_cols),
        columns_detected=list(df.columns),
        inferred_metadata_columns=inferred_mapping,
        analysis_units=analysis_units,
        warnings=warnings
    )
