from pydantic import BaseModel, Field
from typing import List, Dict, Optional, Any

class AnalysisUnit(BaseModel):
    """
    Represents a logically comparable biological analysis partition.
    Never mix different organisms, tissues, or independent contrasts.
    """
    unit_id: str
    organism: str
    tissue: Optional[str] = None
    treatment: Optional[str] = None
    time_point: Optional[str] = None
    contrast: str
    samples: List[str]
    is_raw_counts: bool
    detected_engines: List[str] = Field(default_factory=list) # e.g. ["DESeq2", "edgeR"] if pre-computed

class ProfilingReport(BaseModel):
    """
    Result of profiling an uploaded dataset.
    """
    total_genes: int
    total_samples: int
    columns_detected: List[str]
    inferred_metadata_columns: Dict[str, str] # e.g. {"organism": "species", "condition": "treatment"}
    analysis_units: List[AnalysisUnit]
    warnings: List[str]
