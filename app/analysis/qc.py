import pandas as pd
import numpy as np
from sklearn.decomposition import PCA

def compute_qc_metrics(counts_df: pd.DataFrame, meta_df: pd.DataFrame) -> dict:
    """
    Computes quality control metrics for RNA-seq count data.
    Matches the structure expected by the frontend QC dashboard.
    """
    # counts_df: index = genes, columns = samples
    sample_names = list(counts_df.columns)
    
    # Try to find group column
    possible_group_cols = [c for c in meta_df.columns if c.lower() not in ['sample', 'sample_id', 'id']]
    group_col = possible_group_cols[0] if possible_group_cols else meta_df.columns[0]
    sample_groups = meta_df[group_col].tolist()
    
    total_genes = len(counts_df)
    
    # Library sizes (sum per sample)
    library_sizes = counts_df.sum(axis=0).tolist()
    avg_lib_size = float(np.mean(library_sizes)) if library_sizes else 0.0
    
    # Detection rates (% genes > 0)
    detection_rates = ((counts_df > 0).mean(axis=0) * 100).tolist()
    avg_det_rate = float(np.mean(detection_rates)) if detection_rates else 0.0
    
    # Filter genes with low expression for PCA and correlation (CPM > 1 in at least 2 samples)
    lib_size_array = np.array(library_sizes)
    # Avoid division by zero
    lib_size_array[lib_size_array == 0] = 1
    cpm = (counts_df * 1e6) / lib_size_array
    keep_mask = (cpm >= 1.0).sum(axis=1) >= 2
    
    genes_retained = int(keep_mask.sum())
    filtered_counts = counts_df[keep_mask]
    
    # Correlation Matrix (Spearman)
    if not filtered_counts.empty and len(sample_names) > 1:
        corr_df = filtered_counts.corr(method="spearman")
        corr_matrix = corr_df.values.tolist()
        mean_corrs = corr_df.mean(axis=0).tolist()
    else:
        corr_matrix = [[1.0 for _ in sample_names] for _ in sample_names]
        mean_corrs = [1.0 for _ in sample_names]
        
    # PCA
    if not filtered_counts.empty and len(sample_names) > 1:
        log_cpm = np.log2(cpm[keep_mask] + 1)
        pca = PCA(n_components=min(2, len(sample_names)))
        components = pca.fit_transform(log_cpm.T)
        variance_explained = (pca.explained_variance_ratio_ * 100).tolist()
        
        # Ensure we have 2 components to return even if n_samples == 2
        pc1 = components[:, 0].tolist() if components.shape[1] > 0 else [0] * len(sample_names)
        pc2 = components[:, 1].tolist() if components.shape[1] > 1 else [0] * len(sample_names)
        
        pca_result = {
            "components": [pc1, pc2],
            "varianceExplained": variance_explained
        }
    else:
        pca_result = {
            "components": [[0] * len(sample_names), [0] * len(sample_names)],
            "varianceExplained": [0, 0]
        }
        
    # Outlier detection (simple heuristic)
    outlier_status = []
    outlier_count = 0
    
    median_size = np.median(library_sizes) if library_sizes else 0
    median_rate = np.median(detection_rates) if detection_rates else 0
    
    for i in range(len(sample_names)):
        size = library_sizes[i]
        rate = detection_rates[i]
        avg_corr = mean_corrs[i]
        
        reasons = []
        is_outlier = False
        is_warning = False
        
        if avg_corr < 0.80:
            is_outlier = True
            reasons.append(f"Low correlation (r={avg_corr:.3f})")
        elif avg_corr < 0.85:
            is_warning = True
            reasons.append(f"Reduced correlation (r={avg_corr:.3f})")
            
        if size < 0.1 * median_size:
            is_outlier = True
            reasons.append(f"Severely low library size ({(size/1e6):.2f}M reads)")
        elif size < 0.3 * median_size:
            is_warning = True
            reasons.append(f"Low library size ({(size/1e6):.2f}M reads)")
            
        if rate < 50:
            is_outlier = True
            reasons.append(f"Severely low detection rate ({rate:.1f}%)")
        elif rate < 70:
            is_warning = True
            reasons.append(f"Reduced detection rate ({rate:.1f}%)")
            
        status = "Fail" if is_outlier else ("Warning" if is_warning else "Pass")
        reason = "; ".join(reasons) if reasons else "Normal (Passed all QC criteria)"
        
        if is_outlier:
            outlier_count += 1
            
        outlier_status.append({
            "index": i,
            "avgCorr": float(avg_corr),
            "corrZ": 0,
            "sizeZ": 0,
            "rateZ": 0,
            "isOutlier": is_outlier,
            "isWarning": is_warning,
            "reason": reason,
            "status": status
        })
        
    return {
        "sampleNames": sample_names,
        "sampleGroups": sample_groups,
        "librarySizes": library_sizes,
        "detectionRates": detection_rates,
        "corrMatrix": corr_matrix,
        "outlierStatus": outlier_status,
        "pcaResult": pca_result,
        "avgLibSize": avg_lib_size,
        "avgDetRate": avg_det_rate,
        "outlierCount": outlier_count,
        "genesRetained": genes_retained,
        "totalGenes": total_genes
    }
