import warnings

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def safe_ttest(healthy_values, diseased_values):
    healthy = np.asarray(healthy_values, dtype=float)
    diseased = np.asarray(diseased_values, dtype=float)
    healthy = healthy[np.isfinite(healthy)]
    diseased = diseased[np.isfinite(diseased)]

    if len(healthy) < 2 or len(diseased) < 2:
        return np.nan

    healthy_var = np.var(healthy, ddof=1)
    diseased_var = np.var(diseased, ddof=1)
    if np.isclose(healthy_var, 0.0) and np.isclose(diseased_var, 0.0):
        return 1.0 if np.isclose(np.mean(healthy), np.mean(diseased)) else 0.0

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        result = stats.ttest_ind(healthy, diseased, equal_var=False, nan_policy="omit")

    return float(result.pvalue) if np.isfinite(result.pvalue) else np.nan


def adjust_pvalues(p_values):
    p_values = np.asarray(p_values, dtype=float)
    adjusted = np.full_like(p_values, np.nan, dtype=float)
    valid = np.isfinite(p_values)

    if valid.any():
        adjusted[valid] = multipletests(p_values[valid], method="fdr_bh")[1]

    return adjusted


def classify_gene(log2fc, adjusted_pvalue, log2fc_threshold, fdr_threshold):
    if not np.isfinite(log2fc) or not np.isfinite(adjusted_pvalue):
        return "non significatif"
    if adjusted_pvalue <= fdr_threshold and log2fc >= log2fc_threshold:
        return "sur-exprimé"
    if adjusted_pvalue <= fdr_threshold and log2fc <= -log2fc_threshold:
        return "sous-exprimé"
    return "non significatif"


def differential_expression(
    expression_data,
    gene_column,
    healthy_samples,
    diseased_samples,
    log2fc_threshold=1.0,
    fdr_threshold=0.05,
    pseudocount=1e-9,
):
    healthy_matrix = expression_data[healthy_samples]
    diseased_matrix = expression_data[diseased_samples]

    mean_healthy = healthy_matrix.mean(axis=1)
    mean_diseased = diseased_matrix.mean(axis=1)
    denominator = mean_healthy + pseudocount
    numerator = mean_diseased + pseudocount
    with np.errstate(divide="ignore", invalid="ignore"):
        log2fc = np.where(
            (denominator > 0) & (numerator > 0),
            np.log2(numerator / denominator),
            np.nan,
        )

    p_values = [
        safe_ttest(healthy_matrix.iloc[index].values, diseased_matrix.iloc[index].values)
        for index in range(len(expression_data))
    ]
    adjusted = adjust_pvalues(p_values)

    results = pd.DataFrame(
        {
            "gene": expression_data[gene_column].values,
            "mean_healthy": mean_healthy.values,
            "mean_disease": mean_diseased.values,
            "log2FC": log2fc,
            "p_value": p_values,
            "p_value_adjusted": adjusted,
        }
    )
    results["classification"] = [
        classify_gene(value, pvalue, log2fc_threshold, fdr_threshold)
        for value, pvalue in zip(results["log2FC"], results["p_value_adjusted"])
    ]
    results["abs_log2FC"] = results["log2FC"].abs()
    results = results.sort_values(
        by=["p_value_adjusted", "abs_log2FC"],
        ascending=[True, False],
        na_position="last",
    ).drop(columns=["abs_log2FC"])

    return results.reset_index(drop=True)


def summarize_results(results):
    return {
        "total_genes": int(len(results)),
        "overexpressed": int((results["classification"] == "sur-exprimé").sum()),
        "underexpressed": int((results["classification"] == "sous-exprimé").sum()),
        "not_significant": int((results["classification"] == "non significatif").sum()),
    }
