import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


COLOR_MAP = {
    "sur-exprimé": "#d73027",
    "sous-exprimé": "#4575b4",
    "non significatif": "#9e9e9e",
}


def create_volcano_plot(results, log2fc_threshold, fdr_threshold):
    plot_data = results.copy()
    adjusted = plot_data["p_value_adjusted"].replace(0, np.nextafter(0, 1))
    plot_data["minus_log10_fdr"] = -np.log10(adjusted)
    plot_data["minus_log10_fdr"] = plot_data["minus_log10_fdr"].replace(
        [np.inf, -np.inf],
        np.nan,
    )

    fig = px.scatter(
        plot_data,
        x="log2FC",
        y="minus_log10_fdr",
        color="classification",
        color_discrete_map=COLOR_MAP,
        hover_name="gene",
        hover_data={
            "mean_healthy": ":.3f",
            "mean_disease": ":.3f",
            "p_value_adjusted": ":.3e",
            "classification": True,
        },
        labels={
            "log2FC": "log2 fold change",
            "minus_log10_fdr": "-log10(p-value ajustée)",
            "classification": "Classification",
        },
    )
    fig.add_vline(x=log2fc_threshold, line_dash="dash", line_color="#555555")
    fig.add_vline(x=-log2fc_threshold, line_dash="dash", line_color="#555555")
    if fdr_threshold > 0:
        fig.add_hline(y=-np.log10(fdr_threshold), line_dash="dash", line_color="#555555")
    fig.update_layout(
        template="plotly_white",
        legend_title_text="Classification",
        margin=dict(l=20, r=20, t=30, b=20),
    )
    return fig


def zscore_rows(matrix):
    row_means = matrix.mean(axis=1)
    row_stds = matrix.std(axis=1).replace(0, np.nan)
    zscores = matrix.sub(row_means, axis=0).div(row_stds, axis=0)
    return zscores.fillna(0)


def create_heatmap(expression_data, results, healthy_samples, diseased_samples, top_n=20):
    ordered_genes = results.head(top_n)["gene"].tolist()
    heatmap_data = expression_data[expression_data["gene"].isin(ordered_genes)].copy()
    heatmap_data["gene"] = pd.Categorical(
        heatmap_data["gene"],
        categories=ordered_genes,
        ordered=True,
    )
    heatmap_data = heatmap_data.sort_values("gene")

    samples = healthy_samples + diseased_samples
    normalized = zscore_rows(heatmap_data[samples])

    fig = go.Figure(
        data=go.Heatmap(
            z=normalized.values,
            x=samples,
            y=heatmap_data["gene"].astype(str),
            colorscale="RdBu_r",
            zmid=0,
            colorbar=dict(title="Z-score"),
        )
    )
    fig.update_layout(
        template="plotly_white",
        xaxis_title="Échantillons",
        yaxis_title="Gènes",
        margin=dict(l=20, r=20, t=30, b=20),
        height=max(450, min(900, 24 * len(heatmap_data) + 180)),
    )
    return fig


def create_pca_plot(expression_data, healthy_samples, diseased_samples):
    samples = healthy_samples + diseased_samples
    matrix = expression_data[samples].T

    if matrix.shape[0] < 2 or matrix.shape[1] < 2:
        raise ValueError("La PCA nécessite au moins deux échantillons et deux gènes.")

    scaled = StandardScaler().fit_transform(matrix)
    pca = PCA(n_components=2)
    coordinates = pca.fit_transform(scaled)

    pca_data = pd.DataFrame(
        {
            "sample": samples,
            "group": ["sain"] * len(healthy_samples) + ["malade"] * len(diseased_samples),
            "PC1": coordinates[:, 0],
            "PC2": coordinates[:, 1],
        }
    )

    labels = {
        "PC1": f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f} %)",
        "PC2": f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f} %)",
        "group": "Groupe",
    }
    fig = px.scatter(
        pca_data,
        x="PC1",
        y="PC2",
        color="group",
        text="sample",
        labels=labels,
        color_discrete_map={"sain": "#1b9e77", "malade": "#d95f02"},
    )
    fig.update_traces(textposition="top center")
    fig.update_layout(
        template="plotly_white",
        margin=dict(l=20, r=20, t=30, b=20),
    )
    return fig, pca_data
