from html import escape

import plotly.io as pio


def figure_html(fig, include_plotlyjs):
    return pio.to_html(
        fig,
        full_html=False,
        include_plotlyjs=include_plotlyjs,
        config={"displayModeBar": False},
    )


def build_html_report(
    dataset_name,
    sample_count,
    healthy_samples,
    diseased_samples,
    summary,
    volcano_fig,
    heatmap_fig,
    pca_fig,
    log2fc_threshold,
    fdr_threshold,
):
    healthy_text = escape(", ".join(healthy_samples))
    diseased_text = escape(", ".join(diseased_samples))
    dataset_text = escape(dataset_name)

    volcano_html = figure_html(volcano_fig, include_plotlyjs=True)
    heatmap_html = figure_html(heatmap_fig, include_plotlyjs=False)
    pca_html = figure_html(pca_fig, include_plotlyjs=False)

    return f"""<!doctype html>
<html lang="fr">
<head>
  <meta charset="utf-8">
  <title>Rapport d’analyse d’expression génique</title>
  <style>
    body {{
      font-family: Arial, sans-serif;
      line-height: 1.5;
      margin: 32px;
      color: #222;
    }}
    h1, h2 {{
      color: #16324f;
    }}
    .summary {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
      gap: 12px;
      margin: 20px 0;
    }}
    .metric {{
      border: 1px solid #d9e2ec;
      border-radius: 8px;
      padding: 14px;
      background: #f8fafc;
    }}
    .metric strong {{
      display: block;
      font-size: 1.6rem;
    }}
    section {{
      margin-top: 28px;
    }}
  </style>
</head>
<body>
  <h1>Rapport d’analyse d’expression génique</h1>
  <p>Dataset : <strong>{dataset_text}</strong></p>
  <p>Nombre d’échantillons comparés : <strong>{sample_count}</strong></p>
  <p>Groupe sain : {healthy_text}</p>
  <p>Groupe malade : {diseased_text}</p>
  <p>Critères : |log2FC| ≥ {log2fc_threshold:.2f} et FDR ≤ {fdr_threshold:.3f}</p>

  <div class="summary">
    <div class="metric">Gènes analysés<strong>{summary["total_genes"]}</strong></div>
    <div class="metric">Gènes sur-exprimés<strong>{summary["overexpressed"]}</strong></div>
    <div class="metric">Gènes sous-exprimés<strong>{summary["underexpressed"]}</strong></div>
    <div class="metric">Non significatifs<strong>{summary["not_significant"]}</strong></div>
  </div>

  <section>
    <h2>Volcano plot</h2>
    {volcano_html}
  </section>

  <section>
    <h2>Heatmap</h2>
    {heatmap_html}
  </section>

  <section>
    <h2>PCA</h2>
    {pca_html}
  </section>

  <section>
    <h2>Conclusion prudente</h2>
    <p>
      Cette analyse met en évidence des gènes dont l’expression diffère entre les groupes
      comparés selon les seuils choisis. Ces résultats sont exploratoires, dépendent de la
      qualité du dataset et ne permettent pas d’établir une cause biologique ou une conclusion
      médicale.
    </p>
  </section>
</body>
</html>"""
