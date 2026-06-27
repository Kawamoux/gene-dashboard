from pathlib import Path

import pandas as pd
import streamlit as st

from src.analysis import differential_expression, summarize_results
from src.plots import create_heatmap, create_pca_plot, create_volcano_plot
from src.preprocessing import load_expression_table, prepare_expression_table
from src.report import build_html_report


INTRO_TEXT = (
    "Cette application permet d’explorer des données d’expression génique afin de "
    "comparer deux conditions biologiques, par exemple des échantillons malades et "
    "des échantillons contrôles. Elle calcule les variations d’expression des gènes, "
    "identifie les gènes potentiellement sur-exprimés ou sous-exprimés, puis visualise "
    "les résultats sous forme de volcano plot, heatmap et PCA. Les résultats sont "
    "exploratoires et ne constituent pas une conclusion médicale."
)


def infer_defaults(columns, gene_column):
    samples = [col for col in columns if col != gene_column]
    healthy_keywords = ("healthy", "control", "ctrl", "sain", "temoin", "témoin")
    disease_keywords = ("disease", "malade", "case", "patient")

    healthy = [
        col for col in samples if any(keyword in col.lower() for keyword in healthy_keywords)
    ]
    disease = [
        col for col in samples if any(keyword in col.lower() for keyword in disease_keywords)
    ]

    if not healthy or not disease:
        midpoint = max(1, len(samples) // 2)
        healthy = samples[:midpoint]
        disease = samples[midpoint:]

    return healthy, disease


def csv_bytes(dataframe):
    return dataframe.to_csv(index=False).encode("utf-8")


def main():
    st.set_page_config(
        page_title="Dashboard d’analyse d’expression génique",
        layout="wide",
    )

    st.title("Dashboard d’analyse d’expression génique")
    st.write(INTRO_TEXT)

    uploaded_file = st.file_uploader("Importer un fichier CSV", type=["csv"])
    example_path = Path(__file__).parent / "data" / "example_expression.csv"

    if uploaded_file is None:
        dataset_name = "example_expression.csv"
        raw_data = load_expression_table(example_path)
        st.info("Aucun fichier importé : le dataset synthétique d’exemple est chargé.")
    else:
        dataset_name = uploaded_file.name
        raw_data = load_expression_table(uploaded_file)

    st.subheader("Colonnes et groupes")
    gene_column = st.selectbox("Colonne contenant les noms des gènes", raw_data.columns)
    sample_columns = [col for col in raw_data.columns if col != gene_column]
    default_healthy, default_disease = infer_defaults(raw_data.columns, gene_column)

    col_left, col_right = st.columns(2)
    with col_left:
        healthy_samples = st.multiselect(
            "Échantillons du groupe sain",
            sample_columns,
            default=[col for col in default_healthy if col in sample_columns],
        )
    with col_right:
        diseased_samples = st.multiselect(
            "Échantillons du groupe malade",
            sample_columns,
            default=[col for col in default_disease if col in sample_columns],
        )

    with st.sidebar:
        st.header("Paramètres")
        log2fc_threshold = st.number_input(
            "Seuil |log2FC|",
            min_value=0.0,
            value=1.0,
            step=0.1,
        )
        fdr_threshold = st.number_input(
            "Seuil FDR",
            min_value=0.0,
            max_value=1.0,
            value=0.05,
            step=0.01,
        )
        heatmap_size = st.selectbox("Nombre de gènes dans la heatmap", [20, 50], index=0)

    run_analysis = st.button("Lancer l’analyse", type="primary")

    if run_analysis:
        overlap = set(healthy_samples).intersection(diseased_samples)
        if overlap:
            st.error("Un échantillon ne peut pas appartenir aux deux groupes.")
            return
        if len(healthy_samples) < 2 or len(diseased_samples) < 2:
            st.error("Sélectionner au moins deux échantillons dans chaque groupe.")
            return

        try:
            expression_data = prepare_expression_table(
                raw_data,
                gene_column,
                healthy_samples,
                diseased_samples,
            )
            results = differential_expression(
                expression_data,
                gene_column="gene",
                healthy_samples=healthy_samples,
                diseased_samples=diseased_samples,
                log2fc_threshold=log2fc_threshold,
                fdr_threshold=fdr_threshold,
            )
            summary = summarize_results(results)
            volcano_fig = create_volcano_plot(results, log2fc_threshold, fdr_threshold)
            heatmap_fig = create_heatmap(
                expression_data,
                results,
                healthy_samples,
                diseased_samples,
                top_n=heatmap_size,
            )
            pca_fig, pca_data = create_pca_plot(
                expression_data,
                healthy_samples,
                diseased_samples,
            )

            st.session_state.analysis = {
                "dataset_name": dataset_name,
                "expression_data": expression_data,
                "results": results,
                "summary": summary,
                "volcano_fig": volcano_fig,
                "heatmap_fig": heatmap_fig,
                "pca_fig": pca_fig,
                "pca_data": pca_data,
                "healthy_samples": healthy_samples,
                "diseased_samples": diseased_samples,
                "log2fc_threshold": log2fc_threshold,
                "fdr_threshold": fdr_threshold,
            }
        except ValueError as exc:
            st.error(str(exc))
            return

    if "analysis" not in st.session_state:
        st.write("Importer un CSV, choisir les groupes, puis lancer l’analyse.")
        return

    analysis = st.session_state.analysis
    results = analysis["results"]
    summary = analysis["summary"]

    tabs = st.tabs(
        [
            "Résumé",
            "Tableau des gènes",
            "Volcano plot",
            "Heatmap",
            "PCA",
            "Export",
        ]
    )

    with tabs[0]:
        st.subheader("Résumé")
        metric_a, metric_b, metric_c, metric_d = st.columns(4)
        metric_a.metric("Gènes analysés", summary["total_genes"])
        metric_b.metric("Sur-exprimés", summary["overexpressed"])
        metric_c.metric("Sous-exprimés", summary["underexpressed"])
        metric_d.metric("Non significatifs", summary["not_significant"])

        st.write("Groupes comparés")
        st.write(
            pd.DataFrame(
                {
                    "groupe": ["sain", "malade"],
                    "échantillons": [
                        ", ".join(analysis["healthy_samples"]),
                        ", ".join(analysis["diseased_samples"]),
                    ],
                }
            )
        )
        st.write(
            "Les résultats indiquent des associations statistiques exploratoires, "
            "pas une cause biologique ou médicale."
        )

    with tabs[1]:
        st.subheader("Tableau des gènes différentiellement exprimés")
        st.dataframe(results, use_container_width=True)

    with tabs[2]:
        st.subheader("Volcano plot")
        st.plotly_chart(analysis["volcano_fig"], use_container_width=True)

    with tabs[3]:
        st.subheader("Heatmap")
        st.plotly_chart(analysis["heatmap_fig"], use_container_width=True)

    with tabs[4]:
        st.subheader("PCA")
        st.plotly_chart(analysis["pca_fig"], use_container_width=True)
        st.dataframe(analysis["pca_data"], use_container_width=True)

    with tabs[5]:
        st.subheader("Export")
        st.download_button(
            "Télécharger les résultats CSV",
            data=csv_bytes(results),
            file_name="gene_expression_results.csv",
            mime="text/csv",
        )

        report_html = build_html_report(
            dataset_name=analysis["dataset_name"],
            sample_count=len(analysis["healthy_samples"]) + len(analysis["diseased_samples"]),
            healthy_samples=analysis["healthy_samples"],
            diseased_samples=analysis["diseased_samples"],
            summary=summary,
            volcano_fig=analysis["volcano_fig"],
            heatmap_fig=analysis["heatmap_fig"],
            pca_fig=analysis["pca_fig"],
            log2fc_threshold=analysis["log2fc_threshold"],
            fdr_threshold=analysis["fdr_threshold"],
        )
        st.download_button(
            "Télécharger le rapport HTML",
            data=report_html.encode("utf-8"),
            file_name="gene_expression_report.html",
            mime="text/html",
        )


if __name__ == "__main__":
    main()
