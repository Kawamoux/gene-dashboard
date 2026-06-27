import pandas as pd


def load_expression_table(source):
    try:
        return pd.read_csv(source)
    except Exception as exc:
        raise ValueError("Le fichier CSV ne peut pas être lu.") from exc


def prepare_expression_table(dataframe, gene_column, healthy_samples, diseased_samples):
    selected_columns = [gene_column] + healthy_samples + diseased_samples
    missing_columns = [col for col in selected_columns if col not in dataframe.columns]

    if missing_columns:
        raise ValueError("Certaines colonnes sélectionnées sont absentes du CSV.")

    prepared = dataframe[selected_columns].copy()
    prepared = prepared.rename(columns={gene_column: "gene"})
    prepared = prepared[prepared["gene"].notna()]
    prepared["gene"] = prepared["gene"].astype(str).str.strip()
    prepared = prepared[prepared["gene"] != ""]

    for sample in healthy_samples + diseased_samples:
        prepared[sample] = pd.to_numeric(prepared[sample], errors="coerce")

    if prepared[healthy_samples + diseased_samples].isna().all(axis=1).any():
        prepared = prepared.dropna(subset=healthy_samples + diseased_samples, how="all")

    if prepared.empty:
        raise ValueError("Aucun gène exploitable n’a été trouvé après nettoyage.")

    values = prepared[healthy_samples + diseased_samples]
    if (values < 0).any().any():
        raise ValueError("Les valeurs d’expression doivent être positives ou nulles.")

    prepared = (
        prepared.groupby("gene", as_index=False)[healthy_samples + diseased_samples]
        .mean(numeric_only=True)
        .reset_index(drop=True)
    )

    if prepared[healthy_samples + diseased_samples].isna().any().any():
        raise ValueError("Certaines valeurs d’expression sélectionnées ne sont pas numériques.")

    return prepared
