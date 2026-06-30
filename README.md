# Dashboard d’analyse d’expression génique

Ce dépôt contient un dashboard Python développé avec Streamlit pour explorer des données d’expression génique et comparer deux groupes biologiques, par exemple des échantillons malades et des témoins sains.

Le projet complet se trouve dans le dossier [`gene_dashboard`](gene_dashboard/).

## Lancer l’application

```bash
cd gene_dashboard
pip install -r requirements.txt
streamlit run app.py
```

## Contenu du projet

- `gene_dashboard/app.py` : interface Streamlit.
- `gene_dashboard/src/` : fonctions de nettoyage, analyse statistique, graphiques et rapport HTML.
- `gene_dashboard/data/example_expression.csv` : dataset synthétique d’exemple.
- `gene_dashboard/README.md` : documentation détaillée du projet.

Les résultats produits par l’application sont exploratoires et ne constituent pas une conclusion médicale.
