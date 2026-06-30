# Dashboard d’analyse d’expression génique

J’ai développé un dashboard Python d’analyse d’expression génique permettant de comparer des échantillons malades et sains à partir de données publiques. Le projet combine programmation, statistiques et biologie moléculaire afin d’identifier et visualiser les gènes différentiellement exprimés dans une maladie donnée.

## But du projet

Ce projet propose une application Streamlit pour explorer des données d’expression génique et comparer deux groupes biologiques, par exemple des échantillons malades et des témoins sains. L’application calcule des indicateurs statistiques simples, classe les gènes selon leur variation d’expression et génère des visualisations interprétables.

Les résultats sont exploratoires. Ils ne constituent pas une conclusion médicale et ne permettent pas d’identifier une cause de maladie.

## Expression génique

L’expression génique correspond au niveau d’activité d’un gène dans un échantillon biologique. Selon la technologie utilisée, elle peut être estimée par des comptages, des intensités ou des valeurs normalisées. Un gène très exprimé produit généralement davantage d’ARN mesurable qu’un gène peu exprimé.

## Gène différentiellement exprimé

Un gène est dit différentiellement exprimé lorsque son niveau moyen d’expression varie entre deux conditions. Ici, la comparaison repose sur :

- la moyenne d’expression dans le groupe sain ;
- la moyenne d’expression dans le groupe malade ;
- le log2 fold change ;
- une p-value issue d’un test t de Welch ;
- une p-value ajustée par correction FDR.

Par défaut, un gène est classé comme sur-exprimé ou sous-exprimé si `|log2FC| >= 1` et si la p-value ajustée est `<= 0.05`.

## Installation

Créer un environnement Python, puis installer les dépendances :

```bash
pip install -r requirements.txt
```

Lancer l’application :

```bash
streamlit run app.py
```

## Format du CSV

Le fichier CSV doit contenir :

- une colonne avec les noms des gènes ;
- une colonne par échantillon ;
- des valeurs numériques positives ou nulles.

Exemple :

```csv
gene,Healthy_1,Healthy_2,Healthy_3,Disease_1,Disease_2,Disease_3
Gene_001,12.1,11.8,12.4,31.5,30.7,32.1
Gene_002,7.2,7.5,7.1,3.1,3.3,3.0
```

Un fichier synthétique d’exemple est fourni dans le dossier `data`.

## Utilisation

1. Importer un fichier CSV ou utiliser le fichier d’exemple.
2. Choisir la colonne contenant les noms des gènes.
3. Sélectionner les échantillons du groupe sain.
4. Sélectionner les échantillons du groupe malade.
5. Cliquer sur **Lancer l’analyse**.
6. Consulter les onglets de résumé, tableau, graphiques et export.

## Interprétation des graphiques

Le volcano plot affiche chaque gène selon son log2 fold change et sa significativité statistique. Les points à droite correspondent aux gènes plus exprimés dans le groupe malade, tandis que les points à gauche correspondent aux gènes moins exprimés dans le groupe malade. Les couleurs mettent en évidence les classifications selon les seuils choisis.

La heatmap montre les gènes les plus significatifs et compare les profils d’expression entre les échantillons. Les valeurs sont normalisées par gène afin de faciliter la comparaison visuelle.

La PCA résume la variabilité globale des échantillons en deux axes principaux. Une séparation visuelle entre groupes peut suggérer des différences globales d’expression, mais elle doit être interprétée avec prudence.

## Export

L’application permet de télécharger :

- le tableau complet des résultats au format CSV ;
- un rapport HTML contenant le nom du dataset, les groupes comparés, le nombre d’échantillons, le nombre de gènes sur-exprimés et sous-exprimés, les graphiques et une conclusion prudente.

## Limites

Cette analyse utilise un test statistique simple et ne remplace pas une analyse bioinformatique complète. Elle ne modélise pas les effets de lots, l’âge, le sexe, les traitements, la qualité des échantillons ou d’autres facteurs de confusion. Les résultats dépendent fortement de la qualité des données, du nombre d’échantillons et de la normalisation réalisée en amont.

Le dashboard sert à explorer et visualiser des tendances. Toute interprétation biologique doit être validée par des méthodes adaptées, des contrôles expérimentaux et une expertise du domaine.
