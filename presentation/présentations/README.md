# Présentation d'avancement de stage (Beamer)

Diaporama LaTeX/Beamer pour la soutenance d'avancement devant les tuteurs
(D. Legendre, R. Zamansky).

```
présentations/
├── presentation.tex     # le diaporama (Beamer, 16:9)
├── README.md            # ce fichier
└── images/              # images utilisées (copiées + renommées proprement)
```

> **Note sur les images.** Dans `../sections/images/`, plusieurs fichiers ont un
> nom qui **ne correspond pas** à leur contenu (ex. `grace.png` est en réalité la
> bulle 2D en zigzag, `image_recuperee_1.png` est la caractérisation de la
> turbulence, etc.). Les images utilisées ici ont donc été **copiées et renommées**
> dans `images/` pour que le diaporama soit autonome et sans ambiguïté.

## Comment compiler

### Option A — Overleaf (le plus simple, recommandé)

1. Aller sur https://www.overleaf.com → **New Project** → **Upload Project**.
2. Compresser le dossier `présentations/` (avec son sous-dossier `images/`) en
   `.zip` et l'uploader. (Ou créer un projet vide puis glisser-déposer
   `presentation.tex` **et** le dossier `images/`.)
3. Dans Overleaf : menu **Menu** (en haut à gauche) → **Compiler =
   `pdfLaTeX`**, **TeX Live ≥ 2021**.
4. Cliquer **Recompile**. (Overleaf fait automatiquement les 2 passes
   nécessaires à la table des matières et aux bandeaux TikZ.)

Le projet est **autonome** : tout ce qu'il faut est dans `présentations/`, pas
besoin d'uploader le reste du rapport.

### Option B — En local (MiKTeX / TeX Live)

Depuis le dossier `présentations/` :

```bash
pdflatex presentation.tex
pdflatex presentation.tex      # 2e passe : sommaire + bandeaux
```

ou, plus simple, avec `latexmk` :

```bash
latexmk -pdf presentation.tex
```

> Windows : installer **MiKTeX** (https://miktex.org) puis utiliser TeXworks
> (moteur `pdfLaTeX`) ou la ligne de commande ci-dessus. MiKTeX installe les
> paquets manquants à la volée (beamer, babel-french, colortbl, tikz…).

## Modifier le diaporama

- Le moteur est **pdfLaTeX** (pas besoin de XeLaTeX/LuaLaTeX).
- Aucune police exotique : `lmodern` + thème Beamer `Madrid` recoloré en bleu IMFT.
- Le bandeau rouge « Résultats issus de l'ANCIEN code » est posé par la macro
  `\oldcode` (à appeler juste après `\begin{frame}{...}` sur les diapos
  concernées).
- Pour remplacer une figure obsolète par un résultat du **nouveau** code :
  déposer la nouvelle image dans `images/` (même nom) et retirer le `\oldcode`
  de la diapo correspondante.
