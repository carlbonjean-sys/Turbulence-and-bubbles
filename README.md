# 📖 Mode d'Emploi Total -- Code & Simulations Basilisk DNS
## Ascension d'une Bulle en Turbulence Homogène Isotrope (HIT)

**Auteur :** Carl Bonjean Fraser (ENSEEIHT / IMFT)  
**Tuteurs :** Dominique Legendre & Rémi Zamansky (IMFT)  
**Date :** Juillet 2026  

---

## 🎯 1. Guide de Prise en Main Rapide (Quick Start pour les Tuteurs)

Si vous souhaitez **consulter immédiatement les livrables finaux**, **régénérer les figures** ou **lancer une simulation de test**, voici les commandes directes :

### 📄 Accès Direct aux Livrables (PDF)
* 📘 **Rapport complet (25 pages)** : [`rapport_vers_actuelle/main_final.pdf`](rapport_vers_actuelle/main_final.pdf)
* 📺 **Présentation Beamer (9 slides)** : [`rapport_vers_actuelle/presentation_tuteurs.pdf`](rapport_vers_actuelle/presentation_tuteurs.pdf)
* 🎙️ **Script oral slide-par-slide** : [`rapport_vers_actuelle/script_oral_soutenance.pdf`](rapport_vers_actuelle/script_oral_soutenance.pdf)

### 📊 Régénération de Toutes les Figures du Rapport (1 ligne)
Depuis la racine du projet :
```bash
python scripts/plot_courbe_maitresse.py && python scripts/plot_campagne_beta.py && python scripts/plot_weber_turbulent.py && python scripts/plot_conv_turb.py && python scripts/plot_trajectoires_cote_a_cote.py
```
*Toutes les figures PNG générées sont automatiquement enregistrées dans `rapport_vers_actuelle/sections/images/`.*

---

## 📐 2. Cadre Physique & Paramètres Adimensionnels

La configuration physique modélise l'ascension d'une bulle de gaz déformable ($\rho_l/\rho_g = 850$, $\mu_l/\mu_g = 53.9$) dans un liquide soumis à une turbulence homogène et isotrope forcée en boîte triplement périodique.

| Paramètre Physique | Symbole / Formule | Valeur dans la campagne | Rôle physique |
| :--- | :--- | :--- | :--- |
| **Nombre de Bond** | $Bo = \frac{\rho_l g D^2}{\sigma}$ | **$1.0$** | Fixe la déformabilité initiale de la bulle |
| **Nombre de Galilée** | $Ga = \frac{\sqrt{gD}D}{\nu}$ | **$70.0$** | Fixe le régime d'ascension (oblique/zig-zag) |
| **Vitesse laminaire de référence** | $u_\infty$ | **$12.27$** | Vitesse terminale en liquide au repos ($level=8$) |
| **Intensité turbulente relative** | $\beta = u'/u_\infty$ | **$0.05 \to 0.65$** | **Paramètre de contrôle principal balayé** |
| **Nombre de Froude turbulent** | $Fr' = u'/\sqrt{gD}$ | **$0.077 \to 1.00$** | $Fr' = 1.53375 \times \beta$ |
| **Accélération de la pesanteur** | $g$ | **$4.0$** | Permet d'élever $\sigma = 1022.8$ pour avoir $We_t \le 1.0 < 3.0$ |
| **Taille du domaine** | $L$ | **$16D = 256$** | Évite les interactions périodiques latérales |

---

## 📂 3. Organisation Complète des Dossiers

```text
code stage IMFT/
├── src/                        # Code source C des solveurs Basilisk
│   ├── main.c                  # Solveur DNS principal
│   ├── render_lab.c            # Rendu 3D/2D hors-ligne
│   ├── mechanism.c             # Mesure de vitesse dans la coquille liquide
│   ├── static_bubble.c         # Validation des courants parasites
│   ├── compare_wall.c          # Cas témoin en domaine borné avec fond solide
│   └── diag_volume.c           # Diagnostic de conservation du volume
│
├── campagne_lvl7/              # Scripts de gestion HPC pour cluster (SLURM)
│   ├── deploy_beta_full.sh     # Lancement de la campagne complète
│   ├── deploy_ensemble.sh      # Lancement des réalisations d'ensemble (m0 à m4)
│   ├── fetch_results.sh        # Rapatriement des données .dat du cluster vers le PC local
│   ├── check_volume.sh         # Surveillance de la conservation du volume en temps réel
│   └── push.sh                 # Synchronisation du code vers le cluster
│
├── scripts/                    # Scripts Python de post-traitement & tracés
│   ├── plot_courbe_maitresse.py# Figure 4 & 5 du rapport (ralentissement vs beta et Fr')
│   ├── plot_campagne_beta.py   # Cartographie des points acquis (beta vs Froude)
│   ├── plot_weber_turbulent.py # Carte de sécurité du nombre de Weber turbulent
│   ├── plot_conv_turb.py       # Figure de convergence en maillage en turbulence
│   ├── plot_trajectoires_cote_a_cote.py # Tracé des trajectoires 3D
│   ├── plot_mechanism.py       # Analyse de l'échantillonnage vs traînée non linéaire
│   ├── postprocess.py          # Calcul des moyennes d'ensemble et incertitudes (SEM)
│   └── dataio.py               # Chargeur de fichiers de données Basilisk
│
├── simulations/                # Base de données brutes (.dat) de toutes les séries
│   ├── betasweep/              # Points d'extension (lo050, lo075, lo100, hi050, hi065)
│   └── weber_ensemble/         # Points initiaux (wt05, wt11, wt21, wt25, wt32)
│
├── rapport_vers_actuelle/      # Sources LaTeX et documents compilés
│   ├── main_final.tex          # Rapport complet (main_final.pdf)
│   ├── presentation_tuteurs.tex# Présentation Beamer (presentation_tuteurs.pdf)
│   ├── script_oral_soutenance.tex # Script oral (script_oral_soutenance.pdf)
│   ├── references.bib          # Bibliographie
│   └── sections/images/        # Dossier récepteur des figures générées
│
├── Makefile                    # Makefile de compilation pour Basilisk / MPI
└── README.md                   # Ce mode d'emploi
```

---

## 🛠️ 4. Guide de Compilation & Exécution des Solveurs C (`src/`)

### Compilation
Basilisk doit être installé sur la machine avec `$BASILISK` défini.
```bash
make            # Compile l'exécutable MPI ./bubble
make serial     # Compile la version séquentielle ./bubble_serial
make clean      # Supprime les binaire et fichiers temporaires
```

### Argumentation du Solveur `./bubble`
Le binaire `./bubble` s'exécute via la syntaxe suivante :
```bash
mpirun -np <NPROC> ./bubble MAXLEVEL MAXTIME R0 FORCED KE_TARGET INJECT [BOND] [OMECO] [GRAV] [NU] [MEMBER]
```

**Signification exacte des arguments positionnels :**
1. `MAXLEVEL` : Niveau octree max (`7` pour grille $128^3$, `8` pour $256^3$).
2. `MAXTIME` : Temps adimensionnel final de simulation (ex: `280`).
3. `R0` : Rayon initial de la bulle (défaut: `8.0`, soit diamètre $D=16$).
4. `FORCED` : `1` pour turbulence HIT forcée, `0` pour fluide au repos.
5. `KE_TARGET` : Énergie cinétique cible du forçage de Rosales & Meneveau.
6. `INJECT` : `1` pour injecter la bulle à $t=0$, `0` pour tourner en monophasique (précurseur).
7. `BOND` : Nombre de Bond $Bo$ (défaut: `1.0`).
8. `OMECO` : Paramètre d'ajustement du forçage (défaut: `0.20`).
9. `GRAV` : Pesanteur $g$ (défaut: `4.0`).
10. `NU` : Viscosité cinématique (défaut: `-1` pour utiliser $Ga=70$).
11. `MEMBER` : Indice du membre d'ensemble (`0` à `4`, décale le centre d'injection en $x,y$).

**Exemple d'exécution locale (Test rapide sur 4 cœurs) :**
```bash
mpirun -np 4 ./bubble 7 50 8.0 1 21 1 1.0 0.20 4.0 -1 0
```

---

## 💻 5. Lancement des Campagnes sur Cluster HPC (SLURM / Kairos)

Le dossier `campagne_lvl7/` contient tout le nécessaire pour déployer des calculs sur un supercalculateur SLURM :

### Étapes de déploiement pour un tuteur / chercheur :
1. **Envoyer le code source à jour vers Kairos :**
   ```bash
   bash campagne_lvl7/push.sh
   ```
2. **Lancer une campagne complète ou de nouveaux points :**
   ```bash
   bash campagne_lvl7/deploy_beta_full.sh
   ```
3. **Suivre la vitesse et la conservation du volume en direct :**
   ```bash
   bash campagne_lvl7/check_volume.sh
   ```
4. **Rapatrier automatiquement les résultats `.dat` vers le PC local :**
   ```bash
   bash campagne_lvl7/fetch_results.sh
   ```

---

## 📈 6. Traitement des Données & Ajout de Nouveaux Points (`scripts/`)

Les scripts Python lisent directement les dossiers de `simulations/`.

### Comment ajouter un nouveau point de simulation ?
1. Placez le dossier du run (contenant `frame.dat`) dans `simulations/betasweep/` ou `simulations/weber_ensemble/`.
2. Ouvrez `scripts/plot_courbe_maitresse.py` et ajoutez le nom du dossier dans le dictionnaire `RUNS`.
3. Lancez `python scripts/plot_courbe_maitresse.py`.

Le script calcule automatiquement :
- La moyenne d'ensemble sur les $n$ membres.
- L'erreur standard sur la moyenne (SEM : Standard Error of the Mean).
- La vitesse relative $\bar{u}_{\infty,\text{turb}}/u_\infty$.
- La superposition avec les lois de Spelt & Biesheuvel (1997) et Liu & Deike (2024).

---

## 📄 7. Compilation de la Documentation LaTeX (`rapport_vers_actuelle/`)

Se placer dans le dossier `rapport_vers_actuelle/` :

```bash
cd rapport_vers_actuelle

# 1. Compilation du Rapport complet (main_final.pdf)
pdflatex -interaction=nonstopmode main_final.tex
bibtex main_final
pdflatex -interaction=nonstopmode main_final.tex
pdflatex -interaction=nonstopmode main_final.tex

# 2. Compilation des Slides Beamer (presentation_tuteurs.pdf)
pdflatex -interaction=nonstopmode presentation_tuteurs.tex

# 3. Compilation du Script Oral (script_oral_soutenance.pdf)
pdflatex -interaction=nonstopmode script_oral_soutenance.tex
```

---

## ❓ 8. FAQ / Remarques Techniques Utiles pour les Tuteurs

### Q1 : Pourquoi la pesanteur vaut-elle $g=4.0$ au lieu de $g=1.0$ ?
* **Réponse** : En fixant $g=4.0$, la tension de surface nécessaire pour conserver $Bo=1$ vaut $\sigma = 1022.8$ (au lieu de $255.7$ à $g=1$). Cela réduit le nombre de Weber turbulent $We_t = \rho u'^2 D/\sigma$ d'un facteur 4 et permet de pousser l'intensité turbulente jusqu'à $\beta=0.65$ sans risquer la fragmentation artificielle de la bulle ($We_t \le 1.0 < We_c=3.0$).

### Q2 : Pourquoi ne pas utiliser le module standard `reduced.h` de Basilisk ?
* **Réponse** : Dans une boîte triplement périodique, le potentiel de gravité $g \cdot z$ n'est pas périodique en $z$. L'utilisation de `reduced.h` crée un saut de pression aux frontières verticales qui bloque la bulle (effet d'épinglage). La solution mise en œuvre dans `src/main.c` consiste à appliquer la poussée sous forme d'une force volumique moyenne périodique $\rho a = (\rho - \rho_m)g$.

### Q3 : Comment la vitesse d'ascension est-elle mesurée sans effet de grille ?
* **Réponse** : Le repère mobile (`event move_frame`) déplace le repère de calcul à la vitesse exacte de la bulle. La vitesse réelle d'ascension est enregistrée dans la variable `frame_uz` du fichier `frame.dat`.

---

*Pour toute question complémentaire, l'intégralité du code est commentée dans `src/main.c` et la chaîne de validation quantitative est détaillée dans la Section 4 de `rapport_vers_actuelle/main_final.pdf`.*
