# 📖 Mode d'Emploi Total -- Code & Simulations Basilisk DNS
## Ascension d'une Bulle en Turbulence Homogène Isotrope (HIT)

**Auteur :** Carl Bonjean Fraser (ENSEEIHT / IMFT)  
**Tuteurs :** Dominique Legendre & Rémi Zamansky (IMFT)  
**Date :** Juillet 2026  

---

## 🎯 1. Guide de Prise en Main Rapide

### 📄 Livrables (PDF)
* 📘 **Rapport complet (25 pages)** : [`rapport_vers_actuelle/main_final.pdf`](rapport_vers_actuelle/main_final.pdf)
* 📺 **Présentation Beamer (9 slides)** : [`rapport_vers_actuelle/presentation_tuteurs.pdf`](rapport_vers_actuelle/presentation_tuteurs.pdf)
* 🎙️ **Script oral slide-par-slide** : [`rapport_vers_actuelle/script_oral_soutenance.pdf`](rapport_vers_actuelle/script_oral_soutenance.pdf)

### 📊 Régénération des Figures du Rapport
Depuis la racine du projet :
```bash
python scripts/plot_courbe_maitresse.py && python scripts/plot_campagne_beta.py && python scripts/plot_weber_turbulent.py && python scripts/plot_conv_turb.py && python scripts/plot_trajectoires_cote_a_cote.py
```

---

## 📐 2. Cadre Physique & Paramètres Adimensionnels

La configuration physique modélise l'ascension d'une bulle de gaz déformable ($\rho_l/\rho_g = 850$, $\mu_l/\mu_g = 53.9$) dans un liquide soumis à une turbulence homogène et isotrope forcée en boîte triplement périodique.

| Paramètre Physique | Symbole / Formule | Valeur dans la campagne | Rôle physique |
| :--- | :--- | :--- | :--- |
| **Nombre de Bond** | $Bo = \frac{\rho_l g D^2}{\sigma}$ | **$1.0$** | Déformabilité initiale de la bulle |
| **Nombre de Galilée** | $Ga = \frac{\sqrt{gD}D}{\nu}$ | **$70.0$** | Régime d'ascension |
| **Vitesse laminaire de référence** | $u_\infty$ | **$12.27$** | Vitesse terminale en liquide au repos |
| **Intensité turbulente relative** | $\beta = u'/u_\infty$ | **$0.05 \to 0.65$** | **Paramètre de contrôle principal ($\beta$-sweep)** |
| **Nombre de Froude turbulent** | $Fr' = u'/\sqrt{gD}$ | **$0.077 \to 1.00$** | $Fr' = 1.53375 \times \beta$ |
| **Accélération de la pesanteur** | $g$ | **$4.0$** | Pesanteur |
| **Taille du domaine** | $L$ | **$16D = 256$** | Boîte triplement périodique |

---

## 📂 3. Organisation des Dossiers

```text
code stage IMFT/
├── src/                        # Code source C des solveurs Basilisk
│   ├── main.c                  # Solveur DNS principal
│   ├── render_lab.c            # Rendu 3D/2D hors-ligne
│   ├── mechanism.c             # Échantillonnage préférentiel & traînée
│   ├── static_bubble.c         # Validation des courants parasites
│   ├── compare_wall.c          # Cas témoin domaine borné
│   └── diag_volume.c           # Diagnostic de conservation du volume
│
├── campagne_lvl7/              # Scripts de gestion HPC pour cluster (SLURM)
│   ├── deploy_beta_full.sh     # Lancement de la campagne complète
│   ├── deploy_ensemble.sh      # Lancement des réalisations d'ensemble (m0 à m4)
│   ├── fetch_results.sh        # Rapatriement des données .dat du cluster
│   ├── check_volume.sh         # Surveillance du volume en temps réel
│   └── push.sh                 # Synchronisation vers le cluster
│
├── scripts/                    # Scripts Python de post-traitement & tracés
│   ├── plot_courbe_maitresse.py# Ralentissement vs beta et Fr'
│   ├── plot_campagne_beta.py   # Cartographie des points acquis
│   ├── plot_weber_turbulent.py # Nombre de Weber turbulent
│   ├── plot_conv_turb.py       # Convergence en maillage sous turbulence
│   ├── plot_trajectoires_cote_a_cote.py # Trajectoires 3D
│   ├── plot_mechanism.py       # Analyse du mécanisme de ralentissement
│   ├── postprocess.py          # Calcul des moyennes d'ensemble et SEM
│   └── dataio.py               # Chargeur de données Basilisk
│
├── simulations/                # Données brutes de simulation (.dat)
│   ├── betasweep/              # Balayage en intensité turbulente (lo050..hi065, wt05..wt32)
│   ├── convturb/               # Étude de convergence en maillage (lvl7, lvl8)
│   ├── lam7_frame/             # Run laminaire de référence
│   ├── precursor/              # Précurseur turbulent monophasique
│   └── validation/             # Validation courants parasites et convergence laminaire
│
├── rapport_vers_actuelle/      # Livrables finaux
│   ├── main_final.pdf          # Rapport complet
│   ├── presentation_tuteurs.pdf# Présentation Beamer
│   └── script_oral_soutenance.pdf # Script oral
│
├── Makefile                    # Makefile de compilation Basilisk
└── README.md                   # Ce document
```

---

## 🛠️ 4. Compilation & Exécution des Solveurs C (`src/`)

```bash
make            # Compile le binaire MPI ./bubble
make serial     # Compile la version séquentielle ./bubble_serial
make clean      # Nettoie les binaires
```

### Syntaxe d'exécution :
```bash
mpirun -np <NPROC> ./bubble MAXLEVEL MAXTIME R0 FORCED KE_TARGET INJECT [BOND] [OMECO] [GRAV] [NU] [MEMBER]
```

**Exemple d'exécution rapide (4 cœurs) :**
```bash
mpirun -np 4 ./bubble 7 50 8.0 1 21 1 1.0 0.20 4.0 -1 0
```

---

## 📈 5. Traitement des Données (`scripts/`)

Les scripts Python lisent la base de données standardisée `simulations/betasweep/`.

Pour ajouter une nouvelle série de simulations :
1. Placez le dossier du run (contenant `frame.dat`) dans `simulations/betasweep/`.
2. Ajoutez le tag dans le dictionnaire `POINTS` de `scripts/plot_courbe_maitresse.py`.
3. Lancez `python scripts/plot_courbe_maitresse.py`.
