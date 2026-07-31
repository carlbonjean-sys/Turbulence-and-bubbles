# Films — Précurseur turbulent (HIT)

Films de la **boîte turbulente du précurseur** (turbulence homogène isotrope forcée),
avant injection de la bulle. Champ level-8-AMR (`simulations/precurseur8vortices`),
turbulence développée et statistiquement stationnaire.

**Paramètres de la turbulence filmée :**
- Reynolds de Taylor : **Re_λ ≈ 45**
- Fluctuation rms : u' = 10,7 · Dissipation : ε = 52,8 · Kolmogorov : η = 0,58
- Résolution : k_max·η = 3,9 (toutes les échelles résolues, DNS)
- **12 images, t = 19 → 30** (unités code), ~4 img/s

Les coupes 2D sont prises dans le plan horizontal (z-normal) à mi-domaine.

| Fichier | Ce que montre le film |
|---|---|
| **`precurseur_vitesse.mp4`** | Champ de **vitesse verticale u_z** (coupe z). Grandes structures énergétiques de l'écoulement, montées/descentes de fluide. |
| **`precurseur_vorticite.mp4`** | **Norme de la vorticité \|ω\| = \|rot(u)\|** (coupe z). Les filaments rouges = zones de fort cisaillement/rotation ; visualise le cœur des tourbillons. |
| **`precurseur_enstrophie.mp4`** | **Enstrophie \|ω\|²** (coupe z). Comme la vorticité mais accentue les structures les plus intenses (∝ dissipation locale). |
| **`precurseur_vortex3d.mp4`** | **Isosurface 3D de \|ω\|** (seuil 1,5·ω_rms), colorée par u_z (bleu = descendant, rouge = montant). Les **tubes tourbillonnaires en « vermicelles »** caractéristiques de la turbulence 3D. Le rendu 3D lisse confirme que les structures sont bien résolues par le maillage. |

**Génération :** rendus offline via `render3d` (WSL/OSMesa) sur les snapshots
`snapshot-019..030.00`, assemblés avec ffmpeg. Voir la section "Offline rendering"
de `README.md`.
