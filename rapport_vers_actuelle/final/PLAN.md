# Plan du rapport final de stage

DNS d'une bulle de gaz (Bo=1) montant dans une turbulence homogène isotrope :
la turbulence accélère-t-elle ou ralentit-elle la remontée ?

Passage du **journal de bord chronologique** actuel à un **rapport final scientifique**.
Le fond existe déjà (sec9, ~43 p., audité) ; il s'agit de le réorganiser et de rédiger
les sections manquantes (introduction, état de l'art, méthodes, conclusion).

## Structure cible

| # | Section | Source | Rédacteur |
|---|---|---|---|
| 1 | Introduction, contexte, état de l'art | neuf (+ refs.bib) | Sonnet 5 |
| 2 | Cadre physique et adimensionnel | sec4 étendu | Sonnet 5 |
| 3 | Méthodes numériques | consolidé (sec6/sec9) | Sonnet 5 |
| 4 | Validation du modèle | sec9 (réorganisé) | assemblage |
| 5 | Résultats : effet de la turbulence | sec9 (réorganisé) | assemblage |
| 6 | Mécanisme du ralentissement | sec9 (déjà rédigé) | assemblage |
| 7 | Conclusion et perspectives | neuf | Sonnet 5 |
| A | Annexe : exploration 2D | sec5 condensé | assemblage |

## Fiche de faits (source unique de vérité pour toutes les sections)

- **Objet** : bulle unique gaz/liquide, Bo=1, Ga=70, dans THI triplement périodique.
  Solveur Basilisk (VOF conservant la qdm), octree AMR, MPI, cluster CALMIP/Kairos.
- **Rapports** : ρ_l/ρ_g = 850, μ_l/μ_g = 100. Domaine L0=120, R0=8, D=16 (code units).
- **Adimensionnels** : Bo=1 ; Ga=70 (> Ga_c≈44 ⇒ trajectoire oblique) ; We_t ; Re_λ ;
  Fr' = u'/√(gD) = 1.53·β ; β = u'/u_∞.
- **Référence laminaire** : u_∞ = 12.27 (mesuré, validé en maillage lvl7↔lvl8 à +0.32 %).
  Re_b = 108, C_D = 0.566 (entre bulle propre et sphère solide), χ = 1.38. Accélération
  initiale mesurée 8.17 vs masse ajoutée 7.97 (+2.4 %, sans ajustement). St = 0.080.
- **Résultat central (n=5)** : β=0.15→−6 %, 0.22→−12 %, 0.31→−31 %, 0.33→−32 %,
  0.38→−38 %. Loi 1−cβ². Rejoint la courbe maîtresse de Liu, Farsoiya, Perrard & Deike
  (JFM 2024) : parabole petit-Fr' puis 0.37/Fr'. On attrape la transition vers Fr'≈0.45.
  ⚠ barres d'erreur ±3–6 % à β≳0.3 (u'/u_∞≈0.5) — la magnitude haute reste incertaine.
- **Artefact** : reduced.h ⇒ potentiel φ=[ρ]g(z−Z) linéaire, non périodique ⇒ épingle la
  bulle à z=L0−R0=112. Démasqué par le run laminaire. Correctif : force volumique
  périodique ρa=(ρ−ρ_m)g + repère mobile + garde-fou de qdm.
- **Volume** : conservé à 99.05–99.88 % (les « 93-97 % » étaient un artefact de mesure du
  diagnostic, corrigé). d_eq = 15.99 vs 16.
- **Convergence maillage (laminaire)** : u_∞ = 10.93 / 12.27 / 12.31 (lvl6/7/8), +0.32 %
  lvl7→8. δ/Δx = 0.87 / 1.64 / 3.28. Le paramètre pilote = couche limite, pas D/Δx.
- **Convergence turbulente** : trajectoires divergent (chaos, |Δh|→2.6 diamètres),
  statistiques convergent (Re_λ identique) ⇒ la convergence porte sur les statistiques.
- **Courants parasites** : u_max = 0.0093 lvl7 = 0.08 % de u_∞ ; Laplace +1.8 %.
- **Biais de sillage** (domaine périodique) : −0.8 %, borné (le sillage sature). 2 méthodes.
- **Résolution turbulente** : k_max·η ≥ 3.92 partout (seuil 1.5).
- **Mécanisme** : traînée non linéaire (échantillonnage préférentiel ≤12 %, ~0 en moyenne ;
  fluctuations de glissement 8→50 %). Accord Liu-Deike pour bulle super-Kolmogorov.

## Clés bibliographiques (references.bib, \cite{} uniquement, JAMAIS \citep/\citet)

legendre2025 (revue RMP, tuteur), liu2024 (LE papier du sujet), riviere2021 (We_c=3),
ruth2021, spelt1997, poorte2002, canolozano2016 (Ga_c), mei1994, loisy2017, bassenne2016
(forçage contrôlé), popinet2009 (Basilisk).

## Figures disponibles (scripts/figures/, copiées dans sections/images/)

courbe_maitresse, weber_vs_Fr_litterature, mechanism, convergence_maillage,
convergence_turbulente, trajectoires_cote_a_cote, trajectoires_3d, trajectoires_xz,
epinglage, diagnostic_bubble, ascension_laminaire, biais_sillage, campagne_beta,
pourquoi_meme_vitesse, precurseur_caracterisation, pourquoi_meme_vitesse.

## Règles de rédaction (toutes sections)

- Français, style rapport scientifique, sobre. Pas de « je » sauf conclusion.
- N'inventer AUCUN nombre : n'utiliser que ceux de la fiche de faits. Incertain ⇒ `% TODO`.
- `\cite{cle}` uniquement (le préambule n'a pas natbib).
- Le préambule fournit les packages ; écrire du LaTeX qui compile seul (pas de \begin{document}).
