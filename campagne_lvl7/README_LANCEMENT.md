# Campagne lvl7 — lancement (préparée le 2026-07-08, télétravail)

Tout est prêt et testé localement (compilation MPI + smoke test serial OK en WSL,
`[PARAMS]` et `levels.dat` validés). **Rien n'a été soumis** : c'est toi qui
lances demain au labo.

## Checklist demain matin (4 commandes)

```bash
# 1. Depuis WSL (labo) — pousser main.c + deploy.sh sur Kairos :
cd "/mnt/c/Users/carlb/OneDrive/Desktop/code stage IMFT" && bash campagne_lvl7/push.sh

# 2. Sur Kairos — dry-run (compile, prépare les 6 dossiers, montre les job.slurm) :
ssh kairoslogin
cd /work/p0910/bonjeanf/bulles/campagne_lvl7 && bash deploy.sh

# 3. Vérifier un job généré (partition/account corrects ?) puis soumettre :
cat /work/p0910/bonjeanf/bulles/lvl7/job.slurm
LAUNCH=1 bash deploy.sh

# 4. Suivi :
squeue -u $USER
tail -f /work/p0910/bonjeanf/bulles/lvl7/slurm-*.out
```

Si le dry-run ne trouve pas d'ancien `job.slurm` pour recopier l'en-tête SLURM :
`PARTITION=<ta_partition> ACCOUNT=<ton_compte> LAUNCH=1 bash deploy.sh`.

## Les 3 nombres sans dimension (Buckingham)

Le problème a **exactement 3 groupes indépendants** = **(Bo, We_t, Re_λ)**, pilotés
par 3 arguments (le 4ᵉ groupe Ga/Mo/Re_b est alors *dérivé*) :

| Arg | Bouton | Effet | Précurseur |
|---|---|---|---|
| `BOND` (7) | **Bo** | Bo = BOND exactement | même `end_we1` |
| `GRAV` (9) | **We_t à Bo fixe** | We_t ∝ 1/GRAV | même `end_we1` |
| `NU` (10) | **Re_λ** | ν↓ → Re_λ↑ (≤~55 à lvl7) | **nouveau** précurseur |

g et σ n'agissent qu'à l'interface → le précurseur monophasique est insensible à
`BOND`/`GRAV` : **les axes Bo et We_t se font tous avec `end_we1`**. Seul l'axe Re_λ
(via `NU`) demande un nouveau précurseur (non lancé ici).

## Matrice des runs (7 jobs × 8 cœurs = 56 ≤ garde-fou 64)

Ordre des args : `maxlevel MAXTIME R0 FORCED KE INJECT BOND OMECO GRAV [NU]`.

| Run | Args `./bubble` | Durée | Répond à | Attendu |
|---|---|---|---|---|
| `lvl7` (extension) | `7 226 8 1 15 1 1.0 0.20 1.0` | t=126→226 | P2 (≥10 τ), P5 | Bo=1, We_t≈0.93 ; auto-resume ; `dump_t126.bak` sauvé avant |
| `lvl7_we2` | `7 226 8 1 15 1 1.0 0.20 0.5` | t=76→226 | P3, P5 | **Bo=1**, We_t≈1.9, σ=127.9, dt ×1.41 |
| `lvl7_we3` | `7 226 8 1 15 1 1.0 0.20 0.333` | t=76→226 | P3, P5 | **Bo=1**, We_t≈2.8 — **fragmentation possible = résultat (We_c!), pas un bug** |
| `sens_ome02` | `7 86 8 1 15 1 1.0 0.20 1.0` | t=76→86 | P1 + non-rég. | doit reproduire le début du `lvl7` d'origine* |
| `sens_ome01` | `7 86 8 1 15 1 1.0 0.10 1.0` | t=76→86 | P1 | plus de cellules (seuil = tolérance : plus bas = plus fin) |
| `sens_ome008` | `7 86 8 1 15 1 1.0 0.08 1.0` | t=76→86 | P1 | encore plus de cellules ; borne « fin » du sweep |
| `sens_ome04` | `7 86 8 1 15 1 1.0 0.40 1.0` | t=76→86 | P1 | moins de cellules |

**Tous à Bo=1, Re_λ≈25 (même `end_we1`) ; seul We_t (via GRAV) ou OMECO varie.**
C'est la correction du sweep BOND initial qui couplait Bo et We_t.

\*Non-régression nuancée : le nouveau binaire a `MINLEVEL=5` (base 2⁵ + plancher
adapt) alors que le lvl7 d'origine tournait à base 2⁴/plancher 1. Le fond turbulent
aura donc un peu plus de cellules ; l'interface (à maxlevel 7) et donc `u_z`/volume
ne doivent PAS bouger. Si `u_z(t)` se superpose → binaire validé.

## Validations à faire une fois les jobs partis (5 min)

1. **`[PARAMS]` dans chaque `slurm-*.out`** : `Bo`, `GRAV`, `sigma` corrects ?
   (we2 → `Bo=1 GRAV=0.5 sigma=127.849`, we3 → `Bo=1 GRAV=0.333 sigma=85.23`,
   les autres → `Bo=1 GRAV=1 sigma=255.698`). Tous doivent afficher `MINLEVEL=5`.
2. **Non-régression** : quand `sens_ome02` est fini (≈1 h), comparer son
   `bubble.dat` au `lvl7` d'origine sur t∈[76,86] — superposition attendue à
   l'epsilon près (même binaire logique, mêmes défauts BOND=1/OMECO=0.2).
   Si écart visible → STOP, ne pas exploiter bo2/bo3, investiguer.
3. `levels.dat` apparaît bien dans chaque dossier (nouvel event).

## Rapatriement + post-traitement (le soir même pour les sens, J+2 pour le reste)

```bash
# WSL :
bash campagne_lvl7/fetch_results.sh
# Windows ou WSL :
python scripts/postprocess.py --root . --levels 7 8
python scripts/plot_levels.py simulations/sens_ome01 simulations/sens_ome02 \
    simulations/sens_ome04 --out scripts/figures/amr_sensibilite.png
```

Figures P6/P7 déjà générées avec les données actuelles (dans `scripts/figures/`) :
`cumulative_mean.png`, `t0_robustness.png`, `divergence_lvl7_lvl8.png`,
`convergence.png` + `convergence_table.csv`.
**Résultat clé déjà acquis (P6)** : en excluant le transitoire (t0=injection+8),
la pente de régression donne **U_b ≈ 0** (lvl7 : −0.03, lvl8 : −0.05) sur
t∈[84,126] — le 1.08 « fenêtre complète » était dominé par la montée initiale.
À I≈0.8, pas de montée nette sur 3.5 τ → c'est exactement l'argument pour les
runs longs (P2) et l'ensemble (P5).

## Notes techniques

- **TOLERANCE = 1e-3 partout** : c'est le réglage du run lvl7 validé (le 0.1
  n'était nécessaire qu'au lvl8, multigrid raide). Campagne 100 % lvl7 ⇒ cohérent.
- **Reprise après walltime** : resoumettre le même `job.slurm` tel quel
  (auto-resume depuis `dump`) ; **REPASSER LE MÊME `BOND`** — `f.sigma` est
  refixé par `event init` à chaque démarrage à partir de l'argument.
- Wallclocks : sens ≈ ≲1 h ; bo2/bo3 : ~6600/4400 pas (dt_σ×√Bo) ; extension :
  ~4400 pas. Les walltimes SLURM (48 h) sont très larges — mesurer le vrai
  coût dans `perfs` (colonne wallclock) pour le chiffrage P5.
- Tout vit sur `/work` (jamais `/users`, quota 50 GB).
