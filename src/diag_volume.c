/**
 * diag_volume.c -- RECALCUL a posteriori du volume et de la vitesse de la bulle
 *                  sur un snapshot, avec les sommes CORRECTES.
 *
 * A quoi ca sert
 * --------------
 * Jusqu'au 2026-07-22, `event bubble` de main.c ne sommait que sur les mailles
 * taguees (f < 0.5) et jetait donc la moitie exterieure de la bande d'interface,
 * qui contient pourtant du gaz. Resultat : volume sous-estime de ~4.5 % et
 * vitesse verticale biaisee de ~+0.35 en repere mobile. main.c est corrige, mais
 * TOUS LES RUNS DEJA FAITS ont un bubble.dat entache. Ce programme les rattrape
 * sans avoir a relancer quoi que ce soit : il relit les snapshots et recalcule.
 *
 * 🚨 PIEGE MORTEL, LA RAISON D'ETRE DES 4 LIGNES CI-DESSOUS
 * ---------------------------------------------------------
 * `restore()` NE RESTAURE NI LA TAILLE DU DOMAINE NI LES CONDITIONS AUX LIMITES.
 * Un post-traitement qui se contente de restore() travaille donc sur une boite
 * NON PERIODIQUE. Consequence constatee le 2026-07-22 : tag() ne reconnectait
 * plus la bulle a travers la frontiere des qu'elle chevauchait un bord, le
 * "volume de la bulle mere" s'effondrait a 50 %, et j'en ai conclu a tort a une
 * fausse fragmentation systematique du solveur. Il n'en etait rien : c'etait ce
 * post-traitement qui etait faux. => TOUJOURS re-declarer size(), origin() et
 * periodic() a l'identique de main().
 *
 * Sortie (une ligne par snapshot, sur stdout) :
 *   tag  V_complet  V_tronque  w_complet  w_tronque  xc  yc  zc  n_regions
 *     V_complet = Sum (1-f) dV sur TOUTES les mailles          <- la bonne valeur
 *     V_tronque = idem restreint aux mailles taguees f<0.5     <- l'ancienne
 *     w_*       = vitesse verticale ponderee correspondante
 *     xc,yc,zc  = centre de masse en MOYENNE CIRCULAIRE (periodicite)
 *
 * usage : ./diag_volume <snapshot> <etiquette>
 */
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"

#define WIDTH 120.0            // identique a main.c

int main (int argc, char ** argv)
{
  if (argc < 3) {
    fprintf (ferr, "usage: %s <snapshot> <etiquette>\n", argv[0]);
    return 1;
  }

  /* OBLIGATOIRE AVANT restore() -- cf. l'avertissement en tete. */
  size (WIDTH);
  origin (0., 0., 0.);
  foreach_dimension()
    periodic (right);

  if (!restore (file = argv[1])) {
    fprintf (ferr, "diag_volume: lecture impossible de '%s'\n", argv[1]);
    return 1;
  }

  /* --- somme COMPLETE + centre en moyenne circulaire --- */
  double kk = 2.*pi/L0;
  double Vfull = 0., Wfull = 0.;
  coord cs = {0,0,0}, ss = {0,0,0};
  foreach (reduction(+:Vfull) reduction(+:Wfull)
           reduction(+:cs) reduction(+:ss)) {
    double w = dv()*(1. - clamp(f[], 0., 1.));
    if (w > 0.) {
      Vfull += w;
      Wfull += w*u.z[];
      coord p = {x, y, z};
      foreach_dimension() {
        cs.x += w*cos(kk*p.x);
        ss.x += w*sin(kk*p.x);
      }
    }
  }
  if (Vfull <= 0.) { fprintf(fout, "%s nan nan nan nan nan nan nan 0\n", argv[2]); return 0; }
  coord com;
  foreach_dimension() {
    com.x = atan2(ss.x, cs.x)/kk;
    if (com.x < 0.) com.x += L0;
  }

  /* --- somme TRONQUEE (l'ancienne definition), pour quantifier l'ecart --- */
  scalar m[];
  foreach()
    m[] = (f[] < 0.5);
  int n = tag(m);
  double Vtag = 0., Wtag = 0.;
  if (n >= 1) {
    double * vol = calloc(n, sizeof(double));
    foreach_leaf()
      if (m[] > 0) vol[(int)m[]-1] += dv()*(1. - f[]);
#if _MPI
    MPI_Allreduce(MPI_IN_PLACE, vol, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    int jm = 0;
    for (int j = 1; j < n; j++) if (vol[j] > vol[jm]) jm = j;
    Vtag = vol[jm]; free(vol);
    int mother = jm + 1;
    foreach_leaf()
      if ((int)m[] == mother) Wtag += dv()*(1. - f[])*u.z[];
    if (Vtag > 0.) Wtag /= Vtag;
  }

  fprintf (fout, "%s %.4f %.4f %.6f %.6f %.4f %.4f %.4f %d\n",
           argv[2], Vfull, Vtag, Wfull/Vfull, Wtag, com.x, com.y, com.z, n);
  return 0;
}
