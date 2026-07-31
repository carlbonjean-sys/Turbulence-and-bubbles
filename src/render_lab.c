/**
 * render_lab.c -- Rendu offline d'un snapshot : COUPE VERTICALE (x-z).
 *
 *   usage : ./render_lab <fichier_dump> <basename> [frame_uz]
 *
 * POURQUOI CE FICHIER EXISTE (2026-07-21) -----------------------------------
 * `render.c` coupe en `n = {0,0,1}` (plan HORIZONTAL x-y) et regarde avec
 * `camera = "front"` (= quaternion identite => ecran = plan x-y, l'axe z RENTRE
 * dans l'ecran). Consequence : dans les films produits, **la montee de la bulle
 * est geometriquement invisible** -- on la voit de dessus. Il ne reste a l'ecran
 * que la derive laterale (|u_h| oscillant 0.14 -> 2.57, periode 16 u.t.), ce qui
 * donne l'impression que « la bulle rampe ». Et comme le champ affiche est u.z
 * dans un plan horizontal, le liquide y vaut ~ -frame_uz UNIFORMEMENT (repere
 * mobile) : un aplat de couleur constant => « le fond ne bouge pas ».
 * Ici on coupe donc en `n = {0,1,0}` (plan VERTICAL x-z) avec `camera = "top"`
 * (rotation -pi/2 autour de x => ecran y = z monde => la gravite est vers le bas
 * a l'ecran, comme il se doit).
 *
 * REPERE MOBILE ------------------------------------------------------------
 * Le champ stocke est RELATIF au repere qui suit la bulle : u_labo = u + frame_uz.
 * - sans argument [frame_uz] : on affiche u.z tel quel = **vue relative**, le
 *   liquide defile vers le bas devant la bulle (vue "soufflerie"). C'est la plus
 *   lisible : on voit la bulle grimper dans le fluide et le sillage se former.
 * - avec [frame_uz] : on affiche aussi uz_lab = u.z + frame_uz = **vue labo**,
 *   ou le liquide loin de la bulle est au repos. A relire dans frame.dat (col 1)
 *   a l'instant du snapshot.
 */
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"

scalar omega[], uz_lab[];

int main (int argc, char ** argv)
{
  if (argc < 3) {
    fprintf (ferr, "usage: %s <dump> <basename> [frame_uz]\n", argv[0]);
    return 1;
  }
  /* OBLIGATOIRE AVANT restore() : celui-ci ne restaure NI la taille du domaine NI
     les conditions aux limites. Sans ces trois lignes, les mailles fantomes des
     bords sont en symetrie au lieu d'etre periodiques -> vorticite FAUSSE sur la
     couche de bord, et tout tag()/connexite coupe au franchissement d'un bord.
     (Bug commis le 2026-07-22 dans un diagnostic autonome : il avait fait conclure
     a tort a une "fausse fragmentation". Re-declarer a l'identique de main.c.) */
  size (120.);                 // = WIDTH de main.c
  origin (0., 0., 0.);
  foreach_dimension()
    periodic (right);

  if (!restore (file = argv[1])) {
    fprintf (ferr, "render_lab: impossible de lire '%s'\n", argv[1]);
    return 1;
  }
  char * base = argv[2];
  double frame_uz = (argc > 3 ? atof(argv[3]) : 0.);
  char name[256];

  foreach() {
    double wx = (u.z[0,1,0] - u.z[0,-1,0] - (u.y[0,0,1] - u.y[0,0,-1]))/(2.*Delta);
    double wy = (u.x[0,0,1] - u.x[0,0,-1] - (u.z[1]     - u.z[-1]    ))/(2.*Delta);
    double wz = (u.y[1]     - u.y[-1]     - (u.x[0,1,0] - u.x[0,-1,0]))/(2.*Delta);
    omega[]   = sqrt(wx*wx + wy*wy + wz*wz);
    uz_lab[]  = u.z[] + frame_uz;               // vitesse verticale LABO
  }

  /* Centre de masse du gaz -> plan de coupe vertical y = yc.
     MOYENNE CIRCULAIRE OBLIGATOIRE (comme main.c event bubble) : le domaine est
     periodique, et quand la bulle est A CHEVAL sur un bord (moitie a y~119,
     moitie a y~1) une moyenne arithmetique donne ~60, soit le cote OPPOSE du
     domaine -> on coupait dans du liquide vide pendant ~11 u.t. par traversee.
     Bug constate le 2026-07-21 sur snapshot-228 : arithmetique 62.9 vs vraie 119.7. */
  double kk = 2.*pi/L0;
  double cy = 0., sy = 0., cz = 0., sz = 0., sv = 0.;
  foreach (reduction(+:cy) reduction(+:sy) reduction(+:cz) reduction(+:sz)
           reduction(+:sv)) {
    double w = dv()*(1. - f[]);
    sv += w;
    cy += w*cos(kk*y); sy += w*sin(kk*y);
    cz += w*cos(kk*z); sz += w*sin(kk*z);
  }
  double yc = L0/2., zc = L0/2.;
  if (sv > 1e-10) {
    yc = atan2(sy, cy)/kk; if (yc < 0.) yc += L0;
    zc = atan2(sz, cz)/kk; if (zc < 0.) zc += L0;
  }
  fprintf (ferr, "render_lab: coupe verticale y=%g (bulle a z=%g), frame_uz=%g\n",
           yc, zc, frame_uz);

  /* CADRAGE -- calibre empiriquement le 2026-07-21, NE PAS "corriger" de tete :
     seul (camera="bottom", tx=-0.5, ty=-0.5) donne z VERS LE HAUT avec x correct.
     Verifie sur snapshot-200.00 (bulle en xc=71.5, zc=73.2) : la bulle tombe a
     +11.7 unites a droite du centre et a z=73.0 -> les deux axes sont bons.
     Les autres combinaisons echouent :
       camera="top",    ty=+0.5 -> z INVERSE (le sillage apparait AU-DESSUS) ;
       camera="top",    ty=-0.5 -> domaine hors champ (image quasi vide) ;
       camera="bottom", ty=+0.5 -> domaine hors champ.
     NB : les etiquettes d'axes de box() sont vues "de dos" a cet angle et
     paraissent en miroir -- c'est le rendu des labels, PAS les donnees.

     CAMERA QUASI-ORTHOGRAPHIQUE (tz/near/far + fov etroit) -- INDISPENSABLE :
     la camera regarde LE LONG DE y, et le plan de coupe est en y = yc, qui suit
     la bulle. En perspective, ce plan s'eloigne/se rapproche de l'oeil quand la
     bulle derive en y => l'echelle apparente CHANGE : mesure sur snapshot-200,
     une coupe a y=10 fait 427 px de large, la meme a y=110 en fait 41 (-90 %).
     Sur le film ca se voit comme un dezoom progressif, remis a zero a chaque
     franchissement du bord periodique. En reculant la camera a 30 (unites de
     domaine) avec fov 2.23 deg, la variation tombe a -2.7 %. */
  view (width = 1024, height = 1024, fov = 2.23, camera = "bottom",
        tx = -0.5, ty = -0.5, tz = -30., near = 29., far = 31.,
        bg = {1,1,1}, samples = 4);

  const char * field[]  = {"u.z", "uz_lab", "omega"};
  const char * suffix[] = {"vitesse_rel", "vitesse_lab", "vorticite"};
  int nf = (frame_uz != 0. ? 3 : 3);            // uz_lab == u.z si frame_uz=0

  /* 1. Bulle + maillage, coupe VERTICALE. */
  clear();
  cells    (n = {0,1,0}, alpha = yc);
  draw_vof ("f", lw = 3.);
  box();
  sprintf (name, "%s_bulle_maillage.png", base); save (name);

  /* 2-4. Champs, coupe verticale, interface en surimpression. */
  for (int k = 0; k < nf; k++) {
    clear();
    squares  (field[k], n = {0,1,0}, alpha = yc, linear = true);
    draw_vof ("f", lw = 3.);
    box();
    sprintf (name, "%s_%s.png", base, suffix[k]); save (name);
  }

  return 0;
}
