/**
 * render_lab.c -- Rendu 2D hors-ligne d'un snapshot (coupe verticale x-z).
 *
 * Usage : ./render_lab <fichier_dump> <basename> [frame_uz]
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
  size (120.);                 // WIDTH
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

  /* Centre de masse du gaz par moyenne circulaire (boîte périodique) */
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

  /* Vue caméra quasi-orthographique */
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
