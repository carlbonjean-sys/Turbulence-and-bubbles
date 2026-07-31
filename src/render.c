/**
 * render.c -- Rendu hors-ligne d'un snapshot Basilisk (coupe x-y).
 *
 * Usage : ./render <fichier_dump> <basename>
 */

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"

scalar omega[], omega2[], S2[];

int main (int argc, char ** argv)
{
  if (argc < 3) {
    fprintf (ferr, "usage: %s <dump> <basename>\n", argv[0]);
    return 1;
  }
  if (!restore (file = argv[1])) {
    fprintf (ferr, "render: impossible de lire '%s'\n", argv[1]);
    return 1;
  }
  char * base = argv[2];
  char name[256];

  foreach() {
    double uxx = (u.x[1]    - u.x[-1]   )/(2.*Delta);
    double uxy = (u.x[0,1]  - u.x[0,-1] )/(2.*Delta);
    double uxz = (u.x[0,0,1]- u.x[0,0,-1])/(2.*Delta);
    double uyx = (u.y[1]    - u.y[-1]   )/(2.*Delta);
    double uyy = (u.y[0,1]  - u.y[0,-1] )/(2.*Delta);
    double uyz = (u.y[0,0,1]- u.y[0,0,-1])/(2.*Delta);
    double uzx = (u.z[1]    - u.z[-1]   )/(2.*Delta);
    double uzy = (u.z[0,1]  - u.z[0,-1] )/(2.*Delta);
    double uzz = (u.z[0,0,1]- u.z[0,0,-1])/(2.*Delta);

    double wx = uzy - uyz, wy = uxz - uzx, wz = uyx - uxy;
    omega2[] = wx*wx + wy*wy + wz*wz;
    omega[]  = sqrt(omega2[]);

    double sxy = 0.5*(uxy + uyx), sxz = 0.5*(uxz + uzx), syz = 0.5*(uyz + uzy);
    S2[] = uxx*uxx + uyy*uyy + uzz*uzz + 2.*(sxy*sxy + sxz*sxz + syz*syz);
  }

  double sz = 0., sv = 0.;
  foreach (reduction(+:sz) reduction(+:sv)) {
    double w = dv()*(1. - f[]);
    sz += w*z; sv += w;
  }
  double zc = (sv > 1e-10 ? sz/sv : L0/2.);

  view (width = 1024, height = 1024, fov = 22, camera = "front",
        tx = -0.5, ty = -0.5, bg = {1,1,1}, samples = 4);

  const char * field[]  = {"u.z", "omega", "omega2", "S2"};
  const char * suffix[] = {"vitesse", "vorticite", "enstrophie", "cisaillement"};

  clear();
  cells   (n = {0,0,1}, alpha = zc);
  draw_vof ("f", lw = 3.);
  box();
  sprintf (name, "%s_bulle_maillage.png", base); save (name);

  for (int k = 0; k < 4; k++) {
    clear();
    squares (field[k], n = {0,0,1}, alpha = zc, linear = true);
    draw_vof ("f", lw = 3.);
    box();
    sprintf (name, "%s_%s.png", base, suffix[k]); save (name);
  }

  return 0;
}
