/**
 * diag_volume.c -- Recalcul a posteriori du volume, du centre de masse
 *                  et de la vitesse de la bulle sur snapshot Basilisk.
 *
 * Usage : ./diag_volume <snapshot> <etiquette>
 */
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"

#define WIDTH 120.0

int main (int argc, char ** argv)
{
  if (argc < 3) {
    fprintf (ferr, "usage: %s <snapshot> <etiquette>\n", argv[0]);
    return 1;
  }

  size (WIDTH);
  origin (0., 0., 0.);
  foreach_dimension()
    periodic (right);

  if (!restore (file = argv[1])) {
    fprintf (ferr, "diag_volume: lecture impossible de '%s'\n", argv[1]);
    return 1;
  }

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
