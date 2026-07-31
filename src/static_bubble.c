/**
 * static_bubble.c -- Test de validation des courants parasites (bulle au repos).
 *
 * Usage : ./static_bubble [maxlevel] [MAXTIME]
 */
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "maxruntime.h"

#define RHOR  850.0
#define MUR   100.0
#define WIDTH 120.0
#define R0    8.0
#define GRAV  4.0
#define BOND  1.0

int    maxlevel = 7;
double MAXTIME  = 5.0;

int main (int argc, char * argv[])
{
  maxruntime (&argc, argv);
  if (argc > 1) maxlevel = atoi(argv[1]);
  if (argc > 2) MAXTIME  = atof(argv[2]);

  size (WIDTH);
  origin (0., 0., 0.);
  foreach_dimension()
    periodic (right);

  rho1 = 1.0;  rho2 = rho1/RHOR;
  mu1  = 0.01*sq(WIDTH/(2.0*pi))/2.0;
  mu2  = mu1/MUR;

  N = 1 << (maxlevel - 3);
  TOLERANCE = 1e-3;
  CFL = 0.5;
  run();
}

event init (i = 0)
{
  f.sigma = (rho1 - rho2)*GRAV*sq(2.0*R0)/BOND;
  double xc = L0/2., yc = L0/2., zc = L0/2.;
  refine (sq(x-xc) + sq(y-yc) + sq(z-zc) < sq(1.30*R0) &&
          sq(x-xc) + sq(y-yc) + sq(z-zc) > sq(0.70*R0) &&
          level < maxlevel);
  fraction (f, sq(x-xc) + sq(y-yc) + sq(z-zc) - sq(R0));
  foreach()
    foreach_dimension()
      u.x[] = 0.;
  if (pid() == 0)
    fprintf (ferr, "maxlevel=%d  sigma=%g  mu1=%g  dp_theo=%g  D/dx=%g\n",
             maxlevel, f.sigma, mu1, 2.*f.sigma/R0,
             2.*R0/(WIDTH/(1 << maxlevel)));
}

event spurious (i += 5)
{
  double umax = 0., u2 = 0., vol = 0.;
  double pg = 0., wg = 0., pl = 0., wl = 0.;
  foreach (reduction(max:umax) reduction(+:u2) reduction(+:vol)
           reduction(+:pg) reduction(+:wg) reduction(+:pl) reduction(+:wl)) {
    double m = sqrt(sq(u.x[]) + sq(u.y[]) + sq(u.z[]));
    if (m > umax) umax = m;
    u2  += dv()*sq(m);
    vol += dv();
    if (f[] < 0.01) { pg += dv()*p[]; wg += dv(); }
    if (f[] > 0.99) { pl += dv()*p[]; wl += dv(); }
  }
  double urms   = sqrt(u2/vol);
  double dp_num = (wg > 0. && wl > 0. ? pg/wg - pl/wl : 0.);

  if (pid() == 0) {
    static FILE * fp = NULL;
    if (!fp) {
      fp = fopen("spurious.dat", "a");
      fprintf(fp, "# t umax urms Ca_max dp_num dp_theo\n");
    }
    fprintf(fp, "%g %g %g %g %g %g\n", t, umax, urms,
            mu1*umax/f.sigma, dp_num, 2.*f.sigma/R0);
    fflush(fp);
  }
}

event adapt (i++)
{
  adapt_wavelet ((scalar *){f}, (double[]){1e-3}, maxlevel, 5);
}

event end (t = MAXTIME)
{
  dump (file = "end");
}
