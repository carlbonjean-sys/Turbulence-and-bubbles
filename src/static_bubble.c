/**
 * static_bubble.c -- TEST V1 : COURANTS PARASITES (spurious currents).
 *
 * Une bulle SPHERIQUE au repos dans un liquide au repos, SANS gravite et SANS
 * forcage, est une solution exacte des equations : la tension de surface est
 * exactement equilibree par le saut de pression de Laplace, dp = 2 sigma / R,
 * et la vitesse doit rester RIGOUREUSEMENT NULLE pour toujours.
 *
 * Tout ce qu'on mesure comme vitesse est donc de l'erreur numerique pure : les
 * "courants parasites" du schema CSF (tension de surface). C'est le test de
 * validation le plus propre qui soit, parce que la reponse exacte est connue.
 *
 * On les compare a la vitesse d'ascension mesuree en production (u_inf = 12.27) :
 * si u_parasite / u_inf est de l'ordre du pourcent, la mesure de vitesse est
 * saine ; si c'est 10 %, tout le reste est a relativiser.
 *
 * Parametres RIGOUREUSEMENT identiques a la production (main.c) : meme domaine,
 * meme R0, meme RHOR/MUR, meme viscosite, et sigma calcule avec la MEME formule
 * (Bo=1, GRAV=4) => sigma = 1022.8. Seule la gravite est retiree.
 *
 * Sorties -> spurious.dat :  t  umax  urms  Ca_max  dp_num  dp_theo
 *   umax   = max |u| dans tout le domaine          (doit rester ~0)
 *   urms   = sqrt(<|u|^2>)
 *   Ca_max = mu1 * umax / sigma                    (nombre capillaire parasite)
 *   dp_num = <p>_gaz - <p>_liquide  mesure         (doit valoir 2 sigma / R)
 *
 * Args : maxlevel MAXTIME     (defauts 7, 5)
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
#define GRAV  4.0        // sert UNIQUEMENT a calculer sigma comme en production
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
    periodic (right);              // meme domaine triplement periodique

  rho1 = 1.0;  rho2 = rho1/RHOR;
  mu1  = 0.01*sq(WIDTH/(2.0*pi))/2.0;
  mu2  = mu1/MUR;

  N = 1 << (maxlevel - 3);         // base + adaptation, comme le mode bulle
  TOLERANCE = 1e-3;
  CFL = 0.5;
  /* PAS de champ d'acceleration : ni gravite ni forcage. C'est tout l'interet. */
  run();
}

event init (i = 0)
{
  f.sigma = (rho1 - rho2)*GRAV*sq(2.0*R0)/BOND;      // = 1022.8, comme en prod
  double xc = L0/2., yc = L0/2., zc = L0/2.;
  refine (sq(x-xc) + sq(y-yc) + sq(z-zc) < sq(1.30*R0) &&
          sq(x-xc) + sq(y-yc) + sq(z-zc) > sq(0.70*R0) &&
          level < maxlevel);
  fraction (f, sq(x-xc) + sq(y-yc) + sq(z-zc) - sq(R0));
  foreach()
    foreach_dimension()
      u.x[] = 0.;                                     // liquide ET gaz au repos
  if (pid() == 0)
    fprintf (ferr, "[V1] maxlevel=%d  sigma=%g  mu1=%g  dp_theo=2sigma/R=%g  "
             "D/dx=%g\n", maxlevel, f.sigma, mu1, 2.*f.sigma/R0,
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
    /* pression moyenne loin de l'interface, de part et d'autre */
    if (f[] < 0.01) { pg += dv()*p[]; wg += dv(); }     // coeur gazeux
    if (f[] > 0.99) { pl += dv()*p[]; wl += dv(); }     // liquide
  }
  double urms   = sqrt(u2/vol);
  double dp_num = (wg > 0. && wl > 0. ? pg/wg - pl/wl : 0.);

  if (pid() == 0) {
    static FILE * fp = NULL;
    if (!fp) {
      fp = fopen("spurious.dat", "a");
      fprintf(fp, "# t umax urms Ca_max dp_num dp_theo\n"
                  "# bulle spherique au repos, sans gravite : u DOIT rester 0.\n");
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
