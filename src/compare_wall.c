/**
 * compare_wall.c -- CAS DE COMPARAISON : bulle montant dans un domaine BORNE
 * (fond solide en bas, haut ouvert / surface libre) avec GRAVITE REDUITE (reduced.h).
 *
 * But : comparer le profil de vitesse d'ascension a notre cas de production
 * (domaine triplement periodique + force volumique + repere mobile). reduced.h
 * est VALIDE ici car le domaine n'est PAS periodique selon la gravite (z borne) :
 * c'est son usage standard, bien equilibre. On garde EXACTEMENT les memes
 * parametres physiques (RHOR, MUR, R0, GRAV, Bo=1, meme nu) et le meme maxlevel 7.
 *
 * Question : la vitesse est-elle plus "uniforme" (moins d'a-coups) qu'en periodique ?
 *
 * Args : maxlevel MAXTIME   (defauts 7, 12)
 */
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"                 // gravite reduite : OK ici (z NON periodique)
#include "navier-stokes/perfs.h"
#include "maxruntime.h"

#define RHOR  850.0
#define MUR   100.0
#define WIDTH 180.0            // domaine PLUS HAUT (vs 120 en prod) : ~12 u.t. de
                              // montee (z=24 -> ~172) => atteint ET soutient la
                              // vitesse terminale (~1 periode d'oscillation).
                              // Contrepartie a lvl7 : D/dx = 16/(180/128) ~ 11
                              // cellules (vs 17 en prod) -> resolution qualitative.
#define R0    8.0
#define GRAV  4.0
#define BOND  1.0
int    maxlevel = 7;
double MAXTIME  = 18.0;

/* x,y periodiques ; z BORNE : fond solide (back, z=0) + haut ouvert (front, z=L0).
   Fond = non-glissement ; haut = sortie libre (p=0), approximation de surface libre. */
u.n[back]  = dirichlet(0.);
u.t[back]  = dirichlet(0.);
u.r[back]  = dirichlet(0.);
u.n[front] = neumann(0.);
p[front]   = dirichlet(0.);
pf[front]  = dirichlet(0.);

int main (int argc, char * argv[])
{
  maxruntime (&argc, argv);
  if (argc > 1) maxlevel = atoi(argv[1]);
  if (argc > 2) MAXTIME  = atof(argv[2]);

  size (WIDTH); origin (0., 0., 0.);
  periodic (right);                  // x periodique
  periodic (top);                    // y periodique  (z reste borne : back/front)

  rho1 = 1.0; rho2 = rho1/RHOR;
  /* 🚨 BUG CORRIGE le 2026-07-22 : la ligne etait
         mu1 = 0.01*sq(WIDTH/(2.0*pi))/2.0;   // "MEME viscosite que main.c"
     ... mais WIDTH vaut 180 ICI et 120 dans main.c, donc mu1 sortait a 4.104
     au lieu de 1.824, soit Ga = 31 au lieu de 70 -- SOUS le seuil d'instabilite
     de trajectoire (Ga_c ~ 44). Ce cas de comparaison aurait donc tourne dans un
     regime physique DIFFERENT : la bulle n'aurait pas devie, et on en aurait
     conclu que c'est le cas periodique qui est fautif. Le commentaire affirmait
     l'inverse de ce que faisait le code.
     -> on fige la constante de main.c (120), independamment du WIDTH local. */
  mu1  = 0.01*sq(120.0/(2.0*pi))/2.0;   // = 1.82378, identique a main.c
  mu2  = mu1/MUR;
  f.sigma = (rho1 - rho2)*GRAV*sq(2.0*R0)/BOND;   // Bo = BOND (= 1)

  G.z = -GRAV; Z.z = 0.;             // gravite reduite (bien equilibree, domaine borne)

  init_grid (1 << maxlevel);
  TOLERANCE = 1e-3; CFL = 0.5;
  run();
}

event init (i = 0)
{
  if (restore (file = "dump")) return 0;
  double xc = L0/2., yc = L0/2., zc = 24.;    // bulle a 3 rayons du fond
  refine (sq(x-xc)+sq(y-yc)+sq(z-zc) < sq(1.30*R0) &&
          sq(x-xc)+sq(y-yc)+sq(z-zc) > sq(0.70*R0) && level < maxlevel);
  fraction (f, sq(x-xc)+sq(y-yc)+sq(z-zc) - sq(R0));   // f=1 liquide, f=0 gaz
  foreach()
    foreach_dimension()
      u.x[] = 0.;                    // liquide AU REPOS (laminaire)
}

event adapt (i++)
{
  adapt_wavelet ((scalar *){f, u}, (double[]){1e-3, 3e-2, 3e-2, 3e-2}, maxlevel, 5);
}

/* Suivi de la bulle (gaz, f<0.5) : centre de masse + vitesse (ascension uz + lateral). */
event logbubble (i += 5)
{
  double w = 0., xc = 0., yc = 0., zc = 0., uz = 0., ux = 0., uy = 0.;
  foreach (reduction(+:w) reduction(+:xc) reduction(+:yc) reduction(+:zc)
           reduction(+:uz) reduction(+:ux) reduction(+:uy)) {
    double g = (1. - f[])*dv();
    w += g; xc += g*x; yc += g*y; zc += g*z;
    uz += g*u.z[]; ux += g*u.x[]; uy += g*u.y[];
  }
  if (pid() == 0 && w > 1e-10) {
    static FILE * fp = NULL;
    if (!fp) { fp = fopen("bubble.dat","a"); fprintf(fp,"# t xc yc zc uz ux uy\n"); }
    fprintf (fp, "%g %g %g %g %g %g %g\n", t, xc/w, yc/w, zc/w, uz/w, ux/w, uy/w);
    fflush (fp);
  }
}

event snapshots (t += 0.5; t <= MAXTIME)
{
  char name[80]; sprintf (name, "snapshot-%06.2f", t); dump (file = name);
}
event checkpoint (i += 200) { dump (file = "dump"); }
event end (t = MAXTIME)     { dump (file = "end"); }
