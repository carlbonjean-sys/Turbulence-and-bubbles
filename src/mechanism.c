/**
 * mechanism.c -- POURQUOI la turbulence ralentit la bulle ?
 *                Echantillonnage preferentiel vs trainee non lineaire.
 *
 * Deux mecanismes candidats, SIGNATURES OPPOSEES et mesurables :
 *  (A) ECHANTILLONNAGE PREFERENTIEL : la bulle passe plus de temps dans du
 *      fluide qui DESCEND. Signature : <u_z fluide vu par la bulle> - <u_z
 *      domaine> NEGATIF.
 *  (B) TRAINEE NON LINEAIRE (Liu, Farsoiya, Perrard & Deike 2024) : fluide
 *      ambiant a vitesse moyenne nulle, mais C_D(Re) non lineaire -> les
 *      fluctuations de glissement augmentent la trainee moyenne. Signature :
 *      <u_z fluide> ~ 0 mais forte VARIANCE du glissement.
 *
 * ASTUCE REPERE MOBILE : u.z stocke = u_z^lab - frame_uz PARTOUT. La difference
 * (proche bulle) - (moyenne domaine) est donc INDEPENDANTE de frame_uz (il
 * s'annule) -> calculable sur le snapshot seul. Verifie : uz_domaine mesure
 * vaut bien ~ -frame_uz (le liquide "defile" a -12.2 en repere mobile).
 *
 * DIFFICULTE : la bulle induit son propre ecoulement (sillage sous elle,
 * potentiel autour). PARADE : cone SUPERIEUR (evite sillage et equateur) a
 * plusieurs rayons ; si <u_z> -> valeur negative quand r grandit => (A), si
 * -> 0 => (B). Coquille pleine et hemisphere en controle.
 *
 * SORTIE (1 ligne/snapshot) :
 *   label z_c u_bulle uz_domaine  4x[cone] 4x[hemi] 4x[coquille]  Vb
 *   bandes r/D = [1,1.5] [1.5,2] [2,2.5] [2.5,3] ; toutes STOCKEES (repere mobile).
 *
 * PIEGES qcc appris ici (2026-07-26) :
 *  - periodic() AVANT restore() (ne restaure ni domaine ni CL) ;
 *  - NE PAS nommer une variable 'sf'/'wf' : 'sf' est un scalaire Basilisk
 *    (two-phase-generic.h: scalar sf, fraction lissee) ;
 *  - pas de reduction sur 'coord' hors d'un event -> reductions scalaires ;
 *  - f[]/u.z[] uniquement DANS une expression (clamp(f[],..), w*u.z[]),
 *    jamais en RHS nu 'double x = f[];' ; pas de 'continue' dans foreach ;
 *    tableaux en taille litterale. Compile SERIE.
 *
 * usage : ./mechanism <snapshot> <label>
 */
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

#define WIDTH 120.0
#define R0    8.0

int main (int argc, char ** argv)
{
  if (argc < 3) { fprintf (ferr, "usage: %s <snapshot> <label>\n", argv[0]); return 1; }
  size (WIDTH); origin (0., 0., 0.);
  foreach_dimension() periodic (right);              // OBLIGATOIRE avant restore
  if (!restore (file = argv[1])) {
    fprintf (ferr, "mechanism: lecture impossible de '%s'\n", argv[1]);
    return 1;
  }
  double D = 2.*R0, k = 2.*pi/L0;

  /* --- centre de la bulle (moyenne circulaire) + sa vitesse (somme complete) --- */
  double Vb=0., Wb=0., csx=0.,csy=0.,csz=0., ssx=0.,ssy=0.,ssz=0.;
  foreach (reduction(+:Vb) reduction(+:Wb)
           reduction(+:csx) reduction(+:csy) reduction(+:csz)
           reduction(+:ssx) reduction(+:ssy) reduction(+:ssz)) {
    double g = (1. - clamp(f[],0.,1.))*dv();
    Vb += g; Wb += g*u.z[];
    csx += g*cos(k*x); ssx += g*sin(k*x);
    csy += g*cos(k*y); ssy += g*sin(k*y);
    csz += g*cos(k*z); ssz += g*sin(k*z);
  }
  if (Vb <= 0.) { fprintf (fout, "%s nan\n", argv[2]); return 0; }
  double cx=atan2(ssx,csx)/k; if(cx<0.) cx+=L0;
  double cy=atan2(ssy,csy)/k; if(cy<0.) cy+=L0;
  double cz=atan2(ssz,csz)/k; if(cz<0.) cz+=L0;
  double ub = Wb/Vb;                                  // vitesse bulle (stockee)

  /* --- moyenne de u_z sur tout le liquide (champ lointain), poids frac. liquide --- */
  double accU=0., accW=0.;
  foreach (reduction(+:accU) reduction(+:accW)) {
    double w = dv()*clamp(f[],0.,1.);
    accW += w; accU += w*u.z[];
  }
  double uz_far = accU/accW;

  /* --- coquilles autour de la bulle (image minimale, poids frac. liquide) --- */
  double r1[4]; r1[0]=1.0; r1[1]=1.5; r1[2]=2.0; r1[3]=2.5;
  double r2[4]; r2[0]=1.5; r2[1]=2.0; r2[2]=2.5; r2[3]=3.0;
  double su_cone[4], sw_cone[4], su_hemi[4], sw_hemi[4], su_full[4], sw_full[4];
  for (int b=0; b<4; b++) { su_cone[b]=sw_cone[b]=su_hemi[b]=sw_hemi[b]=su_full[b]=sw_full[b]=0.; }
  foreach() {
    double w = dv()*clamp(f[],0.,1.);                 // poids = fraction liquide
    double dx = x - cx, dy = y - cy, dz = z - cz;
    if (dx >  L0/2.) dx -= L0; if (dx < -L0/2.) dx += L0;
    if (dy >  L0/2.) dy -= L0; if (dy < -L0/2.) dy += L0;
    if (dz >  L0/2.) dz -= L0; if (dz < -L0/2.) dz += L0;
    double rr = sqrt(dx*dx + dy*dy + dz*dz);
    double rD = rr/D;
    for (int b=0; b<4; b++)
      if (rD >= r1[b] && rD < r2[b]) {
        su_full[b] += w*u.z[]; sw_full[b] += w;
        if (dz > 0.)     { su_hemi[b] += w*u.z[]; sw_hemi[b] += w; }
        if (dz > 0.7*rr) { su_cone[b] += w*u.z[]; sw_cone[b] += w; }
      }
  }

  fprintf (fout, "%s %g %g %g", argv[2], cz, ub, uz_far);
  for (int b=0; b<4; b++) fprintf (fout, " %g", sw_cone[b]>0. ? su_cone[b]/sw_cone[b] : 0.);
  for (int b=0; b<4; b++) fprintf (fout, " %g", sw_hemi[b]>0. ? su_hemi[b]/sw_hemi[b] : 0.);
  for (int b=0; b<4; b++) fprintf (fout, " %g", sw_full[b]>0. ? su_full[b]/sw_full[b] : 0.);
  fprintf (fout, " %g\n", Vb);
  return 0;
}
