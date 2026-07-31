/**
 * main.c  --  Bulle de gaz en HIT (Homogeneous Isotropic Turbulence)
 * --------------------------------------------------------------------
 * DNS 3D (Octree) avec Basilisk.
 *
 * Objectif : mesurer la vitesse de remontee d'une bulle (Bo = 1) dans un
 * ecoulement turbulent triplement periodique, et la comparer au cas laminaire.
 *
 * Pipeline en 2 phases :
 *   1) PRECURSEUR  (INJECT=0) : champ ABC -> turbulence forcee/decroissante,
 *      monophasique (f = 1 partout). On dump le champ stationnaire ("end").
 *   2) BULLE       (INJECT=1) : on restaure le precurseur ("restart"), on
 *      raffine autour de la sphere, on injecte la bulle (fraction VOF), et on
 *      laisse remonter sous l'effet de la gravite reduite.
 *
 * Physique imposee :
 *   - Gravite : force volumique PERIODIQUE a_z = -GRAV (1 - rho_m/rho) (PAS reduced.h,
 *     dont le potentiel lineaire epingle la bulle a la frontiere ; cf. avertissement).
 *   - Tension de surface DYNAMIQUE pour garantir Bo = 1 (convention diametre):
 *         Bo = (rho1 - rho2) |G_Z| (2 R0)^2 / sigma = 1
 *     =>  f.sigma = (rho1 - rho2) * fabs(G_Z) * sq(2*R0).
 *   - Adaptation sur l'interface VOF 'f' ET la vitesse 'u'.
 *   - Toutes les ecritures fichiers sont protegees par if (pid() == 0).
 *
 * Compilation : voir Makefile (qcc -D_MPI=1 -O3). Option graphique : -DMOVIE=1.
 *
 * Arguments CLI (positionnels) :
 *   argv[1] maxlevel   (7, 8, 9 ...)
 *   argv[2] MAXTIME    (temps physique final)
 *   argv[3] R0         (rayon de la bulle)
 *   argv[4] FORCED     (1 = forcage lineaire de la turbulence, 0 = libre)
 *   argv[5] KE_TARGET  (energie cinetique cible du forcage controle)
 *   argv[6] INJECT     (0 = precurseur, 1 = injection de la bulle)
 *   argv[7] BOND       (= Bo vise, defaut 1 ; Bo = BOND exactement.
 *                       A REPASSER IDENTIQUE lors d'une reprise depuis "dump" :
 *                       f.sigma est refixe par event init a chaque demarrage)
 *   argv[8] OMECO      (coefficient du seuil AMR vorticite, defaut 0.20 :
 *                       omemax = OMECO * <|omega|> ; PLUS BAS = PLUS RAFFINE)
 *   argv[9]  GRAV      (= |g|, defaut 1. Bouton We_t A Bo FIXE :
 *                       We_t ∝ 1/GRAV. N'affecte pas le precurseur monophasique.)
 *   argv[10] NU        (= viscosite mu1, defaut = formule si <0. Bouton Re_lambda ;
 *                       tout changement de NU exige de REFAIRE le precurseur.)
 *   argv[11] MEMBER    (membre d'ensemble, defaut 0 : decale la position
 *                       d'injection de la bulle -> realisation turbulente
 *                       independante, pour moyenner plusieurs mesures par point.)
 *
 *  Probleme = 3 nombres sans dimension independants (Buckingham) : (Bo, We_t, Re_λ)
 *  pilotes par (BOND, GRAV, NU). Le 4e groupe (Ga, Mo, Re_b) est alors DERIVE.
 */

/* ATTENTION : PAS de "reduced.h" ici -- voir event acceleration.
   reduced.h exprime la poussee via le potentiel phi = [rho] G.(x - Z), LINEAIRE en z
   donc NON PERIODIQUE : il saute de [rho]*g*L0 a la frontiere z=L0. Des que la bulle
   traverse le haut du domaine, ce saut cree une force parasite qui l'EPINGLE sur la
   frontiere (constate : en fluide au repos la bulle monte a uz~12 puis se fige a
   zc = L0 - R). Incompatible avec un domaine periodique selon la gravite. */
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"   // advection VOF conservant la qdm (haut RHOR)
#include "tension.h"
/* (reduced.h RETIRE : non periodique selon g -- cf. avertissement en tete) */
#include "tag.h"
#include "navier-stokes/perfs.h"
#include "maxruntime.h"

#if MOVIE
# include "lambda2.h"
# include "view.h"
#endif

/* ----------------------- Parametres physiques ------------------------ */
#define RHOR   850.0      // rho_liquide / rho_gaz
#define MUR    100.0      // mu_liquide  / mu_gaz
#define WIDTH  120.0      // taille du domaine cubique periodique
#define SNAP_DT 1.0       // periode des snapshots (dump) pour le rendu
#define T_AMR_START 18.0  // precurseur : passe en AMR (vorticite) apres cet instant
                          // (avant, uniforme -> le bruit initial ne fausse pas l'AMR)

/* Champ d'acceleration pour le forcage lineaire de la turbulence. */
face vector av[];

/* ----------------------- Parametres pilotables ----------------------- */
int    maxlevel  = 8;
double MAXTIME   = 60.0;
double R0        = 8.0;
int    FORCED    = 1;
double KE_TARGET = 24.0;   // energie cinetique cible (forcage controle)
int    INJECT    = 0;      // 0 = precurseur ; 1 = bulle
double OMECO     = 0.20;   // seuil AMR vorticite = OMECO * <|omega|>

/* --- Les 3 boutons sans dimension (probleme = 3 groupes independants, Buckingham) --- */
double BOND      = 1.0;    // -> Bo (Bo = BOND exactement, par construction de f.sigma)
double GRAV      = 1.0;    // -> |g| (G.z = -GRAV). DECOUPLE We_t de Bo : a Bo fixe et
                           // turbulence fixe, We_t ∝ 1/GRAV. g n'affecte PAS le
                           // precurseur monophasique -> meme end_we1 pour l'axe We_t.
double NU        = -1.0;   // -> viscosite liquide mu1 (rho1=1). <0 => valeur par defaut
                           // (formule). Pilote Re_lambda (ν↓ => Re_λ↑) ; borne par la
                           // resolution (η) -> Re_λ ≲ 55 a lvl7. Axe Re = NOUVEAU precurseur.
int    MINLEVEL  = 5;      // niveau de base (2^5=32^3) ET plancher des adapt_wavelet
int    MEMBER    = 0;      // membre d'ensemble : decale la position d'injection de la
                          // bulle (turbulence homogene -> realisation independante).
                          // 0 = centre du domaine (defaut). >0 -> table d'offsets.

/* Stats de turbulence (calculees dans event logfile), partagees avec le forcage
   controle et l'adaptation sur vorticite. */
coord  global_ubar = {0., 0., 0.};
double global_ke   = 0.;
double global_eps  = 0.;
/* <f> = fraction volumique moyenne de liquide -> densite moyenne rho_m du domaine,
   necessaire a la poussee periodique-compatible (event acceleration). 1 = monophasique. */
double global_fbar = 1.0;

/* ---------------- REPERE MOBILE (suit la bulle) ---------------------------
   Pourquoi : avec la poussee volumique la bulle monte librement (bien), mais elle
   TRAVERSE le maillage sur ~1900 unites en 150 u.t. -> l'erreur d'advection VOF
   s'accumule (mesure : +1.2 % de volume par traversee de domaine).
   Remede : invariance galileenne. On retranche a TOUT le champ la vitesse de montee
   de la bulle -> la bulle reste au centre, le liquide defile devant elle. Son
   interface ne se deplace plus dans la grille -> plus d'erreur d'advection.
   On ne perd rien : la vraie trajectoire lab est integree dans frame_z, et la
   vitesse terminale est frame_uz a l'equilibre.
   Subtilite : un repere qui ACCELERE est non-inertiel (pseudo-force). A la vitesse
   terminale frame_uz est constant -> repere inertiel -> exact. On relaxe donc
   doucement (FRAME_TAU) au lieu de recaler brutalement a chaque pas : pas de
   pseudo-force violente, et les fluctuations turbulentes de la bulle sont preservees. */
#define FRAME_TAU 1.0      // temps de relaxation du repere (u.t.) ; 0 = repere fixe
double frame_uz = 0.0;     // vitesse du repere = vitesse LAB moyenne de la bulle
double frame_z  = 0.0;     // position LAB de la bulle (integree) = vraie trajectoire

int main (int argc, char * argv[])
{
  maxruntime (&argc, argv);

  if (argc > 1) maxlevel  = atoi(argv[1]);
  if (argc > 2) MAXTIME   = atof(argv[2]);
  if (argc > 3) R0        = atof(argv[3]);
  if (argc > 4) FORCED    = atoi(argv[4]);
  if (argc > 5) KE_TARGET = atof(argv[5]);
  if (argc > 6) INJECT    = atoi(argv[6]);
  if (argc > 7) BOND      = atof(argv[7]);
  if (argc > 8) OMECO     = atof(argv[8]);
  if (argc > 9) GRAV      = atof(argv[9]);
  if (argc > 10) NU       = atof(argv[10]);
  if (argc > 11) MEMBER   = atoi(argv[11]);

  size (WIDTH);
  origin (0., 0., 0.);
  foreach_dimension()
    periodic (right);                    // domaine triplement periodique (HIT)

  /* Proprietes des deux phases. */
  rho1 = 1.0;            // liquide  (f = 1)
  rho2 = rho1/RHOR;      // gaz      (f = 0)  -> la bulle
  mu1  = (NU > 0.0 ? NU : 0.01*sq(WIDTH/(2.0*pi))/2.0);   // ν = mu1 (rho1=1) -> Re_lambda
  mu2  = mu1/MUR;

  /* Gravite : PAS de G/reduced.h. La poussee est appliquee comme force volumique
     periodique dans event acceleration (a_z = -GRAV (1 - rho_m/rho)). GRAV = |g|. */

  /* Champ de forcage de la turbulence. */
  a = av;

  /* Precurseur : grille UNIFORME au niveau 'maxlevel' (turbulence homogene, rien
     a gagner a adapter). Bulle : base au niveau MINLEVEL (2^5) + adaptation jusqu'a
     maxlevel (l'adapt a aussi MINLEVEL comme plancher). */
  if (INJECT == 0)
    N = 1 << maxlevel;
  else
    N = 1 << MINLEVEL;

  /* Defaut Basilisk (1e-3), pas 1e-4 : au niveau 8, le saut rho1/rho2=RHOR,
     mu1/mu2=MUR rend le solveur visqueux raide (mg_solve incremente nrelax
     jusqu'a son plafond 100 des que la reduction du residu est <1.2x/cycle,
     nrelax se propageant d'un pas a l'autre) -> mgu.nrelax reste bloque a
     98-100 au niveau 8 (~3 au niveau 7), cout par pas ~70x celui de lvl7. */
  TOLERANCE = 1e-3;
  CFL = 0.5;
  run();
}

/* ============================ INITIALISATION ========================= */
event init (i = 0)
{
  /* Tension de surface dynamique -> Bo = BOND (convention diametre) :
       Bo = (rho1 - rho2) g (2 R0)^2 / sigma ,  g = GRAV
     => f.sigma = (rho1 - rho2) GRAV (2 R0)^2 / BOND  => Bo = BOND EXACTEMENT,
     independamment de GRAV. C'est GRAV qui deplace We_t a Bo fixe (We_t ∝ 1/GRAV). */
  f.sigma = (rho1 - rho2)*GRAV*sq(2.0*R0)/BOND;
  if (pid() == 0)
    fprintf (ferr, "[PARAMS] maxlevel=%d MINLEVEL=%d R0=%g | Bo=BOND=%g GRAV=%g "
             "sigma=%g nu=mu1=%g | OMECO=%g KE_TARGET=%g INJECT=%d MEMBER=%d\n",
             maxlevel, MINLEVEL, R0, BOND, GRAV, f.sigma, mu1, OMECO,
             KE_TARGET, INJECT, MEMBER);

  /* 1) Reprise depuis un checkpoint glissant. */
  if (restore (file = "dump")) {
    /* L'etat du repere mobile (frame_uz, frame_z) n'est PAS dans le dump : sans
       le restaurer, une reprise apres walltime reinitialiserait frame_uz a 0 et
       toute la trajectoire lab serait faussee (silencieusement). On relit la
       derniere ligne de frame.dat anterieure au t restaure (le fichier peut
       contenir des lignes posterieures au dump : on les ignore). */
    if (INJECT && FRAME_TAU > 0.) {
      double fu = 0., fz = 0.;
      if (pid() == 0) {
        FILE * fp = fopen("frame.dat", "r");
        if (fp) {
          char line[256];
          double tt, u1, z1, ub1;
          while (fgets(line, sizeof line, fp))
            if (line[0] != '#' &&
                sscanf(line, "%lf %lf %lf %lf", &tt, &u1, &z1, &ub1) == 4 &&
                tt <= t + 1e-9) {
              fu = u1; fz = z1;
            }
          fclose(fp);
        }
      }
#if _MPI
      MPI_Bcast(&fu, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&fz, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
      frame_uz = fu;
      frame_z  = fz;
      if (pid() == 0)
        fprintf (ferr, "[FRAME] reprise t=%g : frame_uz=%g frame_z=%g "
                 "(relus de frame.dat)\n", t, frame_uz, frame_z);
    }
    return 0;
  }

  /* 2) Mode bulle : restaurer le precurseur turbulent puis injecter la bulle.
     NB: restore() segfault au-dela de ~128^3 (bug Basilisk a grande echelle sur
     ce build) -> le precurseur DOIT etre en niveau 7 (128^3), pas 8 (256^3). */
  if (INJECT && restore (file = "restart")) {
    /* Position d'injection : centre du domaine, decalee en (x,y) pour les
       membres d'ensemble (MEMBER>0). La turbulence etant homogene, injecter a
       des positions bien separees (~L0/2, > echelle integrale) donne des
       realisations INDEPENDANTES du meme champ -> moyenne d'ensemble pour des
       points fiables. Offsets choisis dans [R0, L0-R0] (pas de bulle a cheval sur
       une frontiere, ou fraction() la couperait) ; z inchange (meme altitude de
       depart). Table cyclique -> autant de membres que voulu. */
    double ox[6] = {  30., -30.,  30., -30.,  40.,   0. };
    double oy[6] = {  30.,  30., -30., -30.,   0.,  40. };
    double dx = 0., dy = 0.;
    if (MEMBER > 0) { int k = (MEMBER - 1) % 6; dx = ox[k]; dy = oy[k]; }
    double xc = L0/2. + dx, yc = L0/2. + dy, zc = L0/2.;
    if (pid() == 0)
      fprintf (ferr, "[MEMBER] %d : injection en (%g, %g, %g)\n", MEMBER, xc, yc, zc);

    /* Raffiner la coquille autour de la future interface AVANT fraction(). */
    refine (sq(x-xc) + sq(y-yc) + sq(z-zc) < sq(1.30*R0) &&
            sq(x-xc) + sq(y-yc) + sq(z-zc) > sq(0.70*R0) &&
            level < maxlevel);

    /* f = 1 a l'exterieur (liquide), f = 0 a l'interieur (gaz = bulle). */
    fraction (f, sq(x-xc) + sq(y-yc) + sq(z-zc) - sq(R0));

    /* La bulle demarre au repos. FORCED=1 : la turbulence est conservee dans le
       liquide (u*=f -> vitesse annulee seulement dans le gaz). FORCED=0 : REFERENCE
       LAMINAIRE -> on annule la vitesse PARTOUT (liquide au repos), et l'absence de
       forcage (event acceleration inactif) laisse la bulle remonter en fluide
       quiescent -> vitesse terminale laminaire V_lam. */
    foreach()
      foreach_dimension()
        u.x[] = (FORCED ? u.x[] : 0.0)*f[];
    return 0;
  }

  /* 3) Demarrage a froid : champ ABC monophasique (precurseur).
     Amplitude calee sur KE_TARGET (ke_ABC = 3/2 amp^2 -> amp = sqrt(2 ke/3)) pour
     que le forcage controle n'ait qu'a MAINTENIR l'energie, pas a la corriger
     d'un facteur 20 (transitoire interminable sinon). */
  /* IC = SUPERPOSITION de modes ABC (k=1,2,3) + un peu de bruit. Le mode ABC k=1
     seul est une solution de Beltrami quasi-stationnaire qui ne transitionne PAS
     (le bruit blanc se dissipe avant d'amorcer). Plusieurs modes incommensurables
     interagissent non-lineairement -> chaos immediat -> transition rapide vers la
     turbulence. Chaque mode ABC est a divergence nulle. Normalise pour ke ~ KE_TARGET
     (le forcage controle ajuste ensuite). */
  double waveno = WIDTH/(2.0*pi);
  double wk[3]  = {1.0, 0.7, 0.5};               // poids des modes k=1,2,3
  double norm   = sqrt(wk[0]*wk[0] + wk[1]*wk[1] + wk[2]*wk[2]);
  double amp    = sqrt(2.0*KE_TARGET/3.0)/norm;
  foreach() {
    f[] = 1.0;
    double xw = x/waveno, yw = y/waveno, zw = z/waveno;
    double ax = 0., ay = 0., az = 0.;
    for (int kk = 1; kk <= 3; kk++) {
      double c = wk[kk-1];
      ax += c*( cos(kk*yw) + sin(kk*zw) );
      ay += c*( sin(kk*xw) + cos(kk*zw) );
      az += c*( cos(kk*xw) + sin(kk*yw) );
    }
    u.x[] = amp*ax + 0.10*amp*noise();
    u.y[] = amp*ay + 0.10*amp*noise();
    u.z[] = amp*az + 0.10*amp*noise();
  }
}

/* ===================== FORCAGE CONTROLE A ENERGIE CIBLE ============= */
/* Forcage lineaire controle (Bassenne et al. 2016) : on ajuste l'amplitude A a
   chaque pas pour tenir ke = KE_TARGET. A = [eps + (KE_TARGET - ke)/tau]/(2 ke),
   avec tau = ke/eps (temps de retournement). Au stationnaire A = eps/(2 ke) -> la
   puissance injectee (2 A ke) egale la dissipation -> energie BORNEE, pas
   d'emballement, et l'intensite turbulente = KE_TARGET (choisie). On force le
   liquide seul (facteur f). ke/eps/ubar viennent de l'event logfile (decalage
   d'un pas, sans incidence pour un controleur a relaxation). */
event acceleration (i++)
{
  if (FORCED && global_ke > 0.) {
    double tau = global_ke/(global_eps + 1e-30);
    double A   = (global_eps + (KE_TARGET - global_ke)/tau)/(2.*global_ke);
    foreach_face()
      av.x[] += A*(0.5*(f[] + f[-1]))*
                ((u.x[] + u.x[-1])/2. - global_ubar.x);
  }

  /* ---- POUSSEE D'ARCHIMEDE COMPATIBLE AVEC LA PERIODICITE ----
     On remplace rho*g par (rho - rho_m)*g, soit l'acceleration
         a_z = g (1 - rho_m/rho)   avec g = -GRAV (vers le bas),
     rho_m = densite moyenne du domaine. Proprietes :
       - la force volumique (rho - rho_m) g est PERIODIQUE (rho l'est) : aucun
         potentiel lineaire, donc la bulle traverse la frontiere et CONTINUE de
         monter (contrairement a reduced.h qui l'epinglait, cf. en-tete) ;
       - sa moyenne est nulle (int (rho - rho_m) g dV = 0) : pas d'emballement de
         la quantite de mouvement totale, ce qui est requis en periodique ;
       - en monophasique (f=1) : rho = rho_m = rho1 -> a_z = 0, donc le PRECURSEUR
         n'est pas affecte par GRAV (propriete conservee, meme end_we1 reutilisable).
     La bulle (rho2 << rho_m) recoit a_z ~ +GRAV*rho_m/rho2 (tres grand) ; c'est la
     formulation standard non-well-balanced : la projection de pression equilibre et
     l'acceleration NETTE de la bulle vaut ~2g (masse ajoutee). */
  double rhom = rho1*global_fbar + rho2*(1.0 - global_fbar);
  foreach_face(z) {
    double ff   = clamp((f[] + f[0,0,-1])/2., 0., 1.);
    double rhof = ff*(rho1 - rho2) + rho2;
    av.z[] += -GRAV*(1.0 - rhom/rhof);
  }
}

/* ===================== REPERE MOBILE : SUIVI DE LA BULLE ============ */
/* Retranche a tout le champ une fraction (dt/FRAME_TAU) de la vitesse verticale de
   la bulle, pour que celle-ci reste centree dans le maillage. Transformation de
   Galilee : la physique est inchangee (a la pseudo-force pres, nulle a l'equilibre).
   -> frame_uz converge vers u_inf ; frame_z est la vraie trajectoire lab.
   Ecrit frame.dat : t frame_uz frame_z ub_dans_repere zc_dans_maillage           */
event move_frame (i++)
{
  if (!INJECT || FRAME_TAU <= 0.) return 0;   // precurseur : pas de repere mobile

  /* vitesse verticale de la bulle DANS le repere courant (moyenne sur le gaz) */
  double vg = 0., wg = 0.;
  foreach (reduction(+:vg) reduction(+:wg)) {
    double gas = (1. - f[])*dv();
    wg += gas;
    vg += gas*u.z[];
  }
  if (wg <= 1e-30) return 0;                  // pas de gaz -> rien a suivre
  double ub = vg/wg;

  /* relaxation douce : shift = ub * dt/tau. A l'equilibre ub->0 (le repere a
     rattrape la bulle) donc shift->0 : pas de derive parasite. */
  double shift = ub*dt/FRAME_TAU;
  foreach()
    u.z[] -= shift;
  frame_uz += shift;                          // vitesse lab cumulee de la bulle
  frame_z  += frame_uz*dt;                    // position lab (trajectoire vraie)

  if (pid() == 0) {
    static FILE * ff = NULL;
    if (!ff) {
      ff = fopen("frame.dat", "a");
      fprintf(ff, "# t frame_uz frame_z ub_dans_repere\n");
    }
    if (i % 5 == 0) {                         // meme cadence que bubble.dat
      fprintf(ff, "%g %g %g %g\n", t, frame_uz, frame_z, ub);
      fflush(ff);
    }
  }
}

/* ===================== STATISTIQUES DE TURBULENCE =================== */
/* Calcule ubar, ke, dissipation eps -> stocke en global (forcage + adapt) et
   ecrit stats.dat (+ Re_lambda et echelle de Kolmogorov eta). */
event logfile (i++)
{
  coord ubar;
  foreach_dimension() {
    stats s = statsf(u.x);
    ubar.x = s.sum/s.volume;
  }
  double ke = 0., vd = 0., vol = 0., fsum = 0.;
  foreach (reduction(+:ke) reduction(+:vd) reduction(+:vol) reduction(+:fsum)) {
    vol  += dv();
    fsum += dv()*f[];                    // -> <f> pour rho_m (poussee periodique)
    foreach_dimension() {
      ke += dv()*sq(u.x[] - ubar.x);
      vd += dv()*(sq(u.x[1] - u.x[-1]) +
                  sq(u.x[0,1] - u.x[0,-1]) +
                  sq(u.x[0,0,1] - u.x[0,0,-1]))/sq(2.*Delta);
    }
  }
  ke /= 2.*vol;
  vd *= mu1/vol;
  global_fbar = fsum/vol;                // densite moyenne rho_m = rho1*<f>+rho2*(1-<f>)

  /* ---- GARDE-FOU : QUANTITE DE MOUVEMENT TOTALE (repere LABORATOIRE) -------
     En domaine triplement periodique il n'existe AUCUNE force nette exterieure :
     la poussee est a moyenne nulle par construction (rho-rho_m)g, la tension de
     surface est une force interne, le forcage turbulent est centre. Donc
     P_lab = int rho u_lab dV doit rester CONSTANT (~0).
     ATTENTION au REPERE MOBILE : u stocke = vitesse DANS LE REPERE de la bulle
     (u_lab = u + frame_uz e_z). Le pz "brut" tend donc vers -frame_uz PAR
     CONSTRUCTION (le liquide defile vers le bas devant la bulle maintenue au
     centre) : ce n'est PAS une derive. Verifie sur le run lvl6 (job 38746) :
     pz = -10.654 vs frame_uz = +10.649 a t=190 -> pz_lab ~ -0.005 ~ px,py.
     L'INVARIANT a surveiller est le trio (px, py, pz_lab = pz + frame_uz) ~ 0.
     S'il derive -> une force fantome injecte de la qdm -> BUG, et toute vitesse
     mesuree est suspecte. */
  double px = 0., py = 0., pz = 0., mass = 0.;
  foreach (reduction(+:px) reduction(+:py) reduction(+:pz) reduction(+:mass)) {
    double rhoc = clamp(f[],0.,1.)*(rho1 - rho2) + rho2;
    double m = rhoc*dv();
    mass += m;
    px += m*u.x[]; py += m*u.y[]; pz += m*u.z[];
  }
  if (pid() == 0 && mass > 0.) {
    static FILE * fm = NULL;
    if (!fm) {
      fm = fopen("momentum.dat", "a");
      fprintf(fm, "# t px py pz frame_uz pz_lab  (vitesses barycentriques ; "
                  "invariant = px, py, pz_lab ~ 0 ; pz seul -> -frame_uz, normal)\n");
    }
    fprintf(fm, "%g %g %g %g %g %g\n", t, px/mass, py/mass, pz/mass,
            frame_uz, pz/mass + frame_uz);
    fflush(fm);
  }
  double Re  = (vd > 0. ? 2./3.*ke/mu1*sqrt(15.*mu1/vd) : 0.);
  double eta = (vd > 0. ? pow(mu1*mu1*mu1/vd, 0.25) : 0.);   // nu=mu1, rho1=1

  /* Export pour le forcage controle et l'adaptation sur vorticite. */
  global_ubar = ubar;
  global_ke   = ke;
  global_eps  = vd;

  if (pid() == 0) {
    static FILE * fd = NULL;
    if (!fd) {
      fd = fopen("stats.dat", "a");
      fprintf(fd, "# t dissipation energy Reynolds eta\n");
    }
    fprintf(fd, "%g %g %g %g %g\n", t, vd, ke, Re, eta);
    fflush(fd);
    fprintf(ferr, "i=%d t=%g ke=%g eps=%g Re_l=%g eta=%g\n", i, t, ke, vd, Re, eta);
  }
}

/* ===================== SUIVI DE LA BULLE ============================ */
/* Bo = 1 : pas de fragmentation attendue. On mesure volume, centre, vitesse,
   boite englobante, diametre equivalent et rapport d'aspect chi.
 *
 * 🚨 REECRIT le 2026-07-22 -- l'ancienne version etait FAUSSE de deux facons,
 *    toutes deux demontrees sur le run LAMINAIRE (ou la reponse est connue) :
 *
 *  (1) SOMME TRONQUEE. Elle construisait le marqueur binaire m = (f < 0.5),
 *      tagguait, puis ne sommait QUE sur la region "mere". Or une maille a
 *      f = 0.6 contient 40 % de gaz mais n'est pas taguee -> son gaz etait
 *      JETE. Comme l'interface fait 1-2 mailles, cela revenait a jeter la
 *      moitie exterieure de la peau de la bulle. Mesure sur snapshot-200.00 :
 *      806 mailles jetees, 4.5 % du volume. Consequences :
 *        - volume rapporte 95.0 % du theorique alors que le VRAI est 99.5 %
 *          (=> les "3-7 % de perte de volume" de tout l'historique du projet
 *           etaient un artefact de MESURE, pas une perte physique : le VOF
 *           conserve tres bien, d_eq = 15.994 vs 16.0 theorique) ;
 *        - vitesse verticale biaisee de +0.35 a +0.48 (3-4 % de u_inf). Sens
 *          explique : en repere mobile le liquide defile vers le bas, donc les
 *          mailles jetees (les plus riches en liquide) sont les plus negatives ;
 *          les retirer REMONTE la moyenne. Pas de biais en x,y (pas de vitesse
 *          moyenne dans ces directions) -> accord a 3 % la-bas.
 *
 *  (2) [NON, FAUSSE ALERTE -- garde en trace pour ne pas la refaire] J'avais
 *      conclu le 2026-07-22 que tag() ne gerait pas la periodicite, sur la foi
 *      d'un programme de diagnostic autonome qui voyait la bulle coupee en deux
 *      des qu'elle chevauchait un bord. C'ETAIT UN BUG DE CE PROGRAMME : il
 *      appelait restore() sans declarer periodic(), or restore() ne restaure NI
 *      la taille du domaine NI les conditions aux limites. Verification : la
 *      vraie simulation n'a jamais rapporte la moindre chute (bubble.dat reste a
 *      94.8-95.5 % sur toute la fenetre incriminee), et le meme diagnostic avec
 *      periodic() declare donne n_regions = 1 partout. => tag() est CORRECT ici,
 *      et rien n'infirme la fragmentation observee au run KE=200.
 *      LECON : tout post-traitement qui fait restore() DOIT re-declarer size(),
 *      origin() et periodic() a l'identique de main().
 *
 * CORRECTIF : volume / centre / vitesse / boite sont calcules par des sommes
 * COMPLETES ponderees par (1-f) sur TOUTES les mailles, avec moyenne circulaire
 * pour le centre. tag() est conserve pour la detection de fragmentation, mais
 * relegue a des colonnes SUPPLEMENTAIRES (n_regions, V_tag) : la mesure
 * principale ne depend plus de lui. Le drapeau "straddle" (bulle a cheval sur un
 * bord) est ecrit a titre informatif.
 *
 * Les colonnes 1-13 gardent leur ordre et leur signification (les valeurs, elles,
 * sont desormais correctes) -> les scripts existants continuent de lire juste.  */
event bubble (i += 5)
{
  /* --- Passe 1 : sommes COMPLETES sur tout le gaz (toutes les mailles).
     Centre par moyenne circulaire => robuste a la periodicite. --- */
  double k = 2.*pi/L0;
  double Vb = 0.;
  coord csum = {0,0,0}, ssum = {0,0,0}, vsum = {0,0,0};
  foreach_leaf() {
    double w = dv()*(1. - clamp(f[], 0., 1.));
    if (w > 0.) {
      Vb += w;
      coord p = {x, y, z};
      foreach_dimension() {
        csum.x += w*cos(k*p.x);
        ssum.x += w*sin(k*p.x);
        vsum.x += w*u.x[];
      }
    }
  }
#if _MPI
  MPI_Allreduce(MPI_IN_PLACE, &Vb,   1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &csum, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &ssum, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &vsum, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (Vb <= 1e-30) return 0;              // pas de gaz (precurseur) -> rien

  coord com, ub;
  foreach_dimension() {
    com.x = atan2(ssum.x, csum.x)/k;      // dans [-L/2, L/2]
    if (com.x < 0.) com.x += L0;          // -> [0, L0)
    ub.x  = vsum.x/Vb;
  }

  /* --- Passe 2 : boite englobante via image minimale autour du COM.
     Critere f < 0.5 = l'isosurface qui DEFINIT la surface de la bulle. Ne PAS
     prendre "toute maille contenant du gaz" (f < 1-eps) : l'advection VOF laisse
     des residus microscopiques (f ~ 0.9999) loin de la bulle, sans effet sur le
     volume mais qui font exploser la boite (mesure : depth 35.6 au lieu de 16.9,
     chi 1.73 au lieu de 1.38). Pas besoin de tag() ici : l'image minimale autour
     du COM gere deja la periodicite. --- */
  coord lo = { 1e30,  1e30,  1e30};
  coord hi = {-1e30, -1e30, -1e30};
  foreach_leaf()
    if (f[] < 0.5) {                      // maille majoritairement gaz
      coord p = {x, y, z};
      foreach_dimension() {
        double d = p.x - com.x;
        if (d >  L0/2.) d -= L0;
        if (d < -L0/2.) d += L0;
        if (d < lo.x) lo.x = d;
        if (d > hi.x) hi.x = d;
      }
    }
#if _MPI
  MPI_Allreduce(MPI_IN_PLACE, &lo, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &hi, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
  double width  = hi.x - lo.x;           // extension selon x
  double depth  = hi.y - lo.y;           // extension selon y
  double height = hi.z - lo.z;           // extension selon z (direction de montee)
  double d_eq   = 2.0*pow(3.0*Vb/(4.0*pi), 1.0/3.0);
  /* Aspect ratio = plus grande longueur / plus petite longueur de la bulle. */
  double lmax   = max(width, max(depth, height));
  double lmin   = min(width, min(depth, height));
  double chi    = (lmin > 0. ? lmax/lmin : 1.);

  /* --- Passe 3 : DIAGNOSTIC de fragmentation (colonnes supplementaires).
     tag() ne recollant pas les morceaux a travers une frontiere periodique, on
     signale explicitement le chevauchement : la boite englobante deborde-t-elle
     du domaine ? Si oui -> n_regions et V_tag sont ininterpretables. --- */
  int straddle = 0;
  foreach_dimension()
    if (com.x + hi.x > L0 || com.x + lo.x < 0.) straddle = 1;

  scalar m[];
  foreach()
    m[] = (f[] < 0.5);
  int n = tag(m);
  double Vtag = 0.;
  if (n >= 1) {
    double * vol = calloc(n, sizeof(double));
    foreach_leaf()
      if (m[] > 0)
        vol[(int)m[] - 1] += dv()*(1. - f[]);
#if _MPI
    MPI_Allreduce(MPI_IN_PLACE, vol, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    for (int j = 0; j < n; j++)
      if (vol[j] > Vtag) Vtag = vol[j];
    free(vol);
  }

  if (pid() == 0) {
    static FILE * fb = NULL;
    if (!fb) {
      fb = fopen("bubble.dat", "a");
      fprintf(fb, "# t volume xc yc zc ub_x ub_y ub_z width depth height d_eq chi"
                  " n_regions V_tag straddle\n"
                  "# col 1-13 : sommes COMPLETES ponderees (1-f) sur toutes les"
                  " mailles, centre en moyenne circulaire.\n"
                  "#            (l'ancienne version ne sommait que sur les mailles"
                  " taguees f<0.5 -> volume sous-estime de ~4.5 %%,\n"
                  "#             vitesse verticale biaisee de ~+0.35 en repere"
                  " mobile. Cf. event bubble dans main.c.)\n"
                  "# col 14-16 : diagnostic fragmentation par tag() (V_tag ="
                  " plus grosse region ; straddle = bulle a cheval sur un bord).\n");
    }
    fprintf(fb, "%g %g %g %g %g %g %g %g %g %g %g %g %g %d %g %d\n",
            t, Vb, com.x, com.y, com.z, ub.x, ub.y, ub.z,
            width, depth, height, d_eq, chi, n, Vtag, straddle);
    fflush(fb);
  }
}

/* =========================== ADAPTATION ============================= */
#if TREE
event adapt (i++)
{
  /* Precurseur : UNIFORME tant que le bruit/transitoire initial n'est pas dissipe
     (adapter pendant le demarrage bruite raffinerait PARTOUT -> inutile + explosion
     de mailles). Apres T_AMR_START : AMR sur la vorticite pour alleger le dump
     (contourne le bug de restore a 256^3 uniforme si l'AMR reduit assez). */
  if (INJECT == 0 && t < T_AMR_START)
    return 0;

  /* Vorticite |omega| = |rot(u)| : cible les petites echelles de la turbulence
     (pas la vitesse, qui suit aussi les gros tourbillons -- remarque des tuteurs). */
  scalar omega[];
  foreach() {
    double wx = (u.z[0,1,0] - u.z[0,-1,0] - (u.y[0,0,1] - u.y[0,0,-1]))/(2.*Delta);
    double wy = (u.x[0,0,1] - u.x[0,0,-1] - (u.z[1]    - u.z[-1]   ))/(2.*Delta);
    double wz = (u.y[1]     - u.y[-1]     - (u.x[0,1,0] - u.x[0,-1,0]))/(2.*Delta);
    omega[] = sqrt(wx*wx + wy*wy + wz*wz);
  }
  /* Seuil cale sur l'intensite reelle (fraction de la vorticite moyenne).
     NB: c'est une TOLERANCE D'ERREUR d'ondelette (plus bas = plus raffine),
     pas un detecteur "raffine ou omega > seuil". adapt_wavelet deraffine
     implicitement quand l'erreur passe sous omemax/1.5 (hysteresis). */
  double omemax = OMECO*normf(omega).avg + 1e-30;

  if (INJECT == 0)
    /* Precurseur (apres transitoire) : AMR sur la vorticite seule.
       Plancher MINLEVEL (4e arg) : pas de deraffinement sous 2^5. */
    adapt_wavelet ((scalar *){omega}, (double[]){omemax}, maxlevel, MINLEVEL);
  else {
    /* Bulle : interface (f) jusqu'au maxlevel + vorticite pour la turbulence.
       Plancher MINLEVEL (4e arg) : le fond turbulent ne deraffine pas sous 2^5. */
    double femax = 1e-3;
    adapt_wavelet ((scalar *){f, omega}, (double[]){femax, omemax}, maxlevel,
                   MINLEVEL);
  }
}
#endif

/* ===================== COMPTAGE DES CELLULES (MPI) ================== */
/* grid->n est LOCAL a chaque rang ; on somme sur tous les rangs pour le total
   global, ce qui permet de dimensionner le nombre de coeurs (~100k cell/coeur).
   Le compte grimpe avec le developpement de la turbulence puis plafonne. */
event cells (i += 20)
{
  long nc = 0;
  foreach (reduction(+:nc))
    nc++;
  if (pid() == 0) {
    static FILE * fp = NULL;
    if (!fp) {
      fp = fopen("cells.dat", "a");
      fprintf(fp, "# i t cells coeurs_suggere(~100k/coeur)\n");
    }
    fprintf(fp, "%d %g %ld %ld\n", i, t, nc, nc/100000 + 1);
    fflush(fp);
    fprintf(ferr, "[CELLS] i=%d t=%g cells=%ld -> ~%ld coeurs\n",
            i, t, nc, nc/100000 + 1);
  }
}

/* ================ HISTOGRAMME DES NIVEAUX DE RAFFINEMENT ============ */
/* Ou vit la grille ? Compte les feuilles par niveau -> levels.dat. Montre (1) que
   l'AMR raffine bien la turbulence/interface vers maxlevel et (2) le deraffinement
   implicite d'adapt_wavelet (les niveaux bas se repeuplent quand |omega| retombe). */
#define NLEVELS 16
event levels (i += 100)
{
  long nl[NLEVELS] = {0};
  foreach_leaf()
    if (level < NLEVELS)
      nl[level]++;
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, nl, NLEVELS, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (pid() == 0) {
    static FILE * fl = NULL;
    if (!fl) {
      fl = fopen("levels.dat", "a");
      fprintf(fl, "# t cellules_par_niveau_0_a_%d\n", NLEVELS - 1);
    }
    fprintf(fl, "%g", t);
    for (int l = 0; l < NLEVELS; l++)
      fprintf(fl, " %ld", nl[l]);
    fputc('\n', fl);
    fflush(fl);
  }
}

/* ===================== CHECKPOINTS & SNAPSHOTS ====================== */
event checkpoint (i += 200)
{
  dump (file = "dump");                  // checkpoint glissant (reprise)
}

event snapshots (t += SNAP_DT; t <= MAXTIME)
{
  char name[80];
  sprintf(name, "snapshot-%06.2f", t);   // pour le rendu offline
  dump (file = name);
}

/* =========================== FILM (option) ========================== */
#if MOVIE
event movie (t += 0.5)
{
  view (fov = 40, camera = "iso", width = 800, height = 800,
        bg = {1,1,1}, samples = 4);
  clear();
  squares ("u.z", linear = true, n = {0,0,1}, alpha = L0/2.);
  draw_vof ("f");
  save ("movie.mp4");
}
#endif

event end (t = MAXTIME)
{
  dump (file = "end");                   // champ final (= "restart" du precurseur)
}
