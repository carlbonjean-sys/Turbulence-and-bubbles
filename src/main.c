/**
 * main.c  --  DNS d'une bulle de gaz en turbulence homogene et isotrope (HIT)
 * ---------------------------------------------------------------------------
 * Solveur DNS 3D (Octree / Basilisk).
 *
 * Configuration physique :
 *   - Poussee d'Archimede : force volumique periodique a_z = -GRAV (1 - rho_m/rho)
 *   - Tension de surface : f.sigma ajustee pour imposer Bo = BOND
 *   - Convergence & adaptation : AMR sur VOF 'f' et vorticite |omega|
 */

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "tag.h"
#include "navier-stokes/perfs.h"
#include "maxruntime.h"

#if MOVIE
# include "lambda2.h"
# include "view.h"
#endif

#define RHOR   850.0
#define MUR    100.0
#define WIDTH  120.0
#define SNAP_DT 1.0
#define T_AMR_START 18.0

face vector av[];

int    maxlevel  = 8;
double MAXTIME   = 60.0;
double R0        = 8.0;
int    FORCED    = 1;
double KE_TARGET = 24.0;
int    INJECT    = 0;
double OMECO     = 0.20;

double BOND      = 1.0;
double GRAV      = 1.0;
double NU        = -1.0;
int    MINLEVEL  = 5;
int    MEMBER    = 0;

coord  global_ubar = {0., 0., 0.};
double global_ke   = 0.;
double global_eps  = 0.;
double global_fbar = 1.0;

#define FRAME_TAU 1.0
double frame_uz = 0.0;
double frame_z  = 0.0;

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
    periodic (right);

  rho1 = 1.0;
  rho2 = rho1/RHOR;
  mu1  = (NU > 0.0 ? NU : 0.01*sq(WIDTH/(2.0*pi))/2.0);
  mu2  = mu1/MUR;

  a = av;

  if (INJECT == 0)
    N = 1 << maxlevel;
  else
    N = 1 << MINLEVEL;

  TOLERANCE = 1e-3;
  CFL = 0.5;
  run();
}

event init (i = 0)
{
  f.sigma = (rho1 - rho2)*GRAV*sq(2.0*R0)/BOND;
  if (pid() == 0)
    fprintf (ferr, "[PARAMS] maxlevel=%d MINLEVEL=%d R0=%g | Bo=%g GRAV=%g "
             "sigma=%g nu=%g | OMECO=%g KE_TARGET=%g INJECT=%d MEMBER=%d\n",
             maxlevel, MINLEVEL, R0, BOND, GRAV, f.sigma, mu1, OMECO,
             KE_TARGET, INJECT, MEMBER);

  if (restore (file = "dump")) {
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
    }
    return 0;
  }

  if (INJECT && restore (file = "restart")) {
    double ox[6] = {  30., -30.,  30., -30.,  40.,   0. };
    double oy[6] = {  30.,  30., -30., -30.,   0.,  40. };
    double dx = 0., dy = 0.;
    if (MEMBER > 0) { int k = (MEMBER - 1) % 6; dx = ox[k]; dy = oy[k]; }
    double xc = L0/2. + dx, yc = L0/2. + dy, zc = L0/2.;

    refine (sq(x-xc) + sq(y-yc) + sq(z-zc) < sq(1.30*R0) &&
            sq(x-xc) + sq(y-yc) + sq(z-zc) > sq(0.70*R0) &&
            level < maxlevel);

    fraction (f, sq(x-xc) + sq(y-yc) + sq(z-zc) - sq(R0));

    foreach()
      foreach_dimension()
        u.x[] = (FORCED ? u.x[] : 0.0)*f[];
    return 0;
  }

  double waveno = WIDTH/(2.0*pi);
  double wk[3]  = {1.0, 0.7, 0.5};
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

event acceleration (i++)
{
  if (FORCED && global_ke > 0.) {
    double tau = global_ke/(global_eps + 1e-30);
    double A   = (global_eps + (KE_TARGET - global_ke)/tau)/(2.*global_ke);
    foreach_face()
      av.x[] += A*(0.5*(f[] + f[-1]))*
                ((u.x[] + u.x[-1])/2. - global_ubar.x);
  }

  double rhom = rho1*global_fbar + rho2*(1.0 - global_fbar);
  foreach_face(z) {
    double ff   = clamp((f[] + f[0,0,-1])/2., 0., 1.);
    double rhof = ff*(rho1 - rho2) + rho2;
    av.z[] += -GRAV*(1.0 - rhom/rhof);
  }
}

event move_frame (i++)
{
  if (!INJECT || FRAME_TAU <= 0.) return 0;

  double vg = 0., wg = 0.;
  foreach (reduction(+:vg) reduction(+:wg)) {
    double gas = (1. - f[])*dv();
    wg += gas;
    vg += gas*u.z[];
  }
  if (wg <= 1e-30) return 0;
  double ub = vg/wg;

  double shift = ub*dt/FRAME_TAU;
  foreach()
    u.z[] -= shift;
  frame_uz += shift;
  frame_z  += frame_uz*dt;

  if (pid() == 0) {
    static FILE * ff = NULL;
    if (!ff) {
      ff = fopen("frame.dat", "a");
      fprintf(ff, "# t frame_uz frame_z ub_dans_repere\n");
    }
    if (i % 5 == 0) {
      fprintf(ff, "%g %g %g %g\n", t, frame_uz, frame_z, ub);
      fflush(ff);
    }
  }
}

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
    fsum += dv()*f[];
    foreach_dimension() {
      ke += dv()*sq(u.x[] - ubar.x);
      vd += dv()*(sq(u.x[1] - u.x[-1]) +
                  sq(u.x[0,1] - u.x[0,-1]) +
                  sq(u.x[0,0,1] - u.x[0,0,-1]))/sq(2.*Delta);
    }
  }
  ke /= 2.*vol;
  vd *= mu1/vol;
  global_fbar = fsum/vol;

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
      fprintf(fm, "# t px py pz frame_uz pz_lab\n");
    }
    fprintf(fm, "%g %g %g %g %g %g\n", t, px/mass, py/mass, pz/mass,
            frame_uz, pz/mass + frame_uz);
    fflush(fm);
  }
  double Re  = (vd > 0. ? 2./3.*ke/mu1*sqrt(15.*mu1/vd) : 0.);
  double eta = (vd > 0. ? pow(mu1*mu1*mu1/vd, 0.25) : 0.);

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
  }
}

event bubble (i += 5)
{
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
  if (Vb <= 1e-30) return 0;

  coord com, ub;
  foreach_dimension() {
    com.x = atan2(ssum.x, csum.x)/k;
    if (com.x < 0.) com.x += L0;
    ub.x  = vsum.x/Vb;
  }

  coord lo = { 1e30,  1e30,  1e30};
  coord hi = {-1e30, -1e30, -1e30};
  foreach_leaf()
    if (f[] < 0.5) {
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
  double width  = hi.x - lo.x;
  double depth  = hi.y - lo.y;
  double height = hi.z - lo.z;
  double d_eq   = 2.0*pow(3.0*Vb/(4.0*pi), 1.0/3.0);
  double lmax   = max(width, max(depth, height));
  double lmin   = min(width, min(depth, height));
  double chi    = (lmin > 0. ? lmax/lmin : 1.);

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
      fprintf(fb, "# t volume xc yc zc ub_x ub_y ub_z width depth height d_eq chi n_regions V_tag straddle\n");
    }
    fprintf(fb, "%g %g %g %g %g %g %g %g %g %g %g %g %g %d %g %d\n",
            t, Vb, com.x, com.y, com.z, ub.x, ub.y, ub.z,
            width, depth, height, d_eq, chi, n, Vtag, straddle);
    fflush(fb);
  }
}

#if TREE
event adapt (i++)
{
  if (INJECT == 0 && t < T_AMR_START)
    return 0;

  scalar omega[];
  foreach() {
    double wx = (u.z[0,1,0] - u.z[0,-1,0] - (u.y[0,0,1] - u.y[0,0,-1]))/(2.*Delta);
    double wy = (u.x[0,0,1] - u.x[0,0,-1] - (u.z[1]    - u.z[-1]   ))/(2.*Delta);
    double wz = (u.y[1]     - u.y[-1]     - (u.x[0,1,0] - u.x[0,-1,0]))/(2.*Delta);
    omega[] = sqrt(wx*wx + wy*wy + wz*wz);
  }
  double omemax = OMECO*normf(omega).avg + 1e-30;

  if (INJECT == 0)
    adapt_wavelet ((scalar *){omega}, (double[]){omemax}, maxlevel, MINLEVEL);
  else {
    double femax = 1e-3;
    adapt_wavelet ((scalar *){f, omega}, (double[]){femax, omemax}, maxlevel,
                   MINLEVEL);
  }
}
#endif

event cells (i += 20)
{
  long nc = 0;
  foreach (reduction(+:nc))
    nc++;
  if (pid() == 0) {
    static FILE * fp = NULL;
    if (!fp) {
      fp = fopen("cells.dat", "a");
      fprintf(fp, "# i t cells\n");
    }
    fprintf(fp, "%d %g %ld\n", i, t, nc);
    fflush(fp);
  }
}

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
      fprintf(fl, "# t cellules_par_niveau\n");
    }
    fprintf(fl, "%g", t);
    for (int l = 0; l < NLEVELS; l++)
      fprintf(fl, " %ld", nl[l]);
    fputc('\n', fl);
    fflush(fl);
  }
}

event checkpoint (i += 200)
{
  dump (file = "dump");
}

event snapshots (t += SNAP_DT; t <= MAXTIME)
{
  char name[80];
  sprintf(name, "snapshot-%06.2f", t);
  dump (file = name);
}

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
  dump (file = "end");
}
