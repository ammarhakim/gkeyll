#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_gkgeom.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

#include <math.h>

// Cerfon equilibrium
struct cerfon_ctx {
  double R0, psi_prefactor;
};

static inline double sq(double x) { return x*x; }
static inline double cub(double x) { return x*x*x; }
static inline double qad(double x) { return x*x*x*x; }
static inline double pen(double x) { return x*x*x*x*x; }
static inline double hex(double x) { return x*x*x*x*x*x; }

void
psi_cerfon(double t, const double *xn, double *fout, void *ctx)
{
  struct cerfon_ctx *s = ctx;
  double R0 = s->R0, psi_prefactor = s->psi_prefactor;
  double R = xn[0], Z = xn[1];
  double x = R/R0, y = Z/R0;

  fout[0] = psi_prefactor*(0.00373804283369699*hex(x)*log(x) - 0.00574955335438162*hex(x) - 0.0448565140043639*qad(x)*sq(y)*log(x) + 0.0503044260840946*qad(x)*sq(y) + 0.017623348727471*qad(x)*log(x) + 0.0956643504553683*qad(x) + 0.0299043426695759*sq(x)*qad(y)*log(x) - 0.0160920841654771*sq(x)*qad(y) - 0.0704933949098842*sq(x)*sq(y)*log(x) + 0.0644725519961135*sq(x)*sq(y) - 7.00898484784405e-5*sq(x)*log(x) - 0.303766642191745*sq(x) - 0.00199362284463839*hex(y) + 0.0117488991516474*qad(y) + 7.00898484784405e-5*sq(y) + 0.0145368720253975);
}

void
cerforn_rt(void)
{
  fprintf(stdout, "---- Cerforn Double-Null Configuration\n");
  
  struct cerfon_ctx ctx = {  .R0 = 2.5, .psi_prefactor = 1.0 };
  
  // create RZ grid
  double lower[] = { 0.01, -6.0 }, upper[] = { 6.0, 6.0 };
  int cells[] = { 64, 128 };

  struct gkyl_rect_grid rzgrid;
  gkyl_rect_grid_init(&rzgrid, 2, lower, upper, cells);

  // RZ ranges
  struct gkyl_range rzlocal, rzlocal_ext;
  int nghost[GKYL_MAX_CDIM] = { 0, 0 };
  gkyl_create_grid_ranges(&rzgrid, nghost, &rzlocal_ext, &rzlocal);

  // RZ basis function
  int rz_poly_order = 2;
  struct gkyl_basis rzbasis;
  gkyl_cart_modal_serendip(&rzbasis, 2, rz_poly_order);

  // allocate psiRZ array, initialize and write it to file
  struct gkyl_array *psiRZ = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  
  gkyl_eval_on_nodes *eon = gkyl_eval_on_nodes_new(&rzgrid,
    &rzbasis, 1, &psi_cerfon, &ctx);
  gkyl_eval_on_nodes_advance(eon, 0.0, &rzlocal, psiRZ);
  gkyl_eval_on_nodes_release(eon);

  gkyl_grid_sub_array_write(&rzgrid, &rzlocal, 0, psiRZ, "cerfon_psi.gkyl");

  gkyl_gkgeom *geo = gkyl_gkgeom_new(&(struct gkyl_gkgeom_inp) {
      // psiRZ and related inputs
      .rzgrid = &rzgrid,
      .rzbasis = &rzbasis,
      .psiRZ = psiRZ,
      .rzlocal = &rzlocal
    }
  );

  int cum_nroots = 0;

  // compute R for various psi, Z
  do {
    double psi = 0.060095, Z = -5.611532889;
    double R[2] = { 0.0 }, dR[2] = { 0.0 };
    int nr = gkyl_gkgeom_R_psiZ(geo, psi, Z, 2, R, dR);
    for (int i=0; i<nr; ++i)
      fprintf(stdout, " >> r[%d] R(%g,%g) = %g\n", i, psi, Z, R[i]);
  } while (0);

  do {
    double psi = 0.060095, Z = 5.611532879;
    double R[2] = { 0.0 }, dR[2] = { 0.0 };
    int nr = gkyl_gkgeom_R_psiZ(geo, psi, Z, 2, R, dR);
    for (int i=0; i<nr; ++i)
      fprintf(stdout, " >> r[%d] R(%g,%g) = %g\n", i, psi, Z, R[i]);
  } while (0);  

  // compute arc-length of various flux-surfaces
  do {
    double psi_ref = 1e-4; // close to the seperatrix
    double arcL = gkyl_gkgeom_integrate_psi_contour(geo, psi_ref,
      lower[1], upper[1], upper[0]);

    struct gkyl_gkgeom_stat stat = gkyl_gkgeom_get_stat(geo);
    fprintf(stdout, "psi = %lg has arc-length %.10lg. Calls to contour func = %ld\n",
      psi_ref, arcL, stat.nquad_cont_calls-cum_nroots);
    
    cum_nroots += stat.nquad_cont_calls;
  } while(0);

  do {
    double psi_ref = 1.2;
    double arcL = gkyl_gkgeom_integrate_psi_contour(geo, psi_ref,
      lower[1], upper[1], upper[0]);

    struct gkyl_gkgeom_stat stat = gkyl_gkgeom_get_stat(geo);
    fprintf(stdout, "psi = %lg has arc-length %.10lg. Calls to contour func = %ld\n",
      psi_ref, arcL, stat.nquad_cont_calls-cum_nroots);

    cum_nroots += stat.nquad_cont_calls;
  } while(0);

  do {
    // compute outboard SOL geometry
    int npsi = 10, ntheta = 16;
    double psi_min = 0.0001, psi_max = 1.2;
    double dpsi = (psi_max-psi_min)/npsi;
    double dtheta = M_PI/ntheta;
  
    // Computational grid: theta X psi X alpha (only 2D for now)
    double clower[] = { -M_PI/2, psi_min };
    double cupper[] = { M_PI/2, psi_max };
    int ccells[] = { 16, 10 };
    
    struct gkyl_rect_grid cgrid;
    gkyl_rect_grid_init(&cgrid, 2, clower, cupper, ccells);

    // create mpc2p DG array
    struct gkyl_range clocal, clocal_ext;
    gkyl_create_grid_ranges(&cgrid, (int[]) { 0, 0, 0 },
      &clocal_ext, &clocal);

    int cpoly_order = 2;
    struct gkyl_basis cbasis;
    gkyl_cart_modal_serendip(&cbasis, 2, cpoly_order);
    struct gkyl_array *mapc2p = gkyl_array_new(GKYL_DOUBLE, 2*cbasis.num_basis, clocal_ext.volume);
    
    struct gkyl_gkgeom_geo_inp ginp = {
      .cgrid = &cgrid,
      .cbasis = &cbasis,
      .ftype = GKYL_GEOM_SOL_DN,
      .rclose = upper[0],
      .zmin = lower[1],
      .zmax = upper[1],
    
      .write_node_coord_array = true,
      .node_file_nm = "cerfon_out_sol_nod.gkyl"
    };

    gkyl_gkgeom_calcgeom(geo, &ginp, mapc2p);
    
    struct gkyl_gkgeom_stat stat = gkyl_gkgeom_get_stat(geo);
    fprintf(stdout, "Total number of contour funcs called = %ld. Total calls from root-finder = %ld\n",
      stat.nquad_cont_calls-cum_nroots, stat.nroot_cont_calls);

    gkyl_array_release(mapc2p);
  } while(0);

  do {
    // compute inboard SOL geometry
    int npsi = 2;
    double psi_min = 0.0001, psi_max = 0.01;
    double dpsi = (psi_max-psi_min)/npsi;
  
    // Computational grid: theta X psi X alpha (only 2D for now)
    double clower[] = { -M_PI/2, psi_min };
    double cupper[] = { M_PI/2, psi_max };
    int ccells[] = { 16, npsi };
    
    struct gkyl_rect_grid cgrid;
    gkyl_rect_grid_init(&cgrid, 2, clower, cupper, ccells);

    // create mpc2p DG array
    struct gkyl_range clocal, clocal_ext;
    gkyl_create_grid_ranges(&cgrid, (int[]) { 0, 0, 0 },
      &clocal_ext, &clocal);

    int cpoly_order = 2;
    struct gkyl_basis cbasis;
    gkyl_cart_modal_serendip(&cbasis, 2, cpoly_order);
    struct gkyl_array *mapc2p = gkyl_array_new(GKYL_DOUBLE, 2*cbasis.num_basis, clocal_ext.volume);
    
    struct gkyl_gkgeom_geo_inp ginp = {
      .cgrid = &cgrid,
      .cbasis = &cbasis,
      .ftype = GKYL_GEOM_SOL_DN,
      .rclose = lower[0],
      .zmin = lower[1],
      .zmax = upper[1],
    
      .write_node_coord_array = true,
      .node_file_nm = "cerfon_in_sol_nod.gkyl"
    };

    gkyl_gkgeom_calcgeom(geo, &ginp, mapc2p);
    
    struct gkyl_gkgeom_stat stat = gkyl_gkgeom_get_stat(geo);
    fprintf(stdout, "Total number of contour funcs called = %ld. Total calls from root-finder = %ld\n",
      stat.nquad_cont_calls-cum_nroots, stat.nroot_cont_calls);

    gkyl_array_release(mapc2p);
  } while(0);

  gkyl_gkgeom_release(geo);
  gkyl_array_release(psiRZ);
}

// WHAM equilibrium
struct wham_ctx {
  double B, gamma, Zm;
};

void
psi_wham(double t, const double *xn, double *fout, void *ctx)
{
  struct wham_ctx *s = ctx;
  double B = s->B, gamma = s->gamma, Zm = s->Zm;
  double R = xn[0], Z = xn[1];

  // double Lorentzian: See Francisquez PoP 2023.
  double psi = sq(R)*B/(2*M_PI*gamma)*
    (
      1/(1+sq((Z-Zm)/gamma))
      +
      1/(1+sq((Z+Zm)/gamma))
    );

  fout[0] = psi;
}

void
wham_2l_rt(void)
{
  fprintf(stdout, "---- WHAM Configuration\n");
  
  struct wham_ctx ctx = {
    .B = 6.51292,
    .gamma = 0.124904,
    .Zm = 0.98
  };
  
  // create RZ grid
  double lower[] = { 0.01, -2.0 };
  double upper[] = { 0.4, 2.0 };
  int cells[] = { 64, 128 };

  struct gkyl_rect_grid rzgrid;
  gkyl_rect_grid_init(&rzgrid, 2, lower, upper, cells);

  // RZ ranges
  struct gkyl_range rzlocal, rzlocal_ext;
  int nghost[GKYL_MAX_CDIM] = { 0, 0 };
  gkyl_create_grid_ranges(&rzgrid, nghost, &rzlocal_ext, &rzlocal);

  // RZ basis function
  int rz_poly_order = 2;
  struct gkyl_basis rzbasis;
  gkyl_cart_modal_serendip(&rzbasis, 2, rz_poly_order);

  // allocate psiRZ array, initialize and write it to file
  struct gkyl_array *psiRZ = gkyl_array_new(GKYL_DOUBLE, rzbasis.num_basis, rzlocal_ext.volume);
  
  gkyl_eval_on_nodes *eon = gkyl_eval_on_nodes_new(&rzgrid,
    &rzbasis, 1, &psi_wham, &ctx);
  gkyl_eval_on_nodes_advance(eon, 0.0, &rzlocal, psiRZ);
  gkyl_eval_on_nodes_release(eon);

  gkyl_grid_sub_array_write(&rzgrid, &rzlocal, 0, psiRZ, "wham_psi.gkyl");

  gkyl_gkgeom *geo = gkyl_gkgeom_new(&(struct gkyl_gkgeom_inp) {
      // psiRZ and related inputs
      .rzgrid = &rzgrid,
      .rzbasis = &rzbasis,
      .psiRZ = psiRZ,
      .rzlocal = &rzlocal
    }
  );

  int cum_nroots = 0;

  do {
    // compute outboard SOL geometry
    int npsi = 10, ntheta = 16;
    double psi_min = 0.001, psi_max = 0.02;
    double dpsi = (psi_max-psi_min)/npsi;
    double dtheta = M_PI/ntheta;
  
    // Computational grid: theta X psi X alpha (only 2D for now)
    double clower[] = { -M_PI/2, psi_min };
    double cupper[] = { M_PI/2, psi_max };
    int ccells[] = { 16, 10 };
    
    struct gkyl_rect_grid cgrid;
    gkyl_rect_grid_init(&cgrid, 2, clower, cupper, ccells);

    // create mpc2p DG array
    struct gkyl_range clocal, clocal_ext;
    gkyl_create_grid_ranges(&cgrid, (int[]) { 0, 0, 0 },
      &clocal_ext, &clocal);

    int cpoly_order = 2;
    struct gkyl_basis cbasis;
    gkyl_cart_modal_serendip(&cbasis, 2, cpoly_order);
    struct gkyl_array *mapc2p = gkyl_array_new(GKYL_DOUBLE, 2*cbasis.num_basis, clocal_ext.volume);
    
    struct gkyl_gkgeom_geo_inp ginp = {
      .cgrid = &cgrid,
      .cbasis = &cbasis,
      .ftype = GKYL_GEOM_SOL_DN,
      .rclose = upper[0],
      .zmin = lower[1],
      .zmax = upper[1],
    
      .write_node_coord_array = true,
      .node_file_nm = "wham_out_sol_nod.gkyl"
    };

    gkyl_gkgeom_calcgeom(geo, &ginp, mapc2p);
    
    struct gkyl_gkgeom_stat stat = gkyl_gkgeom_get_stat(geo);
    fprintf(stdout, "Total number of contour funcs called = %ld. Total calls from root-finder = %ld\n",
      stat.nquad_cont_calls-cum_nroots, stat.nroot_cont_calls);

    gkyl_array_release(mapc2p);
  } while(0);

  gkyl_gkgeom_release(geo);
  gkyl_array_release(psiRZ);
}

void
wham_beta0_rt(void)
{
  fprintf(stdout, "---- WHAM beta-0 Configuration\n");
  
  struct gkyl_rect_grid rzgrid;
  struct gkyl_array *psiRZ =
    gkyl_grid_array_new_from_file(&rzgrid, "psi_dg.gkyl");

  // RZ ranges
  struct gkyl_range rzlocal, rzlocal_ext;
  int nghost[GKYL_MAX_CDIM] = { 0, 0 };
  gkyl_create_grid_ranges(&rzgrid, nghost, &rzlocal_ext, &rzlocal);  
  
  // RZ basis function
  int rz_poly_order = 2;
  struct gkyl_basis rzbasis;
  gkyl_cart_modal_serendip(&rzbasis, 2, rz_poly_order);  
  
  gkyl_gkgeom *geo = gkyl_gkgeom_new(&(struct gkyl_gkgeom_inp) {
      // psiRZ and related inputs
      .rzgrid = &rzgrid,
      .rzbasis = &rzbasis,
      .psiRZ = psiRZ,
      .rzlocal = &rzlocal
    }
  );

  int cum_nroots = 0;

  do {
    // compute outboard SOL geometry
    int npsi = 10, ntheta = 16;
    //double psi_min = 1.0e-4, psi_max = 3.16726875e-03;
    double psi_min = 0.0e-6, psi_max = 3.16726875e-03;    
    double dpsi = (psi_max-psi_min)/npsi;
    double dtheta = M_PI/ntheta;
  
    // Computational grid: theta X psi X alpha (only 2D for now)
    double clower[] = { -M_PI/2, psi_min };
    double cupper[] = { M_PI/2, psi_max };
    int ccells[] = { 16, 10 };
    
    struct gkyl_rect_grid cgrid;
    gkyl_rect_grid_init(&cgrid, 2, clower, cupper, ccells);

    // create mpc2p DG array
    struct gkyl_range clocal, clocal_ext;
    gkyl_create_grid_ranges(&cgrid, (int[]) { 0, 0, 0 },
      &clocal_ext, &clocal);

    int cpoly_order = 2;
    struct gkyl_basis cbasis;
    gkyl_cart_modal_serendip(&cbasis, 2, cpoly_order);
    struct gkyl_array *mapc2p = gkyl_array_new(GKYL_DOUBLE, 2*cbasis.num_basis, clocal_ext.volume);
    
    struct gkyl_gkgeom_geo_inp ginp = {
      .cgrid = &cgrid,
      .cbasis = &cbasis,
      .ftype = GKYL_GEOM_SOL_DN,
      .rclose = rzgrid.upper[0],
      .zmin = -2.0, //rzgrid.lower[1],
      .zmax = 2.0, //rzgrid.upper[1],
    
      .write_node_coord_array = true,
      .node_file_nm = "wham_out_sol_nod.gkyl"
    };

    gkyl_gkgeom_calcgeom(geo, &ginp, mapc2p);
    
    struct gkyl_gkgeom_stat stat = gkyl_gkgeom_get_stat(geo);
    fprintf(stdout, "Total number of contour funcs called = %ld. Total calls from root-finder = %ld\n",
      stat.nquad_cont_calls-cum_nroots, stat.nroot_cont_calls);

    gkyl_array_release(mapc2p);
  } while(0);

  gkyl_gkgeom_release(geo);
  gkyl_array_release(psiRZ);  
}

int
main(int argc, char **argcv)
{
  cerforn_rt();
  wham_2l_rt();
  //wham_beta0_rt();
  
  return 0;
}
