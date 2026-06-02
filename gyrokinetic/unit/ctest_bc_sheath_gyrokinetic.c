/* This is the unit test for the sheath BC updater.
  * The test does the following:
  * 1. Initializes the sheath BC updater with a given configuration space dimension, edge.
  * 2. Project a maxwellian distribution function with input parameters to the edge of the domain (skin cells).
  * 3. Apply the sheath BC updater to the distribution function.
  * 4. Check the ghost values of the distribution function (where the result of the sheath BC is stored).
  *    We check that the reflected function is zero in cells that are fully outside the cutoff 
  *    and non zero in cells that are fully inside. You can check these numbers by setting pars.verbose to true.
  * 5. Deallocate memory used by the updater.
  * 
  * This is done for each configuration space dimension (1x2v, 2x2v, 3x2v), for both lower and upper edge, and
  * positive, negative potential, for electrons and ions.
  * 
  * The test is also able to verify the interface of the surrogate model for the sheath BC by setting 
  * pars.use_surrogate to true. In this case the check uses the surrogate directly to obtain the vcut.
  * 
  * Finally, one can set write_fields to true to write the distribution function and visualize it with postgkyl.
  * 
*/

#include <acutest.h>

#include <gkyl_array_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>
#include <gkyl_bc_sheath_gyrokinetic.h>
#include <gkyl_velocity_map.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_const.h>
#include <float.h>

// Path to the KANN surrogate model used in surrogate tests.
// Set to NULL (or set srgrz_test_enabled = false) to skip surrogate tests.
static const char *srgrz_model_path = "/Users/ahoffman/gkeyll_sheath_ai/model/nn_model_conv_MPE.kann";
static const bool srgrz_test_enabled = true;

static struct gkyl_array*
mkarr(bool use_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (use_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct test_sheath_ctx {
  int cdim; // Configuration space dimension.
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM]; // Grid extents.
  int cells[GKYL_MAX_DIM]; // Number of cells.
  double mass; // Species mass.
  double charge; // Species charge.
  double dens; // Reference density.
  double upar; // Mean parallel drift speed.
  double temp; // Reference temperature.
  double B0; // Reference magnetic field.
  double phi_mpe; // ES potential at the sheath entrance.
  double phi_wall; // ES potential at the wall.
  double impact_angle; // Angle [rad] of incidence if B to the plate (for surrogate sheath).
  double x0; // x-center
  double y0; // y-center
  double z0; // z-center
  double sigmax; // Width of distribution in x.
  double sigmay; // Width of distribution in y.
  double sigmaz; // Width of distribution in z.
  bool use_surrogate; // Whether to use surrogate model for vcut factor.
  const char *surrogate_model_path; // Path to the .kann model file (NULL if not using surrogate).
  bool verbose; // Whether to print detailed output.
};

void
eval_func_vcut_fact(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct test_sheath_ctx *pars = ctx;
  double mu = xn[pars->cdim-1]; // Mu is the last dimension.

  if (pars->use_surrogate) {
    gkyl_bc_sheath_gyrokinetic_evaluate_vcut_fact_surrogate(
      &mu, 1, pars->phi_mpe, pars->phi_wall, pars->dens, 
      pars->temp, 2*pars->charge/pars->mass, pars->B0, pars->impact_angle, fout);
      fout[0] = fout[0];
  } else {
    double Lmu = 1.0; // Characteristic scale length in mu direction.
    double vcut_fact_0 = 1.0;
  
    double vcut_fact = vcut_fact_0 * (1 + exp(-mu/Lmu));
    fout[0] =  1;//pow(vcut_fact, 2);
  }
}

void
eval_func_1x2v(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double z = xn[0];
  double vpar = xn[1], mu = xn[2];

  struct test_sheath_ctx *pars = ctx;

  double vt = sqrt(pars->temp/pars->mass); // Thermal speed.
  double B0 = pars->B0;
  double m = pars->mass;
  double upar = pars->upar;
  double z0 = pars->z0;
  double sigmaz = pars->sigmaz;
  double T = pars->temp;

  double envelope = 1;//exp( -pow(z-z0,2)/(2.0*pow(sigmaz,2)) );
  fout[0] = exp( -(m*pow(vpar-upar,2)/2.0 + mu*B0) / T ) * envelope;
}

void
eval_func_2x2v(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], z = xn[1];
  double vpar = xn[2], mu = xn[3];

  struct test_sheath_ctx *pars = ctx;

  double vt = sqrt(pars->temp/pars->mass); // Thermal speed.
  double B0 = pars->B0;
  double m = pars->mass;
  double upar = pars->upar;
  double x0 = pars->x0;
  double z0 = pars->z0;
  double sigmax = pars->sigmax;
  double sigmaz = pars->sigmaz;
  double T = pars->temp;

  double envelope = 1;//exp( -pow(x-x0,2)/(2.0*pow(sigmax,2)) ) * exp( -pow(z-z0,2)/(2.0*pow(sigmaz,2)) );
  fout[0] = exp( -(m*pow(vpar-upar,2)/2.0 + mu*B0) / T ) * envelope;
}

void
eval_func_3x2v(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double vpar = xn[3], mu = xn[4];
  double x = xn[0], y = xn[1], z = xn[2];

  struct test_sheath_ctx *pars = ctx;

  double vt = sqrt(pars->temp/pars->mass); // Thermal speed.
  double B0 = pars->B0;
  double m = pars->mass;
  double upar = pars->upar;
  double x0 = pars->x0;
  double y0 = pars->y0;
  double z0 = pars->z0;
  double sigmax = pars->sigmax;
  double sigmay = pars->sigmay;
  double sigmaz = pars->sigmaz;
  double T = pars->temp;

  double envelope =1;//exp( -pow(x-x0,2)/(2.0*pow(sigmax,2)) ) * exp( -pow(y-y0,2)/(2.0*pow(sigmay,2)) ) * exp( -pow(z-z0,2)/(2.0*pow(sigmaz,2)) );
  fout[0] = exp( -(m*pow(vpar-upar,2)/2.0 + mu*B0) / T ) * envelope;
}

// Checks that the distribution function values in the ghost cells are set to 
// 0 for velocities beyond the cutoff velocity and >0 otherwise.
void
check_function(struct test_sheath_ctx *pars, struct gkyl_array *distf_ho, struct gkyl_rect_grid grid, struct gkyl_range ghost_r, enum gkyl_edge_loc edge)
{
  double delta_phi = pars->phi_mpe - pars->phi_wall;
  int cdim = pars->cdim;
  double vpar_cut_radius;
  // This cut is if we have reflection. If q*dphi > 0, we expect 0 everywhere.
  vpar_cut_radius = sqrt(fabs(pars->charge*delta_phi)/(0.5*pars->mass));

  int num_cells = 0;
  int num_zero_cells = 0;
  int num_zero_cells_expected = 0;
  int num_uknown_cells = 0;
  // We need the sign of charge*delta_phi to determine 
  // which part of the distribution function is cut off by the sheath BC.
  double qphi_sign = pars->charge * delta_phi > 0 ? 1.0 : -1.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &ghost_r);
  while (gkyl_range_iter_next(&iter)) {
    int *idx_g = iter.idx;
    long linidx_g = gkyl_range_idx(&ghost_r, idx_g);
    const double *distf_c = gkyl_array_cfetch(distf_ho, linidx_g);
    double tol = 1e-12;
    double ref_val = 0.0;

    // Get cell boundaries along vpar.
    double xc[GKYL_MAX_DIM] = {0};
    gkyl_rect_grid_cell_center(&grid, idx_g, xc);
    double cell_lower_vpar = xc[cdim] - 0.5*grid.dx[cdim];
    double cell_upper_vpar = xc[cdim] + 0.5*grid.dx[cdim];
    // Compute the maximal "velocity radius" of the cell
    double vr_min = fmin(fabs(cell_lower_vpar), fabs(cell_upper_vpar));
    double vr_max = fmax(fabs(cell_lower_vpar), fabs(cell_upper_vpar));

    // Compute the vpar cut radius at the lower and upper mu boundaries (assumes vcut_fact is monotonic in mu).
    double cell_lower_mu = xc[cdim+1] - 0.5*grid.dx[cdim+1];
    double cell_upper_mu = xc[cdim+1] + 0.5*grid.dx[cdim+1];
    double vcut_coord_lower[cdim], vcut_coord_upper[cdim];
    vcut_coord_lower[cdim-1] = cell_lower_mu;
    vcut_coord_upper[cdim-1] = cell_upper_mu;
    double vcut_fact_lower[1], vcut_fact_upper[1];
    eval_func_vcut_fact(0, vcut_coord_lower, vcut_fact_lower, pars);
    eval_func_vcut_fact(0, vcut_coord_upper, vcut_fact_upper, pars);
    double vpcr_lower = sqrt(vcut_fact_lower[0]) * vpar_cut_radius;
    double vpcr_upper = sqrt(vcut_fact_upper[0]) * vpar_cut_radius;
    double vpcr_min = fmin(vpcr_lower, vpcr_upper);
    double vpcr_max = fmax(vpcr_lower, vpcr_upper);

    num_cells++;
    num_zero_cells += distf_c[0] > tol ? 0 : 1;
    if ( (vr_min >= vpcr_max) || (qphi_sign > 0) ) { 
      // Check that all coefficients are 0 for cells that are fully outside the cutoff or if we have full absorption (q*delta_phi > 0).
      num_zero_cells_expected++;
      for (int k=0; k<distf_ho->ncomp; k++) {
        TEST_CHECK( gkyl_compare(distf_c[k], ref_val, tol) );
        TEST_MSG( "Expected %.9e | Got: %.9e at idx=%d,%d,%d\n", ref_val, distf_c[k], idx_g[0], idx_g[1], idx_g[2]);
      }
    } else if ( ((vpcr_min < vr_max)  && (vpcr_min > vr_min)) || ((vpcr_max > vr_min)  && (vpcr_max < vr_max)) ) {
      // For cells that are partially cut off by the sheath BC we don't know.
      num_uknown_cells++;
    } else {
      // These must be non zero cells that are fully within the cutoff, so we check that the average value is >0.
      double cell_avg = distf_c[0];
      TEST_CHECK( cell_avg > ref_val );
      TEST_MSG( "Expected > %.9e | Got: %.9e at idx=%d,%d,%d\n", ref_val, cell_avg, idx_g[0], idx_g[1], idx_g[2]);
    }
  }
  // Uncomment for debugging.
  // printf("\n");
  // printf("Edge: %s | Charge: %.1e | Delta phi: %.2e\n", 
  //   edge == GKYL_LOWER_EDGE? "lower" : "upper", charge, delta_phi);
  // printf("resulting zero cells = %d | fully outside cells = %d | crossing cells = %d | total cells = %d\n", num_zero_cells, num_zero_cells_expected, num_uknown_cells, num_cells);
}

void write_out_fields(struct test_sheath_ctx *pars, int cdim, int vdim, enum gkyl_edge_loc edge, bool use_gpu,
  struct gkyl_array *distf_ho, struct gkyl_array *phi_ho, struct gkyl_array *phiw_ho,
  struct gkyl_rect_grid grid_ext, struct gkyl_range local_ext,
  struct gkyl_rect_grid grid_conf, struct gkyl_range local_conf, struct gkyl_basis *basis, struct gkyl_basis basis_conf) {

  struct gkyl_rect_grid grid_norm_ext = grid_ext;
  // Normalize vparallel by thermal speed.
  double vt = sqrt(pars->temp/pars->mass);
  grid_norm_ext.dx[cdim] /= vt;
  grid_norm_ext.lower[cdim] /= vt;
  grid_norm_ext.upper[cdim] /= vt;
  // Normalize mu by B*vt^2.
  grid_norm_ext.dx[cdim+1] *= pars->B0/pars->temp;
  grid_norm_ext.lower[cdim+1] *= pars->B0/pars->temp;
  grid_norm_ext.upper[cdim+1] *= pars->B0/pars->temp;
  
  struct gkyl_msgpack_map_elem io_meta[] = {
    { .key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = basis->poly_order },
    { .key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = basis->id }
  };
  int io_meta_len = sizeof(io_meta)/sizeof(io_meta[0]);
  struct gkyl_msgpack_data *mt = gkyl_msgpack_create(io_meta_len, io_meta);

  struct gkyl_msgpack_map_elem io_meta_conf[] = {
    { .key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = basis_conf.poly_order },
    { .key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = basis_conf.id }
  };
  int io_meta_conf_len = sizeof(io_meta_conf)/sizeof(io_meta_conf[0]);
  struct gkyl_msgpack_data *mt_conf = gkyl_msgpack_create(io_meta_conf_len, io_meta_conf);

  char fname[256];
  const char *fmt = "bc_sheath_%dx%dv_%s_%s_%s.gkyl";
  snprintf(fname, sizeof(fname), fmt, cdim, vdim, edge == GKYL_LOWER_EDGE? "lower" : "upper", use_gpu? "gpu" : "cpu", "distf_out");
  gkyl_grid_sub_array_write(&grid_norm_ext, &local_ext, mt, distf_ho, fname);
  snprintf(fname, sizeof(fname), fmt, cdim, vdim, edge == GKYL_LOWER_EDGE? "lower" : "upper", use_gpu? "gpu" : "cpu", "phi_mpe");
  gkyl_grid_sub_array_write(&grid_conf, &local_conf, mt_conf, phi_ho, fname);
  snprintf(fname, sizeof(fname), fmt, cdim, vdim, edge == GKYL_LOWER_EDGE? "lower" : "upper", use_gpu? "gpu" : "cpu", "phi_wall");
  gkyl_grid_sub_array_write(&grid_conf, &local_conf, mt_conf, phiw_ho, fname);

  gkyl_msgpack_data_release(mt);
  gkyl_msgpack_data_release(mt_conf);
}

void
test_bc_sheath_gyrokinetic_1x2v(struct test_sheath_ctx *pars, enum gkyl_edge_loc edge, bool write_fields, bool use_gpu)
{
  /*
  This test applies the sheath BC to a Maxwellian distribution function on both upper and lower edges in the z direction,
  according to a potential drop between the magnetic presheath entrance and the wall. 
  The folloing files are written out if requested:
  - "bc_sheath_1x2v_distf_in_${edge}_${cpugpu}.gkyl": Input distribution function.
  - "bc_sheath_1x2v_phi_mpe.gkyl": Potential at the magnetic presheath entrance (constant everywhere).
  - "bc_sheath_1x2v_phi_wall.gkyl": Potential at the wall (constant everywhere).
  - "bc_sheath_1x2v_distf_out_${edge}_${cpugpu}.gkyl": Distribution after the sheath BC has been applied (includes ghost cells).
  */
  double ms = pars->mass;
  double qs = pars->charge;
  double ns = pars->dens;
  double Ts = pars->temp;
  double B0 = pars->B0;
  double phi_mpe = pars->phi_mpe;
  double phi_wall = pars->phi_wall;
  double impact_angle = pars->impact_angle;
  double vt = sqrt(Ts/ms);
  double q2Dm = 2*qs/ms;
  double mu0 = Ts/B0;

  double phi_surr = -qs * (phi_mpe - phi_wall) / Ts;
  double gamma_surr = sqrt(GKYL_ELECTRON_MASS * ns / GKYL_EPSILON0) / B0;
  double alpha_surr = impact_angle * 180.0 / GKYL_PI;

  // if (pars->use_surrogate) {
  //   printf("\n surrogate params: alpha = %.2f | gamma = %.2f | phi = %.2f\n", alpha_surr, gamma_surr, phi_surr);
  // }

  int poly_order = 1;
  double lower[] = {-1.0, -3.0*vt, 0.};
  double upper[] = { 1.0,  3.0*vt, 5.0*mu0};
  int vdim = 2;
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int cdim = ndim - vdim;
  int dir = cdim-1; // Always applied in the parallel direction.

  double lower_conf[cdim], upper_conf[cdim];
  int cells_conf[cdim];
  for (int d=0; d<cdim; d++) {
    lower_conf[d] = lower[d];
    upper_conf[d] = upper[d];
    cells_conf[d] = pars->cells[d];
  }
  double lower_vel[vdim], upper_vel[vdim];
  int cells_vel[vdim];
  for (int d=0; d<vdim; d++) {
    lower_vel[d] = lower[cdim+d];
    upper_vel[d] = upper[cdim+d];
    cells_vel[d] = pars->cells[cdim+d];
  }

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, pars->cells);
  struct gkyl_rect_grid grid_conf;
  gkyl_rect_grid_init(&grid_conf, cdim, lower_conf, upper_conf, cells_conf);
  struct gkyl_rect_grid grid_vel;
  gkyl_rect_grid_init(&grid_vel, vdim, lower_vel, upper_vel, cells_vel);

  // Basis functions.
  struct gkyl_basis *basis;
  basis = use_gpu? gkyl_cart_modal_gkhybrid_cu_dev_new(cdim, vdim)
                  : gkyl_cart_modal_gkhybrid_new(cdim, vdim);

  struct gkyl_basis basis_ho;
  gkyl_cart_modal_gkhybrid(&basis_ho, cdim, vdim);

  struct gkyl_basis basis_conf;
  gkyl_cart_modal_serendip(&basis_conf, cdim, poly_order);

  int ghost_conf[] = { 1, 1, 1 };
  struct gkyl_range local_conf, local_conf_ext; // local, local-ext position-space ranges
  gkyl_create_grid_ranges(&grid_conf, ghost_conf, &local_conf_ext, &local_conf);

  int ghost_vel[] = { 0, 0 };
  struct gkyl_range local_vel, local_vel_ext; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&grid_vel, ghost_vel, &local_vel_ext, &local_vel);

  int ghost[GKYL_MAX_DIM];
  for (int d=0; d<cdim; d++)
    ghost[d] = ghost_conf[d];
  for (int d=0; d<vdim; d++)
    ghost[cdim+d] = ghost_vel[d];
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Initialize velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, grid_vel,
    local, local_ext, local_vel, local_vel_ext, use_gpu);

  // Extended grid for the distribution function, which includes ghost cells.
  double lower_ext[ndim], upper_ext[ndim];
  int cells_ext[ndim];
  for (int d=0; d<ndim; d++) {
    double dx = (upper[d]-lower[d])/pars->cells[d];
    lower_ext[d] = lower[d]-dx*ghost[d];
    upper_ext[d] = upper[d]+dx*ghost[d];
    cells_ext[d] = pars->cells[d]+2*ghost[d];
  }
  struct gkyl_rect_grid grid_ext;
  gkyl_rect_grid_init(&grid_ext, ndim, lower_ext, upper_ext, cells_ext);

  // Create the skin/ghost ranges.
  struct gkyl_range skin_r, ghost_r;
  gkyl_skin_ghost_ranges(&skin_r, &ghost_r, dir, edge, &local_ext, ghost);

  // Initialize the distribution
  struct gkyl_array *distf = mkarr(use_gpu, basis_ho.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu? mkarr(false, distf->ncomp, distf->size) : gkyl_array_acquire(distf);
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = &grid,
      .basis = &basis_ho,
      .num_ret_vals = 1,
      .eval = eval_func_1x2v,
      .ctx = pars,
    }
  );
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf_ho);
  gkyl_array_copy(distf, distf_ho);

  // Initialize the electrostatic potential at MPE and WALL
  double dgnormc = pow(sqrt(2.0), cdim);
  struct gkyl_array *phi = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
  struct gkyl_array *phi_ho = use_gpu? mkarr(false, phi->ncomp, phi->size) : gkyl_array_acquire(phi);
  gkyl_array_shiftc(phi_ho, phi_mpe * dgnormc, 0 * basis_conf.num_basis);
  gkyl_array_copy(phi, phi_ho);

  struct gkyl_array *phiw = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
  struct gkyl_array *phiw_ho = use_gpu? mkarr(false, phiw->ncomp, phiw->size) : gkyl_array_acquire(phiw);
  gkyl_array_shiftc(phiw_ho, phi_wall * dgnormc, 0 * basis_conf.num_basis);
  gkyl_array_copy(phiw, phiw_ho);

  // Create the BC updater.
  struct gkyl_bc_sheath_gyrokinetic *bcsheath = gkyl_bc_sheath_gyrokinetic_new(dir, edge,
    basis, &skin_r, &ghost_r, gvm, cdim, 2.*qs/ms, pars->use_surrogate, pars->surrogate_model_path, use_gpu);

  // Build the vcut_fact DG array to make vpar cut vary.
  struct gkyl_basis vcut_fact_basis = gkyl_bc_sheath_gyrokinetic_get_vcut_fact_basis(bcsheath);
  struct gkyl_range *vcut_fact_local = gkyl_bc_sheath_gyrokinetic_get_vcut_fact_range(bcsheath);
  struct gkyl_array *vcut_fact = mkarr(use_gpu, vcut_fact_basis.num_basis, vcut_fact_local->volume);
  struct gkyl_array *vcut_fact_ho = use_gpu? mkarr(false, vcut_fact->ncomp, vcut_fact->size) : gkyl_array_acquire(vcut_fact);

  struct gkyl_rect_grid vcut_grid;
  double vcut_lower[cdim], vcut_upper[cdim];
  for (int d=0; d<cdim-1; d++) {
    vcut_lower[d] = lower[d];
    vcut_upper[d] = upper[d];
  }
  vcut_lower[cdim-1] = lower[cdim+vdim-1];
  vcut_upper[cdim-1] = upper[cdim+vdim-1];
  gkyl_rect_grid_init(&vcut_grid, cdim, vcut_lower, vcut_upper, pars->cells);
  gkyl_proj_on_basis *projVcut = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = &vcut_grid,
      .basis = &vcut_fact_basis,
      .num_ret_vals = 1,
      .eval = eval_func_vcut_fact,
      .ctx = pars,
    }
  );
  gkyl_proj_on_basis_advance(projVcut, 0.0, vcut_fact_local, vcut_fact_ho);
  gkyl_array_copy(vcut_fact, vcut_fact_ho);

  // Update the BC updater with vcut_fact array.
  gkyl_bc_sheath_gyrokinetic_set_vcut_fact(bcsheath, vcut_fact);

  if (pars->use_surrogate) {
    struct gkyl_array *density = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
    struct gkyl_array *density_ho = use_gpu? mkarr(false, density->ncomp, density->size) : gkyl_array_acquire(density);
    gkyl_array_shiftc(density_ho, ns * dgnormc, 0 * basis_conf.num_basis);
    gkyl_array_copy(density, density_ho);

    struct gkyl_array *temperature = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
    struct gkyl_array *temperature_ho = use_gpu? mkarr(false, temperature->ncomp, temperature->size) : gkyl_array_acquire(temperature);
    gkyl_array_shiftc(temperature_ho, Ts * dgnormc, 0 * basis_conf.num_basis);
    gkyl_array_copy(temperature, temperature_ho);

    struct gkyl_array *bmag = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
    struct gkyl_array *bmag_ho = use_gpu? mkarr(false, bmag->ncomp, bmag->size) : gkyl_array_acquire(bmag);
    gkyl_array_shiftc(bmag_ho, B0 * dgnormc, 0 * basis_conf.num_basis);
    gkyl_array_copy(bmag, bmag_ho);

    struct gkyl_array *bimpact_angle = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
    struct gkyl_array *bimpact_angle_ho = use_gpu? mkarr(false, bimpact_angle->ncomp, bimpact_angle->size) : gkyl_array_acquire(bimpact_angle);
    gkyl_array_shiftc(bimpact_angle_ho, impact_angle * dgnormc, 0 * basis_conf.num_basis);
    gkyl_array_copy(bimpact_angle, bimpact_angle_ho);

    gkyl_bc_sheath_gyrokinetic_update_vcut_fact_surrogate(bcsheath, phi, phiw, density, temperature, bmag, bimpact_angle, &local_conf);

    gkyl_array_release(density);
    gkyl_array_release(density_ho);
    gkyl_array_release(temperature);
    gkyl_array_release(temperature_ho);
    gkyl_array_release(bmag);
    gkyl_array_release(bmag_ho);
    gkyl_array_release(bimpact_angle);
    gkyl_array_release(bimpact_angle_ho);
  }

  // Advance the BC updater.
  gkyl_bc_sheath_gyrokinetic_advance(bcsheath, phi, phiw, distf, &local_conf);

  // Copy back to host.
  gkyl_array_copy(distf_ho, distf);

  // Check the result.
  check_function(pars, distf_ho, grid, ghost_r, edge);

  // Write out the distribution function after applying BC if requested.
  if (write_fields)
    write_out_fields(pars, cdim, vdim, edge, use_gpu, distf_ho, phi_ho, phiw_ho, 
      grid_ext, local_ext, grid_conf, local_conf, &basis_ho, basis_conf);

  // Clean up.
  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
  gkyl_array_release(distf_ho);
  gkyl_array_release(phi);
  gkyl_array_release(phi_ho);
  gkyl_array_release(phiw);
  gkyl_array_release(phiw_ho);
  if (use_gpu)
    gkyl_cart_modal_basis_release_cu(basis);
  else
    gkyl_cart_modal_basis_release(basis);

  gkyl_velocity_map_release(gvm);
  gkyl_bc_sheath_gyrokinetic_release(bcsheath);
  gkyl_proj_on_basis_release(projVcut);
  gkyl_array_release(vcut_fact);
  gkyl_array_release(vcut_fact_ho);
}

void
test_bc_sheath_gyrokinetic_2x2v(struct test_sheath_ctx *pars, enum gkyl_edge_loc edge, bool write_fields, bool use_gpu)
{
  /*
  This test applies the sheath BC to a Maxwellian distribution function on both upper and lower edges in the z direction,
  according to a potential drop between the magnetic presheath entrance and the wall.
  The following files are written out if requested:
  - "bc_sheath_2x2v_distf_in.gkyl": Input distribution function.
  - "bc_sheath_2x2v_phi_mpe.gkyl": Potential at the magnetic presheath entrance (constant everywhere).
  - "bc_sheath_2x2v_phi_wall.gkyl": Potential at the wall (constant everywhere).
  - "bc_sheath_2x2v_distf_out.gkyl": Distribution after the sheath BC has been applied (includes ghost cells).
  */

  double ms = pars->mass;
  double qs = pars->charge;
  double ns = pars->dens;
  double Ts = pars->temp;
  double B0 = pars->B0;
  double phi_mpe = pars->phi_mpe;
  double phi_wall = pars->phi_wall;
  double impact_angle = pars->impact_angle;
  double vt = sqrt(Ts/ms);
  double q2Dm = 2*qs/ms;
  double mu0 = Ts/B0;

  double phi_surr = -qs * (phi_mpe - phi_wall) / Ts;
  double gamma_surr = sqrt(ms * ns / GKYL_EPSILON0) / B0;
  double alpha_surr = impact_angle * 180.0 / GKYL_PI;

  // if (pars->use_surrogate)
    // printf("\n surrogate params: alpha = %.2f | gamma = %.2f | phi = %.2f\n", alpha_surr, gamma_surr, phi_surr);

  int poly_order = 1;
  double lower[] = { 0.0, -1.0, -3.0*vt, 0.};
  double upper[] = { 1.0,  1.0,  3.0*vt, 5*mu0};
  int vdim = 2;
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int cdim = ndim - vdim;
  int dir = cdim-1; // Always applied in the parallel direction.

  double lower_conf[cdim], upper_conf[cdim];
  int cells_conf[cdim];
  for (int d=0; d<cdim; d++) {
    lower_conf[d] = lower[d];
    upper_conf[d] = upper[d];
    cells_conf[d] = pars->cells[d];
  }
  double lower_vel[vdim], upper_vel[vdim];
  int cells_vel[vdim];
  for (int d=0; d<vdim; d++) {
    lower_vel[d] = lower[cdim+d];
    upper_vel[d] = upper[cdim+d];
    cells_vel[d] = pars->cells[cdim+d];
  }

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, pars->cells);
  struct gkyl_rect_grid grid_conf;
  gkyl_rect_grid_init(&grid_conf, cdim, lower_conf, upper_conf, cells_conf);
  struct gkyl_rect_grid grid_vel;
  gkyl_rect_grid_init(&grid_vel, vdim, lower_vel, upper_vel, cells_vel);

  // Basis functions.
  struct gkyl_basis *basis;
  basis = use_gpu? gkyl_cart_modal_gkhybrid_cu_dev_new(cdim, vdim)
                  : gkyl_cart_modal_gkhybrid_new(cdim, vdim);

  struct gkyl_basis basis_ho;
  gkyl_cart_modal_gkhybrid(&basis_ho, cdim, vdim);

  struct gkyl_basis basis_conf;
  gkyl_cart_modal_serendip(&basis_conf, cdim, poly_order);

  int ghost_conf[] = { 1, 1, 1 };
  struct gkyl_range local_conf, local_conf_ext; // local, local-ext position-space ranges
  gkyl_create_grid_ranges(&grid_conf, ghost_conf, &local_conf_ext, &local_conf);

  int ghost_vel[] = { 0, 0 };
  struct gkyl_range local_vel, local_vel_ext; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&grid_vel, ghost_vel, &local_vel_ext, &local_vel);

  int ghost[GKYL_MAX_DIM];
  for (int d=0; d<cdim; d++)
    ghost[d] = ghost_conf[d];
  for (int d=0; d<vdim; d++)
    ghost[cdim+d] = ghost_vel[d];
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Initialize velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, grid_vel,
    local, local_ext, local_vel, local_vel_ext, use_gpu);

  // Extended grid for the distribution function, which includes ghost cells.
  double lower_ext[ndim], upper_ext[ndim];
  int cells_ext[ndim];
  for (int d=0; d<ndim; d++) {
    double dx = (upper[d]-lower[d])/pars->cells[d];
    lower_ext[d] = lower[d]-dx*ghost[d];
    upper_ext[d] = upper[d]+dx*ghost[d];
    cells_ext[d] = pars->cells[d]+2*ghost[d];
  }
  struct gkyl_rect_grid grid_ext;
  gkyl_rect_grid_init(&grid_ext, ndim, lower_ext, upper_ext, cells_ext);

  // Create the skin/ghost ranges.
  struct gkyl_range skin_r, ghost_r;
  gkyl_skin_ghost_ranges(&skin_r, &ghost_r, dir, edge, &local_ext, ghost);

  // Initialize the distribution.
  struct gkyl_array *distf = mkarr(use_gpu, basis_ho.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu? mkarr(false, distf->ncomp, distf->size) : gkyl_array_acquire(distf);
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = &grid,
      .basis = &basis_ho,
      .num_ret_vals = 1,
      .eval = eval_func_2x2v,
      .ctx = pars,
    }
  );
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf_ho);
  gkyl_array_copy(distf, distf_ho);

  // Initialize the electrostatic potential at MPE and WALL.
  double dgnormc = pow(sqrt(2.0), cdim);
  struct gkyl_array *phi = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
  struct gkyl_array *phi_ho = use_gpu? mkarr(false, phi->ncomp, phi->size) : gkyl_array_acquire(phi);
  gkyl_array_shiftc(phi_ho, phi_mpe * dgnormc, 0 * basis_conf.num_basis);
  gkyl_array_copy(phi, phi_ho);

  struct gkyl_array *phiw = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
  struct gkyl_array *phiw_ho = use_gpu? mkarr(false, phiw->ncomp, phiw->size) : gkyl_array_acquire(phiw);
  gkyl_array_shiftc(phiw_ho, phi_wall * dgnormc, 0 * basis_conf.num_basis);
  gkyl_array_copy(phiw, phiw_ho);

  // Create the BC updater.
  struct gkyl_bc_sheath_gyrokinetic *bcsheath = gkyl_bc_sheath_gyrokinetic_new(dir, edge,
    basis, &skin_r, &ghost_r, gvm, cdim, 2.*qs/ms, pars->use_surrogate, pars->surrogate_model_path, use_gpu);

  // Build the vcut_fact DG array to make vpar cut vary.
  struct gkyl_basis vcut_fact_basis = gkyl_bc_sheath_gyrokinetic_get_vcut_fact_basis(bcsheath);
  struct gkyl_range *vcut_fact_local = gkyl_bc_sheath_gyrokinetic_get_vcut_fact_range(bcsheath);
  struct gkyl_array *vcut_fact = mkarr(use_gpu, vcut_fact_basis.num_basis, vcut_fact_local->volume);
  struct gkyl_array *vcut_fact_ho = use_gpu? mkarr(false, vcut_fact->ncomp, vcut_fact->size) : gkyl_array_acquire(vcut_fact);

  struct gkyl_rect_grid vcut_grid;
  double vcut_lower[cdim], vcut_upper[cdim];
  for (int d=0; d<cdim-1; d++) {
    vcut_lower[d] = lower[d];
    vcut_upper[d] = upper[d];
  }
  vcut_lower[cdim-1] = lower[cdim+vdim-1];
  vcut_upper[cdim-1] = upper[cdim+vdim-1];
  gkyl_rect_grid_init(&vcut_grid, cdim, vcut_lower, vcut_upper, pars->cells);
  gkyl_proj_on_basis *projVcut = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = &vcut_grid,
      .basis = &vcut_fact_basis,
      .num_ret_vals = 1,
      .eval = eval_func_vcut_fact,
      .ctx = pars,
    }
  );
  gkyl_proj_on_basis_advance(projVcut, 0.0, vcut_fact_local, vcut_fact_ho);
  gkyl_array_copy(vcut_fact, vcut_fact_ho);

  // Update the BC updater with vcut_fact array.
  gkyl_bc_sheath_gyrokinetic_set_vcut_fact(bcsheath, vcut_fact);

  if (pars->use_surrogate) {
    struct gkyl_array *density = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
    struct gkyl_array *density_ho = use_gpu? mkarr(false, density->ncomp, density->size) : gkyl_array_acquire(density);
    gkyl_array_shiftc(density_ho, ns * dgnormc, 0 * basis_conf.num_basis);
    gkyl_array_copy(density, density_ho);

    struct gkyl_array *temperature = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
    struct gkyl_array *temperature_ho = use_gpu? mkarr(false, temperature->ncomp, temperature->size) : gkyl_array_acquire(temperature);
    gkyl_array_shiftc(temperature_ho, Ts * dgnormc, 0 * basis_conf.num_basis);
    gkyl_array_copy(temperature, temperature_ho);

    struct gkyl_array *bmag = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
    struct gkyl_array *bmag_ho = use_gpu? mkarr(false, bmag->ncomp, bmag->size) : gkyl_array_acquire(bmag);
    gkyl_array_shiftc(bmag_ho, B0 * dgnormc, 0 * basis_conf.num_basis);
    gkyl_array_copy(bmag, bmag_ho);

    struct gkyl_array *bimpact_angle = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
    struct gkyl_array *bimpact_angle_ho = use_gpu? mkarr(false, bimpact_angle->ncomp, bimpact_angle->size) : gkyl_array_acquire(bimpact_angle);
    gkyl_array_shiftc(bimpact_angle_ho, impact_angle * dgnormc, 0 * basis_conf.num_basis);
    gkyl_array_copy(bimpact_angle, bimpact_angle_ho);

    gkyl_bc_sheath_gyrokinetic_update_vcut_fact_surrogate(bcsheath, phi, phiw, density, temperature, bmag, bimpact_angle, &local_conf);

    gkyl_array_release(density);
    gkyl_array_release(density_ho);
    gkyl_array_release(temperature);
    gkyl_array_release(temperature_ho);
    gkyl_array_release(bmag);
    gkyl_array_release(bmag_ho);
    gkyl_array_release(bimpact_angle);
    gkyl_array_release(bimpact_angle_ho);
  }

  // Advance the BC updater.
  gkyl_bc_sheath_gyrokinetic_advance(bcsheath, phi, phiw, distf, &local_conf);

  // Copy back to host.
  gkyl_array_copy(distf_ho, distf);

  // Check the result.
  check_function(pars, distf_ho, grid, ghost_r, edge);

  // Write out the distribution function after applying BC if requested.
  if (write_fields)
    write_out_fields(pars, cdim, vdim, edge, use_gpu, distf_ho, phi_ho, phiw_ho, 
      grid_ext, local_ext, grid_conf, local_conf, &basis_ho, basis_conf);

  // Clean up.
  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
  gkyl_array_release(distf_ho);
  gkyl_array_release(phi);
  gkyl_array_release(phi_ho);
  gkyl_array_release(phiw);
  gkyl_array_release(phiw_ho);
  if (use_gpu)
    gkyl_cart_modal_basis_release_cu(basis);
  else
    gkyl_cart_modal_basis_release(basis);

  gkyl_velocity_map_release(gvm);
  gkyl_bc_sheath_gyrokinetic_release(bcsheath);
  gkyl_proj_on_basis_release(projVcut);
  gkyl_array_release(vcut_fact);
  gkyl_array_release(vcut_fact_ho);
}

void
test_bc_sheath_gyrokinetic_3x2v(struct test_sheath_ctx *pars, enum gkyl_edge_loc edge, bool write_fields, bool use_gpu)
{
  /*
  This test applies the sheath BC to a Maxwellian distribution function on both upper and lower edges in the z direction,
  according to a potential drop between the magnetic presheath entrance and the wall.
  The following files are written out if requested:
  - "bc_sheath_3x2v_distf_in.gkyl": Input distribution function.
  - "bc_sheath_3x2v_phi_mpe.gkyl": Potential at the magnetic presheath entrance (constant everywhere).
  - "bc_sheath_3x2v_phi_wall.gkyl": Potential at the wall (constant everywhere).
  - "bc_sheath_3x2v_distf_out.gkyl": Distribution after the sheath BC has been applied (includes ghost cells).
  */

  double ms = pars->mass;
  double qs = pars->charge;
  double ns = pars->dens;
  double Ts = pars->temp;
  double B0 = pars->B0;
  double phi_mpe = pars->phi_mpe;
  double phi_wall = pars->phi_wall;
  double impact_angle = pars->impact_angle;
  double vt = sqrt(Ts/ms);
  double mu0 = Ts/B0;
  double q2Dm = 2*qs/ms;

  double phi_surr = -qs * (phi_mpe - phi_wall) / Ts;
  double gamma_surr = sqrt(ms * ns / GKYL_EPSILON0) / B0;
  double alpha_surr = impact_angle * 180.0 / GKYL_PI;

  // if (pars->use_surrogate)
  //  printf("\n surrogate params: alpha = %.2f | gamma = %.2f | phi = %.2f\n", alpha_surr, gamma_surr, phi_surr);

  int poly_order = 1;
  double lower[] = { 0.0, -2.0, -1.0, -3.0*vt, 0.};
  double upper[] = { 1.0,  2.0,  1.0,  3.0*vt, 5*mu0};
  int vdim = 2;
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int cdim = ndim - vdim;
  int dir = cdim-1; // Always applied in the parallel direction.
  double x0 = 0.0; // x-center of distribution function
  double y0 = 0.0; // y-center of distribution function
  double z0 = 0.2; // z-center of distribution function
  double sigmax = 0.5; // Width of distribution function in x.
  double sigmay = 0.5; // Width of distribution function in y.
  double sigmaz = 1.0; // Width of distribution function in z.

  double lower_conf[cdim], upper_conf[cdim];
  int cells_conf[cdim];
  for (int d=0; d<cdim; d++) {
    lower_conf[d] = lower[d];
    upper_conf[d] = upper[d];
    cells_conf[d] = pars->cells[d];
  }
  double lower_vel[vdim], upper_vel[vdim];
  int cells_vel[vdim];
  for (int d=0; d<vdim; d++) {
    lower_vel[d] = lower[cdim+d];
    upper_vel[d] = upper[cdim+d];
    cells_vel[d] = pars->cells[cdim+d];
  }

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, pars->cells);
  struct gkyl_rect_grid grid_conf;
  gkyl_rect_grid_init(&grid_conf, cdim, lower_conf, upper_conf, cells_conf);
  struct gkyl_rect_grid grid_vel;
  gkyl_rect_grid_init(&grid_vel, vdim, lower_vel, upper_vel, cells_vel);

  // Basis functions.
  struct gkyl_basis *basis;
  basis = use_gpu? gkyl_cart_modal_gkhybrid_cu_dev_new(cdim, vdim)
                  : gkyl_cart_modal_gkhybrid_new(cdim, vdim);

  struct gkyl_basis basis_ho;
  gkyl_cart_modal_gkhybrid(&basis_ho, cdim, vdim);

  struct gkyl_basis basis_conf;
  gkyl_cart_modal_serendip(&basis_conf, cdim, poly_order);

  int ghost_conf[] = { 1, 1, 1 };
  struct gkyl_range local_conf, local_conf_ext; // local, local-ext position-space ranges
  gkyl_create_grid_ranges(&grid_conf, ghost_conf, &local_conf_ext, &local_conf);

  int ghost_vel[] = { 0, 0 };
  struct gkyl_range local_vel, local_vel_ext; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&grid_vel, ghost_vel, &local_vel_ext, &local_vel);

  int ghost[GKYL_MAX_DIM];
  for (int d=0; d<cdim; d++)
    ghost[d] = ghost_conf[d];
  for (int d=0; d<vdim; d++)
    ghost[cdim+d] = ghost_vel[d];
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Initialize velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, grid_vel,
    local, local_ext, local_vel, local_vel_ext, use_gpu);

  // Extended grid for the distribution function, which includes ghost cells.
  double lower_ext[ndim], upper_ext[ndim];
  int cells_ext[ndim];
  for (int d=0; d<ndim; d++) {
    double dx = (upper[d]-lower[d])/pars->cells[d];
    lower_ext[d] = lower[d]-dx*ghost[d];
    upper_ext[d] = upper[d]+dx*ghost[d];
    cells_ext[d] = pars->cells[d]+2*ghost[d];
  }
  struct gkyl_rect_grid grid_ext;
  gkyl_rect_grid_init(&grid_ext, ndim, lower_ext, upper_ext, cells_ext);

  // Create the skin/ghost ranges.
  struct gkyl_range skin_r, ghost_r;
  gkyl_skin_ghost_ranges(&skin_r, &ghost_r, dir, edge, &local_ext, ghost);

  // Initialize the distribution.
  struct gkyl_array *distf = mkarr(use_gpu, basis_ho.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu? mkarr(false, distf->ncomp, distf->size) : gkyl_array_acquire(distf);
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = &grid,
      .basis = &basis_ho,
      .num_ret_vals = 1,
      .eval = eval_func_3x2v,
      .ctx = pars,
    }
  );
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf_ho);
  gkyl_array_copy(distf, distf_ho);

  // Initialize the electrostatic potential at MPE and WALL.
  double dgnormc = pow(sqrt(2.0), cdim);
  struct gkyl_array *phi = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
  struct gkyl_array *phi_ho = use_gpu? mkarr(false, phi->ncomp, phi->size) : gkyl_array_acquire(phi);
  gkyl_array_shiftc(phi_ho, phi_mpe * dgnormc, 0 * basis_conf.num_basis);
  gkyl_array_copy(phi, phi_ho);

  struct gkyl_array *phiw = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
  struct gkyl_array *phiw_ho = use_gpu? mkarr(false, phiw->ncomp, phiw->size) : gkyl_array_acquire(phiw);
  gkyl_array_shiftc(phiw_ho, phi_wall * dgnormc, 0 * basis_conf.num_basis);
  gkyl_array_copy(phiw, phiw_ho);

  // Create the BC updater.
  struct gkyl_bc_sheath_gyrokinetic *bcsheath = gkyl_bc_sheath_gyrokinetic_new(dir, edge,
    basis, &skin_r, &ghost_r, gvm, cdim, 2.*qs/ms, pars->use_surrogate, pars->surrogate_model_path, use_gpu);

  // Build the vcut_fact DG array to make vpar cut vary.
  struct gkyl_basis vcut_fact_basis = gkyl_bc_sheath_gyrokinetic_get_vcut_fact_basis(bcsheath);
  struct gkyl_range *vcut_fact_local = gkyl_bc_sheath_gyrokinetic_get_vcut_fact_range(bcsheath);
  struct gkyl_array *vcut_fact = mkarr(use_gpu, vcut_fact_basis.num_basis, vcut_fact_local->volume);
  struct gkyl_array *vcut_fact_ho = use_gpu? mkarr(false, vcut_fact->ncomp, vcut_fact->size) : gkyl_array_acquire(vcut_fact);

  struct gkyl_rect_grid vcut_grid;
  double vcut_lower[cdim], vcut_upper[cdim];
  for (int d=0; d<cdim-1; d++) {
    vcut_lower[d] = lower[d];
    vcut_upper[d] = upper[d];
  }
  vcut_lower[cdim-1] = lower[cdim+vdim-1];
  vcut_upper[cdim-1] = upper[cdim+vdim-1];
  gkyl_rect_grid_init(&vcut_grid, cdim, vcut_lower, vcut_upper, pars->cells);
  gkyl_proj_on_basis *projVcut = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = &vcut_grid,
      .basis = &vcut_fact_basis,
      .num_ret_vals = 1,
      .eval = eval_func_vcut_fact,
      .ctx = pars,
    }
  );
  gkyl_proj_on_basis_advance(projVcut, 0.0, vcut_fact_local, vcut_fact_ho);
  gkyl_array_copy(vcut_fact, vcut_fact_ho);

  // Update the BC updater with vcut_fact array.
  gkyl_bc_sheath_gyrokinetic_set_vcut_fact(bcsheath, vcut_fact);

  if (pars->use_surrogate) {
    struct gkyl_array *density = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
    struct gkyl_array *density_ho = use_gpu? mkarr(false, density->ncomp, density->size) : gkyl_array_acquire(density);
    gkyl_array_shiftc(density_ho, ns * dgnormc, 0 * basis_conf.num_basis);
    gkyl_array_copy(density, density_ho);

    struct gkyl_array *temperature = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
    struct gkyl_array *temperature_ho = use_gpu? mkarr(false, temperature->ncomp, temperature->size) : gkyl_array_acquire(temperature);
    gkyl_array_shiftc(temperature_ho, Ts * dgnormc, 0 * basis_conf.num_basis);
    gkyl_array_copy(temperature, temperature_ho);

    struct gkyl_array *bmag = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
    struct gkyl_array *bmag_ho = use_gpu? mkarr(false, bmag->ncomp, bmag->size) : gkyl_array_acquire(bmag);
    gkyl_array_shiftc(bmag_ho, B0 * dgnormc, 0 * basis_conf.num_basis);
    gkyl_array_copy(bmag, bmag_ho);

    struct gkyl_array *bimpact_angle = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
    struct gkyl_array *bimpact_angle_ho = use_gpu? mkarr(false, bimpact_angle->ncomp, bimpact_angle->size) : gkyl_array_acquire(bimpact_angle);
    gkyl_array_shiftc(bimpact_angle_ho, impact_angle * dgnormc, 0 * basis_conf.num_basis);
    gkyl_array_copy(bimpact_angle, bimpact_angle_ho);

    gkyl_bc_sheath_gyrokinetic_update_vcut_fact_surrogate(bcsheath, phi, phiw, density, temperature, bmag, bimpact_angle, &local_conf);

    gkyl_array_release(density);
    gkyl_array_release(density_ho);
    gkyl_array_release(temperature);
    gkyl_array_release(temperature_ho);
    gkyl_array_release(bmag);
    gkyl_array_release(bmag_ho);
    gkyl_array_release(bimpact_angle);
    gkyl_array_release(bimpact_angle_ho);
  }

  // Advance the BC updater.
  gkyl_bc_sheath_gyrokinetic_advance(bcsheath, phi, phiw, distf, &local_conf);

  // Copy back to host.
  gkyl_array_copy(distf_ho, distf);

  // Check the result.
  check_function(pars, distf_ho, grid, ghost_r, edge);

  // Write out the distribution function after applying BC if requested.
  if (write_fields)
    write_out_fields(pars, cdim, vdim, edge, use_gpu, distf_ho, phi_ho, phiw_ho, 
      grid_ext, local_ext, grid_conf, local_conf, &basis_ho, basis_conf);

  // Clean up.
  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
  gkyl_array_release(distf_ho);
  gkyl_array_release(phi);
  gkyl_array_release(phi_ho);
  gkyl_array_release(phiw);
  gkyl_array_release(phiw_ho);
  if (use_gpu)
    gkyl_cart_modal_basis_release_cu(basis);
  else
    gkyl_cart_modal_basis_release(basis);

  gkyl_velocity_map_release(gvm);
  gkyl_bc_sheath_gyrokinetic_release(bcsheath);
  gkyl_proj_on_basis_release(projVcut);
  gkyl_array_release(vcut_fact);
  gkyl_array_release(vcut_fact_ho);
}

void test_bc_sheath_gk_1x2v_ho()
{
  if (!srgrz_test_enabled || !srgrz_model_path) return;
  bool write_fields;

  struct test_sheath_ctx pars = {
    .cdim = 1,
    .cells = {4, 16, 12},
    .mass = GKYL_ELECTRON_MASS,
    .charge = -GKYL_ELEMENTARY_CHARGE,
    .dens = 1.0e19,
    .upar = 0.0,
    .temp = 5.0 * GKYL_ELEMENTARY_CHARGE,
    .B0 = 1.0,
    .phi_mpe = 15.0,
    .phi_wall = 0.0,
    .impact_angle = 5.0 * GKYL_PI / 180.0, // 5 degrees
    .z0 = 0.0,
    .sigmaz = 1.0,
    .use_surrogate = true,
    .surrogate_model_path = srgrz_model_path,
    .verbose = false,
  };

  // Electrons with positive sheath entrance potential.
  write_fields = false;
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_LOWER_EDGE, write_fields, false);
  pars.use_surrogate = false;
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_UPPER_EDGE, write_fields, false);

  // // // Reverse potential sign.
  write_fields = false;
  pars.phi_mpe *= -1.0;
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_LOWER_EDGE, write_fields, false);
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_UPPER_EDGE, write_fields, false);

  // // // Ions with negative sheath entrance potential.
  pars.use_surrogate = false; // the surrogate is for electrons.
  write_fields = false;
  pars.charge = GKYL_ELEMENTARY_CHARGE;
  pars.mass = GKYL_PROTON_MASS;
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_LOWER_EDGE, write_fields, false);
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_UPPER_EDGE, write_fields, false);

  // // // // Reverse potential sign.
  write_fields = false;
  pars.phi_mpe *= -1.0;
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_LOWER_EDGE, write_fields, false);
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_UPPER_EDGE, write_fields, false);
}

void test_bc_sheath_gk_2x2v_ho()
{
  bool write_fields;

  struct test_sheath_ctx pars = {
    .cdim = 2,
    .cells = {4, 4, 16, 12},
    .mass = GKYL_ELECTRON_MASS,
    .charge = -GKYL_ELEMENTARY_CHARGE,
    .dens = 1.0e19,
    .upar = 0.0,
    .temp = 5.0 * GKYL_ELEMENTARY_CHARGE,
    .B0 = 1.0,
    .phi_mpe = 15.0,
    .phi_wall = 0.0,
    .impact_angle = 5.0 * GKYL_PI / 180.0, // 5 degrees
    .x0 = 0.0,
    .sigmax = 1.0,
    .z0 = 0.0,
    .sigmaz = 1.0,
    .use_surrogate = false,
    .verbose = false,
  };

  // Electrons with positive sheath entrance potential.
  write_fields = false;
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_LOWER_EDGE, write_fields, false);
  pars.use_surrogate = false;
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_UPPER_EDGE, write_fields, false);

  // // Electrons with negative sheath entrance potential.
  pars.phi_mpe *= -1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_LOWER_EDGE, write_fields, false);
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_UPPER_EDGE, write_fields, false);

  // // Ions with negative sheath entrance potential.
  pars.use_surrogate = false; // the surrogate is for electrons.
  write_fields = false;
  pars.charge = GKYL_ELEMENTARY_CHARGE;
  pars.mass = GKYL_PROTON_MASS;
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_LOWER_EDGE, write_fields, false);
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_UPPER_EDGE, write_fields, false);

  // // // Reverse potential sign.
  write_fields = false;
  pars.phi_mpe *= -1.0;
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_LOWER_EDGE, write_fields, false);
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_UPPER_EDGE, write_fields, false);
}
void test_bc_sheath_gk_3x2v_ho()
{
  if (!srgrz_test_enabled || !srgrz_model_path) return;
  bool write_fields;

  struct test_sheath_ctx pars = {
    .cdim = 3,
    .cells = {4, 4, 4, 16, 12},
    .mass = GKYL_ELECTRON_MASS,
    .charge = -GKYL_ELEMENTARY_CHARGE,
    .dens = 1.0e19,
    .upar = 0.0,
    .temp = 5.0 * GKYL_ELEMENTARY_CHARGE,
    .B0 = 1.0,
    .phi_mpe = 15.0,
    .phi_wall = 0.0,
    .impact_angle = 5.0 * GKYL_PI / 180.0, // 5 degrees
    .x0 = 0.0,
    .sigmax = 1.0,
    .y0 = 0.0,
    .sigmay = 1.0,
    .z0 = 0.0,
    .sigmaz = 1.0,
    .use_surrogate = true,
    .surrogate_model_path = srgrz_model_path,
    .verbose = false,
  };

  // Electrons with positive sheath entrance potential.
  write_fields = false;
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_LOWER_EDGE, write_fields, false);
  pars.use_surrogate = false;
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_UPPER_EDGE, write_fields, false);

  // // Electrons with negative sheath entrance potential.
  write_fields = false;
  pars.phi_mpe *= -1.0;
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_LOWER_EDGE, write_fields, false);
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_UPPER_EDGE, write_fields, false);

  // // Ions with negative sheath entrance potential.
  pars.use_surrogate = false; // the surrogate is for electrons.
  write_fields = false;
  pars.charge = GKYL_ELEMENTARY_CHARGE;
  pars.mass = GKYL_PROTON_MASS;
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_LOWER_EDGE, write_fields, false);
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_UPPER_EDGE, write_fields, false);

  // // // Reverse potential sign.
  write_fields = false;
  pars.phi_mpe *= -1.0;
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_LOWER_EDGE, write_fields, false);
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_UPPER_EDGE, write_fields, false);
}

#ifdef GKYL_HAVE_CUDA
void test_bc_sheath_gk_1x2v_dev()
{
  bool write_fields;
  
  struct test_sheath_ctx pars = {
    .cdim = 1,
    .cells = {4, 16, 12},
    .mass = GKYL_ELECTRON_MASS,
    .charge = -GKYL_ELEMENTARY_CHARGE,
    .dens = 1.0e19,
    .upar = 0.0,
    .temp = 5.0 * GKYL_ELEMENTARY_CHARGE,
    .B0 = 1.0,
    .phi_mpe = 15.0,
    .phi_wall = 0.0,
    .impact_angle = 5.0 * GKYL_PI / 180.0, // 5 degrees
    .z0 = 0.0,
    .sigmaz = 1.0,
    .use_surrogate = true,
    .verbose = false,
  };

  // Electrons with positive sheath entrance potential.
  write_fields = false;
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_LOWER_EDGE, write_fields, true);
  pars.use_surrogate = false;
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_UPPER_EDGE, write_fields, true);

  // Electrons with negative sheath entrance potential.
  pars.phi_mpe *= -1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_LOWER_EDGE, write_fields, true);
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_UPPER_EDGE, write_fields, true);

  // Ions with negative sheath entrance potential.
  pars.charge = GKYL_ELEMENTARY_CHARGE;
  pars.mass = GKYL_PROTON_MASS;
  write_fields = false;
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_LOWER_EDGE, write_fields, true);
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_UPPER_EDGE, write_fields, true);

  // Ions with positive sheath entrance potential.
  pars.phi_mpe *= -1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_LOWER_EDGE, write_fields, true);
  test_bc_sheath_gyrokinetic_1x2v(&pars, GKYL_UPPER_EDGE, write_fields, true);
}

void test_bc_sheath_gk_2x2v_dev()
{
  bool write_fields;

  struct test_sheath_ctx pars = {
    .cdim = 2,
    .cells = {4, 4, 16, 12},
    .mass = GKYL_ELECTRON_MASS,
    .charge = -GKYL_ELEMENTARY_CHARGE,
    .dens = 1.0e19,
    .upar = 0.0,
    .temp = 5.0 * GKYL_ELEMENTARY_CHARGE,
    .B0 = 1.0,
    .phi_mpe = 15.0,
    .phi_wall = 0.0,
    .impact_angle = 5.0 * GKYL_PI / 180.0, // 5 degrees
    .x0 = 0.0,
    .sigmax = 1.0,
    .z0 = 0.0,
    .sigmaz = 1.0,
    .use_surrogate = true,
    .verbose = false,
  };

  // Electrons with positive sheath entrance potential.
  write_fields = false;
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_LOWER_EDGE, write_fields, true);
  pars.use_surrogate = false;
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_UPPER_EDGE, write_fields, true);

  // Electrons with negative sheath entrance potential.
  pars.phi_mpe *= -1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_LOWER_EDGE, write_fields, true);
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_UPPER_EDGE, write_fields, true);

  // Ions with negative sheath entrance potential.
  pars.charge = GKYL_ELEMENTARY_CHARGE;
  pars.mass = GKYL_PROTON_MASS;
  write_fields = false;
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_LOWER_EDGE, write_fields, true);
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_UPPER_EDGE, write_fields, true);

  // // Ions with positive sheath entrance potential.
  pars.phi_mpe *= -1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_LOWER_EDGE, write_fields, true);
  test_bc_sheath_gyrokinetic_2x2v(&pars, GKYL_UPPER_EDGE, write_fields, true);
}

void test_bc_sheath_gk_3x2v_dev()
{
  bool write_fields;

  struct test_sheath_ctx pars = {
    .cdim = 3,
    .cells = {4, 4, 4, 16, 12},
    .mass = GKYL_ELECTRON_MASS,
    .charge = -GKYL_ELEMENTARY_CHARGE,
    .dens = 1.0e19,
    .upar = 0.0,
    .temp = 5.0 * GKYL_ELEMENTARY_CHARGE,
    .B0 = 1.0,
    .phi_mpe = 15.0,
    .phi_wall = 0.0,
    .impact_angle = 5.0 * GKYL_PI / 180.0, // 5 degrees
    .x0 = 0.0,
    .sigmax = 1.0,
    .y0 = 0.0,
    .sigmay = 1.0,
    .z0 = 0.0,
    .sigmaz = 1.0,
    .use_surrogate = true,
    .verbose = false,
  };

  // Electrons with positive sheath entrance potential.
  write_fields = false;
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_LOWER_EDGE, write_fields, true);
  pars.use_surrogate = false;
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_UPPER_EDGE, write_fields, true);

  // Electrons with negative sheath entrance potential.
  pars.phi_mpe *= -1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_LOWER_EDGE, write_fields, true);
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_UPPER_EDGE, write_fields, true);

  // Ions with negative sheath entrance potential.
  pars.charge = GKYL_ELEMENTARY_CHARGE;
  pars.mass = GKYL_PROTON_MASS;
  write_fields = false;
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_LOWER_EDGE, write_fields, true);
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_UPPER_EDGE, write_fields, true);

  // Ions with positive sheath entrance potential.
  pars.phi_mpe *= -1.0;
  write_fields = false;
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_LOWER_EDGE, write_fields, true);
  test_bc_sheath_gyrokinetic_3x2v(&pars, GKYL_UPPER_EDGE, write_fields, true);
}
#endif

TEST_LIST = {
  { "test_bc_sheath_gk_1x2v_ho", test_bc_sheath_gk_1x2v_ho },
  { "test_bc_sheath_gk_2x2v_ho", test_bc_sheath_gk_2x2v_ho },
  { "test_bc_sheath_gk_3x2v_ho", test_bc_sheath_gk_3x2v_ho },
#ifdef GKYL_HAVE_CUDA
  { "test_bc_sheath_gk_1x2v_dev", test_bc_sheath_gk_1x2v_dev },
  { "test_bc_sheath_gk_2x2v_dev", test_bc_sheath_gk_2x2v_dev },
  { "test_bc_sheath_gk_3x2v_dev", test_bc_sheath_gk_3x2v_dev },
#endif
  { NULL, NULL },
};