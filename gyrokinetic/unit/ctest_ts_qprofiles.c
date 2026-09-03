// Test creation and deallocation of updater that applies the
// twist shift BCs.
//
#include <acutest.h>

#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_basis.h>
#include <gkyl_proj_on_basis.h>
#include <mpack.h>
#include <gkyl_array_rio.h>
#include <gkyl_bc_twistshift.h>
#include <gkyl_velocity_map.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_dg_updater_moment_gyrokinetic.h>

// test function
void evalFunc(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  double n = 2.0;
  fout[0] = 1.0 + 0*sin(2.0*M_PI*n*y/0.2);
}

// Meta-data for IO
struct test_bc_twistshift_output_meta {
  int poly_order; // polynomial order
  const char *basis_type; // name of basis functions
};

// returned gkyl_array_meta must be freed using gyrokinetic_array_meta_release
static struct gkyl_msgpack_data*
test_bc_twistshift_array_meta_new(struct test_bc_twistshift_output_meta meta)
{
  struct gkyl_msgpack_data *mt = gkyl_malloc(sizeof(*mt));

  mt->meta_sz = 0;
  mpack_writer_t writer;
  mpack_writer_init_growable(&writer, &mt->meta, &mt->meta_sz);

  // add some data to mpack
  mpack_build_map(&writer);

  mpack_write_cstr(&writer, "polyOrder");
  mpack_write_i64(&writer, meta.poly_order);

  mpack_write_cstr(&writer, "basisType");
  mpack_write_cstr(&writer, meta.basis_type);

  mpack_complete_map(&writer);

  int status = mpack_writer_destroy(&writer);

  if (status != mpack_ok) {
    free(mt->meta); // we need to use free here as mpack does its own malloc
    gkyl_free(mt);
    mt = 0;
  }

  return mt;
}

static void
test_bc_twistshift_array_meta_release(struct gkyl_msgpack_data *mt)
{
  if (!mt) return;
  MPACK_FREE(mt->meta);
  gkyl_free(mt);
}

struct skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

// Create ghost and skin sub-ranges given a parent range
static void
skin_ghost_ranges_init(struct skin_ghost_ranges *sgr,
  const struct gkyl_range *parent, const int *ghost)
{
  int ndim = parent->ndim;

  for (int d=0; d<ndim; ++d) {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
      d, GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
      d, GKYL_UPPER_EDGE, parent, ghost);
  }
}
// Apply periodic BCs along parallel direction
void
apply_periodic_bc(struct gkyl_array *buff, struct gkyl_array *fld, const int dir, const struct skin_ghost_ranges sgr)
{
  gkyl_array_copy_to_buffer(buff->data, fld, &(sgr.lower_skin[dir]));
  gkyl_array_copy_from_buffer(fld, buff->data, &(sgr.upper_ghost[dir]));

  gkyl_array_copy_to_buffer(buff->data, fld, &(sgr.upper_skin[dir]));
  gkyl_array_copy_from_buffer(fld, buff->data, &(sgr.lower_ghost[dir]));
}

static struct gkyl_array*
mkarr(bool on_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (on_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct test_bc_twistshift_ctx {
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  int cells[GKYL_MAX_DIM];
  double B0;
  double vt;
  double mass;
};

void
mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; xp[1] = xc[1]; xp[2] = xc[2];
}

double r_x(double x, double a_mid, double x_inner)
{
  return x+a_mid-x_inner;
}

// /*
double qprofile_TCV(double r, double R_axis) 
{
  // Magnetic safety factor as a function of minor radius r.
  double a[4] = {484.0615913225881, -1378.25993228584, 1309.3099150729233, 
    -414.13270311478726};
  return a[0]*pow(r+R_axis,3.0) + a[1]*pow(r+R_axis,2.0) + a[2]*(r+R_axis) + a[3];
}
// */

// /*
double qprofile_lowshear(double r, double R_axis) 
{
  return 5*pow(r+R_axis,1.0) + 5;
}
// */

// /*
double qprofile_highshear(double r, double R_axis) 
{
  double a[4] = {484.0615913225881, -1378.25993228584, 1309.3099150729233, 
  -414.13270311478726};
  return (a[0]*pow(r+R_axis,3.0) + a[1]*pow(r+R_axis,2.0) + a[2]*(r+R_axis) + a[3]) 
    * (20*a[0]*pow(r+R_axis,2.0) + 2*a[1]*pow(r+R_axis,1.0) + a[2])/30000;
}
// */

double qprofile_linear(double r, double R_axis) 
{
  double a = 44.61713037036269;
  double b = -46.18401975122209;
  return a*(r+R_axis) + b;
}

double qprofile_quadratic(double r, double R_axis) 
{
  double a = 247.89266869551582;
  double b = -510.5619559135587;
  double c = 264.34987390804827;
  return a*pow(r+R_axis,2.0) + b*(r+R_axis) + c;
}

double qprofile_cubic(double r, double R_axis) 
{
  double a = 484.0615913225881;
  double b = -1378.25993228584;
  double c = 1309.3099150729233;
  double d = -414.13270311478726;
  return a*pow(r+R_axis,3.0) + b*pow(r+R_axis,2.0) + c*(r+R_axis) + d;
}


double qprofile_pwlin8(double r, double R_axis) {
  double R = r + R_axis;
  double q = 0.0;
  if (R <= 1.0748) q = 21.5288 * R + -21.172;
  if (R >= 1.0748 && R <= 1.0898) q = 27.0051 * R + -27.0579;
  if (R >= 1.0898 && R <= 1.1048) q = 33.1349 * R + -33.7382;
  if (R >= 1.1048 && R <= 1.1198) q = 39.9182 * R + -41.2323;
  if (R >= 1.1198 && R <= 1.1348) q = 47.355 * R + -49.56;
  if (R >= 1.1348 && R <= 1.1498) q = 55.4453 * R + -58.7408;
  if (R >= 1.1498 && R <= 1.1648) q = 64.189 * R + -68.7944;
  if (R >= 1.1648) q = 73.5862 * R + -79.7402;
  return q;
}

double qprofile_pwlin48(double r, double R_axis) {
  double R = r + R_axis;
  double q = 0.0;
  if (R <= 1.0623) q = 19.4134 * R + -18.93;
  if (R >= 1.0623 && R <= 1.0648) q = 20.2353 * R + -19.8032;
  if (R >= 1.0648 && R <= 1.0673) q = 21.0754 * R + -20.6978;
  if (R >= 1.0673 && R <= 1.0698) q = 21.9337 * R + -21.6138;
  if (R >= 1.0698 && R <= 1.0723) q = 22.8101 * R + -22.5514;
  if (R >= 1.0723 && R <= 1.0748) q = 23.7047 * R + -23.5106;
  if (R >= 1.0748 && R <= 1.0773) q = 24.6174 * R + -24.4916;
  if (R >= 1.0773 && R <= 1.0798) q = 25.5483 * R + -25.4944;
  if (R >= 1.0798 && R <= 1.0823) q = 26.4973 * R + -26.5192;
  if (R >= 1.0823 && R <= 1.0848) q = 27.4645 * R + -27.566;
  if (R >= 1.0848 && R <= 1.0873) q = 28.4498 * R + -28.6349;
  if (R >= 1.0873 && R <= 1.0898) q = 29.4533 * R + -29.7259;
  if (R >= 1.0898 && R <= 1.0923) q = 30.4749 * R + -30.8393;
  if (R >= 1.0923 && R <= 1.0948) q = 31.5147 * R + -31.9751;
  if (R >= 1.0948 && R <= 1.0973) q = 32.5727 * R + -33.1333;
  if (R >= 1.0973 && R <= 1.0998) q = 33.6488 * R + -34.3141;
  if (R >= 1.0998 && R <= 1.1023) q = 34.743 * R + -35.5175;
  if (R >= 1.1023 && R <= 1.1048) q = 35.8554 * R + -36.7437;
  if (R >= 1.1048 && R <= 1.1073) q = 36.986 * R + -37.9928;
  if (R >= 1.1073 && R <= 1.1098) q = 38.1347 * R + -39.2647;
  if (R >= 1.1098 && R <= 1.1123) q = 39.3015 * R + -40.5597;
  if (R >= 1.1123 && R <= 1.1148) q = 40.4865 * R + -41.8778;
  if (R >= 1.1148 && R <= 1.1173) q = 41.6897 * R + -43.2191;
  if (R >= 1.1173 && R <= 1.1198) q = 42.911 * R + -44.5836;
  if (R >= 1.1198 && R <= 1.1223) q = 44.1505 * R + -45.9716;
  if (R >= 1.1223 && R <= 1.1248) q = 45.4081 * R + -47.383;
  if (R >= 1.1248 && R <= 1.1273) q = 46.6838 * R + -48.818;
  if (R >= 1.1273 && R <= 1.1298) q = 47.9778 * R + -50.2766;
  if (R >= 1.1298 && R <= 1.1323) q = 49.2898 * R + -51.759;
  if (R >= 1.1323 && R <= 1.1348) q = 50.6201 * R + -53.2652;
  if (R >= 1.1348 && R <= 1.1373) q = 51.9684 * R + -54.7953;
  if (R >= 1.1373 && R <= 1.1398) q = 53.335 * R + -56.3495;
  if (R >= 1.1398 && R <= 1.1423) q = 54.7196 * R + -57.9277;
  if (R >= 1.1423 && R <= 1.1448) q = 56.1225 * R + -59.5302;
  if (R >= 1.1448 && R <= 1.1473) q = 57.5435 * R + -61.1569;
  if (R >= 1.1473 && R <= 1.1498) q = 58.9826 * R + -62.808;
  if (R >= 1.1498 && R <= 1.1523) q = 60.4399 * R + -64.4836;
  if (R >= 1.1523 && R <= 1.1548) q = 61.9153 * R + -66.1838;
  if (R >= 1.1548 && R <= 1.1573) q = 63.4089 * R + -67.9086;
  if (R >= 1.1573 && R <= 1.1598) q = 64.9207 * R + -69.6581;
  if (R >= 1.1598 && R <= 1.1623) q = 66.4506 * R + -71.4325;
  if (R >= 1.1623 && R <= 1.1648) q = 67.9986 * R + -73.2318;
  if (R >= 1.1648 && R <= 1.1673) q = 69.5648 * R + -75.0561;
  if (R >= 1.1673 && R <= 1.1698) q = 71.1492 * R + -76.9055;
  if (R >= 1.1698 && R <= 1.1723) q = 72.7517 * R + -78.7801;
  if (R >= 1.1723 && R <= 1.1748) q = 74.3724 * R + -80.68;
  if (R >= 1.1748 && R <= 1.1773) q = 76.0112 * R + -82.6053;
  if (R >= 1.1773) q = 77.6681 * R + -84.556;
  return q;
}

double qprofile(double r, double R_axis) {
  double eta = 0.000;
  // return qprofile_TCV(r, R_axis);
  // return qprofile_lowshear(r, R_axis);
  // return qprofile_highshear(r, R_axis);
  // return (1-eta)*qprofile_lowshear(r, R_axis) + eta*qprofile_highshear(r, R_axis);
  // return qprofile_linear(r, R_axis);
  // return qprofile_quadratic(r, R_axis);
  // return qprofile_cubic(r, R_axis);
  return qprofile_pwlin8(r, R_axis);
  // return qprofile_pwlin48(r, R_axis);
}

void bc_shift_func_lo(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xc[0];
  double a_shift   = 0.5;               // Parameter in Shafranov shift.
  double Z_axis    = 0.1414361745;       // Magnetic axis height [m].
  double R_axis    = 0.8867856264;       // Magnetic axis major radius [m].
  double B_axis    = 1.4;                // Magnetic field at the magnetic axis [T].
  double R_LCFSmid = 1.0870056099999; // Major radius of the LCFS at the outboard midplane [m
  double x_inner   = 0.04;               // Radial extent inside LCFS    
  double x_outer   = 0.08;               // Radial extent outside LCFS
  double Rmid_min  = R_LCFSmid - x_inner;      // Minimum midplane major radius of simulation box [m].
  double Rmid_max  = R_LCFSmid + x_outer;      // Maximum midplane major radius of simulation box [m].
  double R0        = 0.5*(Rmid_min+Rmid_max);  // Major radius of the simulation box [m].
  double a_mid     = R_LCFSmid-R_axis;   // Minor radius at outboard midplane [m].
  // Redefine a_mid with Shafranov shift, to ensure LCFS radial location.
  a_mid = R_axis/a_shift - sqrt(R_axis*(R_axis - 2*a_shift*R_LCFSmid + 2*a_shift*R_axis))/a_shift;
  double r0        = R0-R_axis;           // Minor radius of the simulation box [m].
  double Lx        = Rmid_max-Rmid_min;
  double x_min     = 0.0;
  double x_max     = Lx;
  double q0        = qprofile(r_x(0.5*(x_min+x_max),a_mid,x_inner),R_axis);   
  double Lz = 2.*M_PI-1e-10;
  double r = r_x(x,a_mid,x_inner);
  fout[0] = -r0/q0*qprofile(r,R_axis)*Lz;
}

void bc_shift_func_up(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  bc_shift_func_lo(t, xc, fout, ctx);
  fout[0] *= -1;
}

void
test_bc_twistshift_3x_load_M0(const int *cells, enum gkyl_edge_loc edge,
  bool check_distf, bool use_gpu, bool write_f)
{
  /* 
  The file used in this test comes from a NT TCV simulation:
    /pscratch/sd/a/ah1032/gkyl_main/production/tcv_nt_48x32x16x12x6_src1.3/wk/gk_tcv_NT_iwl_3x2v-ion_M0_56.gkyl
  It contains the M0 ion moment and presents smooth sheared fluctuation at z=+pi but spurious kx sturctures at z=-pi.
  */
  const char *fname = "ctest_bc_twistshift_3x_M0_NT.gkyl";
  double eV = 1.6021766208e-19;
  double mp = 1.67262192e-27;
  double AMU = 2.01410177811;
  double qi = eV;
  double mass = mp*AMU;
  double T0 = 100*eV;
  double n0  = 2.0e19;   // [1/m^3]
  double vt = sqrt(T0/mass); // Thermal speeds.
  double B0 = 1.1214937537309428; // Magnetic field magnitude.
  int bc_dir = 2; // Direction in which to apply TS.
  double c_s = sqrt(T0/mass);
  double omega_ci = fabs(qi*B0/mass);
  double rho_s = c_s/omega_ci;
  double Lx = 0.12, Ly = 150*rho_s, Lz = 2*M_PI - 1e-10;
  double x_min = 0, x_max = Lx;
  double y_min = -0.5*Ly, y_max = 0.5*Ly;
  double z_min = 0, z_max = Lz;

  const int poly_order = 1;
  const double lower[] = {x_min, y_min, z_min};
  const double upper[] = {x_max, y_max, z_max};
  const int vdim = 0;
  const int ndim = sizeof(lower)/sizeof(lower[0]);
  const int cdim = ndim - vdim;

  double lower_conf[cdim], upper_conf[cdim];
  int cells_conf[cdim];
  for (int d=0; d<cdim; d++) {
    lower_conf[d] = lower[d];
    upper_conf[d] = upper[d];
    cells_conf[d] = cells[d];
  }

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid grid_conf;
  gkyl_rect_grid_init(&grid_conf, cdim, lower_conf, upper_conf, cells_conf);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  struct gkyl_basis basis_conf;
  gkyl_cart_modal_serendip(&basis_conf, cdim, poly_order);

  // Ranges.
  int ghost_conf[cdim];
  for (int d=0; d<cdim; d++) ghost_conf[d] = 1;
  struct gkyl_range local_conf, local_ext_conf; // local, local-ext position-space ranges
  gkyl_create_grid_ranges(&grid_conf, ghost_conf, &local_ext_conf, &local_conf);

  int ghost[ndim];
  for (int d=0; d<cdim; d++) ghost[d] = ghost_conf[d];
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct skin_ghost_ranges skin_ghost; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);
  struct skin_ghost_ranges skin_ghost_conf; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost_conf, &local_ext_conf, ghost_conf);

  // Pick skin and ghost ranges based on 'edge'.
  struct gkyl_range skin_rng, ghost_rng;
  struct gkyl_range skin_rng_conf, ghost_rng_conf;
  if (edge == GKYL_LOWER_EDGE) {
    skin_rng = skin_ghost.upper_skin[bc_dir];
    ghost_rng = skin_ghost.lower_ghost[bc_dir];
    skin_rng_conf = skin_ghost_conf.upper_skin[bc_dir];
    ghost_rng_conf = skin_ghost_conf.lower_ghost[bc_dir];
  }
  else {
    skin_rng = skin_ghost.lower_skin[bc_dir];
    ghost_rng = skin_ghost.upper_ghost[bc_dir];
    skin_rng_conf = skin_ghost_conf.lower_skin[bc_dir];
    ghost_rng_conf = skin_ghost_conf.upper_ghost[bc_dir];
  }

  struct test_bc_twistshift_ctx proj_ctx = {
    .lower = {lower[0], lower[1], lower[2]},
    .upper = {upper[0], upper[1], upper[2]},
    .cells = {cells[0], cells[1], cells[2]},
    .B0 = B0,
    .vt = vt,
    .mass = mass,
  };

  // declare donor and target arrays
  struct gkyl_array *donor = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *donor_ho = use_gpu? mkarr(false, basis.num_basis, local_ext.volume) : gkyl_array_acquire(donor);
  struct gkyl_array *target = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *target_ho = use_gpu? mkarr(false, basis.num_basis, local_ext.volume) : gkyl_array_acquire(target);

  // Initialize the distribution
//  gkyl_grid_sub_array_read(&grid, &local, donor, fname);
  gkyl_array_clear(donor, 0.0);
  gkyl_array_shiftc(donor, 1.0*pow(sqrt(2.0),cdim+vdim), 0);

  gkyl_array_copy(donor, donor_ho);
  struct gkyl_msgpack_data *mt = test_bc_twistshift_array_meta_new( (struct test_bc_twistshift_output_meta) {
      .poly_order = poly_order,
      .basis_type = basis.id
    }
  );

  // project the target function and weight
  // gkyl_proj_on_basis *projf = gkyl_proj_on_basis_new(
  //   &grid, &basis, poly_order + 1, 1, evalFunc, NULL);
  // gkyl_proj_on_basis_advance(projf, 0.0, &local, donor_ho);
  // gkyl_array_copy(donor, donor_ho);
  // gkyl_proj_on_basis_release(projf);

  // Write out the donor field.
  gkyl_grid_sub_array_write(&grid, &local, mt, donor_ho, "ctest_ts_qprofiles_donor.gkyl");

  // Copy the donor to the target and modify the target.
  gkyl_array_copy(target, donor);
  
  // Create a range only extended in bc_dir.
  struct gkyl_range update_rng;
  int lower_bcdir_ext[ndim], upper_bcdir_ext[ndim];
  for (int d=0; d<ndim; d++) {
    lower_bcdir_ext[d] = local.lower[d];
    upper_bcdir_ext[d] = local.upper[d];
  }
  lower_bcdir_ext[bc_dir] = local.lower[bc_dir] - ghost[bc_dir];
  upper_bcdir_ext[bc_dir] = local.upper[bc_dir] + ghost[bc_dir];
  gkyl_sub_range_init(&update_rng, &local_ext, lower_bcdir_ext, upper_bcdir_ext);

  // Create the twist-shift updater and shift the donor field.
  struct gkyl_bc_twistshift_inp tsinp = {
    .bc_dir = bc_dir,
    .shift_dir = 1, // y shift.
    .shear_dir = 0, // shift varies with x.
    .edge = edge,
    .cdim = cdim,
    .bcdir_ext_update_r = &update_rng,
    .num_ghost = ghost,
    .basis = &basis,
    .grid = &grid,
    .shift_func = bc_shift_func_lo,
    .shift_func_ctx = &proj_ctx,
    .use_gpu = use_gpu,
    .shift_poly_order = 2,
  };

  struct gkyl_bc_twistshift *tsup = gkyl_bc_twistshift_inew(&tsinp);

  // First apply periodicity in z.
  struct gkyl_array *buff_per = mkarr(use_gpu, basis.num_basis, skin_rng.volume);
  apply_periodic_bc(buff_per, target, bc_dir, skin_ghost);

  gkyl_bc_twistshift_advance(tsup, target, target);
  gkyl_array_copy(target_ho, target);

  // Write out the target in the extended range.
  double lower_ext[ndim], upper_ext[ndim];
  int cells_ext[ndim];
  for (int d=0; d<ndim; d++) {
    double dx = (upper[d]-lower[d])/cells[d];
    lower_ext[d] = lower[d]-dx*ghost[d];
    upper_ext[d] = upper[d]+dx*ghost[d];
    cells_ext[d] = cells[d]+2*ghost[d];
  }

  struct gkyl_rect_grid grid_ext;
  gkyl_rect_grid_init(&grid_ext, ndim, lower_ext, upper_ext, cells_ext);
  gkyl_grid_sub_array_write(&grid_ext, &local_ext, mt, target_ho, "ctest_ts_qprofiles_target.gkyl");

  // Build the error
  struct gkyl_array *err_ho = mkarr(false, basis.num_basis, local_ext.volume);
  gkyl_array_copy(err_ho, target_ho);
  gkyl_array_accumulate(err_ho, -1.0, donor_ho);

  // gkyl_grid_sub_array_write(&grid_ext, &local_ext, mt, err_ho, "error.gkyl");

  // Copy the ghost cell back to the skin cell, and apply the negated shift to
  // see if we recover the donor field approximately.
  gkyl_array_copy_range_to_range(target, target, &skin_rng, &ghost_rng);
  tsinp.shift_func = bc_shift_func_up;

  struct gkyl_bc_twistshift *tsup_m = gkyl_bc_twistshift_inew(&tsinp);
  gkyl_bc_twistshift_advance(tsup_m, target, target);
  gkyl_array_copy(target_ho, target);

  for (int d=0; d<ndim; d++) {
    double dx = (upper[d]-lower[d])/cells[d];
    lower_ext[d] = lower[d]-dx*ghost[d];
    upper_ext[d] = upper[d]+dx*ghost[d];
    cells_ext[d] = cells[d]+2*ghost[d];
  }
  gkyl_rect_grid_init(&grid_ext, ndim, lower_ext, upper_ext, cells_ext);
  gkyl_grid_sub_array_write(&grid_ext, &local_ext, mt, target_ho, "ctest_ts_qprofiles_target_unshifted.gkyl");

  gkyl_array_release(buff_per);
  test_bc_twistshift_array_meta_release(mt);
  gkyl_bc_twistshift_release(tsup);
  gkyl_bc_twistshift_release(tsup_m);
  gkyl_array_release(donor_ho);
  gkyl_array_release(donor);
  gkyl_array_release(target_ho);
  gkyl_array_release(target);
  gkyl_array_release(err_ho);
}

void mapc2p_vel_ion(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  double vpar_max_ion = 346299.94613717566;
  double mu_max_ion = 1.1428900033869495e-16;
  double cvpar = vc[0], cmu = vc[1];
  // Linear map up to vpar_max/2, then quadratic.
  if (fabs(cvpar) <= 0.5)
    vp[0] = vpar_max_ion*cvpar;
  else if (cvpar < -0.5)
    vp[0] = -vpar_max_ion*2.0*pow(cvpar,2);
  else
    vp[0] =  vpar_max_ion*2.0*pow(cvpar,2);
  // Quadratic map in mu.
  vp[1] = mu_max_ion*pow(cmu,2);
}

void
test_bc_twistshift_3x2v_load_ion(const int *cells, enum gkyl_edge_loc edge,
  bool check_distf, bool use_gpu, bool write_f)
{
  /* 
  The file used in this test comes from a NT TCV simulation:
    /pscratch/sd/a/ah1032/gkyl_main/production/tcv_nt_48x32x16x12x6_src1.3/wk/gk_tcv_NT_iwl_3x2v-ion_M0_56.gkyl
  It contains the M0 ion moment and presents smooth sheared fluctuation at z=+pi but spurious kx sturctures at z=-pi.
  */
  const char *fname = "ctest_bc_twistshift_3x2v_ion_NT.gkyl";
  double eV = 1.6021766208e-19;
  double mp = 1.67262192e-27;
  double AMU = 2.01410177811;
  double qi = eV;
  double mass = mp*AMU;
  double T0 = 100*eV;
  double n0  = 2.0e19;   // [1/m^3]
  double vt = sqrt(T0/mass); // Thermal speeds.
  double B0 = 1.1214937537309428; // Magnetic field magnitude.
  int bc_dir = 2; // Direction in which to apply TS.
  double c_s = sqrt(T0/mass);
  double omega_ci = fabs(qi*B0/mass);
  double rho_s = c_s/omega_ci;
  double Lx = 0.12, Ly = 150*rho_s, Lz = 2*M_PI - 1e-10;
  double x_min = 0, x_max = Lx;
  double y_min = -0.5*Ly, y_max = 0.5*Ly;
  double z_min = 0, z_max = Lz;
  double vtmax = 5*vt;
  double mumax = mp*pow(4*vt,2)/(2*B0);

  const int poly_order = 1;
  const double lower[] = {x_min, y_min, z_min, -5*vtmax, 0.};
  const double upper[] = {x_max, y_max, z_max,  5*vtmax, mumax};
  const int vdim = 2;
  const int ndim = sizeof(lower)/sizeof(lower[0]);
  const int cdim = ndim - vdim;

  double lower_conf[cdim], upper_conf[cdim];
  int cells_conf[cdim];
  for (int d=0; d<cdim; d++) {
    lower_conf[d] = lower[d];
    upper_conf[d] = upper[d];
    cells_conf[d] = cells[d];
  }
  double lower_vel[vdim], upper_vel[vdim];
  int cells_vel[vdim];
  for (int d=0; d<vdim; d++) {
    lower_vel[d] = lower[cdim+d];
    upper_vel[d] = upper[cdim+d];
    cells_vel[d] = cells[cdim+d];
  }
  

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid grid_conf;
  gkyl_rect_grid_init(&grid_conf, cdim, lower_conf, upper_conf, cells_conf);
  struct gkyl_rect_grid grid_vel;
  gkyl_rect_grid_init(&grid_vel, vdim, lower_vel, upper_vel, cells_vel);

	// Basis functions.
	struct gkyl_basis basis;
	if (poly_order == 1) 
			gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
	else
			gkyl_cart_modal_serendip(&basis, ndim, poly_order);
	struct gkyl_basis basis_conf;
	gkyl_cart_modal_serendip(&basis_conf, cdim, poly_order);

  // Ranges.
  int ghost_conf[cdim];
  for (int d=0; d<cdim; d++) ghost_conf[d] = 1;
  struct gkyl_range local_conf, local_ext_conf; // local, local-ext position-space ranges
  gkyl_create_grid_ranges(&grid_conf, ghost_conf, &local_ext_conf, &local_conf);

  int ghost_vel[vdim];
  for (int d=0; d<vdim; d++) ghost_vel[d] = 0;
  struct gkyl_range local_vel, local_ext_vel; // local, local-ext position-space ranges
  gkyl_create_grid_ranges(&grid_vel, ghost_vel, &local_ext_vel, &local_vel);

  int ghost[ndim];
  for (int d=0; d<cdim; d++) ghost[d] = ghost_conf[d];
  for (int d=cdim; d<ndim; d++) ghost[d] = ghost_vel[d-cdim];
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct skin_ghost_ranges skin_ghost; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);
  struct skin_ghost_ranges skin_ghost_conf; // skin/ghost.
  skin_ghost_ranges_init(&skin_ghost_conf, &local_ext_conf, ghost_conf);

  // Pick skin and ghost ranges based on 'edge'.
  struct gkyl_range skin_rng, ghost_rng;
  struct gkyl_range skin_rng_conf, ghost_rng_conf;
  if (edge == GKYL_LOWER_EDGE) {
    skin_rng = skin_ghost.upper_skin[bc_dir];
    ghost_rng = skin_ghost.lower_ghost[bc_dir];
    skin_rng_conf = skin_ghost_conf.upper_skin[bc_dir];
    ghost_rng_conf = skin_ghost_conf.lower_ghost[bc_dir];
  }
  else {
    skin_rng = skin_ghost.lower_skin[bc_dir];
    ghost_rng = skin_ghost.upper_ghost[bc_dir];
    skin_rng_conf = skin_ghost_conf.lower_skin[bc_dir];
    ghost_rng_conf = skin_ghost_conf.upper_ghost[bc_dir];
  }

  struct test_bc_twistshift_ctx proj_ctx = {
    .lower = {lower[0], lower[1], lower[2], lower[3], lower[4]},
    .upper = {upper[0], upper[1], upper[2], upper[3], upper[4]},
    .cells = {cells[0], cells[1], cells[2], cells[3], cells[4]},
    .B0 = B0,
    .vt = vt,
    .mass = mass,
  };

  // Initialize the distribution
  struct gkyl_array *distf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu? mkarr(false, basis.num_basis, local_ext.volume) : gkyl_array_acquire(distf);
//  gkyl_grid_sub_array_read(&grid, &local, distf, fname);
  gkyl_array_clear(distf, 0.0);
  gkyl_array_shiftc(distf, 1.0*pow(sqrt(2.0),cdim+vdim), 0);

  gkyl_array_copy(distf, distf_ho);
  struct gkyl_msgpack_data *mt = test_bc_twistshift_array_meta_new( (struct test_bc_twistshift_output_meta) {
      .poly_order = poly_order,
      .basis_type = basis.id
    }
  );
  // Write out the donor field.
  gkyl_grid_sub_array_write(&grid, &local, mt, distf_ho, "ctest_ts_qprofiles_3x2v_ion_NT_do.gkyl");

  // Create a range only extended in bc_dir.
  struct gkyl_range update_rng;
  int lower_bcdir_ext[ndim], upper_bcdir_ext[ndim];
  for (int d=0; d<ndim; d++) {
    lower_bcdir_ext[d] = local.lower[d];
    upper_bcdir_ext[d] = local.upper[d];
  }
  lower_bcdir_ext[bc_dir] = local.lower[bc_dir] - ghost[bc_dir];
  upper_bcdir_ext[bc_dir] = local.upper[bc_dir] + ghost[bc_dir];
  gkyl_sub_range_init(&update_rng, &local_ext, lower_bcdir_ext, upper_bcdir_ext);

  // Create the twist-shift updater and shift the donor field.
  struct gkyl_bc_twistshift_inp tsinp = {
    .bc_dir = bc_dir,
    .shift_dir = 1, // y shift.
    .shear_dir = 0, // shift varies with x.
    .edge = edge,
    .cdim = cdim,
    .bcdir_ext_update_r = &update_rng,
    .num_ghost = ghost,
    .basis = &basis,
    .grid = &grid,
    .shift_func = bc_shift_func_lo,
    .shift_func_ctx = &proj_ctx,
    .use_gpu = use_gpu,
  };

  struct gkyl_bc_twistshift *tsup = gkyl_bc_twistshift_inew(&tsinp);

  // First apply periodicity in z.
  struct gkyl_array *buff_per = mkarr(use_gpu, basis.num_basis, skin_rng.volume);
  apply_periodic_bc(buff_per, distf, bc_dir, skin_ghost);

  gkyl_bc_twistshift_advance(tsup, distf, distf);
  gkyl_array_copy(distf_ho, distf);

  // Write out the target in the extended range.
  double lower_ext[ndim], upper_ext[ndim];
  int cells_ext[ndim];
  for (int d=0; d<ndim; d++) {
    double dx = (upper[d]-lower[d])/cells[d];
    lower_ext[d] = lower[d]-dx*ghost[d];
    upper_ext[d] = upper[d]+dx*ghost[d];
    cells_ext[d] = cells[d]+2*ghost[d];
  }
  struct gkyl_rect_grid grid_ext;
  gkyl_rect_grid_init(&grid_ext, ndim, lower_ext, upper_ext, cells_ext);
  gkyl_grid_sub_array_write(&grid_ext, &local_ext, mt, distf_ho, "ctest_ts_qprofiles_3x2v_ion_NT_tar.gkyl");

  // Copy the ghost cell back to the skin cell, and apply the negated shift to
  // see if we recover the donor field approximately.
  gkyl_array_copy_range_to_range(distf, distf, &skin_rng, &ghost_rng);
  tsinp.shift_func = bc_shift_func_up;

  struct gkyl_bc_twistshift *tsup_m = gkyl_bc_twistshift_inew(&tsinp);
  gkyl_bc_twistshift_advance(tsup_m, distf, distf);
  gkyl_array_copy(distf_ho, distf);

  for (int d=0; d<ndim; d++) {
    double dx = (upper[d]-lower[d])/cells[d];
    lower_ext[d] = lower[d]-dx*ghost[d];
    upper_ext[d] = upper[d]+dx*ghost[d];
    cells_ext[d] = cells[d]+2*ghost[d];
  }
  gkyl_rect_grid_init(&grid_ext, ndim, lower_ext, upper_ext, cells_ext);
  gkyl_grid_sub_array_write(&grid_ext, &local_ext, mt, distf_ho, "ctest_ts_qprofiles_3x_M0_NT_tar_shifted.gkyl");

  gkyl_array_release(buff_per);
  test_bc_twistshift_array_meta_release(mt);
  gkyl_bc_twistshift_release(tsup);
  gkyl_bc_twistshift_release(tsup_m);
  gkyl_array_release(distf_ho);
  gkyl_array_release(distf);
}


void
test_3x(bool use_gpu)
{
  const int cells0[] = {48, 32, 16};

  enum gkyl_edge_loc edgelo = GKYL_LOWER_EDGE; // Lower edge.
  test_bc_twistshift_3x_load_M0(cells0, edgelo, true, use_gpu, true);
}

void test_ts_qprofiles_3x_ho(){ test_3x(false); }


void
test_3x2v(bool use_gpu)
{
  const int cells0[] = {48, 32, 16, 12, 6};

  enum gkyl_edge_loc edgelo = GKYL_LOWER_EDGE; // Lower edge.
  test_bc_twistshift_3x2v_load_ion(cells0, edgelo, true, use_gpu, false);
}

void test_ts_qprofiles_3x2v_ho(){ test_3x2v(false); }


#ifdef GKYL_HAVE_CUDA
void test_dev(){ test_3x(true); }
#endif

TEST_LIST = {
  { "test_ts_qprofiles_3x_ho", test_ts_qprofiles_3x_ho },
	{ "test_ts_qprofiles_3x2v_ho", test_ts_qprofiles_3x2v_ho },
#ifdef GKYL_HAVE_CUDA
  { "test_ts_qprofiles_3x_dev", test_ts_qprofiles_3x_dev },
#endif
  { NULL, NULL },
};

/*
pgkyl ctest_bc_twistshift_3x_M0_NT.gkyl interp sel --z0 :32 --z2 0 pl -x 'x' -y 'y' --clabel '$M_0(z=z_{max}-\Delta z/4)$' --zmin 0 --zmax 5e19 &

pgkyl ctest_bc_twistshift_3x_M0_NT_tar.gkyl interp sel --z0 2:34 --z1 2:66 --z2 1 pl -x 'x' -y 'y' --clabel '$M_0(z=z_{min}-\Delta z/4)$' &


pgkyl ctest_bc_twistshift_3x2v_ion_NT.gkyl interp sel --z0 :32 --z2 0 --z3 0 --z4 0 pl -x 'x' -y 'y' --clabel '$f_i(z=-\pi, v_\parallel=0, \mu=0)$' &

pgkyl ctest_bc_twistshift_3x2v_ion_NT.gkyl interp sel --z0 :32 --z2 -1 --z3 0 --z4 0 pl -x 'x' -y 'y' --clabel '$f_i(z= \pi, v_\parallel=0, \mu=0)$' &

pgkyl ctest_bc_twistshift_3x2v_ion_NT_do.gkyl interp sel --z0 :32 --z2 0 --z3 0 --z4 0 pl -x 'x' -y 'y' --clabel '$f_i(z=-\pi, v_\parallel=0, \mu=0)$' &

pgkyl ctest_bc_twistshift_3x2v_ion_NT_do.gkyl interp sel --z0 :32 --z2 -1 --z3 0 --z4 0 pl -x 'x' -y 'y' --clabel '$f_i(z= \pi, v_\parallel=0, \mu=0)$' &

pgkyl ctest_bc_twistshift_3x2v_ion_NT_tar.gkyl interp sel --z0 2:34 --z1 2:66 --z2 1 --z3 0 --z4 0 pl -x 'x' -y 'y' --clabel '$f_i(z=-\pi, v_\parallel=0, \mu=0)$' &

*/
