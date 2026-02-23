// Test creation and deallocation of updater that applies the
// conducting sheath BC for gyrokinetics.
//
#include <acutest.h>

#include <gkyl_array_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>
#include <gkyl_bc_sheath_gyrokinetic.h>
#include <gkyl_velocity_map.h>
#include <gkyl_proj_on_basis.h>
#include <mpack.h>

// Meta-data for IO
struct test_bc_sheath_output_meta {
  int poly_order; // polynomial order
  const char *basis_type; // name of basis functions
};

// returned gkyl_array_meta must be freed using gyrokinetic_array_meta_release
static struct gkyl_msgpack_data*
test_bc_sheath_array_meta_new(struct test_bc_sheath_output_meta meta)
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
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM]; // Grid extents.
  int cells[GKYL_MAX_DIM]; // Number of cells.
  double upar; // Mean parallel drift speed.
  double B0; // Reference magnetic field.
  double vt; // Reference thermal speed.
  double mass; // Species mass.
};

void
eval_func_1x2v(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double vpar = xn[1], mu = xn[2];

  struct test_sheath_ctx *pars = ctx;

  double vt = pars->vt;
  double B0 = pars->B0;
  double m = pars->mass;
  double upar = pars->upar;

  fout[0] = exp( -(pow(vpar-upar,2) + 2.0*mu*B0/m) /(2.0*pow(vt,2)) );
}

void
test_bc_sheath_gyrokinetic_1x2v(const int *cells, enum gkyl_edge_loc edge, int dir, bool write_fields, bool use_gpu)
{

  double mass = 1.; // Species mass.
  double charge = -1.; // Species charge (+ or - for electrons/ions).
  double vt = 5.0; // Reference thermal speed (for grid extents).
  double B0 = 1.0; // Magnetic field magnitude.
  double upar = 0.0; // Parallel flow speed in distribution function.
  double vt = 1.5*vt; // Thermal speed in distribution function.
  double phi_wall = 0.0; // Potential at wall.
  double phi_mpe = 10.0; // Potential at the magnetic presheath entrance.
  
  int poly_order = 1;
  double lower[] = {-2.0, -5.0*vt, 0.},
  double upper[] = { 2.0,  5.0*vt, mass*(pow(5.0*vt,2))/(2.0*B0)};
  int vdim = 2;
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int cdim = ndim - vdim;

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

  int ghost_conf[] = { 1, 1, 1 };
  struct gkyl_range local_conf, local_conf_ext; // local, local-ext position-space ranges
  gkyl_create_grid_ranges(&grid_conf, ghost, &local_conf_ext, &local_conf);

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
    double dx = (upper[d]-lower[d])/cells[d];
    lower_ext[d] = lower[d]-dx*ghost[d];
    upper_ext[d] = upper[d]+dx*ghost[d];
    cells_ext[d] = cells[d]+2*ghost[d];
  }
  struct gkyl_rect_grid grid_ext;
  gkyl_rect_grid_init(&grid_ext, ndim, lower_ext, upper_ext, cells_ext);

  // Create the skin/ghost ranges.
  struct gkyl_range skin_r, ghost_r;
  gkyl_skin_ghost_ranges(&skin_r, &ghost_r, dir, edge, &local_ext, ghost);

  // Initialize the distribution
  struct gkyl_array *distf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu? mkarr(false, distf->ncomp, distf->size) : gkyl_array_acquire(distf);
  struct test_sheath_ctx proj_ctx = {
    .B0 = B0,
    .vt = vt,
    .mass = mass,
    .upar = upar,
  };
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = &grid,
      .basis = &basis,
      .num_ret_vals = 1,
      .eval = eval_func_1x2v,
      .ctx = &proj_ctx,
    }
  );
  gkyl_proj_on_basis_advance(projDistf, 0.0, &local, distf_ho);
  gkyl_array_copy(distf, distf_ho);

  // Initialize the electrostatic potential at MPE and WALL
  double dgnormc = pow(sqrt(2.), cdim);
  struct gkyl_array *phi = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
  struct gkyl_array *phi_ho = use_gpu? mkarr(false, phi->ncomp, phi->size) : gkyl_array_acquire(phi);
  gkyl_array_shiftc(phi_ho, phi_mpe * dgnormc, 0 * basis_conf.num_basis);
  gkyl_array_copy(phi, phi_ho);

  struct gkyl_array *phiw = mkarr(use_gpu, basis_conf.num_basis, local_conf_ext.volume);
  struct gkyl_array *phiw_ho = use_gpu? mkarr(false, phiw_ho->ncomp, phiw_ho->size) : gkyl_array_acquire(phiw);
  gkyl_array_shiftc(phiw_ho, phi_wall * dgnormc, 0 * basis_conf.num_basis);
  gkyl_array_copy(phiw, phiw_ho);

  // Write out input fields if requested.
  struct gkyl_msgpack_data *mt = test_bc_sheath_array_meta_new( (struct test_bc_sheath_output_meta) {
      .poly_order = poly_order,
      .basis_type = basis.id
    }
  );
  if (write_fields) {
    gkyl_grid_sub_array_write(&grid_ext, &local_ext, mt, distf_ho, "bc_sheath_1x2v_distf_skin.gkyl");
    gkyl_grid_sub_array_write(&grid_conf, &local_conf, mt, phi_ho, "bc_sheath_1x2v_phi_mpe.gkyl");
    gkyl_grid_sub_array_write(&grid_conf, &local_conf, mt, phiw_ho, "bc_sheath_1x2v_phi_wall.gkyl");
  }

  // Create the BC updater.
  struct gkyl_bc_sheath_gyrokinetic *bcsheath = gkyl_bc_sheath_gyrokinetic_new(dir, edge,
    &basis, &skin_r, &ghost_r, gvm, cdim, 2.*charge/mass, use_gpu);

  // Advance the BC updater.
  gkyl_bc_sheath_gyrokinetic_advance(bcsheath, phi, phiw, distf, &local_conf);

  // Write out the distribution function after applying BC if requested.
  if (write_fields) {
    gkyl_grid_sub_array_write(&grid_ext, &local_ext, mt, distf_ho, "bc_sheath_1x2v_distf_ghost.gkyl");
  }

  // Clean up.
  gkyl_msgpack_data_release(mt);
  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
  gkyl_array_release(distf_ho);
  gkyl_array_release(phi);
  gkyl_array_release(phi_ho);
  gkyl_array_release(phiw);
  gkyl_array_release(phiw_ho);

  gkyl_velocity_map_release(gvm);
  gkyl_bc_sheath_gyrokinetic_release(bcsheath);
}

void test_bc_gksheath(){
  test_bc_sheath_gyrokinetic_1x2v((int[]){4, 16, 12}, GKYL_LOWER_EDGE, 0, true, false);
}

#ifdef GKYL_HAVE_CUDA
void test_bc_gksheath_cu(){
  test_bc_sheath_gyrokinetic_1x2v((int[]){4, 16, 12}, GKYL_LOWER_EDGE, 0, true, true);
}
#endif

TEST_LIST = {
  { "test_bc_gksheath", test_bc_gksheath },
#ifdef GKYL_HAVE_CUDA
  { "test_bc_gksheath_cu", test_bc_gksheath_cu },
#endif
  { NULL, NULL },
};
