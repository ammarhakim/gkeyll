#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_numerical_diffusion_gyrokinetic.h>
#include <gkyl_dg_updater_gk_numerical_diffusion.h>
#include <gkyl_dg_updater_gk_numerical_diffusion_priv.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

struct gkyl_dg_eqn*
gkyl_dg_updater_gk_numerical_diffusion_acquire_eqn(const struct gkyl_dg_updater_gk_numerical_diffusion *up)
{
  return gkyl_dg_eqn_acquire(up->dgeqn);
}

struct gkyl_dg_updater_gk_numerical_diffusion*
gkyl_dg_updater_gk_numerical_diffusion_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, const struct gkyl_basis *cbasis, const struct gkyl_range *conf_range,
  const enum gkyl_gyrokinetic_bc_type *bc_lower, const enum gkyl_gyrokinetic_bc_type *bc_upper,
  const bool *has_diff_in_dir, int diff_order,
  const struct gkyl_array *nu, const struct gkyl_array *jacobgeo_inv, bool use_gpu)
{
  struct gkyl_dg_updater_gk_numerical_diffusion *up = gkyl_malloc(sizeof(struct gkyl_dg_updater_gk_numerical_diffusion));

  int pdim = basis->ndim;
  int cdim = cbasis->ndim;
  up->use_gpu = use_gpu;

  // Choice of boundary_surf and boundary_diag kernels:
  //   boundary_surf: zero_flux or local
  //   boundary_diag: local or recovery
  // MF 2025/09/10: as of now these options are meant for (here
  // N/A means not applicable):
  //            bound_surf  bound_diag  hyper_dg-zero_flux
  // SKIP:       N/A         recovery    no
  // ABSORB:     N/A         N/A         no
  // PERIODIC:   N/A         N/A         no
  // FIXED_FUNC: N/A         N/A         no
  // ZERO_FLUX:  zero_flux   N/A         yes
  // ELSE:       local       local       yes

  up->dgeqn = gkyl_dg_numerical_diffusion_gyrokinetic_new(basis, cbasis, conf_range,
    has_diff_in_dir, diff_order, bc_lower, bc_upper, up->use_gpu);

  gkyl_dg_numerical_diffusion_gyrokinetic_set_auxfields(up->dgeqn, (struct gkyl_dg_numerical_diffusion_gyrokinetic_auxfields) {
    .nu = nu, .jacobgeo_inv = jacobgeo_inv });

  int num_up_dirs = 0;
  for (int d=0; d<cdim; d++) num_up_dirs += has_diff_in_dir[d]? 1 : 0;

  int up_dirs[GKYL_MAX_DIM], zero_flux_flags[2*GKYL_MAX_DIM];
  int linc = 0;
  for (int d=0; d<cdim; ++d) {
    if (has_diff_in_dir[d]) up_dirs[linc] = d;
    linc += 1;

    zero_flux_flags[d] = is_zero_flux_bc[d]? 1 : 0;
    zero_flux_flags[d+pdim] = is_zero_flux_bc[d+pdim]? 1 : 0;
  }

  // Determine if hyper_dg should apply boundary_surf kernels.
  int use_boundary_surf[2*GKYL_MAX_DIM] = {0};
  for (int d=0; d<cdim; ++d) {
    if ( !((bc_lower[d] == GKYL_BC_GK_SKIP) ||
           (bc_lower[d] == GKYL_BC_GK_SPECIES_ABSORB) ||
           (bc_lower[d] == GKYL_BC_GK_SPECIES_PERIODIC) ||
           (bc_lower[d] == GKYL_BC_GK_SPECIES_FIXED_FUNC)) )
      use_boundary_surf[d]      = 1;
    if ( !((bc_upper[d] == GKYL_BC_GK_SKIP) ||
           (bc_upper[d] == GKYL_BC_GK_SPECIES_ABSORB) ||
           (bc_upper[d] == GKYL_BC_GK_SPECIES_PERIODIC) ||
           (bc_upper[d] == GKYL_BC_GK_SPECIES_FIXED_FUNC)) )
      use_boundary_surf[d+pdim] = 1;
  }

  up->hyperdg = gkyl_hyper_dg_new(grid, basis, up->dgeqn, num_up_dirs, up_dirs, use_boundary_surf, 1, up->use_gpu);

  up->diffusion_tm = 0.0;

  return up;
}

void
gkyl_dg_updater_gk_numerical_diffusion_advance(struct gkyl_dg_updater_gk_numerical_diffusion *up,
  const struct gkyl_range *update_rng, const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs)
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_hyper_dg_advance(up->hyperdg, update_rng, fIn, cflrate, rhs);
  up->diffusion_tm += gkyl_time_diff_now_sec(wst);
}

struct gkyl_dg_updater_gk_numerical_diffusion_tm
gkyl_dg_updater_gk_numerical_diffusion_get_tm(const struct gkyl_dg_updater_gk_numerical_diffusion *up)
{
  return (struct gkyl_dg_updater_gk_numerical_diffusion_tm) {
    .diffusion_tm = up->diffusion_tm,
  };
}

void
gkyl_dg_updater_gk_numerical_diffusion_release(struct gkyl_dg_updater_gk_numerical_diffusion *up)
{
  gkyl_dg_eqn_release(up->dgeqn);
  gkyl_hyper_dg_release(up->hyperdg);
  gkyl_free(up);
}
