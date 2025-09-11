#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_gk_anomalous_diffusion.h>
#include <gkyl_dg_updater_gk_anomalous_diffusion.h>
#include <gkyl_dg_updater_gk_anomalous_diffusion_priv.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

struct gkyl_dg_eqn*
gkyl_dg_updater_gk_anomalous_diffusion_acquire_eqn(const struct gkyl_dg_updater_gk_anomalous_diffusion *up)
{
  return gkyl_dg_eqn_acquire(up->dgeqn);
}

struct gkyl_dg_updater_gk_anomalous_diffusion*
gkyl_dg_updater_gk_anomalous_diffusion_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, const struct gkyl_basis *cbasis, const struct gkyl_range *conf_range,
  const bool *is_zero_flux_bc, const bool *is_absorb_bc, double skip_cell_threshold, const struct gkyl_array *nu,
  const struct gkyl_array *jacobgeo_inv, bool use_gpu)
{
  struct gkyl_dg_updater_gk_anomalous_diffusion *up = gkyl_malloc(sizeof(struct gkyl_dg_updater_gk_anomalous_diffusion));

  int pdim = basis->ndim;
  int cdim = cbasis->ndim;
  up->use_gpu = use_gpu;

  // Choice of boundary_surf and boundary_diag kernels:
  //   boundary_surf: zero_flux or local
  //   boundary_diag: local or recovery
  // MF 2025/09/10: as of now these options are meant for (here
  // N/A means not applicable):
  //            bound_surf  bound_diag  hyper_dg-zero_flux
  // PERIODIC:  N/A         N/A         no
  // COPY:      N/A         recovery    no
  // SKIP:      N/A         recovery    no
  // ABSORB:    local       local       yes
  // FUNC:      N/A         recovery    no
  // ZERO_FLUX: zero_flux   N/A         yes

  up->dgeqn = gkyl_gk_anomalous_diffusion_new(basis, cbasis,
    conf_range, is_zero_flux_bc, is_absorb_bc, skip_cell_threshold, up->use_gpu);

  gkyl_gk_anomalous_diffusion_set_auxfields(up->dgeqn, (struct gkyl_gk_anomalous_diffusion_auxfields) {
    .nu = nu, .jacobgeo_inv = jacobgeo_inv });

  int num_up_dirs = 1;
  int up_dirs[GKYL_MAX_DIM];
  up_dirs[0] = 0;

  // Determine if hyper_dg should apply boundary_surf kernels.
  int use_boundary_surf[2*GKYL_MAX_DIM] = {0};
  use_boundary_surf[0]      = is_zero_flux_bc[0]      || is_absorb_bc[0]?      1 : 0;
  use_boundary_surf[0+pdim] = is_zero_flux_bc[0+pdim] || is_absorb_bc[0+pdim]? 1 : 0;

  up->hyperdg = gkyl_hyper_dg_new(grid, basis, up->dgeqn, num_up_dirs, up_dirs, use_boundary_surf, 1, up->use_gpu);

  up->diffusion_tm = 0.0;

  return up;
}

void
gkyl_dg_updater_gk_anomalous_diffusion_advance(struct gkyl_dg_updater_gk_anomalous_diffusion *up,
  const struct gkyl_range *update_rng, const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs)
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_hyper_dg_advance(up->hyperdg, update_rng, fIn, cflrate, rhs);
  up->diffusion_tm += gkyl_time_diff_now_sec(wst);
}

struct gkyl_dg_updater_gk_anomalous_diffusion_tm
gkyl_dg_updater_gk_anomalous_diffusion_get_tm(const struct gkyl_dg_updater_gk_anomalous_diffusion *up)
{
  return (struct gkyl_dg_updater_gk_anomalous_diffusion_tm) {
    .diffusion_tm = up->diffusion_tm,
  };
}

void
gkyl_dg_updater_gk_anomalous_diffusion_release(struct gkyl_dg_updater_gk_anomalous_diffusion *up)
{
  gkyl_dg_eqn_release(up->dgeqn);
  gkyl_hyper_dg_release(up->hyperdg);
  gkyl_free(up);
}
