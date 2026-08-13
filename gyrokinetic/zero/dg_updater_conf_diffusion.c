#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_conf_diffusion.h>
#include <gkyl_dg_updater_conf_diffusion.h>
#include <gkyl_dg_updater_conf_diffusion_priv.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

struct gkyl_dg_updater_conf_diffusion*
gkyl_dg_updater_conf_diffusion_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, const struct gkyl_range *conf_range,
  const struct gkyl_array *diffusion_tensor,
  const struct gkyl_array *jacobgeo_inv,
  const int zero_flux_flags[2*GKYL_MAX_CDIM], bool use_gpu)
{
  struct gkyl_dg_updater_conf_diffusion *up = gkyl_malloc(sizeof(*up));
  up->dgeqn = gkyl_conf_diffusion_new(basis, conf_range, use_gpu);
  gkyl_conf_diffusion_set_auxfields(up->dgeqn,
    (struct gkyl_conf_diffusion_auxfields) {
      .diffusion_tensor = diffusion_tensor,
      .jacobgeo_inv = jacobgeo_inv,
    });

  int up_dirs[GKYL_MAX_DIM] = { 0 };
  int boundary_flags[2*GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<basis->ndim; ++d)
    up_dirs[d] = d;
  for (int d=0; d<2*basis->ndim; ++d)
    boundary_flags[d] = zero_flux_flags[d];
  up->hyperdg = gkyl_hyper_dg_new(grid, basis, up->dgeqn, basis->ndim, up_dirs,
    boundary_flags, 1, use_gpu);
  up->diffusion_tm = 0.0;
  return up;
}

void
gkyl_dg_updater_conf_diffusion_advance(
  struct gkyl_dg_updater_conf_diffusion *up,
  const struct gkyl_range *update_rng,
  const struct gkyl_array *GKYL_RESTRICT fIn,
  struct gkyl_array *GKYL_RESTRICT cflrate,
  struct gkyl_array *GKYL_RESTRICT rhs)
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_hyper_dg_advance(up->hyperdg, update_rng, fIn, cflrate, rhs);
  up->diffusion_tm += gkyl_time_diff_now_sec(wst);
}

struct gkyl_dg_updater_conf_diffusion_tm
gkyl_dg_updater_conf_diffusion_get_tm(
  const struct gkyl_dg_updater_conf_diffusion *up)
{
  return (struct gkyl_dg_updater_conf_diffusion_tm) {
    .diffusion_tm = up->diffusion_tm,
  };
}

void
gkyl_dg_updater_conf_diffusion_release(
  struct gkyl_dg_updater_conf_diffusion *up)
{
  gkyl_dg_eqn_release(up->dgeqn);
  gkyl_hyper_dg_release(up->hyperdg);
  gkyl_free(up);
}
