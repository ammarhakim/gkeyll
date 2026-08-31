/* -*- c++ -*- */
extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_cart_modal_serendip_priv.h>
#include <gkyl_conf_diffusion_priv.h>
}

__global__ static void
conf_diffusion_set_auxfields_cu_kernel(const struct gkyl_dg_eqn *eqn,
  const struct gkyl_array *diffusion_tensor,
  const struct gkyl_array *jacobgeo_inv)
{
  struct conf_diffusion *diffusion =
    container_of(eqn, struct conf_diffusion, eqn);
  diffusion->auxfields.diffusion_tensor = diffusion_tensor;
  diffusion->auxfields.jacobgeo_inv = jacobgeo_inv;
}

void
gkyl_conf_diffusion_set_auxfields_cu(const struct gkyl_dg_eqn *eqn,
  struct gkyl_conf_diffusion_auxfields auxin)
{
  conf_diffusion_set_auxfields_cu_kernel<<<1,1>>>(eqn,
    auxin.diffusion_tensor->on_dev, auxin.jacobgeo_inv->on_dev);
}

__global__ static void
conf_diffusion_set_cu_dev_ptrs(struct conf_diffusion *diffusion, int cdim)
{
  diffusion->basis.flip_odd_sign = fos_list[cdim].fs[1];
  diffusion->basis.flip_even_sign = fes_list[cdim].fs[1];
  diffusion->eqn.num_equations = 1;
  diffusion->eqn.vol_term = conf_diffusion_vol;
  diffusion->eqn.surf_term = conf_diffusion_surf;
  diffusion->eqn.boundary_surf_term = conf_diffusion_boundary_surf;
  diffusion->eqn.boundary_diag_term = conf_diffusion_boundary_surf;
  diffusion->eqn.gen_surf_term = 0;
  diffusion->eqn.gen_boundary_surf_term = 0;

  if (cdim == 1) {
    diffusion->vol = conf_diffusion_vol_1x_ser_p1;
    diffusion->surf[0] = conf_diffusion_surfx_1x_ser_p1;
    diffusion->boundary_surf[0] = conf_diffusion_boundary_surfx_1x_ser_p1;
  }
  else if (cdim == 2) {
    diffusion->vol = conf_diffusion_vol_2x_ser_p1;
    diffusion->surf[0] = conf_diffusion_surfx_2x_ser_p1;
    diffusion->surf[1] = conf_diffusion_surfy_2x_ser_p1;
    diffusion->boundary_surf[0] = conf_diffusion_boundary_surfx_2x_ser_p1;
    diffusion->boundary_surf[1] = conf_diffusion_boundary_surfy_2x_ser_p1;
  }
  else {
    diffusion->vol = conf_diffusion_vol_3x_ser_p1;
    diffusion->surf[0] = conf_diffusion_surfx_3x_ser_p1;
    diffusion->surf[1] = conf_diffusion_surfy_3x_ser_p1;
    diffusion->surf[2] = conf_diffusion_surfz_3x_ser_p1;
    diffusion->boundary_surf[0] = conf_diffusion_boundary_surfx_3x_ser_p1;
    diffusion->boundary_surf[1] = conf_diffusion_boundary_surfy_3x_ser_p1;
    diffusion->boundary_surf[2] = conf_diffusion_boundary_surfz_3x_ser_p1;
  }
  for (int d=cdim; d<GKYL_MAX_CDIM; ++d) {
    diffusion->surf[d] = 0;
    diffusion->boundary_surf[d] = 0;
  }
  diffusion->auxfields.diffusion_tensor = 0;
  diffusion->auxfields.jacobgeo_inv = 0;
}

struct gkyl_dg_eqn*
gkyl_conf_diffusion_cu_dev_new(const struct gkyl_basis *basis,
  const struct gkyl_range *conf_range)
{
  struct conf_diffusion *diffusion =
    (struct conf_diffusion*) gkyl_malloc(sizeof(*diffusion));
  diffusion->basis = *basis;
  diffusion->conf_range = *conf_range;
  diffusion->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(diffusion->eqn.flags);
  diffusion->eqn.ref_count = gkyl_ref_count_init(gkyl_conf_diffusion_free);

  struct conf_diffusion *diffusion_cu =
    (struct conf_diffusion*) gkyl_cu_malloc(sizeof(*diffusion_cu));
  gkyl_cu_memcpy(diffusion_cu, diffusion, sizeof(*diffusion_cu),
    GKYL_CU_MEMCPY_H2D);
  conf_diffusion_set_cu_dev_ptrs<<<1,1>>>(diffusion_cu, basis->ndim);
  diffusion->eqn.on_dev = &diffusion_cu->eqn;
  return &diffusion->eqn;
}
