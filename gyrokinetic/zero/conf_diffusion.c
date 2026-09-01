#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_conf_diffusion.h>
#include <gkyl_conf_diffusion_priv.h>
#include <gkyl_ref_count.h>
#include <gkyl_util.h>

static const conf_diffusion_vol_kern_t ser_vol_kernels[] = {
  conf_diffusion_vol_1x_ser_p1,
  conf_diffusion_vol_2x_ser_p1,
  conf_diffusion_vol_3x_ser_p1,
};

static const conf_diffusion_surf_kern_t
ser_surf_kernels[][GKYL_MAX_CDIM][GKYL_MAX_CDIM] = {
  { { conf_diffusion_surfxx_1x_ser_p1, 0, 0 }, { 0 }, { 0 } },
  { { conf_diffusion_surfxx_2x_ser_p1,
      conf_diffusion_surfxy_2x_ser_p1, 0 },
    { conf_diffusion_surfyx_2x_ser_p1,
      conf_diffusion_surfyy_2x_ser_p1, 0 }, { 0 } },
  { { conf_diffusion_surfxx_3x_ser_p1,
      conf_diffusion_surfxy_3x_ser_p1,
      conf_diffusion_surfxz_3x_ser_p1 },
    { conf_diffusion_surfyx_3x_ser_p1,
      conf_diffusion_surfyy_3x_ser_p1,
      conf_diffusion_surfyz_3x_ser_p1 },
    { conf_diffusion_surfzx_3x_ser_p1,
      conf_diffusion_surfzy_3x_ser_p1,
      conf_diffusion_surfzz_3x_ser_p1 } },
};

static const conf_diffusion_boundary_surf_kern_t
ser_boundary_surf_kernels[][GKYL_MAX_CDIM][GKYL_MAX_CDIM] = {
  { { conf_diffusion_boundary_surfxx_1x_ser_p1, 0, 0 }, { 0 }, { 0 } },
  { { conf_diffusion_boundary_surfxx_2x_ser_p1,
      conf_diffusion_boundary_surfxy_2x_ser_p1, 0 },
    { conf_diffusion_boundary_surfyx_2x_ser_p1,
      conf_diffusion_boundary_surfyy_2x_ser_p1, 0 }, { 0 } },
  { { conf_diffusion_boundary_surfxx_3x_ser_p1,
      conf_diffusion_boundary_surfxy_3x_ser_p1,
      conf_diffusion_boundary_surfxz_3x_ser_p1 },
    { conf_diffusion_boundary_surfyx_3x_ser_p1,
      conf_diffusion_boundary_surfyy_3x_ser_p1,
      conf_diffusion_boundary_surfyz_3x_ser_p1 },
    { conf_diffusion_boundary_surfzx_3x_ser_p1,
      conf_diffusion_boundary_surfzy_3x_ser_p1,
      conf_diffusion_boundary_surfzz_3x_ser_p1 } },
};

void
gkyl_conf_diffusion_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *eqn = container_of(ref, struct gkyl_dg_eqn, ref_count);
  if (gkyl_dg_eqn_is_cu_dev(eqn)) {
    struct conf_diffusion *diffusion_cu =
      container_of(eqn->on_dev, struct conf_diffusion, eqn);
    gkyl_cu_free(diffusion_cu);
  }
  struct conf_diffusion *diffusion =
    container_of(eqn, struct conf_diffusion, eqn);
  gkyl_free(diffusion);
}

void
gkyl_conf_diffusion_set_auxfields(const struct gkyl_dg_eqn *eqn,
  struct gkyl_conf_diffusion_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(auxin.diffusion_tensor)
    && gkyl_array_is_cu_dev(auxin.jacobgeo_inv)) {
    gkyl_conf_diffusion_set_auxfields_cu(eqn->on_dev, auxin);
    return;
  }
#endif
  struct conf_diffusion *diffusion =
    container_of(eqn, struct conf_diffusion, eqn);
  diffusion->auxfields = auxin;
}

struct gkyl_dg_eqn*
gkyl_conf_diffusion_new(const struct gkyl_basis *basis,
  const struct gkyl_range *conf_range, bool use_gpu)
{
  assert(basis->ndim > 0 && basis->ndim <= 3);
  assert(basis->poly_order == 1);
  assert(basis->b_type == GKYL_BASIS_MODAL_SERENDIPITY);

#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    return gkyl_conf_diffusion_cu_dev_new(basis, conf_range);
#endif
  assert(!use_gpu);

  struct conf_diffusion *diffusion = gkyl_malloc(sizeof(*diffusion));
  diffusion->basis = *basis;
  diffusion->eqn.num_equations = 1;
  diffusion->eqn.vol_term = conf_diffusion_vol;
  diffusion->eqn.surf_term = conf_diffusion_surf;
  diffusion->eqn.boundary_surf_term = conf_diffusion_boundary_surf;
  diffusion->eqn.boundary_diag_term = conf_diffusion_boundary_surf;
  diffusion->eqn.gen_surf_term = 0;
  diffusion->eqn.gen_boundary_surf_term = 0;
  diffusion->vol = ser_vol_kernels[basis->ndim-1];
  for (int i=0; i<GKYL_MAX_CDIM; ++i)
    for (int j=0; j<GKYL_MAX_CDIM; ++j) {
      diffusion->surf[i][j] = ser_surf_kernels[basis->ndim-1][i][j];
      diffusion->boundary_surf[i][j] =
        ser_boundary_surf_kernels[basis->ndim-1][i][j];
    }
  diffusion->eqn.flags = 0;
  diffusion->eqn.ref_count = gkyl_ref_count_init(gkyl_conf_diffusion_free);
  diffusion->eqn.on_dev = &diffusion->eqn;
  diffusion->conf_range = *conf_range;
  diffusion->auxfields.diffusion_tensor = 0;
  diffusion->auxfields.jacobgeo_inv = 0;
  return &diffusion->eqn;
}
