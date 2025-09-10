#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_cx.h>
#include <gkyl_dg_cx_priv.h>
#include <gkyl_util.h>

gkyl_dg_cx*
gkyl_dg_cx_new(struct gkyl_dg_cx_inp *inp, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_dg_cx_cu_dev_new(inp);
  } 
#endif  
  gkyl_dg_cx *up = gkyl_malloc(sizeof(struct gkyl_dg_cx));

  up->cbasis = inp->cbasis;
  up->pbasis_gk = inp->pbasis_gk;
  up->pbasis_vl = inp->pbasis_vl;
  up->conf_rng = inp->conf_rng;
  up->conf_rng_ext = inp->conf_rng_ext;
  up->phase_rng = inp->phase_rng;
  up->grid = inp->grid;
  up->mass_ion = inp->mass_ion;
  up->mass_neut = inp->mass_neut;
  up->type_ion = inp->type_ion;
  up->vt_sq_ion_min = inp->vt_sq_ion_min;
  up->vt_sq_neut_min = inp->vt_sq_neut_min;

  int cdim = up->cbasis->ndim;
  int poly_order = up->cbasis->poly_order;
  int vdim_vl = up->pbasis_vl->ndim - cdim;
  enum gkyl_basis_type b_type = up->pbasis_vl->b_type;

  // get parameters for fitting function
  fit_param(up->type_ion, &up->a, &up->b);

  int tblidx = cv_index[cdim].vdim[vdim_vl];
  assert(tblidx != -1);
  up->react_rate = choose_kern(b_type, tblidx, poly_order);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up;
  
  return up;
}

void gkyl_dg_cx_coll(const struct gkyl_dg_cx *up, 
  struct gkyl_array *prim_vars_ion, struct gkyl_array *prim_vars_neut,
  struct gkyl_array *upar_b_i, struct gkyl_array *coef_cx, struct gkyl_array *cflrate)
{
#ifdef GKYL_HAVE_CUDA
  if(gkyl_array_is_cu_dev(coef_cx)) {
    return gkyl_dg_cx_coll_cu(up, prim_vars_ion, prim_vars_neut, 
      upar_b_i, coef_cx, cflrate);
  }
#endif
  
  struct gkyl_range vel_rng;
  struct gkyl_range_iter conf_iter, vel_iter;

  gkyl_range_iter_init(&conf_iter, up->conf_rng);
  while (gkyl_range_iter_next(&conf_iter)) {
    long loc = gkyl_range_idx(up->conf_rng, conf_iter.idx);

    const double *prim_vars_ion_d = gkyl_array_cfetch(prim_vars_ion, loc);
    const double *prim_vars_neut_d = gkyl_array_cfetch(prim_vars_neut, loc);
    const double *upar_b_i_d = gkyl_array_cfetch(upar_b_i, loc);

    double *coef_cx_d = gkyl_array_fetch(coef_cx, loc);
    
    double cflr = up->react_rate(up->a, up->b, up->vt_sq_ion_min, up->vt_sq_neut_min,
      prim_vars_ion_d, prim_vars_neut_d, upar_b_i_d, coef_cx_d);
  }
}

void
gkyl_dg_cx_release(gkyl_dg_cx* cx)
{
  free(cx);
}
