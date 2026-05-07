#include <gkyl_positivity_fdot_restrict.h>
#include <gkyl_positivity_fdot_restrict_priv.h>

#include <gkyl_alloc.h>
#include <gkyl_range.h>

#include <float.h>
#include <math.h>

struct gkyl_positivity_fdot_restrict*
gkyl_positivity_fdot_restrict_new(struct gkyl_positivity_fdot_restrict_inp inp)
{
  struct gkyl_positivity_fdot_restrict *up = gkyl_malloc(sizeof(*up));
  
  up->basis = inp.basis; // Pass device basis
  up->type = inp.type;
  up->safety_factor = inp.safety_factor;
  up->use_gpu = inp.use_gpu;

  switch (inp.type) {
    case GKYL_POSITIVITY_FDOT_RESTRICT_QUAD:
      up->fdot_restriction_func = restrict_fdot_quad;
      break;
    case GKYL_POSITIVITY_FDOT_RESTRICT_AVG:
      up->fdot_restriction_func = restrict_fdot_avg;
      break;
    case GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD:
      up->fdot_restriction_func = restrict_fdot_diode_quad;
      break;
    case GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG:
      up->fdot_restriction_func = restrict_fdot_diode_avg;
      break;
  }

#ifdef GKYL_HAVE_CUDA
  up->on_dev = NULL;
  if (up->use_gpu) {
    up->on_dev = gkyl_cu_malloc(sizeof(*up));
    gkyl_cu_memcpy(up->on_dev, up, sizeof(*up), GKYL_CU_MEMCPY_H2D);
  }
#endif
  
  return up;
}

void
gkyl_positivity_fdot_restrict_advance(gkyl_positivity_fdot_restrict* up,
  const struct gkyl_range *range, const struct gkyl_array *f,
  struct gkyl_array *dfdt, double dt)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    return gkyl_positivity_fdot_restrict_advance_cu(up, range, f, dfdt, dt);
  }
#endif
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(range, iter.idx);
    const double *f_ptr = gkyl_array_cfetch(f, idx);
    double *dfdt_ptr = gkyl_array_fetch(dfdt, idx);
    up->fdot_restriction_func(up, f_ptr, dfdt_ptr, dt);
  }
}

void
gkyl_positivity_fdot_restrict_release(gkyl_positivity_fdot_restrict* up)
{
  if (up) {
#ifdef GKYL_HAVE_CUDA
    if (up->on_dev)
      gkyl_cu_free(up->on_dev);
#endif
    gkyl_free(up);
  }
}
