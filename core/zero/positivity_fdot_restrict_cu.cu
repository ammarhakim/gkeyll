extern "C" {
#include <gkyl_positivity_fdot_restrict.h>
#include <gkyl_positivity_fdot_restrict_priv.h>

#include <gkyl_alloc.h>
#include <gkyl_range.h>
}

__global__ static void
gkyl_positivity_fdot_restrict_set_func_cu_kern(gkyl_positivity_fdot_restrict *up)
{
  switch (up->type) {
    case GKYL_POSITIVITY_FDOT_RESTRICT_AVG:
      up->fdot_restriction_func = restrict_fdot_avg;
      break;
    case GKYL_POSITIVITY_FDOT_RESTRICT_QUAD:
      up->fdot_restriction_func = restrict_fdot_quad;
      break;
    case GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG:
      up->fdot_restriction_func = restrict_fdot_diode_avg;
      break;
    case GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_QUAD:
      up->fdot_restriction_func = restrict_fdot_diode_quad;
      break;
    case GKYL_POSITIVITY_FDOT_RESTRICT_AVG_ALL_COMP:
      up->fdot_restriction_func = restrict_fdot_avg_all_comp;
      break;
    case GKYL_POSITIVITY_FDOT_RESTRICT_DIODE_AVG_ALL_COMP:
      up->fdot_restriction_func = restrict_fdot_diode_avg_all_comp;
      break;
  }
}

// Upper bound on quadrature nodes. GPUs cannot dynamically allocate memory
#define GKYL_POS_FDOT_MAX_QUAD 243

__global__ static void
gkyl_positivity_fdot_restrict_cu_kern(gkyl_positivity_fdot_restrict *up,
  struct gkyl_range range, const struct gkyl_array *f, struct gkyl_array *dfdt, double dt)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= range.volume) return;
  double fq[GKYL_POS_FDOT_MAX_QUAD], dfdtq[GKYL_POS_FDOT_MAX_QUAD];
  up->fdot_restriction_func(up, (const double*)gkyl_array_cfetch(f, idx),
    (double*)gkyl_array_fetch(dfdt, idx), dt, fq, dfdtq);
}

void
gkyl_positivity_fdot_restrict_cu_set_ptrs(struct gkyl_positivity_fdot_restrict *up)
{
  gkyl_positivity_fdot_restrict_set_func_cu_kern<<<1,1>>>(up->on_dev);
}

void
gkyl_positivity_fdot_restrict_advance_cu(gkyl_positivity_fdot_restrict *up,
  const struct gkyl_range *range, const struct gkyl_array *f,
  struct gkyl_array *dfdt, double dt)
{
  gkyl_positivity_fdot_restrict_cu_kern<<<range->nblocks, range->nthreads>>>(
    up->on_dev, *range, f->on_dev, dfdt->on_dev, dt);
}
