#include <gkyl_positivity_fdot_restrict.h>
#include <gkyl_positivity_fdot_restrict_priv.h>

#include <gkyl_alloc.h>
#include <gkyl_range.h>

#include <float.h>
#include <math.h>

__global__ static void
gkyl_positivity_fdot_restrict_advance_cu_kern(const gkyl_positivity_fdot_restrict* up,
  const struct gkyl_range *range, const struct gkyl_array *f,
  struct gkyl_array *dfdt, double dt)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= range->volume) return;

  const double *fptr = (const double*) gkyl_array_cfetch(f, idx);
  double *dfdtptr = (double*) gkyl_array_fetch(dfdt, idx);

  up->fdot_restrict_func(up, fptr, dfdtptr, dt);
}

void
gkyl_positivity_fdot_restrict_advance_cu(gkyl_positivity_fdot_restrict* up,
  const struct gkyl_range *range, const struct gkyl_array *f,
  struct gkyl_array *dfdt, double dt)
{
  gkyl_positivity_fdot_restrict_advance_cu_kern<<<range->nblocks, range->nthreads>>>(
    up->on_dev, range, f->on_dev, dfdt->on_dev, dt);
}
