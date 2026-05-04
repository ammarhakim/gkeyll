#include <gkyl_positivity_shift_zhang_shu.h>
#include <gkyl_positivity_shift_zhang_shu_priv.h>

#include <gkyl_alloc.h>

#include <float.h>
#include <math.h>

struct gkyl_positivity_shift_zhang_shu*
gkyl_positivity_shift_zhang_shu_new(struct gkyl_basis basis)
{
  assert(basis.modal_to_quad_nodal != NULL);
  assert(basis.quad_nodal_to_modal != NULL);

  struct gkyl_positivity_shift_zhang_shu *up = gkyl_malloc(sizeof(*up));
  up->basis = basis;
  up->cellav_fac = 1./pow(sqrt(2.), basis.ndim);

  return up;
}

static bool
limit_cell(const struct gkyl_positivity_shift_zhang_shu *up, double *f)
{
  const int nbasis = up->basis.num_basis;
  double fquad[nbasis];

  double fmin = DBL_MAX;
  for (int k=0; k<nbasis; ++k) {
    fquad[k] = 0.0;
    up->basis.modal_to_quad_nodal(f, fquad, k);
    fmin = fmin < fquad[k] ? fmin : fquad[k];
  }

  if (fmin >= 0.0)
    return true;

  double fbar = up->cellav_fac * f[0];
  if (fbar <= 0.0)
    return false;

  double theta = fbar / (fbar - fmin);
  if (theta < 1.0) {
    for (int k=1; k<nbasis; ++k)
      f[k] *= theta;
  }

  return true;
}

bool
gkyl_positivity_shift_zhang_shu_advance(gkyl_positivity_shift_zhang_shu* up,
  const struct gkyl_range *range, struct gkyl_array *f)
{
  assert(f->type == GKYL_DOUBLE);
  assert(f->ncomp == up->basis.num_basis);

  bool success = true;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(range, iter.idx);
    double *fc = gkyl_array_fetch(f, lidx);
    success = limit_cell(up, fc) && success;
  }

  return success;
}

void
gkyl_positivity_shift_zhang_shu_release(gkyl_positivity_shift_zhang_shu* up)
{
  gkyl_free(up);
}