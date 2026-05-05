#include <gkyl_limiter.h>
#include <gkyl_limiter_priv.h>

#include <gkyl_alloc.h>

#include <float.h>
#include <math.h>
#include <stdio.h>

static double
cell_average(const struct gkyl_limiter *up, const double *f)
{
  return up->cellav_fac * f[0];
}

// Finds the minimum and maximum values of the quadrature points of the cell. Used in MRS limiter
static void
min_max_quad(const struct gkyl_limiter *up, const double *f, double *fmin, double *fmax)
{
  const int nbasis = up->basis.num_basis;
  double fquad[nbasis];

  *fmin = DBL_MAX;
  *fmax = -DBL_MAX;
  for (int k = 0; k < nbasis; ++k) {
    fquad[k] = 0.0;
    up->basis.modal_to_quad_nodal(f, fquad, k);
    if (fquad[k] < *fmin) *fmin = fquad[k];
    if (fquad[k] > *fmax) *fmax = fquad[k];
  }
}

static void
min_max_quad_neighbor(const struct gkyl_limiter *up, const struct gkyl_range *range,
  struct gkyl_array *f, const int *f_idx, double *fmin, double *fmax)
{
  *fmin = DBL_MAX;
  *fmax = -DBL_MAX;

  int idx[GKYL_MAX_DIM];
  for (int d = 0; d < up->basis.ndim; ++d) {
    idx[d] = f_idx[d];
  }

  for (int d = 0; d < up->basis.ndim; ++d) {
    // check lower neighbor
    if (idx[d] > 0) {
      idx[d] -= 1;
      const double *f_neighbor = gkyl_array_cfetch(f, gkyl_range_idx(range, idx));
      double fmin_neighbor, fmax_neighbor;
      min_max_quad(up, f_neighbor, &fmin_neighbor, &fmax_neighbor);
      if (fmin_neighbor < *fmin) *fmin = fmin_neighbor;
      if (fmax_neighbor > *fmax) *fmax = fmax_neighbor;
      idx[d] += 1;
    }

    // check upper neighbor
    if (idx[d] < range->upper[d]) {
      idx[d] += 1;
      const double *f_neighbor = gkyl_array_cfetch(f, gkyl_range_idx(range, idx));
      double fmin_neighbor, fmax_neighbor;
      min_max_quad(up, f_neighbor, &fmin_neighbor, &fmax_neighbor);
      if (fmin_neighbor < *fmin) *fmin = fmin_neighbor;
      if (fmax_neighbor > *fmax) *fmax = fmax_neighbor;
      idx[d] -= 1;
    }
  }
}

static void
limit_cell_zs(const struct gkyl_limiter *up, const struct gkyl_range *range,
  struct gkyl_array *f, const int *idx, double *fc)
{
  (void) range;
  (void) f;
  (void) idx;

  const int nbasis = up->basis.num_basis;
  double fquad[nbasis];

  double fmin = DBL_MAX;
  for (int k = 0; k < nbasis; ++k) {
    fquad[k] = 0.0;
    up->basis.modal_to_quad_nodal(fc, fquad, k);
    fmin = fmin < fquad[k] ? fmin : fquad[k];
  }

  if (fmin >= 0.0) {
    return;
  }

  const double fbar = cell_average(up, fc);

  const double theta = fbar / (fbar - fmin);
  if (theta < 1.0) {
    for (int k = 1; k < nbasis; ++k) {
      fc[k] *= theta;
    }
  }
}

static void
limit_cell_mrs(const struct gkyl_limiter *up, const struct gkyl_range *range,
  struct gkyl_array *f, const int *idx, double *fc)
{
  const int nbasis = up->basis.num_basis;
  double fquad[nbasis];

  double fmin, fmax;
  min_max_quad(up, fc, &fmin, &fmax);
  double fmin_neighbor, fmax_neighbor;
  min_max_quad_neighbor(up, range, f, idx, &fmin_neighbor, &fmax_neighbor);
  const double fbar = cell_average(up, fc);

  double theta = 1.0;
  if (fmin < 0.0) {
    theta = fbar / (fbar - fmin);
  }

  if (fmax > fbar && fmax_neighbor > fbar) {
    const double theta_max = (fmax_neighbor - fbar) / (fmax - fbar);
    if (theta_max < theta) {
      theta = theta_max;
    }
  }

  if (theta < 1.0) {
    if (theta < 0.0) {
      theta = 0.0;
    }
    for (int k = 1; k < nbasis; ++k) {
      fc[k] *= theta;
    }
  }
}

static void
per_cell_limiter(const struct gkyl_limiter *up, const struct gkyl_range *range,
  struct gkyl_array *f)
{
  assert(f->ncomp == up->basis.num_basis);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(range, iter.idx);
    double *fc = gkyl_array_fetch(f, lidx);
    up->limit_cell_func(up, range, f, iter.idx, fc);
  }
}

static void
per_cell_limiter_timestep(const struct gkyl_limiter *up, const struct gkyl_range *range,
  struct gkyl_array *f, struct gkyl_array *dfdt, double *dt)
{
  assert(f->ncomp == up->basis.num_basis);
  assert(dfdt->ncomp == up->basis.num_basis);

  double dt_bound = DBL_MAX;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(range, iter.idx);
    const double *fc = gkyl_array_cfetch(f, lidx);
    const double *dfc = gkyl_array_cfetch(dfdt, lidx);

    const double fbar = cell_average(up, fc);
    const double dfbar = cell_average(up, dfc);

    if (dfbar < 0.0) {
      const double cell_dt = up->dt_factor * fbar / (-dfbar);
      if (cell_dt < dt_bound) {
        dt_bound = cell_dt;
      }
    }
    else if (fbar <= 0.0) {
      dt_bound = 0.0;
      break;
    }
  }

  if (dt_bound < 0.0) {
    dt_bound = 0.0;
  }

  *dt = dt_bound;
}

struct gkyl_limiter*
gkyl_limiter_new(struct gkyl_limiter_inp inp)
{
  assert(inp.basis.modal_to_quad_nodal != NULL);
  assert(inp.basis.quad_nodal_to_modal != NULL);

  struct gkyl_limiter *up = gkyl_malloc(sizeof(*up));
  up->basis = inp.basis;
  up->cellav_fac = 1. / pow(sqrt(2.), inp.basis.ndim);
  up->dt_factor = inp.dt_factor;

  switch (inp.type) {
    case GKYL_LIMITER_ZS:
      up->limiter_func = per_cell_limiter;
      up->limiter_timestep_func = per_cell_limiter_timestep;
      up->limit_cell_func = limit_cell_zs;
      break;
    case GKYL_LIMITER_MRS:
      up->limiter_func = per_cell_limiter;
      up->limiter_timestep_func = per_cell_limiter_timestep;
      up->limit_cell_func = limit_cell_mrs;
      break;
    default:
      assert(false);
  }

  return up;
}

void
gkyl_limiter_advance(gkyl_limiter *up, const struct gkyl_range *range, struct gkyl_array *f)
{
  up->limiter_func(up, range, f);
}

void
gkyl_limiter_advance_timestep(gkyl_limiter *up, const struct gkyl_range *range, struct gkyl_array *f, struct gkyl_array *dfdt, double *dt)
{
  up->limiter_timestep_func(up, range, f, dfdt, dt);
}

void
gkyl_limiter_release(gkyl_limiter *up)
{
  gkyl_free(up);
}