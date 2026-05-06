#include <gkyl_positivity.h>
#include <gkyl_positivity_priv.h>

#include <gkyl_alloc.h>

#include <float.h>
#include <math.h>
#include <stdio.h>

static double
cell_average(const struct gkyl_positivity *up, const double *f)
{
  return up->cellav_fac * f[0];
}

// Finds the minimum and maximum values of the quadrature points of the cell.
static void
min_max_quad(const struct gkyl_positivity *up, const double *f, double *fmin, double *fmax)
{
  const int nbasis = up->basis.num_basis;

  *fmin = DBL_MAX;
  *fmax = -DBL_MAX;
  for (int k = 0; k < nbasis; ++k) {
    up->fquad[k] = 0.0;
    up->basis.modal_to_quad_nodal(f, up->fquad, k);
    if (up->fquad[k] < *fmin) *fmin = up->fquad[k];
    if (up->fquad[k] > *fmax) *fmax = up->fquad[k];
  }
}

static void
min_max_quad_neighbor(const struct gkyl_positivity *up, const struct gkyl_range *range,
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
limit_cell_zs(const struct gkyl_positivity *up, const struct gkyl_range *range,
  struct gkyl_array *f, const int *idx, double *fc)
{
  const int nbasis = up->basis.num_basis;

  double fmin = DBL_MAX;
  double fmax = -DBL_MAX;
  min_max_quad(up, fc, &fmin, &fmax);

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
limit_cell_mrs(const struct gkyl_positivity *up, const struct gkyl_range *range,
  struct gkyl_array *f, const int *idx, double *fc)
{
  const int nbasis = up->basis.num_basis;

  double qm, qM; // minimum and Maximum of quadrature points in cell
  min_max_quad(up, fc, &qm, &qM);

  if (qm >= 0.0) {
    return;
  }

  double qm_neighbor, qM_neighbor;
  min_max_quad_neighbor(up, range, f, idx, &qm_neighbor, &qM_neighbor);
  const double qbar = cell_average(up, fc);

  double alpha = (fmax(fmax(qM,qbar),qM_neighbor) - fmin(fmin(qm,qbar),qm_neighbor))*0.0;

  double Mi = fmax(qbar + alpha, qM_neighbor);
  double mi = fmax(fmin(qbar - alpha, qm_neighbor), 0.0);

  double theta_Mi = fmin(((Mi - qbar) / (qM - qbar))/1.1, 1.0);
  double theta_mi = fmin(((mi - qbar) / (qm - qbar))/1.1, 1.0);
  double theta = fmin(fmin(theta_Mi, theta_mi), 1.0);

  if (theta < 1.0) {
    for (int k = 1; k < nbasis; ++k) {
      fc[k] *= theta;
    }
  }
}

static void
limit_time_step(const struct gkyl_positivity *up, const double *fc,
  struct gkyl_array *dfdt, double *dt, long lidx, double *dt_bound)
{
  const double *dfc = gkyl_array_cfetch(dfdt, lidx);
  const double fbar = cell_average(up, fc);
  const double dfbar = cell_average(up, dfc);

  if (dfbar < 0.0) {
    const double cell_dt = up->dt_factor * fbar / (-dfbar);
    if (cell_dt < *dt_bound) *dt_bound = cell_dt;
  }
}

static void
per_cell_limiter(const struct gkyl_positivity *up, const struct gkyl_range *range,
  struct gkyl_array *f, struct gkyl_array *dfdt, double *dt)
{
  double dt_bound = DBL_MAX;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(range, iter.idx);
    double *fc = gkyl_array_fetch(f, lidx);

    up->limit_cell_func(up, range, f, iter.idx, fc);

    up->limit_timestep_func(up, fc, dfdt, &dt_bound, lidx, &dt_bound);
  }
  *dt = dt_bound;
}

void
limiter_func_diabled(const struct gkyl_positivity *up, const struct gkyl_range *range,
  struct gkyl_array *f, struct gkyl_array *dfdt, double *dt)
{}

void
limit_cell_disabled(const struct gkyl_positivity *up, const struct gkyl_range *range,
  struct gkyl_array *f, const int *idx, double *fc)
{}

struct gkyl_positivity*
gkyl_positivity_new(struct gkyl_positivity_inp inp)
{
  struct gkyl_positivity *up = gkyl_malloc(sizeof(*up));
  up->type = inp.type;
  up->basis = inp.basis;
  up->dt_factor = inp.dt_factor;
  up->cellav_fac = 1. / pow(sqrt(2.), inp.basis.ndim);
  up->fquad = gkyl_malloc(up->basis.num_basis * sizeof(double));

  switch (inp.type) {
    case GKYL_POSITIVITY_ZS:
      up->positivity_func = per_cell_limiter;
      up->limit_cell_func = limit_cell_zs;
      up->limit_timestep_func = limit_time_step;
      break;
    case GKYL_POSITIVITY_MRS:
      up->positivity_func = per_cell_limiter;
      up->limit_cell_func = limit_cell_mrs;
      up->limit_timestep_func = limit_time_step;
      break;
    default:
      up->positivity_func = limiter_func_diabled;
      up->limit_cell_func = limit_cell_disabled;
      up->limit_timestep_func = limit_time_step;
  }

  return up;
}

void
gkyl_positivity_advance(gkyl_positivity *up, const struct gkyl_range *range, 
  struct gkyl_array *f, struct gkyl_array *dfdt, double *dt)
{
  up->positivity_func(up, range, f, dfdt, dt);
}

void
gkyl_positivity_release(gkyl_positivity *up)
{
  gkyl_free(up->fquad);
  gkyl_free(up);
}