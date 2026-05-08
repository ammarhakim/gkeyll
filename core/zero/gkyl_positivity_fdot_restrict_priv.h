#pragma once

#include <gkyl_positivity_fdot_restrict.h>

struct gkyl_positivity_fdot_restrict {
  const struct gkyl_basis *basis;
  enum gkyl_positivity_fdot_restrict_type type;
  double safety_factor;
  bool use_gpu;
  struct gkyl_positivity_fdot_restrict *on_dev;

  // Arrays for quadrature point evaluations (CPU)
  double *fquad;
  double *dfdt_quad;

  void (*fdot_restriction_func)(struct gkyl_positivity_fdot_restrict *up,
    const double *f, double *dfdt, double dt, double *fq, double *dfdtq);
};

#ifdef GKYL_HAVE_CUDA
void gkyl_positivity_fdot_restrict_cu_set_ptrs(struct gkyl_positivity_fdot_restrict *up);
void gkyl_positivity_fdot_restrict_advance_cu(gkyl_positivity_fdot_restrict *up,
  const struct gkyl_range *range, const struct gkyl_array *f,
  struct gkyl_array *dfdt, double dt);
#endif

GKYL_CU_DH
static void
restrict_fdot_avg(struct gkyl_positivity_fdot_restrict *up,
  const double *f, double *dfdt, double dt, double *fq, double *dfdtq)
{
  double f_pred = f[0] + dt * dfdt[0];
  if (f_pred < 0.0)
    dfdt[0] = (up->safety_factor - 1.0) * f[0] / dt;
}

GKYL_CU_DH
static void
restrict_fdot_quad(struct gkyl_positivity_fdot_restrict *up,
  const double *f, double *dfdt, double dt, double *fq, double *dfdtq)
{
  const int nquad = up->basis->num_quad;
  const double safety_factor = up->safety_factor;
  bool reproject = false;

  for (int k = 0; k < nquad; ++k) {
    up->basis->modal_to_quad_nodal(dfdt, dfdtq, k);
    up->basis->modal_to_quad_nodal(f, fq, k);

    if (fq[k] + dt * dfdtq[k] < 0.0) {
      dfdtq[k] = (safety_factor - 1.0) * fq[k] / dt;
      reproject = true;
    }
  }

  if (reproject) {
    for (int k = 0; k < nquad; ++k)
      up->basis->quad_nodal_to_modal(dfdtq, dfdt, k);
  }
}

GKYL_CU_DH
static void
restrict_fdot_diode_avg(struct gkyl_positivity_fdot_restrict *up,
  const double *f, double *dfdt, double dt, double *fq, double *dfdtq)
{
  if (f[0] < 0.0 && dfdt[0] < 0.0)
    dfdt[0] = 0.0;
}

GKYL_CU_DH
static void
restrict_fdot_diode_quad(struct gkyl_positivity_fdot_restrict *up,
  const double *f, double *dfdt, double dt, double *fq, double *dfdtq)
{
  const int nquad = up->basis->num_quad;
  bool reproject = false;

  for (int k = 0; k < nquad; ++k) {
    up->basis->modal_to_quad_nodal(f, fq, k);
    up->basis->modal_to_quad_nodal(dfdt, dfdtq, k);

    if (fq[k] < 0.0 && dfdtq[k] < 0.0) {
      dfdtq[k] = 0.0;
      reproject = true;
    }
  }

  if (reproject) {
    for (int k = 0; k < nquad; ++k)
      up->basis->quad_nodal_to_modal(dfdtq, dfdt, k);
  }
}
