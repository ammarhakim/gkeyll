#pragma once

#include <gkyl_conf_diffusion.h>
#include <gkyl_conf_diffusion_kernels.h>
#include <gkyl_ref_count.h>

typedef double (*conf_diffusion_vol_kern_t)(const double*, const double*,
  const double*, const double*, const double*, double*);
typedef double (*conf_diffusion_surf_kern_t)(const double*, const double*,
  const double*, const double*, const double*, const double*, const double*,
  const double*, double*);
typedef double (*conf_diffusion_boundary_surf_kern_t)(const double*,
  const double*, const double*, const double*, int, const double*,
  const double*, double*);

struct conf_diffusion {
  struct gkyl_dg_eqn eqn;
  struct gkyl_basis basis;
  conf_diffusion_vol_kern_t vol;
  // One surface kernel for each flux/derivative direction pair (i,j),
  // following the general tensor-diffusion operator's decomposition.
  conf_diffusion_surf_kern_t surf[GKYL_MAX_CDIM][GKYL_MAX_CDIM];
  conf_diffusion_boundary_surf_kern_t
    boundary_surf[GKYL_MAX_CDIM][GKYL_MAX_CDIM];
  struct gkyl_range conf_range;
  struct gkyl_conf_diffusion_auxfields auxfields;
};

GKYL_CU_DH static inline const double*
conf_diffusion_fetch_tensor(const struct conf_diffusion *diffusion,
  const int *idx)
{
  return (const double*) gkyl_array_cfetch(diffusion->auxfields.diffusion_tensor,
    gkyl_range_idx(&diffusion->conf_range, idx));
}

GKYL_CU_DH static inline const double*
conf_diffusion_fetch_jacobgeo_inv(const struct conf_diffusion *diffusion,
  const int *idx)
{
  return (const double*) gkyl_array_cfetch(diffusion->auxfields.jacobgeo_inv,
    gkyl_range_idx(&diffusion->conf_range, idx));
}

GKYL_CU_DH static double
conf_diffusion_vol(const struct gkyl_dg_eqn *eqn,
  const double *xc, const double *dx, const int *idx,
  const double *qIn, double *GKYL_RESTRICT qRhsOut)
{
  const struct conf_diffusion *diffusion =
    container_of(eqn, struct conf_diffusion, eqn);
  return diffusion->vol(xc, dx, conf_diffusion_fetch_tensor(diffusion, idx),
    conf_diffusion_fetch_jacobgeo_inv(diffusion, idx), qIn, qRhsOut);
}

GKYL_CU_DH static double
conf_diffusion_surf(const struct gkyl_dg_eqn *eqn, int dir,
  const double *xcL, const double *xcC, const double *xcR,
  const double *dxL, const double *dxC, const double *dxR,
  const int *idxL, const int *idxC, const int *idxR,
  const double *qInL, const double *qInC, const double *qInR,
  double *GKYL_RESTRICT qRhsOut)
{
  const struct conf_diffusion *diffusion =
    container_of(eqn, struct conf_diffusion, eqn);
  if (dir < 0 || dir >= diffusion->basis.ndim)
    return 0.0;
  const double *Kl = conf_diffusion_fetch_tensor(diffusion, idxL);
  const double *Kc = conf_diffusion_fetch_tensor(diffusion, idxC);
  const double *Kr = conf_diffusion_fetch_tensor(diffusion, idxR);
  double cflFreq = 0.0;
  for (int deriv_dir=0; deriv_dir<diffusion->basis.ndim; ++deriv_dir)
    cflFreq += diffusion->surf[dir][deriv_dir](xcC, dxC, Kl, Kc, Kr,
      qInL, qInC, qInR, qRhsOut);
  return cflFreq;
}

GKYL_CU_DH static double
conf_diffusion_boundary_surf(const struct gkyl_dg_eqn *eqn, int dir,
  const double *xcEdge, const double *xcSkin,
  const double *dxEdge, const double *dxSkin,
  const int *idxEdge, const int *idxSkin, int edge,
  const double *qInEdge, const double *qInSkin,
  double *GKYL_RESTRICT qRhsOut)
{
  const struct conf_diffusion *diffusion =
    container_of(eqn, struct conf_diffusion, eqn);
  if (dir < 0 || dir >= diffusion->basis.ndim)
    return 0.0;

  const double *Kedge = conf_diffusion_fetch_tensor(diffusion, idxEdge);
  const double *Kskin = conf_diffusion_fetch_tensor(diffusion, idxSkin);
  double cflFreq = 0.0;
  for (int deriv_dir=0; deriv_dir<diffusion->basis.ndim; ++deriv_dir)
    cflFreq += diffusion->boundary_surf[dir][deriv_dir](xcSkin, dxSkin,
      Kedge, Kskin, edge, qInEdge, qInSkin, qRhsOut);
  return cflFreq;
}

void gkyl_conf_diffusion_free(const struct gkyl_ref_count *ref);

#ifdef GKYL_HAVE_CUDA
struct gkyl_dg_eqn* gkyl_conf_diffusion_cu_dev_new(
  const struct gkyl_basis *basis, const struct gkyl_range *conf_range);
void gkyl_conf_diffusion_set_auxfields_cu(const struct gkyl_dg_eqn *eqn,
  struct gkyl_conf_diffusion_auxfields auxin);
#endif
