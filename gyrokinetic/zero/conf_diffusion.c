#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_conf_diffusion.h>
#include <gkyl_conf_diffusion_kernels.h>
#include <gkyl_ref_count.h>
#include <gkyl_util.h>

struct conf_diffusion {
  struct gkyl_dg_eqn eqn;
  struct gkyl_basis basis;
  double (*vol)(const double *w, const double *dx, const double *K,
    const double *jacobgeo_inv, const double *f, double *out);
  double (*surf[GKYL_MAX_CDIM])(const double *wc, const double *dxc,
    const double *Dl, const double *Dc, const double *Dr,
    const double *fl, const double *fc, const double *fr, double *out);
  struct gkyl_range conf_range;
  struct gkyl_conf_diffusion_auxfields auxfields;
};

typedef double (*conf_diffusion_vol_kern_t)(const double*, const double*,
  const double*, const double*, const double*, double*);
typedef double (*conf_diffusion_surf_kern_t)(const double*, const double*,
  const double*, const double*, const double*, const double*, const double*,
  const double*, double*);

static const conf_diffusion_vol_kern_t ser_vol_kernels[] = {
  conf_diffusion_vol_1x_ser_p1,
  conf_diffusion_vol_2x_ser_p1,
  conf_diffusion_vol_3x_ser_p1,
};

static const conf_diffusion_surf_kern_t ser_surf_kernels[][GKYL_MAX_CDIM] = {
  { conf_diffusion_surfx_1x_ser_p1, 0, 0 },
  { conf_diffusion_surfx_2x_ser_p1, conf_diffusion_surfy_2x_ser_p1, 0 },
  { conf_diffusion_surfx_3x_ser_p1, conf_diffusion_surfy_3x_ser_p1,
    conf_diffusion_surfz_3x_ser_p1 },
};

static inline const double*
fetch_tensor(const struct conf_diffusion *diffusion, const int *idx)
{
  return gkyl_array_cfetch(diffusion->auxfields.diffusion_tensor,
    gkyl_range_idx(&diffusion->conf_range, idx));
}

static inline const double*
fetch_jacobgeo_inv(const struct conf_diffusion *diffusion, const int *idx)
{
  return gkyl_array_cfetch(diffusion->auxfields.jacobgeo_inv,
    gkyl_range_idx(&diffusion->conf_range, idx));
}

static double
conf_diffusion_vol(const struct gkyl_dg_eqn *eqn,
  const double *xc, const double *dx, const int *idx,
  const double *qIn, double *GKYL_RESTRICT qRhsOut)
{
  const struct conf_diffusion *diffusion =
    container_of(eqn, struct conf_diffusion, eqn);
  return diffusion->vol(xc, dx, fetch_tensor(diffusion, idx),
    fetch_jacobgeo_inv(diffusion, idx), qIn, qRhsOut);
}

static double
conf_diffusion_surf(const struct gkyl_dg_eqn *eqn, int dir,
  const double *xcL, const double *xcC, const double *xcR,
  const double *dxL, const double *dxC, const double *dxR,
  const int *idxL, const int *idxC, const int *idxR,
  const double *qInL, const double *qInC, const double *qInR,
  double *GKYL_RESTRICT qRhsOut)
{
  const struct conf_diffusion *diffusion =
    container_of(eqn, struct conf_diffusion, eqn);
  if (dir < 0 || dir >= GKYL_MAX_CDIM || diffusion->surf[dir] == 0)
    return 0.0;

  return diffusion->surf[dir](xcC, dxC,
    fetch_tensor(diffusion, idxL), fetch_tensor(diffusion, idxC),
    fetch_tensor(diffusion, idxR), qInL, qInC, qInR, qRhsOut);
}

static double
conf_diffusion_boundary_surf(const struct gkyl_dg_eqn *eqn, int dir,
  const double *xcEdge, const double *xcSkin,
  const double *dxEdge, const double *dxSkin,
  const int *idxEdge, const int *idxSkin, int edge,
  const double *qInEdge, const double *qInSkin,
  double *GKYL_RESTRICT qRhsOut)
{
  const struct conf_diffusion *diffusion =
    container_of(eqn, struct conf_diffusion, eqn);
  if (dir < 0 || dir >= GKYL_MAX_CDIM || diffusion->surf[dir] == 0)
    return 0.0;

  const int nb = diffusion->basis.num_basis;
  const int nK = diffusion->basis.ndim*diffusion->basis.ndim;
  const double *Kedge = fetch_tensor(diffusion, idxEdge);
  const double *Kskin = fetch_tensor(diffusion, idxSkin);
  double qghost[nb], Kghost[nK*nb];

  // Homogeneous natural boundary condition n_i K^{ij} d_j f = 0. The
  // interior surface kernel already contains the contribution from both cell
  // faces. Reflect f and use a sign-reversed reflected K on the exterior face;
  // its arithmetic face value is then exactly zero, leaving only the interior
  // interface contribution in the skin-cell update.
  diffusion->basis.flip_odd_sign(dir, qInSkin, qghost);
  for (int k=0; k<nK; ++k)
    diffusion->basis.flip_even_sign(dir, &Kskin[k*nb], &Kghost[k*nb]);

  if (edge == -1)
    return diffusion->surf[dir](xcSkin, dxSkin, Kghost, Kskin, Kedge,
      qghost, qInSkin, qInEdge, qRhsOut);
  return diffusion->surf[dir](xcSkin, dxSkin, Kedge, Kskin, Kghost,
    qInEdge, qInSkin, qghost, qRhsOut);
}

static void
conf_diffusion_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *eqn = container_of(ref, struct gkyl_dg_eqn, ref_count);
  struct conf_diffusion *diffusion =
    container_of(eqn, struct conf_diffusion, eqn);
  gkyl_free(diffusion);
}

void
gkyl_conf_diffusion_set_auxfields(const struct gkyl_dg_eqn *eqn,
  struct gkyl_conf_diffusion_auxfields auxin)
{
  struct conf_diffusion *diffusion =
    container_of(eqn, struct conf_diffusion, eqn);
  diffusion->auxfields = auxin;
}

struct gkyl_dg_eqn*
gkyl_conf_diffusion_new(const struct gkyl_basis *basis,
  const struct gkyl_range *conf_range, bool use_gpu)
{
  assert(!use_gpu);
  assert(basis->ndim > 0 && basis->ndim <= 3);
  assert(basis->poly_order == 1);
  assert(basis->b_type == GKYL_BASIS_MODAL_SERENDIPITY);

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
  for (int d=0; d<GKYL_MAX_CDIM; ++d)
    diffusion->surf[d] = ser_surf_kernels[basis->ndim-1][d];
  diffusion->eqn.flags = 0;
  diffusion->eqn.ref_count = gkyl_ref_count_init(conf_diffusion_free);
  diffusion->eqn.on_dev = &diffusion->eqn;
  diffusion->conf_range = *conf_range;
  diffusion->auxfields.diffusion_tensor = 0;
  diffusion->auxfields.jacobgeo_inv = 0;
  return &diffusion->eqn;
}
