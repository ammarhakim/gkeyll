#pragma once

#include <gkyl_gk_numerical_diffusion.h>
#include <gkyl_gk_numerical_diffusion_kernels.h>
#include <gkyl_ref_count.h>

// private header for use in diffusion DG equation object creation
// functions

static inline int diffdirs_linidx(const bool *has_diff_in_dir, int cdim) {
  // Compute the linear index into the array of volume kernels (one
  // kernel for each combination of diffusive directions).

  // Linear index into list of volume kernels.
  int dirs_bin_key[] = {1,2,4,8,16,32}; // Binary: 000001, 000010, 000100, 001000, 010000, 100000.
  int dirs_linidx = 0; // Binary 000000.
  for (int d=0; d<cdim; d++) {
     if (has_diff_in_dir[d]) dirs_linidx = dirs_linidx | dirs_bin_key[d];
  }
  dirs_linidx -= 1;
  return dirs_linidx;
}

// Types for various kernels
typedef double (*gk_nume_diff_surf_t)(const double *w, const double *dx,
  const double *nul, const double *nuc, const double *nur,
  const double *jacobgeo_invl, const double *jacobgeo_invc, const double *jacobgeo_invr,
  const double *Jfl, const double *Jfc, const double *Jfr, double* GKYL_RESTRICT out);

typedef double (*gk_nume_diff_boundary_surf_t)(const double *wSkin, const double *dxSkin,
  const double *nuEdge, const double *nuSkin, const double *jacobgeo_invEdge, const double *jacobgeo_invSkin,
  int edge, const double *JfEdge, const double *JfSkin, double* GKYL_RESTRICT out);

struct gk_numerical_diffusion {
  struct gkyl_dg_eqn eqn;
  gk_nume_diff_surf_t surf[GKYL_MAX_CDIM];
  gk_nume_diff_boundary_surf_t boundary_surf_lower[GKYL_MAX_CDIM], boundary_surf_upper[GKYL_MAX_CDIM];
  gk_nume_diff_boundary_surf_t boundary_diag_lower[GKYL_MAX_CDIM], boundary_diag_upper[GKYL_MAX_CDIM];
  struct gkyl_range conf_range;
  struct gkyl_gk_numerical_diffusion_auxfields auxfields;
  bool has_diff_in_dir[GKYL_MAX_CDIM];
};

#define _cfnu(idx) (const double *) gkyl_array_cfetch(gknd->auxfields.nu, gkyl_range_idx(&gknd->conf_range, idx))

#define _cfJacInv(idx) (const double *) gkyl_array_cfetch(gknd->auxfields.jacobgeo_inv, gkyl_range_idx(&gknd->conf_range, idx))

// for use in kernel tables
typedef struct { vol_termf_t kernels[7]; } gkyl_gk_numerical_diffusion_vol_kern_list_diffdir;
typedef struct { gkyl_gk_numerical_diffusion_vol_kern_list_diffdir list[2]; } gkyl_gk_numerical_diffusion_vol_kern_list_polyOrder;
typedef struct { gkyl_gk_numerical_diffusion_vol_kern_list_polyOrder list[2]; } gkyl_gk_numerical_diffusion_vol_kern_list;

// for use in kernel tables
typedef struct { gk_nume_diff_surf_t kernels[2]; } gkyl_gk_numerical_diffusion_surf_kernels_polyOrder;
typedef struct { gkyl_gk_numerical_diffusion_surf_kernels_polyOrder list[4]; } gkyl_gk_numerical_diffusion_surf_kernels_dir;
typedef struct { gkyl_gk_numerical_diffusion_surf_kernels_dir list[2]; } gkyl_gk_numerical_diffusion_surf_kern_list;

typedef struct { gk_nume_diff_boundary_surf_t kernels[2]; } gkyl_gk_numerical_diffusion_boundary_surf_kernels_polyOrder;
typedef struct { gkyl_gk_numerical_diffusion_boundary_surf_kernels_polyOrder list[4]; } gkyl_gk_numerical_diffusion_boundary_surf_kernels_dir;
typedef struct { gkyl_gk_numerical_diffusion_boundary_surf_kernels_dir list[2]; } gkyl_gk_numerical_diffusion_boundary_surf_kern_list;

// ............... Inhomogeneous (spatially varying) diffusivity ............... //

// Serendipity volume kernels
// Need to be separated like this for GPU build

// 1x 2nd order diffusion.
GKYL_CU_DH static double ker_gk_numerical_diffusion_order2_vol_1x1v_ser_p1_varnu_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct gk_numerical_diffusion* gknd = container_of(eqn, struct gk_numerical_diffusion, eqn);
  return gk_numerical_diffusion_order2_vol_1x1v_ser_p1_varnu_diffdirsx(xc, dx, _cfnu(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_gk_numerical_diffusion_order2_vol_1x2v_ser_p1_varnu_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct gk_numerical_diffusion* gknd = container_of(eqn, struct gk_numerical_diffusion, eqn);
  return gk_numerical_diffusion_order2_vol_1x2v_ser_p1_varnu_diffdirsx(xc, dx, _cfnu(idx), _cfJacInv(idx), qIn, qRhsOut);
}
// 2x 2nd order diffusion.
GKYL_CU_DH static double ker_gk_numerical_diffusion_order2_vol_2x2v_ser_p1_varnu_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct gk_numerical_diffusion* gknd = container_of(eqn, struct gk_numerical_diffusion, eqn);
  return gk_numerical_diffusion_order2_vol_2x2v_ser_p1_varnu_diffdirsx(xc, dx, _cfnu(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_gk_numerical_diffusion_order2_vol_2x2v_ser_p1_varnu_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct gk_numerical_diffusion* gknd = container_of(eqn, struct gk_numerical_diffusion, eqn);
  return gk_numerical_diffusion_order2_vol_2x2v_ser_p1_varnu_diffdirsz(xc, dx, _cfnu(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_gk_numerical_diffusion_order2_vol_2x2v_ser_p1_varnu_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct gk_numerical_diffusion* gknd = container_of(eqn, struct gk_numerical_diffusion, eqn);
  return gk_numerical_diffusion_order2_vol_2x2v_ser_p1_varnu_diffdirsxz(xc, dx, _cfnu(idx), _cfJacInv(idx), qIn, qRhsOut);
}
// 3x 2nd order diffusion.
GKYL_CU_DH static double ker_gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct gk_numerical_diffusion* gknd = container_of(eqn, struct gk_numerical_diffusion, eqn);
  return gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffdirsx(xc, dx, _cfnu(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct gk_numerical_diffusion* gknd = container_of(eqn, struct gk_numerical_diffusion, eqn);
  return gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffdirsy(xc, dx, _cfnu(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct gk_numerical_diffusion* gknd = container_of(eqn, struct gk_numerical_diffusion, eqn);
  return gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffdirsxy(xc, dx, _cfnu(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct gk_numerical_diffusion* gknd = container_of(eqn, struct gk_numerical_diffusion, eqn);
  return gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffdirsz(xc, dx, _cfnu(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct gk_numerical_diffusion* gknd = container_of(eqn, struct gk_numerical_diffusion, eqn);
  return gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffdirsxz(xc, dx, _cfnu(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffdirsyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct gk_numerical_diffusion* gknd = container_of(eqn, struct gk_numerical_diffusion, eqn);
  return gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffdirsyz(xc, dx, _cfnu(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffdirsxyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct gk_numerical_diffusion* gknd = container_of(eqn, struct gk_numerical_diffusion, eqn);
  return gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffdirsxyz(xc, dx, _cfnu(idx), _cfJacInv(idx), qIn, qRhsOut);
}


// 1x 4th order diffusion.
GKYL_CU_DH static double ker_gk_numerical_diffusion_order4_vol_1x1v_ser_p1_varnu_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct gk_numerical_diffusion* gknd = container_of(eqn, struct gk_numerical_diffusion, eqn);
  return gk_numerical_diffusion_order4_vol_1x1v_ser_p1_varnu_diffdirsx(xc, dx, _cfnu(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_gk_numerical_diffusion_order4_vol_1x2v_ser_p1_varnu_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct gk_numerical_diffusion* gknd = container_of(eqn, struct gk_numerical_diffusion, eqn);
  return gk_numerical_diffusion_order4_vol_1x2v_ser_p1_varnu_diffdirsx(xc, dx, _cfnu(idx), _cfJacInv(idx), qIn, qRhsOut);
}
// 2x 4th order diffusion.
GKYL_CU_DH static double ker_gk_numerical_diffusion_order4_vol_2x2v_ser_p1_varnu_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct gk_numerical_diffusion* gknd = container_of(eqn, struct gk_numerical_diffusion, eqn);
  return gk_numerical_diffusion_order4_vol_2x2v_ser_p1_varnu_diffdirsx(xc, dx, _cfnu(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_gk_numerical_diffusion_order4_vol_2x2v_ser_p1_varnu_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct gk_numerical_diffusion* gknd = container_of(eqn, struct gk_numerical_diffusion, eqn);
  return gk_numerical_diffusion_order4_vol_2x2v_ser_p1_varnu_diffdirsz(xc, dx, _cfnu(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_gk_numerical_diffusion_order4_vol_2x2v_ser_p1_varnu_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct gk_numerical_diffusion* gknd = container_of(eqn, struct gk_numerical_diffusion, eqn);
  return gk_numerical_diffusion_order4_vol_2x2v_ser_p1_varnu_diffdirsxz(xc, dx, _cfnu(idx), _cfJacInv(idx), qIn, qRhsOut);
}
// 3x 4th order diffusion.
GKYL_CU_DH static double ker_gk_numerical_diffusion_order4_vol_3x2v_ser_p1_varnu_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct gk_numerical_diffusion* gknd = container_of(eqn, struct gk_numerical_diffusion, eqn);
  return gk_numerical_diffusion_order4_vol_3x2v_ser_p1_varnu_diffdirsx(xc, dx, _cfnu(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_gk_numerical_diffusion_order4_vol_3x2v_ser_p1_varnu_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct gk_numerical_diffusion* gknd = container_of(eqn, struct gk_numerical_diffusion, eqn);
  return gk_numerical_diffusion_order4_vol_3x2v_ser_p1_varnu_diffdirsy(xc, dx, _cfnu(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_gk_numerical_diffusion_order4_vol_3x2v_ser_p1_varnu_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct gk_numerical_diffusion* gknd = container_of(eqn, struct gk_numerical_diffusion, eqn);
  return gk_numerical_diffusion_order4_vol_3x2v_ser_p1_varnu_diffdirsxy(xc, dx, _cfnu(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_gk_numerical_diffusion_order4_vol_3x2v_ser_p1_varnu_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct gk_numerical_diffusion* gknd = container_of(eqn, struct gk_numerical_diffusion, eqn);
  return gk_numerical_diffusion_order4_vol_3x2v_ser_p1_varnu_diffdirsz(xc, dx, _cfnu(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_gk_numerical_diffusion_order4_vol_3x2v_ser_p1_varnu_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct gk_numerical_diffusion* gknd = container_of(eqn, struct gk_numerical_diffusion, eqn);
  return gk_numerical_diffusion_order4_vol_3x2v_ser_p1_varnu_diffdirsxz(xc, dx, _cfnu(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_gk_numerical_diffusion_order4_vol_3x2v_ser_p1_varnu_diffdirsyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct gk_numerical_diffusion* gknd = container_of(eqn, struct gk_numerical_diffusion, eqn);
  return gk_numerical_diffusion_order4_vol_3x2v_ser_p1_varnu_diffdirsyz(xc, dx, _cfnu(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_gk_numerical_diffusion_order4_vol_3x2v_ser_p1_varnu_diffdirsxyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct gk_numerical_diffusion* gknd = container_of(eqn, struct gk_numerical_diffusion, eqn);
  return gk_numerical_diffusion_order4_vol_3x2v_ser_p1_varnu_diffdirsxyz(xc, dx, _cfnu(idx), _cfJacInv(idx), qIn, qRhsOut);
}

// Volume kernel list.
GKYL_CU_D
static const gkyl_gk_numerical_diffusion_vol_kern_list ser_vol_kernels_varnu[] = {
  // 1x1v
  {.list={
      // 2nd order diffusion.
      {.list={
          {ker_gk_numerical_diffusion_order2_vol_1x1v_ser_p1_varnu_diffdirsx,NULL,NULL,NULL,NULL,NULL,NULL},
          {NULL,NULL,NULL,NULL,NULL,NULL,NULL},},
      },
      // 4th order diffusion.
      {.list={
          {ker_gk_numerical_diffusion_order4_vol_1x1v_ser_p1_varnu_diffdirsx,NULL,NULL,NULL,NULL,NULL,NULL},
          {NULL,NULL,NULL,NULL,NULL,NULL,NULL},},
      },
    },
  },
  // 1x2v
  {.list={
      // 2nd order diffusion.
      {.list={
          {ker_gk_numerical_diffusion_order2_vol_1x2v_ser_p1_varnu_diffdirsx,NULL,NULL,NULL,NULL,NULL,NULL},
          {NULL,NULL,NULL,NULL,NULL,NULL,NULL},},
      },
      // 4th order diffusion.
      {.list={
          {ker_gk_numerical_diffusion_order4_vol_1x2v_ser_p1_varnu_diffdirsx,NULL,NULL,NULL,NULL,NULL,NULL},
          {NULL,NULL,NULL,NULL,NULL,NULL,NULL},},
      },
    },
  },
  // 2x2v
  {.list={
      // 2nd order diffusion.
      {.list={
          {ker_gk_numerical_diffusion_order2_vol_2x2v_ser_p1_varnu_diffdirsx,ker_gk_numerical_diffusion_order2_vol_2x2v_ser_p1_varnu_diffdirsz,ker_gk_numerical_diffusion_order2_vol_2x2v_ser_p1_varnu_diffdirsxz,NULL,NULL,NULL,NULL},
          {NULL,NULL,NULL,NULL,NULL,NULL,NULL},},
      },
      // 4th order diffusion.
      {.list={
          {ker_gk_numerical_diffusion_order4_vol_2x2v_ser_p1_varnu_diffdirsx,ker_gk_numerical_diffusion_order4_vol_2x2v_ser_p1_varnu_diffdirsz,ker_gk_numerical_diffusion_order4_vol_2x2v_ser_p1_varnu_diffdirsxz,NULL,NULL,NULL,NULL},
          {NULL,NULL,NULL,NULL,NULL,NULL,NULL},},
      },
    },
  },
  // 3x2v
  {.list={
      // 2nd order diffusion.
      {.list={
          {ker_gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffdirsx,ker_gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffdirsy,ker_gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffdirsxy,ker_gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffdirsz,ker_gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffdirsxz,ker_gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffdirsyz,ker_gk_numerical_diffusion_order2_vol_3x2v_ser_p1_varnu_diffdirsxyz},
          {NULL,NULL,NULL,NULL,NULL,NULL,NULL},},
      },
      // 4th order diffusion.
      {.list={
          {ker_gk_numerical_diffusion_order4_vol_3x2v_ser_p1_varnu_diffdirsx,ker_gk_numerical_diffusion_order4_vol_3x2v_ser_p1_varnu_diffdirsy,ker_gk_numerical_diffusion_order4_vol_3x2v_ser_p1_varnu_diffdirsxy,ker_gk_numerical_diffusion_order4_vol_3x2v_ser_p1_varnu_diffdirsz,ker_gk_numerical_diffusion_order4_vol_3x2v_ser_p1_varnu_diffdirsxz,ker_gk_numerical_diffusion_order4_vol_3x2v_ser_p1_varnu_diffdirsyz,ker_gk_numerical_diffusion_order4_vol_3x2v_ser_p1_varnu_diffdirsxyz},
          {NULL,NULL,NULL,NULL,NULL,NULL,NULL},},
      },
    },
  },
};

// Surface kernel list.
GKYL_CU_D static const gkyl_gk_numerical_diffusion_surf_kern_list ser_gyrokinetic_surf_kernels_varnu[] = {
  // Along x.
  {.list={
      // 2nd order diffusion.
      {.list = {
          { gk_numerical_diffusion_order2_surfx_1x1v_ser_p1_varnu, NULL },
          { gk_numerical_diffusion_order2_surfx_1x2v_ser_p1_varnu, NULL },
          { gk_numerical_diffusion_order2_surfx_2x2v_ser_p1_varnu, NULL },
          { gk_numerical_diffusion_order2_surfx_3x2v_ser_p1_varnu, NULL },},
      },
      // 4th order diffusion.
      {.list = {
          { gk_numerical_diffusion_order2_surfx_1x1v_ser_p1_varnu, NULL },
          { gk_numerical_diffusion_order2_surfx_1x2v_ser_p1_varnu, NULL },
          { gk_numerical_diffusion_order2_surfx_2x2v_ser_p1_varnu, NULL },
          { gk_numerical_diffusion_order2_surfx_3x2v_ser_p1_varnu, NULL },},
      },
    },
  },
  // Along y
  {.list={
      // 2nd order diffusion.
      {.list = {
          { NULL, NULL },
          { NULL, NULL },
          { gk_numerical_diffusion_order2_surfz_2x2v_ser_p1_varnu, NULL },
          { gk_numerical_diffusion_order2_surfy_3x2v_ser_p1_varnu, NULL },},
      },
      // 4th order diffusion.
      {.list = {
          { NULL, NULL },
          { NULL, NULL },
          { gk_numerical_diffusion_order2_surfz_2x2v_ser_p1_varnu, NULL },
          { gk_numerical_diffusion_order2_surfy_3x2v_ser_p1_varnu, NULL },},
      },
    },
  },
  // Along z.
  {.list={
      // 2nd order diffusion.
      {.list = {
          { NULL, NULL },
          { NULL, NULL },
          { NULL, NULL },
          { gk_numerical_diffusion_order2_surfz_3x2v_ser_p1_varnu, NULL },},
      },
      // 4th order diffusion.
      {.list = {
          { NULL, NULL },
          { NULL, NULL },
          { NULL, NULL },
          { gk_numerical_diffusion_order2_surfz_3x2v_ser_p1_varnu, NULL },},
      },
    },
  },
};

// Zero-flux boundary surface kernel list.
GKYL_CU_D static const gkyl_gk_numerical_diffusion_boundary_surf_kern_list ser_gyrokinetic_boundary_surf_lower_zeroflux_kernels_varnu[] = {
  // Along x.
  {.list= {
      // 2nd order diffusion.
      {.list= {
          { gk_numerical_diffusion_order2_boundary_surfx_lower_zero_flux_1x1v_ser_p1_varnu, NULL }, // 1x1v
          { gk_numerical_diffusion_order2_boundary_surfx_lower_zero_flux_1x2v_ser_p1_varnu, NULL }, // 1x2v
          { gk_numerical_diffusion_order2_boundary_surfx_lower_zero_flux_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order2_boundary_surfx_lower_zero_flux_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
      // 4th order diffusion.
      {.list= {
          { gk_numerical_diffusion_order4_boundary_surfx_lower_zero_flux_1x1v_ser_p1_varnu, NULL }, // 1x1v
          { gk_numerical_diffusion_order4_boundary_surfx_lower_zero_flux_1x2v_ser_p1_varnu, NULL }, // 1x2v
          { gk_numerical_diffusion_order4_boundary_surfx_lower_zero_flux_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order4_boundary_surfx_lower_zero_flux_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
    },
  },
  // Along y.
  {.list= {
      // 2nd order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { gk_numerical_diffusion_order2_boundary_surfz_lower_zero_flux_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order2_boundary_surfy_lower_zero_flux_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
      // 4th order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { gk_numerical_diffusion_order4_boundary_surfz_lower_zero_flux_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order4_boundary_surfy_lower_zero_flux_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
    },
  },
  // Along z.
  {.list= {
      // 2nd order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { NULL, NULL }, // 2x2v
          { gk_numerical_diffusion_order2_boundary_surfz_lower_zero_flux_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
      // 4th order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { NULL, NULL }, // 2x2v
          { gk_numerical_diffusion_order4_boundary_surfz_lower_zero_flux_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
    },
  },
};

GKYL_CU_D static const gkyl_gk_numerical_diffusion_boundary_surf_kern_list ser_gyrokinetic_boundary_surf_upper_zeroflux_kernels_varnu[] = {
  // Along x.
  {.list= {
      // 2nd order diffusion.
      {.list= {
          { gk_numerical_diffusion_order2_boundary_surfx_upper_zero_flux_1x1v_ser_p1_varnu, NULL }, // 1x1v
          { gk_numerical_diffusion_order2_boundary_surfx_upper_zero_flux_1x2v_ser_p1_varnu, NULL }, // 1x2v
          { gk_numerical_diffusion_order2_boundary_surfx_upper_zero_flux_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order2_boundary_surfx_upper_zero_flux_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
      // 4th order diffusion.
      {.list= {
          { gk_numerical_diffusion_order4_boundary_surfx_upper_zero_flux_1x1v_ser_p1_varnu, NULL }, // 1x1v
          { gk_numerical_diffusion_order4_boundary_surfx_upper_zero_flux_1x2v_ser_p1_varnu, NULL }, // 1x2v
          { gk_numerical_diffusion_order4_boundary_surfx_upper_zero_flux_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order4_boundary_surfx_upper_zero_flux_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
    },
  },
  // Along y.
  {.list= {
      // 2nd order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { gk_numerical_diffusion_order2_boundary_surfz_upper_zero_flux_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order2_boundary_surfy_upper_zero_flux_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
      // 4th order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { gk_numerical_diffusion_order4_boundary_surfz_upper_zero_flux_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order4_boundary_surfy_upper_zero_flux_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
    },
  },
  // Along z.
  {.list= {
      // 2nd order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { NULL, NULL }, // 2x2v
          { gk_numerical_diffusion_order2_boundary_surfz_upper_zero_flux_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
      // 4th order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { NULL, NULL }, // 2x2v
          { gk_numerical_diffusion_order4_boundary_surfz_upper_zero_flux_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
    },
  },
};

// Boundary-local boundary surface kernel list.
GKYL_CU_D static const gkyl_gk_numerical_diffusion_boundary_surf_kern_list ser_gyrokinetic_boundary_surf_lower_boundlocal_kernels_varnu[] = {
  // Along x.
  {.list= {
      // 2nd order diffusion.
      {.list= {
          { gk_numerical_diffusion_order2_boundary_surfx_lower_bound_local_1x1v_ser_p1_varnu, NULL }, // 1x1v
          { gk_numerical_diffusion_order2_boundary_surfx_lower_bound_local_1x2v_ser_p1_varnu, NULL }, // 1x2v
          { gk_numerical_diffusion_order2_boundary_surfx_lower_bound_local_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order2_boundary_surfx_lower_bound_local_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
      // 4th order diffusion.
      {.list= {
          { gk_numerical_diffusion_order4_boundary_surfx_lower_bound_local_1x1v_ser_p1_varnu, NULL }, // 1x1v
          { gk_numerical_diffusion_order4_boundary_surfx_lower_bound_local_1x2v_ser_p1_varnu, NULL }, // 1x2v
          { gk_numerical_diffusion_order4_boundary_surfx_lower_bound_local_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order4_boundary_surfx_lower_bound_local_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
    },
  },
  // Along y.
  {.list= {
      // 2nd order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { gk_numerical_diffusion_order2_boundary_surfz_lower_bound_local_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order2_boundary_surfy_lower_bound_local_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
      // 4th order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { gk_numerical_diffusion_order4_boundary_surfz_lower_bound_local_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order4_boundary_surfy_lower_bound_local_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
    },
  },
  // Along z.
  {.list= {
      // 2nd order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { NULL, NULL }, // 2x2v
          { gk_numerical_diffusion_order2_boundary_surfz_lower_bound_local_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
      // 4th order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { NULL, NULL }, // 2x2v
          { gk_numerical_diffusion_order4_boundary_surfz_lower_bound_local_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
    },
  },
};

GKYL_CU_D static const gkyl_gk_numerical_diffusion_boundary_surf_kern_list ser_gyrokinetic_boundary_surf_upper_boundlocal_kernels_varnu[] = {
  // Along x.
  {.list= {
      // 2nd order diffusion.
      {.list= {
          { gk_numerical_diffusion_order2_boundary_surfx_upper_bound_local_1x1v_ser_p1_varnu, NULL }, // 1x1v
          { gk_numerical_diffusion_order2_boundary_surfx_upper_bound_local_1x2v_ser_p1_varnu, NULL }, // 1x2v
          { gk_numerical_diffusion_order2_boundary_surfx_upper_bound_local_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order2_boundary_surfx_upper_bound_local_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
      // 4th order diffusion.
      {.list= {
          { gk_numerical_diffusion_order4_boundary_surfx_upper_bound_local_1x1v_ser_p1_varnu, NULL }, // 1x1v
          { gk_numerical_diffusion_order4_boundary_surfx_upper_bound_local_1x2v_ser_p1_varnu, NULL }, // 1x2v
          { gk_numerical_diffusion_order4_boundary_surfx_upper_bound_local_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order4_boundary_surfx_upper_bound_local_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
    },
  },
  // Along y.
  {.list= {
      // 2nd order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { gk_numerical_diffusion_order2_boundary_surfz_upper_bound_local_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order2_boundary_surfy_upper_bound_local_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
      // 4th order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { gk_numerical_diffusion_order4_boundary_surfz_upper_bound_local_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order4_boundary_surfy_upper_bound_local_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
    },
  },
  // Along z.
  {.list= {
      // 2nd order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { NULL, NULL }, // 2x2v
          { gk_numerical_diffusion_order2_boundary_surfz_upper_bound_local_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
      // 4th order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { NULL, NULL }, // 2x2v
          { gk_numerical_diffusion_order4_boundary_surfz_upper_bound_local_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
    },
  },
};

// Boundary-local boundary diagnostic kernel list.
GKYL_CU_D static const gkyl_gk_numerical_diffusion_boundary_surf_kern_list ser_gyrokinetic_boundary_diag_lower_boundlocal_kernels_varnu[] = {
  // Along x.
  {.list= {
      // 2nd order diffusion.
      {.list= {
          { gk_numerical_diffusion_order2_boundary_diagx_lower_bound_local_1x1v_ser_p1_varnu, NULL }, // 1x1v
          { gk_numerical_diffusion_order2_boundary_diagx_lower_bound_local_1x2v_ser_p1_varnu, NULL }, // 1x2v
          { gk_numerical_diffusion_order2_boundary_diagx_lower_bound_local_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order2_boundary_diagx_lower_bound_local_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
      // 4th order diffusion.
      {.list= {
          { gk_numerical_diffusion_order4_boundary_diagx_lower_bound_local_1x1v_ser_p1_varnu, NULL }, // 1x1v
          { gk_numerical_diffusion_order4_boundary_diagx_lower_bound_local_1x2v_ser_p1_varnu, NULL }, // 1x2v
          { gk_numerical_diffusion_order4_boundary_diagx_lower_bound_local_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order4_boundary_diagx_lower_bound_local_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
    },
  },
  // Along y.
  {.list= {
      // 2nd order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { gk_numerical_diffusion_order2_boundary_diagz_lower_bound_local_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order2_boundary_diagy_lower_bound_local_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
      // 4th order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { gk_numerical_diffusion_order4_boundary_diagz_lower_bound_local_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order4_boundary_diagy_lower_bound_local_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
    },
  },
  // Along z.
  {.list= {
      // 2nd order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { NULL, NULL }, // 2x2v
          { gk_numerical_diffusion_order2_boundary_diagz_lower_bound_local_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
      // 4th order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { NULL, NULL }, // 2x2v
          { gk_numerical_diffusion_order4_boundary_diagz_lower_bound_local_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
    },
  },
};

GKYL_CU_D static const gkyl_gk_numerical_diffusion_boundary_surf_kern_list ser_gyrokinetic_boundary_diag_upper_boundlocal_kernels_varnu[] = {
  // Along x.
  {.list= {
      // 2nd order diffusion.
      {.list= {
          { gk_numerical_diffusion_order2_boundary_diagx_upper_bound_local_1x1v_ser_p1_varnu, NULL }, // 1x1v
          { gk_numerical_diffusion_order2_boundary_diagx_upper_bound_local_1x2v_ser_p1_varnu, NULL }, // 1x2v
          { gk_numerical_diffusion_order2_boundary_diagx_upper_bound_local_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order2_boundary_diagx_upper_bound_local_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
      // 4th order diffusion.
      {.list= {
          { gk_numerical_diffusion_order4_boundary_diagx_upper_bound_local_1x1v_ser_p1_varnu, NULL }, // 1x1v
          { gk_numerical_diffusion_order4_boundary_diagx_upper_bound_local_1x2v_ser_p1_varnu, NULL }, // 1x2v
          { gk_numerical_diffusion_order4_boundary_diagx_upper_bound_local_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order4_boundary_diagx_upper_bound_local_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
    },
  },
  // Along y.
  {.list= {
      // 2nd order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { gk_numerical_diffusion_order2_boundary_diagz_upper_bound_local_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order2_boundary_diagy_upper_bound_local_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
      // 4th order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { gk_numerical_diffusion_order4_boundary_diagz_upper_bound_local_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order4_boundary_diagy_upper_bound_local_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
    },
  },
  // Along z.
  {.list= {
      // 2nd order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { NULL, NULL }, // 2x2v
          { gk_numerical_diffusion_order2_boundary_diagz_upper_bound_local_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
      // 4th order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { NULL, NULL }, // 2x2v
          { gk_numerical_diffusion_order4_boundary_diagz_upper_bound_local_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
    },
  },
};

// Boundary-recovery boundary diagnostic kernel list.
GKYL_CU_D static const gkyl_gk_numerical_diffusion_boundary_surf_kern_list ser_gyrokinetic_boundary_diag_lower_boundrecovery_kernels_varnu[] = {
  // Along x.
  {.list= {
      // 2nd order diffusion.
      {.list= {
          { gk_numerical_diffusion_order2_boundary_diagx_lower_bound_recovery_1x1v_ser_p1_varnu, NULL }, // 1x1v
          { gk_numerical_diffusion_order2_boundary_diagx_lower_bound_recovery_1x2v_ser_p1_varnu, NULL }, // 1x2v
          { gk_numerical_diffusion_order2_boundary_diagx_lower_bound_recovery_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order2_boundary_diagx_lower_bound_recovery_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
      // 4th order diffusion.
      {.list= {
          { gk_numerical_diffusion_order4_boundary_diagx_lower_bound_recovery_1x1v_ser_p1_varnu, NULL }, // 1x1v
          { gk_numerical_diffusion_order4_boundary_diagx_lower_bound_recovery_1x2v_ser_p1_varnu, NULL }, // 1x2v
          { gk_numerical_diffusion_order4_boundary_diagx_lower_bound_recovery_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order4_boundary_diagx_lower_bound_recovery_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
    },
  },
  // Along y.
  {.list= {
      // 2nd order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { gk_numerical_diffusion_order2_boundary_diagz_lower_bound_recovery_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order2_boundary_diagy_lower_bound_recovery_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
      // 4th order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { gk_numerical_diffusion_order4_boundary_diagz_lower_bound_recovery_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order4_boundary_diagy_lower_bound_recovery_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
    },
  },
  // Along z.
  {.list= {
      // 2nd order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { NULL, NULL }, // 2x2v
          { gk_numerical_diffusion_order2_boundary_diagz_lower_bound_recovery_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
      // 4th order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { NULL, NULL }, // 2x2v
          { gk_numerical_diffusion_order4_boundary_diagz_lower_bound_recovery_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
    },
  },
};

GKYL_CU_D static const gkyl_gk_numerical_diffusion_boundary_surf_kern_list ser_gyrokinetic_boundary_diag_upper_boundrecovery_kernels_varnu[] = {
  // Along x.
  {.list= {
      // 2nd order diffusion.
      {.list= {
          { gk_numerical_diffusion_order2_boundary_diagx_upper_bound_recovery_1x1v_ser_p1_varnu, NULL }, // 1x1v
          { gk_numerical_diffusion_order2_boundary_diagx_upper_bound_recovery_1x2v_ser_p1_varnu, NULL }, // 1x2v
          { gk_numerical_diffusion_order2_boundary_diagx_upper_bound_recovery_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order2_boundary_diagx_upper_bound_recovery_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
      // 4th order diffusion.
      {.list= {
          { gk_numerical_diffusion_order4_boundary_diagx_upper_bound_recovery_1x1v_ser_p1_varnu, NULL }, // 1x1v
          { gk_numerical_diffusion_order4_boundary_diagx_upper_bound_recovery_1x2v_ser_p1_varnu, NULL }, // 1x2v
          { gk_numerical_diffusion_order4_boundary_diagx_upper_bound_recovery_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order4_boundary_diagx_upper_bound_recovery_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
    },
  },
  // Along y.
  {.list= {
      // 2nd order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { gk_numerical_diffusion_order2_boundary_diagz_upper_bound_recovery_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order2_boundary_diagy_upper_bound_recovery_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
      // 4th order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { gk_numerical_diffusion_order4_boundary_diagz_upper_bound_recovery_2x2v_ser_p1_varnu, NULL }, // 2x2v
          { gk_numerical_diffusion_order4_boundary_diagy_upper_bound_recovery_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
    },
  },
  // Along z.
  {.list= {
      // 2nd order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { NULL, NULL }, // 2x2v
          { gk_numerical_diffusion_order2_boundary_diagz_upper_bound_recovery_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
      // 4th order diffusion.
      {.list= {
          { NULL, NULL }, // 1x1v
          { NULL, NULL }, // 1x2v
          { NULL, NULL }, // 2x2v
          { gk_numerical_diffusion_order4_boundary_diagz_upper_bound_recovery_3x2v_ser_p1_varnu, NULL }, // 3x2v
        },
      },
    },
  },
};

// Macro for choosing volume and surface kernels.
#define CKVOL(lst,cdim,vdim,diff_order,poly_order,diffdir_linidx) lst[cdim+vdim-2].list[diff_order/2-1].list[poly_order-1].kernels[diffdir_linidx]
#define CKSURF(lst,cdim,vdim,diff_order,poly_order,dir) lst[dir].list[diff_order/2-1].list[cdim+vdim-2].kernels[poly_order-1]
#define CKBSURF(lst,cdim,vdim,diff_order,poly_order) lst->list[diff_order/2-1].list[cdim+vdim-2].kernels[poly_order-1]

GKYL_CU_D static double surf(const struct gkyl_dg_eqn* eqn, int dir,
  const double* xcL, const double* xcC, const double* xcR, 
  const double* dxL, const double* dxC, const double* dxR,
  const int* idxL, const int* idxC, const int* idxR,
  const double* qInL, const double* qInC, const double* qInR,
  double* GKYL_RESTRICT qRhsOut)
{
  struct gk_numerical_diffusion* gknd = container_of(eqn, struct gk_numerical_diffusion, eqn);

  if (gknd->has_diff_in_dir[dir]) {
    gknd->surf[dir](xcC, dxC, _cfnu(idxL), _cfnu(idxC), _cfnu(idxR), _cfJacInv(idxL), _cfJacInv(idxC), _cfJacInv(idxR), qInL, qInC, qInR, qRhsOut);
  }
  return 0.;  // CFL frequency computed in volume term.
}

GKYL_CU_D static double boundary_surf(const struct gkyl_dg_eqn* eqn, int dir,
  const double* xcEdge, const double* xcSkin, const double* dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut)
{ 
  struct gk_numerical_diffusion* gknd = container_of(eqn, struct gk_numerical_diffusion, eqn);
  
  if (gknd->has_diff_in_dir[dir]) {
    if (edge == -1) {
      gknd->boundary_surf_lower[dir](xcSkin, dxSkin, _cfnu(idxEdge), _cfnu(idxSkin), _cfJacInv(idxEdge), _cfJacInv(idxSkin), edge, qInEdge, qInSkin, qRhsOut);
    }
    else {
      gknd->boundary_surf_upper[dir](xcSkin, dxSkin, _cfnu(idxEdge), _cfnu(idxSkin), _cfJacInv(idxEdge), _cfJacInv(idxSkin), edge, qInEdge, qInSkin, qRhsOut);
    }
  }
  return 0.;  // CFL frequency computed in volume term.
}

GKYL_CU_D static double boundary_diag(const struct gkyl_dg_eqn* eqn, int dir,
  const double* xcSkin, const double* xcGhost, const double* dxSkin, const double* dxGhost,
  const int* idxSkin, const int* idxGhost, const int edge,
  const double* qInSkin, const double* qInGhost, double* GKYL_RESTRICT qRhsGhost)
{ 
  // This function is based on boundary_surf above, but notice we use Edge
  // where the boundary_surf used Skin, because we assume this kernel is called
  // in the ghost range (e.g. by the boundary_flux updater).
  struct gk_numerical_diffusion* gknd = container_of(eqn, struct gk_numerical_diffusion, eqn);
  
  if (gknd->has_diff_in_dir[dir]) {
    if (edge == -1) {
      gknd->boundary_diag_lower[dir](xcSkin, dxSkin, _cfnu(idxSkin), _cfnu(idxGhost), _cfJacInv(idxSkin), _cfJacInv(idxGhost), edge, qInSkin, qInGhost, qRhsGhost);
    }
    else {
      gknd->boundary_diag_upper[dir](xcSkin, dxSkin, _cfnu(idxSkin), _cfnu(idxGhost), _cfJacInv(idxSkin), _cfJacInv(idxGhost), edge, qInSkin, qInGhost, qRhsGhost);
    }
  }
  return 0.;  // CFL frequency computed in volume term.
}

#undef _cfnu
#undef _cfJacInv

/**
 * Free diffusion equation object
 *
 * @param ref Reference counter for diffusion equation.
 */
void gkyl_gk_numerical_diffusion_free(const struct gkyl_ref_count* ref);

#ifdef GKYL_HAVE_CUDA
/**
 * Create a new gyrokinetic diffusion equation object on the device.
 *
 * @param basis Basis functions of the equation system.
 * @param cbasis Configuration space basis.
 * @param conf_range Conf-space range object.
 * @param bc_lower Boundary condition at lower boundaries.
 * @param bc_upper Boundary condition at upper boundaries.
 * @param has_diff_in_dir Whether to apply diffusion in each direction.
 * @param diff_order Diffusion order.
 * @return Pointer to diffusion equation object
 */
struct gkyl_dg_eqn*
gkyl_gk_numerical_diffusion_cu_dev_new(const struct gkyl_basis *basis, const struct gkyl_basis *cbasis,
  const struct gkyl_range *conf_range, const enum gkyl_gyrokinetic_bc_type *bc_lower, const enum gkyl_gyrokinetic_bc_type *bc_upper,
  const bool *has_diff_in_dir, int diff_order);

/**
 * CUDA device function to set auxiliary fields (e.g. diffusion tensor D) needed in updating diffusion equation.
 *
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_gk_numerical_diffusion_set_auxfields_cu(const struct gkyl_dg_eqn *eqn,
  struct gkyl_gk_numerical_diffusion_auxfields auxin);
#endif
