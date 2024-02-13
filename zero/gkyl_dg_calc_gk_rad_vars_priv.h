// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_rad_gyrokinetic_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

typedef void (*rad_gyrokinetic_nu_t)(const double *w, const double *dxv, 
  double charge, double mass, double a, double alpha, double beta, double gamma, double v0, 
  const double *bmag, double* GKYL_RESTRICT drag_rad_surf, double* GKYL_RESTRICT drag_rad); 

typedef void (*rad_gyrokinetic_nI_nu_t)(const double *vnu_surf, const double *vnu,
  const double *vsqnu_surf, const double *vsqnu, const double *nI, 
  double* GKYL_RESTRICT nvnu_surf, double* GKYL_RESTRICT nvnu, 
  double* GKYL_RESTRICT nvsqnu_surf, double* GKYL_RESTRICT nvsqnu); 

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below.
GKYL_CU_D
static struct { int vdim[3]; } cv_index[] = {
  {-1, -1, -1}, // 0x makes no sense.
  {-1,  0,  1}, // 1x kernel indices.
  {-1, -1,  2}, // 2x kernel indices.
  {-1, -1,  3}, // 3x kernel indices.
};

// for use in kernel tables
typedef struct { rad_gyrokinetic_nu_t kernels[3]; } gkyl_dg_rad_gyrokinetic_nu_kern_list;
typedef struct { rad_gyrokinetic_nI_nu_t kernels[3]; } gkyl_dg_rad_gyrokinetic_nI_nu_kern_list;

struct gkyl_dg_calc_gk_rad_vars {
  struct gkyl_rect_grid phase_grid; // Phase space grid for cell spacing and cell center
  int cdim; // Configuration space dimensionality
  int pdim; // Phase space dimensionality
  rad_gyrokinetic_nu_t rad_nu_vpar; // kernel for computing surface and volume expansions of vpar radiation drag
  rad_gyrokinetic_nu_t rad_nu_mu; // kernel for computing surface and volume expansions of mu radiation drag
  rad_gyrokinetic_nI_nu_t rad_nI_nu; // kernel for density weighted radiation drag

  double charge, mass;
  double a, alpha, beta, gamma, v0; // Fitting parameters for radiation drag coefficient
  const struct gk_geometry *gk_geom; // Pointer to geometry struct

  uint32_t flags;
  struct gkyl_dg_calc_gk_rad_vars *on_dev; // pointer to itself or device data
};

// Radiation drag surface and volume expansions in vpar (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_rad_gyrokinetic_nu_kern_list ser_rad_gyrokinetic_nu_vpar_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, rad_gyrokinetic_drag_nuvpar_1x2v_ser_p1, rad_gyrokinetic_drag_nuvpar_1x2v_ser_p2 }, // 1
  { NULL, rad_gyrokinetic_drag_nuvpar_2x2v_ser_p1, rad_gyrokinetic_drag_nuvpar_2x2v_ser_p2 }, // 2
  { NULL, rad_gyrokinetic_drag_nuvpar_3x2v_ser_p1, NULL }, // 3
};

// Radiation drag surface and volume expansions in mu (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_rad_gyrokinetic_nu_kern_list ser_rad_gyrokinetic_nu_mu_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, rad_gyrokinetic_drag_numu_1x2v_ser_p1, rad_gyrokinetic_drag_numu_1x2v_ser_p2 }, // 1
  { NULL, rad_gyrokinetic_drag_numu_2x2v_ser_p1, rad_gyrokinetic_drag_numu_2x2v_ser_p2 }, // 2
  { NULL, rad_gyrokinetic_drag_numu_3x2v_ser_p1, NULL }, // 3
};

// Density-weighted radiation drag surface and volume expansions (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_rad_gyrokinetic_nI_nu_kern_list ser_rad_gyrokinetic_nI_nu_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, rad_gyrokinetic_drag_nI_nu_1x2v_ser_p1, rad_gyrokinetic_drag_nI_nu_1x2v_ser_p2 }, // 1
  { NULL, rad_gyrokinetic_drag_nI_nu_2x2v_ser_p1, rad_gyrokinetic_drag_nI_nu_2x2v_ser_p2 }, // 2
  { NULL, rad_gyrokinetic_drag_nI_nu_3x2v_ser_p1, NULL }, // 3
};

GKYL_CU_D
static rad_gyrokinetic_nu_t
choose_rad_gyrokinetic_nu_vpar_kern(int cdim, int vdim, int poly_order)
{
  return ser_rad_gyrokinetic_nu_vpar_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
}

GKYL_CU_D
static rad_gyrokinetic_nu_t
choose_rad_gyrokinetic_nu_mu_kern(int cdim, int vdim, int poly_order)
{
  return ser_rad_gyrokinetic_nu_mu_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
}

GKYL_CU_D
static rad_gyrokinetic_nI_nu_t
choose_rad_gyrokinetic_nI_nu_kern(int cdim, int vdim, int poly_order)
{
  return ser_rad_gyrokinetic_nI_nu_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
}
