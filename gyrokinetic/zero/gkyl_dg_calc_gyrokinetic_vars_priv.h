// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

typedef double (*gyrokinetic_flux_surf_t)(const double *w, const double *dxv, 
  const double *vmap, const double *vmapSq, const double q_, const double m_, 
  const struct gkyl_dg_surf_geom *dgs, const struct gkyl_gk_dg_surf_geom *gkdgs, 
  const double *bmag, const double *phi, const double *JfL, const double *JfR, double* GKYL_RESTRICT flux_surf); 

typedef double (*gyrokinetic_flux_surfvpar_t)( 
  const double *w, const double *dxv, 
  const double *vmap_prime_l, const double *vmap_prime_r,
  const double *vmap, const double *vmapSq, const double q_, const double m_, 
  const struct gkyl_dg_vol_geom *dgv, const struct gkyl_gk_dg_vol_geom *gkdgv, 
  const double *bmag, const double *phi, const double *JfL, const double *JfR, double* GKYL_RESTRICT flux_surf); 

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
typedef struct { gyrokinetic_flux_surf_t kernels[3]; } gkyl_dg_gyrokinetic_flux_surf_kern_list;
typedef struct { gyrokinetic_flux_surfvpar_t kernels[3]; } gkyl_dg_gyrokinetic_flux_surfvpar_kern_list;

struct gkyl_dg_calc_gyrokinetic_vars {
  struct gkyl_rect_grid phase_grid; // Phase space grid for cell spacing and cell center
  int cdim; // Configuration space dimensionality
  int pdim; // Phase space dimensionality
  gyrokinetic_flux_surf_t flux_surf[3]; // kernel for computing surface expansion of phase space flux alpha
  gyrokinetic_flux_surfvpar_t flux_surfvpar[1]; // kernel for computing surface expansion of phase space flux alpha
  gyrokinetic_flux_surf_t flux_edge_surf[3]; // kernel for computing surface expansion of phase space flux alpha
                                               // at upper configuration space edge
  double charge, mass;
  const struct gk_geometry *gk_geom; // Pointer to geometry struct.
  const struct gkyl_dg_geom *dg_geom; // Pointer to vol dg geometry struct.
  const struct gkyl_gk_dg_geom *gk_dg_geom; // Pointer to vol gk dg geometry struct.
  const struct gkyl_velocity_map *vel_map; // Velocity space mapping object.

  uint32_t flags;
  struct gkyl_dg_calc_gyrokinetic_vars *on_dev; // pointer to itself or device data
};

//
// Serendipity surface kernels general geometry
//
// Gyrokinetic phase space flux alpha surface expansions in x (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_gyrokinetic_flux_surf_kern_list ser_gyrokinetic_flux_surfx_kernels[] = {
  { NULL, gyrokinetic_flux_surfx_1x1v_ser_p1, NULL }, // 0
  { NULL, gyrokinetic_flux_surfx_1x2v_ser_p1, NULL }, // 1
  { NULL, gyrokinetic_flux_surfx_2x2v_ser_p1, NULL }, // 2
  { NULL, gyrokinetic_flux_surfx_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux alpha edge surface expansions in x (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_gyrokinetic_flux_surf_kern_list ser_gyrokinetic_flux_edge_surfx_kernels[] = {
  { NULL, gyrokinetic_flux_edge_surfx_1x1v_ser_p1, NULL }, // 0
  { NULL, gyrokinetic_flux_edge_surfx_1x2v_ser_p1, NULL }, // 1
  { NULL, gyrokinetic_flux_edge_surfx_2x2v_ser_p1, NULL }, // 2
  { NULL, gyrokinetic_flux_edge_surfx_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux flux surface expansions in y (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_gyrokinetic_flux_surf_kern_list ser_gyrokinetic_flux_surfy_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gyrokinetic_flux_surfy_2x2v_ser_p1, NULL }, // 2
  { NULL, gyrokinetic_flux_surfy_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux flux edge surface expansions in y (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_gyrokinetic_flux_surf_kern_list ser_gyrokinetic_flux_edge_surfy_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gyrokinetic_flux_edge_surfy_2x2v_ser_p1, NULL }, // 2
  { NULL, gyrokinetic_flux_edge_surfy_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux flux surface expansions in z (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_gyrokinetic_flux_surf_kern_list ser_gyrokinetic_flux_surfz_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gyrokinetic_flux_surfz_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux flux edge surface expansions in z (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_gyrokinetic_flux_surf_kern_list ser_gyrokinetic_flux_edge_surfz_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gyrokinetic_flux_edge_surfz_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux alpha surface expansions in vpar (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_gyrokinetic_flux_surfvpar_kern_list ser_gyrokinetic_flux_surfvpar_kernels[] = {
  { NULL, gyrokinetic_flux_surfvpar_1x1v_ser_p1, NULL }, // 0
  { NULL, gyrokinetic_flux_surfvpar_1x2v_ser_p1, NULL }, // 1
  { NULL, gyrokinetic_flux_surfvpar_2x2v_ser_p1, NULL }, // 2
  { NULL, gyrokinetic_flux_surfvpar_3x2v_ser_p1, NULL }, // 3
};

//
// Serendipity surface kernels general geometry, no toroidal field (by=0)
//
// Gyrokinetic phase space flux alpha surface expansions in x (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_gyrokinetic_flux_surf_kern_list ser_gyrokinetic_flux_no_by_surfx_kernels[] = {
  { NULL, gyrokinetic_flux_surfx_1x1v_ser_p1, NULL }, // 0
  { NULL, gyrokinetic_flux_surfx_1x2v_ser_p1, NULL }, // 1
  { NULL, gyrokinetic_flux_no_by_surfx_2x2v_ser_p1, NULL }, // 2
  { NULL, gyrokinetic_flux_no_by_surfx_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux alpha edge surface expansions in x (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_gyrokinetic_flux_surf_kern_list ser_gyrokinetic_flux_no_by_edge_surfx_kernels[] = {
  { NULL, gyrokinetic_flux_edge_surfx_1x1v_ser_p1, NULL }, // 0
  { NULL, gyrokinetic_flux_edge_surfx_1x2v_ser_p1, NULL }, // 1
  { NULL, gyrokinetic_flux_no_by_edge_surfx_2x2v_ser_p1, NULL }, // 2
  { NULL, gyrokinetic_flux_no_by_edge_surfx_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux alpha surface expansions in y (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_gyrokinetic_flux_surf_kern_list ser_gyrokinetic_flux_no_by_surfy_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gyrokinetic_flux_no_by_surfy_2x2v_ser_p1, NULL }, // 2
  { NULL, gyrokinetic_flux_no_by_surfy_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux alpha edge surface expansions in y (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_gyrokinetic_flux_surf_kern_list ser_gyrokinetic_flux_no_by_edge_surfy_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gyrokinetic_flux_no_by_edge_surfy_2x2v_ser_p1, NULL }, // 2
  { NULL, gyrokinetic_flux_no_by_edge_surfy_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux alpha surface expansions in z (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_gyrokinetic_flux_surf_kern_list ser_gyrokinetic_flux_no_by_surfz_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gyrokinetic_flux_no_by_surfz_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux alpha edge surface expansions in z (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_gyrokinetic_flux_surf_kern_list ser_gyrokinetic_flux_no_by_edge_surfz_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gyrokinetic_flux_no_by_edge_surfz_3x2v_ser_p1, NULL }, // 3
};

// Gyrokinetic phase space flux alpha surface expansions in vpar (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_gyrokinetic_flux_surfvpar_kern_list ser_gyrokinetic_flux_no_by_surfvpar_kernels[] = {
  { NULL, gyrokinetic_flux_surfvpar_1x1v_ser_p1, NULL }, // 0
  { NULL, gyrokinetic_flux_surfvpar_1x2v_ser_p1, NULL }, // 1
  { NULL, gyrokinetic_flux_no_by_surfvpar_2x2v_ser_p1, NULL }, // 2
  { NULL, gyrokinetic_flux_no_by_surfvpar_3x2v_ser_p1, NULL }, // 3
};

GKYL_CU_D
static gyrokinetic_flux_surf_t
choose_gyrokinetic_flux_surf_conf_kern(int dir, int cdim, int vdim, int poly_order)
{
  if (dir == 0)
    return ser_gyrokinetic_flux_surfx_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  else if (dir == 1)
    return ser_gyrokinetic_flux_surfy_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  else if (dir == 2)
    return ser_gyrokinetic_flux_surfz_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  else
    return NULL;
}

GKYL_CU_D
static gyrokinetic_flux_surf_t
choose_gyrokinetic_flux_edge_surf_conf_kern(int dir, int cdim, int vdim, int poly_order)
{
  if (dir == 0)
    return ser_gyrokinetic_flux_edge_surfx_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  else if (dir == 1)
    return ser_gyrokinetic_flux_edge_surfy_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  else if (dir == 2)
    return ser_gyrokinetic_flux_edge_surfz_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  else
    return NULL;
}

GKYL_CU_D
static gyrokinetic_flux_surfvpar_t
choose_gyrokinetic_flux_surf_vpar_kern(int cdim, int vdim, int poly_order)
{
  return ser_gyrokinetic_flux_surfvpar_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
}

GKYL_CU_D
static gyrokinetic_flux_surf_t
choose_gyrokinetic_flux_no_by_surf_conf_kern(int dir, int cdim, int vdim, int poly_order)
{
  if (dir == 0)
    return ser_gyrokinetic_flux_no_by_surfx_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  else if (dir == 1)
    return ser_gyrokinetic_flux_no_by_surfy_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  else if (dir == 2)
    return ser_gyrokinetic_flux_no_by_surfz_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  else
    return NULL;
}

GKYL_CU_D
static gyrokinetic_flux_surf_t
choose_gyrokinetic_flux_no_by_edge_surf_conf_kern(int dir, int cdim, int vdim, int poly_order)
{
  if (dir == 0)
    return ser_gyrokinetic_flux_no_by_edge_surfx_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  else if (dir == 1)
    return ser_gyrokinetic_flux_no_by_edge_surfy_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  else if (dir == 2)
    return ser_gyrokinetic_flux_no_by_edge_surfz_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  else
    return NULL;
}

GKYL_CU_D
static gyrokinetic_flux_surfvpar_t
choose_gyrokinetic_flux_no_by_surf_vpar_kern(int cdim, int vdim, int poly_order)
{
  return ser_gyrokinetic_flux_no_by_surfvpar_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
}

#ifdef GKYL_HAVE_CUDA
/**
 * Create new updater to compute gyrokinetic variables on
 * NV-GPU. See new() method for documentation.
 */
struct gkyl_dg_calc_gyrokinetic_vars* 
gkyl_dg_calc_gyrokinetic_vars_cu_dev_new(const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  double charge, double mass, enum gkyl_gkmodel_id gkmodel_id, 
  const struct gk_geometry *gk_geom, const struct gkyl_dg_geom *dg_geom, 
  const struct gkyl_gk_dg_geom *gk_dg_geom, const struct gkyl_velocity_map *vel_map);

/**
 * Host-side wrappers for gyrokinetic vars operations on device
 */
void gkyl_dg_calc_gyrokinetic_vars_flux_surf_cu(struct gkyl_dg_calc_gyrokinetic_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_range *conf_ext_range, const struct gkyl_range *phase_ext_range, const struct gkyl_array *phi, 
  const struct gkyl_array* fin, struct gkyl_array* flux_surf, struct gkyl_array* cflrate);
#endif
