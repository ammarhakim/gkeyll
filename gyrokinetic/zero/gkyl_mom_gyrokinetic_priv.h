#pragma once

// Private header, not for direct use in user code

#include <gkyl_util.h>
#include <gkyl_array.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_velocity_map.h>
#include <gkyl_mom_type.h>
#include <gkyl_range.h>
#include <gkyl_ref_count.h>
#include <gkyl_mom_gyrokinetic_kernels.h>


struct mom_type_gyrokinetic {
  struct gkyl_mom_type momt;
  double mass; // Mass of species.
  double charge; // Charge of species.
  struct gkyl_range conf_range; // Configuration space range.
  const struct gk_geometry *gk_geom; // Pointer to geometry struct.
  const struct gkyl_velocity_map *vel_map; // Velocity space mapping.
  struct gkyl_array *phi; // Electrostatic potential.
};

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
static struct { int vdim[3]; } cv_index[] = {
  {-1, -1, -1}, // 0x makes no sense
  {-1,  0,  1}, // 1x kernel indices
  {-1, -1,  2}, // 2x kernel indices
  {-1, -1,  3}, // 3x kernel indices  
};

// for use in kernel tables
typedef struct {
  momf_t kernels[3];
} gkyl_gyrokinetic_mom_kern_list;

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0_1x1v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0_1x1v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0_1x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0_1x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0_2x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0_2x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0_3x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M1_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M1_1x1v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M1_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M1_1x1v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M1_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M1_1x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M1_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M1_1x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M1_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M1_2x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M1_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M1_2x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M1_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M1_3x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M2_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M2_1x1v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M2_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M2_1x1v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M2_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M2_1x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M2_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M2_1x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M2_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M2_2x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M2_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M2_2x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M2_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M2_3x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M3_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M3_1x1v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M3_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M3_1x1v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M3_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M3_1x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M3_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M3_1x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M3_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M3_2x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M3_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M3_2x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M3_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M3_3x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M2PAR_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M2PAR_1x1v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M2PAR_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M2PAR_1x1v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M2PAR_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M2PAR_1x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M2PAR_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M2PAR_1x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M2PAR_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M2PAR_2x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M2PAR_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M2PAR_2x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M2PAR_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M2PAR_3x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M2PERP_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M2PERP_1x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M2PERP_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M2PERP_1x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M2PERP_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M2PERP_2x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M2PERP_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M2PERP_2x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M2PERP_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M2PERP_3x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M3PAR_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M3PAR_1x1v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M3PAR_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M3PAR_1x1v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M3PAR_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M3PAR_1x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M3PAR_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M3PAR_1x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M3PAR_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M3PAR_2x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M3PAR_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M3PAR_2x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M3PAR_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M3PAR_3x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M3PERP_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M3PERP_1x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M3PERP_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M3PERP_1x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M3PERP_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M3PERP_2x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M3PERP_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M3PERP_2x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M3PERP_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M3PERP_3x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0M1M2_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0M1M2_1x1v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0M1M2_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0M1M2_1x1v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0M1M2_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0M1M2_1x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0M1M2_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0M1M2_1x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0M1M2_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0M1M2_2x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0M1M2_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0M1M2_2x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0M1M2_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0M1M2_3x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0M1M2M3_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0M1M2M3_1x1v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0M1M2M3_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0M1M2M3_1x1v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0M1M2M3_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0M1M2M3_1x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0M1M2M3_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0M1M2M3_1x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0M1M2M3_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0M1M2M3_2x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0M1M2M3_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0M1M2M3_2x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0M1M2M3_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0M1M2M3_3x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0M1M2PARM2PERP_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0M1M2PARM2PERP_1x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0M1M2PARM2PERP_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0M1M2PARM2PERP_1x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0M1M2PARM2PERP_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0M1M2PARM2PERP_2x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0M1M2PARM2PERP_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0M1M2PARM2PERP_2x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0M1M2PARM2PERP_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0M1M2PARM2PERP_3x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_1x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_1x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_2x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_2x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_3x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), f, out);
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_HAMILTONIAN_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_HAMILTONIAN_1x1v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass, mom_gk->charge,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    (const double*) gkyl_array_cfetch(mom_gk->phi, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_HAMILTONIAN_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_HAMILTONIAN_1x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass, mom_gk->charge,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    (const double*) gkyl_array_cfetch(mom_gk->phi, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_HAMILTONIAN_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_HAMILTONIAN_2x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass, mom_gk->charge,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    (const double*) gkyl_array_cfetch(mom_gk->phi, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_HAMILTONIAN_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_HAMILTONIAN_3x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass, mom_gk->charge,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    (const double*) gkyl_array_cfetch(mom_gk->phi, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_HAMILTONIAN_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_HAMILTONIAN_1x1v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass, mom_gk->charge,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    (const double*) gkyl_array_cfetch(mom_gk->phi, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_HAMILTONIAN_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_HAMILTONIAN_1x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass, mom_gk->charge,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    (const double*) gkyl_array_cfetch(mom_gk->phi, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_mom_gyrokinetic_HAMILTONIAN_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return mom_gyrokinetic_HAMILTONIAN_2x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass, mom_gk->charge,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    (const double*) gkyl_array_cfetch(mom_gk->phi, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0_1x1v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0_1x1v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0_1x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0_1x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0_2x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);   
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0_2x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0_3x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M1_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M1_1x1v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M1_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M1_1x1v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M1_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M1_1x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M1_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M1_1x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M1_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M1_2x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);   
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M1_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M1_2x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M1_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M1_3x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M2_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M2_1x1v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M2_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M2_1x1v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M2_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M2_1x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M2_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M2_1x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M2_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M2_2x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);   
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M2_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M2_2x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M2_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M2_3x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M3_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M3_1x1v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M3_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M3_1x1v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M3_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M3_1x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M3_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M3_1x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M3_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M3_2x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);   
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M3_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M3_2x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M3_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M3_3x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M2PAR_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M2PAR_1x1v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M2PAR_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M2PAR_1x1v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M2PAR_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M2PAR_1x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M2PAR_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M2PAR_1x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M2PAR_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M2PAR_2x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);   
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M2PAR_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M2PAR_2x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M2PAR_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M2PAR_3x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M2PERP_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M2PERP_1x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M2PERP_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M2PERP_1x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M2PERP_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M2PERP_2x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);   
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M2PERP_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M2PERP_2x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M2PERP_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M2PERP_3x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M3PAR_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M3PAR_1x1v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M3PAR_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M3PAR_1x1v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M3PAR_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M3PAR_1x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M3PAR_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M3PAR_1x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M3PAR_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M3PAR_2x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);   
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M3PAR_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M3PAR_2x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M3PAR_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M3PAR_3x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M3PERP_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M3PERP_1x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M3PERP_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M3PERP_1x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M3PERP_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M3PERP_2x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);   
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M3PERP_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M3PERP_2x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M3PERP_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M3PERP_3x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0M1M2_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0M1M2_1x1v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0M1M2_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0M1M2_1x1v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0M1M2_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0M1M2_1x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0M1M2_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0M1M2_1x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0M1M2_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0M1M2_2x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);   
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0M1M2_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0M1M2_2x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0M1M2_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0M1M2_3x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0M1M2M3_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0M1M2M3_1x1v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0M1M2M3_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0M1M2M3_1x1v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0M1M2M3_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0M1M2M3_1x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0M1M2M3_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0M1M2M3_1x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0M1M2M3_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0M1M2M3_2x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);   
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0M1M2M3_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0M1M2M3_2x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0M1M2M3_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0M1M2M3_3x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0M1M2PARM2PERP_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0M1M2PARM2PERP_1x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0M1M2PARM2PERP_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0M1M2PARM2PERP_1x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0M1M2PARM2PERP_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0M1M2PARM2PERP_2x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);   
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0M1M2PARM2PERP_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0M1M2PARM2PERP_2x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0M1M2PARM2PERP_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0M1M2PARM2PERP_3x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_1x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_1x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_2x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);   
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_2x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_3x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_HAMILTONIAN_1x1v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_HAMILTONIAN_1x1v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass, mom_gk->charge,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    (const double*) gkyl_array_cfetch(mom_gk->phi, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_HAMILTONIAN_1x1v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_HAMILTONIAN_1x1v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass, mom_gk->charge,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    (const double*) gkyl_array_cfetch(mom_gk->phi, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_HAMILTONIAN_1x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_HAMILTONIAN_1x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass, mom_gk->charge,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    (const double*) gkyl_array_cfetch(mom_gk->phi, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_HAMILTONIAN_1x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_HAMILTONIAN_1x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass, mom_gk->charge,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    (const double*) gkyl_array_cfetch(mom_gk->phi, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_HAMILTONIAN_2x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_HAMILTONIAN_2x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass, mom_gk->charge,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    (const double*) gkyl_array_cfetch(mom_gk->phi, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_HAMILTONIAN_2x2v_ser_p2(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_HAMILTONIAN_2x2v_ser_p2(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass, mom_gk->charge,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    (const double*) gkyl_array_cfetch(mom_gk->phi, cidx), f, out);  
}

GKYL_CU_DH
static void
kernel_int_mom_gyrokinetic_HAMILTONIAN_3x2v_ser_p1(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_gyrokinetic *mom_gk = container_of(momt, struct mom_type_gyrokinetic, momt);

  int idx_vel[2];
  for (int d=momt->cdim; d<momt->pdim; d++) idx_vel[d-momt->cdim] = idx[d];

  long cidx = gkyl_range_idx(&mom_gk->conf_range, idx);
  long vidx = gkyl_range_idx(&mom_gk->vel_map->local_vel, idx_vel);
  return int_mom_gyrokinetic_HAMILTONIAN_3x2v_ser_p1(dx,
    (const double*) gkyl_array_cfetch(mom_gk->vel_map->vmap, vidx), mom_gk->mass, mom_gk->charge,
    (const double*) gkyl_array_cfetch(mom_gk->gk_geom->geo_corn.bmag, cidx), 
    (const double*) gkyl_array_cfetch(mom_gk->phi, cidx), f, out);  
}

//
// Serendipity basis kernels
//

// M0 kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m0_kernels[] = {
  // 1x kernels
  { NULL, kernel_mom_gyrokinetic_M0_1x1v_ser_p1, kernel_mom_gyrokinetic_M0_1x1v_ser_p2 }, // 0
  { NULL, kernel_mom_gyrokinetic_M0_1x2v_ser_p1, kernel_mom_gyrokinetic_M0_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_mom_gyrokinetic_M0_2x2v_ser_p1, kernel_mom_gyrokinetic_M0_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_mom_gyrokinetic_M0_3x2v_ser_p1, NULL }, // 3
};

// M1 kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m1_kernels[] = {
  // 1x kernels
  { NULL, kernel_mom_gyrokinetic_M1_1x1v_ser_p1, kernel_mom_gyrokinetic_M1_1x1v_ser_p2 }, // 0
  { NULL, kernel_mom_gyrokinetic_M1_1x2v_ser_p1, kernel_mom_gyrokinetic_M1_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_mom_gyrokinetic_M1_2x2v_ser_p1, kernel_mom_gyrokinetic_M1_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_mom_gyrokinetic_M1_3x2v_ser_p1, NULL }, // 3
};

// M2 kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m2_kernels[] = {
  // 1x kernels
  { NULL, kernel_mom_gyrokinetic_M2_1x1v_ser_p1, kernel_mom_gyrokinetic_M2_1x1v_ser_p2 }, // 0
  { NULL, kernel_mom_gyrokinetic_M2_1x2v_ser_p1, kernel_mom_gyrokinetic_M2_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_mom_gyrokinetic_M2_2x2v_ser_p1, kernel_mom_gyrokinetic_M2_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_mom_gyrokinetic_M2_3x2v_ser_p1, NULL }, // 3
};

// M3 kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m3_kernels[] = {
  // 1x kernels
  { NULL, kernel_mom_gyrokinetic_M3_1x1v_ser_p1, kernel_mom_gyrokinetic_M3_1x1v_ser_p2 }, // 0
  { NULL, kernel_mom_gyrokinetic_M3_1x2v_ser_p1, kernel_mom_gyrokinetic_M3_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_mom_gyrokinetic_M3_2x2v_ser_p1, kernel_mom_gyrokinetic_M3_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_mom_gyrokinetic_M3_3x2v_ser_p1, NULL }, // 3
};

// M2 parallel kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m2par_kernels[] = {
  // 1x kernels
  { NULL, kernel_mom_gyrokinetic_M2PAR_1x1v_ser_p1, kernel_mom_gyrokinetic_M2PAR_1x1v_ser_p2 }, // 0
  { NULL, kernel_mom_gyrokinetic_M2PAR_1x2v_ser_p1, kernel_mom_gyrokinetic_M2PAR_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_mom_gyrokinetic_M2PAR_2x2v_ser_p1, kernel_mom_gyrokinetic_M2PAR_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_mom_gyrokinetic_M2PAR_3x2v_ser_p1, NULL }, // 3
};

// M2 perpendicular kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m2perp_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, kernel_mom_gyrokinetic_M2PERP_1x2v_ser_p1, kernel_mom_gyrokinetic_M2PERP_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_mom_gyrokinetic_M2PERP_2x2v_ser_p1, kernel_mom_gyrokinetic_M2PERP_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_mom_gyrokinetic_M2PERP_3x2v_ser_p1, NULL }, // 3
};

// M3 parallel kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m3par_kernels[] = {
  // 1x kernels
  { NULL, kernel_mom_gyrokinetic_M3PAR_1x1v_ser_p1, kernel_mom_gyrokinetic_M3PAR_1x1v_ser_p2 }, // 0
  { NULL, kernel_mom_gyrokinetic_M3PAR_1x2v_ser_p1, kernel_mom_gyrokinetic_M3PAR_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_mom_gyrokinetic_M3PAR_2x2v_ser_p1, kernel_mom_gyrokinetic_M3PAR_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_mom_gyrokinetic_M3PAR_3x2v_ser_p1, NULL }, // 3
};

// M3 perpendicular kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m3perp_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, kernel_mom_gyrokinetic_M3PERP_1x2v_ser_p1, kernel_mom_gyrokinetic_M3PERP_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_mom_gyrokinetic_M3PERP_2x2v_ser_p1, kernel_mom_gyrokinetic_M3PERP_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_mom_gyrokinetic_M3PERP_3x2v_ser_p1, NULL }, // 3
};

// Density, parallel momentum and kinetic energy computed together.
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m0m1m2_kernels[] = {
  // 1x kernels
  { NULL, kernel_mom_gyrokinetic_M0M1M2_1x1v_ser_p1, kernel_mom_gyrokinetic_M0M1M2_1x1v_ser_p2 }, // 0
  { NULL, kernel_mom_gyrokinetic_M0M1M2_1x2v_ser_p1, kernel_mom_gyrokinetic_M0M1M2_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_mom_gyrokinetic_M0M1M2_2x2v_ser_p1, kernel_mom_gyrokinetic_M0M1M2_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_mom_gyrokinetic_M0M1M2_3x2v_ser_p1, NULL }, // 3
};

// Density, parallel momentum, kinetic energy and heat flux computed together.
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m0m1m2m3_kernels[] = {
  // 1x kernels
  { NULL, kernel_mom_gyrokinetic_M0M1M2M3_1x1v_ser_p1, kernel_mom_gyrokinetic_M0M1M2M3_1x1v_ser_p2 }, // 0
  { NULL, kernel_mom_gyrokinetic_M0M1M2M3_1x2v_ser_p1, kernel_mom_gyrokinetic_M0M1M2M3_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_mom_gyrokinetic_M0M1M2M3_2x2v_ser_p1, kernel_mom_gyrokinetic_M0M1M2M3_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_mom_gyrokinetic_M0M1M2M3_3x2v_ser_p1, NULL }, // 3
};

// Density, parallel momentum, and parallel and perpendicular kinetic energy computed together.
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m0m1m2parm2perp_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, kernel_mom_gyrokinetic_M0M1M2PARM2PERP_1x2v_ser_p1, kernel_mom_gyrokinetic_M0M1M2PARM2PERP_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_mom_gyrokinetic_M0M1M2PARM2PERP_2x2v_ser_p1, kernel_mom_gyrokinetic_M0M1M2PARM2PERP_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_mom_gyrokinetic_M0M1M2PARM2PERP_3x2v_ser_p1, NULL }, // 3
};

// Density, parallel momentum, and parallel and perpendicular kinetic energy,
// and parallel and perpendicular heat flux computed together.
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_m0m1m2parm2perpm3parm3perp_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, kernel_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_1x2v_ser_p1, kernel_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_2x2v_ser_p1, kernel_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_3x2v_ser_p1, NULL }, // 3
};

// Hamiltonian moment kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_hamiltonian_kernels[] = {
  // 1x kernels
  { NULL, kernel_mom_gyrokinetic_HAMILTONIAN_1x1v_ser_p1, kernel_mom_gyrokinetic_HAMILTONIAN_1x1v_ser_p2 }, // 0
  { NULL, kernel_mom_gyrokinetic_HAMILTONIAN_1x2v_ser_p1, kernel_mom_gyrokinetic_HAMILTONIAN_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_mom_gyrokinetic_HAMILTONIAN_2x2v_ser_p1, kernel_mom_gyrokinetic_HAMILTONIAN_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_mom_gyrokinetic_HAMILTONIAN_3x2v_ser_p1, NULL }, // 3
};

// Integrated M0 kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_int_m0_kernels[] = {
  // 1x kernels
  { NULL, kernel_int_mom_gyrokinetic_M0_1x1v_ser_p1, kernel_int_mom_gyrokinetic_M0_1x1v_ser_p2 }, // 0
  { NULL, kernel_int_mom_gyrokinetic_M0_1x2v_ser_p1, kernel_int_mom_gyrokinetic_M0_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_int_mom_gyrokinetic_M0_2x2v_ser_p1, kernel_int_mom_gyrokinetic_M0_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_int_mom_gyrokinetic_M0_3x2v_ser_p1, NULL }, // 3
};

// Integrated M1 kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_int_m1_kernels[] = {
  // 1x kernels
  { NULL, kernel_int_mom_gyrokinetic_M1_1x1v_ser_p1, kernel_int_mom_gyrokinetic_M1_1x1v_ser_p2 }, // 0
  { NULL, kernel_int_mom_gyrokinetic_M1_1x2v_ser_p1, kernel_int_mom_gyrokinetic_M1_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_int_mom_gyrokinetic_M1_2x2v_ser_p1, kernel_int_mom_gyrokinetic_M1_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_int_mom_gyrokinetic_M1_3x2v_ser_p1, NULL }, // 3
};

// Integrated M2 kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_int_m2_kernels[] = {
  // 1x kernels
  { NULL, kernel_int_mom_gyrokinetic_M2_1x1v_ser_p1, kernel_int_mom_gyrokinetic_M2_1x1v_ser_p2 }, // 0
  { NULL, kernel_int_mom_gyrokinetic_M2_1x2v_ser_p1, kernel_int_mom_gyrokinetic_M2_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_int_mom_gyrokinetic_M2_2x2v_ser_p1, kernel_int_mom_gyrokinetic_M2_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_int_mom_gyrokinetic_M2_3x2v_ser_p1, NULL }, // 3
};

// Integrated M3 kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_int_m3_kernels[] = {
  // 1x kernels
  { NULL, kernel_int_mom_gyrokinetic_M3_1x1v_ser_p1, kernel_int_mom_gyrokinetic_M3_1x1v_ser_p2 }, // 0
  { NULL, kernel_int_mom_gyrokinetic_M3_1x2v_ser_p1, kernel_int_mom_gyrokinetic_M3_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_int_mom_gyrokinetic_M3_2x2v_ser_p1, kernel_int_mom_gyrokinetic_M3_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_int_mom_gyrokinetic_M3_3x2v_ser_p1, NULL }, // 3
};

// Integrated M2 parallel kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_int_m2par_kernels[] = {
  // 1x kernels
  { NULL, kernel_int_mom_gyrokinetic_M2PAR_1x1v_ser_p1, kernel_int_mom_gyrokinetic_M2PAR_1x1v_ser_p2 }, // 0
  { NULL, kernel_int_mom_gyrokinetic_M2PAR_1x2v_ser_p1, kernel_int_mom_gyrokinetic_M2PAR_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_int_mom_gyrokinetic_M2PAR_2x2v_ser_p1, kernel_int_mom_gyrokinetic_M2PAR_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_int_mom_gyrokinetic_M2PAR_3x2v_ser_p1, NULL }, // 3
};

// Integrated M2 perpendicular kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_int_m2perp_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, kernel_int_mom_gyrokinetic_M2PERP_1x2v_ser_p1, kernel_int_mom_gyrokinetic_M2PERP_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_int_mom_gyrokinetic_M2PERP_2x2v_ser_p1, kernel_int_mom_gyrokinetic_M2PERP_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_int_mom_gyrokinetic_M2PERP_3x2v_ser_p1, NULL }, // 3
};

// Integrated M3 parallel kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_int_m3par_kernels[] = {
  // 1x kernels
  { NULL, kernel_int_mom_gyrokinetic_M3PAR_1x1v_ser_p1, kernel_int_mom_gyrokinetic_M3PAR_1x1v_ser_p2 }, // 0
  { NULL, kernel_int_mom_gyrokinetic_M3PAR_1x2v_ser_p1, kernel_int_mom_gyrokinetic_M3PAR_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_int_mom_gyrokinetic_M3PAR_2x2v_ser_p1, kernel_int_mom_gyrokinetic_M3PAR_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_int_mom_gyrokinetic_M3PAR_3x2v_ser_p1, NULL }, // 3
};

// Integrated M3 perpendicular kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_int_m3perp_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, kernel_int_mom_gyrokinetic_M3PERP_1x2v_ser_p1, kernel_int_mom_gyrokinetic_M3PERP_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_int_mom_gyrokinetic_M3PERP_2x2v_ser_p1, kernel_int_mom_gyrokinetic_M3PERP_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_int_mom_gyrokinetic_M3PERP_3x2v_ser_p1, NULL }, // 3
};

// Integrated M0M1M2_moments kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_int_m0m1m2_kernels[] = {
  // 1x kernels
  { NULL, kernel_int_mom_gyrokinetic_M0M1M2_1x1v_ser_p1, kernel_int_mom_gyrokinetic_M0M1M2_1x1v_ser_p2 }, // 0
  { NULL, kernel_int_mom_gyrokinetic_M0M1M2_1x2v_ser_p1, kernel_int_mom_gyrokinetic_M0M1M2_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_int_mom_gyrokinetic_M0M1M2_2x2v_ser_p1, kernel_int_mom_gyrokinetic_M0M1M2_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_int_mom_gyrokinetic_M0M1M2_3x2v_ser_p1, NULL }, // 3
};

// Integrated M0M1M2M3_moments kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_int_m0m1m2m3_kernels[] = {
  // 1x kernels
  { NULL, kernel_int_mom_gyrokinetic_M0M1M2M3_1x1v_ser_p1, kernel_int_mom_gyrokinetic_M0M1M2M3_1x1v_ser_p2 }, // 0
  { NULL, kernel_int_mom_gyrokinetic_M0M1M2M3_1x2v_ser_p1, kernel_int_mom_gyrokinetic_M0M1M2M3_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_int_mom_gyrokinetic_M0M1M2M3_2x2v_ser_p1, kernel_int_mom_gyrokinetic_M0M1M2M3_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_int_mom_gyrokinetic_M0M1M2M3_3x2v_ser_p1, NULL }, // 3
};

// Integrated M0M1M2PARM2PERP_moments kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_int_m0m1m2parm2perp_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, kernel_int_mom_gyrokinetic_M0M1M2PARM2PERP_1x2v_ser_p1, kernel_int_mom_gyrokinetic_M0M1M2PARM2PERP_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_int_mom_gyrokinetic_M0M1M2PARM2PERP_2x2v_ser_p1, kernel_int_mom_gyrokinetic_M0M1M2PARM2PERP_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_int_mom_gyrokinetic_M0M1M2PARM2PERP_3x2v_ser_p1, NULL }, // 3
};

// Integrated M0M1M2PARM2PERPM3PARM3PERP_moments kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_int_m0m1m2parm2perpm3parm3perp_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, kernel_int_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_1x2v_ser_p1, kernel_int_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_int_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_2x2v_ser_p1, kernel_int_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_int_mom_gyrokinetic_M0M1M2PARM2PERPM3PARM3PERP_3x2v_ser_p1, NULL }, // 3
};

// Integrated Hamiltonian kernel list
GKYL_CU_D
static const gkyl_gyrokinetic_mom_kern_list ser_int_hamiltonian_kernels[] = {
  // 1x kernels
  { NULL, kernel_int_mom_gyrokinetic_HAMILTONIAN_1x1v_ser_p1, kernel_int_mom_gyrokinetic_HAMILTONIAN_1x1v_ser_p2 }, // 0
  { NULL, kernel_int_mom_gyrokinetic_HAMILTONIAN_1x2v_ser_p1, kernel_int_mom_gyrokinetic_HAMILTONIAN_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_int_mom_gyrokinetic_HAMILTONIAN_2x2v_ser_p1, kernel_int_mom_gyrokinetic_HAMILTONIAN_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_int_mom_gyrokinetic_HAMILTONIAN_3x2v_ser_p1, NULL }, // 3
};

/**
 * Free moment object.
 *
 * @param ref Reference counter for moment to free
 */
void gkyl_gk_mom_free(const struct gkyl_ref_count *ref);

#ifdef GKYL_HAVE_CUDA
/**
 * Create new gyrokinetic moment type object on NV-GPU:
 * see new() method above for documentation.
 */
struct gkyl_mom_type*
gkyl_mom_gyrokinetic_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range,
  double mass, double charge, const struct gkyl_velocity_map *vel_map,
  const struct gk_geometry *gk_geom, struct gkyl_array *phi, enum gkyl_distribution_moments mom_type);

/**
 * Create new integrated Gyrokinetic moment type object on NV-GPU:
 * see new() method above for documentation.
 */
struct gkyl_mom_type*
gkyl_int_mom_gyrokinetic_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range,
  double mass, double charge, const struct gkyl_velocity_map *vel_map,
  const struct gk_geometry *gk_geom, struct gkyl_array *phi, enum gkyl_distribution_moments mom_type);
#endif    
