#pragma once

// Private header, not for direct use in user code

#include <gkyl_dg_cx_kernels.h>
#include <gkyl_util.h>
#include <assert.h>

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1, -1,  3,  4}, // 2x kernel indices
  {-1, -1, -1,  5}, // 3x kernel indices  
};

typedef double (*dg_cx_react_ratef_t)(const double a, const double b, double vt_sq_ion_min, double vt_sq_neut_min, 
  const double *prim_vars_ion, const double *prim_vars_neut, const double *u_ion, double* GKYL_RESTRICT v_sigma_cx) ;

// for use in kernel tables
typedef struct { dg_cx_react_ratef_t kernels[3]; } gkyl_cx_react_rate_kern_list;

//
// Serendipity basis kernels
// 

// CX reaction rate kernel list
GKYL_CU_D
static const gkyl_cx_react_rate_kern_list ser_cx_react_rate_kernels[] = {
  { NULL, sigma_cx_1x1v_ser_p1, sigma_cx_1x1v_ser_p2 }, // 0
  { NULL, sigma_cx_1x2v_ser_p1, sigma_cx_1x2v_ser_p2 }, // 1
  { NULL, sigma_cx_1x3v_ser_p1, sigma_cx_1x3v_ser_p2 }, // 2
  { NULL, sigma_cx_2x2v_ser_p1, sigma_cx_2x2v_ser_p2 }, // 3
  { NULL, sigma_cx_2x3v_ser_p1, sigma_cx_2x3v_ser_p2 }, // 4
  { NULL, sigma_cx_3x3v_ser_p1, NULL }, // 5
};

struct gkyl_dg_cx {
  const struct gkyl_rect_grid *grid; // grid object

  struct gkyl_basis *cbasis;
  struct gkyl_basis *pbasis_gk;
  struct gkyl_basis *pbasis_vl;
  const struct gkyl_range *conf_rng;
  const struct gkyl_range *conf_rng_ext;
  const struct gkyl_range *phase_rng; 
  
  double a; // Fitting function coefficient.
  double b; // Fitting function coefficient.
  double mass_ion;
  double mass_neut;
  double vt_sq_ion_min;
  double vt_sq_neut_min; 
  
  enum gkyl_ion_type type_ion;
  enum gkyl_react_self_type type_self;

  uint32_t flags;
  dg_cx_react_ratef_t react_rate; // pointer to reaction rate kernel

  struct gkyl_dg_cx *on_dev; // pointer to itself or device data

};

GKYL_CU_D
static dg_cx_react_ratef_t
choose_kern(enum gkyl_basis_type b_type, int tblidx, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_HYBRID:
    case GKYL_BASIS_MODAL_TENSOR:
      return ser_cx_react_rate_kernels[tblidx].kernels[poly_order];
      break;
    default:
      assert(false);
      break; 
  }
  return 0;
}

static void
fit_param(enum gkyl_ion_type type_ion, double *a, double *b)
{
  // These values are from E. Meier's PhD Thesis
  if (type_ion == GKYL_ION_H) {
    a[0] = 1.12e-18;
    b[0] = 7.15e-20;
  }
  else if (type_ion == GKYL_ION_D) {
    a[0] = 1.09e-18;
    b[0] = 7.15e-20;
  }
  else if (type_ion == GKYL_ION_HE) {
    a[0] = 6.484e-19;
    b[0] = 4.350e-20;
  } 
  else if (type_ion == GKYL_ION_NE) {
    a[0] = 7.95e-19;
    b[0] = 5.65e-20;
  }
}

