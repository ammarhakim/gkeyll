#include <gkyl_dg_gyrokinetic_kernels.h> 
GKYL_CU_DH double dg_gyrokinetic_add_em_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
    const double q_, const double m_, const double *bmag, const double *phi,
    const double *apar, const double *rtg33inv, const double *bioverJB,
    const double *fin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // vmap: velocity space mapping.
  // vmapSq: velocity space mapping squared.
  // q_,m_: species charge and mass.
  // bmag: magnetic field amplitude.
  // phi: electrostatic potential .
  // apar: parallel component of vector potential.
  // fin: Distribution function.
  // out: output increment.

  return 0.; 
} 
