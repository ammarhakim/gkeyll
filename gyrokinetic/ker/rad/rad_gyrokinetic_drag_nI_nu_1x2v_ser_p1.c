#include <gkyl_rad_gyrokinetic_kernels.h> 
GKYL_CU_DH void rad_gyrokinetic_drag_nI_nu_1x2v_ser_p1(const double *vnu_surf, const double *vnu,
    const double *vsqnu_surf, const double *vsqnu, const double *nI, 
    double* GKYL_RESTRICT nvnu_surf, double* GKYL_RESTRICT nvnu, 
    double* GKYL_RESTRICT nvsqnu_surf, double* GKYL_RESTRICT nvsqnu) 
{ 
  // vnu_surf/vnu: Input radiation drag in vparallel direction, surface and volume expansions.
  // vsqnu_surf/vsqnu: Input radiation drag in mu direction, surface and volume expansions.
  // nI: Input ion density.
  // nvnu_surf/nvnu: Accumulated output density-weighted radiation drag in vparallel direction, surface and volume expansions.
  // nvsqnu_surf/nvsqnu: Accumulated output density-weighted radiation drag in mu direction, surface and volume expansions.

  nvnu_surf[0] += 0.7071067811865475*(nI[1]*vnu_surf[1]+nI[0]*vnu_surf[0]); 
  nvnu_surf[1] += 0.7071067811865475*(nI[0]*vnu_surf[1]+vnu_surf[0]*nI[1]); 
  nvnu_surf[2] += 0.7071067811865475*(nI[1]*vnu_surf[3]+nI[0]*vnu_surf[2]); 
  nvnu_surf[3] += 0.7071067811865475*(nI[0]*vnu_surf[3]+nI[1]*vnu_surf[2]); 

  nvnu[0] += 0.7071067811865475*(nI[1]*vnu[1]+nI[0]*vnu[0]); 
  nvnu[1] += 0.7071067811865475*(nI[0]*vnu[1]+vnu[0]*nI[1]); 
  nvnu[2] += 0.7071067811865475*(nI[1]*vnu[4]+nI[0]*vnu[2]); 
  nvnu[3] += 0.7071067811865475*(nI[1]*vnu[5]+nI[0]*vnu[3]); 
  nvnu[4] += 0.7071067811865475*(nI[0]*vnu[4]+nI[1]*vnu[2]); 
  nvnu[5] += 0.7071067811865475*(nI[0]*vnu[5]+nI[1]*vnu[3]); 
  nvnu[6] += 0.7071067811865475*(nI[1]*vnu[7]+nI[0]*vnu[6]); 
  nvnu[7] += 0.7071067811865475*(nI[0]*vnu[7]+nI[1]*vnu[6]); 
  nvnu[8] += 0.7071067811865475*(nI[1]*vnu[9]+nI[0]*vnu[8]); 
  nvnu[9] += 0.7071067811865475*(nI[0]*vnu[9]+nI[1]*vnu[8]); 
  nvnu[10] += 0.7071067811865475*(nI[1]*vnu[11]+nI[0]*vnu[10]); 
  nvnu[11] += 0.7071067811865475*(nI[0]*vnu[11]+nI[1]*vnu[10]); 

  nvsqnu_surf[0] += 0.7071067811865475*(nI[1]*vsqnu_surf[1]+nI[0]*vsqnu_surf[0]); 
  nvsqnu_surf[1] += 0.7071067811865475*(nI[0]*vsqnu_surf[1]+vsqnu_surf[0]*nI[1]); 
  nvsqnu_surf[2] += 0.7071067811865475*(nI[1]*vsqnu_surf[3]+nI[0]*vsqnu_surf[2]); 
  nvsqnu_surf[3] += 0.7071067811865475*(nI[0]*vsqnu_surf[3]+nI[1]*vsqnu_surf[2]); 
  nvsqnu_surf[4] += 0.7071067811865475*(nI[1]*vsqnu_surf[5]+nI[0]*vsqnu_surf[4]); 
  nvsqnu_surf[5] += 0.7071067811865475*(nI[0]*vsqnu_surf[5]+nI[1]*vsqnu_surf[4]); 

  nvsqnu[0] += 0.7071067811865475*(nI[1]*vsqnu[1]+nI[0]*vsqnu[0]); 
  nvsqnu[1] += 0.7071067811865475*(nI[0]*vsqnu[1]+vsqnu[0]*nI[1]); 
  nvsqnu[2] += 0.7071067811865475*(nI[1]*vsqnu[4]+nI[0]*vsqnu[2]); 
  nvsqnu[3] += 0.7071067811865475*(nI[1]*vsqnu[5]+nI[0]*vsqnu[3]); 
  nvsqnu[4] += 0.7071067811865475*(nI[0]*vsqnu[4]+nI[1]*vsqnu[2]); 
  nvsqnu[5] += 0.7071067811865475*(nI[0]*vsqnu[5]+nI[1]*vsqnu[3]); 
  nvsqnu[6] += 0.7071067811865475*(nI[1]*vsqnu[7]+nI[0]*vsqnu[6]); 
  nvsqnu[7] += 0.7071067811865475*(nI[0]*vsqnu[7]+nI[1]*vsqnu[6]); 
  nvsqnu[8] += 0.7071067811865475*(nI[1]*vsqnu[9]+nI[0]*vsqnu[8]); 
  nvsqnu[9] += 0.7071067811865475*(nI[0]*vsqnu[9]+nI[1]*vsqnu[8]); 
  nvsqnu[10] += 0.7071067811865475*(nI[1]*vsqnu[11]+nI[0]*vsqnu[10]); 
  nvsqnu[11] += 0.7071067811865475*(nI[0]*vsqnu[11]+nI[1]*vsqnu[10]); 

} 
