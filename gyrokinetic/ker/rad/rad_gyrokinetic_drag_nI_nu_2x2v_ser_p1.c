#include <gkyl_rad_gyrokinetic_kernels.h> 
GKYL_CU_DH void rad_gyrokinetic_drag_nI_nu_2x2v_ser_p1(const double *vnu_surf, const double *vnu,
    const double *vsqnu_surf, const double *vsqnu, const double *nI, 
    double* GKYL_RESTRICT nvnu_surf, double* GKYL_RESTRICT nvnu, 
    double* GKYL_RESTRICT nvsqnu_surf, double* GKYL_RESTRICT nvsqnu) 
{ 
  // vnu_surf/vnu: Input radiation drag in vparallel direction, surface and volume expansions.
  // vsqnu_surf/vsqnu: Input radiation drag in mu direction, surface and volume expansions.
  // nI: Input ion density.
  // nvnu_surf/nvnu: Accumulated output density-weighted radiation drag in vparallel direction, surface and volume expansions.
  // nvsqnu_surf/nvsqnu: Accumulated output density-weighted radiation drag in mu direction, surface and volume expansions.

  nvnu_surf[0] += 0.5*(nI[3]*vnu_surf[4]+nI[2]*vnu_surf[2]+nI[1]*vnu_surf[1]+nI[0]*vnu_surf[0]); 
  nvnu_surf[1] += 0.5*(nI[2]*vnu_surf[4]+vnu_surf[2]*nI[3]+nI[0]*vnu_surf[1]+vnu_surf[0]*nI[1]); 
  nvnu_surf[2] += 0.5*(nI[1]*vnu_surf[4]+vnu_surf[1]*nI[3]+nI[0]*vnu_surf[2]+vnu_surf[0]*nI[2]); 
  nvnu_surf[3] += 0.5*(nI[3]*vnu_surf[7]+nI[2]*vnu_surf[6]+nI[1]*vnu_surf[5]+nI[0]*vnu_surf[3]); 
  nvnu_surf[4] += 0.5*(nI[0]*vnu_surf[4]+vnu_surf[0]*nI[3]+nI[1]*vnu_surf[2]+vnu_surf[1]*nI[2]); 
  nvnu_surf[5] += 0.5*(nI[2]*vnu_surf[7]+nI[3]*vnu_surf[6]+nI[0]*vnu_surf[5]+nI[1]*vnu_surf[3]); 
  nvnu_surf[6] += 0.5*(nI[1]*vnu_surf[7]+nI[0]*vnu_surf[6]+nI[3]*vnu_surf[5]+nI[2]*vnu_surf[3]); 
  nvnu_surf[7] += 0.5*(nI[0]*vnu_surf[7]+nI[1]*vnu_surf[6]+nI[2]*vnu_surf[5]+nI[3]*vnu_surf[3]); 

  nvnu[0] += 0.5*(nI[3]*vnu[5]+nI[2]*vnu[2]+nI[1]*vnu[1]+nI[0]*vnu[0]); 
  nvnu[1] += 0.5*(nI[2]*vnu[5]+vnu[2]*nI[3]+nI[0]*vnu[1]+vnu[0]*nI[1]); 
  nvnu[2] += 0.5*(nI[1]*vnu[5]+vnu[1]*nI[3]+nI[0]*vnu[2]+vnu[0]*nI[2]); 
  nvnu[3] += 0.5*(nI[3]*vnu[11]+nI[2]*vnu[7]+nI[1]*vnu[6]+nI[0]*vnu[3]); 
  nvnu[4] += 0.5*(nI[3]*vnu[12]+nI[2]*vnu[9]+nI[1]*vnu[8]+nI[0]*vnu[4]); 
  nvnu[5] += 0.5*(nI[0]*vnu[5]+vnu[0]*nI[3]+nI[1]*vnu[2]+vnu[1]*nI[2]); 
  nvnu[6] += 0.5*(nI[2]*vnu[11]+nI[3]*vnu[7]+nI[0]*vnu[6]+nI[1]*vnu[3]); 
  nvnu[7] += 0.5*(nI[1]*vnu[11]+nI[0]*vnu[7]+nI[3]*vnu[6]+nI[2]*vnu[3]); 
  nvnu[8] += 0.5*(nI[2]*vnu[12]+nI[3]*vnu[9]+nI[0]*vnu[8]+nI[1]*vnu[4]); 
  nvnu[9] += 0.5*(nI[1]*vnu[12]+nI[0]*vnu[9]+nI[3]*vnu[8]+nI[2]*vnu[4]); 
  nvnu[10] += 0.5*(nI[3]*vnu[15]+nI[2]*vnu[14]+nI[1]*vnu[13]+nI[0]*vnu[10]); 
  nvnu[11] += 0.5*(nI[0]*vnu[11]+nI[1]*vnu[7]+nI[2]*vnu[6]+nI[3]*vnu[3]); 
  nvnu[12] += 0.5*(nI[0]*vnu[12]+nI[1]*vnu[9]+nI[2]*vnu[8]+nI[3]*vnu[4]); 
  nvnu[13] += 0.5*(nI[2]*vnu[15]+nI[3]*vnu[14]+nI[0]*vnu[13]+nI[1]*vnu[10]); 
  nvnu[14] += 0.5*(nI[1]*vnu[15]+nI[0]*vnu[14]+nI[3]*vnu[13]+nI[2]*vnu[10]); 
  nvnu[15] += 0.5*(nI[0]*vnu[15]+nI[1]*vnu[14]+nI[2]*vnu[13]+nI[3]*vnu[10]); 
  nvnu[16] += 0.5*nI[3]*vnu[20]+0.5000000000000001*(nI[2]*vnu[18]+nI[1]*vnu[17])+0.5*nI[0]*vnu[16]; 
  nvnu[17] += 0.5000000000000001*nI[2]*vnu[20]+0.5*(nI[3]*vnu[18]+nI[0]*vnu[17])+0.5000000000000001*nI[1]*vnu[16]; 
  nvnu[18] += 0.5000000000000001*nI[1]*vnu[20]+0.5*(nI[0]*vnu[18]+nI[3]*vnu[17])+0.5000000000000001*nI[2]*vnu[16]; 
  nvnu[19] += 0.5*nI[3]*vnu[23]+0.5000000000000001*(nI[2]*vnu[22]+nI[1]*vnu[21])+0.5*nI[0]*vnu[19]; 
  nvnu[20] += 0.5*nI[0]*vnu[20]+0.5000000000000001*(nI[1]*vnu[18]+nI[2]*vnu[17])+0.5*nI[3]*vnu[16]; 
  nvnu[21] += 0.5000000000000001*nI[2]*vnu[23]+0.5*(nI[3]*vnu[22]+nI[0]*vnu[21])+0.5000000000000001*nI[1]*vnu[19]; 
  nvnu[22] += 0.5000000000000001*nI[1]*vnu[23]+0.5*(nI[0]*vnu[22]+nI[3]*vnu[21])+0.5000000000000001*nI[2]*vnu[19]; 
  nvnu[23] += 0.5*nI[0]*vnu[23]+0.5000000000000001*(nI[1]*vnu[22]+nI[2]*vnu[21])+0.5*nI[3]*vnu[19]; 

  nvsqnu_surf[0] += 0.5*(nI[3]*vsqnu_surf[4]+nI[2]*vsqnu_surf[2]+nI[1]*vsqnu_surf[1]+nI[0]*vsqnu_surf[0]); 
  nvsqnu_surf[1] += 0.5*(nI[2]*vsqnu_surf[4]+vsqnu_surf[2]*nI[3]+nI[0]*vsqnu_surf[1]+vsqnu_surf[0]*nI[1]); 
  nvsqnu_surf[2] += 0.5*(nI[1]*vsqnu_surf[4]+vsqnu_surf[1]*nI[3]+nI[0]*vsqnu_surf[2]+vsqnu_surf[0]*nI[2]); 
  nvsqnu_surf[3] += 0.5*(nI[3]*vsqnu_surf[7]+nI[2]*vsqnu_surf[6]+nI[1]*vsqnu_surf[5]+nI[0]*vsqnu_surf[3]); 
  nvsqnu_surf[4] += 0.5*(nI[0]*vsqnu_surf[4]+vsqnu_surf[0]*nI[3]+nI[1]*vsqnu_surf[2]+vsqnu_surf[1]*nI[2]); 
  nvsqnu_surf[5] += 0.5*(nI[2]*vsqnu_surf[7]+nI[3]*vsqnu_surf[6]+nI[0]*vsqnu_surf[5]+nI[1]*vsqnu_surf[3]); 
  nvsqnu_surf[6] += 0.5*(nI[1]*vsqnu_surf[7]+nI[0]*vsqnu_surf[6]+nI[3]*vsqnu_surf[5]+nI[2]*vsqnu_surf[3]); 
  nvsqnu_surf[7] += 0.5*(nI[0]*vsqnu_surf[7]+nI[1]*vsqnu_surf[6]+nI[2]*vsqnu_surf[5]+nI[3]*vsqnu_surf[3]); 
  nvsqnu_surf[8] += 0.5*nI[3]*vsqnu_surf[11]+0.5000000000000001*(nI[2]*vsqnu_surf[10]+nI[1]*vsqnu_surf[9])+0.5*nI[0]*vsqnu_surf[8]; 
  nvsqnu_surf[9] += 0.5000000000000001*nI[2]*vsqnu_surf[11]+0.5*(nI[3]*vsqnu_surf[10]+nI[0]*vsqnu_surf[9])+0.5000000000000001*nI[1]*vsqnu_surf[8]; 
  nvsqnu_surf[10] += 0.5000000000000001*nI[1]*vsqnu_surf[11]+0.5*(nI[0]*vsqnu_surf[10]+nI[3]*vsqnu_surf[9])+0.5000000000000001*nI[2]*vsqnu_surf[8]; 
  nvsqnu_surf[11] += 0.5*nI[0]*vsqnu_surf[11]+0.5000000000000001*(nI[1]*vsqnu_surf[10]+nI[2]*vsqnu_surf[9])+0.5*nI[3]*vsqnu_surf[8]; 

  nvsqnu[0] += 0.5*(nI[3]*vsqnu[5]+nI[2]*vsqnu[2]+nI[1]*vsqnu[1]+nI[0]*vsqnu[0]); 
  nvsqnu[1] += 0.5*(nI[2]*vsqnu[5]+vsqnu[2]*nI[3]+nI[0]*vsqnu[1]+vsqnu[0]*nI[1]); 
  nvsqnu[2] += 0.5*(nI[1]*vsqnu[5]+vsqnu[1]*nI[3]+nI[0]*vsqnu[2]+vsqnu[0]*nI[2]); 
  nvsqnu[3] += 0.5*(nI[3]*vsqnu[11]+nI[2]*vsqnu[7]+nI[1]*vsqnu[6]+nI[0]*vsqnu[3]); 
  nvsqnu[4] += 0.5*(nI[3]*vsqnu[12]+nI[2]*vsqnu[9]+nI[1]*vsqnu[8]+nI[0]*vsqnu[4]); 
  nvsqnu[5] += 0.5*(nI[0]*vsqnu[5]+vsqnu[0]*nI[3]+nI[1]*vsqnu[2]+vsqnu[1]*nI[2]); 
  nvsqnu[6] += 0.5*(nI[2]*vsqnu[11]+nI[3]*vsqnu[7]+nI[0]*vsqnu[6]+nI[1]*vsqnu[3]); 
  nvsqnu[7] += 0.5*(nI[1]*vsqnu[11]+nI[0]*vsqnu[7]+nI[3]*vsqnu[6]+nI[2]*vsqnu[3]); 
  nvsqnu[8] += 0.5*(nI[2]*vsqnu[12]+nI[3]*vsqnu[9]+nI[0]*vsqnu[8]+nI[1]*vsqnu[4]); 
  nvsqnu[9] += 0.5*(nI[1]*vsqnu[12]+nI[0]*vsqnu[9]+nI[3]*vsqnu[8]+nI[2]*vsqnu[4]); 
  nvsqnu[10] += 0.5*(nI[3]*vsqnu[15]+nI[2]*vsqnu[14]+nI[1]*vsqnu[13]+nI[0]*vsqnu[10]); 
  nvsqnu[11] += 0.5*(nI[0]*vsqnu[11]+nI[1]*vsqnu[7]+nI[2]*vsqnu[6]+nI[3]*vsqnu[3]); 
  nvsqnu[12] += 0.5*(nI[0]*vsqnu[12]+nI[1]*vsqnu[9]+nI[2]*vsqnu[8]+nI[3]*vsqnu[4]); 
  nvsqnu[13] += 0.5*(nI[2]*vsqnu[15]+nI[3]*vsqnu[14]+nI[0]*vsqnu[13]+nI[1]*vsqnu[10]); 
  nvsqnu[14] += 0.5*(nI[1]*vsqnu[15]+nI[0]*vsqnu[14]+nI[3]*vsqnu[13]+nI[2]*vsqnu[10]); 
  nvsqnu[15] += 0.5*(nI[0]*vsqnu[15]+nI[1]*vsqnu[14]+nI[2]*vsqnu[13]+nI[3]*vsqnu[10]); 
  nvsqnu[16] += 0.5*nI[3]*vsqnu[20]+0.5000000000000001*(nI[2]*vsqnu[18]+nI[1]*vsqnu[17])+0.5*nI[0]*vsqnu[16]; 
  nvsqnu[17] += 0.5000000000000001*nI[2]*vsqnu[20]+0.5*(nI[3]*vsqnu[18]+nI[0]*vsqnu[17])+0.5000000000000001*nI[1]*vsqnu[16]; 
  nvsqnu[18] += 0.5000000000000001*nI[1]*vsqnu[20]+0.5*(nI[0]*vsqnu[18]+nI[3]*vsqnu[17])+0.5000000000000001*nI[2]*vsqnu[16]; 
  nvsqnu[19] += 0.5*nI[3]*vsqnu[23]+0.5000000000000001*(nI[2]*vsqnu[22]+nI[1]*vsqnu[21])+0.5*nI[0]*vsqnu[19]; 
  nvsqnu[20] += 0.5*nI[0]*vsqnu[20]+0.5000000000000001*(nI[1]*vsqnu[18]+nI[2]*vsqnu[17])+0.5*nI[3]*vsqnu[16]; 
  nvsqnu[21] += 0.5000000000000001*nI[2]*vsqnu[23]+0.5*(nI[3]*vsqnu[22]+nI[0]*vsqnu[21])+0.5000000000000001*nI[1]*vsqnu[19]; 
  nvsqnu[22] += 0.5000000000000001*nI[1]*vsqnu[23]+0.5*(nI[0]*vsqnu[22]+nI[3]*vsqnu[21])+0.5000000000000001*nI[2]*vsqnu[19]; 
  nvsqnu[23] += 0.5*nI[0]*vsqnu[23]+0.5000000000000001*(nI[1]*vsqnu[22]+nI[2]*vsqnu[21])+0.5*nI[3]*vsqnu[19]; 

} 
