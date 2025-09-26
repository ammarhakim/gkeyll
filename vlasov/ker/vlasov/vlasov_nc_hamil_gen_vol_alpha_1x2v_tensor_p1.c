#include <gkyl_nc_hamil_vol_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_alpha_1x2v_tensor_p1(const double *w, const double *dxv, const int dir,
   const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dv0 = dxv[1]; 
  const double dv1 = dxv[2]; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  const double wx1 = w[1]; 
  const double wx2 = w[2]; 
  
  const double *pt_conf_0 = &poisson_tensor_conf[0]; 
  if (dir == 0) { 
  out[0] += 1.732050807568877*pt_conf_0[0]*hamil[1]*dv10; 
  out[1] += 1.732050807568877*hamil[1]*pt_conf_0[1]*dv10; 
  out[2] += 3.872983346207417*pt_conf_0[0]*hamil[4]*dv10; 
  out[3] += 1.732050807568877*pt_conf_0[0]*hamil[3]*dv10; 
  out[4] += 3.872983346207417*pt_conf_0[1]*hamil[4]*dv10; 
  out[5] += 1.732050807568877*pt_conf_0[1]*hamil[3]*dv10; 
  out[6] += 3.872983346207417*pt_conf_0[0]*hamil[6]*dv10; 
  out[8] += 1.732050807568877*pt_conf_0[0]*hamil[7]*dv10; 
  out[9] += 3.872983346207417*pt_conf_0[1]*hamil[6]*dv10; 
  out[12] += 1.732050807568877*pt_conf_0[1]*hamil[7]*dv10; 
  out[13] += 3.872983346207417*pt_conf_0[0]*hamil[8]*dv10; 
  out[15] += 3.872983346207417*pt_conf_0[1]*hamil[8]*dv10; 
  } 
  if (dir == 1) { 
  } 

  const double *pt_conf_1 = &poisson_tensor_conf[2]; 
  if (dir == 0) { 
  out[0] += 1.732050807568877*pt_conf_1[0]*hamil[2]*dv11; 
  out[1] += 1.732050807568877*pt_conf_1[1]*hamil[2]*dv11; 
  out[2] += 1.732050807568877*pt_conf_1[0]*hamil[3]*dv11; 
  out[3] += 3.872983346207417*pt_conf_1[0]*hamil[5]*dv11; 
  out[4] += 1.732050807568877*pt_conf_1[1]*hamil[3]*dv11; 
  out[5] += 3.872983346207417*pt_conf_1[1]*hamil[5]*dv11; 
  out[6] += 3.872983346207417*pt_conf_1[0]*hamil[7]*dv11; 
  out[7] += 1.732050807568877*pt_conf_1[0]*hamil[6]*dv11; 
  out[9] += 3.872983346207417*pt_conf_1[1]*hamil[7]*dv11; 
  out[10] += 1.732050807568877*pt_conf_1[1]*hamil[6]*dv11; 
  out[11] += 3.872983346207417*pt_conf_1[0]*hamil[8]*dv11; 
  out[14] += 3.872983346207417*pt_conf_1[1]*hamil[8]*dv11; 
  } 
  if (dir == 2) { 
  } 

  const double *pt_conf_4 = &poisson_tensor_conf[8]; 
  if (dir == 1) { 
  out[0] += dv11*(1.732050807568877*pt_conf_4[0]*hamil[2]*wx1+0.5*pt_conf_4[0]*hamil[3]*dv0); 
  out[1] += dv11*(1.732050807568877*pt_conf_4[1]*hamil[2]*wx1+0.5*pt_conf_4[1]*hamil[3]*dv0); 
  out[2] += dv11*(1.732050807568877*pt_conf_4[0]*hamil[3]*wx1+0.447213595499958*pt_conf_4[0]*hamil[6]*dv0+0.5*pt_conf_4[0]*hamil[2]*dv0); 
  out[3] += dv11*(3.872983346207417*pt_conf_4[0]*hamil[5]*wx1+1.118033988749895*pt_conf_4[0]*hamil[7]*dv0); 
  out[4] += dv11*(1.732050807568877*pt_conf_4[1]*hamil[3]*wx1+0.447213595499958*pt_conf_4[1]*hamil[6]*dv0+0.5*pt_conf_4[1]*hamil[2]*dv0); 
  out[5] += dv11*(3.872983346207417*pt_conf_4[1]*hamil[5]*wx1+1.118033988749895*pt_conf_4[1]*hamil[7]*dv0); 
  out[6] += dv11*(3.872983346207417*pt_conf_4[0]*hamil[7]*wx1+pt_conf_4[0]*hamil[8]*dv0+1.118033988749895*pt_conf_4[0]*hamil[5]*dv0); 
  out[7] += dv11*(1.732050807568877*pt_conf_4[0]*hamil[6]*wx1+0.4472135954999579*pt_conf_4[0]*hamil[3]*dv0); 
  out[9] += dv11*(3.872983346207417*pt_conf_4[1]*hamil[7]*wx1+pt_conf_4[1]*hamil[8]*dv0+1.118033988749895*pt_conf_4[1]*hamil[5]*dv0); 
  out[10] += dv11*(1.732050807568877*pt_conf_4[1]*hamil[6]*wx1+0.447213595499958*pt_conf_4[1]*hamil[3]*dv0); 
  out[11] += dv11*(3.872983346207417*pt_conf_4[0]*hamil[8]*wx1+pt_conf_4[0]*hamil[7]*dv0); 
  out[14] += dv11*(3.872983346207417*pt_conf_4[1]*hamil[8]*wx1+1.0*pt_conf_4[1]*hamil[7]*dv0); 
  } 
  if (dir == 2) { 
  out[0] += dv11*((-1.732050807568877*pt_conf_4[0]*hamil[1]*wx1)-1.118033988749895*pt_conf_4[0]*hamil[4]*dv0); 
  out[1] += dv11*((-1.732050807568877*hamil[1]*pt_conf_4[1]*wx1)-1.118033988749895*pt_conf_4[1]*hamil[4]*dv0); 
  out[2] += dv11*((-3.872983346207417*pt_conf_4[0]*hamil[4]*wx1)-0.5*pt_conf_4[0]*hamil[1]*dv0); 
  out[3] += dv11*((-1.732050807568877*pt_conf_4[0]*hamil[3]*wx1)-1.118033988749895*pt_conf_4[0]*hamil[6]*dv0); 
  out[4] += dv11*((-3.872983346207417*pt_conf_4[1]*hamil[4]*wx1)-0.5*hamil[1]*pt_conf_4[1]*dv0); 
  out[5] += dv11*((-1.732050807568877*pt_conf_4[1]*hamil[3]*wx1)-1.118033988749895*pt_conf_4[1]*hamil[6]*dv0); 
  out[6] += dv11*((-3.872983346207417*pt_conf_4[0]*hamil[6]*wx1)-0.5*pt_conf_4[0]*hamil[3]*dv0); 
  out[7] += -1.0*pt_conf_4[0]*hamil[4]*dv0*dv11; 
  out[8] += dv11*((-1.732050807568877*pt_conf_4[0]*hamil[7]*wx1)-1.118033988749895*pt_conf_4[0]*hamil[8]*dv0); 
  out[9] += dv11*((-3.872983346207417*pt_conf_4[1]*hamil[6]*wx1)-0.5*pt_conf_4[1]*hamil[3]*dv0); 
  out[10] += -1.0*pt_conf_4[1]*hamil[4]*dv0*dv11; 
  out[11] += -1.0*pt_conf_4[0]*hamil[6]*dv0*dv11; 
  out[12] += dv11*((-1.732050807568877*pt_conf_4[1]*hamil[7]*wx1)-1.118033988749895*pt_conf_4[1]*hamil[8]*dv0); 
  out[13] += dv11*((-3.872983346207417*pt_conf_4[0]*hamil[8]*wx1)-0.5*pt_conf_4[0]*hamil[7]*dv0); 
  out[14] += -1.0*pt_conf_4[1]*hamil[6]*dv0*dv11; 
  out[15] += dv11*((-3.872983346207417*pt_conf_4[1]*hamil[8]*wx1)-0.5000000000000001*pt_conf_4[1]*hamil[7]*dv0); 
  out[16] += -1.0*pt_conf_4[0]*hamil[8]*dv0*dv11; 
  out[17] += -1.0*pt_conf_4[1]*hamil[8]*dv0*dv11; 
  } 

  const double *pt_conf_5 = &poisson_tensor_conf[10]; 
  if (dir == 1) { 
  out[0] += dv11*(1.732050807568877*pt_conf_5[0]*hamil[2]*wx2+1.118033988749895*pt_conf_5[0]*hamil[5]*dv1); 
  out[1] += dv11*(1.732050807568877*pt_conf_5[1]*hamil[2]*wx2+1.118033988749895*pt_conf_5[1]*hamil[5]*dv1); 
  out[2] += dv11*(1.732050807568877*pt_conf_5[0]*hamil[3]*wx2+1.118033988749895*pt_conf_5[0]*hamil[7]*dv1); 
  out[3] += dv11*(3.872983346207417*pt_conf_5[0]*hamil[5]*wx2+0.5*pt_conf_5[0]*hamil[2]*dv1); 
  out[4] += dv11*(1.732050807568877*pt_conf_5[1]*hamil[3]*wx2+1.118033988749895*pt_conf_5[1]*hamil[7]*dv1); 
  out[5] += dv11*(3.872983346207417*pt_conf_5[1]*hamil[5]*wx2+0.5*pt_conf_5[1]*hamil[2]*dv1); 
  out[6] += dv11*(3.872983346207417*pt_conf_5[0]*hamil[7]*wx2+0.5*pt_conf_5[0]*hamil[3]*dv1); 
  out[7] += dv11*(1.732050807568877*pt_conf_5[0]*hamil[6]*wx2+1.118033988749895*pt_conf_5[0]*hamil[8]*dv1); 
  out[8] += pt_conf_5[0]*hamil[5]*dv1*dv11; 
  out[9] += dv11*(3.872983346207417*pt_conf_5[1]*hamil[7]*wx2+0.5*pt_conf_5[1]*hamil[3]*dv1); 
  out[10] += dv11*(1.732050807568877*pt_conf_5[1]*hamil[6]*wx2+1.118033988749895*pt_conf_5[1]*hamil[8]*dv1); 
  out[11] += dv11*(3.872983346207417*pt_conf_5[0]*hamil[8]*wx2+0.5*pt_conf_5[0]*hamil[6]*dv1); 
  out[12] += 1.0*pt_conf_5[1]*hamil[5]*dv1*dv11; 
  out[13] += pt_conf_5[0]*hamil[7]*dv1*dv11; 
  out[14] += dv11*(3.872983346207417*pt_conf_5[1]*hamil[8]*wx2+0.5000000000000001*pt_conf_5[1]*hamil[6]*dv1); 
  out[15] += 1.0*pt_conf_5[1]*hamil[7]*dv1*dv11; 
  out[16] += pt_conf_5[0]*hamil[8]*dv1*dv11; 
  out[17] += pt_conf_5[1]*hamil[8]*dv1*dv11; 
  } 
  if (dir == 2) { 
  out[0] += dv11*((-1.732050807568877*pt_conf_5[0]*hamil[1]*wx2)-0.5*pt_conf_5[0]*hamil[3]*dv1); 
  out[1] += dv11*((-1.732050807568877*hamil[1]*pt_conf_5[1]*wx2)-0.5*pt_conf_5[1]*hamil[3]*dv1); 
  out[2] += dv11*((-3.872983346207417*pt_conf_5[0]*hamil[4]*wx2)-1.118033988749895*pt_conf_5[0]*hamil[6]*dv1); 
  out[3] += dv11*((-1.732050807568877*pt_conf_5[0]*hamil[3]*wx2)-0.447213595499958*pt_conf_5[0]*hamil[7]*dv1-0.5*pt_conf_5[0]*hamil[1]*dv1); 
  out[4] += dv11*((-3.872983346207417*pt_conf_5[1]*hamil[4]*wx2)-1.118033988749895*pt_conf_5[1]*hamil[6]*dv1); 
  out[5] += dv11*((-1.732050807568877*pt_conf_5[1]*hamil[3]*wx2)-0.447213595499958*pt_conf_5[1]*hamil[7]*dv1-0.5*hamil[1]*pt_conf_5[1]*dv1); 
  out[6] += dv11*((-3.872983346207417*pt_conf_5[0]*hamil[6]*wx2)-1.0*pt_conf_5[0]*hamil[8]*dv1-1.118033988749895*pt_conf_5[0]*hamil[4]*dv1); 
  out[8] += dv11*((-1.732050807568877*pt_conf_5[0]*hamil[7]*wx2)-0.4472135954999579*pt_conf_5[0]*hamil[3]*dv1); 
  out[9] += dv11*((-3.872983346207417*pt_conf_5[1]*hamil[6]*wx2)-1.0*pt_conf_5[1]*hamil[8]*dv1-1.118033988749895*pt_conf_5[1]*hamil[4]*dv1); 
  out[12] += dv11*((-1.732050807568877*pt_conf_5[1]*hamil[7]*wx2)-0.447213595499958*pt_conf_5[1]*hamil[3]*dv1); 
  out[13] += dv11*((-3.872983346207417*pt_conf_5[0]*hamil[8]*wx2)-1.0*pt_conf_5[0]*hamil[6]*dv1); 
  out[15] += dv11*((-3.872983346207417*pt_conf_5[1]*hamil[8]*wx2)-1.0*pt_conf_5[1]*hamil[6]*dv1); 
  } 

} 
