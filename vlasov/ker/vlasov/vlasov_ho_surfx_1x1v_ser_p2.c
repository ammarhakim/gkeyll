#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_ho_surfx_1x1v_ser_p2(const double *w, const double *dxv,
  const double *Fhat_l_nodal, const double *Fhat_r_nodal, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 

  const double *Fhat_l_nodal_c = &Fhat_l_nodal[0]; 
  const double *Fhat_r_nodal_c = &Fhat_r_nodal[0]; 
  double G1[3] = {0.0}; 
  double Ghat[3] = {0.0}; 
  G1[0] = 0.24597051986529*Fhat_l_nodal_c[3]+0.46113626132125757*Fhat_l_nodal_c[2]+0.46113626132125757*Fhat_l_nodal_c[1]+0.24597051986529*Fhat_l_nodal_c[0]; 
  G1[1] = 0.3668728630454641*Fhat_l_nodal_c[3]+0.271546746793545*Fhat_l_nodal_c[2]-0.271546746793545*Fhat_l_nodal_c[1]-0.3668728630454641*Fhat_l_nodal_c[0]; 
  G1[2] = 0.3367876570272824*Fhat_l_nodal_c[3]-0.3367876570272824*Fhat_l_nodal_c[2]-0.3367876570272824*Fhat_l_nodal_c[1]+0.3367876570272824*Fhat_l_nodal_c[0]; 
  Ghat[0] = G1[0]; 
  Ghat[1] = G1[1]; 
  Ghat[2] = G1[2]; 
  out[0] += 0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -(1.224744871391589*Ghat[0]*dx10); 
  out[2] += 0.7071067811865475*Ghat[1]*dx10; 
  out[3] += -(1.224744871391589*Ghat[1]*dx10); 
  out[4] += 1.5811388300841895*Ghat[0]*dx10; 
  out[5] += 0.7071067811865475*Ghat[2]*dx10; 
  out[6] += 1.5811388300841898*Ghat[1]*dx10; 
  out[7] += -(1.224744871391589*Ghat[2]*dx10); 

  G1[0] = 0.24597051986529*Fhat_r_nodal_c[3]+0.46113626132125757*Fhat_r_nodal_c[2]+0.46113626132125757*Fhat_r_nodal_c[1]+0.24597051986529*Fhat_r_nodal_c[0]; 
  G1[1] = 0.3668728630454641*Fhat_r_nodal_c[3]+0.271546746793545*Fhat_r_nodal_c[2]-0.271546746793545*Fhat_r_nodal_c[1]-0.3668728630454641*Fhat_r_nodal_c[0]; 
  G1[2] = 0.3367876570272824*Fhat_r_nodal_c[3]-0.3367876570272824*Fhat_r_nodal_c[2]-0.3367876570272824*Fhat_r_nodal_c[1]+0.3367876570272824*Fhat_r_nodal_c[0]; 
  Ghat[0] = G1[0]; 
  Ghat[1] = G1[1]; 
  Ghat[2] = G1[2]; 
  out[0] += -(0.7071067811865475*Ghat[0]*dx10); 
  out[1] += -(1.224744871391589*Ghat[0]*dx10); 
  out[2] += -(0.7071067811865475*Ghat[1]*dx10); 
  out[3] += -(1.224744871391589*Ghat[1]*dx10); 
  out[4] += -(1.5811388300841895*Ghat[0]*dx10); 
  out[5] += -(0.7071067811865475*Ghat[2]*dx10); 
  out[6] += -(1.5811388300841898*Ghat[1]*dx10); 
  out[7] += -(1.224744871391589*Ghat[2]*dx10); 

  return 0.0;

} 
