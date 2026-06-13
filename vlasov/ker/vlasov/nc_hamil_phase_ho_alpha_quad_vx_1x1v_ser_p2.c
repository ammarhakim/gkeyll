#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void nc_hamil_phase_ho_alpha_quad_vx_1x1v_ser_p2(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 

  const double wx1 = w[1]; 
  const double dv0 = dxv[1]; 
  const double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 
  const double *poisson_tensor_conf_x0 = &poisson_tensor_conf[0]; 
  alpha_quad[0] += -1.0*(0.9681844646844028*poisson_tensor_conf_x0[2]-1.054672281193885*poisson_tensor_conf_x0[1]+0.7071067811865475*poisson_tensor_conf_x0[0])*(1.936491673103709*hamil[7]+5.002749890427371*hamil[6]-2.888338995926613*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1])*dx10; 
  alpha_quad[1] += -1.0*((-0.5164305132317774*poisson_tensor_conf_x0[2])-0.416390039500913*poisson_tensor_conf_x0[1]+0.7071067811865475*poisson_tensor_conf_x0[0])*(1.936491673103709*hamil[7]+1.97511137974555*hamil[6]-1.140331086775586*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1])*dx10; 
  alpha_quad[2] += -1.0*((-0.5164305132317774*poisson_tensor_conf_x0[2])+0.416390039500913*poisson_tensor_conf_x0[1]+0.7071067811865475*poisson_tensor_conf_x0[0])*(1.936491673103709*hamil[7]-1.97511137974555*hamil[6]+1.140331086775586*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1])*dx10; 
  alpha_quad[3] += -1.0*(0.9681844646844028*poisson_tensor_conf_x0[2]+1.054672281193885*poisson_tensor_conf_x0[1]+0.7071067811865475*poisson_tensor_conf_x0[0])*(1.936491673103709*hamil[7]-5.002749890427371*hamil[6]+2.888338995926613*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1])*dx10; 

} 
