#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void nc_hamil_alpha_quad_vx_1x1v_ser_p2(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx10 = 2.0/dxv[0]; 

  const double wx1 = w[1]; 
  const double dv10 = 2.0/dxv[1]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[3]; 
  alpha_quad[0] += -dx10*((1.224744871391589*hamil[1]-4.743416490252569*hamil[2])*(0.9681844646844028*poisson_tensor_conf_0[2]-1.054672281193885*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(wx1-0.5*dv10)); 
  alpha_quad[1] += -dx10*((1.224744871391589*hamil[1]-4.743416490252569*hamil[2])*((-0.5164305132317774*poisson_tensor_conf_0[2])-0.416390039500913*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(wx1-0.5*dv10)); 
  alpha_quad[2] += -dx10*((1.224744871391589*hamil[1]-4.743416490252569*hamil[2])*((-0.5164305132317774*poisson_tensor_conf_0[2])+0.416390039500913*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(wx1-0.5*dv10)); 
  alpha_quad[3] += -dx10*((1.224744871391589*hamil[1]-4.743416490252569*hamil[2])*(0.9681844646844028*poisson_tensor_conf_0[2]+1.054672281193885*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(wx1-0.5*dv10)); 

} 
