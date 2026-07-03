#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_1x1v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge, 
    const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  const double dv10 = 2.0/dxv[1]; 

  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  if (hamil_pt_edge == -1) {
    alpha_quad[0] += (1.224744871391589*hamil[1]-3.674234614174766*hamil[2])*(1.58113883008419*poisson_tensor_conf_0[2]-1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*dv10;
    alpha_quad[1] += 1.224744871391589*hamil[1]*(1.58113883008419*poisson_tensor_conf_0[2]-1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*dv10;
    alpha_quad[2] += (3.674234614174766*hamil[2]+1.224744871391589*hamil[1])*(1.58113883008419*poisson_tensor_conf_0[2]-1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*dv10;
  }
  else if (hamil_pt_edge == 1) {
    alpha_quad[0] += (1.224744871391589*hamil[1]-3.674234614174766*hamil[2])*(1.58113883008419*poisson_tensor_conf_0[2]+1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*dv10;
    alpha_quad[1] += 1.224744871391589*hamil[1]*(1.58113883008419*poisson_tensor_conf_0[2]+1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*dv10;
    alpha_quad[2] += (3.674234614174766*hamil[2]+1.224744871391589*hamil[1])*(1.58113883008419*poisson_tensor_conf_0[2]+1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*dv10;
  }

} 
