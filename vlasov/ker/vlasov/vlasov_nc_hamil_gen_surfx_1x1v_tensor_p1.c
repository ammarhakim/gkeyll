#include <gkyl_nc_hamil_surf_comps_kernels.h> 
#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_nc_hamil_gen_surfx_1x1v_tensor_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 

  double Ghat_r[3] = {0.0}; 
  double Ghat_l[3] = {0.0}; 
  double alpha_max = 0.0; 
  vlasov_nc_hamil_gen_surf_node0_x_1x1v_tensor_p1(w, dxv, poisson_tensor_conf, hamil, fl, fc, fr, &alpha_max, Ghat_l, Ghat_r); 
  vlasov_nc_hamil_gen_surf_node1_x_1x1v_tensor_p1(w, dxv, poisson_tensor_conf, hamil, fl, fc, fr, &alpha_max, Ghat_l, Ghat_r); 
  vlasov_nc_hamil_gen_surf_node2_x_1x1v_tensor_p1(w, dxv, poisson_tensor_conf, hamil, fl, fc, fr, &alpha_max, Ghat_l, Ghat_r); 
  vlasov_nc_hamil_gen_surf_node3_x_1x1v_tensor_p1(w, dxv, poisson_tensor_conf, hamil, fl, fc, fr, &alpha_max, Ghat_l, Ghat_r); 
  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx10; 
  out[1] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx10; 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx10; 
  out[3] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx10; 
  out[4] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dx10; 
  out[5] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dx10; 

  return fabs(1.5*dx10*alpha_max);

} 
