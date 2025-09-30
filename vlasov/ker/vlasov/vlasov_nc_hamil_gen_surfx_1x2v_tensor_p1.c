#include <gkyl_nc_hamil_surf_comps_kernels.h> 
#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_nc_hamil_gen_surfx_1x2v_tensor_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 

  double Ghat_r[9] = {0.0}; 
  double Ghat_l[9] = {0.0}; 
  double alpha_max = 0.0; 
  vlasov_nc_hamil_gen_surf_node0_x_1x2v_tensor_p1(w, dxv, poisson_tensor_conf, hamil, fl, fc, fr, &alpha_max, Ghat_l, Ghat_r); 
  vlasov_nc_hamil_gen_surf_node1_x_1x2v_tensor_p1(w, dxv, poisson_tensor_conf, hamil, fl, fc, fr, &alpha_max, Ghat_l, Ghat_r); 
  vlasov_nc_hamil_gen_surf_node2_x_1x2v_tensor_p1(w, dxv, poisson_tensor_conf, hamil, fl, fc, fr, &alpha_max, Ghat_l, Ghat_r); 
  vlasov_nc_hamil_gen_surf_node3_x_1x2v_tensor_p1(w, dxv, poisson_tensor_conf, hamil, fl, fc, fr, &alpha_max, Ghat_l, Ghat_r); 
  vlasov_nc_hamil_gen_surf_node4_x_1x2v_tensor_p1(w, dxv, poisson_tensor_conf, hamil, fl, fc, fr, &alpha_max, Ghat_l, Ghat_r); 
  vlasov_nc_hamil_gen_surf_node5_x_1x2v_tensor_p1(w, dxv, poisson_tensor_conf, hamil, fl, fc, fr, &alpha_max, Ghat_l, Ghat_r); 
  vlasov_nc_hamil_gen_surf_node6_x_1x2v_tensor_p1(w, dxv, poisson_tensor_conf, hamil, fl, fc, fr, &alpha_max, Ghat_l, Ghat_r); 
  vlasov_nc_hamil_gen_surf_node7_x_1x2v_tensor_p1(w, dxv, poisson_tensor_conf, hamil, fl, fc, fr, &alpha_max, Ghat_l, Ghat_r); 
  vlasov_nc_hamil_gen_surf_node8_x_1x2v_tensor_p1(w, dxv, poisson_tensor_conf, hamil, fl, fc, fr, &alpha_max, Ghat_l, Ghat_r); 
  vlasov_nc_hamil_gen_surf_node9_x_1x2v_tensor_p1(w, dxv, poisson_tensor_conf, hamil, fl, fc, fr, &alpha_max, Ghat_l, Ghat_r); 
  vlasov_nc_hamil_gen_surf_node10_x_1x2v_tensor_p1(w, dxv, poisson_tensor_conf, hamil, fl, fc, fr, &alpha_max, Ghat_l, Ghat_r); 
  vlasov_nc_hamil_gen_surf_node11_x_1x2v_tensor_p1(w, dxv, poisson_tensor_conf, hamil, fl, fc, fr, &alpha_max, Ghat_l, Ghat_r); 
  vlasov_nc_hamil_gen_surf_node12_x_1x2v_tensor_p1(w, dxv, poisson_tensor_conf, hamil, fl, fc, fr, &alpha_max, Ghat_l, Ghat_r); 
  vlasov_nc_hamil_gen_surf_node13_x_1x2v_tensor_p1(w, dxv, poisson_tensor_conf, hamil, fl, fc, fr, &alpha_max, Ghat_l, Ghat_r); 
  vlasov_nc_hamil_gen_surf_node14_x_1x2v_tensor_p1(w, dxv, poisson_tensor_conf, hamil, fl, fc, fr, &alpha_max, Ghat_l, Ghat_r); 
  vlasov_nc_hamil_gen_surf_node15_x_1x2v_tensor_p1(w, dxv, poisson_tensor_conf, hamil, fl, fc, fr, &alpha_max, Ghat_l, Ghat_r); 
  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx10; 
  out[1] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx10; 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx10; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dx10; 
  out[4] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx10; 
  out[5] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dx10; 
  out[6] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dx10; 
  out[7] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dx10; 
  out[8] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dx10; 
  out[9] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dx10; 
  out[10] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dx10; 
  out[11] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dx10; 
  out[12] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dx10; 
  out[13] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dx10; 
  out[14] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dx10; 
  out[15] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dx10; 
  out[16] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dx10; 
  out[17] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dx10; 

  return fabs(1.5*dx10*alpha_max);

} 
