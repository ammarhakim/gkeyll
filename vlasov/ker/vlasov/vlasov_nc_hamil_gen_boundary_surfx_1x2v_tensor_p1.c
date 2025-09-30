#include <gkyl_nc_hamil_surf_comps_kernels.h> 
#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_nc_hamil_gen_boundary_surfx_1x2v_tensor_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  double Ghat[9]; 

  if (edge == -1) { 

  vlasov_nc_hamil_gen_boundary_surf_node0_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fskin, fedge, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node1_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fskin, fedge, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node2_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fskin, fedge, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node3_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fskin, fedge, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node4_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fskin, fedge, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node5_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fskin, fedge, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node6_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fskin, fedge, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node7_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fskin, fedge, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node8_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fskin, fedge, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node9_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fskin, fedge, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node10_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fskin, fedge, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node11_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fskin, fedge, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node12_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fskin, fedge, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node13_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fskin, fedge, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node14_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fskin, fedge, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node15_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fskin, fedge, Ghat); 
  out[0] += -0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += -0.7071067811865475*Ghat[1]*dx10; 
  out[3] += -0.7071067811865475*Ghat[2]*dx10; 
  out[4] += -1.224744871391589*Ghat[1]*dx10; 
  out[5] += -1.224744871391589*Ghat[2]*dx10; 
  out[6] += -0.7071067811865475*Ghat[3]*dx10; 
  out[7] += -0.7071067811865475*Ghat[4]*dx10; 
  out[8] += -0.7071067811865475*Ghat[5]*dx10; 
  out[9] += -1.224744871391589*Ghat[3]*dx10; 
  out[10] += -1.224744871391589*Ghat[4]*dx10; 
  out[11] += -0.7071067811865475*Ghat[6]*dx10; 
  out[12] += -1.224744871391589*Ghat[5]*dx10; 
  out[13] += -0.7071067811865475*Ghat[7]*dx10; 
  out[14] += -1.224744871391589*Ghat[6]*dx10; 
  out[15] += -1.224744871391589*Ghat[7]*dx10; 
  out[16] += -0.7071067811865475*Ghat[8]*dx10; 
  out[17] += -1.224744871391589*Ghat[8]*dx10; 

  } else { 

  vlasov_nc_hamil_gen_boundary_surf_node0_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fedge, fskin, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node1_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fedge, fskin, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node2_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fedge, fskin, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node3_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fedge, fskin, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node4_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fedge, fskin, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node5_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fedge, fskin, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node6_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fedge, fskin, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node7_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fedge, fskin, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node8_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fedge, fskin, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node9_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fedge, fskin, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node10_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fedge, fskin, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node11_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fedge, fskin, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node12_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fedge, fskin, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node13_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fedge, fskin, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node14_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fedge, fskin, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node15_x_1x2v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fedge, fskin, Ghat); 
  out[0] += 0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += 0.7071067811865475*Ghat[1]*dx10; 
  out[3] += 0.7071067811865475*Ghat[2]*dx10; 
  out[4] += -1.224744871391589*Ghat[1]*dx10; 
  out[5] += -1.224744871391589*Ghat[2]*dx10; 
  out[6] += 0.7071067811865475*Ghat[3]*dx10; 
  out[7] += 0.7071067811865475*Ghat[4]*dx10; 
  out[8] += 0.7071067811865475*Ghat[5]*dx10; 
  out[9] += -1.224744871391589*Ghat[3]*dx10; 
  out[10] += -1.224744871391589*Ghat[4]*dx10; 
  out[11] += 0.7071067811865475*Ghat[6]*dx10; 
  out[12] += -1.224744871391589*Ghat[5]*dx10; 
  out[13] += 0.7071067811865475*Ghat[7]*dx10; 
  out[14] += -1.224744871391589*Ghat[6]*dx10; 
  out[15] += -1.224744871391589*Ghat[7]*dx10; 
  out[16] += 0.7071067811865475*Ghat[8]*dx10; 
  out[17] += -1.224744871391589*Ghat[8]*dx10; 

  } 
  return 0.0;

} 
