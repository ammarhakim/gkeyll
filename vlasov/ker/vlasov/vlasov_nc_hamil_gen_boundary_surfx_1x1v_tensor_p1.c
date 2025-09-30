#include <gkyl_nc_hamil_surf_comps_kernels.h> 
#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_nc_hamil_gen_boundary_surfx_1x1v_tensor_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, 
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  double Ghat[3]; 

  if (edge == -1) { 

  vlasov_nc_hamil_gen_boundary_surf_node0_x_1x1v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fskin, fedge, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node1_x_1x1v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fskin, fedge, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node2_x_1x1v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fskin, fedge, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node3_x_1x1v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fskin, fedge, Ghat); 
  out[0] += -0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += -0.7071067811865475*Ghat[1]*dx10; 
  out[3] += -1.224744871391589*Ghat[1]*dx10; 
  out[4] += -0.7071067811865475*Ghat[2]*dx10; 
  out[5] += -1.224744871391589*Ghat[2]*dx10; 

  } else { 

  vlasov_nc_hamil_gen_boundary_surf_node0_x_1x1v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fedge, fskin, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node1_x_1x1v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fedge, fskin, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node2_x_1x1v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fedge, fskin, Ghat); 
  vlasov_nc_hamil_gen_boundary_surf_node3_x_1x1v_tensor_p1(w, dxv, edge, poisson_tensor_conf, hamil, fedge, fskin, Ghat); 
  out[0] += 0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += 0.7071067811865475*Ghat[1]*dx10; 
  out[3] += -1.224744871391589*Ghat[1]*dx10; 
  out[4] += 0.7071067811865475*Ghat[2]*dx10; 
  out[5] += -1.224744871391589*Ghat[2]*dx10; 

  } 
  return 0.0;

} 
