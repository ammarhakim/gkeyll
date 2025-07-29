#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_gen_boundary_surfx_1x1v_tensor_p1(const double *w, const double *dxv, 
  const double *hamil, 
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[1]; 

  double Ghat[2]; 

  double Ghat_nodal[2]; 
  if (edge == -1) { 

  double alpha_quad = 0.0; 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  alpha_quad = 1.5*hamil[3]+0.8660254037844386*hamil[2]; 
  f_l_quad = -(0.8660254037844386*fskin[3])-0.5*fskin[2]+0.8660254037844386*fskin[1]+0.5*fskin[0]; 
  f_r_quad = 0.8660254037844386*fedge[3]-0.5*fedge[2]-0.8660254037844386*fedge[1]+0.5*fedge[0]; 
  Ghat_nodal[0] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = 1.5*hamil[3]+0.8660254037844386*hamil[2]; 
  f_l_quad = 0.8660254037844386*fskin[3]+0.5*fskin[2]+0.8660254037844386*fskin[1]+0.5*fskin[0]; 
  f_r_quad = -(0.8660254037844386*fedge[3])+0.5*fedge[2]-0.8660254037844386*fedge[1]+0.5*fedge[0]; 
  Ghat_nodal[1] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  Ghat[0] = 0.7071067811865475*Ghat_nodal[1]+0.7071067811865475*Ghat_nodal[0]; 
  Ghat[1] = 0.7071067811865475*Ghat_nodal[1]-0.7071067811865475*Ghat_nodal[0]; 

  out[0] += -(0.7071067811865475*Ghat[0]*dv10*dx10); 
  out[1] += -(1.224744871391589*Ghat[0]*dv10*dx10); 
  out[2] += -(0.7071067811865475*Ghat[1]*dv10*dx10); 
  out[3] += -(1.224744871391589*Ghat[1]*dv10*dx10); 

  } else { 

  double alpha_quad = 0.0; 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  alpha_quad = 0.8660254037844386*hamil[2]-1.5*hamil[3]; 
  f_l_quad = -(0.8660254037844386*fedge[3])-0.5*fedge[2]+0.8660254037844386*fedge[1]+0.5*fedge[0]; 
  f_r_quad = 0.8660254037844386*fskin[3]-0.5*fskin[2]-0.8660254037844386*fskin[1]+0.5*fskin[0]; 
  Ghat_nodal[0] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  alpha_quad = 0.8660254037844386*hamil[2]-1.5*hamil[3]; 
  f_l_quad = 0.8660254037844386*fedge[3]+0.5*fedge[2]+0.8660254037844386*fedge[1]+0.5*fedge[0]; 
  f_r_quad = -(0.8660254037844386*fskin[3])+0.5*fskin[2]-0.8660254037844386*fskin[1]+0.5*fskin[0]; 
  Ghat_nodal[1] = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 

  Ghat[0] = 0.7071067811865475*Ghat_nodal[1]+0.7071067811865475*Ghat_nodal[0]; 
  Ghat[1] = 0.7071067811865475*Ghat_nodal[1]-0.7071067811865475*Ghat_nodal[0]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv10*dx10; 
  out[1] += -(1.224744871391589*Ghat[0]*dv10*dx10); 
  out[2] += 0.7071067811865475*Ghat[1]*dv10*dx10; 
  out[3] += -(1.224744871391589*Ghat[1]*dv10*dx10); 

  } 
  return 0.0;

} 
