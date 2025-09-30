#include <gkyl_nc_hamil_surf_comps_kernels.h> 
GKYL_CU_DH void vlasov_nc_hamil_gen_boundary_surf_node4_x_1x2v_tensor_p1(const double *w, const double *dxv, const int edge,
   const double *poisson_tensor_conf, const double *hamil, const double *fl, 
   const double *fr, double* GKYL_RESTRICT Ghat) 
{ 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[2]; 

  double alpha_quad = 0.0; 
  if (edge == -1) { 
    alpha_quad = (1.224744871391589*poisson_tensor_conf_1[1]+0.7071067811865475*poisson_tensor_conf_1[0])*(2.109478270807616*hamil[8]+1.700840128541523*hamil[7]-0.6324956225107456*hamil[6]-2.888338995926613*hamil[5]-0.5099715653772844*hamil[3]+0.8660254037844386*hamil[2])*dv11+(1.224744871391589*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*((-1.561363675455314*hamil[8])+1.185778957683234*hamil[7]+1.700840128541523*hamil[6]-1.140331086775586*hamil[4]-1.291704467391079*hamil[3]+0.8660254037844386*hamil[1])*dv10; 
  } else { 

    alpha_quad = (0.7071067811865475*poisson_tensor_conf_1[0]-1.224744871391589*poisson_tensor_conf_1[1])*(2.109478270807616*hamil[8]+1.700840128541523*hamil[7]-0.6324956225107456*hamil[6]-2.888338995926613*hamil[5]-0.5099715653772844*hamil[3]+0.8660254037844386*hamil[2])*dv11+(0.7071067811865475*poisson_tensor_conf_0[0]-1.224744871391589*poisson_tensor_conf_0[1])*((-1.561363675455314*hamil[8])+1.185778957683234*hamil[7]+1.700840128541523*hamil[6]-1.140331086775586*hamil[4]-1.291704467391079*hamil[3]+0.8660254037844386*hamil[1])*dv10; 
  } 
 
  f_l_quad = (-0.6123724356957941*fl[17])-0.3535533905932736*fl[16]-0.4937465470290729*fl[15]+0.6670756010385543*fl[14]-0.2850647018386835*fl[13]+0.8384723419661304*fl[12]+0.3851362778294408*fl[11]-0.4472419437481548*fl[10]+0.537852874200477*fl[9]+0.4840922323422013*fl[8]-0.2582152566158886*fl[7]+0.3105295017040593*fl[6]-0.9133729881811893*fl[5]-0.3606043520905964*fl[4]-0.5273361405969424*fl[3]-0.2081950197504564*fl[2]+0.6123724356957944*fl[1]+0.3535533905932737*fl[0]; 
  f_r_quad = 0.6123724356957941*fr[17]-0.3535533905932736*fr[16]+0.4937465470290729*fr[15]-0.6670756010385543*fr[14]-0.2850647018386835*fr[13]-0.8384723419661304*fr[12]+0.3851362778294408*fr[11]+0.4472419437481548*fr[10]-0.537852874200477*fr[9]+0.4840922323422013*fr[8]-0.2582152566158886*fr[7]+0.3105295017040593*fr[6]+0.9133729881811893*fr[5]+0.3606043520905964*fr[4]-0.5273361405969424*fr[3]-0.2081950197504564*fr[2]-0.6123724356957944*fr[1]+0.3535533905932737*fr[0]; 
  double Ghat_nodal = 0.5*(alpha_quad*(f_r_quad + f_l_quad) - fabs(alpha_quad)*(f_r_quad - f_l_quad)); 
  Ghat[0] = 0.1134259259259259*Ghat_nodal; 
  Ghat[1] = -0.06679249447653647*Ghat_nodal; 
  Ghat[2] = -0.1691783804450111*Ghat_nodal; 
  Ghat[3] = 0.09962313244682938*Ghat_nodal; 
  Ghat[4] = -0.08283983508321352*Ghat_nodal; 
  Ghat[5] = 0.1553050010207066*Ghat_nodal; 
  Ghat[6] = 0.1235582519719727*Ghat_nodal; 
  Ghat[7] = -0.09145359262597844*Ghat_nodal; 
  Ghat[8] = -0.1134259259259259*Ghat_nodal; 

} 
