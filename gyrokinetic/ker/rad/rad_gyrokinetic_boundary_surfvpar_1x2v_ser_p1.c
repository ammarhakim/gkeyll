#include <gkyl_rad_gyrokinetic_kernels.h> 
GKYL_CU_DH double rad_gyrokinetic_boundary_surfvpar_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap,
    const double *vmap_prime_edge, const double *vmap_prime_skin,
    const double *nvnu_edge, const double *nvnu_skin, const double *nvsqnu_edge, const double *nvsqnu_skin, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  // w[3]:     cell-center coordinates. 
  // dxv[3]:   cell spacing. 
  // vmap: velocity space mapping.
  // vmap_prime_edge,vmap_prime_skin: velocity space mapping derivative in edge and skin cells.
  // nvnu_edge: Surface expansion sum_s n_s*nu_s(v) in vparallel direction on the lower edges of the edge cell.
  // nvnu_skin: Surface expansion sum_s n_s*nu_s(v) in vparallel direction on the lower edges of the skin cell.
  // nvsqnu_edge: Surface expansion sum_s n_s*nu_s(v) in mu direction on the lower edges of the edge cell.
  // nvsqnu_skin: Surface expansion sum_s n_s*nu_s(v) in mu direction on the lower edges of the edge skin.
  // edge: determines if the update is for the left edge (-1) or right edge (+1).
  // fskin,fedge: distribution function in skin cell/last edge cell.
  // out: output increment in center cell.

  double rdv2 = 2.0/dxv[1]; 

  if (edge == -1) { 

  double Ghat_r[4] = {0.0}; 
  if (0.7071067811865475*vmap[0]>0) {

  Ghat_r[0] = (0.7905694150420947*nvnu_edge[3]*fedge[11]+0.7905694150420948*(nvnu_edge[2]*fedge[10]+nvnu_edge[1]*fedge[9])+0.7905694150420947*nvnu_edge[0]*fedge[8]-0.6123724356957944*(nvnu_edge[3]*fedge[7]+nvnu_edge[2]*fedge[6])+0.3535533905932737*nvnu_edge[3]*fedge[5]-0.6123724356957944*nvnu_edge[1]*fedge[4]+0.3535533905932737*nvnu_edge[2]*fedge[3]-0.6123724356957944*nvnu_edge[0]*fedge[2]+0.3535533905932737*(fedge[1]*nvnu_edge[1]+fedge[0]*nvnu_edge[0]))/vmap_prime_edge[0]; 
  Ghat_r[1] = (0.7905694150420947*nvnu_edge[2]*fedge[11]+0.7905694150420948*(nvnu_edge[3]*fedge[10]+nvnu_edge[0]*fedge[9])+0.7905694150420947*nvnu_edge[1]*fedge[8]-0.6123724356957944*(nvnu_edge[2]*fedge[7]+nvnu_edge[3]*fedge[6])+0.3535533905932737*nvnu_edge[2]*fedge[5]-0.6123724356957944*nvnu_edge[0]*fedge[4]+0.3535533905932737*fedge[3]*nvnu_edge[3]-0.6123724356957944*nvnu_edge[1]*fedge[2]+0.3535533905932737*(fedge[0]*nvnu_edge[1]+nvnu_edge[0]*fedge[1]))/vmap_prime_edge[0]; 
  Ghat_r[2] = (0.7905694150420947*nvnu_edge[1]*fedge[11]+0.7905694150420948*(nvnu_edge[0]*fedge[10]+nvnu_edge[3]*fedge[9])+0.7905694150420947*nvnu_edge[2]*fedge[8]-0.6123724356957944*(nvnu_edge[1]*fedge[7]+nvnu_edge[0]*fedge[6])+0.3535533905932737*nvnu_edge[1]*fedge[5]-0.6123724356957944*nvnu_edge[3]*fedge[4]+0.3535533905932737*(fedge[1]*nvnu_edge[3]+nvnu_edge[0]*fedge[3])+(0.3535533905932737*fedge[0]-0.6123724356957944*fedge[2])*nvnu_edge[2])/vmap_prime_edge[0]; 
  Ghat_r[3] = (0.7905694150420947*nvnu_edge[0]*fedge[11]+0.7905694150420948*(nvnu_edge[1]*fedge[10]+nvnu_edge[2]*fedge[9])+0.7905694150420947*nvnu_edge[3]*fedge[8]-0.6123724356957944*(nvnu_edge[0]*fedge[7]+nvnu_edge[1]*fedge[6])+0.3535533905932737*nvnu_edge[0]*fedge[5]-0.6123724356957944*(nvnu_edge[2]*fedge[4]+fedge[2]*nvnu_edge[3])+0.3535533905932737*(fedge[0]*nvnu_edge[3]+nvnu_edge[1]*fedge[3]+fedge[1]*nvnu_edge[2]))/vmap_prime_edge[0]; 

  } else { 

  Ghat_r[0] = (0.7905694150420947*nvnu_edge[3]*fskin[11]+0.7905694150420948*(nvnu_edge[2]*fskin[10]+nvnu_edge[1]*fskin[9])+0.7905694150420947*nvnu_edge[0]*fskin[8]+0.6123724356957944*(nvnu_edge[3]*fskin[7]+nvnu_edge[2]*fskin[6])+0.3535533905932737*nvnu_edge[3]*fskin[5]+0.6123724356957944*nvnu_edge[1]*fskin[4]+0.3535533905932737*nvnu_edge[2]*fskin[3]+0.6123724356957944*nvnu_edge[0]*fskin[2]+0.3535533905932737*(fskin[1]*nvnu_edge[1]+fskin[0]*nvnu_edge[0]))/vmap_prime_skin[0]; 
  Ghat_r[1] = (0.7905694150420947*nvnu_edge[2]*fskin[11]+0.7905694150420948*(nvnu_edge[3]*fskin[10]+nvnu_edge[0]*fskin[9])+0.7905694150420947*nvnu_edge[1]*fskin[8]+0.6123724356957944*(nvnu_edge[2]*fskin[7]+nvnu_edge[3]*fskin[6])+0.3535533905932737*nvnu_edge[2]*fskin[5]+0.6123724356957944*nvnu_edge[0]*fskin[4]+0.3535533905932737*fskin[3]*nvnu_edge[3]+0.6123724356957944*nvnu_edge[1]*fskin[2]+0.3535533905932737*(fskin[0]*nvnu_edge[1]+nvnu_edge[0]*fskin[1]))/vmap_prime_skin[0]; 
  Ghat_r[2] = (0.7905694150420947*nvnu_edge[1]*fskin[11]+0.7905694150420948*(nvnu_edge[0]*fskin[10]+nvnu_edge[3]*fskin[9])+0.7905694150420947*nvnu_edge[2]*fskin[8]+0.6123724356957944*(nvnu_edge[1]*fskin[7]+nvnu_edge[0]*fskin[6])+0.3535533905932737*nvnu_edge[1]*fskin[5]+0.6123724356957944*nvnu_edge[3]*fskin[4]+0.3535533905932737*(fskin[1]*nvnu_edge[3]+nvnu_edge[0]*fskin[3])+(0.6123724356957944*fskin[2]+0.3535533905932737*fskin[0])*nvnu_edge[2])/vmap_prime_skin[0]; 
  Ghat_r[3] = (0.7905694150420947*nvnu_edge[0]*fskin[11]+0.7905694150420948*(nvnu_edge[1]*fskin[10]+nvnu_edge[2]*fskin[9])+0.7905694150420947*nvnu_edge[3]*fskin[8]+0.6123724356957944*(nvnu_edge[0]*fskin[7]+nvnu_edge[1]*fskin[6])+0.3535533905932737*nvnu_edge[0]*fskin[5]+0.6123724356957944*(nvnu_edge[2]*fskin[4]+fskin[2]*nvnu_edge[3])+0.3535533905932737*(fskin[0]*nvnu_edge[3]+nvnu_edge[1]*fskin[3]+fskin[1]*nvnu_edge[2]))/vmap_prime_skin[0]; 

  } 
  out[0] += 0.7071067811865475*Ghat_r[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_r[1]*rdv2; 
  out[2] += 1.224744871391589*Ghat_r[0]*rdv2; 
  out[3] += 0.7071067811865475*Ghat_r[2]*rdv2; 
  out[4] += 1.224744871391589*Ghat_r[1]*rdv2; 
  out[5] += 0.7071067811865475*Ghat_r[3]*rdv2; 
  out[6] += 1.224744871391589*Ghat_r[2]*rdv2; 
  out[7] += 1.224744871391589*Ghat_r[3]*rdv2; 
  out[8] += 1.58113883008419*Ghat_r[0]*rdv2; 
  out[9] += 1.58113883008419*Ghat_r[1]*rdv2; 
  out[10] += 1.58113883008419*Ghat_r[2]*rdv2; 
  out[11] += 1.58113883008419*Ghat_r[3]*rdv2; 

  } else { 

  double Ghat_l[4] = {0.0}; 
  if (w[1]>0) {

  Ghat_l[0] = (0.7905694150420947*nvnu_skin[3]*fskin[11]+0.7905694150420948*(nvnu_skin[2]*fskin[10]+nvnu_skin[1]*fskin[9])+0.7905694150420947*nvnu_skin[0]*fskin[8]-0.6123724356957944*(nvnu_skin[3]*fskin[7]+nvnu_skin[2]*fskin[6])+0.3535533905932737*nvnu_skin[3]*fskin[5]-0.6123724356957944*nvnu_skin[1]*fskin[4]+0.3535533905932737*nvnu_skin[2]*fskin[3]-0.6123724356957944*nvnu_skin[0]*fskin[2]+0.3535533905932737*(fskin[1]*nvnu_skin[1]+fskin[0]*nvnu_skin[0]))/vmap_prime_skin[0]; 
  Ghat_l[1] = (0.7905694150420947*nvnu_skin[2]*fskin[11]+0.7905694150420948*(nvnu_skin[3]*fskin[10]+nvnu_skin[0]*fskin[9])+0.7905694150420947*nvnu_skin[1]*fskin[8]-0.6123724356957944*(nvnu_skin[2]*fskin[7]+nvnu_skin[3]*fskin[6])+0.3535533905932737*nvnu_skin[2]*fskin[5]-0.6123724356957944*nvnu_skin[0]*fskin[4]+0.3535533905932737*fskin[3]*nvnu_skin[3]-0.6123724356957944*nvnu_skin[1]*fskin[2]+0.3535533905932737*(fskin[0]*nvnu_skin[1]+nvnu_skin[0]*fskin[1]))/vmap_prime_skin[0]; 
  Ghat_l[2] = (0.7905694150420947*nvnu_skin[1]*fskin[11]+0.7905694150420948*(nvnu_skin[0]*fskin[10]+nvnu_skin[3]*fskin[9])+0.7905694150420947*nvnu_skin[2]*fskin[8]-0.6123724356957944*(nvnu_skin[1]*fskin[7]+nvnu_skin[0]*fskin[6])+0.3535533905932737*nvnu_skin[1]*fskin[5]-0.6123724356957944*nvnu_skin[3]*fskin[4]+0.3535533905932737*(fskin[1]*nvnu_skin[3]+nvnu_skin[0]*fskin[3])+(0.3535533905932737*fskin[0]-0.6123724356957944*fskin[2])*nvnu_skin[2])/vmap_prime_skin[0]; 
  Ghat_l[3] = (0.7905694150420947*nvnu_skin[0]*fskin[11]+0.7905694150420948*(nvnu_skin[1]*fskin[10]+nvnu_skin[2]*fskin[9])+0.7905694150420947*nvnu_skin[3]*fskin[8]-0.6123724356957944*(nvnu_skin[0]*fskin[7]+nvnu_skin[1]*fskin[6])+0.3535533905932737*nvnu_skin[0]*fskin[5]-0.6123724356957944*(nvnu_skin[2]*fskin[4]+fskin[2]*nvnu_skin[3])+0.3535533905932737*(fskin[0]*nvnu_skin[3]+nvnu_skin[1]*fskin[3]+fskin[1]*nvnu_skin[2]))/vmap_prime_skin[0]; 

  } else { 

  Ghat_l[0] = (0.7905694150420947*nvnu_skin[3]*fedge[11]+0.7905694150420948*(nvnu_skin[2]*fedge[10]+nvnu_skin[1]*fedge[9])+0.7905694150420947*nvnu_skin[0]*fedge[8]+0.6123724356957944*(nvnu_skin[3]*fedge[7]+nvnu_skin[2]*fedge[6])+0.3535533905932737*nvnu_skin[3]*fedge[5]+0.6123724356957944*nvnu_skin[1]*fedge[4]+0.3535533905932737*nvnu_skin[2]*fedge[3]+0.6123724356957944*nvnu_skin[0]*fedge[2]+0.3535533905932737*(fedge[1]*nvnu_skin[1]+fedge[0]*nvnu_skin[0]))/vmap_prime_edge[0]; 
  Ghat_l[1] = (0.7905694150420947*nvnu_skin[2]*fedge[11]+0.7905694150420948*(nvnu_skin[3]*fedge[10]+nvnu_skin[0]*fedge[9])+0.7905694150420947*nvnu_skin[1]*fedge[8]+0.6123724356957944*(nvnu_skin[2]*fedge[7]+nvnu_skin[3]*fedge[6])+0.3535533905932737*nvnu_skin[2]*fedge[5]+0.6123724356957944*nvnu_skin[0]*fedge[4]+0.3535533905932737*fedge[3]*nvnu_skin[3]+0.6123724356957944*nvnu_skin[1]*fedge[2]+0.3535533905932737*(fedge[0]*nvnu_skin[1]+nvnu_skin[0]*fedge[1]))/vmap_prime_edge[0]; 
  Ghat_l[2] = (0.7905694150420947*nvnu_skin[1]*fedge[11]+0.7905694150420948*(nvnu_skin[0]*fedge[10]+nvnu_skin[3]*fedge[9])+0.7905694150420947*nvnu_skin[2]*fedge[8]+0.6123724356957944*(nvnu_skin[1]*fedge[7]+nvnu_skin[0]*fedge[6])+0.3535533905932737*nvnu_skin[1]*fedge[5]+0.6123724356957944*nvnu_skin[3]*fedge[4]+0.3535533905932737*(fedge[1]*nvnu_skin[3]+nvnu_skin[0]*fedge[3])+(0.6123724356957944*fedge[2]+0.3535533905932737*fedge[0])*nvnu_skin[2])/vmap_prime_edge[0]; 
  Ghat_l[3] = (0.7905694150420947*nvnu_skin[0]*fedge[11]+0.7905694150420948*(nvnu_skin[1]*fedge[10]+nvnu_skin[2]*fedge[9])+0.7905694150420947*nvnu_skin[3]*fedge[8]+0.6123724356957944*(nvnu_skin[0]*fedge[7]+nvnu_skin[1]*fedge[6])+0.3535533905932737*nvnu_skin[0]*fedge[5]+0.6123724356957944*(nvnu_skin[2]*fedge[4]+fedge[2]*nvnu_skin[3])+0.3535533905932737*(fedge[0]*nvnu_skin[3]+nvnu_skin[1]*fedge[3]+fedge[1]*nvnu_skin[2]))/vmap_prime_edge[0]; 

  } 
  out[0] += -0.7071067811865475*Ghat_l[0]*rdv2; 
  out[1] += -0.7071067811865475*Ghat_l[1]*rdv2; 
  out[2] += 1.224744871391589*Ghat_l[0]*rdv2; 
  out[3] += -0.7071067811865475*Ghat_l[2]*rdv2; 
  out[4] += 1.224744871391589*Ghat_l[1]*rdv2; 
  out[5] += -0.7071067811865475*Ghat_l[3]*rdv2; 
  out[6] += 1.224744871391589*Ghat_l[2]*rdv2; 
  out[7] += 1.224744871391589*Ghat_l[3]*rdv2; 
  out[8] += -1.58113883008419*Ghat_l[0]*rdv2; 
  out[9] += -1.58113883008419*Ghat_l[1]*rdv2; 
  out[10] += -1.58113883008419*Ghat_l[2]*rdv2; 
  out[11] += -1.58113883008419*Ghat_l[3]*rdv2; 

  } 

  double vmap_prime_min = fmin(fabs(vmap_prime_edge[0]),fabs(vmap_prime_skin[0]));
  double cflFreq = fmax(fabs(nvnu_edge[0]/vmap_prime_min), fabs(nvnu_skin[0]/vmap_prime_min)); 
  return 1.25*rdv2*cflFreq; 

} 
