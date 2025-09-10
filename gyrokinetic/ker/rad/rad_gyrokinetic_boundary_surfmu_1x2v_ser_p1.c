#include <gkyl_rad_gyrokinetic_kernels.h> 
GKYL_CU_DH double rad_gyrokinetic_boundary_surfmu_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap,
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

  double rdv2 = 2.0/dxv[2]; 

  if (edge == -1) { 

  double Ghat_r[6] = {0.0}; 
  Ghat_r[0] = ((-0.6123724356957944*(nvsqnu_edge[5]*fedge[11]+nvsqnu_edge[4]*fedge[10]))+0.3535533905932737*(nvsqnu_edge[5]*fedge[9]+nvsqnu_edge[4]*fedge[8])-0.6123724356957944*(nvsqnu_edge[3]*fedge[7]+nvsqnu_edge[2]*fedge[6]+nvsqnu_edge[1]*fedge[5])+0.3535533905932737*nvsqnu_edge[3]*fedge[4]-0.6123724356957944*nvsqnu_edge[0]*fedge[3]+0.3535533905932737*(fedge[2]*nvsqnu_edge[2]+fedge[1]*nvsqnu_edge[1]+fedge[0]*nvsqnu_edge[0]))/vmap_prime_edge[1]; 
  Ghat_r[1] = ((-0.6123724356957944*(nvsqnu_edge[4]*fedge[11]+nvsqnu_edge[5]*fedge[10]))+0.3535533905932737*(nvsqnu_edge[4]*fedge[9]+nvsqnu_edge[5]*fedge[8])-0.6123724356957944*(nvsqnu_edge[2]*fedge[7]+nvsqnu_edge[3]*fedge[6]+nvsqnu_edge[0]*fedge[5])+0.3535533905932737*(nvsqnu_edge[2]*fedge[4]+fedge[2]*nvsqnu_edge[3])-0.6123724356957944*nvsqnu_edge[1]*fedge[3]+0.3535533905932737*(fedge[0]*nvsqnu_edge[1]+nvsqnu_edge[0]*fedge[1]))/vmap_prime_edge[1]; 
  Ghat_r[2] = ((-0.5477225575051661*(nvsqnu_edge[3]*fedge[11]+nvsqnu_edge[2]*fedge[10]))+0.3162277660168379*nvsqnu_edge[3]*fedge[9]+0.3162277660168379*nvsqnu_edge[2]*fedge[8]+((-0.5477225575051661*nvsqnu_edge[5])-0.6123724356957944*nvsqnu_edge[1])*fedge[7]+((-0.5477225575051661*nvsqnu_edge[4])-0.6123724356957944*nvsqnu_edge[0])*fedge[6]+0.3162277660168379*fedge[4]*nvsqnu_edge[5]-0.6123724356957944*nvsqnu_edge[3]*fedge[5]+0.3162277660168379*fedge[2]*nvsqnu_edge[4]+0.3535533905932737*(nvsqnu_edge[1]*fedge[4]+fedge[1]*nvsqnu_edge[3])-0.6123724356957944*nvsqnu_edge[2]*fedge[3]+0.3535533905932737*(fedge[0]*nvsqnu_edge[2]+nvsqnu_edge[0]*fedge[2]))/vmap_prime_edge[1]; 
  Ghat_r[3] = ((-0.5477225575051661*(nvsqnu_edge[2]*fedge[11]+nvsqnu_edge[3]*fedge[10]))+0.3162277660168379*nvsqnu_edge[2]*fedge[9]+0.3162277660168379*nvsqnu_edge[3]*fedge[8]+((-0.5477225575051661*nvsqnu_edge[4])-0.6123724356957944*nvsqnu_edge[0])*fedge[7]+((-0.5477225575051661*nvsqnu_edge[5])-0.6123724356957944*nvsqnu_edge[1])*fedge[6]+0.3162277660168379*fedge[2]*nvsqnu_edge[5]-0.6123724356957944*nvsqnu_edge[2]*fedge[5]+fedge[4]*(0.3162277660168379*nvsqnu_edge[4]+0.3535533905932737*nvsqnu_edge[0])-0.6123724356957944*fedge[3]*nvsqnu_edge[3]+0.3535533905932737*(fedge[0]*nvsqnu_edge[3]+fedge[1]*nvsqnu_edge[2]+nvsqnu_edge[1]*fedge[2]))/vmap_prime_edge[1]; 
  Ghat_r[4] = (((-0.3912303982179757*nvsqnu_edge[5])-0.6123724356957944*nvsqnu_edge[1])*fedge[11]+((-0.3912303982179757*nvsqnu_edge[4])-0.6123724356957944*nvsqnu_edge[0])*fedge[10]+(0.2258769757263128*nvsqnu_edge[5]+0.3535533905932737*nvsqnu_edge[1])*fedge[9]+(0.2258769757263128*nvsqnu_edge[4]+0.3535533905932737*nvsqnu_edge[0])*fedge[8]-0.5477225575051661*(nvsqnu_edge[3]*fedge[7]+nvsqnu_edge[2]*fedge[6])+(0.3535533905932737*fedge[1]-0.6123724356957944*fedge[5])*nvsqnu_edge[5]+(0.3535533905932737*fedge[0]-0.6123724356957944*fedge[3])*nvsqnu_edge[4]+0.3162277660168379*(nvsqnu_edge[3]*fedge[4]+fedge[2]*nvsqnu_edge[2]))/vmap_prime_edge[1]; 
  Ghat_r[5] = (((-0.3912303982179757*nvsqnu_edge[4])-0.6123724356957944*nvsqnu_edge[0])*fedge[11]+((-0.3912303982179757*nvsqnu_edge[5])-0.6123724356957944*nvsqnu_edge[1])*fedge[10]+(0.2258769757263128*nvsqnu_edge[4]+0.3535533905932737*nvsqnu_edge[0])*fedge[9]+(0.2258769757263128*nvsqnu_edge[5]+0.3535533905932737*nvsqnu_edge[1])*fedge[8]-0.5477225575051661*(nvsqnu_edge[2]*fedge[7]+nvsqnu_edge[3]*fedge[6])+(0.3535533905932737*fedge[0]-0.6123724356957944*fedge[3])*nvsqnu_edge[5]+nvsqnu_edge[4]*(0.3535533905932737*fedge[1]-0.6123724356957944*fedge[5])+0.3162277660168379*(nvsqnu_edge[2]*fedge[4]+fedge[2]*nvsqnu_edge[3]))/vmap_prime_edge[1]; 

  out[0] += 0.7071067811865475*Ghat_r[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_r[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat_r[2]*rdv2; 
  out[3] += 1.224744871391589*Ghat_r[0]*rdv2; 
  out[4] += 0.7071067811865475*Ghat_r[3]*rdv2; 
  out[5] += 1.224744871391589*Ghat_r[1]*rdv2; 
  out[6] += 1.224744871391589*Ghat_r[2]*rdv2; 
  out[7] += 1.224744871391589*Ghat_r[3]*rdv2; 
  out[8] += 0.7071067811865475*Ghat_r[4]*rdv2; 
  out[9] += 0.7071067811865475*Ghat_r[5]*rdv2; 
  out[10] += 1.224744871391589*Ghat_r[4]*rdv2; 
  out[11] += 1.224744871391589*Ghat_r[5]*rdv2; 

  } else { 

  double Ghat_l[6] = {0.0}; 
  Ghat_l[0] = ((-0.6123724356957944*(nvsqnu_skin[5]*fskin[11]+nvsqnu_skin[4]*fskin[10]))+0.3535533905932737*(nvsqnu_skin[5]*fskin[9]+nvsqnu_skin[4]*fskin[8])-0.6123724356957944*(nvsqnu_skin[3]*fskin[7]+nvsqnu_skin[2]*fskin[6]+nvsqnu_skin[1]*fskin[5])+0.3535533905932737*nvsqnu_skin[3]*fskin[4]-0.6123724356957944*nvsqnu_skin[0]*fskin[3]+0.3535533905932737*(fskin[2]*nvsqnu_skin[2]+fskin[1]*nvsqnu_skin[1]+fskin[0]*nvsqnu_skin[0]))/vmap_prime_skin[1]; 
  Ghat_l[1] = ((-0.6123724356957944*(nvsqnu_skin[4]*fskin[11]+nvsqnu_skin[5]*fskin[10]))+0.3535533905932737*(nvsqnu_skin[4]*fskin[9]+nvsqnu_skin[5]*fskin[8])-0.6123724356957944*(nvsqnu_skin[2]*fskin[7]+nvsqnu_skin[3]*fskin[6]+nvsqnu_skin[0]*fskin[5])+0.3535533905932737*(nvsqnu_skin[2]*fskin[4]+fskin[2]*nvsqnu_skin[3])-0.6123724356957944*nvsqnu_skin[1]*fskin[3]+0.3535533905932737*(fskin[0]*nvsqnu_skin[1]+nvsqnu_skin[0]*fskin[1]))/vmap_prime_skin[1]; 
  Ghat_l[2] = ((-0.5477225575051661*(nvsqnu_skin[3]*fskin[11]+nvsqnu_skin[2]*fskin[10]))+0.3162277660168379*nvsqnu_skin[3]*fskin[9]+0.3162277660168379*nvsqnu_skin[2]*fskin[8]+((-0.5477225575051661*nvsqnu_skin[5])-0.6123724356957944*nvsqnu_skin[1])*fskin[7]+((-0.5477225575051661*nvsqnu_skin[4])-0.6123724356957944*nvsqnu_skin[0])*fskin[6]+0.3162277660168379*fskin[4]*nvsqnu_skin[5]-0.6123724356957944*nvsqnu_skin[3]*fskin[5]+0.3162277660168379*fskin[2]*nvsqnu_skin[4]+0.3535533905932737*(nvsqnu_skin[1]*fskin[4]+fskin[1]*nvsqnu_skin[3])-0.6123724356957944*nvsqnu_skin[2]*fskin[3]+0.3535533905932737*(fskin[0]*nvsqnu_skin[2]+nvsqnu_skin[0]*fskin[2]))/vmap_prime_skin[1]; 
  Ghat_l[3] = ((-0.5477225575051661*(nvsqnu_skin[2]*fskin[11]+nvsqnu_skin[3]*fskin[10]))+0.3162277660168379*nvsqnu_skin[2]*fskin[9]+0.3162277660168379*nvsqnu_skin[3]*fskin[8]+((-0.5477225575051661*nvsqnu_skin[4])-0.6123724356957944*nvsqnu_skin[0])*fskin[7]+((-0.5477225575051661*nvsqnu_skin[5])-0.6123724356957944*nvsqnu_skin[1])*fskin[6]+0.3162277660168379*fskin[2]*nvsqnu_skin[5]-0.6123724356957944*nvsqnu_skin[2]*fskin[5]+fskin[4]*(0.3162277660168379*nvsqnu_skin[4]+0.3535533905932737*nvsqnu_skin[0])-0.6123724356957944*fskin[3]*nvsqnu_skin[3]+0.3535533905932737*(fskin[0]*nvsqnu_skin[3]+fskin[1]*nvsqnu_skin[2]+nvsqnu_skin[1]*fskin[2]))/vmap_prime_skin[1]; 
  Ghat_l[4] = (((-0.3912303982179757*nvsqnu_skin[5])-0.6123724356957944*nvsqnu_skin[1])*fskin[11]+((-0.3912303982179757*nvsqnu_skin[4])-0.6123724356957944*nvsqnu_skin[0])*fskin[10]+(0.2258769757263128*nvsqnu_skin[5]+0.3535533905932737*nvsqnu_skin[1])*fskin[9]+(0.2258769757263128*nvsqnu_skin[4]+0.3535533905932737*nvsqnu_skin[0])*fskin[8]-0.5477225575051661*(nvsqnu_skin[3]*fskin[7]+nvsqnu_skin[2]*fskin[6])+(0.3535533905932737*fskin[1]-0.6123724356957944*fskin[5])*nvsqnu_skin[5]+(0.3535533905932737*fskin[0]-0.6123724356957944*fskin[3])*nvsqnu_skin[4]+0.3162277660168379*(nvsqnu_skin[3]*fskin[4]+fskin[2]*nvsqnu_skin[2]))/vmap_prime_skin[1]; 
  Ghat_l[5] = (((-0.3912303982179757*nvsqnu_skin[4])-0.6123724356957944*nvsqnu_skin[0])*fskin[11]+((-0.3912303982179757*nvsqnu_skin[5])-0.6123724356957944*nvsqnu_skin[1])*fskin[10]+(0.2258769757263128*nvsqnu_skin[4]+0.3535533905932737*nvsqnu_skin[0])*fskin[9]+(0.2258769757263128*nvsqnu_skin[5]+0.3535533905932737*nvsqnu_skin[1])*fskin[8]-0.5477225575051661*(nvsqnu_skin[2]*fskin[7]+nvsqnu_skin[3]*fskin[6])+(0.3535533905932737*fskin[0]-0.6123724356957944*fskin[3])*nvsqnu_skin[5]+nvsqnu_skin[4]*(0.3535533905932737*fskin[1]-0.6123724356957944*fskin[5])+0.3162277660168379*(nvsqnu_skin[2]*fskin[4]+fskin[2]*nvsqnu_skin[3]))/vmap_prime_skin[1]; 

  out[0] += -0.7071067811865475*Ghat_l[0]*rdv2; 
  out[1] += -0.7071067811865475*Ghat_l[1]*rdv2; 
  out[2] += -0.7071067811865475*Ghat_l[2]*rdv2; 
  out[3] += 1.224744871391589*Ghat_l[0]*rdv2; 
  out[4] += -0.7071067811865475*Ghat_l[3]*rdv2; 
  out[5] += 1.224744871391589*Ghat_l[1]*rdv2; 
  out[6] += 1.224744871391589*Ghat_l[2]*rdv2; 
  out[7] += 1.224744871391589*Ghat_l[3]*rdv2; 
  out[8] += -0.7071067811865475*Ghat_l[4]*rdv2; 
  out[9] += -0.7071067811865475*Ghat_l[5]*rdv2; 
  out[10] += 1.224744871391589*Ghat_l[4]*rdv2; 
  out[11] += 1.224744871391589*Ghat_l[5]*rdv2; 

  } 

  double vmap_prime_min = fmin(fabs(vmap_prime_edge[1]),fabs(vmap_prime_skin[1]));
  double cflFreq = fmax(fabs(nvsqnu_edge[0]/vmap_prime_min), fabs(nvsqnu_skin[0]/vmap_prime_min)); 
  return 1.25*rdv2*cflFreq; 

} 
