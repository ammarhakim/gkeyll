#include <gkyl_rad_gyrokinetic_kernels.h> 
GKYL_CU_DH double rad_gyrokinetic_surfvpar_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap,
    const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r,
    const double *nvnu_l, const double *nvnu_r, const double *nvsqnu_l, const double *nvsqnu_r, 
    const double *fl, const double *fc, const double *fr, 
    double* GKYL_RESTRICT out) 
{ 
  // w[3]:     cell-center coordinates. 
  // dxv[3]:   cell spacing. 
  // vmap: velocity space mapping.
  // vmap_prime_l,vmap_prime_c,vmap_prime_r: velocity space mapping derivative in left, center and right cells.
  // nvnu_l: Surface expansion sum_s n_s*nu_s(v) in vparallel direction on the left.
  // nvnu_r: Surface expansion sum_s n_s*nu_s(v) in vparallel direction on the right.
  // nvsqnu_l: Surface expansion sum_s n_s*nu_s(v) in mu direction on the left.
  // nvsqnu_r: Surface expansion sum_s n_s*nu_s(v) in mu direction on the right.
  // fl/fc/fr:  distribution function in cells 
  // out:       incremented distribution function in cell 

  double rdv2 = 2.0/dxv[1]; 

  double Ghat_r[4] = {0.0}; 
  double Ghat_l[4] = {0.0}; 
  if (0.7071067811865475*vmap[0]>0) {

  Ghat_l[0] = (0.7905694150420947*nvnu_l[3]*fc[11]+0.7905694150420948*(nvnu_l[2]*fc[10]+nvnu_l[1]*fc[9])+0.7905694150420947*nvnu_l[0]*fc[8]-0.6123724356957944*(nvnu_l[3]*fc[7]+nvnu_l[2]*fc[6])+0.3535533905932737*nvnu_l[3]*fc[5]-0.6123724356957944*nvnu_l[1]*fc[4]+0.3535533905932737*nvnu_l[2]*fc[3]-0.6123724356957944*nvnu_l[0]*fc[2]+0.3535533905932737*(fc[1]*nvnu_l[1]+fc[0]*nvnu_l[0]))/vmap_prime_c[0]; 
  Ghat_l[1] = (0.7905694150420947*nvnu_l[2]*fc[11]+0.7905694150420948*(nvnu_l[3]*fc[10]+nvnu_l[0]*fc[9])+0.7905694150420947*nvnu_l[1]*fc[8]-0.6123724356957944*(nvnu_l[2]*fc[7]+nvnu_l[3]*fc[6])+0.3535533905932737*nvnu_l[2]*fc[5]-0.6123724356957944*nvnu_l[0]*fc[4]+0.3535533905932737*fc[3]*nvnu_l[3]-0.6123724356957944*nvnu_l[1]*fc[2]+0.3535533905932737*(fc[0]*nvnu_l[1]+nvnu_l[0]*fc[1]))/vmap_prime_c[0]; 
  Ghat_l[2] = (0.7905694150420947*nvnu_l[1]*fc[11]+0.7905694150420948*(nvnu_l[0]*fc[10]+nvnu_l[3]*fc[9])+0.7905694150420947*nvnu_l[2]*fc[8]-0.6123724356957944*(nvnu_l[1]*fc[7]+nvnu_l[0]*fc[6])+0.3535533905932737*nvnu_l[1]*fc[5]-0.6123724356957944*nvnu_l[3]*fc[4]+0.3535533905932737*(fc[1]*nvnu_l[3]+nvnu_l[0]*fc[3])+(0.3535533905932737*fc[0]-0.6123724356957944*fc[2])*nvnu_l[2])/vmap_prime_c[0]; 
  Ghat_l[3] = (0.7905694150420947*nvnu_l[0]*fc[11]+0.7905694150420948*(nvnu_l[1]*fc[10]+nvnu_l[2]*fc[9])+0.7905694150420947*nvnu_l[3]*fc[8]-0.6123724356957944*(nvnu_l[0]*fc[7]+nvnu_l[1]*fc[6])+0.3535533905932737*nvnu_l[0]*fc[5]-0.6123724356957944*(nvnu_l[2]*fc[4]+fc[2]*nvnu_l[3])+0.3535533905932737*(fc[0]*nvnu_l[3]+nvnu_l[1]*fc[3]+fc[1]*nvnu_l[2]))/vmap_prime_c[0]; 

  Ghat_r[0] = (0.7905694150420947*nvnu_r[3]*fr[11]+0.7905694150420948*(nvnu_r[2]*fr[10]+nvnu_r[1]*fr[9])+0.7905694150420947*nvnu_r[0]*fr[8]-0.6123724356957944*(nvnu_r[3]*fr[7]+nvnu_r[2]*fr[6])+0.3535533905932737*nvnu_r[3]*fr[5]-0.6123724356957944*nvnu_r[1]*fr[4]+0.3535533905932737*nvnu_r[2]*fr[3]-0.6123724356957944*nvnu_r[0]*fr[2]+0.3535533905932737*(fr[1]*nvnu_r[1]+fr[0]*nvnu_r[0]))/vmap_prime_r[0]; 
  Ghat_r[1] = (0.7905694150420947*nvnu_r[2]*fr[11]+0.7905694150420948*(nvnu_r[3]*fr[10]+nvnu_r[0]*fr[9])+0.7905694150420947*nvnu_r[1]*fr[8]-0.6123724356957944*(nvnu_r[2]*fr[7]+nvnu_r[3]*fr[6])+0.3535533905932737*nvnu_r[2]*fr[5]-0.6123724356957944*nvnu_r[0]*fr[4]+0.3535533905932737*fr[3]*nvnu_r[3]-0.6123724356957944*nvnu_r[1]*fr[2]+0.3535533905932737*(fr[0]*nvnu_r[1]+nvnu_r[0]*fr[1]))/vmap_prime_r[0]; 
  Ghat_r[2] = (0.7905694150420947*nvnu_r[1]*fr[11]+0.7905694150420948*(nvnu_r[0]*fr[10]+nvnu_r[3]*fr[9])+0.7905694150420947*nvnu_r[2]*fr[8]-0.6123724356957944*(nvnu_r[1]*fr[7]+nvnu_r[0]*fr[6])+0.3535533905932737*nvnu_r[1]*fr[5]-0.6123724356957944*nvnu_r[3]*fr[4]+0.3535533905932737*(fr[1]*nvnu_r[3]+nvnu_r[0]*fr[3])+(0.3535533905932737*fr[0]-0.6123724356957944*fr[2])*nvnu_r[2])/vmap_prime_r[0]; 
  Ghat_r[3] = (0.7905694150420947*nvnu_r[0]*fr[11]+0.7905694150420948*(nvnu_r[1]*fr[10]+nvnu_r[2]*fr[9])+0.7905694150420947*nvnu_r[3]*fr[8]-0.6123724356957944*(nvnu_r[0]*fr[7]+nvnu_r[1]*fr[6])+0.3535533905932737*nvnu_r[0]*fr[5]-0.6123724356957944*(nvnu_r[2]*fr[4]+fr[2]*nvnu_r[3])+0.3535533905932737*(fr[0]*nvnu_r[3]+nvnu_r[1]*fr[3]+fr[1]*nvnu_r[2]))/vmap_prime_r[0]; 

  } else { 

  Ghat_l[0] = (0.7905694150420947*nvnu_l[3]*fl[11]+0.7905694150420948*(nvnu_l[2]*fl[10]+nvnu_l[1]*fl[9])+0.7905694150420947*nvnu_l[0]*fl[8]+0.6123724356957944*(nvnu_l[3]*fl[7]+nvnu_l[2]*fl[6])+0.3535533905932737*nvnu_l[3]*fl[5]+0.6123724356957944*nvnu_l[1]*fl[4]+0.3535533905932737*nvnu_l[2]*fl[3]+0.6123724356957944*nvnu_l[0]*fl[2]+0.3535533905932737*(fl[1]*nvnu_l[1]+fl[0]*nvnu_l[0]))/vmap_prime_l[0]; 
  Ghat_l[1] = (0.7905694150420947*nvnu_l[2]*fl[11]+0.7905694150420948*(nvnu_l[3]*fl[10]+nvnu_l[0]*fl[9])+0.7905694150420947*nvnu_l[1]*fl[8]+0.6123724356957944*(nvnu_l[2]*fl[7]+nvnu_l[3]*fl[6])+0.3535533905932737*nvnu_l[2]*fl[5]+0.6123724356957944*nvnu_l[0]*fl[4]+0.3535533905932737*fl[3]*nvnu_l[3]+0.6123724356957944*nvnu_l[1]*fl[2]+0.3535533905932737*(fl[0]*nvnu_l[1]+nvnu_l[0]*fl[1]))/vmap_prime_l[0]; 
  Ghat_l[2] = (0.7905694150420947*nvnu_l[1]*fl[11]+0.7905694150420948*(nvnu_l[0]*fl[10]+nvnu_l[3]*fl[9])+0.7905694150420947*nvnu_l[2]*fl[8]+0.6123724356957944*(nvnu_l[1]*fl[7]+nvnu_l[0]*fl[6])+0.3535533905932737*nvnu_l[1]*fl[5]+0.6123724356957944*nvnu_l[3]*fl[4]+0.3535533905932737*(fl[1]*nvnu_l[3]+nvnu_l[0]*fl[3])+(0.6123724356957944*fl[2]+0.3535533905932737*fl[0])*nvnu_l[2])/vmap_prime_l[0]; 
  Ghat_l[3] = (0.7905694150420947*nvnu_l[0]*fl[11]+0.7905694150420948*(nvnu_l[1]*fl[10]+nvnu_l[2]*fl[9])+0.7905694150420947*nvnu_l[3]*fl[8]+0.6123724356957944*(nvnu_l[0]*fl[7]+nvnu_l[1]*fl[6])+0.3535533905932737*nvnu_l[0]*fl[5]+0.6123724356957944*(nvnu_l[2]*fl[4]+fl[2]*nvnu_l[3])+0.3535533905932737*(fl[0]*nvnu_l[3]+nvnu_l[1]*fl[3]+fl[1]*nvnu_l[2]))/vmap_prime_l[0]; 

  Ghat_r[0] = (0.7905694150420947*nvnu_r[3]*fc[11]+0.7905694150420948*(nvnu_r[2]*fc[10]+nvnu_r[1]*fc[9])+0.7905694150420947*nvnu_r[0]*fc[8]+0.6123724356957944*(nvnu_r[3]*fc[7]+nvnu_r[2]*fc[6])+0.3535533905932737*nvnu_r[3]*fc[5]+0.6123724356957944*nvnu_r[1]*fc[4]+0.3535533905932737*nvnu_r[2]*fc[3]+0.6123724356957944*nvnu_r[0]*fc[2]+0.3535533905932737*(fc[1]*nvnu_r[1]+fc[0]*nvnu_r[0]))/vmap_prime_c[0]; 
  Ghat_r[1] = (0.7905694150420947*nvnu_r[2]*fc[11]+0.7905694150420948*(nvnu_r[3]*fc[10]+nvnu_r[0]*fc[9])+0.7905694150420947*nvnu_r[1]*fc[8]+0.6123724356957944*(nvnu_r[2]*fc[7]+nvnu_r[3]*fc[6])+0.3535533905932737*nvnu_r[2]*fc[5]+0.6123724356957944*nvnu_r[0]*fc[4]+0.3535533905932737*fc[3]*nvnu_r[3]+0.6123724356957944*nvnu_r[1]*fc[2]+0.3535533905932737*(fc[0]*nvnu_r[1]+nvnu_r[0]*fc[1]))/vmap_prime_c[0]; 
  Ghat_r[2] = (0.7905694150420947*nvnu_r[1]*fc[11]+0.7905694150420948*(nvnu_r[0]*fc[10]+nvnu_r[3]*fc[9])+0.7905694150420947*nvnu_r[2]*fc[8]+0.6123724356957944*(nvnu_r[1]*fc[7]+nvnu_r[0]*fc[6])+0.3535533905932737*nvnu_r[1]*fc[5]+0.6123724356957944*nvnu_r[3]*fc[4]+0.3535533905932737*(fc[1]*nvnu_r[3]+nvnu_r[0]*fc[3])+(0.6123724356957944*fc[2]+0.3535533905932737*fc[0])*nvnu_r[2])/vmap_prime_c[0]; 
  Ghat_r[3] = (0.7905694150420947*nvnu_r[0]*fc[11]+0.7905694150420948*(nvnu_r[1]*fc[10]+nvnu_r[2]*fc[9])+0.7905694150420947*nvnu_r[3]*fc[8]+0.6123724356957944*(nvnu_r[0]*fc[7]+nvnu_r[1]*fc[6])+0.3535533905932737*nvnu_r[0]*fc[5]+0.6123724356957944*(nvnu_r[2]*fc[4]+fc[2]*nvnu_r[3])+0.3535533905932737*(fc[0]*nvnu_r[3]+nvnu_r[1]*fc[3]+fc[1]*nvnu_r[2]))/vmap_prime_c[0]; 

  } 
  out[0] += (0.7071067811865475*Ghat_r[0]-0.7071067811865475*Ghat_l[0])*rdv2; 
  out[1] += (0.7071067811865475*Ghat_r[1]-0.7071067811865475*Ghat_l[1])*rdv2; 
  out[2] += 1.224744871391589*(Ghat_r[0]+Ghat_l[0])*rdv2; 
  out[3] += (0.7071067811865475*Ghat_r[2]-0.7071067811865475*Ghat_l[2])*rdv2; 
  out[4] += 1.224744871391589*(Ghat_r[1]+Ghat_l[1])*rdv2; 
  out[5] += (0.7071067811865475*Ghat_r[3]-0.7071067811865475*Ghat_l[3])*rdv2; 
  out[6] += 1.224744871391589*(Ghat_r[2]+Ghat_l[2])*rdv2; 
  out[7] += 1.224744871391589*(Ghat_r[3]+Ghat_l[3])*rdv2; 
  out[8] += (1.58113883008419*Ghat_r[0]-1.58113883008419*Ghat_l[0])*rdv2; 
  out[9] += (1.58113883008419*Ghat_r[1]-1.58113883008419*Ghat_l[1])*rdv2; 
  out[10] += (1.58113883008419*Ghat_r[2]-1.58113883008419*Ghat_l[2])*rdv2; 
  out[11] += (1.58113883008419*Ghat_r[3]-1.58113883008419*Ghat_l[3])*rdv2; 

  double vmap_prime_min = fmin(fmin(fabs(vmap_prime_l[0]),fabs(vmap_prime_c[0])),fabs(vmap_prime_r[0]));
  double cflFreq = fmax(fabs(nvnu_l[0]/vmap_prime_min), fabs(nvnu_r[0]/vmap_prime_min)); 
  return 1.25*rdv2*cflFreq; 

} 
