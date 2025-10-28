#include <gkyl_rad_gyrokinetic_kernels.h> 
GKYL_CU_DH double rad_gyrokinetic_surfmu_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap,
    const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r,
    const double *nvnu_l, const double *nvnu_r, const double *nvsqnu_l, const double *nvsqnu_r, 
    const double *fl, const double *fc, const double *fr, 
    double* GKYL_RESTRICT out) 
{ 
  // w[3]:     cell-center coordinates. 
  // dxv[3]:   cell spacing. 
  // vmap_prime: velocity space mapping.
  // vmap_prime_l,vmap_prime_c,vmap_prime_r: velocity space mapping derivative in left, center and right cells.
  // nvnu_l: Surface expansion sum_s n_s*nu_s(v) in vparallel direction on the left.
  // nvnu_r: Surface expansion sum_s n_s*nu_s(v) in vparallel direction on the right.
  // nvsqnu_l: Surface expansion sum_s n_s*nu_s(v) in mu direction on the left.
  // nvsqnu_r: Surface expansion sum_s n_s*nu_s(v) in mu direction on the right.
  // fl/fc/fr:  distribution function in cells 
  // out:       incremented distribution function in cell 

  double rdv2 = 2.0/dxv[2]; 

  double Ghat_r[6] = {0.0}; 
  double Ghat_l[6] = {0.0}; 
  Ghat_l[0] = ((-0.6123724356957944*(nvsqnu_l[5]*fc[11]+nvsqnu_l[4]*fc[10]))+0.3535533905932737*(nvsqnu_l[5]*fc[9]+nvsqnu_l[4]*fc[8])-0.6123724356957944*(nvsqnu_l[3]*fc[7]+nvsqnu_l[2]*fc[6]+nvsqnu_l[1]*fc[5])+0.3535533905932737*nvsqnu_l[3]*fc[4]-0.6123724356957944*nvsqnu_l[0]*fc[3]+0.3535533905932737*(fc[2]*nvsqnu_l[2]+fc[1]*nvsqnu_l[1]+fc[0]*nvsqnu_l[0]))/vmap_prime_c[1]; 
  Ghat_l[1] = ((-0.6123724356957944*(nvsqnu_l[4]*fc[11]+nvsqnu_l[5]*fc[10]))+0.3535533905932737*(nvsqnu_l[4]*fc[9]+nvsqnu_l[5]*fc[8])-0.6123724356957944*(nvsqnu_l[2]*fc[7]+nvsqnu_l[3]*fc[6]+nvsqnu_l[0]*fc[5])+0.3535533905932737*(nvsqnu_l[2]*fc[4]+fc[2]*nvsqnu_l[3])-0.6123724356957944*nvsqnu_l[1]*fc[3]+0.3535533905932737*(fc[0]*nvsqnu_l[1]+nvsqnu_l[0]*fc[1]))/vmap_prime_c[1]; 
  Ghat_l[2] = ((-0.5477225575051661*(nvsqnu_l[3]*fc[11]+nvsqnu_l[2]*fc[10]))+0.3162277660168379*nvsqnu_l[3]*fc[9]+0.3162277660168379*nvsqnu_l[2]*fc[8]+((-0.5477225575051661*nvsqnu_l[5])-0.6123724356957944*nvsqnu_l[1])*fc[7]+((-0.5477225575051661*nvsqnu_l[4])-0.6123724356957944*nvsqnu_l[0])*fc[6]+0.3162277660168379*fc[4]*nvsqnu_l[5]-0.6123724356957944*nvsqnu_l[3]*fc[5]+0.3162277660168379*fc[2]*nvsqnu_l[4]+0.3535533905932737*(nvsqnu_l[1]*fc[4]+fc[1]*nvsqnu_l[3])-0.6123724356957944*nvsqnu_l[2]*fc[3]+0.3535533905932737*(fc[0]*nvsqnu_l[2]+nvsqnu_l[0]*fc[2]))/vmap_prime_c[1]; 
  Ghat_l[3] = ((-0.5477225575051661*(nvsqnu_l[2]*fc[11]+nvsqnu_l[3]*fc[10]))+0.3162277660168379*nvsqnu_l[2]*fc[9]+0.3162277660168379*nvsqnu_l[3]*fc[8]+((-0.5477225575051661*nvsqnu_l[4])-0.6123724356957944*nvsqnu_l[0])*fc[7]+((-0.5477225575051661*nvsqnu_l[5])-0.6123724356957944*nvsqnu_l[1])*fc[6]+0.3162277660168379*fc[2]*nvsqnu_l[5]-0.6123724356957944*nvsqnu_l[2]*fc[5]+fc[4]*(0.3162277660168379*nvsqnu_l[4]+0.3535533905932737*nvsqnu_l[0])-0.6123724356957944*fc[3]*nvsqnu_l[3]+0.3535533905932737*(fc[0]*nvsqnu_l[3]+fc[1]*nvsqnu_l[2]+nvsqnu_l[1]*fc[2]))/vmap_prime_c[1]; 
  Ghat_l[4] = (((-0.3912303982179757*nvsqnu_l[5])-0.6123724356957944*nvsqnu_l[1])*fc[11]+((-0.3912303982179757*nvsqnu_l[4])-0.6123724356957944*nvsqnu_l[0])*fc[10]+(0.2258769757263128*nvsqnu_l[5]+0.3535533905932737*nvsqnu_l[1])*fc[9]+(0.2258769757263128*nvsqnu_l[4]+0.3535533905932737*nvsqnu_l[0])*fc[8]-0.5477225575051661*(nvsqnu_l[3]*fc[7]+nvsqnu_l[2]*fc[6])+(0.3535533905932737*fc[1]-0.6123724356957944*fc[5])*nvsqnu_l[5]+(0.3535533905932737*fc[0]-0.6123724356957944*fc[3])*nvsqnu_l[4]+0.3162277660168379*(nvsqnu_l[3]*fc[4]+fc[2]*nvsqnu_l[2]))/vmap_prime_c[1]; 
  Ghat_l[5] = (((-0.3912303982179757*nvsqnu_l[4])-0.6123724356957944*nvsqnu_l[0])*fc[11]+((-0.3912303982179757*nvsqnu_l[5])-0.6123724356957944*nvsqnu_l[1])*fc[10]+(0.2258769757263128*nvsqnu_l[4]+0.3535533905932737*nvsqnu_l[0])*fc[9]+(0.2258769757263128*nvsqnu_l[5]+0.3535533905932737*nvsqnu_l[1])*fc[8]-0.5477225575051661*(nvsqnu_l[2]*fc[7]+nvsqnu_l[3]*fc[6])+(0.3535533905932737*fc[0]-0.6123724356957944*fc[3])*nvsqnu_l[5]+nvsqnu_l[4]*(0.3535533905932737*fc[1]-0.6123724356957944*fc[5])+0.3162277660168379*(nvsqnu_l[2]*fc[4]+fc[2]*nvsqnu_l[3]))/vmap_prime_c[1]; 

  Ghat_r[0] = ((-0.6123724356957944*(nvsqnu_r[5]*fr[11]+nvsqnu_r[4]*fr[10]))+0.3535533905932737*(nvsqnu_r[5]*fr[9]+nvsqnu_r[4]*fr[8])-0.6123724356957944*(nvsqnu_r[3]*fr[7]+nvsqnu_r[2]*fr[6]+nvsqnu_r[1]*fr[5])+0.3535533905932737*nvsqnu_r[3]*fr[4]-0.6123724356957944*nvsqnu_r[0]*fr[3]+0.3535533905932737*(fr[2]*nvsqnu_r[2]+fr[1]*nvsqnu_r[1]+fr[0]*nvsqnu_r[0]))/vmap_prime_r[1]; 
  Ghat_r[1] = ((-0.6123724356957944*(nvsqnu_r[4]*fr[11]+nvsqnu_r[5]*fr[10]))+0.3535533905932737*(nvsqnu_r[4]*fr[9]+nvsqnu_r[5]*fr[8])-0.6123724356957944*(nvsqnu_r[2]*fr[7]+nvsqnu_r[3]*fr[6]+nvsqnu_r[0]*fr[5])+0.3535533905932737*(nvsqnu_r[2]*fr[4]+fr[2]*nvsqnu_r[3])-0.6123724356957944*nvsqnu_r[1]*fr[3]+0.3535533905932737*(fr[0]*nvsqnu_r[1]+nvsqnu_r[0]*fr[1]))/vmap_prime_r[1]; 
  Ghat_r[2] = ((-0.5477225575051661*(nvsqnu_r[3]*fr[11]+nvsqnu_r[2]*fr[10]))+0.3162277660168379*nvsqnu_r[3]*fr[9]+0.3162277660168379*nvsqnu_r[2]*fr[8]+((-0.5477225575051661*nvsqnu_r[5])-0.6123724356957944*nvsqnu_r[1])*fr[7]+((-0.5477225575051661*nvsqnu_r[4])-0.6123724356957944*nvsqnu_r[0])*fr[6]+0.3162277660168379*fr[4]*nvsqnu_r[5]-0.6123724356957944*nvsqnu_r[3]*fr[5]+0.3162277660168379*fr[2]*nvsqnu_r[4]+0.3535533905932737*(nvsqnu_r[1]*fr[4]+fr[1]*nvsqnu_r[3])-0.6123724356957944*nvsqnu_r[2]*fr[3]+0.3535533905932737*(fr[0]*nvsqnu_r[2]+nvsqnu_r[0]*fr[2]))/vmap_prime_r[1]; 
  Ghat_r[3] = ((-0.5477225575051661*(nvsqnu_r[2]*fr[11]+nvsqnu_r[3]*fr[10]))+0.3162277660168379*nvsqnu_r[2]*fr[9]+0.3162277660168379*nvsqnu_r[3]*fr[8]+((-0.5477225575051661*nvsqnu_r[4])-0.6123724356957944*nvsqnu_r[0])*fr[7]+((-0.5477225575051661*nvsqnu_r[5])-0.6123724356957944*nvsqnu_r[1])*fr[6]+0.3162277660168379*fr[2]*nvsqnu_r[5]-0.6123724356957944*nvsqnu_r[2]*fr[5]+fr[4]*(0.3162277660168379*nvsqnu_r[4]+0.3535533905932737*nvsqnu_r[0])-0.6123724356957944*fr[3]*nvsqnu_r[3]+0.3535533905932737*(fr[0]*nvsqnu_r[3]+fr[1]*nvsqnu_r[2]+nvsqnu_r[1]*fr[2]))/vmap_prime_r[1]; 
  Ghat_r[4] = (((-0.3912303982179757*nvsqnu_r[5])-0.6123724356957944*nvsqnu_r[1])*fr[11]+((-0.3912303982179757*nvsqnu_r[4])-0.6123724356957944*nvsqnu_r[0])*fr[10]+(0.2258769757263128*nvsqnu_r[5]+0.3535533905932737*nvsqnu_r[1])*fr[9]+(0.2258769757263128*nvsqnu_r[4]+0.3535533905932737*nvsqnu_r[0])*fr[8]-0.5477225575051661*(nvsqnu_r[3]*fr[7]+nvsqnu_r[2]*fr[6])+(0.3535533905932737*fr[1]-0.6123724356957944*fr[5])*nvsqnu_r[5]+(0.3535533905932737*fr[0]-0.6123724356957944*fr[3])*nvsqnu_r[4]+0.3162277660168379*(nvsqnu_r[3]*fr[4]+fr[2]*nvsqnu_r[2]))/vmap_prime_r[1]; 
  Ghat_r[5] = (((-0.3912303982179757*nvsqnu_r[4])-0.6123724356957944*nvsqnu_r[0])*fr[11]+((-0.3912303982179757*nvsqnu_r[5])-0.6123724356957944*nvsqnu_r[1])*fr[10]+(0.2258769757263128*nvsqnu_r[4]+0.3535533905932737*nvsqnu_r[0])*fr[9]+(0.2258769757263128*nvsqnu_r[5]+0.3535533905932737*nvsqnu_r[1])*fr[8]-0.5477225575051661*(nvsqnu_r[2]*fr[7]+nvsqnu_r[3]*fr[6])+(0.3535533905932737*fr[0]-0.6123724356957944*fr[3])*nvsqnu_r[5]+nvsqnu_r[4]*(0.3535533905932737*fr[1]-0.6123724356957944*fr[5])+0.3162277660168379*(nvsqnu_r[2]*fr[4]+fr[2]*nvsqnu_r[3]))/vmap_prime_r[1]; 

  out[0] += (0.7071067811865475*Ghat_r[0]-0.7071067811865475*Ghat_l[0])*rdv2; 
  out[1] += (0.7071067811865475*Ghat_r[1]-0.7071067811865475*Ghat_l[1])*rdv2; 
  out[2] += (0.7071067811865475*Ghat_r[2]-0.7071067811865475*Ghat_l[2])*rdv2; 
  out[3] += 1.224744871391589*(Ghat_r[0]+Ghat_l[0])*rdv2; 
  out[4] += (0.7071067811865475*Ghat_r[3]-0.7071067811865475*Ghat_l[3])*rdv2; 
  out[5] += 1.224744871391589*(Ghat_r[1]+Ghat_l[1])*rdv2; 
  out[6] += 1.224744871391589*(Ghat_r[2]+Ghat_l[2])*rdv2; 
  out[7] += 1.224744871391589*(Ghat_r[3]+Ghat_l[3])*rdv2; 
  out[8] += (0.7071067811865475*Ghat_r[4]-0.7071067811865475*Ghat_l[4])*rdv2; 
  out[9] += (0.7071067811865475*Ghat_r[5]-0.7071067811865475*Ghat_l[5])*rdv2; 
  out[10] += 1.224744871391589*(Ghat_r[4]+Ghat_l[4])*rdv2; 
  out[11] += 1.224744871391589*(Ghat_r[5]+Ghat_l[5])*rdv2; 

  double vmap_prime_min = fmin(fmin(fabs(vmap_prime_l[1]),fabs(vmap_prime_c[1])),fabs(vmap_prime_r[1]));
  double cflFreq = fmax(fabs(nvsqnu_l[0]/vmap_prime_min), fabs(nvsqnu_r[0]/vmap_prime_min)); 
  return 0.75*rdv2*cflFreq; 

} 
