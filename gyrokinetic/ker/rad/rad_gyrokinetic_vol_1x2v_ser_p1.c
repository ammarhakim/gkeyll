#include <gkyl_rad_gyrokinetic_kernels.h> 
#include <gkyl_rad_gyrokinetic_kernels.h> 
GKYL_CU_DH double rad_gyrokinetic_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap_prime, 
    const double *nvnu, const double *nvsqnu, 
    const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[3]: cell-center coordinates. 
  // dxv[3]: cell spacing. 
  // vmap_prime: velocity space mapping derivative.
  // nvnu: Volume expansion of sum_s n_s*nu_s(v) in vparallel direction.
  // nvsqnu: Volume expansion of sum_s n_s*nu_s(v) in mu direction.
  // f: input distribution function.
  // out: incremented output 

  double rdv2[2]; 
  rdv2[0] = 2.0/dxv[1]; 
  rdv2[1] = 2.0/dxv[2]; 

  double nvnu_f[12] = {0.0}; 
  double nvsqnu_f[12] = {0.0}; 
  nvnu_f[0] = (0.3535533905932737*(f[11]*nvnu[11]+f[10]*nvnu[10]+f[9]*nvnu[9]+f[8]*nvnu[8]+f[7]*nvnu[7]+f[6]*nvnu[6]+f[5]*nvnu[5]+f[4]*nvnu[4]+f[3]*nvnu[3]+f[2]*nvnu[2]+f[1]*nvnu[1]+f[0]*nvnu[0]))/vmap_prime[0]; 
  nvnu_f[1] = (0.3535533905932737*(f[10]*nvnu[11]+nvnu[10]*f[11]+f[8]*nvnu[9]+nvnu[8]*f[9]+f[6]*nvnu[7]+nvnu[6]*f[7]+f[3]*nvnu[5]+nvnu[3]*f[5]+f[2]*nvnu[4]+nvnu[2]*f[4]+f[0]*nvnu[1]+nvnu[0]*f[1]))/vmap_prime[0]; 
  nvnu_f[2] = (0.3162277660168379*(f[7]*nvnu[11]+nvnu[7]*f[11])+0.3162277660168379*(f[6]*nvnu[10]+nvnu[6]*f[10]+f[4]*nvnu[9]+nvnu[4]*f[9])+0.3162277660168379*(f[2]*nvnu[8]+nvnu[2]*f[8])+0.3535533905932737*(f[5]*nvnu[7]+nvnu[5]*f[7]+f[3]*nvnu[6]+nvnu[3]*f[6]+f[1]*nvnu[4]+nvnu[1]*f[4]+f[0]*nvnu[2]+nvnu[0]*f[2]))/vmap_prime[0]; 
  nvnu_f[3] = (0.3535533905932737*(f[9]*nvnu[11]+nvnu[9]*f[11]+f[8]*nvnu[10]+nvnu[8]*f[10]+f[4]*nvnu[7]+nvnu[4]*f[7]+f[2]*nvnu[6]+nvnu[2]*f[6]+f[1]*nvnu[5]+nvnu[1]*f[5]+f[0]*nvnu[3]+nvnu[0]*f[3]))/vmap_prime[0]; 
  nvnu_f[4] = (0.3162277660168379*(f[6]*nvnu[11]+nvnu[6]*f[11])+0.3162277660168379*(f[7]*nvnu[10]+nvnu[7]*f[10]+f[2]*nvnu[9]+nvnu[2]*f[9])+0.3162277660168379*(f[4]*nvnu[8]+nvnu[4]*f[8])+0.3535533905932737*(f[3]*nvnu[7]+nvnu[3]*f[7]+f[5]*nvnu[6]+nvnu[5]*f[6]+f[0]*nvnu[4]+nvnu[0]*f[4]+f[1]*nvnu[2]+nvnu[1]*f[2]))/vmap_prime[0]; 
  nvnu_f[5] = (0.3535533905932737*(f[8]*nvnu[11]+nvnu[8]*f[11]+f[9]*nvnu[10]+nvnu[9]*f[10]+f[2]*nvnu[7]+nvnu[2]*f[7]+f[4]*nvnu[6]+nvnu[4]*f[6]+f[0]*nvnu[5]+nvnu[0]*f[5]+f[1]*nvnu[3]+nvnu[1]*f[3]))/vmap_prime[0]; 
  nvnu_f[6] = (0.3162277660168379*(f[4]*nvnu[11]+nvnu[4]*f[11])+0.3162277660168379*(f[2]*nvnu[10]+nvnu[2]*f[10]+f[7]*nvnu[9]+nvnu[7]*f[9])+0.3162277660168379*(f[6]*nvnu[8]+nvnu[6]*f[8])+0.3535533905932737*(f[1]*nvnu[7]+nvnu[1]*f[7]+f[0]*nvnu[6]+nvnu[0]*f[6]+f[4]*nvnu[5]+nvnu[4]*f[5]+f[2]*nvnu[3]+nvnu[2]*f[3]))/vmap_prime[0]; 
  nvnu_f[7] = (0.3162277660168379*(f[2]*nvnu[11]+nvnu[2]*f[11])+0.3162277660168379*(f[4]*nvnu[10]+nvnu[4]*f[10]+f[6]*nvnu[9]+nvnu[6]*f[9])+0.3162277660168379*(f[7]*nvnu[8]+nvnu[7]*f[8])+0.3535533905932737*(f[0]*nvnu[7]+nvnu[0]*f[7]+f[1]*nvnu[6]+nvnu[1]*f[6]+f[2]*nvnu[5]+nvnu[2]*f[5]+f[3]*nvnu[4]+nvnu[3]*f[4]))/vmap_prime[0]; 
  nvnu_f[8] = (0.2258769757263128*f[11]*nvnu[11]+0.3535533905932737*(f[5]*nvnu[11]+nvnu[5]*f[11])+0.2258769757263128*f[10]*nvnu[10]+0.3535533905932737*(f[3]*nvnu[10]+nvnu[3]*f[10])+0.2258769757263128*f[9]*nvnu[9]+0.3535533905932737*(f[1]*nvnu[9]+nvnu[1]*f[9])+0.2258769757263128*f[8]*nvnu[8]+0.3535533905932737*(f[0]*nvnu[8]+nvnu[0]*f[8])+0.3162277660168379*(f[7]*nvnu[7]+f[6]*nvnu[6]+f[4]*nvnu[4]+f[2]*nvnu[2]))/vmap_prime[0]; 
  nvnu_f[9] = ((0.2258769757263128*f[10]+0.3535533905932737*f[3])*nvnu[11]+0.2258769757263128*nvnu[10]*f[11]+0.3535533905932737*(nvnu[3]*f[11]+f[5]*nvnu[10]+nvnu[5]*f[10])+(0.2258769757263128*f[8]+0.3535533905932737*f[0])*nvnu[9]+0.2258769757263128*nvnu[8]*f[9]+0.3535533905932737*(nvnu[0]*f[9]+f[1]*nvnu[8]+nvnu[1]*f[8])+0.3162277660168379*(f[6]*nvnu[7]+nvnu[6]*f[7]+f[2]*nvnu[4]+nvnu[2]*f[4]))/vmap_prime[0]; 
  nvnu_f[10] = ((0.2258769757263128*f[9]+0.3535533905932737*f[1])*nvnu[11]+(0.2258769757263128*nvnu[9]+0.3535533905932737*nvnu[1])*f[11]+(0.2258769757263128*f[8]+0.3535533905932737*f[0])*nvnu[10]+0.2258769757263128*nvnu[8]*f[10]+0.3535533905932737*(nvnu[0]*f[10]+f[5]*nvnu[9]+nvnu[5]*f[9]+f[3]*nvnu[8]+nvnu[3]*f[8])+0.3162277660168379*(f[4]*nvnu[7]+nvnu[4]*f[7]+f[2]*nvnu[6]+nvnu[2]*f[6]))/vmap_prime[0]; 
  nvnu_f[11] = ((0.2258769757263128*f[8]+0.3535533905932737*f[0])*nvnu[11]+(0.2258769757263128*nvnu[8]+0.3535533905932737*nvnu[0])*f[11]+(0.2258769757263128*f[9]+0.3535533905932737*f[1])*nvnu[10]+0.2258769757263128*nvnu[9]*f[10]+0.3535533905932737*(nvnu[1]*f[10]+f[3]*nvnu[9]+nvnu[3]*f[9]+f[5]*nvnu[8]+nvnu[5]*f[8])+0.3162277660168379*(f[2]*nvnu[7]+nvnu[2]*f[7]+f[4]*nvnu[6]+nvnu[4]*f[6]))/vmap_prime[0]; 

  nvsqnu_f[0] = (0.3535533905932737*(f[11]*nvsqnu[11]+f[10]*nvsqnu[10]+f[9]*nvsqnu[9]+f[8]*nvsqnu[8]+f[7]*nvsqnu[7]+f[6]*nvsqnu[6]+f[5]*nvsqnu[5]+f[4]*nvsqnu[4]+f[3]*nvsqnu[3]+f[2]*nvsqnu[2]+f[1]*nvsqnu[1]+f[0]*nvsqnu[0]))/vmap_prime[1]; 
  nvsqnu_f[1] = (0.3535533905932737*(f[10]*nvsqnu[11]+nvsqnu[10]*f[11]+f[8]*nvsqnu[9]+nvsqnu[8]*f[9]+f[6]*nvsqnu[7]+nvsqnu[6]*f[7]+f[3]*nvsqnu[5]+nvsqnu[3]*f[5]+f[2]*nvsqnu[4]+nvsqnu[2]*f[4]+f[0]*nvsqnu[1]+nvsqnu[0]*f[1]))/vmap_prime[1]; 
  nvsqnu_f[2] = (0.3162277660168379*(f[7]*nvsqnu[11]+nvsqnu[7]*f[11])+0.3162277660168379*(f[6]*nvsqnu[10]+nvsqnu[6]*f[10]+f[4]*nvsqnu[9]+nvsqnu[4]*f[9])+0.3162277660168379*(f[2]*nvsqnu[8]+nvsqnu[2]*f[8])+0.3535533905932737*(f[5]*nvsqnu[7]+nvsqnu[5]*f[7]+f[3]*nvsqnu[6]+nvsqnu[3]*f[6]+f[1]*nvsqnu[4]+nvsqnu[1]*f[4]+f[0]*nvsqnu[2]+nvsqnu[0]*f[2]))/vmap_prime[1]; 
  nvsqnu_f[3] = (0.3535533905932737*(f[9]*nvsqnu[11]+nvsqnu[9]*f[11]+f[8]*nvsqnu[10]+nvsqnu[8]*f[10]+f[4]*nvsqnu[7]+nvsqnu[4]*f[7]+f[2]*nvsqnu[6]+nvsqnu[2]*f[6]+f[1]*nvsqnu[5]+nvsqnu[1]*f[5]+f[0]*nvsqnu[3]+nvsqnu[0]*f[3]))/vmap_prime[1]; 
  nvsqnu_f[4] = (0.3162277660168379*(f[6]*nvsqnu[11]+nvsqnu[6]*f[11])+0.3162277660168379*(f[7]*nvsqnu[10]+nvsqnu[7]*f[10]+f[2]*nvsqnu[9]+nvsqnu[2]*f[9])+0.3162277660168379*(f[4]*nvsqnu[8]+nvsqnu[4]*f[8])+0.3535533905932737*(f[3]*nvsqnu[7]+nvsqnu[3]*f[7]+f[5]*nvsqnu[6]+nvsqnu[5]*f[6]+f[0]*nvsqnu[4]+nvsqnu[0]*f[4]+f[1]*nvsqnu[2]+nvsqnu[1]*f[2]))/vmap_prime[1]; 
  nvsqnu_f[5] = (0.3535533905932737*(f[8]*nvsqnu[11]+nvsqnu[8]*f[11]+f[9]*nvsqnu[10]+nvsqnu[9]*f[10]+f[2]*nvsqnu[7]+nvsqnu[2]*f[7]+f[4]*nvsqnu[6]+nvsqnu[4]*f[6]+f[0]*nvsqnu[5]+nvsqnu[0]*f[5]+f[1]*nvsqnu[3]+nvsqnu[1]*f[3]))/vmap_prime[1]; 
  nvsqnu_f[6] = (0.3162277660168379*(f[4]*nvsqnu[11]+nvsqnu[4]*f[11])+0.3162277660168379*(f[2]*nvsqnu[10]+nvsqnu[2]*f[10]+f[7]*nvsqnu[9]+nvsqnu[7]*f[9])+0.3162277660168379*(f[6]*nvsqnu[8]+nvsqnu[6]*f[8])+0.3535533905932737*(f[1]*nvsqnu[7]+nvsqnu[1]*f[7]+f[0]*nvsqnu[6]+nvsqnu[0]*f[6]+f[4]*nvsqnu[5]+nvsqnu[4]*f[5]+f[2]*nvsqnu[3]+nvsqnu[2]*f[3]))/vmap_prime[1]; 
  nvsqnu_f[7] = (0.3162277660168379*(f[2]*nvsqnu[11]+nvsqnu[2]*f[11])+0.3162277660168379*(f[4]*nvsqnu[10]+nvsqnu[4]*f[10]+f[6]*nvsqnu[9]+nvsqnu[6]*f[9])+0.3162277660168379*(f[7]*nvsqnu[8]+nvsqnu[7]*f[8])+0.3535533905932737*(f[0]*nvsqnu[7]+nvsqnu[0]*f[7]+f[1]*nvsqnu[6]+nvsqnu[1]*f[6]+f[2]*nvsqnu[5]+nvsqnu[2]*f[5]+f[3]*nvsqnu[4]+nvsqnu[3]*f[4]))/vmap_prime[1]; 
  nvsqnu_f[8] = (0.2258769757263128*f[11]*nvsqnu[11]+0.3535533905932737*(f[5]*nvsqnu[11]+nvsqnu[5]*f[11])+0.2258769757263128*f[10]*nvsqnu[10]+0.3535533905932737*(f[3]*nvsqnu[10]+nvsqnu[3]*f[10])+0.2258769757263128*f[9]*nvsqnu[9]+0.3535533905932737*(f[1]*nvsqnu[9]+nvsqnu[1]*f[9])+0.2258769757263128*f[8]*nvsqnu[8]+0.3535533905932737*(f[0]*nvsqnu[8]+nvsqnu[0]*f[8])+0.3162277660168379*(f[7]*nvsqnu[7]+f[6]*nvsqnu[6]+f[4]*nvsqnu[4]+f[2]*nvsqnu[2]))/vmap_prime[1]; 
  nvsqnu_f[9] = ((0.2258769757263128*f[10]+0.3535533905932737*f[3])*nvsqnu[11]+0.2258769757263128*nvsqnu[10]*f[11]+0.3535533905932737*(nvsqnu[3]*f[11]+f[5]*nvsqnu[10]+nvsqnu[5]*f[10])+(0.2258769757263128*f[8]+0.3535533905932737*f[0])*nvsqnu[9]+0.2258769757263128*nvsqnu[8]*f[9]+0.3535533905932737*(nvsqnu[0]*f[9]+f[1]*nvsqnu[8]+nvsqnu[1]*f[8])+0.3162277660168379*(f[6]*nvsqnu[7]+nvsqnu[6]*f[7]+f[2]*nvsqnu[4]+nvsqnu[2]*f[4]))/vmap_prime[1]; 
  nvsqnu_f[10] = ((0.2258769757263128*f[9]+0.3535533905932737*f[1])*nvsqnu[11]+(0.2258769757263128*nvsqnu[9]+0.3535533905932737*nvsqnu[1])*f[11]+(0.2258769757263128*f[8]+0.3535533905932737*f[0])*nvsqnu[10]+0.2258769757263128*nvsqnu[8]*f[10]+0.3535533905932737*(nvsqnu[0]*f[10]+f[5]*nvsqnu[9]+nvsqnu[5]*f[9]+f[3]*nvsqnu[8]+nvsqnu[3]*f[8])+0.3162277660168379*(f[4]*nvsqnu[7]+nvsqnu[4]*f[7]+f[2]*nvsqnu[6]+nvsqnu[2]*f[6]))/vmap_prime[1]; 
  nvsqnu_f[11] = ((0.2258769757263128*f[8]+0.3535533905932737*f[0])*nvsqnu[11]+(0.2258769757263128*nvsqnu[8]+0.3535533905932737*nvsqnu[0])*f[11]+(0.2258769757263128*f[9]+0.3535533905932737*f[1])*nvsqnu[10]+0.2258769757263128*nvsqnu[9]*f[10]+0.3535533905932737*(nvsqnu[1]*f[10]+f[3]*nvsqnu[9]+nvsqnu[3]*f[9]+f[5]*nvsqnu[8]+nvsqnu[5]*f[8])+0.3162277660168379*(f[2]*nvsqnu[7]+nvsqnu[2]*f[7]+f[4]*nvsqnu[6]+nvsqnu[4]*f[6]))/vmap_prime[1]; 

  out[2] += -1.732050807568877*nvnu_f[0]*rdv2[0]; 
  out[3] += -1.732050807568877*nvsqnu_f[0]*rdv2[1]; 
  out[4] += -1.732050807568877*rdv2[0]*nvnu_f[1]; 
  out[5] += -1.732050807568877*nvsqnu_f[1]*rdv2[1]; 
  out[6] += (-1.732050807568877*rdv2[0]*nvnu_f[3])-1.732050807568877*rdv2[1]*nvsqnu_f[2]; 
  out[7] += (-1.732050807568877*rdv2[0]*nvnu_f[5])-1.732050807568877*rdv2[1]*nvsqnu_f[4]; 
  out[8] += -3.872983346207417*rdv2[0]*nvnu_f[2]; 
  out[9] += -3.872983346207417*rdv2[0]*nvnu_f[4]; 
  out[10] += (-1.732050807568877*rdv2[1]*nvsqnu_f[8])-3.872983346207417*rdv2[0]*nvnu_f[6]; 
  out[11] += (-1.732050807568877*rdv2[1]*nvsqnu_f[9])-3.872983346207417*rdv2[0]*nvnu_f[7]; 

  return 0.0; 
} 

