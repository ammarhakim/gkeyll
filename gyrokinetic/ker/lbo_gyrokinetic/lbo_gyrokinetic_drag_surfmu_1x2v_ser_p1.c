#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_drag_surfmu_1x2v_ser_p1(const double *dxv, const double *vmap, const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // dxv[3]: cell spacing. 
  // vmap: velocity space mapping.
  // vmap_prime_l,vmap_prime_c,vmap_prime_r: velocity space mapping derivative in left, center and right cells.
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[4]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 

  double rdv2 = 2.0/dxv[2]; 

  double Ghat_r[6] = {0.0}; 
  double Ghat_l[6] = {0.0}; 
  double alphaDrSurf_l[6] = {0.0}; 
  alphaDrSurf_l[0] = 2.0*nuSum[0]*vmap[2]-3.4641016151377544*nuSum[0]*vmap[3]; 
  alphaDrSurf_l[1] = 2.0*nuSum[1]*vmap[2]-3.4641016151377544*nuSum[1]*vmap[3]; 

  double alphaDrSurf_r[6] = {0.0}; 
  alphaDrSurf_r[0] = 3.4641016151377544*nuSum[0]*vmap[3]+2.0*nuSum[0]*vmap[2]; 
  alphaDrSurf_r[1] = 3.4641016151377544*nuSum[1]*vmap[3]+2.0*nuSum[1]*vmap[2]; 

  Ghat_l[0] = -((0.25*(2.4494897427831783*(alphaDrSurf_l[1]*fc[5]+alphaDrSurf_l[0]*fc[3])-1.4142135623730951*(alphaDrSurf_l[1]*fc[1]+alphaDrSurf_l[0]*fc[0])))/vmap_prime_c[1]); 
  Ghat_l[1] = -((0.25*(2.4494897427831783*(alphaDrSurf_l[0]*fc[5]+alphaDrSurf_l[1]*fc[3])-1.4142135623730951*(alphaDrSurf_l[0]*fc[1]+fc[0]*alphaDrSurf_l[1])))/vmap_prime_c[1]); 
  Ghat_l[2] = -((0.25*(2.4494897427831783*(alphaDrSurf_l[1]*fc[7]+alphaDrSurf_l[0]*fc[6])-1.4142135623730951*(alphaDrSurf_l[1]*fc[4]+alphaDrSurf_l[0]*fc[2])))/vmap_prime_c[1]); 
  Ghat_l[3] = -((0.25*(2.4494897427831783*(alphaDrSurf_l[0]*fc[7]+alphaDrSurf_l[1]*fc[6])-1.4142135623730951*(alphaDrSurf_l[0]*fc[4]+alphaDrSurf_l[1]*fc[2])))/vmap_prime_c[1]); 
  Ghat_l[4] = -((0.016666666666666666*(36.74234614174767*alphaDrSurf_l[1]*fc[11]+36.74234614174768*alphaDrSurf_l[0]*fc[10]-21.21320343559643*alphaDrSurf_l[1]*fc[9]-21.213203435596427*alphaDrSurf_l[0]*fc[8]))/vmap_prime_c[1]); 
  Ghat_l[5] = -((0.016666666666666666*(36.74234614174768*alphaDrSurf_l[0]*fc[11]+36.74234614174767*alphaDrSurf_l[1]*fc[10]-21.213203435596427*alphaDrSurf_l[0]*fc[9]-21.21320343559643*alphaDrSurf_l[1]*fc[8]))/vmap_prime_c[1]); 

  Ghat_r[0] = -((0.25*(2.4494897427831783*(alphaDrSurf_r[1]*fr[5]+alphaDrSurf_r[0]*fr[3])-1.4142135623730951*(alphaDrSurf_r[1]*fr[1]+alphaDrSurf_r[0]*fr[0])))/vmap_prime_r[1]); 
  Ghat_r[1] = -((0.25*(2.4494897427831783*(alphaDrSurf_r[0]*fr[5]+alphaDrSurf_r[1]*fr[3])-1.4142135623730951*(alphaDrSurf_r[0]*fr[1]+fr[0]*alphaDrSurf_r[1])))/vmap_prime_r[1]); 
  Ghat_r[2] = -((0.25*(2.4494897427831783*(alphaDrSurf_r[1]*fr[7]+alphaDrSurf_r[0]*fr[6])-1.4142135623730951*(alphaDrSurf_r[1]*fr[4]+alphaDrSurf_r[0]*fr[2])))/vmap_prime_r[1]); 
  Ghat_r[3] = -((0.25*(2.4494897427831783*(alphaDrSurf_r[0]*fr[7]+alphaDrSurf_r[1]*fr[6])-1.4142135623730951*(alphaDrSurf_r[0]*fr[4]+alphaDrSurf_r[1]*fr[2])))/vmap_prime_r[1]); 
  Ghat_r[4] = -((0.016666666666666666*(36.74234614174767*alphaDrSurf_r[1]*fr[11]+36.74234614174768*alphaDrSurf_r[0]*fr[10]-21.21320343559643*alphaDrSurf_r[1]*fr[9]-21.213203435596427*alphaDrSurf_r[0]*fr[8]))/vmap_prime_r[1]); 
  Ghat_r[5] = -((0.016666666666666666*(36.74234614174768*alphaDrSurf_r[0]*fr[11]+36.74234614174767*alphaDrSurf_r[1]*fr[10]-21.213203435596427*alphaDrSurf_r[0]*fr[9]-21.21320343559643*alphaDrSurf_r[1]*fr[8]))/vmap_prime_r[1]); 

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

  return 0.;

} 
