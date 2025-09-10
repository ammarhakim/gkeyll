#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_drag_surfmu_2x2v_ser_p1(const double *dxv, const double *vmap, const double *vmap_prime_l, const double *vmap_prime_c, const double *vmap_prime_r, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // dxv[4]: cell spacing. 
  // vmap: velocity space mapping.
  // vmap_prime_l,vmap_prime_c,vmap_prime_r: velocity space mapping derivative in left, center and right cells.
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[8]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 

  double rdv2 = 2.0/dxv[3]; 

  double Ghat_r[12] = {0.0}; 
  double Ghat_l[12] = {0.0}; 
  double alphaDrSurf_l[12] = {0.0}; 
  alphaDrSurf_l[0] = 2.0*nuSum[0]*vmap[2]-3.4641016151377544*nuSum[0]*vmap[3]; 
  alphaDrSurf_l[1] = 2.0*nuSum[1]*vmap[2]-3.4641016151377544*nuSum[1]*vmap[3]; 
  alphaDrSurf_l[2] = 2.0*nuSum[2]*vmap[2]-3.4641016151377544*nuSum[2]*vmap[3]; 
  alphaDrSurf_l[4] = 2.0*vmap[2]*nuSum[3]-3.4641016151377544*nuSum[3]*vmap[3]; 

  double alphaDrSurf_r[12] = {0.0}; 
  alphaDrSurf_r[0] = 3.4641016151377544*nuSum[0]*vmap[3]+2.0*nuSum[0]*vmap[2]; 
  alphaDrSurf_r[1] = 3.4641016151377544*nuSum[1]*vmap[3]+2.0*nuSum[1]*vmap[2]; 
  alphaDrSurf_r[2] = 3.4641016151377544*nuSum[2]*vmap[3]+2.0*nuSum[2]*vmap[2]; 
  alphaDrSurf_r[4] = 3.4641016151377544*nuSum[3]*vmap[3]+2.0*vmap[2]*nuSum[3]; 

  Ghat_l[0] = -((0.25*(1.7320508075688772*(alphaDrSurf_l[4]*fc[12]+alphaDrSurf_l[2]*fc[9]+alphaDrSurf_l[1]*fc[8])-1.0*alphaDrSurf_l[4]*fc[5]+1.7320508075688772*alphaDrSurf_l[0]*fc[4]-1.0*(alphaDrSurf_l[2]*fc[2]+alphaDrSurf_l[1]*fc[1]+alphaDrSurf_l[0]*fc[0])))/vmap_prime_c[1]); 
  Ghat_l[1] = -((0.25*(1.7320508075688772*(alphaDrSurf_l[2]*fc[12]+alphaDrSurf_l[4]*fc[9]+alphaDrSurf_l[0]*fc[8])-1.0*alphaDrSurf_l[2]*fc[5]+1.7320508075688772*alphaDrSurf_l[1]*fc[4]-1.0*(fc[2]*alphaDrSurf_l[4]+alphaDrSurf_l[0]*fc[1]+fc[0]*alphaDrSurf_l[1])))/vmap_prime_c[1]); 
  Ghat_l[2] = -((0.25*(1.7320508075688772*(alphaDrSurf_l[1]*fc[12]+alphaDrSurf_l[0]*fc[9]+alphaDrSurf_l[4]*fc[8])-1.0*alphaDrSurf_l[1]*fc[5]+1.7320508075688772*alphaDrSurf_l[2]*fc[4]-1.0*(fc[1]*alphaDrSurf_l[4]+alphaDrSurf_l[0]*fc[2]+fc[0]*alphaDrSurf_l[2])))/vmap_prime_c[1]); 
  Ghat_l[3] = -((0.25*(1.7320508075688772*(alphaDrSurf_l[4]*fc[15]+alphaDrSurf_l[2]*fc[14]+alphaDrSurf_l[1]*fc[13])-1.0*alphaDrSurf_l[4]*fc[11]+1.7320508075688772*alphaDrSurf_l[0]*fc[10]-1.0*(alphaDrSurf_l[2]*fc[7]+alphaDrSurf_l[1]*fc[6]+alphaDrSurf_l[0]*fc[3])))/vmap_prime_c[1]); 
  Ghat_l[4] = -((0.25*(1.7320508075688772*(alphaDrSurf_l[0]*fc[12]+alphaDrSurf_l[1]*fc[9]+alphaDrSurf_l[2]*fc[8])-1.0*alphaDrSurf_l[0]*fc[5]+1.7320508075688772*alphaDrSurf_l[4]*fc[4]-1.0*(fc[0]*alphaDrSurf_l[4]+alphaDrSurf_l[1]*fc[2]+fc[1]*alphaDrSurf_l[2])))/vmap_prime_c[1]); 
  Ghat_l[5] = -((0.25*(1.7320508075688772*(alphaDrSurf_l[2]*fc[15]+alphaDrSurf_l[4]*fc[14]+alphaDrSurf_l[0]*fc[13])-1.0*alphaDrSurf_l[2]*fc[11]+1.7320508075688772*alphaDrSurf_l[1]*fc[10]-1.0*(alphaDrSurf_l[4]*fc[7]+alphaDrSurf_l[0]*fc[6]+alphaDrSurf_l[1]*fc[3])))/vmap_prime_c[1]); 
  Ghat_l[6] = -((0.25*(1.7320508075688772*(alphaDrSurf_l[1]*fc[15]+alphaDrSurf_l[0]*fc[14]+alphaDrSurf_l[4]*fc[13])-1.0*alphaDrSurf_l[1]*fc[11]+1.7320508075688772*alphaDrSurf_l[2]*fc[10]-1.0*(alphaDrSurf_l[0]*fc[7]+alphaDrSurf_l[4]*fc[6]+alphaDrSurf_l[2]*fc[3])))/vmap_prime_c[1]); 
  Ghat_l[7] = -((0.25*(1.7320508075688772*(alphaDrSurf_l[0]*fc[15]+alphaDrSurf_l[1]*fc[14]+alphaDrSurf_l[2]*fc[13])-1.0*alphaDrSurf_l[0]*fc[11]+1.7320508075688772*alphaDrSurf_l[4]*fc[10]-1.0*(alphaDrSurf_l[1]*fc[7]+alphaDrSurf_l[2]*fc[6]+fc[3]*alphaDrSurf_l[4])))/vmap_prime_c[1]); 
  Ghat_l[8] = -((0.016666666666666666*(25.98076211353316*alphaDrSurf_l[4]*fc[23]+25.980762113533157*(alphaDrSurf_l[2]*fc[22]+alphaDrSurf_l[1]*fc[21])-15.0*alphaDrSurf_l[4]*fc[20]+25.98076211353316*alphaDrSurf_l[0]*fc[19]-15.000000000000002*(alphaDrSurf_l[2]*fc[18]+alphaDrSurf_l[1]*fc[17])-15.0*alphaDrSurf_l[0]*fc[16]))/vmap_prime_c[1]); 
  Ghat_l[9] = -((0.016666666666666666*(25.980762113533157*alphaDrSurf_l[2]*fc[23]+25.98076211353316*(alphaDrSurf_l[4]*fc[22]+alphaDrSurf_l[0]*fc[21])-15.000000000000002*alphaDrSurf_l[2]*fc[20]+25.980762113533157*alphaDrSurf_l[1]*fc[19]-15.0*(alphaDrSurf_l[4]*fc[18]+alphaDrSurf_l[0]*fc[17])-15.000000000000002*alphaDrSurf_l[1]*fc[16]))/vmap_prime_c[1]); 
  Ghat_l[10] = -((0.016666666666666666*(25.980762113533157*alphaDrSurf_l[1]*fc[23]+25.98076211353316*(alphaDrSurf_l[0]*fc[22]+alphaDrSurf_l[4]*fc[21])-15.000000000000002*alphaDrSurf_l[1]*fc[20]+25.980762113533157*alphaDrSurf_l[2]*fc[19]-15.0*(alphaDrSurf_l[0]*fc[18]+alphaDrSurf_l[4]*fc[17])-15.000000000000002*alphaDrSurf_l[2]*fc[16]))/vmap_prime_c[1]); 
  Ghat_l[11] = -((0.016666666666666666*(25.98076211353316*alphaDrSurf_l[0]*fc[23]+25.980762113533157*(alphaDrSurf_l[1]*fc[22]+alphaDrSurf_l[2]*fc[21])-15.0*alphaDrSurf_l[0]*fc[20]+25.98076211353316*alphaDrSurf_l[4]*fc[19]-15.000000000000002*(alphaDrSurf_l[1]*fc[18]+alphaDrSurf_l[2]*fc[17])-15.0*alphaDrSurf_l[4]*fc[16]))/vmap_prime_c[1]); 

  Ghat_r[0] = -((0.25*(1.7320508075688772*(alphaDrSurf_r[4]*fr[12]+alphaDrSurf_r[2]*fr[9]+alphaDrSurf_r[1]*fr[8])-1.0*alphaDrSurf_r[4]*fr[5]+1.7320508075688772*alphaDrSurf_r[0]*fr[4]-1.0*(alphaDrSurf_r[2]*fr[2]+alphaDrSurf_r[1]*fr[1]+alphaDrSurf_r[0]*fr[0])))/vmap_prime_r[1]); 
  Ghat_r[1] = -((0.25*(1.7320508075688772*(alphaDrSurf_r[2]*fr[12]+alphaDrSurf_r[4]*fr[9]+alphaDrSurf_r[0]*fr[8])-1.0*alphaDrSurf_r[2]*fr[5]+1.7320508075688772*alphaDrSurf_r[1]*fr[4]-1.0*(fr[2]*alphaDrSurf_r[4]+alphaDrSurf_r[0]*fr[1]+fr[0]*alphaDrSurf_r[1])))/vmap_prime_r[1]); 
  Ghat_r[2] = -((0.25*(1.7320508075688772*(alphaDrSurf_r[1]*fr[12]+alphaDrSurf_r[0]*fr[9]+alphaDrSurf_r[4]*fr[8])-1.0*alphaDrSurf_r[1]*fr[5]+1.7320508075688772*alphaDrSurf_r[2]*fr[4]-1.0*(fr[1]*alphaDrSurf_r[4]+alphaDrSurf_r[0]*fr[2]+fr[0]*alphaDrSurf_r[2])))/vmap_prime_r[1]); 
  Ghat_r[3] = -((0.25*(1.7320508075688772*(alphaDrSurf_r[4]*fr[15]+alphaDrSurf_r[2]*fr[14]+alphaDrSurf_r[1]*fr[13])-1.0*alphaDrSurf_r[4]*fr[11]+1.7320508075688772*alphaDrSurf_r[0]*fr[10]-1.0*(alphaDrSurf_r[2]*fr[7]+alphaDrSurf_r[1]*fr[6]+alphaDrSurf_r[0]*fr[3])))/vmap_prime_r[1]); 
  Ghat_r[4] = -((0.25*(1.7320508075688772*(alphaDrSurf_r[0]*fr[12]+alphaDrSurf_r[1]*fr[9]+alphaDrSurf_r[2]*fr[8])-1.0*alphaDrSurf_r[0]*fr[5]+1.7320508075688772*alphaDrSurf_r[4]*fr[4]-1.0*(fr[0]*alphaDrSurf_r[4]+alphaDrSurf_r[1]*fr[2]+fr[1]*alphaDrSurf_r[2])))/vmap_prime_r[1]); 
  Ghat_r[5] = -((0.25*(1.7320508075688772*(alphaDrSurf_r[2]*fr[15]+alphaDrSurf_r[4]*fr[14]+alphaDrSurf_r[0]*fr[13])-1.0*alphaDrSurf_r[2]*fr[11]+1.7320508075688772*alphaDrSurf_r[1]*fr[10]-1.0*(alphaDrSurf_r[4]*fr[7]+alphaDrSurf_r[0]*fr[6]+alphaDrSurf_r[1]*fr[3])))/vmap_prime_r[1]); 
  Ghat_r[6] = -((0.25*(1.7320508075688772*(alphaDrSurf_r[1]*fr[15]+alphaDrSurf_r[0]*fr[14]+alphaDrSurf_r[4]*fr[13])-1.0*alphaDrSurf_r[1]*fr[11]+1.7320508075688772*alphaDrSurf_r[2]*fr[10]-1.0*(alphaDrSurf_r[0]*fr[7]+alphaDrSurf_r[4]*fr[6]+alphaDrSurf_r[2]*fr[3])))/vmap_prime_r[1]); 
  Ghat_r[7] = -((0.25*(1.7320508075688772*(alphaDrSurf_r[0]*fr[15]+alphaDrSurf_r[1]*fr[14]+alphaDrSurf_r[2]*fr[13])-1.0*alphaDrSurf_r[0]*fr[11]+1.7320508075688772*alphaDrSurf_r[4]*fr[10]-1.0*(alphaDrSurf_r[1]*fr[7]+alphaDrSurf_r[2]*fr[6]+fr[3]*alphaDrSurf_r[4])))/vmap_prime_r[1]); 
  Ghat_r[8] = -((0.016666666666666666*(25.98076211353316*alphaDrSurf_r[4]*fr[23]+25.980762113533157*(alphaDrSurf_r[2]*fr[22]+alphaDrSurf_r[1]*fr[21])-15.0*alphaDrSurf_r[4]*fr[20]+25.98076211353316*alphaDrSurf_r[0]*fr[19]-15.000000000000002*(alphaDrSurf_r[2]*fr[18]+alphaDrSurf_r[1]*fr[17])-15.0*alphaDrSurf_r[0]*fr[16]))/vmap_prime_r[1]); 
  Ghat_r[9] = -((0.016666666666666666*(25.980762113533157*alphaDrSurf_r[2]*fr[23]+25.98076211353316*(alphaDrSurf_r[4]*fr[22]+alphaDrSurf_r[0]*fr[21])-15.000000000000002*alphaDrSurf_r[2]*fr[20]+25.980762113533157*alphaDrSurf_r[1]*fr[19]-15.0*(alphaDrSurf_r[4]*fr[18]+alphaDrSurf_r[0]*fr[17])-15.000000000000002*alphaDrSurf_r[1]*fr[16]))/vmap_prime_r[1]); 
  Ghat_r[10] = -((0.016666666666666666*(25.980762113533157*alphaDrSurf_r[1]*fr[23]+25.98076211353316*(alphaDrSurf_r[0]*fr[22]+alphaDrSurf_r[4]*fr[21])-15.000000000000002*alphaDrSurf_r[1]*fr[20]+25.980762113533157*alphaDrSurf_r[2]*fr[19]-15.0*(alphaDrSurf_r[0]*fr[18]+alphaDrSurf_r[4]*fr[17])-15.000000000000002*alphaDrSurf_r[2]*fr[16]))/vmap_prime_r[1]); 
  Ghat_r[11] = -((0.016666666666666666*(25.98076211353316*alphaDrSurf_r[0]*fr[23]+25.980762113533157*(alphaDrSurf_r[1]*fr[22]+alphaDrSurf_r[2]*fr[21])-15.0*alphaDrSurf_r[0]*fr[20]+25.98076211353316*alphaDrSurf_r[4]*fr[19]-15.000000000000002*(alphaDrSurf_r[1]*fr[18]+alphaDrSurf_r[2]*fr[17])-15.0*alphaDrSurf_r[4]*fr[16]))/vmap_prime_r[1]); 

  out[0] += (0.7071067811865475*Ghat_r[0]-0.7071067811865475*Ghat_l[0])*rdv2; 
  out[1] += (0.7071067811865475*Ghat_r[1]-0.7071067811865475*Ghat_l[1])*rdv2; 
  out[2] += (0.7071067811865475*Ghat_r[2]-0.7071067811865475*Ghat_l[2])*rdv2; 
  out[3] += (0.7071067811865475*Ghat_r[3]-0.7071067811865475*Ghat_l[3])*rdv2; 
  out[4] += 1.224744871391589*(Ghat_r[0]+Ghat_l[0])*rdv2; 
  out[5] += (0.7071067811865475*Ghat_r[4]-0.7071067811865475*Ghat_l[4])*rdv2; 
  out[6] += (0.7071067811865475*Ghat_r[5]-0.7071067811865475*Ghat_l[5])*rdv2; 
  out[7] += (0.7071067811865475*Ghat_r[6]-0.7071067811865475*Ghat_l[6])*rdv2; 
  out[8] += 1.224744871391589*(Ghat_r[1]+Ghat_l[1])*rdv2; 
  out[9] += 1.224744871391589*(Ghat_r[2]+Ghat_l[2])*rdv2; 
  out[10] += 1.224744871391589*(Ghat_r[3]+Ghat_l[3])*rdv2; 
  out[11] += (0.7071067811865475*Ghat_r[7]-0.7071067811865475*Ghat_l[7])*rdv2; 
  out[12] += 1.224744871391589*(Ghat_r[4]+Ghat_l[4])*rdv2; 
  out[13] += 1.224744871391589*(Ghat_r[5]+Ghat_l[5])*rdv2; 
  out[14] += 1.224744871391589*(Ghat_r[6]+Ghat_l[6])*rdv2; 
  out[15] += 1.224744871391589*(Ghat_r[7]+Ghat_l[7])*rdv2; 
  out[16] += (0.7071067811865475*Ghat_r[8]-0.7071067811865475*Ghat_l[8])*rdv2; 
  out[17] += (0.7071067811865475*Ghat_r[9]-0.7071067811865475*Ghat_l[9])*rdv2; 
  out[18] += (0.7071067811865475*Ghat_r[10]-0.7071067811865475*Ghat_l[10])*rdv2; 
  out[19] += 1.224744871391589*(Ghat_r[8]+Ghat_l[8])*rdv2; 
  out[20] += (0.7071067811865475*Ghat_r[11]-0.7071067811865475*Ghat_l[11])*rdv2; 
  out[21] += 1.224744871391589*(Ghat_r[9]+Ghat_l[9])*rdv2; 
  out[22] += 1.224744871391589*(Ghat_r[10]+Ghat_l[10])*rdv2; 
  out[23] += 1.224744871391589*(Ghat_r[11]+Ghat_l[11])*rdv2; 

  return 0.;

} 
