#include <gkyl_lbo_pkpm_kernels.h> 
GKYL_CU_DH double lbo_pkpm_drag_vol_1x1v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:       Cell-center coordinates. 
  // dxv[NDIM]:     Cell spacing. 
  // nuSum:         Collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: Sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // f:             Input distribution functions [F_0, T_perp/m G_1 = T_perp/m (F_0 - F_1)].
  // out:           Incremented output distribution functions. 
  const double rdvpar = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double *F_0 = &f[0]; 
  const double *G_1 = &f[9]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[9]; 

  const double *sumNuUPar = &nuPrimMomsSum[0]; 

  double alphaDrag[9]; 
  alphaDrag[0] = rdvpar*(1.414213562373095*sumNuUPar[0]-1.414213562373095*nuSum[0]*wvpar); 
  alphaDrag[1] = rdvpar*(1.414213562373095*sumNuUPar[1]-1.414213562373095*nuSum[1]*wvpar); 
  alphaDrag[2] = -0.408248290463863*nuSum[0]*dvpar*rdvpar; 
  alphaDrag[3] = -0.408248290463863*nuSum[1]*dvpar*rdvpar; 
  alphaDrag[4] = rdvpar*(1.414213562373095*sumNuUPar[2]-1.414213562373095*nuSum[2]*wvpar); 
  alphaDrag[6] = -0.408248290463863*nuSum[2]*dvpar*rdvpar; 

  out_F_0[2] += 0.8660254037844386*F_0[6]*alphaDrag[6]+0.8660254037844386*F_0[4]*alphaDrag[4]+0.8660254037844386*F_0[3]*alphaDrag[3]+0.8660254037844386*F_0[2]*alphaDrag[2]+0.8660254037844386*F_0[1]*alphaDrag[1]+0.8660254037844386*F_0[0]*alphaDrag[0]; 
  out_F_0[3] += 0.7745966692414834*F_0[3]*alphaDrag[6]+0.7745966692414834*alphaDrag[3]*F_0[6]+0.7745966692414833*F_0[1]*alphaDrag[4]+0.7745966692414833*alphaDrag[1]*F_0[4]+0.8660254037844386*F_0[2]*alphaDrag[3]+0.8660254037844386*alphaDrag[2]*F_0[3]+0.8660254037844386*F_0[0]*alphaDrag[1]+0.8660254037844386*alphaDrag[0]*F_0[1]; 
  out_F_0[5] += 1.732050807568877*alphaDrag[6]*F_0[8]+1.732050807568877*alphaDrag[3]*F_0[7]+1.936491673103709*F_0[4]*alphaDrag[6]+1.936491673103709*alphaDrag[4]*F_0[6]+1.732050807568877*alphaDrag[2]*F_0[5]+1.936491673103709*F_0[1]*alphaDrag[3]+1.936491673103709*alphaDrag[1]*F_0[3]+1.936491673103709*F_0[0]*alphaDrag[2]+1.936491673103709*alphaDrag[0]*F_0[2]; 
  out_F_0[6] += 0.5532833351724881*F_0[6]*alphaDrag[6]+0.8660254037844386*F_0[2]*alphaDrag[6]+0.8660254037844386*alphaDrag[2]*F_0[6]+0.5532833351724881*F_0[4]*alphaDrag[4]+0.8660254037844387*F_0[0]*alphaDrag[4]+0.8660254037844387*alphaDrag[0]*F_0[4]+0.7745966692414834*F_0[3]*alphaDrag[3]+0.7745966692414834*F_0[1]*alphaDrag[1]; 
  out_F_0[7] += 1.549193338482967*alphaDrag[3]*F_0[8]+1.549193338482967*alphaDrag[6]*F_0[7]+1.732050807568877*alphaDrag[2]*F_0[7]+1.732050807568877*F_0[1]*alphaDrag[6]+1.732050807568877*alphaDrag[1]*F_0[6]+1.732050807568877*alphaDrag[3]*F_0[5]+1.732050807568877*F_0[3]*alphaDrag[4]+1.732050807568877*alphaDrag[3]*F_0[4]+1.936491673103709*F_0[0]*alphaDrag[3]+1.936491673103709*alphaDrag[0]*F_0[3]+1.936491673103709*F_0[1]*alphaDrag[2]+1.936491673103709*alphaDrag[1]*F_0[2]; 
  out_F_0[8] += 1.106566670344976*alphaDrag[6]*F_0[8]+1.732050807568877*alphaDrag[2]*F_0[8]+1.549193338482967*alphaDrag[3]*F_0[7]+1.732050807568877*F_0[5]*alphaDrag[6]+1.237179148263484*F_0[4]*alphaDrag[6]+1.936491673103709*F_0[0]*alphaDrag[6]+1.237179148263484*alphaDrag[4]*F_0[6]+1.936491673103709*alphaDrag[0]*F_0[6]+1.936491673103709*F_0[2]*alphaDrag[4]+1.936491673103709*alphaDrag[2]*F_0[4]+1.732050807568877*F_0[1]*alphaDrag[3]+1.732050807568877*alphaDrag[1]*F_0[3]; 
  out_G_1[2] += 0.8660254037844386*G_1[6]*alphaDrag[6]+0.8660254037844386*G_1[4]*alphaDrag[4]+0.8660254037844386*G_1[3]*alphaDrag[3]+0.8660254037844386*G_1[2]*alphaDrag[2]+0.8660254037844386*G_1[1]*alphaDrag[1]+0.8660254037844386*G_1[0]*alphaDrag[0]; 
  out_G_1[3] += 0.7745966692414834*G_1[3]*alphaDrag[6]+0.7745966692414834*alphaDrag[3]*G_1[6]+0.7745966692414833*G_1[1]*alphaDrag[4]+0.7745966692414833*alphaDrag[1]*G_1[4]+0.8660254037844386*G_1[2]*alphaDrag[3]+0.8660254037844386*alphaDrag[2]*G_1[3]+0.8660254037844386*G_1[0]*alphaDrag[1]+0.8660254037844386*alphaDrag[0]*G_1[1]; 
  out_G_1[5] += 1.732050807568877*alphaDrag[6]*G_1[8]+1.732050807568877*alphaDrag[3]*G_1[7]+1.936491673103709*G_1[4]*alphaDrag[6]+1.936491673103709*alphaDrag[4]*G_1[6]+1.732050807568877*alphaDrag[2]*G_1[5]+1.936491673103709*G_1[1]*alphaDrag[3]+1.936491673103709*alphaDrag[1]*G_1[3]+1.936491673103709*G_1[0]*alphaDrag[2]+1.936491673103709*alphaDrag[0]*G_1[2]; 
  out_G_1[6] += 0.5532833351724881*G_1[6]*alphaDrag[6]+0.8660254037844386*G_1[2]*alphaDrag[6]+0.8660254037844386*alphaDrag[2]*G_1[6]+0.5532833351724881*G_1[4]*alphaDrag[4]+0.8660254037844387*G_1[0]*alphaDrag[4]+0.8660254037844387*alphaDrag[0]*G_1[4]+0.7745966692414834*G_1[3]*alphaDrag[3]+0.7745966692414834*G_1[1]*alphaDrag[1]; 
  out_G_1[7] += 1.549193338482967*alphaDrag[3]*G_1[8]+1.549193338482967*alphaDrag[6]*G_1[7]+1.732050807568877*alphaDrag[2]*G_1[7]+1.732050807568877*G_1[1]*alphaDrag[6]+1.732050807568877*alphaDrag[1]*G_1[6]+1.732050807568877*alphaDrag[3]*G_1[5]+1.732050807568877*G_1[3]*alphaDrag[4]+1.732050807568877*alphaDrag[3]*G_1[4]+1.936491673103709*G_1[0]*alphaDrag[3]+1.936491673103709*alphaDrag[0]*G_1[3]+1.936491673103709*G_1[1]*alphaDrag[2]+1.936491673103709*alphaDrag[1]*G_1[2]; 
  out_G_1[8] += 1.106566670344976*alphaDrag[6]*G_1[8]+1.732050807568877*alphaDrag[2]*G_1[8]+1.549193338482967*alphaDrag[3]*G_1[7]+1.732050807568877*G_1[5]*alphaDrag[6]+1.237179148263484*G_1[4]*alphaDrag[6]+1.936491673103709*G_1[0]*alphaDrag[6]+1.237179148263484*alphaDrag[4]*G_1[6]+1.936491673103709*alphaDrag[0]*G_1[6]+1.936491673103709*G_1[2]*alphaDrag[4]+1.936491673103709*alphaDrag[2]*G_1[4]+1.732050807568877*G_1[1]*alphaDrag[3]+1.732050807568877*alphaDrag[1]*G_1[3]; 

  return fabs(1.25*alphaDrag[0]-1.397542485937369*alphaDrag[4]); 

} 
