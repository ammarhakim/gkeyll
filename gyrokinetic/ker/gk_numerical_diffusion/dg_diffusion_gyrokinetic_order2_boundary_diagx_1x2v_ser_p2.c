#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_diagx_1x2v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qGhost, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qGhost/qSkin: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],2.);

  double boundSurf_incr[20] = {0.0}; 

  if (edge == -1) { 

  boundSurf_incr[1] = 0.9682458365518543*coeff[0]*qSkin[7]-1.25*coeff[0]*qSkin[1]+0.8660254037844386*coeff[0]*qSkin[0]; 
  boundSurf_incr[4] = 0.9682458365518543*coeff[0]*qSkin[11]-1.25*coeff[0]*qSkin[4]+0.8660254037844386*coeff[0]*qSkin[2]; 
  boundSurf_incr[5] = 0.9682458365518543*coeff[0]*qSkin[13]-1.25*coeff[0]*qSkin[5]+0.8660254037844386*coeff[0]*qSkin[3]; 
  boundSurf_incr[7] = -(3.75*coeff[0]*qSkin[7])+4.841229182759272*coeff[0]*qSkin[1]-3.3541019662496847*coeff[0]*qSkin[0]; 
  boundSurf_incr[10] = 0.9682458365518543*coeff[0]*qSkin[17]-1.25*coeff[0]*qSkin[10]+0.8660254037844386*coeff[0]*qSkin[6]; 
  boundSurf_incr[11] = -(3.75*coeff[0]*qSkin[11])+4.841229182759271*coeff[0]*qSkin[4]-3.3541019662496843*coeff[0]*qSkin[2]; 
  boundSurf_incr[12] = 0.8660254037844387*coeff[0]*qSkin[8]-1.25*coeff[0]*qSkin[12]; 
  boundSurf_incr[13] = -(3.75*coeff[0]*qSkin[13])+4.841229182759271*coeff[0]*qSkin[5]-3.3541019662496843*coeff[0]*qSkin[3]; 
  boundSurf_incr[15] = 0.8660254037844387*coeff[0]*qSkin[9]-1.25*coeff[0]*qSkin[15]; 
  boundSurf_incr[17] = -(3.75*coeff[0]*qSkin[17])+4.841229182759272*coeff[0]*qSkin[10]-3.3541019662496847*coeff[0]*qSkin[6]; 
  boundSurf_incr[18] = 0.8660254037844387*coeff[0]*qSkin[14]-1.25*coeff[0]*qSkin[18]; 
  boundSurf_incr[19] = 0.8660254037844387*coeff[0]*qSkin[16]-1.25*coeff[0]*qSkin[19]; 

  } else { 

  boundSurf_incr[1] = -(0.9682458365518543*coeff[0]*qSkin[7])-1.25*coeff[0]*qSkin[1]-0.8660254037844386*coeff[0]*qSkin[0]; 
  boundSurf_incr[4] = -(0.9682458365518543*coeff[0]*qSkin[11])-1.25*coeff[0]*qSkin[4]-0.8660254037844386*coeff[0]*qSkin[2]; 
  boundSurf_incr[5] = -(0.9682458365518543*coeff[0]*qSkin[13])-1.25*coeff[0]*qSkin[5]-0.8660254037844386*coeff[0]*qSkin[3]; 
  boundSurf_incr[7] = -(3.75*coeff[0]*qSkin[7])-4.841229182759272*coeff[0]*qSkin[1]-3.3541019662496847*coeff[0]*qSkin[0]; 
  boundSurf_incr[10] = -(0.9682458365518543*coeff[0]*qSkin[17])-1.25*coeff[0]*qSkin[10]-0.8660254037844386*coeff[0]*qSkin[6]; 
  boundSurf_incr[11] = -(3.75*coeff[0]*qSkin[11])-4.841229182759271*coeff[0]*qSkin[4]-3.3541019662496843*coeff[0]*qSkin[2]; 
  boundSurf_incr[12] = -(1.25*coeff[0]*qSkin[12])-0.8660254037844387*coeff[0]*qSkin[8]; 
  boundSurf_incr[13] = -(3.75*coeff[0]*qSkin[13])-4.841229182759271*coeff[0]*qSkin[5]-3.3541019662496843*coeff[0]*qSkin[3]; 
  boundSurf_incr[15] = -(1.25*coeff[0]*qSkin[15])-0.8660254037844387*coeff[0]*qSkin[9]; 
  boundSurf_incr[17] = -(3.75*coeff[0]*qSkin[17])-4.841229182759272*coeff[0]*qSkin[10]-3.3541019662496847*coeff[0]*qSkin[6]; 
  boundSurf_incr[18] = -(1.25*coeff[0]*qSkin[18])-0.8660254037844387*coeff[0]*qSkin[14]; 
  boundSurf_incr[19] = -(1.25*coeff[0]*qSkin[19])-0.8660254037844387*coeff[0]*qSkin[16]; 

  }

  out[0] += boundSurf_incr[0]*rdx2Sq; 
  out[1] += boundSurf_incr[1]*rdx2Sq; 
  out[2] += boundSurf_incr[2]*rdx2Sq; 
  out[3] += boundSurf_incr[3]*rdx2Sq; 
  out[4] += boundSurf_incr[4]*rdx2Sq; 
  out[5] += boundSurf_incr[5]*rdx2Sq; 
  out[6] += boundSurf_incr[6]*rdx2Sq; 
  out[7] += boundSurf_incr[7]*rdx2Sq; 
  out[8] += boundSurf_incr[8]*rdx2Sq; 
  out[9] += boundSurf_incr[9]*rdx2Sq; 
  out[10] += boundSurf_incr[10]*rdx2Sq; 
  out[11] += boundSurf_incr[11]*rdx2Sq; 
  out[12] += boundSurf_incr[12]*rdx2Sq; 
  out[13] += boundSurf_incr[13]*rdx2Sq; 
  out[14] += boundSurf_incr[14]*rdx2Sq; 
  out[15] += boundSurf_incr[15]*rdx2Sq; 
  out[16] += boundSurf_incr[16]*rdx2Sq; 
  out[17] += boundSurf_incr[17]*rdx2Sq; 
  out[18] += boundSurf_incr[18]*rdx2Sq; 
  out[19] += boundSurf_incr[19]*rdx2Sq; 

  return 0.;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_diagx_1x2v_ser_p2_varcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qGhost, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qGhost/qSkin: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],2.);

  double fGhost[20];
  fGhost[0] = 0.7071067811865475*(jacobgeo_inv[2]*qGhost[7]+jacobgeo_inv[1]*qGhost[1]+jacobgeo_inv[0]*qGhost[0]); 
  fGhost[1] = 0.1414213562373095*(4.47213595499958*(jacobgeo_inv[1]*qGhost[7]+qGhost[1]*jacobgeo_inv[2])+5.0*(jacobgeo_inv[0]*qGhost[1]+qGhost[0]*jacobgeo_inv[1])); 
  fGhost[2] = 0.04714045207910316*(15.000000000000002*jacobgeo_inv[2]*qGhost[11]+15.0*(jacobgeo_inv[1]*qGhost[4]+jacobgeo_inv[0]*qGhost[2])); 
  fGhost[3] = 0.04714045207910316*(15.000000000000002*jacobgeo_inv[2]*qGhost[13]+15.0*(jacobgeo_inv[1]*qGhost[5]+jacobgeo_inv[0]*qGhost[3])); 
  fGhost[4] = 0.04714045207910316*(13.416407864998737*jacobgeo_inv[1]*qGhost[11]+(13.416407864998739*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qGhost[4]+15.0*jacobgeo_inv[1]*qGhost[2]); 
  fGhost[5] = 0.04714045207910316*(13.416407864998737*jacobgeo_inv[1]*qGhost[13]+(13.416407864998739*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qGhost[5]+15.0*jacobgeo_inv[1]*qGhost[3]); 
  fGhost[6] = 0.7071067811865475*(jacobgeo_inv[2]*qGhost[17]+jacobgeo_inv[1]*qGhost[10]+jacobgeo_inv[0]*qGhost[6]); 
  fGhost[7] = 0.020203050891044214*((22.3606797749979*jacobgeo_inv[2]+35.0*jacobgeo_inv[0])*qGhost[7]+35.0*qGhost[0]*jacobgeo_inv[2]+31.304951684997057*jacobgeo_inv[1]*qGhost[1]); 
  fGhost[8] = 0.04714045207910316*(15.000000000000002*jacobgeo_inv[1]*qGhost[12]+15.0*jacobgeo_inv[0]*qGhost[8]); 
  fGhost[9] = 0.04714045207910316*(15.000000000000002*jacobgeo_inv[1]*qGhost[15]+15.0*jacobgeo_inv[0]*qGhost[9]); 
  fGhost[10] = 0.1414213562373095*(4.47213595499958*jacobgeo_inv[1]*qGhost[17]+(4.47213595499958*jacobgeo_inv[2]+5.0*jacobgeo_inv[0])*qGhost[10]+5.0*jacobgeo_inv[1]*qGhost[6]); 
  fGhost[11] = 0.006734350297014738*((67.0820393249937*jacobgeo_inv[2]+105.0*jacobgeo_inv[0])*qGhost[11]+93.91485505499116*jacobgeo_inv[1]*qGhost[4]+105.00000000000001*jacobgeo_inv[2]*qGhost[2]); 
  fGhost[12] = 0.04714045207910316*((13.416407864998739*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qGhost[12]+15.000000000000002*jacobgeo_inv[1]*qGhost[8]); 
  fGhost[13] = 0.006734350297014738*((67.0820393249937*jacobgeo_inv[2]+105.0*jacobgeo_inv[0])*qGhost[13]+93.91485505499116*jacobgeo_inv[1]*qGhost[5]+105.00000000000001*jacobgeo_inv[2]*qGhost[3]); 
  fGhost[14] = 0.04714045207910316*(15.000000000000002*jacobgeo_inv[1]*qGhost[18]+15.0*jacobgeo_inv[0]*qGhost[14]); 
  fGhost[15] = 0.04714045207910316*((13.416407864998739*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qGhost[15]+15.000000000000002*jacobgeo_inv[1]*qGhost[9]); 
  fGhost[16] = 0.04714045207910316*(15.000000000000002*jacobgeo_inv[1]*qGhost[19]+15.0*jacobgeo_inv[0]*qGhost[16]); 
  fGhost[17] = 0.020203050891044214*((22.3606797749979*jacobgeo_inv[2]+35.0*jacobgeo_inv[0])*qGhost[17]+31.304951684997057*jacobgeo_inv[1]*qGhost[10]+35.0*jacobgeo_inv[2]*qGhost[6]); 
  fGhost[18] = 0.04714045207910316*((13.416407864998739*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qGhost[18]+15.000000000000002*jacobgeo_inv[1]*qGhost[14]); 
  fGhost[19] = 0.04714045207910316*((13.416407864998739*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qGhost[19]+15.000000000000002*jacobgeo_inv[1]*qGhost[16]); 

  double fSkin[20];
  fSkin[0] = 0.7071067811865475*(jacobgeo_inv[2]*qSkin[7]+jacobgeo_inv[1]*qSkin[1]+jacobgeo_inv[0]*qSkin[0]); 
  fSkin[1] = 0.1414213562373095*(4.47213595499958*(jacobgeo_inv[1]*qSkin[7]+qSkin[1]*jacobgeo_inv[2])+5.0*(jacobgeo_inv[0]*qSkin[1]+qSkin[0]*jacobgeo_inv[1])); 
  fSkin[2] = 0.04714045207910316*(15.000000000000002*jacobgeo_inv[2]*qSkin[11]+15.0*(jacobgeo_inv[1]*qSkin[4]+jacobgeo_inv[0]*qSkin[2])); 
  fSkin[3] = 0.04714045207910316*(15.000000000000002*jacobgeo_inv[2]*qSkin[13]+15.0*(jacobgeo_inv[1]*qSkin[5]+jacobgeo_inv[0]*qSkin[3])); 
  fSkin[4] = 0.04714045207910316*(13.416407864998737*jacobgeo_inv[1]*qSkin[11]+(13.416407864998739*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qSkin[4]+15.0*jacobgeo_inv[1]*qSkin[2]); 
  fSkin[5] = 0.04714045207910316*(13.416407864998737*jacobgeo_inv[1]*qSkin[13]+(13.416407864998739*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qSkin[5]+15.0*jacobgeo_inv[1]*qSkin[3]); 
  fSkin[6] = 0.7071067811865475*(jacobgeo_inv[2]*qSkin[17]+jacobgeo_inv[1]*qSkin[10]+jacobgeo_inv[0]*qSkin[6]); 
  fSkin[7] = 0.020203050891044214*((22.3606797749979*jacobgeo_inv[2]+35.0*jacobgeo_inv[0])*qSkin[7]+35.0*qSkin[0]*jacobgeo_inv[2]+31.304951684997057*jacobgeo_inv[1]*qSkin[1]); 
  fSkin[8] = 0.04714045207910316*(15.000000000000002*jacobgeo_inv[1]*qSkin[12]+15.0*jacobgeo_inv[0]*qSkin[8]); 
  fSkin[9] = 0.04714045207910316*(15.000000000000002*jacobgeo_inv[1]*qSkin[15]+15.0*jacobgeo_inv[0]*qSkin[9]); 
  fSkin[10] = 0.1414213562373095*(4.47213595499958*jacobgeo_inv[1]*qSkin[17]+(4.47213595499958*jacobgeo_inv[2]+5.0*jacobgeo_inv[0])*qSkin[10]+5.0*jacobgeo_inv[1]*qSkin[6]); 
  fSkin[11] = 0.006734350297014738*((67.0820393249937*jacobgeo_inv[2]+105.0*jacobgeo_inv[0])*qSkin[11]+93.91485505499116*jacobgeo_inv[1]*qSkin[4]+105.00000000000001*jacobgeo_inv[2]*qSkin[2]); 
  fSkin[12] = 0.04714045207910316*((13.416407864998739*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qSkin[12]+15.000000000000002*jacobgeo_inv[1]*qSkin[8]); 
  fSkin[13] = 0.006734350297014738*((67.0820393249937*jacobgeo_inv[2]+105.0*jacobgeo_inv[0])*qSkin[13]+93.91485505499116*jacobgeo_inv[1]*qSkin[5]+105.00000000000001*jacobgeo_inv[2]*qSkin[3]); 
  fSkin[14] = 0.04714045207910316*(15.000000000000002*jacobgeo_inv[1]*qSkin[18]+15.0*jacobgeo_inv[0]*qSkin[14]); 
  fSkin[15] = 0.04714045207910316*((13.416407864998739*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qSkin[15]+15.000000000000002*jacobgeo_inv[1]*qSkin[9]); 
  fSkin[16] = 0.04714045207910316*(15.000000000000002*jacobgeo_inv[1]*qSkin[19]+15.0*jacobgeo_inv[0]*qSkin[16]); 
  fSkin[17] = 0.020203050891044214*((22.3606797749979*jacobgeo_inv[2]+35.0*jacobgeo_inv[0])*qSkin[17]+31.304951684997057*jacobgeo_inv[1]*qSkin[10]+35.0*jacobgeo_inv[2]*qSkin[6]); 
  fSkin[18] = 0.04714045207910316*((13.416407864998739*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qSkin[18]+15.000000000000002*jacobgeo_inv[1]*qSkin[14]); 
  fSkin[19] = 0.04714045207910316*((13.416407864998739*jacobgeo_inv[2]+15.0*jacobgeo_inv[0])*qSkin[19]+15.000000000000002*jacobgeo_inv[1]*qSkin[16]); 

  double boundSurf_incr[20] = {0.0}; 

  if (edge == -1) { 

  boundSurf_incr[1] = 1.5309310892394856*coeff[2]*fSkin[7]-1.185854122563142*coeff[1]*fSkin[7]+0.6846531968814573*coeff[0]*fSkin[7]-1.9764235376052366*fSkin[1]*coeff[2]+1.369306393762915*fSkin[0]*coeff[2]+1.5309310892394856*coeff[1]*fSkin[1]-0.883883476483184*coeff[0]*fSkin[1]-1.060660171779821*fSkin[0]*coeff[1]+0.6123724356957944*coeff[0]*fSkin[0]; 
  boundSurf_incr[4] = 1.5309310892394858*coeff[2]*fSkin[11]-1.1858541225631418*coeff[1]*fSkin[11]+0.6846531968814574*coeff[0]*fSkin[11]-1.9764235376052366*coeff[2]*fSkin[4]+1.5309310892394856*coeff[1]*fSkin[4]-0.883883476483184*coeff[0]*fSkin[4]+1.369306393762915*coeff[2]*fSkin[2]-1.060660171779821*coeff[1]*fSkin[2]+0.6123724356957944*coeff[0]*fSkin[2]; 
  boundSurf_incr[5] = 1.5309310892394858*coeff[2]*fSkin[13]-1.1858541225631418*coeff[1]*fSkin[13]+0.6846531968814574*coeff[0]*fSkin[13]-1.9764235376052366*coeff[2]*fSkin[5]+1.5309310892394856*coeff[1]*fSkin[5]-0.883883476483184*coeff[0]*fSkin[5]+1.369306393762915*coeff[2]*fSkin[3]-1.060660171779821*coeff[1]*fSkin[3]+0.6123724356957944*coeff[0]*fSkin[3]; 
  boundSurf_incr[7] = -(5.929270612815711*coeff[2]*fSkin[7])+4.592793267718456*coeff[1]*fSkin[7]-2.651650429449552*coeff[0]*fSkin[7]+7.654655446197428*fSkin[1]*coeff[2]-5.303300858899105*fSkin[0]*coeff[2]-5.929270612815711*coeff[1]*fSkin[1]+3.4232659844072875*coeff[0]*fSkin[1]+4.107919181288745*fSkin[0]*coeff[1]-2.371708245126284*coeff[0]*fSkin[0]; 
  boundSurf_incr[10] = 1.5309310892394856*coeff[2]*fSkin[17]-1.185854122563142*coeff[1]*fSkin[17]+0.6846531968814573*coeff[0]*fSkin[17]-1.9764235376052366*coeff[2]*fSkin[10]+1.5309310892394856*coeff[1]*fSkin[10]-0.883883476483184*coeff[0]*fSkin[10]+1.369306393762915*coeff[2]*fSkin[6]-1.060660171779821*coeff[1]*fSkin[6]+0.6123724356957944*coeff[0]*fSkin[6]; 
  boundSurf_incr[11] = -(5.929270612815711*coeff[2]*fSkin[11])+4.592793267718456*coeff[1]*fSkin[11]-2.651650429449552*coeff[0]*fSkin[11]+7.65465544619743*coeff[2]*fSkin[4]-5.929270612815709*coeff[1]*fSkin[4]+3.4232659844072866*coeff[0]*fSkin[4]-5.303300858899106*coeff[2]*fSkin[2]+4.107919181288745*coeff[1]*fSkin[2]-2.371708245126284*coeff[0]*fSkin[2]; 
  boundSurf_incr[12] = -(1.9764235376052366*coeff[2]*fSkin[12])+1.5309310892394856*coeff[1]*fSkin[12]-0.883883476483184*coeff[0]*fSkin[12]+1.369306393762915*coeff[2]*fSkin[8]-1.060660171779821*coeff[1]*fSkin[8]+0.6123724356957944*coeff[0]*fSkin[8]; 
  boundSurf_incr[13] = -(5.929270612815711*coeff[2]*fSkin[13])+4.592793267718456*coeff[1]*fSkin[13]-2.651650429449552*coeff[0]*fSkin[13]+7.65465544619743*coeff[2]*fSkin[5]-5.929270612815709*coeff[1]*fSkin[5]+3.4232659844072866*coeff[0]*fSkin[5]-5.303300858899106*coeff[2]*fSkin[3]+4.107919181288745*coeff[1]*fSkin[3]-2.371708245126284*coeff[0]*fSkin[3]; 
  boundSurf_incr[15] = -(1.9764235376052366*coeff[2]*fSkin[15])+1.5309310892394856*coeff[1]*fSkin[15]-0.883883476483184*coeff[0]*fSkin[15]+1.369306393762915*coeff[2]*fSkin[9]-1.060660171779821*coeff[1]*fSkin[9]+0.6123724356957944*coeff[0]*fSkin[9]; 
  boundSurf_incr[17] = -(5.929270612815711*coeff[2]*fSkin[17])+4.592793267718456*coeff[1]*fSkin[17]-2.651650429449552*coeff[0]*fSkin[17]+7.654655446197428*coeff[2]*fSkin[10]-5.929270612815711*coeff[1]*fSkin[10]+3.4232659844072875*coeff[0]*fSkin[10]-5.303300858899105*coeff[2]*fSkin[6]+4.107919181288745*coeff[1]*fSkin[6]-2.371708245126284*coeff[0]*fSkin[6]; 
  boundSurf_incr[18] = -(1.9764235376052366*coeff[2]*fSkin[18])+1.5309310892394856*coeff[1]*fSkin[18]-0.883883476483184*coeff[0]*fSkin[18]+1.369306393762915*coeff[2]*fSkin[14]-1.060660171779821*coeff[1]*fSkin[14]+0.6123724356957944*coeff[0]*fSkin[14]; 
  boundSurf_incr[19] = -(1.9764235376052366*coeff[2]*fSkin[19])+1.5309310892394856*coeff[1]*fSkin[19]-0.883883476483184*coeff[0]*fSkin[19]+1.369306393762915*coeff[2]*fSkin[16]-1.060660171779821*coeff[1]*fSkin[16]+0.6123724356957944*coeff[0]*fSkin[16]; 

  } else { 

  boundSurf_incr[1] = -(1.5309310892394856*coeff[2]*fSkin[7])-1.185854122563142*coeff[1]*fSkin[7]-0.6846531968814573*coeff[0]*fSkin[7]-1.9764235376052366*fSkin[1]*coeff[2]-1.369306393762915*fSkin[0]*coeff[2]-1.5309310892394856*coeff[1]*fSkin[1]-0.883883476483184*coeff[0]*fSkin[1]-1.060660171779821*fSkin[0]*coeff[1]-0.6123724356957944*coeff[0]*fSkin[0]; 
  boundSurf_incr[4] = -(1.5309310892394858*coeff[2]*fSkin[11])-1.1858541225631418*coeff[1]*fSkin[11]-0.6846531968814574*coeff[0]*fSkin[11]-1.9764235376052366*coeff[2]*fSkin[4]-1.5309310892394856*coeff[1]*fSkin[4]-0.883883476483184*coeff[0]*fSkin[4]-1.369306393762915*coeff[2]*fSkin[2]-1.060660171779821*coeff[1]*fSkin[2]-0.6123724356957944*coeff[0]*fSkin[2]; 
  boundSurf_incr[5] = -(1.5309310892394858*coeff[2]*fSkin[13])-1.1858541225631418*coeff[1]*fSkin[13]-0.6846531968814574*coeff[0]*fSkin[13]-1.9764235376052366*coeff[2]*fSkin[5]-1.5309310892394856*coeff[1]*fSkin[5]-0.883883476483184*coeff[0]*fSkin[5]-1.369306393762915*coeff[2]*fSkin[3]-1.060660171779821*coeff[1]*fSkin[3]-0.6123724356957944*coeff[0]*fSkin[3]; 
  boundSurf_incr[7] = -(5.929270612815711*coeff[2]*fSkin[7])-4.592793267718456*coeff[1]*fSkin[7]-2.651650429449552*coeff[0]*fSkin[7]-7.654655446197428*fSkin[1]*coeff[2]-5.303300858899105*fSkin[0]*coeff[2]-5.929270612815711*coeff[1]*fSkin[1]-3.4232659844072875*coeff[0]*fSkin[1]-4.107919181288745*fSkin[0]*coeff[1]-2.371708245126284*coeff[0]*fSkin[0]; 
  boundSurf_incr[10] = -(1.5309310892394856*coeff[2]*fSkin[17])-1.185854122563142*coeff[1]*fSkin[17]-0.6846531968814573*coeff[0]*fSkin[17]-1.9764235376052366*coeff[2]*fSkin[10]-1.5309310892394856*coeff[1]*fSkin[10]-0.883883476483184*coeff[0]*fSkin[10]-1.369306393762915*coeff[2]*fSkin[6]-1.060660171779821*coeff[1]*fSkin[6]-0.6123724356957944*coeff[0]*fSkin[6]; 
  boundSurf_incr[11] = -(5.929270612815711*coeff[2]*fSkin[11])-4.592793267718456*coeff[1]*fSkin[11]-2.651650429449552*coeff[0]*fSkin[11]-7.65465544619743*coeff[2]*fSkin[4]-5.929270612815709*coeff[1]*fSkin[4]-3.4232659844072866*coeff[0]*fSkin[4]-5.303300858899106*coeff[2]*fSkin[2]-4.107919181288745*coeff[1]*fSkin[2]-2.371708245126284*coeff[0]*fSkin[2]; 
  boundSurf_incr[12] = -(1.9764235376052366*coeff[2]*fSkin[12])-1.5309310892394856*coeff[1]*fSkin[12]-0.883883476483184*coeff[0]*fSkin[12]-1.369306393762915*coeff[2]*fSkin[8]-1.060660171779821*coeff[1]*fSkin[8]-0.6123724356957944*coeff[0]*fSkin[8]; 
  boundSurf_incr[13] = -(5.929270612815711*coeff[2]*fSkin[13])-4.592793267718456*coeff[1]*fSkin[13]-2.651650429449552*coeff[0]*fSkin[13]-7.65465544619743*coeff[2]*fSkin[5]-5.929270612815709*coeff[1]*fSkin[5]-3.4232659844072866*coeff[0]*fSkin[5]-5.303300858899106*coeff[2]*fSkin[3]-4.107919181288745*coeff[1]*fSkin[3]-2.371708245126284*coeff[0]*fSkin[3]; 
  boundSurf_incr[15] = -(1.9764235376052366*coeff[2]*fSkin[15])-1.5309310892394856*coeff[1]*fSkin[15]-0.883883476483184*coeff[0]*fSkin[15]-1.369306393762915*coeff[2]*fSkin[9]-1.060660171779821*coeff[1]*fSkin[9]-0.6123724356957944*coeff[0]*fSkin[9]; 
  boundSurf_incr[17] = -(5.929270612815711*coeff[2]*fSkin[17])-4.592793267718456*coeff[1]*fSkin[17]-2.651650429449552*coeff[0]*fSkin[17]-7.654655446197428*coeff[2]*fSkin[10]-5.929270612815711*coeff[1]*fSkin[10]-3.4232659844072875*coeff[0]*fSkin[10]-5.303300858899105*coeff[2]*fSkin[6]-4.107919181288745*coeff[1]*fSkin[6]-2.371708245126284*coeff[0]*fSkin[6]; 
  boundSurf_incr[18] = -(1.9764235376052366*coeff[2]*fSkin[18])-1.5309310892394856*coeff[1]*fSkin[18]-0.883883476483184*coeff[0]*fSkin[18]-1.369306393762915*coeff[2]*fSkin[14]-1.060660171779821*coeff[1]*fSkin[14]-0.6123724356957944*coeff[0]*fSkin[14]; 
  boundSurf_incr[19] = -(1.9764235376052366*coeff[2]*fSkin[19])-1.5309310892394856*coeff[1]*fSkin[19]-0.883883476483184*coeff[0]*fSkin[19]-1.369306393762915*coeff[2]*fSkin[16]-1.060660171779821*coeff[1]*fSkin[16]-0.6123724356957944*coeff[0]*fSkin[16]; 

  }

  out[0] += boundSurf_incr[0]*rdx2Sq; 
  out[1] += boundSurf_incr[1]*rdx2Sq; 
  out[2] += boundSurf_incr[2]*rdx2Sq; 
  out[3] += boundSurf_incr[3]*rdx2Sq; 
  out[4] += boundSurf_incr[4]*rdx2Sq; 
  out[5] += boundSurf_incr[5]*rdx2Sq; 
  out[6] += boundSurf_incr[6]*rdx2Sq; 
  out[7] += boundSurf_incr[7]*rdx2Sq; 
  out[8] += boundSurf_incr[8]*rdx2Sq; 
  out[9] += boundSurf_incr[9]*rdx2Sq; 
  out[10] += boundSurf_incr[10]*rdx2Sq; 
  out[11] += boundSurf_incr[11]*rdx2Sq; 
  out[12] += boundSurf_incr[12]*rdx2Sq; 
  out[13] += boundSurf_incr[13]*rdx2Sq; 
  out[14] += boundSurf_incr[14]*rdx2Sq; 
  out[15] += boundSurf_incr[15]*rdx2Sq; 
  out[16] += boundSurf_incr[16]*rdx2Sq; 
  out[17] += boundSurf_incr[17]*rdx2Sq; 
  out[18] += boundSurf_incr[18]*rdx2Sq; 
  out[19] += boundSurf_incr[19]*rdx2Sq; 

  return 0.;
}

