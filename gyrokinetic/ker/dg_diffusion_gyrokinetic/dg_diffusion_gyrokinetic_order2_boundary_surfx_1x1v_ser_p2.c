#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_surfx_1x1v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qSkin/Edge: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],2.);

  double vol_incr[8] = {0.0}; 
  vol_incr[4] = 6.708203932499369*coeff[0]*qSkin[0]; 
  vol_incr[6] = 6.7082039324993685*coeff[0]*qSkin[2]; 

  double edgeSurf_incr[8] = {0.0}; 
  double boundSurf_incr[8] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(0.6708203932499369*coeff[0]*qSkin[4])+0.6708203932499369*coeff[0]*qEdge[4]-1.190784930203603*coeff[0]*qSkin[1]-1.190784930203603*coeff[0]*qEdge[1]-0.9375*coeff[0]*qSkin[0]+0.9375*coeff[0]*qEdge[0]; 
  edgeSurf_incr[1] = -(1.5855025573536612*coeff[0]*qSkin[4])+0.7382874503707886*coeff[0]*qEdge[4]-2.671875*coeff[0]*qSkin[1]-1.453125*coeff[0]*qEdge[1]-2.0568103339880417*coeff[0]*qSkin[0]+1.1907849302036029*coeff[0]*qEdge[0]; 
  edgeSurf_incr[2] = -(0.6708203932499369*coeff[0]*qSkin[6])+0.6708203932499369*coeff[0]*qEdge[6]-1.190784930203603*coeff[0]*qSkin[3]-1.190784930203603*coeff[0]*qEdge[3]-0.9375*coeff[0]*qSkin[2]+0.9375*coeff[0]*qEdge[2]; 
  edgeSurf_incr[3] = -(1.5855025573536612*coeff[0]*qSkin[6])+0.7382874503707888*coeff[0]*qEdge[6]-2.671875*coeff[0]*qSkin[3]-1.453125*coeff[0]*qEdge[3]-2.0568103339880417*coeff[0]*qSkin[2]+1.1907849302036029*coeff[0]*qEdge[2]; 
  edgeSurf_incr[4] = -(3.140625*coeff[0]*qSkin[4])-0.140625*coeff[0]*qEdge[4]-5.022775277112744*coeff[0]*qSkin[1]-0.3025768239224549*coeff[0]*qEdge[1]-3.7733647120308955*coeff[0]*qSkin[0]+0.4192627457812108*coeff[0]*qEdge[0]; 
  edgeSurf_incr[5] = -(1.190784930203603*coeff[0]*qSkin[7])-1.190784930203603*coeff[0]*qEdge[7]-0.9375*coeff[0]*qSkin[5]+0.9375*coeff[0]*qEdge[5]; 
  edgeSurf_incr[6] = -(3.140625*coeff[0]*qSkin[6])-0.140625*coeff[0]*qEdge[6]-5.022775277112744*coeff[0]*qSkin[3]-0.30257682392245444*coeff[0]*qEdge[3]-3.773364712030894*coeff[0]*qSkin[2]+0.41926274578121053*coeff[0]*qEdge[2]; 
  edgeSurf_incr[7] = -(2.671875*coeff[0]*qSkin[7])-1.453125*coeff[0]*qEdge[7]-2.0568103339880417*coeff[0]*qSkin[5]+1.190784930203603*coeff[0]*qEdge[5]; 

  boundSurf_incr[1] = 0.9682458365518543*coeff[0]*qSkin[4]-1.25*coeff[0]*qSkin[1]+0.8660254037844386*coeff[0]*qSkin[0]; 
  boundSurf_incr[3] = 0.9682458365518543*coeff[0]*qSkin[6]-1.25*coeff[0]*qSkin[3]+0.8660254037844386*coeff[0]*qSkin[2]; 
  boundSurf_incr[4] = -(3.75*coeff[0]*qSkin[4])+4.841229182759272*coeff[0]*qSkin[1]-3.3541019662496847*coeff[0]*qSkin[0]; 
  boundSurf_incr[6] = -(3.75*coeff[0]*qSkin[6])+4.841229182759271*coeff[0]*qSkin[3]-3.3541019662496843*coeff[0]*qSkin[2]; 
  boundSurf_incr[7] = 0.8660254037844387*coeff[0]*qSkin[5]-1.25*coeff[0]*qSkin[7]; 

  } else { 

  edgeSurf_incr[0] = -(0.6708203932499369*coeff[0]*qSkin[4])+0.6708203932499369*coeff[0]*qEdge[4]+1.190784930203603*coeff[0]*qSkin[1]+1.190784930203603*coeff[0]*qEdge[1]-0.9375*coeff[0]*qSkin[0]+0.9375*coeff[0]*qEdge[0]; 
  edgeSurf_incr[1] = 1.5855025573536612*coeff[0]*qSkin[4]-0.7382874503707886*coeff[0]*qEdge[4]-2.671875*coeff[0]*qSkin[1]-1.453125*coeff[0]*qEdge[1]+2.0568103339880417*coeff[0]*qSkin[0]-1.1907849302036029*coeff[0]*qEdge[0]; 
  edgeSurf_incr[2] = -(0.6708203932499369*coeff[0]*qSkin[6])+0.6708203932499369*coeff[0]*qEdge[6]+1.190784930203603*coeff[0]*qSkin[3]+1.190784930203603*coeff[0]*qEdge[3]-0.9375*coeff[0]*qSkin[2]+0.9375*coeff[0]*qEdge[2]; 
  edgeSurf_incr[3] = 1.5855025573536612*coeff[0]*qSkin[6]-0.7382874503707888*coeff[0]*qEdge[6]-2.671875*coeff[0]*qSkin[3]-1.453125*coeff[0]*qEdge[3]+2.0568103339880417*coeff[0]*qSkin[2]-1.1907849302036029*coeff[0]*qEdge[2]; 
  edgeSurf_incr[4] = -(3.140625*coeff[0]*qSkin[4])-0.140625*coeff[0]*qEdge[4]+5.022775277112744*coeff[0]*qSkin[1]+0.3025768239224549*coeff[0]*qEdge[1]-3.7733647120308955*coeff[0]*qSkin[0]+0.4192627457812108*coeff[0]*qEdge[0]; 
  edgeSurf_incr[5] = 1.190784930203603*coeff[0]*qSkin[7]+1.190784930203603*coeff[0]*qEdge[7]-0.9375*coeff[0]*qSkin[5]+0.9375*coeff[0]*qEdge[5]; 
  edgeSurf_incr[6] = -(3.140625*coeff[0]*qSkin[6])-0.140625*coeff[0]*qEdge[6]+5.022775277112744*coeff[0]*qSkin[3]+0.30257682392245444*coeff[0]*qEdge[3]-3.773364712030894*coeff[0]*qSkin[2]+0.41926274578121053*coeff[0]*qEdge[2]; 
  edgeSurf_incr[7] = -(2.671875*coeff[0]*qSkin[7])-1.453125*coeff[0]*qEdge[7]+2.0568103339880417*coeff[0]*qSkin[5]-1.190784930203603*coeff[0]*qEdge[5]; 

  boundSurf_incr[1] = -(0.9682458365518543*coeff[0]*qSkin[4])-1.25*coeff[0]*qSkin[1]-0.8660254037844386*coeff[0]*qSkin[0]; 
  boundSurf_incr[3] = -(0.9682458365518543*coeff[0]*qSkin[6])-1.25*coeff[0]*qSkin[3]-0.8660254037844386*coeff[0]*qSkin[2]; 
  boundSurf_incr[4] = -(3.75*coeff[0]*qSkin[4])-4.841229182759272*coeff[0]*qSkin[1]-3.3541019662496847*coeff[0]*qSkin[0]; 
  boundSurf_incr[6] = -(3.75*coeff[0]*qSkin[6])-4.841229182759271*coeff[0]*qSkin[3]-3.3541019662496843*coeff[0]*qSkin[2]; 
  boundSurf_incr[7] = -(1.25*coeff[0]*qSkin[7])-0.8660254037844387*coeff[0]*qSkin[5]; 

  }

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*rdx2Sq; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*rdx2Sq; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*rdx2Sq; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*rdx2Sq; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*rdx2Sq; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*rdx2Sq; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*rdx2Sq; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*rdx2Sq; 

  return 0.;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_surfx_1x1v_ser_p2_varcoeff(const double *w, const double *dx, const double *coeff, const double *jacobgeo_inv, int edge, const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // jacobgeo_inv: one divided by the configuration space Jacobian.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // qSkin/Edge: scalar field in skin and egde cells.
  // out: Incremented output.

  const double rdx2Sq = pow(2./dx[0],2.);

  double fSkin[8];
  fSkin[0] = 0.7071067811865476*(jacobgeo_inv[2]*qSkin[4]+jacobgeo_inv[1]*qSkin[1]+jacobgeo_inv[0]*qSkin[0]); 
  fSkin[1] = 0.1*(6.324555320336761*(jacobgeo_inv[1]*qSkin[4]+qSkin[1]*jacobgeo_inv[2])+7.0710678118654755*(jacobgeo_inv[0]*qSkin[1]+qSkin[0]*jacobgeo_inv[1])); 
  fSkin[2] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[2]*qSkin[6]+21.213203435596427*(jacobgeo_inv[1]*qSkin[3]+jacobgeo_inv[0]*qSkin[2])); 
  fSkin[3] = 0.03333333333333333*(18.97366596101028*jacobgeo_inv[1]*qSkin[6]+(18.97366596101028*jacobgeo_inv[2]+21.213203435596427*jacobgeo_inv[0])*qSkin[3]+21.213203435596427*jacobgeo_inv[1]*qSkin[2]); 
  fSkin[4] = 0.014285714285714285*((31.622776601683807*jacobgeo_inv[2]+49.49747468305833*jacobgeo_inv[0])*qSkin[4]+49.49747468305833*qSkin[0]*jacobgeo_inv[2]+44.27188724235732*jacobgeo_inv[1]*qSkin[1]); 
  fSkin[5] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[1]*qSkin[7]+21.213203435596427*jacobgeo_inv[0]*qSkin[5]); 
  fSkin[6] = 0.004761904761904762*((94.86832980505142*jacobgeo_inv[2]+148.49242404917499*jacobgeo_inv[0])*qSkin[6]+132.81566172707196*jacobgeo_inv[1]*qSkin[3]+148.49242404917499*jacobgeo_inv[2]*qSkin[2]); 
  fSkin[7] = 0.03333333333333333*((18.97366596101028*jacobgeo_inv[2]+21.213203435596427*jacobgeo_inv[0])*qSkin[7]+21.21320343559643*jacobgeo_inv[1]*qSkin[5]); 

  double fEdge[8];
  fEdge[0] = 0.7071067811865476*(jacobgeo_inv[2]*qEdge[4]+jacobgeo_inv[1]*qEdge[1]+jacobgeo_inv[0]*qEdge[0]); 
  fEdge[1] = 0.1*(6.324555320336761*(jacobgeo_inv[1]*qEdge[4]+qEdge[1]*jacobgeo_inv[2])+7.0710678118654755*(jacobgeo_inv[0]*qEdge[1]+qEdge[0]*jacobgeo_inv[1])); 
  fEdge[2] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[2]*qEdge[6]+21.213203435596427*(jacobgeo_inv[1]*qEdge[3]+jacobgeo_inv[0]*qEdge[2])); 
  fEdge[3] = 0.03333333333333333*(18.97366596101028*jacobgeo_inv[1]*qEdge[6]+(18.97366596101028*jacobgeo_inv[2]+21.213203435596427*jacobgeo_inv[0])*qEdge[3]+21.213203435596427*jacobgeo_inv[1]*qEdge[2]); 
  fEdge[4] = 0.014285714285714285*((31.622776601683807*jacobgeo_inv[2]+49.49747468305833*jacobgeo_inv[0])*qEdge[4]+49.49747468305833*qEdge[0]*jacobgeo_inv[2]+44.27188724235732*jacobgeo_inv[1]*qEdge[1]); 
  fEdge[5] = 0.03333333333333333*(21.21320343559643*jacobgeo_inv[1]*qEdge[7]+21.213203435596427*jacobgeo_inv[0]*qEdge[5]); 
  fEdge[6] = 0.004761904761904762*((94.86832980505142*jacobgeo_inv[2]+148.49242404917499*jacobgeo_inv[0])*qEdge[6]+132.81566172707196*jacobgeo_inv[1]*qEdge[3]+148.49242404917499*jacobgeo_inv[2]*qEdge[2]); 
  fEdge[7] = 0.03333333333333333*((18.97366596101028*jacobgeo_inv[2]+21.213203435596427*jacobgeo_inv[0])*qEdge[7]+21.21320343559643*jacobgeo_inv[1]*qEdge[5]); 

  double vol_incr[8] = {0.0}; 
  vol_incr[1] = 4.743416490252569*fSkin[1]*coeff[2]+2.1213203435596424*fSkin[0]*coeff[1]; 
  vol_incr[3] = 4.743416490252569*coeff[2]*fSkin[3]+2.1213203435596424*coeff[1]*fSkin[2]; 
  vol_incr[4] = 14.230249470757707*coeff[2]*fSkin[4]+10.606601717798211*fSkin[0]*coeff[2]+9.48683298050514*coeff[1]*fSkin[1]+4.743416490252569*coeff[0]*fSkin[0]; 
  vol_incr[6] = 14.230249470757707*coeff[2]*fSkin[6]+9.48683298050514*coeff[1]*fSkin[3]+10.606601717798213*coeff[2]*fSkin[2]+4.743416490252569*coeff[0]*fSkin[2]; 
  vol_incr[7] = 4.743416490252569*coeff[2]*fSkin[7]+2.1213203435596424*coeff[1]*fSkin[5]; 

  double edgeSurf_incr[8] = {0.0}; 
  double boundSurf_incr[8] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 0.5303300858899105*coeff[2]*fSkin[4]-0.4743416490252568*coeff[0]*fSkin[4]-0.5303300858899105*coeff[2]*fEdge[4]+0.4743416490252568*coeff[0]*fEdge[4]+0.9413981457120035*fSkin[1]*coeff[2]+0.9413981457120035*fEdge[1]*coeff[2]+0.7411588266019635*fSkin[0]*coeff[2]-0.7411588266019635*fEdge[0]*coeff[2]-0.8420120990817169*coeff[0]*fSkin[1]-0.8420120990817169*coeff[0]*fEdge[1]-0.6629126073623879*coeff[0]*fSkin[0]+0.6629126073623879*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 0.24877630200141654*coeff[2]*fSkin[4]-0.5188111786213743*coeff[1]*fSkin[4]-1.1211196098933862*coeff[0]*fSkin[4]-1.588341005085966*coeff[2]*fEdge[4]-0.5188111786213743*coeff[1]*fEdge[4]+0.5220480626221115*coeff[0]*fEdge[4]+0.6670429439417671*fSkin[1]*coeff[2]+2.594055893106872*fEdge[1]*coeff[2]+0.5990715472712751*fSkin[0]*coeff[2]-1.9683779410341897*fEdge[0]*coeff[2]-0.7463289060042488*coeff[1]*fSkin[1]-1.8893009309828055*coeff[0]*fSkin[1]+0.7463289060042488*coeff[1]*fEdge[1]-1.0275145414117013*coeff[0]*fEdge[1]-0.5303300858899105*fSkin[0]*coeff[1]-0.5303300858899105*fEdge[0]*coeff[1]-1.454384534777511*coeff[0]*fSkin[0]+0.8420120990817168*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 0.5303300858899104*coeff[2]*fSkin[6]-0.4743416490252568*coeff[0]*fSkin[6]-0.5303300858899104*coeff[2]*fEdge[6]+0.4743416490252568*coeff[0]*fEdge[6]+0.9413981457120035*coeff[2]*fSkin[3]-0.8420120990817169*coeff[0]*fSkin[3]+0.9413981457120035*coeff[2]*fEdge[3]-0.8420120990817169*coeff[0]*fEdge[3]+0.7411588266019635*coeff[2]*fSkin[2]-0.6629126073623879*coeff[0]*fSkin[2]-0.7411588266019635*coeff[2]*fEdge[2]+0.6629126073623879*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 0.24877630200141687*coeff[2]*fSkin[6]-0.5188111786213743*coeff[1]*fSkin[6]-1.1211196098933864*coeff[0]*fSkin[6]-1.5883410050859663*coeff[2]*fEdge[6]-0.5188111786213743*coeff[1]*fEdge[6]+0.5220480626221116*coeff[0]*fEdge[6]+0.6670429439417671*coeff[2]*fSkin[3]-0.7463289060042488*coeff[1]*fSkin[3]-1.8893009309828055*coeff[0]*fSkin[3]+2.594055893106872*coeff[2]*fEdge[3]+0.7463289060042488*coeff[1]*fEdge[3]-1.0275145414117013*coeff[0]*fEdge[3]+0.5990715472712751*coeff[2]*fSkin[2]-0.5303300858899105*coeff[1]*fSkin[2]-1.454384534777511*coeff[0]*fSkin[2]-1.9683779410341897*coeff[2]*fEdge[2]-0.5303300858899105*coeff[1]*fEdge[2]+0.8420120990817168*coeff[0]*fEdge[2]; 
  edgeSurf_incr[4] = -(1.40820177054373*coeff[2]*fSkin[4])-2.009347054626824*coeff[1]*fSkin[4]-2.220757234663999*coeff[0]*fSkin[4]-3.779910015670014*coeff[2]*fEdge[4]-2.009347054626824*coeff[1]*fEdge[4]-0.0994368911043575*coeff[0]*fEdge[4]-1.626614282316952*fSkin[1]*coeff[2]+5.836674777725536*fEdge[1]*coeff[2]-0.9943689110435825*fSkin[0]*coeff[2]-4.3089319478555215*fEdge[0]*coeff[2]-2.890519423747657*coeff[1]*fSkin[1]-3.5516384588225587*coeff[0]*fSkin[1]+2.890519423747657*coeff[1]*fEdge[1]-0.21395412402545588*coeff[0]*fEdge[1]-2.053959590644372*fSkin[0]*coeff[1]-2.053959590644372*fEdge[0]*coeff[1]-2.6681717757670693*coeff[0]*fSkin[0]+0.29646353064078523*coeff[0]*fEdge[0]; 
  edgeSurf_incr[5] = 0.9413981457120036*coeff[2]*fSkin[7]-0.842012099081717*coeff[0]*fSkin[7]+0.9413981457120036*coeff[2]*fEdge[7]-0.842012099081717*coeff[0]*fEdge[7]+0.7411588266019635*coeff[2]*fSkin[5]-0.6629126073623879*coeff[0]*fSkin[5]-0.7411588266019635*coeff[2]*fEdge[5]+0.6629126073623879*coeff[0]*fEdge[5]; 
  edgeSurf_incr[6] = -(1.40820177054373*coeff[2]*fSkin[6])-2.009347054626824*coeff[1]*fSkin[6]-2.220757234663999*coeff[0]*fSkin[6]-3.779910015670014*coeff[2]*fEdge[6]-2.009347054626824*coeff[1]*fEdge[6]-0.0994368911043575*coeff[0]*fEdge[6]-1.6266142823169534*coeff[2]*fSkin[3]-2.8905194237476564*coeff[1]*fSkin[3]-3.551638458822559*coeff[0]*fSkin[3]+5.836674777725538*coeff[2]*fEdge[3]+2.8905194237476564*coeff[1]*fEdge[3]-0.21395412402545588*coeff[0]*fEdge[3]-0.9943689110435823*coeff[2]*fSkin[2]-2.0539595906443724*coeff[1]*fSkin[2]-2.668171775767069*coeff[0]*fSkin[2]-4.3089319478555215*coeff[2]*fEdge[2]-2.0539595906443724*coeff[1]*fEdge[2]+0.29646353064078523*coeff[0]*fEdge[2]; 
  edgeSurf_incr[7] = 0.6670429439417671*coeff[2]*fSkin[7]-0.7463289060042488*coeff[1]*fSkin[7]-1.8893009309828055*coeff[0]*fSkin[7]+2.594055893106872*coeff[2]*fEdge[7]+0.7463289060042488*coeff[1]*fEdge[7]-1.0275145414117013*coeff[0]*fEdge[7]+0.599071547271275*coeff[2]*fSkin[5]-0.5303300858899104*coeff[1]*fSkin[5]-1.4543845347775113*coeff[0]*fSkin[5]-1.96837794103419*coeff[2]*fEdge[5]-0.5303300858899104*coeff[1]*fEdge[5]+0.842012099081717*coeff[0]*fEdge[5]; 

  boundSurf_incr[1] = 1.5309310892394856*coeff[2]*fSkin[4]-1.185854122563142*coeff[1]*fSkin[4]+0.6846531968814573*coeff[0]*fSkin[4]-1.9764235376052366*fSkin[1]*coeff[2]+1.369306393762915*fSkin[0]*coeff[2]+1.5309310892394856*coeff[1]*fSkin[1]-0.883883476483184*coeff[0]*fSkin[1]-1.060660171779821*fSkin[0]*coeff[1]+0.6123724356957944*coeff[0]*fSkin[0]; 
  boundSurf_incr[3] = 1.5309310892394858*coeff[2]*fSkin[6]-1.1858541225631418*coeff[1]*fSkin[6]+0.6846531968814574*coeff[0]*fSkin[6]-1.9764235376052366*coeff[2]*fSkin[3]+1.5309310892394856*coeff[1]*fSkin[3]-0.883883476483184*coeff[0]*fSkin[3]+1.369306393762915*coeff[2]*fSkin[2]-1.060660171779821*coeff[1]*fSkin[2]+0.6123724356957944*coeff[0]*fSkin[2]; 
  boundSurf_incr[4] = -(5.929270612815711*coeff[2]*fSkin[4])+4.592793267718456*coeff[1]*fSkin[4]-2.651650429449552*coeff[0]*fSkin[4]+7.654655446197428*fSkin[1]*coeff[2]-5.303300858899105*fSkin[0]*coeff[2]-5.929270612815711*coeff[1]*fSkin[1]+3.4232659844072875*coeff[0]*fSkin[1]+4.107919181288745*fSkin[0]*coeff[1]-2.371708245126284*coeff[0]*fSkin[0]; 
  boundSurf_incr[6] = -(5.929270612815711*coeff[2]*fSkin[6])+4.592793267718456*coeff[1]*fSkin[6]-2.651650429449552*coeff[0]*fSkin[6]+7.65465544619743*coeff[2]*fSkin[3]-5.929270612815709*coeff[1]*fSkin[3]+3.4232659844072866*coeff[0]*fSkin[3]-5.303300858899106*coeff[2]*fSkin[2]+4.107919181288745*coeff[1]*fSkin[2]-2.371708245126284*coeff[0]*fSkin[2]; 
  boundSurf_incr[7] = -(1.9764235376052366*coeff[2]*fSkin[7])+1.5309310892394856*coeff[1]*fSkin[7]-0.883883476483184*coeff[0]*fSkin[7]+1.369306393762915*coeff[2]*fSkin[5]-1.060660171779821*coeff[1]*fSkin[5]+0.6123724356957944*coeff[0]*fSkin[5]; 

  } else { 

  edgeSurf_incr[0] = 0.5303300858899105*coeff[2]*fSkin[4]-0.4743416490252568*coeff[0]*fSkin[4]-0.5303300858899105*coeff[2]*fEdge[4]+0.4743416490252568*coeff[0]*fEdge[4]-0.9413981457120035*fSkin[1]*coeff[2]-0.9413981457120035*fEdge[1]*coeff[2]+0.7411588266019635*fSkin[0]*coeff[2]-0.7411588266019635*fEdge[0]*coeff[2]+0.8420120990817169*coeff[0]*fSkin[1]+0.8420120990817169*coeff[0]*fEdge[1]-0.6629126073623879*coeff[0]*fSkin[0]+0.6629126073623879*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = -(0.24877630200141654*coeff[2]*fSkin[4])-0.5188111786213743*coeff[1]*fSkin[4]+1.1211196098933862*coeff[0]*fSkin[4]+1.588341005085966*coeff[2]*fEdge[4]-0.5188111786213743*coeff[1]*fEdge[4]-0.5220480626221115*coeff[0]*fEdge[4]+0.6670429439417671*fSkin[1]*coeff[2]+2.594055893106872*fEdge[1]*coeff[2]-0.5990715472712751*fSkin[0]*coeff[2]+1.9683779410341897*fEdge[0]*coeff[2]+0.7463289060042488*coeff[1]*fSkin[1]-1.8893009309828055*coeff[0]*fSkin[1]-0.7463289060042488*coeff[1]*fEdge[1]-1.0275145414117013*coeff[0]*fEdge[1]-0.5303300858899105*fSkin[0]*coeff[1]-0.5303300858899105*fEdge[0]*coeff[1]+1.454384534777511*coeff[0]*fSkin[0]-0.8420120990817168*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 0.5303300858899104*coeff[2]*fSkin[6]-0.4743416490252568*coeff[0]*fSkin[6]-0.5303300858899104*coeff[2]*fEdge[6]+0.4743416490252568*coeff[0]*fEdge[6]-0.9413981457120035*coeff[2]*fSkin[3]+0.8420120990817169*coeff[0]*fSkin[3]-0.9413981457120035*coeff[2]*fEdge[3]+0.8420120990817169*coeff[0]*fEdge[3]+0.7411588266019635*coeff[2]*fSkin[2]-0.6629126073623879*coeff[0]*fSkin[2]-0.7411588266019635*coeff[2]*fEdge[2]+0.6629126073623879*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = -(0.24877630200141687*coeff[2]*fSkin[6])-0.5188111786213743*coeff[1]*fSkin[6]+1.1211196098933864*coeff[0]*fSkin[6]+1.5883410050859663*coeff[2]*fEdge[6]-0.5188111786213743*coeff[1]*fEdge[6]-0.5220480626221116*coeff[0]*fEdge[6]+0.6670429439417671*coeff[2]*fSkin[3]+0.7463289060042488*coeff[1]*fSkin[3]-1.8893009309828055*coeff[0]*fSkin[3]+2.594055893106872*coeff[2]*fEdge[3]-0.7463289060042488*coeff[1]*fEdge[3]-1.0275145414117013*coeff[0]*fEdge[3]-0.5990715472712751*coeff[2]*fSkin[2]-0.5303300858899105*coeff[1]*fSkin[2]+1.454384534777511*coeff[0]*fSkin[2]+1.9683779410341897*coeff[2]*fEdge[2]-0.5303300858899105*coeff[1]*fEdge[2]-0.8420120990817168*coeff[0]*fEdge[2]; 
  edgeSurf_incr[4] = -(1.40820177054373*coeff[2]*fSkin[4])+2.009347054626824*coeff[1]*fSkin[4]-2.220757234663999*coeff[0]*fSkin[4]-3.779910015670014*coeff[2]*fEdge[4]+2.009347054626824*coeff[1]*fEdge[4]-0.0994368911043575*coeff[0]*fEdge[4]+1.626614282316952*fSkin[1]*coeff[2]-5.836674777725536*fEdge[1]*coeff[2]-0.9943689110435825*fSkin[0]*coeff[2]-4.3089319478555215*fEdge[0]*coeff[2]-2.890519423747657*coeff[1]*fSkin[1]+3.5516384588225587*coeff[0]*fSkin[1]+2.890519423747657*coeff[1]*fEdge[1]+0.21395412402545588*coeff[0]*fEdge[1]+2.053959590644372*fSkin[0]*coeff[1]+2.053959590644372*fEdge[0]*coeff[1]-2.6681717757670693*coeff[0]*fSkin[0]+0.29646353064078523*coeff[0]*fEdge[0]; 
  edgeSurf_incr[5] = -(0.9413981457120036*coeff[2]*fSkin[7])+0.842012099081717*coeff[0]*fSkin[7]-0.9413981457120036*coeff[2]*fEdge[7]+0.842012099081717*coeff[0]*fEdge[7]+0.7411588266019635*coeff[2]*fSkin[5]-0.6629126073623879*coeff[0]*fSkin[5]-0.7411588266019635*coeff[2]*fEdge[5]+0.6629126073623879*coeff[0]*fEdge[5]; 
  edgeSurf_incr[6] = -(1.40820177054373*coeff[2]*fSkin[6])+2.009347054626824*coeff[1]*fSkin[6]-2.220757234663999*coeff[0]*fSkin[6]-3.779910015670014*coeff[2]*fEdge[6]+2.009347054626824*coeff[1]*fEdge[6]-0.0994368911043575*coeff[0]*fEdge[6]+1.6266142823169534*coeff[2]*fSkin[3]-2.8905194237476564*coeff[1]*fSkin[3]+3.551638458822559*coeff[0]*fSkin[3]-5.836674777725538*coeff[2]*fEdge[3]+2.8905194237476564*coeff[1]*fEdge[3]+0.21395412402545588*coeff[0]*fEdge[3]-0.9943689110435823*coeff[2]*fSkin[2]+2.0539595906443724*coeff[1]*fSkin[2]-2.668171775767069*coeff[0]*fSkin[2]-4.3089319478555215*coeff[2]*fEdge[2]+2.0539595906443724*coeff[1]*fEdge[2]+0.29646353064078523*coeff[0]*fEdge[2]; 
  edgeSurf_incr[7] = 0.6670429439417671*coeff[2]*fSkin[7]+0.7463289060042488*coeff[1]*fSkin[7]-1.8893009309828055*coeff[0]*fSkin[7]+2.594055893106872*coeff[2]*fEdge[7]-0.7463289060042488*coeff[1]*fEdge[7]-1.0275145414117013*coeff[0]*fEdge[7]-0.599071547271275*coeff[2]*fSkin[5]-0.5303300858899104*coeff[1]*fSkin[5]+1.4543845347775113*coeff[0]*fSkin[5]+1.96837794103419*coeff[2]*fEdge[5]-0.5303300858899104*coeff[1]*fEdge[5]-0.842012099081717*coeff[0]*fEdge[5]; 

  boundSurf_incr[1] = -(1.5309310892394856*coeff[2]*fSkin[4])-1.185854122563142*coeff[1]*fSkin[4]-0.6846531968814573*coeff[0]*fSkin[4]-1.9764235376052366*fSkin[1]*coeff[2]-1.369306393762915*fSkin[0]*coeff[2]-1.5309310892394856*coeff[1]*fSkin[1]-0.883883476483184*coeff[0]*fSkin[1]-1.060660171779821*fSkin[0]*coeff[1]-0.6123724356957944*coeff[0]*fSkin[0]; 
  boundSurf_incr[3] = -(1.5309310892394858*coeff[2]*fSkin[6])-1.1858541225631418*coeff[1]*fSkin[6]-0.6846531968814574*coeff[0]*fSkin[6]-1.9764235376052366*coeff[2]*fSkin[3]-1.5309310892394856*coeff[1]*fSkin[3]-0.883883476483184*coeff[0]*fSkin[3]-1.369306393762915*coeff[2]*fSkin[2]-1.060660171779821*coeff[1]*fSkin[2]-0.6123724356957944*coeff[0]*fSkin[2]; 
  boundSurf_incr[4] = -(5.929270612815711*coeff[2]*fSkin[4])-4.592793267718456*coeff[1]*fSkin[4]-2.651650429449552*coeff[0]*fSkin[4]-7.654655446197428*fSkin[1]*coeff[2]-5.303300858899105*fSkin[0]*coeff[2]-5.929270612815711*coeff[1]*fSkin[1]-3.4232659844072875*coeff[0]*fSkin[1]-4.107919181288745*fSkin[0]*coeff[1]-2.371708245126284*coeff[0]*fSkin[0]; 
  boundSurf_incr[6] = -(5.929270612815711*coeff[2]*fSkin[6])-4.592793267718456*coeff[1]*fSkin[6]-2.651650429449552*coeff[0]*fSkin[6]-7.65465544619743*coeff[2]*fSkin[3]-5.929270612815709*coeff[1]*fSkin[3]-3.4232659844072866*coeff[0]*fSkin[3]-5.303300858899106*coeff[2]*fSkin[2]-4.107919181288745*coeff[1]*fSkin[2]-2.371708245126284*coeff[0]*fSkin[2]; 
  boundSurf_incr[7] = -(1.9764235376052366*coeff[2]*fSkin[7])-1.5309310892394856*coeff[1]*fSkin[7]-0.883883476483184*coeff[0]*fSkin[7]-1.369306393762915*coeff[2]*fSkin[5]-1.060660171779821*coeff[1]*fSkin[5]-0.6123724356957944*coeff[0]*fSkin[5]; 

  }

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*rdx2Sq; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*rdx2Sq; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*rdx2Sq; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*rdx2Sq; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*rdx2Sq; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*rdx2Sq; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*rdx2Sq; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*rdx2Sq; 

  return 0.;
}

