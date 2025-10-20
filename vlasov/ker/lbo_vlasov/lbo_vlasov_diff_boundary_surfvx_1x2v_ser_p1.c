#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[3]: Cell-center coordinates. 
  // dxv[3]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[6]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fSkin/Edge: Distribution function in cells 
  // out: Incremented distribution function in cell 
  const double *nuVtSqSum = &nuPrimMomsSum[4];

  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 

  double vol_incr[16] = {0.0}; 
  vol_incr[8] = 4.743416490252569*fSkin[1]*nuVtSqSum[1]+4.743416490252569*fSkin[0]*nuVtSqSum[0]; 
  vol_incr[9] = 4.743416490252569*fSkin[0]*nuVtSqSum[1]+4.743416490252569*nuVtSqSum[0]*fSkin[1]; 
  vol_incr[10] = 4.743416490252569*nuVtSqSum[1]*fSkin[5]+4.743416490252569*nuVtSqSum[0]*fSkin[3]; 
  vol_incr[11] = 4.743416490252569*nuVtSqSum[0]*fSkin[5]+4.743416490252569*nuVtSqSum[1]*fSkin[3]; 

  double edgeSurf_incr[16] = {0.0}; 
  double boundSurf_incr[16] = {0.0}; 

  if (edge == -1) { 

  double edgeSurf[16] = {0.0}; 
  edgeSurf[0] = (-0.6708203932499369*fSkin[8])+0.6708203932499369*fEdge[8]-1.190784930203603*fSkin[2]-1.190784930203603*fEdge[2]-0.9375*fSkin[0]+0.9375*fEdge[0]; 
  edgeSurf[1] = (-0.6708203932499369*fSkin[9])+0.6708203932499369*fEdge[9]-1.190784930203603*fSkin[4]-1.190784930203603*fEdge[4]-0.9375*fSkin[1]+0.9375*fEdge[1]; 
  edgeSurf[2] = (-1.585502557353661*fSkin[8])+0.7382874503707888*fEdge[8]-2.671875*fSkin[2]-1.453125*fEdge[2]-2.056810333988042*fSkin[0]+1.190784930203603*fEdge[0]; 
  edgeSurf[3] = (-0.6708203932499369*fSkin[10])+0.6708203932499369*fEdge[10]-1.190784930203603*fSkin[6]-1.190784930203603*fEdge[6]-0.9375*fSkin[3]+0.9375*fEdge[3]; 
  edgeSurf[4] = (-1.585502557353661*fSkin[9])+0.7382874503707888*fEdge[9]-2.671875*fSkin[4]-1.453125*fEdge[4]-2.056810333988042*fSkin[1]+1.190784930203603*fEdge[1]; 
  edgeSurf[5] = (-0.6708203932499369*fSkin[11])+0.6708203932499369*fEdge[11]-1.190784930203603*fSkin[7]-1.190784930203603*fEdge[7]-0.9375*fSkin[5]+0.9375*fEdge[5]; 
  edgeSurf[6] = (-1.585502557353661*fSkin[10])+0.7382874503707888*fEdge[10]-2.671875*fSkin[6]-1.453125*fEdge[6]-2.056810333988042*fSkin[3]+1.190784930203603*fEdge[3]; 
  edgeSurf[7] = (-1.585502557353661*fSkin[11])+0.7382874503707888*fEdge[11]-2.671875*fSkin[7]-1.453125*fEdge[7]-2.056810333988042*fSkin[5]+1.190784930203603*fEdge[5]; 
  edgeSurf[8] = (-3.140625*fSkin[8])-0.140625*fEdge[8]-5.022775277112744*fSkin[2]-0.3025768239224545*fEdge[2]-3.773364712030896*fSkin[0]+0.4192627457812106*fEdge[0]; 
  edgeSurf[9] = (-3.140625*fSkin[9])-0.140625*fEdge[9]-5.022775277112744*fSkin[4]-0.3025768239224544*fEdge[4]-3.773364712030894*fSkin[1]+0.4192627457812105*fEdge[1]; 
  edgeSurf[10] = (-3.140625*fSkin[10])-0.140625*fEdge[10]-5.022775277112744*fSkin[6]-0.3025768239224544*fEdge[6]-3.773364712030894*fSkin[3]+0.4192627457812105*fEdge[3]; 
  edgeSurf[11] = (-3.140625*fSkin[11])-0.140625*fEdge[11]-5.022775277112744*fSkin[7]-0.3025768239224545*fEdge[7]-3.773364712030896*fSkin[5]+0.4192627457812106*fEdge[5]; 
  edgeSurf[12] = (-1.190784930203603*fSkin[14])-1.190784930203603*fEdge[14]-0.9375*fSkin[12]+0.9375*fEdge[12]; 
  edgeSurf[13] = (-1.190784930203603*fSkin[15])-1.190784930203603*fEdge[15]-0.9375*fSkin[13]+0.9375*fEdge[13]; 
  edgeSurf[14] = (-2.671875*fSkin[14])-1.453125*fEdge[14]-2.056810333988042*fSkin[12]+1.190784930203603*fEdge[12]; 
  edgeSurf[15] = (-2.671875*fSkin[15])-1.453125*fEdge[15]-2.056810333988042*fSkin[13]+1.190784930203603*fEdge[13]; 

  double boundSurf[16] = {0.0}; 
  boundSurf[2] = 1.936491673103709*fSkin[8]-1.5*fSkin[2]+0.8660254037844386*fSkin[0]; 
  boundSurf[4] = 1.936491673103709*fSkin[9]-1.5*fSkin[4]+0.8660254037844386*fSkin[1]; 
  boundSurf[6] = 1.936491673103709*fSkin[10]-1.5*fSkin[6]+0.8660254037844386*fSkin[3]; 
  boundSurf[7] = 1.936491673103709*fSkin[11]-1.5*fSkin[7]+0.8660254037844386*fSkin[5]; 
  boundSurf[8] = (-7.5*fSkin[8])+5.809475019311125*fSkin[2]-3.354101966249685*fSkin[0]; 
  boundSurf[9] = (-7.5*fSkin[9])+5.809475019311126*fSkin[4]-3.354101966249684*fSkin[1]; 
  boundSurf[10] = (-7.5*fSkin[10])+5.809475019311126*fSkin[6]-3.354101966249684*fSkin[3]; 
  boundSurf[11] = (-7.5*fSkin[11])+5.809475019311125*fSkin[7]-3.354101966249685*fSkin[5]; 
  boundSurf[14] = 0.8660254037844387*fSkin[12]-1.5*fSkin[14]; 
  boundSurf[15] = 0.8660254037844387*fSkin[13]-1.5*fSkin[15]; 

  edgeSurf_incr[0] = 0.7071067811865475*edgeSurf[1]*nuVtSqSum[1]+0.7071067811865475*edgeSurf[0]*nuVtSqSum[0]; 
  edgeSurf_incr[1] = 0.7071067811865475*edgeSurf[0]*nuVtSqSum[1]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[1]; 
  edgeSurf_incr[2] = 0.7071067811865475*nuVtSqSum[1]*edgeSurf[4]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[2]; 
  edgeSurf_incr[3] = 0.7071067811865475*nuVtSqSum[1]*edgeSurf[5]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[3]; 
  edgeSurf_incr[4] = 0.7071067811865475*nuVtSqSum[0]*edgeSurf[4]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[2]; 
  edgeSurf_incr[5] = 0.7071067811865475*nuVtSqSum[0]*edgeSurf[5]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[3]; 
  edgeSurf_incr[6] = 0.7071067811865475*nuVtSqSum[1]*edgeSurf[7]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[6]; 
  edgeSurf_incr[7] = 0.7071067811865475*nuVtSqSum[0]*edgeSurf[7]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[6]; 
  edgeSurf_incr[8] = 0.7071067811865475*nuVtSqSum[1]*edgeSurf[9]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[8]; 
  edgeSurf_incr[9] = 0.7071067811865475*nuVtSqSum[0]*edgeSurf[9]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[8]; 
  edgeSurf_incr[10] = 0.7071067811865475*nuVtSqSum[1]*edgeSurf[11]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[10]; 
  edgeSurf_incr[11] = 0.7071067811865475*nuVtSqSum[0]*edgeSurf[11]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[10]; 
  edgeSurf_incr[12] = 0.7071067811865475*nuVtSqSum[1]*edgeSurf[13]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[12]; 
  edgeSurf_incr[13] = 0.7071067811865475*nuVtSqSum[0]*edgeSurf[13]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[12]; 
  edgeSurf_incr[14] = 0.7071067811865475*nuVtSqSum[1]*edgeSurf[15]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[14]; 
  edgeSurf_incr[15] = 0.7071067811865475*nuVtSqSum[0]*edgeSurf[15]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[14]; 

  boundSurf_incr[0] = 0.7071067811865475*boundSurf[1]*nuVtSqSum[1]+0.7071067811865475*boundSurf[0]*nuVtSqSum[0]; 
  boundSurf_incr[1] = 0.7071067811865475*boundSurf[0]*nuVtSqSum[1]+0.7071067811865475*nuVtSqSum[0]*boundSurf[1]; 
  boundSurf_incr[2] = 0.7071067811865475*nuVtSqSum[1]*boundSurf[4]+0.7071067811865475*nuVtSqSum[0]*boundSurf[2]; 
  boundSurf_incr[3] = 0.7071067811865475*nuVtSqSum[1]*boundSurf[5]+0.7071067811865475*nuVtSqSum[0]*boundSurf[3]; 
  boundSurf_incr[4] = 0.7071067811865475*nuVtSqSum[0]*boundSurf[4]+0.7071067811865475*nuVtSqSum[1]*boundSurf[2]; 
  boundSurf_incr[5] = 0.7071067811865475*nuVtSqSum[0]*boundSurf[5]+0.7071067811865475*nuVtSqSum[1]*boundSurf[3]; 
  boundSurf_incr[6] = 0.7071067811865475*nuVtSqSum[1]*boundSurf[7]+0.7071067811865475*nuVtSqSum[0]*boundSurf[6]; 
  boundSurf_incr[7] = 0.7071067811865475*nuVtSqSum[0]*boundSurf[7]+0.7071067811865475*nuVtSqSum[1]*boundSurf[6]; 
  boundSurf_incr[8] = 0.7071067811865475*nuVtSqSum[1]*boundSurf[9]+0.7071067811865475*nuVtSqSum[0]*boundSurf[8]; 
  boundSurf_incr[9] = 0.7071067811865475*nuVtSqSum[0]*boundSurf[9]+0.7071067811865475*nuVtSqSum[1]*boundSurf[8]; 
  boundSurf_incr[10] = 0.7071067811865475*nuVtSqSum[1]*boundSurf[11]+0.7071067811865475*nuVtSqSum[0]*boundSurf[10]; 
  boundSurf_incr[11] = 0.7071067811865475*nuVtSqSum[0]*boundSurf[11]+0.7071067811865475*nuVtSqSum[1]*boundSurf[10]; 
  boundSurf_incr[12] = 0.7071067811865475*nuVtSqSum[1]*boundSurf[13]+0.7071067811865475*nuVtSqSum[0]*boundSurf[12]; 
  boundSurf_incr[13] = 0.7071067811865475*nuVtSqSum[0]*boundSurf[13]+0.7071067811865475*nuVtSqSum[1]*boundSurf[12]; 
  boundSurf_incr[14] = 0.7071067811865475*nuVtSqSum[1]*boundSurf[15]+0.7071067811865475*nuVtSqSum[0]*boundSurf[14]; 
  boundSurf_incr[15] = 0.7071067811865475*nuVtSqSum[0]*boundSurf[15]+0.7071067811865475*nuVtSqSum[1]*boundSurf[14]; 


  } else { 

  double edgeSurf[16] = {0.0}; 
  edgeSurf[0] = (-0.6708203932499369*fSkin[8])+0.6708203932499369*fEdge[8]+1.190784930203603*fSkin[2]+1.190784930203603*fEdge[2]-0.9375*fSkin[0]+0.9375*fEdge[0]; 
  edgeSurf[1] = (-0.6708203932499369*fSkin[9])+0.6708203932499369*fEdge[9]+1.190784930203603*fSkin[4]+1.190784930203603*fEdge[4]-0.9375*fSkin[1]+0.9375*fEdge[1]; 
  edgeSurf[2] = 1.585502557353661*fSkin[8]-0.7382874503707888*fEdge[8]-2.671875*fSkin[2]-1.453125*fEdge[2]+2.056810333988042*fSkin[0]-1.190784930203603*fEdge[0]; 
  edgeSurf[3] = (-0.6708203932499369*fSkin[10])+0.6708203932499369*fEdge[10]+1.190784930203603*fSkin[6]+1.190784930203603*fEdge[6]-0.9375*fSkin[3]+0.9375*fEdge[3]; 
  edgeSurf[4] = 1.585502557353661*fSkin[9]-0.7382874503707888*fEdge[9]-2.671875*fSkin[4]-1.453125*fEdge[4]+2.056810333988042*fSkin[1]-1.190784930203603*fEdge[1]; 
  edgeSurf[5] = (-0.6708203932499369*fSkin[11])+0.6708203932499369*fEdge[11]+1.190784930203603*fSkin[7]+1.190784930203603*fEdge[7]-0.9375*fSkin[5]+0.9375*fEdge[5]; 
  edgeSurf[6] = 1.585502557353661*fSkin[10]-0.7382874503707888*fEdge[10]-2.671875*fSkin[6]-1.453125*fEdge[6]+2.056810333988042*fSkin[3]-1.190784930203603*fEdge[3]; 
  edgeSurf[7] = 1.585502557353661*fSkin[11]-0.7382874503707888*fEdge[11]-2.671875*fSkin[7]-1.453125*fEdge[7]+2.056810333988042*fSkin[5]-1.190784930203603*fEdge[5]; 
  edgeSurf[8] = (-3.140625*fSkin[8])-0.140625*fEdge[8]+5.022775277112744*fSkin[2]+0.3025768239224545*fEdge[2]-3.773364712030896*fSkin[0]+0.4192627457812106*fEdge[0]; 
  edgeSurf[9] = (-3.140625*fSkin[9])-0.140625*fEdge[9]+5.022775277112744*fSkin[4]+0.3025768239224544*fEdge[4]-3.773364712030894*fSkin[1]+0.4192627457812105*fEdge[1]; 
  edgeSurf[10] = (-3.140625*fSkin[10])-0.140625*fEdge[10]+5.022775277112744*fSkin[6]+0.3025768239224544*fEdge[6]-3.773364712030894*fSkin[3]+0.4192627457812105*fEdge[3]; 
  edgeSurf[11] = (-3.140625*fSkin[11])-0.140625*fEdge[11]+5.022775277112744*fSkin[7]+0.3025768239224545*fEdge[7]-3.773364712030896*fSkin[5]+0.4192627457812106*fEdge[5]; 
  edgeSurf[12] = 1.190784930203603*fSkin[14]+1.190784930203603*fEdge[14]-0.9375*fSkin[12]+0.9375*fEdge[12]; 
  edgeSurf[13] = 1.190784930203603*fSkin[15]+1.190784930203603*fEdge[15]-0.9375*fSkin[13]+0.9375*fEdge[13]; 
  edgeSurf[14] = (-2.671875*fSkin[14])-1.453125*fEdge[14]+2.056810333988042*fSkin[12]-1.190784930203603*fEdge[12]; 
  edgeSurf[15] = (-2.671875*fSkin[15])-1.453125*fEdge[15]+2.056810333988042*fSkin[13]-1.190784930203603*fEdge[13]; 

  double boundSurf[16] = {0.0}; 
  boundSurf[2] = (-1.936491673103709*fSkin[8])-1.5*fSkin[2]-0.8660254037844386*fSkin[0]; 
  boundSurf[4] = (-1.936491673103709*fSkin[9])-1.5*fSkin[4]-0.8660254037844386*fSkin[1]; 
  boundSurf[6] = (-1.936491673103709*fSkin[10])-1.5*fSkin[6]-0.8660254037844386*fSkin[3]; 
  boundSurf[7] = (-1.936491673103709*fSkin[11])-1.5*fSkin[7]-0.8660254037844386*fSkin[5]; 
  boundSurf[8] = (-7.5*fSkin[8])-5.809475019311125*fSkin[2]-3.354101966249685*fSkin[0]; 
  boundSurf[9] = (-7.5*fSkin[9])-5.809475019311126*fSkin[4]-3.354101966249684*fSkin[1]; 
  boundSurf[10] = (-7.5*fSkin[10])-5.809475019311126*fSkin[6]-3.354101966249684*fSkin[3]; 
  boundSurf[11] = (-7.5*fSkin[11])-5.809475019311125*fSkin[7]-3.354101966249685*fSkin[5]; 
  boundSurf[14] = (-1.5*fSkin[14])-0.8660254037844387*fSkin[12]; 
  boundSurf[15] = (-1.5*fSkin[15])-0.8660254037844387*fSkin[13]; 

  edgeSurf_incr[0] = 0.7071067811865475*edgeSurf[1]*nuVtSqSum[1]+0.7071067811865475*edgeSurf[0]*nuVtSqSum[0]; 
  edgeSurf_incr[1] = 0.7071067811865475*edgeSurf[0]*nuVtSqSum[1]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[1]; 
  edgeSurf_incr[2] = 0.7071067811865475*nuVtSqSum[1]*edgeSurf[4]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[2]; 
  edgeSurf_incr[3] = 0.7071067811865475*nuVtSqSum[1]*edgeSurf[5]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[3]; 
  edgeSurf_incr[4] = 0.7071067811865475*nuVtSqSum[0]*edgeSurf[4]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[2]; 
  edgeSurf_incr[5] = 0.7071067811865475*nuVtSqSum[0]*edgeSurf[5]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[3]; 
  edgeSurf_incr[6] = 0.7071067811865475*nuVtSqSum[1]*edgeSurf[7]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[6]; 
  edgeSurf_incr[7] = 0.7071067811865475*nuVtSqSum[0]*edgeSurf[7]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[6]; 
  edgeSurf_incr[8] = 0.7071067811865475*nuVtSqSum[1]*edgeSurf[9]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[8]; 
  edgeSurf_incr[9] = 0.7071067811865475*nuVtSqSum[0]*edgeSurf[9]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[8]; 
  edgeSurf_incr[10] = 0.7071067811865475*nuVtSqSum[1]*edgeSurf[11]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[10]; 
  edgeSurf_incr[11] = 0.7071067811865475*nuVtSqSum[0]*edgeSurf[11]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[10]; 
  edgeSurf_incr[12] = 0.7071067811865475*nuVtSqSum[1]*edgeSurf[13]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[12]; 
  edgeSurf_incr[13] = 0.7071067811865475*nuVtSqSum[0]*edgeSurf[13]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[12]; 
  edgeSurf_incr[14] = 0.7071067811865475*nuVtSqSum[1]*edgeSurf[15]+0.7071067811865475*nuVtSqSum[0]*edgeSurf[14]; 
  edgeSurf_incr[15] = 0.7071067811865475*nuVtSqSum[0]*edgeSurf[15]+0.7071067811865475*nuVtSqSum[1]*edgeSurf[14]; 

  boundSurf_incr[0] = 0.7071067811865475*boundSurf[1]*nuVtSqSum[1]+0.7071067811865475*boundSurf[0]*nuVtSqSum[0]; 
  boundSurf_incr[1] = 0.7071067811865475*boundSurf[0]*nuVtSqSum[1]+0.7071067811865475*nuVtSqSum[0]*boundSurf[1]; 
  boundSurf_incr[2] = 0.7071067811865475*nuVtSqSum[1]*boundSurf[4]+0.7071067811865475*nuVtSqSum[0]*boundSurf[2]; 
  boundSurf_incr[3] = 0.7071067811865475*nuVtSqSum[1]*boundSurf[5]+0.7071067811865475*nuVtSqSum[0]*boundSurf[3]; 
  boundSurf_incr[4] = 0.7071067811865475*nuVtSqSum[0]*boundSurf[4]+0.7071067811865475*nuVtSqSum[1]*boundSurf[2]; 
  boundSurf_incr[5] = 0.7071067811865475*nuVtSqSum[0]*boundSurf[5]+0.7071067811865475*nuVtSqSum[1]*boundSurf[3]; 
  boundSurf_incr[6] = 0.7071067811865475*nuVtSqSum[1]*boundSurf[7]+0.7071067811865475*nuVtSqSum[0]*boundSurf[6]; 
  boundSurf_incr[7] = 0.7071067811865475*nuVtSqSum[0]*boundSurf[7]+0.7071067811865475*nuVtSqSum[1]*boundSurf[6]; 
  boundSurf_incr[8] = 0.7071067811865475*nuVtSqSum[1]*boundSurf[9]+0.7071067811865475*nuVtSqSum[0]*boundSurf[8]; 
  boundSurf_incr[9] = 0.7071067811865475*nuVtSqSum[0]*boundSurf[9]+0.7071067811865475*nuVtSqSum[1]*boundSurf[8]; 
  boundSurf_incr[10] = 0.7071067811865475*nuVtSqSum[1]*boundSurf[11]+0.7071067811865475*nuVtSqSum[0]*boundSurf[10]; 
  boundSurf_incr[11] = 0.7071067811865475*nuVtSqSum[0]*boundSurf[11]+0.7071067811865475*nuVtSqSum[1]*boundSurf[10]; 
  boundSurf_incr[12] = 0.7071067811865475*nuVtSqSum[1]*boundSurf[13]+0.7071067811865475*nuVtSqSum[0]*boundSurf[12]; 
  boundSurf_incr[13] = 0.7071067811865475*nuVtSqSum[0]*boundSurf[13]+0.7071067811865475*nuVtSqSum[1]*boundSurf[12]; 
  boundSurf_incr[14] = 0.7071067811865475*nuVtSqSum[1]*boundSurf[15]+0.7071067811865475*nuVtSqSum[0]*boundSurf[14]; 
  boundSurf_incr[15] = 0.7071067811865475*nuVtSqSum[0]*boundSurf[15]+0.7071067811865475*nuVtSqSum[1]*boundSurf[14]; 

  } 

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*rdvSq4; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*rdvSq4; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*rdvSq4; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*rdvSq4; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*rdvSq4; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*rdvSq4; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*rdvSq4; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*rdvSq4; 
  out[8] += (vol_incr[8]+edgeSurf_incr[8]+boundSurf_incr[8])*rdvSq4; 
  out[9] += (vol_incr[9]+edgeSurf_incr[9]+boundSurf_incr[9])*rdvSq4; 
  out[10] += (vol_incr[10]+edgeSurf_incr[10]+boundSurf_incr[10])*rdvSq4; 
  out[11] += (vol_incr[11]+edgeSurf_incr[11]+boundSurf_incr[11])*rdvSq4; 
  out[12] += (vol_incr[12]+edgeSurf_incr[12]+boundSurf_incr[12])*rdvSq4; 
  out[13] += (vol_incr[13]+edgeSurf_incr[13]+boundSurf_incr[13])*rdvSq4; 
  out[14] += (vol_incr[14]+edgeSurf_incr[14]+boundSurf_incr[14])*rdvSq4; 
  out[15] += (vol_incr[15]+edgeSurf_incr[15]+boundSurf_incr[15])*rdvSq4; 

  return 0.;

} 
