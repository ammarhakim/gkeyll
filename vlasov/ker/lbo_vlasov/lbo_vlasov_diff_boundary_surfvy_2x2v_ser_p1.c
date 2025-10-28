#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_diff_boundary_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[4]: Cell-center coordinates. 
  // dxv[4]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[12]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fSkin/Edge: Distribution function in cells 
  // out: Incremented distribution function in cell 
  const double *nuVtSqSum = &nuPrimMomsSum[8];

  double rdvSq4 = 4.0/(dxv[3]*dxv[3]); 

  double vol_incr[32] = {0.0}; 
  vol_incr[24] = 3.354101966249685*nuVtSqSum[3]*fSkin[5]+3.354101966249685*fSkin[2]*nuVtSqSum[2]+3.354101966249685*fSkin[1]*nuVtSqSum[1]+3.354101966249685*fSkin[0]*nuVtSqSum[0]; 
  vol_incr[25] = 3.354101966249684*nuVtSqSum[2]*fSkin[5]+3.354101966249684*fSkin[2]*nuVtSqSum[3]+3.354101966249684*fSkin[0]*nuVtSqSum[1]+3.354101966249684*nuVtSqSum[0]*fSkin[1]; 
  vol_incr[26] = 3.354101966249684*nuVtSqSum[1]*fSkin[5]+3.354101966249684*fSkin[1]*nuVtSqSum[3]+3.354101966249684*fSkin[0]*nuVtSqSum[2]+3.354101966249684*nuVtSqSum[0]*fSkin[2]; 
  vol_incr[27] = 3.354101966249684*nuVtSqSum[3]*fSkin[11]+3.354101966249684*nuVtSqSum[2]*fSkin[7]+3.354101966249684*nuVtSqSum[1]*fSkin[6]+3.354101966249684*nuVtSqSum[0]*fSkin[3]; 
  vol_incr[28] = 3.354101966249685*nuVtSqSum[0]*fSkin[5]+3.354101966249685*fSkin[0]*nuVtSqSum[3]+3.354101966249685*fSkin[1]*nuVtSqSum[2]+3.354101966249685*nuVtSqSum[1]*fSkin[2]; 
  vol_incr[29] = 3.354101966249685*nuVtSqSum[2]*fSkin[11]+3.354101966249685*nuVtSqSum[3]*fSkin[7]+3.354101966249685*nuVtSqSum[0]*fSkin[6]+3.354101966249685*nuVtSqSum[1]*fSkin[3]; 
  vol_incr[30] = 3.354101966249685*nuVtSqSum[1]*fSkin[11]+3.354101966249685*nuVtSqSum[0]*fSkin[7]+3.354101966249685*nuVtSqSum[3]*fSkin[6]+3.354101966249685*nuVtSqSum[2]*fSkin[3]; 
  vol_incr[31] = 3.354101966249684*nuVtSqSum[0]*fSkin[11]+3.354101966249684*nuVtSqSum[1]*fSkin[7]+3.354101966249684*nuVtSqSum[2]*fSkin[6]+3.354101966249684*fSkin[3]*nuVtSqSum[3]; 

  double edgeSurf_incr[32] = {0.0}; 
  double boundSurf_incr[32] = {0.0}; 

  if (edge == -1) { 

  double edgeSurf[32] = {0.0}; 
  edgeSurf[0] = (-0.6708203932499369*fSkin[24])+0.6708203932499369*fEdge[24]-1.190784930203603*fSkin[4]-1.190784930203603*fEdge[4]-0.9375*fSkin[0]+0.9375*fEdge[0]; 
  edgeSurf[1] = (-0.6708203932499369*fSkin[25])+0.6708203932499369*fEdge[25]-1.190784930203603*fSkin[8]-1.190784930203603*fEdge[8]-0.9375*fSkin[1]+0.9375*fEdge[1]; 
  edgeSurf[2] = (-0.6708203932499369*fSkin[26])+0.6708203932499369*fEdge[26]-1.190784930203603*fSkin[9]-1.190784930203603*fEdge[9]-0.9375*fSkin[2]+0.9375*fEdge[2]; 
  edgeSurf[3] = (-0.6708203932499369*fSkin[27])+0.6708203932499369*fEdge[27]-1.190784930203603*fSkin[10]-1.190784930203603*fEdge[10]-0.9375*fSkin[3]+0.9375*fEdge[3]; 
  edgeSurf[4] = (-1.585502557353661*fSkin[24])+0.7382874503707888*fEdge[24]-2.671875*fSkin[4]-1.453125*fEdge[4]-2.056810333988042*fSkin[0]+1.190784930203603*fEdge[0]; 
  edgeSurf[5] = (-0.6708203932499369*fSkin[28])+0.6708203932499369*fEdge[28]-1.190784930203603*fSkin[12]-1.190784930203603*fEdge[12]-0.9375*fSkin[5]+0.9375*fEdge[5]; 
  edgeSurf[6] = (-0.6708203932499369*fSkin[29])+0.6708203932499369*fEdge[29]-1.190784930203603*fSkin[13]-1.190784930203603*fEdge[13]-0.9375*fSkin[6]+0.9375*fEdge[6]; 
  edgeSurf[7] = (-0.6708203932499369*fSkin[30])+0.6708203932499369*fEdge[30]-1.190784930203603*fSkin[14]-1.190784930203603*fEdge[14]-0.9375*fSkin[7]+0.9375*fEdge[7]; 
  edgeSurf[8] = (-1.585502557353661*fSkin[25])+0.7382874503707888*fEdge[25]-2.671875*fSkin[8]-1.453125*fEdge[8]-2.056810333988042*fSkin[1]+1.190784930203603*fEdge[1]; 
  edgeSurf[9] = (-1.585502557353661*fSkin[26])+0.7382874503707888*fEdge[26]-2.671875*fSkin[9]-1.453125*fEdge[9]-2.056810333988042*fSkin[2]+1.190784930203603*fEdge[2]; 
  edgeSurf[10] = (-1.585502557353661*fSkin[27])+0.7382874503707888*fEdge[27]-2.671875*fSkin[10]-1.453125*fEdge[10]-2.056810333988042*fSkin[3]+1.190784930203603*fEdge[3]; 
  edgeSurf[11] = (-0.6708203932499369*fSkin[31])+0.6708203932499369*fEdge[31]-1.190784930203603*fSkin[15]-1.190784930203603*fEdge[15]-0.9375*fSkin[11]+0.9375*fEdge[11]; 
  edgeSurf[12] = (-1.585502557353661*fSkin[28])+0.7382874503707888*fEdge[28]-2.671875*fSkin[12]-1.453125*fEdge[12]-2.056810333988042*fSkin[5]+1.190784930203603*fEdge[5]; 
  edgeSurf[13] = (-1.585502557353661*fSkin[29])+0.7382874503707888*fEdge[29]-2.671875*fSkin[13]-1.453125*fEdge[13]-2.056810333988042*fSkin[6]+1.190784930203603*fEdge[6]; 
  edgeSurf[14] = (-1.585502557353661*fSkin[30])+0.7382874503707888*fEdge[30]-2.671875*fSkin[14]-1.453125*fEdge[14]-2.056810333988042*fSkin[7]+1.190784930203603*fEdge[7]; 
  edgeSurf[15] = (-1.585502557353661*fSkin[31])+0.7382874503707888*fEdge[31]-2.671875*fSkin[15]-1.453125*fEdge[15]-2.056810333988042*fSkin[11]+1.190784930203603*fEdge[11]; 
  edgeSurf[16] = (-1.190784930203603*fSkin[19])-1.190784930203603*fEdge[19]-0.9375*fSkin[16]+0.9375*fEdge[16]; 
  edgeSurf[17] = (-1.190784930203603*fSkin[21])-1.190784930203603*fEdge[21]-0.9375*fSkin[17]+0.9375*fEdge[17]; 
  edgeSurf[18] = (-1.190784930203603*fSkin[22])-1.190784930203603*fEdge[22]-0.9375*fSkin[18]+0.9375*fEdge[18]; 
  edgeSurf[19] = (-2.671875*fSkin[19])-1.453125*fEdge[19]-2.056810333988042*fSkin[16]+1.190784930203603*fEdge[16]; 
  edgeSurf[20] = (-1.190784930203603*fSkin[23])-1.190784930203603*fEdge[23]-0.9375*fSkin[20]+0.9375*fEdge[20]; 
  edgeSurf[21] = (-2.671875*fSkin[21])-1.453125*fEdge[21]-2.056810333988042*fSkin[17]+1.190784930203603*fEdge[17]; 
  edgeSurf[22] = (-2.671875*fSkin[22])-1.453125*fEdge[22]-2.056810333988042*fSkin[18]+1.190784930203603*fEdge[18]; 
  edgeSurf[23] = (-2.671875*fSkin[23])-1.453125*fEdge[23]-2.056810333988042*fSkin[20]+1.190784930203603*fEdge[20]; 
  edgeSurf[24] = (-3.140625*fSkin[24])-0.140625*fEdge[24]-5.022775277112744*fSkin[4]-0.3025768239224545*fEdge[4]-3.773364712030896*fSkin[0]+0.4192627457812106*fEdge[0]; 
  edgeSurf[25] = (-3.140625*fSkin[25])-0.140625*fEdge[25]-5.022775277112744*fSkin[8]-0.3025768239224544*fEdge[8]-3.773364712030894*fSkin[1]+0.4192627457812105*fEdge[1]; 
  edgeSurf[26] = (-3.140625*fSkin[26])-0.140625*fEdge[26]-5.022775277112744*fSkin[9]-0.3025768239224544*fEdge[9]-3.773364712030894*fSkin[2]+0.4192627457812105*fEdge[2]; 
  edgeSurf[27] = (-3.140625*fSkin[27])-0.140625*fEdge[27]-5.022775277112744*fSkin[10]-0.3025768239224544*fEdge[10]-3.773364712030894*fSkin[3]+0.4192627457812105*fEdge[3]; 
  edgeSurf[28] = (-3.140625*fSkin[28])-0.140625*fEdge[28]-5.022775277112744*fSkin[12]-0.3025768239224545*fEdge[12]-3.773364712030896*fSkin[5]+0.4192627457812106*fEdge[5]; 
  edgeSurf[29] = (-3.140625*fSkin[29])-0.140625*fEdge[29]-5.022775277112744*fSkin[13]-0.3025768239224545*fEdge[13]-3.773364712030896*fSkin[6]+0.4192627457812106*fEdge[6]; 
  edgeSurf[30] = (-3.140625*fSkin[30])-0.140625*fEdge[30]-5.022775277112744*fSkin[14]-0.3025768239224545*fEdge[14]-3.773364712030896*fSkin[7]+0.4192627457812106*fEdge[7]; 
  edgeSurf[31] = (-3.140625*fSkin[31])-0.140625*fEdge[31]-5.022775277112744*fSkin[15]-0.3025768239224544*fEdge[15]-3.773364712030894*fSkin[11]+0.4192627457812105*fEdge[11]; 

  double boundSurf[32] = {0.0}; 
  boundSurf[4] = 1.936491673103709*fSkin[24]-1.5*fSkin[4]+0.8660254037844386*fSkin[0]; 
  boundSurf[8] = 1.936491673103709*fSkin[25]-1.5*fSkin[8]+0.8660254037844386*fSkin[1]; 
  boundSurf[9] = 1.936491673103709*fSkin[26]-1.5*fSkin[9]+0.8660254037844386*fSkin[2]; 
  boundSurf[10] = 1.936491673103709*fSkin[27]-1.5*fSkin[10]+0.8660254037844386*fSkin[3]; 
  boundSurf[12] = 1.936491673103709*fSkin[28]-1.5*fSkin[12]+0.8660254037844386*fSkin[5]; 
  boundSurf[13] = 1.936491673103709*fSkin[29]-1.5*fSkin[13]+0.8660254037844386*fSkin[6]; 
  boundSurf[14] = 1.936491673103709*fSkin[30]-1.5*fSkin[14]+0.8660254037844386*fSkin[7]; 
  boundSurf[15] = 1.936491673103709*fSkin[31]-1.5*fSkin[15]+0.8660254037844386*fSkin[11]; 
  boundSurf[19] = 0.8660254037844387*fSkin[16]-1.5*fSkin[19]; 
  boundSurf[21] = 0.8660254037844387*fSkin[17]-1.5*fSkin[21]; 
  boundSurf[22] = 0.8660254037844387*fSkin[18]-1.5*fSkin[22]; 
  boundSurf[23] = 0.8660254037844387*fSkin[20]-1.5*fSkin[23]; 
  boundSurf[24] = (-7.5*fSkin[24])+5.809475019311125*fSkin[4]-3.354101966249685*fSkin[0]; 
  boundSurf[25] = (-7.5*fSkin[25])+5.809475019311126*fSkin[8]-3.354101966249684*fSkin[1]; 
  boundSurf[26] = (-7.5*fSkin[26])+5.809475019311126*fSkin[9]-3.354101966249684*fSkin[2]; 
  boundSurf[27] = (-7.5*fSkin[27])+5.809475019311126*fSkin[10]-3.354101966249684*fSkin[3]; 
  boundSurf[28] = (-7.5*fSkin[28])+5.809475019311125*fSkin[12]-3.354101966249685*fSkin[5]; 
  boundSurf[29] = (-7.5*fSkin[29])+5.809475019311125*fSkin[13]-3.354101966249685*fSkin[6]; 
  boundSurf[30] = (-7.5*fSkin[30])+5.809475019311125*fSkin[14]-3.354101966249685*fSkin[7]; 
  boundSurf[31] = (-7.5*fSkin[31])+5.809475019311126*fSkin[15]-3.354101966249684*fSkin[11]; 

  edgeSurf_incr[0] = 0.5*nuVtSqSum[3]*edgeSurf[5]+0.5*edgeSurf[2]*nuVtSqSum[2]+0.5*edgeSurf[1]*nuVtSqSum[1]+0.5*edgeSurf[0]*nuVtSqSum[0]; 
  edgeSurf_incr[1] = 0.5*nuVtSqSum[2]*edgeSurf[5]+0.5*edgeSurf[2]*nuVtSqSum[3]+0.5*edgeSurf[0]*nuVtSqSum[1]+0.5*nuVtSqSum[0]*edgeSurf[1]; 
  edgeSurf_incr[2] = 0.5*nuVtSqSum[1]*edgeSurf[5]+0.5*edgeSurf[1]*nuVtSqSum[3]+0.5*edgeSurf[0]*nuVtSqSum[2]+0.5*nuVtSqSum[0]*edgeSurf[2]; 
  edgeSurf_incr[3] = 0.5*nuVtSqSum[3]*edgeSurf[11]+0.5*nuVtSqSum[2]*edgeSurf[7]+0.5*nuVtSqSum[1]*edgeSurf[6]+0.5*nuVtSqSum[0]*edgeSurf[3]; 
  edgeSurf_incr[4] = 0.5*nuVtSqSum[3]*edgeSurf[12]+0.5*nuVtSqSum[2]*edgeSurf[9]+0.5*nuVtSqSum[1]*edgeSurf[8]+0.5*nuVtSqSum[0]*edgeSurf[4]; 
  edgeSurf_incr[5] = 0.5*nuVtSqSum[0]*edgeSurf[5]+0.5*edgeSurf[0]*nuVtSqSum[3]+0.5*edgeSurf[1]*nuVtSqSum[2]+0.5*nuVtSqSum[1]*edgeSurf[2]; 
  edgeSurf_incr[6] = 0.5*nuVtSqSum[2]*edgeSurf[11]+0.5*nuVtSqSum[3]*edgeSurf[7]+0.5*nuVtSqSum[0]*edgeSurf[6]+0.5*nuVtSqSum[1]*edgeSurf[3]; 
  edgeSurf_incr[7] = 0.5*nuVtSqSum[1]*edgeSurf[11]+0.5*nuVtSqSum[0]*edgeSurf[7]+0.5*nuVtSqSum[3]*edgeSurf[6]+0.5*nuVtSqSum[2]*edgeSurf[3]; 
  edgeSurf_incr[8] = 0.5*nuVtSqSum[2]*edgeSurf[12]+0.5*nuVtSqSum[3]*edgeSurf[9]+0.5*nuVtSqSum[0]*edgeSurf[8]+0.5*nuVtSqSum[1]*edgeSurf[4]; 
  edgeSurf_incr[9] = 0.5*nuVtSqSum[1]*edgeSurf[12]+0.5*nuVtSqSum[0]*edgeSurf[9]+0.5*nuVtSqSum[3]*edgeSurf[8]+0.5*nuVtSqSum[2]*edgeSurf[4]; 
  edgeSurf_incr[10] = 0.5*nuVtSqSum[3]*edgeSurf[15]+0.5*nuVtSqSum[2]*edgeSurf[14]+0.5*nuVtSqSum[1]*edgeSurf[13]+0.5*nuVtSqSum[0]*edgeSurf[10]; 
  edgeSurf_incr[11] = 0.5*nuVtSqSum[0]*edgeSurf[11]+0.5*nuVtSqSum[1]*edgeSurf[7]+0.5*nuVtSqSum[2]*edgeSurf[6]+0.5*edgeSurf[3]*nuVtSqSum[3]; 
  edgeSurf_incr[12] = 0.5*nuVtSqSum[0]*edgeSurf[12]+0.5*nuVtSqSum[1]*edgeSurf[9]+0.5*nuVtSqSum[2]*edgeSurf[8]+0.5*nuVtSqSum[3]*edgeSurf[4]; 
  edgeSurf_incr[13] = 0.5*nuVtSqSum[2]*edgeSurf[15]+0.5*nuVtSqSum[3]*edgeSurf[14]+0.5*nuVtSqSum[0]*edgeSurf[13]+0.5*nuVtSqSum[1]*edgeSurf[10]; 
  edgeSurf_incr[14] = 0.5*nuVtSqSum[1]*edgeSurf[15]+0.5*nuVtSqSum[0]*edgeSurf[14]+0.5*nuVtSqSum[3]*edgeSurf[13]+0.5*nuVtSqSum[2]*edgeSurf[10]; 
  edgeSurf_incr[15] = 0.5*nuVtSqSum[0]*edgeSurf[15]+0.5*nuVtSqSum[1]*edgeSurf[14]+0.5*nuVtSqSum[2]*edgeSurf[13]+0.5*nuVtSqSum[3]*edgeSurf[10]; 
  edgeSurf_incr[16] = 0.5*nuVtSqSum[3]*edgeSurf[20]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[18]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[17]+0.5*nuVtSqSum[0]*edgeSurf[16]; 
  edgeSurf_incr[17] = 0.5000000000000001*nuVtSqSum[2]*edgeSurf[20]+0.5*nuVtSqSum[3]*edgeSurf[18]+0.5*nuVtSqSum[0]*edgeSurf[17]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[16]; 
  edgeSurf_incr[18] = 0.5000000000000001*nuVtSqSum[1]*edgeSurf[20]+0.5*nuVtSqSum[0]*edgeSurf[18]+0.5*nuVtSqSum[3]*edgeSurf[17]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[16]; 
  edgeSurf_incr[19] = 0.5*nuVtSqSum[3]*edgeSurf[23]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[22]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[21]+0.5*nuVtSqSum[0]*edgeSurf[19]; 
  edgeSurf_incr[20] = 0.5*nuVtSqSum[0]*edgeSurf[20]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[18]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[17]+0.5*nuVtSqSum[3]*edgeSurf[16]; 
  edgeSurf_incr[21] = 0.5000000000000001*nuVtSqSum[2]*edgeSurf[23]+0.5*nuVtSqSum[3]*edgeSurf[22]+0.5*nuVtSqSum[0]*edgeSurf[21]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[19]; 
  edgeSurf_incr[22] = 0.5000000000000001*nuVtSqSum[1]*edgeSurf[23]+0.5*nuVtSqSum[0]*edgeSurf[22]+0.5*nuVtSqSum[3]*edgeSurf[21]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[19]; 
  edgeSurf_incr[23] = 0.5*nuVtSqSum[0]*edgeSurf[23]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[22]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[21]+0.5*nuVtSqSum[3]*edgeSurf[19]; 
  edgeSurf_incr[24] = 0.5*nuVtSqSum[3]*edgeSurf[28]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[26]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[25]+0.5*nuVtSqSum[0]*edgeSurf[24]; 
  edgeSurf_incr[25] = 0.5000000000000001*nuVtSqSum[2]*edgeSurf[28]+0.5*nuVtSqSum[3]*edgeSurf[26]+0.5*nuVtSqSum[0]*edgeSurf[25]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[24]; 
  edgeSurf_incr[26] = 0.5000000000000001*nuVtSqSum[1]*edgeSurf[28]+0.5*nuVtSqSum[0]*edgeSurf[26]+0.5*nuVtSqSum[3]*edgeSurf[25]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[24]; 
  edgeSurf_incr[27] = 0.5*nuVtSqSum[3]*edgeSurf[31]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[30]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[29]+0.5*nuVtSqSum[0]*edgeSurf[27]; 
  edgeSurf_incr[28] = 0.5*nuVtSqSum[0]*edgeSurf[28]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[26]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[25]+0.5*nuVtSqSum[3]*edgeSurf[24]; 
  edgeSurf_incr[29] = 0.5000000000000001*nuVtSqSum[2]*edgeSurf[31]+0.5*nuVtSqSum[3]*edgeSurf[30]+0.5*nuVtSqSum[0]*edgeSurf[29]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[27]; 
  edgeSurf_incr[30] = 0.5000000000000001*nuVtSqSum[1]*edgeSurf[31]+0.5*nuVtSqSum[0]*edgeSurf[30]+0.5*nuVtSqSum[3]*edgeSurf[29]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[27]; 
  edgeSurf_incr[31] = 0.5*nuVtSqSum[0]*edgeSurf[31]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[30]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[29]+0.5*nuVtSqSum[3]*edgeSurf[27]; 

  boundSurf_incr[0] = 0.5*nuVtSqSum[3]*boundSurf[5]+0.5*boundSurf[2]*nuVtSqSum[2]+0.5*boundSurf[1]*nuVtSqSum[1]+0.5*boundSurf[0]*nuVtSqSum[0]; 
  boundSurf_incr[1] = 0.5*nuVtSqSum[2]*boundSurf[5]+0.5*boundSurf[2]*nuVtSqSum[3]+0.5*boundSurf[0]*nuVtSqSum[1]+0.5*nuVtSqSum[0]*boundSurf[1]; 
  boundSurf_incr[2] = 0.5*nuVtSqSum[1]*boundSurf[5]+0.5*boundSurf[1]*nuVtSqSum[3]+0.5*boundSurf[0]*nuVtSqSum[2]+0.5*nuVtSqSum[0]*boundSurf[2]; 
  boundSurf_incr[3] = 0.5*nuVtSqSum[3]*boundSurf[11]+0.5*nuVtSqSum[2]*boundSurf[7]+0.5*nuVtSqSum[1]*boundSurf[6]+0.5*nuVtSqSum[0]*boundSurf[3]; 
  boundSurf_incr[4] = 0.5*nuVtSqSum[3]*boundSurf[12]+0.5*nuVtSqSum[2]*boundSurf[9]+0.5*nuVtSqSum[1]*boundSurf[8]+0.5*nuVtSqSum[0]*boundSurf[4]; 
  boundSurf_incr[5] = 0.5*nuVtSqSum[0]*boundSurf[5]+0.5*boundSurf[0]*nuVtSqSum[3]+0.5*boundSurf[1]*nuVtSqSum[2]+0.5*nuVtSqSum[1]*boundSurf[2]; 
  boundSurf_incr[6] = 0.5*nuVtSqSum[2]*boundSurf[11]+0.5*nuVtSqSum[3]*boundSurf[7]+0.5*nuVtSqSum[0]*boundSurf[6]+0.5*nuVtSqSum[1]*boundSurf[3]; 
  boundSurf_incr[7] = 0.5*nuVtSqSum[1]*boundSurf[11]+0.5*nuVtSqSum[0]*boundSurf[7]+0.5*nuVtSqSum[3]*boundSurf[6]+0.5*nuVtSqSum[2]*boundSurf[3]; 
  boundSurf_incr[8] = 0.5*nuVtSqSum[2]*boundSurf[12]+0.5*nuVtSqSum[3]*boundSurf[9]+0.5*nuVtSqSum[0]*boundSurf[8]+0.5*nuVtSqSum[1]*boundSurf[4]; 
  boundSurf_incr[9] = 0.5*nuVtSqSum[1]*boundSurf[12]+0.5*nuVtSqSum[0]*boundSurf[9]+0.5*nuVtSqSum[3]*boundSurf[8]+0.5*nuVtSqSum[2]*boundSurf[4]; 
  boundSurf_incr[10] = 0.5*nuVtSqSum[3]*boundSurf[15]+0.5*nuVtSqSum[2]*boundSurf[14]+0.5*nuVtSqSum[1]*boundSurf[13]+0.5*nuVtSqSum[0]*boundSurf[10]; 
  boundSurf_incr[11] = 0.5*nuVtSqSum[0]*boundSurf[11]+0.5*nuVtSqSum[1]*boundSurf[7]+0.5*nuVtSqSum[2]*boundSurf[6]+0.5*boundSurf[3]*nuVtSqSum[3]; 
  boundSurf_incr[12] = 0.5*nuVtSqSum[0]*boundSurf[12]+0.5*nuVtSqSum[1]*boundSurf[9]+0.5*nuVtSqSum[2]*boundSurf[8]+0.5*nuVtSqSum[3]*boundSurf[4]; 
  boundSurf_incr[13] = 0.5*nuVtSqSum[2]*boundSurf[15]+0.5*nuVtSqSum[3]*boundSurf[14]+0.5*nuVtSqSum[0]*boundSurf[13]+0.5*nuVtSqSum[1]*boundSurf[10]; 
  boundSurf_incr[14] = 0.5*nuVtSqSum[1]*boundSurf[15]+0.5*nuVtSqSum[0]*boundSurf[14]+0.5*nuVtSqSum[3]*boundSurf[13]+0.5*nuVtSqSum[2]*boundSurf[10]; 
  boundSurf_incr[15] = 0.5*nuVtSqSum[0]*boundSurf[15]+0.5*nuVtSqSum[1]*boundSurf[14]+0.5*nuVtSqSum[2]*boundSurf[13]+0.5*nuVtSqSum[3]*boundSurf[10]; 
  boundSurf_incr[16] = 0.5*nuVtSqSum[3]*boundSurf[20]+0.5000000000000001*nuVtSqSum[2]*boundSurf[18]+0.5000000000000001*nuVtSqSum[1]*boundSurf[17]+0.5*nuVtSqSum[0]*boundSurf[16]; 
  boundSurf_incr[17] = 0.5000000000000001*nuVtSqSum[2]*boundSurf[20]+0.5*nuVtSqSum[3]*boundSurf[18]+0.5*nuVtSqSum[0]*boundSurf[17]+0.5000000000000001*nuVtSqSum[1]*boundSurf[16]; 
  boundSurf_incr[18] = 0.5000000000000001*nuVtSqSum[1]*boundSurf[20]+0.5*nuVtSqSum[0]*boundSurf[18]+0.5*nuVtSqSum[3]*boundSurf[17]+0.5000000000000001*nuVtSqSum[2]*boundSurf[16]; 
  boundSurf_incr[19] = 0.5*nuVtSqSum[3]*boundSurf[23]+0.5000000000000001*nuVtSqSum[2]*boundSurf[22]+0.5000000000000001*nuVtSqSum[1]*boundSurf[21]+0.5*nuVtSqSum[0]*boundSurf[19]; 
  boundSurf_incr[20] = 0.5*nuVtSqSum[0]*boundSurf[20]+0.5000000000000001*nuVtSqSum[1]*boundSurf[18]+0.5000000000000001*nuVtSqSum[2]*boundSurf[17]+0.5*nuVtSqSum[3]*boundSurf[16]; 
  boundSurf_incr[21] = 0.5000000000000001*nuVtSqSum[2]*boundSurf[23]+0.5*nuVtSqSum[3]*boundSurf[22]+0.5*nuVtSqSum[0]*boundSurf[21]+0.5000000000000001*nuVtSqSum[1]*boundSurf[19]; 
  boundSurf_incr[22] = 0.5000000000000001*nuVtSqSum[1]*boundSurf[23]+0.5*nuVtSqSum[0]*boundSurf[22]+0.5*nuVtSqSum[3]*boundSurf[21]+0.5000000000000001*nuVtSqSum[2]*boundSurf[19]; 
  boundSurf_incr[23] = 0.5*nuVtSqSum[0]*boundSurf[23]+0.5000000000000001*nuVtSqSum[1]*boundSurf[22]+0.5000000000000001*nuVtSqSum[2]*boundSurf[21]+0.5*nuVtSqSum[3]*boundSurf[19]; 
  boundSurf_incr[24] = 0.5*nuVtSqSum[3]*boundSurf[28]+0.5000000000000001*nuVtSqSum[2]*boundSurf[26]+0.5000000000000001*nuVtSqSum[1]*boundSurf[25]+0.5*nuVtSqSum[0]*boundSurf[24]; 
  boundSurf_incr[25] = 0.5000000000000001*nuVtSqSum[2]*boundSurf[28]+0.5*nuVtSqSum[3]*boundSurf[26]+0.5*nuVtSqSum[0]*boundSurf[25]+0.5000000000000001*nuVtSqSum[1]*boundSurf[24]; 
  boundSurf_incr[26] = 0.5000000000000001*nuVtSqSum[1]*boundSurf[28]+0.5*nuVtSqSum[0]*boundSurf[26]+0.5*nuVtSqSum[3]*boundSurf[25]+0.5000000000000001*nuVtSqSum[2]*boundSurf[24]; 
  boundSurf_incr[27] = 0.5*nuVtSqSum[3]*boundSurf[31]+0.5000000000000001*nuVtSqSum[2]*boundSurf[30]+0.5000000000000001*nuVtSqSum[1]*boundSurf[29]+0.5*nuVtSqSum[0]*boundSurf[27]; 
  boundSurf_incr[28] = 0.5*nuVtSqSum[0]*boundSurf[28]+0.5000000000000001*nuVtSqSum[1]*boundSurf[26]+0.5000000000000001*nuVtSqSum[2]*boundSurf[25]+0.5*nuVtSqSum[3]*boundSurf[24]; 
  boundSurf_incr[29] = 0.5000000000000001*nuVtSqSum[2]*boundSurf[31]+0.5*nuVtSqSum[3]*boundSurf[30]+0.5*nuVtSqSum[0]*boundSurf[29]+0.5000000000000001*nuVtSqSum[1]*boundSurf[27]; 
  boundSurf_incr[30] = 0.5000000000000001*nuVtSqSum[1]*boundSurf[31]+0.5*nuVtSqSum[0]*boundSurf[30]+0.5*nuVtSqSum[3]*boundSurf[29]+0.5000000000000001*nuVtSqSum[2]*boundSurf[27]; 
  boundSurf_incr[31] = 0.5*nuVtSqSum[0]*boundSurf[31]+0.5000000000000001*nuVtSqSum[1]*boundSurf[30]+0.5000000000000001*nuVtSqSum[2]*boundSurf[29]+0.5*nuVtSqSum[3]*boundSurf[27]; 


  } else { 

  double edgeSurf[32] = {0.0}; 
  edgeSurf[0] = (-0.6708203932499369*fSkin[24])+0.6708203932499369*fEdge[24]+1.190784930203603*fSkin[4]+1.190784930203603*fEdge[4]-0.9375*fSkin[0]+0.9375*fEdge[0]; 
  edgeSurf[1] = (-0.6708203932499369*fSkin[25])+0.6708203932499369*fEdge[25]+1.190784930203603*fSkin[8]+1.190784930203603*fEdge[8]-0.9375*fSkin[1]+0.9375*fEdge[1]; 
  edgeSurf[2] = (-0.6708203932499369*fSkin[26])+0.6708203932499369*fEdge[26]+1.190784930203603*fSkin[9]+1.190784930203603*fEdge[9]-0.9375*fSkin[2]+0.9375*fEdge[2]; 
  edgeSurf[3] = (-0.6708203932499369*fSkin[27])+0.6708203932499369*fEdge[27]+1.190784930203603*fSkin[10]+1.190784930203603*fEdge[10]-0.9375*fSkin[3]+0.9375*fEdge[3]; 
  edgeSurf[4] = 1.585502557353661*fSkin[24]-0.7382874503707888*fEdge[24]-2.671875*fSkin[4]-1.453125*fEdge[4]+2.056810333988042*fSkin[0]-1.190784930203603*fEdge[0]; 
  edgeSurf[5] = (-0.6708203932499369*fSkin[28])+0.6708203932499369*fEdge[28]+1.190784930203603*fSkin[12]+1.190784930203603*fEdge[12]-0.9375*fSkin[5]+0.9375*fEdge[5]; 
  edgeSurf[6] = (-0.6708203932499369*fSkin[29])+0.6708203932499369*fEdge[29]+1.190784930203603*fSkin[13]+1.190784930203603*fEdge[13]-0.9375*fSkin[6]+0.9375*fEdge[6]; 
  edgeSurf[7] = (-0.6708203932499369*fSkin[30])+0.6708203932499369*fEdge[30]+1.190784930203603*fSkin[14]+1.190784930203603*fEdge[14]-0.9375*fSkin[7]+0.9375*fEdge[7]; 
  edgeSurf[8] = 1.585502557353661*fSkin[25]-0.7382874503707888*fEdge[25]-2.671875*fSkin[8]-1.453125*fEdge[8]+2.056810333988042*fSkin[1]-1.190784930203603*fEdge[1]; 
  edgeSurf[9] = 1.585502557353661*fSkin[26]-0.7382874503707888*fEdge[26]-2.671875*fSkin[9]-1.453125*fEdge[9]+2.056810333988042*fSkin[2]-1.190784930203603*fEdge[2]; 
  edgeSurf[10] = 1.585502557353661*fSkin[27]-0.7382874503707888*fEdge[27]-2.671875*fSkin[10]-1.453125*fEdge[10]+2.056810333988042*fSkin[3]-1.190784930203603*fEdge[3]; 
  edgeSurf[11] = (-0.6708203932499369*fSkin[31])+0.6708203932499369*fEdge[31]+1.190784930203603*fSkin[15]+1.190784930203603*fEdge[15]-0.9375*fSkin[11]+0.9375*fEdge[11]; 
  edgeSurf[12] = 1.585502557353661*fSkin[28]-0.7382874503707888*fEdge[28]-2.671875*fSkin[12]-1.453125*fEdge[12]+2.056810333988042*fSkin[5]-1.190784930203603*fEdge[5]; 
  edgeSurf[13] = 1.585502557353661*fSkin[29]-0.7382874503707888*fEdge[29]-2.671875*fSkin[13]-1.453125*fEdge[13]+2.056810333988042*fSkin[6]-1.190784930203603*fEdge[6]; 
  edgeSurf[14] = 1.585502557353661*fSkin[30]-0.7382874503707888*fEdge[30]-2.671875*fSkin[14]-1.453125*fEdge[14]+2.056810333988042*fSkin[7]-1.190784930203603*fEdge[7]; 
  edgeSurf[15] = 1.585502557353661*fSkin[31]-0.7382874503707888*fEdge[31]-2.671875*fSkin[15]-1.453125*fEdge[15]+2.056810333988042*fSkin[11]-1.190784930203603*fEdge[11]; 
  edgeSurf[16] = 1.190784930203603*fSkin[19]+1.190784930203603*fEdge[19]-0.9375*fSkin[16]+0.9375*fEdge[16]; 
  edgeSurf[17] = 1.190784930203603*fSkin[21]+1.190784930203603*fEdge[21]-0.9375*fSkin[17]+0.9375*fEdge[17]; 
  edgeSurf[18] = 1.190784930203603*fSkin[22]+1.190784930203603*fEdge[22]-0.9375*fSkin[18]+0.9375*fEdge[18]; 
  edgeSurf[19] = (-2.671875*fSkin[19])-1.453125*fEdge[19]+2.056810333988042*fSkin[16]-1.190784930203603*fEdge[16]; 
  edgeSurf[20] = 1.190784930203603*fSkin[23]+1.190784930203603*fEdge[23]-0.9375*fSkin[20]+0.9375*fEdge[20]; 
  edgeSurf[21] = (-2.671875*fSkin[21])-1.453125*fEdge[21]+2.056810333988042*fSkin[17]-1.190784930203603*fEdge[17]; 
  edgeSurf[22] = (-2.671875*fSkin[22])-1.453125*fEdge[22]+2.056810333988042*fSkin[18]-1.190784930203603*fEdge[18]; 
  edgeSurf[23] = (-2.671875*fSkin[23])-1.453125*fEdge[23]+2.056810333988042*fSkin[20]-1.190784930203603*fEdge[20]; 
  edgeSurf[24] = (-3.140625*fSkin[24])-0.140625*fEdge[24]+5.022775277112744*fSkin[4]+0.3025768239224545*fEdge[4]-3.773364712030896*fSkin[0]+0.4192627457812106*fEdge[0]; 
  edgeSurf[25] = (-3.140625*fSkin[25])-0.140625*fEdge[25]+5.022775277112744*fSkin[8]+0.3025768239224544*fEdge[8]-3.773364712030894*fSkin[1]+0.4192627457812105*fEdge[1]; 
  edgeSurf[26] = (-3.140625*fSkin[26])-0.140625*fEdge[26]+5.022775277112744*fSkin[9]+0.3025768239224544*fEdge[9]-3.773364712030894*fSkin[2]+0.4192627457812105*fEdge[2]; 
  edgeSurf[27] = (-3.140625*fSkin[27])-0.140625*fEdge[27]+5.022775277112744*fSkin[10]+0.3025768239224544*fEdge[10]-3.773364712030894*fSkin[3]+0.4192627457812105*fEdge[3]; 
  edgeSurf[28] = (-3.140625*fSkin[28])-0.140625*fEdge[28]+5.022775277112744*fSkin[12]+0.3025768239224545*fEdge[12]-3.773364712030896*fSkin[5]+0.4192627457812106*fEdge[5]; 
  edgeSurf[29] = (-3.140625*fSkin[29])-0.140625*fEdge[29]+5.022775277112744*fSkin[13]+0.3025768239224545*fEdge[13]-3.773364712030896*fSkin[6]+0.4192627457812106*fEdge[6]; 
  edgeSurf[30] = (-3.140625*fSkin[30])-0.140625*fEdge[30]+5.022775277112744*fSkin[14]+0.3025768239224545*fEdge[14]-3.773364712030896*fSkin[7]+0.4192627457812106*fEdge[7]; 
  edgeSurf[31] = (-3.140625*fSkin[31])-0.140625*fEdge[31]+5.022775277112744*fSkin[15]+0.3025768239224544*fEdge[15]-3.773364712030894*fSkin[11]+0.4192627457812105*fEdge[11]; 

  double boundSurf[32] = {0.0}; 
  boundSurf[4] = (-1.936491673103709*fSkin[24])-1.5*fSkin[4]-0.8660254037844386*fSkin[0]; 
  boundSurf[8] = (-1.936491673103709*fSkin[25])-1.5*fSkin[8]-0.8660254037844386*fSkin[1]; 
  boundSurf[9] = (-1.936491673103709*fSkin[26])-1.5*fSkin[9]-0.8660254037844386*fSkin[2]; 
  boundSurf[10] = (-1.936491673103709*fSkin[27])-1.5*fSkin[10]-0.8660254037844386*fSkin[3]; 
  boundSurf[12] = (-1.936491673103709*fSkin[28])-1.5*fSkin[12]-0.8660254037844386*fSkin[5]; 
  boundSurf[13] = (-1.936491673103709*fSkin[29])-1.5*fSkin[13]-0.8660254037844386*fSkin[6]; 
  boundSurf[14] = (-1.936491673103709*fSkin[30])-1.5*fSkin[14]-0.8660254037844386*fSkin[7]; 
  boundSurf[15] = (-1.936491673103709*fSkin[31])-1.5*fSkin[15]-0.8660254037844386*fSkin[11]; 
  boundSurf[19] = (-1.5*fSkin[19])-0.8660254037844387*fSkin[16]; 
  boundSurf[21] = (-1.5*fSkin[21])-0.8660254037844387*fSkin[17]; 
  boundSurf[22] = (-1.5*fSkin[22])-0.8660254037844387*fSkin[18]; 
  boundSurf[23] = (-1.5*fSkin[23])-0.8660254037844387*fSkin[20]; 
  boundSurf[24] = (-7.5*fSkin[24])-5.809475019311125*fSkin[4]-3.354101966249685*fSkin[0]; 
  boundSurf[25] = (-7.5*fSkin[25])-5.809475019311126*fSkin[8]-3.354101966249684*fSkin[1]; 
  boundSurf[26] = (-7.5*fSkin[26])-5.809475019311126*fSkin[9]-3.354101966249684*fSkin[2]; 
  boundSurf[27] = (-7.5*fSkin[27])-5.809475019311126*fSkin[10]-3.354101966249684*fSkin[3]; 
  boundSurf[28] = (-7.5*fSkin[28])-5.809475019311125*fSkin[12]-3.354101966249685*fSkin[5]; 
  boundSurf[29] = (-7.5*fSkin[29])-5.809475019311125*fSkin[13]-3.354101966249685*fSkin[6]; 
  boundSurf[30] = (-7.5*fSkin[30])-5.809475019311125*fSkin[14]-3.354101966249685*fSkin[7]; 
  boundSurf[31] = (-7.5*fSkin[31])-5.809475019311126*fSkin[15]-3.354101966249684*fSkin[11]; 

  edgeSurf_incr[0] = 0.5*nuVtSqSum[3]*edgeSurf[5]+0.5*edgeSurf[2]*nuVtSqSum[2]+0.5*edgeSurf[1]*nuVtSqSum[1]+0.5*edgeSurf[0]*nuVtSqSum[0]; 
  edgeSurf_incr[1] = 0.5*nuVtSqSum[2]*edgeSurf[5]+0.5*edgeSurf[2]*nuVtSqSum[3]+0.5*edgeSurf[0]*nuVtSqSum[1]+0.5*nuVtSqSum[0]*edgeSurf[1]; 
  edgeSurf_incr[2] = 0.5*nuVtSqSum[1]*edgeSurf[5]+0.5*edgeSurf[1]*nuVtSqSum[3]+0.5*edgeSurf[0]*nuVtSqSum[2]+0.5*nuVtSqSum[0]*edgeSurf[2]; 
  edgeSurf_incr[3] = 0.5*nuVtSqSum[3]*edgeSurf[11]+0.5*nuVtSqSum[2]*edgeSurf[7]+0.5*nuVtSqSum[1]*edgeSurf[6]+0.5*nuVtSqSum[0]*edgeSurf[3]; 
  edgeSurf_incr[4] = 0.5*nuVtSqSum[3]*edgeSurf[12]+0.5*nuVtSqSum[2]*edgeSurf[9]+0.5*nuVtSqSum[1]*edgeSurf[8]+0.5*nuVtSqSum[0]*edgeSurf[4]; 
  edgeSurf_incr[5] = 0.5*nuVtSqSum[0]*edgeSurf[5]+0.5*edgeSurf[0]*nuVtSqSum[3]+0.5*edgeSurf[1]*nuVtSqSum[2]+0.5*nuVtSqSum[1]*edgeSurf[2]; 
  edgeSurf_incr[6] = 0.5*nuVtSqSum[2]*edgeSurf[11]+0.5*nuVtSqSum[3]*edgeSurf[7]+0.5*nuVtSqSum[0]*edgeSurf[6]+0.5*nuVtSqSum[1]*edgeSurf[3]; 
  edgeSurf_incr[7] = 0.5*nuVtSqSum[1]*edgeSurf[11]+0.5*nuVtSqSum[0]*edgeSurf[7]+0.5*nuVtSqSum[3]*edgeSurf[6]+0.5*nuVtSqSum[2]*edgeSurf[3]; 
  edgeSurf_incr[8] = 0.5*nuVtSqSum[2]*edgeSurf[12]+0.5*nuVtSqSum[3]*edgeSurf[9]+0.5*nuVtSqSum[0]*edgeSurf[8]+0.5*nuVtSqSum[1]*edgeSurf[4]; 
  edgeSurf_incr[9] = 0.5*nuVtSqSum[1]*edgeSurf[12]+0.5*nuVtSqSum[0]*edgeSurf[9]+0.5*nuVtSqSum[3]*edgeSurf[8]+0.5*nuVtSqSum[2]*edgeSurf[4]; 
  edgeSurf_incr[10] = 0.5*nuVtSqSum[3]*edgeSurf[15]+0.5*nuVtSqSum[2]*edgeSurf[14]+0.5*nuVtSqSum[1]*edgeSurf[13]+0.5*nuVtSqSum[0]*edgeSurf[10]; 
  edgeSurf_incr[11] = 0.5*nuVtSqSum[0]*edgeSurf[11]+0.5*nuVtSqSum[1]*edgeSurf[7]+0.5*nuVtSqSum[2]*edgeSurf[6]+0.5*edgeSurf[3]*nuVtSqSum[3]; 
  edgeSurf_incr[12] = 0.5*nuVtSqSum[0]*edgeSurf[12]+0.5*nuVtSqSum[1]*edgeSurf[9]+0.5*nuVtSqSum[2]*edgeSurf[8]+0.5*nuVtSqSum[3]*edgeSurf[4]; 
  edgeSurf_incr[13] = 0.5*nuVtSqSum[2]*edgeSurf[15]+0.5*nuVtSqSum[3]*edgeSurf[14]+0.5*nuVtSqSum[0]*edgeSurf[13]+0.5*nuVtSqSum[1]*edgeSurf[10]; 
  edgeSurf_incr[14] = 0.5*nuVtSqSum[1]*edgeSurf[15]+0.5*nuVtSqSum[0]*edgeSurf[14]+0.5*nuVtSqSum[3]*edgeSurf[13]+0.5*nuVtSqSum[2]*edgeSurf[10]; 
  edgeSurf_incr[15] = 0.5*nuVtSqSum[0]*edgeSurf[15]+0.5*nuVtSqSum[1]*edgeSurf[14]+0.5*nuVtSqSum[2]*edgeSurf[13]+0.5*nuVtSqSum[3]*edgeSurf[10]; 
  edgeSurf_incr[16] = 0.5*nuVtSqSum[3]*edgeSurf[20]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[18]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[17]+0.5*nuVtSqSum[0]*edgeSurf[16]; 
  edgeSurf_incr[17] = 0.5000000000000001*nuVtSqSum[2]*edgeSurf[20]+0.5*nuVtSqSum[3]*edgeSurf[18]+0.5*nuVtSqSum[0]*edgeSurf[17]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[16]; 
  edgeSurf_incr[18] = 0.5000000000000001*nuVtSqSum[1]*edgeSurf[20]+0.5*nuVtSqSum[0]*edgeSurf[18]+0.5*nuVtSqSum[3]*edgeSurf[17]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[16]; 
  edgeSurf_incr[19] = 0.5*nuVtSqSum[3]*edgeSurf[23]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[22]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[21]+0.5*nuVtSqSum[0]*edgeSurf[19]; 
  edgeSurf_incr[20] = 0.5*nuVtSqSum[0]*edgeSurf[20]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[18]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[17]+0.5*nuVtSqSum[3]*edgeSurf[16]; 
  edgeSurf_incr[21] = 0.5000000000000001*nuVtSqSum[2]*edgeSurf[23]+0.5*nuVtSqSum[3]*edgeSurf[22]+0.5*nuVtSqSum[0]*edgeSurf[21]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[19]; 
  edgeSurf_incr[22] = 0.5000000000000001*nuVtSqSum[1]*edgeSurf[23]+0.5*nuVtSqSum[0]*edgeSurf[22]+0.5*nuVtSqSum[3]*edgeSurf[21]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[19]; 
  edgeSurf_incr[23] = 0.5*nuVtSqSum[0]*edgeSurf[23]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[22]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[21]+0.5*nuVtSqSum[3]*edgeSurf[19]; 
  edgeSurf_incr[24] = 0.5*nuVtSqSum[3]*edgeSurf[28]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[26]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[25]+0.5*nuVtSqSum[0]*edgeSurf[24]; 
  edgeSurf_incr[25] = 0.5000000000000001*nuVtSqSum[2]*edgeSurf[28]+0.5*nuVtSqSum[3]*edgeSurf[26]+0.5*nuVtSqSum[0]*edgeSurf[25]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[24]; 
  edgeSurf_incr[26] = 0.5000000000000001*nuVtSqSum[1]*edgeSurf[28]+0.5*nuVtSqSum[0]*edgeSurf[26]+0.5*nuVtSqSum[3]*edgeSurf[25]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[24]; 
  edgeSurf_incr[27] = 0.5*nuVtSqSum[3]*edgeSurf[31]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[30]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[29]+0.5*nuVtSqSum[0]*edgeSurf[27]; 
  edgeSurf_incr[28] = 0.5*nuVtSqSum[0]*edgeSurf[28]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[26]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[25]+0.5*nuVtSqSum[3]*edgeSurf[24]; 
  edgeSurf_incr[29] = 0.5000000000000001*nuVtSqSum[2]*edgeSurf[31]+0.5*nuVtSqSum[3]*edgeSurf[30]+0.5*nuVtSqSum[0]*edgeSurf[29]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[27]; 
  edgeSurf_incr[30] = 0.5000000000000001*nuVtSqSum[1]*edgeSurf[31]+0.5*nuVtSqSum[0]*edgeSurf[30]+0.5*nuVtSqSum[3]*edgeSurf[29]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[27]; 
  edgeSurf_incr[31] = 0.5*nuVtSqSum[0]*edgeSurf[31]+0.5000000000000001*nuVtSqSum[1]*edgeSurf[30]+0.5000000000000001*nuVtSqSum[2]*edgeSurf[29]+0.5*nuVtSqSum[3]*edgeSurf[27]; 

  boundSurf_incr[0] = 0.5*nuVtSqSum[3]*boundSurf[5]+0.5*boundSurf[2]*nuVtSqSum[2]+0.5*boundSurf[1]*nuVtSqSum[1]+0.5*boundSurf[0]*nuVtSqSum[0]; 
  boundSurf_incr[1] = 0.5*nuVtSqSum[2]*boundSurf[5]+0.5*boundSurf[2]*nuVtSqSum[3]+0.5*boundSurf[0]*nuVtSqSum[1]+0.5*nuVtSqSum[0]*boundSurf[1]; 
  boundSurf_incr[2] = 0.5*nuVtSqSum[1]*boundSurf[5]+0.5*boundSurf[1]*nuVtSqSum[3]+0.5*boundSurf[0]*nuVtSqSum[2]+0.5*nuVtSqSum[0]*boundSurf[2]; 
  boundSurf_incr[3] = 0.5*nuVtSqSum[3]*boundSurf[11]+0.5*nuVtSqSum[2]*boundSurf[7]+0.5*nuVtSqSum[1]*boundSurf[6]+0.5*nuVtSqSum[0]*boundSurf[3]; 
  boundSurf_incr[4] = 0.5*nuVtSqSum[3]*boundSurf[12]+0.5*nuVtSqSum[2]*boundSurf[9]+0.5*nuVtSqSum[1]*boundSurf[8]+0.5*nuVtSqSum[0]*boundSurf[4]; 
  boundSurf_incr[5] = 0.5*nuVtSqSum[0]*boundSurf[5]+0.5*boundSurf[0]*nuVtSqSum[3]+0.5*boundSurf[1]*nuVtSqSum[2]+0.5*nuVtSqSum[1]*boundSurf[2]; 
  boundSurf_incr[6] = 0.5*nuVtSqSum[2]*boundSurf[11]+0.5*nuVtSqSum[3]*boundSurf[7]+0.5*nuVtSqSum[0]*boundSurf[6]+0.5*nuVtSqSum[1]*boundSurf[3]; 
  boundSurf_incr[7] = 0.5*nuVtSqSum[1]*boundSurf[11]+0.5*nuVtSqSum[0]*boundSurf[7]+0.5*nuVtSqSum[3]*boundSurf[6]+0.5*nuVtSqSum[2]*boundSurf[3]; 
  boundSurf_incr[8] = 0.5*nuVtSqSum[2]*boundSurf[12]+0.5*nuVtSqSum[3]*boundSurf[9]+0.5*nuVtSqSum[0]*boundSurf[8]+0.5*nuVtSqSum[1]*boundSurf[4]; 
  boundSurf_incr[9] = 0.5*nuVtSqSum[1]*boundSurf[12]+0.5*nuVtSqSum[0]*boundSurf[9]+0.5*nuVtSqSum[3]*boundSurf[8]+0.5*nuVtSqSum[2]*boundSurf[4]; 
  boundSurf_incr[10] = 0.5*nuVtSqSum[3]*boundSurf[15]+0.5*nuVtSqSum[2]*boundSurf[14]+0.5*nuVtSqSum[1]*boundSurf[13]+0.5*nuVtSqSum[0]*boundSurf[10]; 
  boundSurf_incr[11] = 0.5*nuVtSqSum[0]*boundSurf[11]+0.5*nuVtSqSum[1]*boundSurf[7]+0.5*nuVtSqSum[2]*boundSurf[6]+0.5*boundSurf[3]*nuVtSqSum[3]; 
  boundSurf_incr[12] = 0.5*nuVtSqSum[0]*boundSurf[12]+0.5*nuVtSqSum[1]*boundSurf[9]+0.5*nuVtSqSum[2]*boundSurf[8]+0.5*nuVtSqSum[3]*boundSurf[4]; 
  boundSurf_incr[13] = 0.5*nuVtSqSum[2]*boundSurf[15]+0.5*nuVtSqSum[3]*boundSurf[14]+0.5*nuVtSqSum[0]*boundSurf[13]+0.5*nuVtSqSum[1]*boundSurf[10]; 
  boundSurf_incr[14] = 0.5*nuVtSqSum[1]*boundSurf[15]+0.5*nuVtSqSum[0]*boundSurf[14]+0.5*nuVtSqSum[3]*boundSurf[13]+0.5*nuVtSqSum[2]*boundSurf[10]; 
  boundSurf_incr[15] = 0.5*nuVtSqSum[0]*boundSurf[15]+0.5*nuVtSqSum[1]*boundSurf[14]+0.5*nuVtSqSum[2]*boundSurf[13]+0.5*nuVtSqSum[3]*boundSurf[10]; 
  boundSurf_incr[16] = 0.5*nuVtSqSum[3]*boundSurf[20]+0.5000000000000001*nuVtSqSum[2]*boundSurf[18]+0.5000000000000001*nuVtSqSum[1]*boundSurf[17]+0.5*nuVtSqSum[0]*boundSurf[16]; 
  boundSurf_incr[17] = 0.5000000000000001*nuVtSqSum[2]*boundSurf[20]+0.5*nuVtSqSum[3]*boundSurf[18]+0.5*nuVtSqSum[0]*boundSurf[17]+0.5000000000000001*nuVtSqSum[1]*boundSurf[16]; 
  boundSurf_incr[18] = 0.5000000000000001*nuVtSqSum[1]*boundSurf[20]+0.5*nuVtSqSum[0]*boundSurf[18]+0.5*nuVtSqSum[3]*boundSurf[17]+0.5000000000000001*nuVtSqSum[2]*boundSurf[16]; 
  boundSurf_incr[19] = 0.5*nuVtSqSum[3]*boundSurf[23]+0.5000000000000001*nuVtSqSum[2]*boundSurf[22]+0.5000000000000001*nuVtSqSum[1]*boundSurf[21]+0.5*nuVtSqSum[0]*boundSurf[19]; 
  boundSurf_incr[20] = 0.5*nuVtSqSum[0]*boundSurf[20]+0.5000000000000001*nuVtSqSum[1]*boundSurf[18]+0.5000000000000001*nuVtSqSum[2]*boundSurf[17]+0.5*nuVtSqSum[3]*boundSurf[16]; 
  boundSurf_incr[21] = 0.5000000000000001*nuVtSqSum[2]*boundSurf[23]+0.5*nuVtSqSum[3]*boundSurf[22]+0.5*nuVtSqSum[0]*boundSurf[21]+0.5000000000000001*nuVtSqSum[1]*boundSurf[19]; 
  boundSurf_incr[22] = 0.5000000000000001*nuVtSqSum[1]*boundSurf[23]+0.5*nuVtSqSum[0]*boundSurf[22]+0.5*nuVtSqSum[3]*boundSurf[21]+0.5000000000000001*nuVtSqSum[2]*boundSurf[19]; 
  boundSurf_incr[23] = 0.5*nuVtSqSum[0]*boundSurf[23]+0.5000000000000001*nuVtSqSum[1]*boundSurf[22]+0.5000000000000001*nuVtSqSum[2]*boundSurf[21]+0.5*nuVtSqSum[3]*boundSurf[19]; 
  boundSurf_incr[24] = 0.5*nuVtSqSum[3]*boundSurf[28]+0.5000000000000001*nuVtSqSum[2]*boundSurf[26]+0.5000000000000001*nuVtSqSum[1]*boundSurf[25]+0.5*nuVtSqSum[0]*boundSurf[24]; 
  boundSurf_incr[25] = 0.5000000000000001*nuVtSqSum[2]*boundSurf[28]+0.5*nuVtSqSum[3]*boundSurf[26]+0.5*nuVtSqSum[0]*boundSurf[25]+0.5000000000000001*nuVtSqSum[1]*boundSurf[24]; 
  boundSurf_incr[26] = 0.5000000000000001*nuVtSqSum[1]*boundSurf[28]+0.5*nuVtSqSum[0]*boundSurf[26]+0.5*nuVtSqSum[3]*boundSurf[25]+0.5000000000000001*nuVtSqSum[2]*boundSurf[24]; 
  boundSurf_incr[27] = 0.5*nuVtSqSum[3]*boundSurf[31]+0.5000000000000001*nuVtSqSum[2]*boundSurf[30]+0.5000000000000001*nuVtSqSum[1]*boundSurf[29]+0.5*nuVtSqSum[0]*boundSurf[27]; 
  boundSurf_incr[28] = 0.5*nuVtSqSum[0]*boundSurf[28]+0.5000000000000001*nuVtSqSum[1]*boundSurf[26]+0.5000000000000001*nuVtSqSum[2]*boundSurf[25]+0.5*nuVtSqSum[3]*boundSurf[24]; 
  boundSurf_incr[29] = 0.5000000000000001*nuVtSqSum[2]*boundSurf[31]+0.5*nuVtSqSum[3]*boundSurf[30]+0.5*nuVtSqSum[0]*boundSurf[29]+0.5000000000000001*nuVtSqSum[1]*boundSurf[27]; 
  boundSurf_incr[30] = 0.5000000000000001*nuVtSqSum[1]*boundSurf[31]+0.5*nuVtSqSum[0]*boundSurf[30]+0.5*nuVtSqSum[3]*boundSurf[29]+0.5000000000000001*nuVtSqSum[2]*boundSurf[27]; 
  boundSurf_incr[31] = 0.5*nuVtSqSum[0]*boundSurf[31]+0.5000000000000001*nuVtSqSum[1]*boundSurf[30]+0.5000000000000001*nuVtSqSum[2]*boundSurf[29]+0.5*nuVtSqSum[3]*boundSurf[27]; 

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
  out[16] += (vol_incr[16]+edgeSurf_incr[16]+boundSurf_incr[16])*rdvSq4; 
  out[17] += (vol_incr[17]+edgeSurf_incr[17]+boundSurf_incr[17])*rdvSq4; 
  out[18] += (vol_incr[18]+edgeSurf_incr[18]+boundSurf_incr[18])*rdvSq4; 
  out[19] += (vol_incr[19]+edgeSurf_incr[19]+boundSurf_incr[19])*rdvSq4; 
  out[20] += (vol_incr[20]+edgeSurf_incr[20]+boundSurf_incr[20])*rdvSq4; 
  out[21] += (vol_incr[21]+edgeSurf_incr[21]+boundSurf_incr[21])*rdvSq4; 
  out[22] += (vol_incr[22]+edgeSurf_incr[22]+boundSurf_incr[22])*rdvSq4; 
  out[23] += (vol_incr[23]+edgeSurf_incr[23]+boundSurf_incr[23])*rdvSq4; 
  out[24] += (vol_incr[24]+edgeSurf_incr[24]+boundSurf_incr[24])*rdvSq4; 
  out[25] += (vol_incr[25]+edgeSurf_incr[25]+boundSurf_incr[25])*rdvSq4; 
  out[26] += (vol_incr[26]+edgeSurf_incr[26]+boundSurf_incr[26])*rdvSq4; 
  out[27] += (vol_incr[27]+edgeSurf_incr[27]+boundSurf_incr[27])*rdvSq4; 
  out[28] += (vol_incr[28]+edgeSurf_incr[28]+boundSurf_incr[28])*rdvSq4; 
  out[29] += (vol_incr[29]+edgeSurf_incr[29]+boundSurf_incr[29])*rdvSq4; 
  out[30] += (vol_incr[30]+edgeSurf_incr[30]+boundSurf_incr[30])*rdvSq4; 
  out[31] += (vol_incr[31]+edgeSurf_incr[31]+boundSurf_incr[31])*rdvSq4; 

  return 0.;

} 
