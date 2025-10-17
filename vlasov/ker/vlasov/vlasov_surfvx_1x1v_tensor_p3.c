#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_surfvx_1x1v_tensor_p3(const double *w, const double *dxv,
  const double *Fhat_l_nodal, const double *Fhat_r_nodal, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[1]; 

  const double *Fhat_l_nodal_c = &Fhat_l_nodal[0]; 
  const double *Fhat_r_nodal_c = &Fhat_r_nodal[0]; 
  out[0] += ((-0.1184634425280946*Fhat_r_nodal_c[4])+0.1184634425280946*Fhat_l_nodal_c[4]-0.2393143352496831*Fhat_r_nodal_c[3]+0.2393143352496831*Fhat_l_nodal_c[3]-0.2844444444444445*Fhat_r_nodal_c[2]+0.2844444444444445*Fhat_l_nodal_c[2]-0.2393143352496831*Fhat_r_nodal_c[1]+0.2393143352496831*Fhat_l_nodal_c[1]-0.1184634425280946*Fhat_r_nodal_c[0]+0.1184634425280946*Fhat_l_nodal_c[0])*dv10; 
  out[1] += ((-0.1859342410113515*Fhat_r_nodal_c[4])+0.1859342410113515*Fhat_l_nodal_c[4]-0.2231979993378559*Fhat_r_nodal_c[3]+0.2231979993378559*(Fhat_l_nodal_c[3]+Fhat_r_nodal_c[1])-0.2231979993378559*Fhat_l_nodal_c[1]+0.1859342410113515*Fhat_r_nodal_c[0]-0.1859342410113515*Fhat_l_nodal_c[0])*dv10; 
  out[2] += ((-0.2051847012981755*(Fhat_r_nodal_c[4]+Fhat_l_nodal_c[4]))-0.4145045876320226*(Fhat_r_nodal_c[3]+Fhat_l_nodal_c[3])-0.4926722297084807*(Fhat_r_nodal_c[2]+Fhat_l_nodal_c[2])-0.4145045876320226*(Fhat_r_nodal_c[1]+Fhat_l_nodal_c[1])-0.2051847012981755*(Fhat_r_nodal_c[0]+Fhat_l_nodal_c[0]))*dv10; 
  out[3] += ((-0.3220475522984177*(Fhat_r_nodal_c[4]+Fhat_l_nodal_c[4]))-0.3865902750008912*(Fhat_r_nodal_c[3]+Fhat_l_nodal_c[3])+0.3865902750008912*(Fhat_r_nodal_c[1]+Fhat_l_nodal_c[1])+0.3220475522984177*(Fhat_r_nodal_c[0]+Fhat_l_nodal_c[0]))*dv10; 
  out[4] += ((-0.1938330593514186*Fhat_r_nodal_c[4])+0.1938330593514186*Fhat_l_nodal_c[4]+0.03482378095143357*Fhat_r_nodal_c[3]-0.03482378095143357*Fhat_l_nodal_c[3]+0.3180185567999701*Fhat_r_nodal_c[2]-0.3180185567999701*Fhat_l_nodal_c[2]+0.03482378095143357*Fhat_r_nodal_c[1]-0.03482378095143357*Fhat_l_nodal_c[1]-0.1938330593514186*Fhat_r_nodal_c[0]+0.1938330593514186*Fhat_l_nodal_c[0])*dv10; 
  out[5] += ((-0.264892310341459*Fhat_r_nodal_c[4])+0.264892310341459*Fhat_l_nodal_c[4]-0.5351231216084656*Fhat_r_nodal_c[3]+0.5351231216084656*Fhat_l_nodal_c[3]-0.6360371135999403*Fhat_r_nodal_c[2]+0.6360371135999403*Fhat_l_nodal_c[2]-0.5351231216084656*Fhat_r_nodal_c[1]+0.5351231216084656*Fhat_l_nodal_c[1]-0.264892310341459*Fhat_r_nodal_c[0]+0.264892310341459*Fhat_l_nodal_c[0])*dv10; 
  out[6] += ((-0.3357287069831708*(Fhat_r_nodal_c[4]+Fhat_l_nodal_c[4]))+0.06031655791953221*(Fhat_r_nodal_c[3]+Fhat_l_nodal_c[3])+0.5508242981272772*(Fhat_r_nodal_c[2]+Fhat_l_nodal_c[2])+0.06031655791953221*(Fhat_r_nodal_c[1]+Fhat_l_nodal_c[1])-0.3357287069831708*(Fhat_r_nodal_c[0]+Fhat_l_nodal_c[0]))*dv10; 
  out[7] += ((-0.4157616022462114*Fhat_r_nodal_c[4])+0.4157616022462114*Fhat_l_nodal_c[4]-0.4990858989613992*Fhat_r_nodal_c[3]+0.4990858989613992*(Fhat_l_nodal_c[3]+Fhat_r_nodal_c[1])-0.4990858989613992*Fhat_l_nodal_c[1]+0.4157616022462114*Fhat_r_nodal_c[0]-0.4157616022462114*Fhat_l_nodal_c[0])*dv10; 
  out[8] += ((-0.1570355987780474*Fhat_r_nodal_c[4])+0.1570355987780474*Fhat_l_nodal_c[4]+0.2642722473443246*Fhat_r_nodal_c[3]-0.2642722473443246*(Fhat_l_nodal_c[3]+Fhat_r_nodal_c[1])+0.2642722473443246*Fhat_l_nodal_c[1]+0.1570355987780474*Fhat_r_nodal_c[0]-0.1570355987780474*Fhat_l_nodal_c[0])*dv10; 
  out[9] += ((-0.3134248083819311*(Fhat_r_nodal_c[4]+Fhat_l_nodal_c[4]))-0.6331662162434001*(Fhat_r_nodal_c[3]+Fhat_l_nodal_c[3])-0.752569261813928*(Fhat_r_nodal_c[2]+Fhat_l_nodal_c[2])-0.6331662162434001*(Fhat_r_nodal_c[1]+Fhat_l_nodal_c[1])-0.3134248083819311*(Fhat_r_nodal_c[0]+Fhat_l_nodal_c[0]))*dv10; 
  out[10] += ((-0.4334238969965234*Fhat_r_nodal_c[4])+0.4334238969965234*Fhat_l_nodal_c[4]+0.07786834144096777*Fhat_r_nodal_c[3]-0.07786834144096777*Fhat_l_nodal_c[3]+0.7111111111111112*Fhat_r_nodal_c[2]-0.7111111111111112*Fhat_l_nodal_c[2]+0.07786834144096777*Fhat_r_nodal_c[1]-0.07786834144096777*Fhat_l_nodal_c[1]-0.4334238969965234*Fhat_r_nodal_c[0]+0.4334238969965234*Fhat_l_nodal_c[0])*dv10; 
  out[11] += ((-0.2719936356805793*(Fhat_r_nodal_c[4]+Fhat_l_nodal_c[4]))+0.4577329594307796*(Fhat_r_nodal_c[3]+Fhat_l_nodal_c[3])-0.4577329594307796*(Fhat_r_nodal_c[1]+Fhat_l_nodal_c[1])+0.2719936356805793*(Fhat_r_nodal_c[0]+Fhat_l_nodal_c[0]))*dv10; 
  out[12] += ((-0.4919357619275829*(Fhat_r_nodal_c[4]+Fhat_l_nodal_c[4]))-0.5905263993751261*(Fhat_r_nodal_c[3]+Fhat_l_nodal_c[3])+0.5905263993751261*(Fhat_r_nodal_c[1]+Fhat_l_nodal_c[1])+0.4919357619275829*(Fhat_r_nodal_c[0]+Fhat_l_nodal_c[0]))*dv10; 
  out[13] += ((-0.351142273755097*Fhat_r_nodal_c[4])+0.351142273755097*Fhat_l_nodal_c[4]+0.5909307096285482*Fhat_r_nodal_c[3]-0.5909307096285482*(Fhat_l_nodal_c[3]+Fhat_r_nodal_c[1])+0.5909307096285482*Fhat_l_nodal_c[1]+0.351142273755097*Fhat_r_nodal_c[0]-0.351142273755097*Fhat_l_nodal_c[0])*dv10; 
  out[14] += ((-0.5128340709066764*(Fhat_r_nodal_c[4]+Fhat_l_nodal_c[4]))+0.09213506410848148*(Fhat_r_nodal_c[3]+Fhat_l_nodal_c[3])+0.84139801359639*(Fhat_r_nodal_c[2]+Fhat_l_nodal_c[2])+0.09213506410848148*(Fhat_r_nodal_c[1]+Fhat_l_nodal_c[1])-0.5128340709066764*(Fhat_r_nodal_c[0]+Fhat_l_nodal_c[0]))*dv10; 
  out[15] += ((-0.415477141350832*(Fhat_r_nodal_c[4]+Fhat_l_nodal_c[4]))+0.6991986448892327*(Fhat_r_nodal_c[3]+Fhat_l_nodal_c[3])-0.6991986448892327*(Fhat_r_nodal_c[1]+Fhat_l_nodal_c[1])+0.415477141350832*(Fhat_r_nodal_c[0]+Fhat_l_nodal_c[0]))*dv10; 

  return 0.0;

} 
