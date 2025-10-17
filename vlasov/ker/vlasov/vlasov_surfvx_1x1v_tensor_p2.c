#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_surfvx_1x1v_tensor_p2(const double *w, const double *dxv,
  const double *Fhat_l_nodal, const double *Fhat_r_nodal, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[1]; 

  const double *Fhat_l_nodal_c = &Fhat_l_nodal[0]; 
  const double *Fhat_r_nodal_c = &Fhat_r_nodal[0]; 
  out[0] += ((-0.1739274225687269*Fhat_r_nodal_c[3])+0.1739274225687269*Fhat_l_nodal_c[3]-0.3260725774312731*Fhat_r_nodal_c[2]+0.3260725774312731*Fhat_l_nodal_c[2]-0.3260725774312731*Fhat_r_nodal_c[1]+0.3260725774312731*Fhat_l_nodal_c[1]-0.1739274225687269*Fhat_r_nodal_c[0]+0.1739274225687269*Fhat_l_nodal_c[0])*dv10; 
  out[1] += ((-0.2594182892927712*Fhat_r_nodal_c[3])+0.2594182892927712*Fhat_l_nodal_c[3]-0.192012546066862*Fhat_r_nodal_c[2]+0.192012546066862*(Fhat_l_nodal_c[2]+Fhat_r_nodal_c[1])-0.192012546066862*Fhat_l_nodal_c[1]+0.2594182892927712*Fhat_r_nodal_c[0]-0.2594182892927712*Fhat_l_nodal_c[0])*dv10; 
  out[2] += ((-0.3012511327185366*(Fhat_r_nodal_c[3]+Fhat_l_nodal_c[3]))-0.5647742710659018*(Fhat_r_nodal_c[2]+Fhat_l_nodal_c[2]+Fhat_r_nodal_c[1]+Fhat_l_nodal_c[1])-0.3012511327185366*(Fhat_r_nodal_c[0]+Fhat_l_nodal_c[0]))*dv10; 
  out[3] += ((-0.449325657467681*(Fhat_r_nodal_c[3]+Fhat_l_nodal_c[3]))-0.3325754854784645*(Fhat_r_nodal_c[2]+Fhat_l_nodal_c[2])+0.3325754854784645*(Fhat_r_nodal_c[1]+Fhat_l_nodal_c[1])+0.449325657467681*(Fhat_r_nodal_c[0]+Fhat_l_nodal_c[0]))*dv10; 
  out[4] += ((-0.2381448361039205*Fhat_r_nodal_c[3])+0.2381448361039205*(Fhat_l_nodal_c[3]+Fhat_r_nodal_c[2])-0.2381448361039205*Fhat_l_nodal_c[2]+0.2381448361039205*Fhat_r_nodal_c[1]-0.2381448361039205*(Fhat_l_nodal_c[1]+Fhat_r_nodal_c[0])+0.2381448361039205*Fhat_l_nodal_c[0])*dv10; 
  out[5] += ((-0.3889135400150045*Fhat_r_nodal_c[3])+0.3889135400150045*Fhat_l_nodal_c[3]-0.7291204487348905*Fhat_r_nodal_c[2]+0.7291204487348905*Fhat_l_nodal_c[2]-0.7291204487348905*Fhat_r_nodal_c[1]+0.7291204487348905*Fhat_l_nodal_c[1]-0.3889135400150045*Fhat_r_nodal_c[0]+0.3889135400150045*Fhat_l_nodal_c[0])*dv10; 
  out[6] += ((-0.4124789556921534*(Fhat_r_nodal_c[3]+Fhat_l_nodal_c[3]))+0.4124789556921534*(Fhat_r_nodal_c[2]+Fhat_l_nodal_c[2]+Fhat_r_nodal_c[1]+Fhat_l_nodal_c[1])-0.4124789556921534*(Fhat_r_nodal_c[0]+Fhat_l_nodal_c[0]))*dv10; 
  out[7] += ((-0.5800769294653423*Fhat_r_nodal_c[3])+0.5800769294653423*Fhat_l_nodal_c[3]-0.4293531055383132*Fhat_r_nodal_c[2]+0.4293531055383132*(Fhat_l_nodal_c[2]+Fhat_r_nodal_c[1])-0.4293531055383132*Fhat_l_nodal_c[1]+0.5800769294653423*Fhat_r_nodal_c[0]-0.5800769294653423*Fhat_l_nodal_c[0])*dv10; 
  out[8] += ((-0.5325080420189123*Fhat_r_nodal_c[3])+0.5325080420189123*(Fhat_l_nodal_c[3]+Fhat_r_nodal_c[2])-0.5325080420189123*Fhat_l_nodal_c[2]+0.5325080420189123*Fhat_r_nodal_c[1]-0.5325080420189123*(Fhat_l_nodal_c[1]+Fhat_r_nodal_c[0])+0.5325080420189123*Fhat_l_nodal_c[0])*dv10; 

  return 0.0;

} 
