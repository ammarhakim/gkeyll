#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_boundary_surfvx_1x1v_tensor_p2(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[1]; 

  const double *Fhat_nodal = &flux[0]; 
  if (edge == -1) { 

  out[0] += ((-0.1739274225687269*Fhat_nodal[3])-0.3260725774312731*(Fhat_nodal[2]+Fhat_nodal[1])-0.1739274225687269*Fhat_nodal[0])*dv10; 
  out[1] += ((-0.2594182892927712*Fhat_nodal[3])-0.192012546066862*Fhat_nodal[2]+0.192012546066862*Fhat_nodal[1]+0.2594182892927712*Fhat_nodal[0])*dv10; 
  out[2] += ((-0.3012511327185366*Fhat_nodal[3])-0.5647742710659018*(Fhat_nodal[2]+Fhat_nodal[1])-0.3012511327185366*Fhat_nodal[0])*dv10; 
  out[3] += ((-0.449325657467681*Fhat_nodal[3])-0.3325754854784645*Fhat_nodal[2]+0.3325754854784645*Fhat_nodal[1]+0.449325657467681*Fhat_nodal[0])*dv10; 
  out[4] += ((-0.2381448361039205*Fhat_nodal[3])+0.2381448361039205*(Fhat_nodal[2]+Fhat_nodal[1])-0.2381448361039205*Fhat_nodal[0])*dv10; 
  out[5] += ((-0.3889135400150045*Fhat_nodal[3])-0.7291204487348905*(Fhat_nodal[2]+Fhat_nodal[1])-0.3889135400150045*Fhat_nodal[0])*dv10; 
  out[6] += ((-0.4124789556921534*Fhat_nodal[3])+0.4124789556921534*(Fhat_nodal[2]+Fhat_nodal[1])-0.4124789556921534*Fhat_nodal[0])*dv10; 
  out[7] += ((-0.5800769294653423*Fhat_nodal[3])-0.4293531055383132*Fhat_nodal[2]+0.4293531055383132*Fhat_nodal[1]+0.5800769294653423*Fhat_nodal[0])*dv10; 
  out[8] += ((-0.5325080420189123*Fhat_nodal[3])+0.5325080420189123*(Fhat_nodal[2]+Fhat_nodal[1])-0.5325080420189123*Fhat_nodal[0])*dv10; 

  } else { 

  out[0] += (0.1739274225687269*Fhat_nodal[3]+0.3260725774312731*(Fhat_nodal[2]+Fhat_nodal[1])+0.1739274225687269*Fhat_nodal[0])*dv10; 
  out[1] += (0.2594182892927712*Fhat_nodal[3]+0.192012546066862*Fhat_nodal[2]-0.192012546066862*Fhat_nodal[1]-0.2594182892927712*Fhat_nodal[0])*dv10; 
  out[2] += ((-0.3012511327185366*Fhat_nodal[3])-0.5647742710659018*(Fhat_nodal[2]+Fhat_nodal[1])-0.3012511327185366*Fhat_nodal[0])*dv10; 
  out[3] += ((-0.449325657467681*Fhat_nodal[3])-0.3325754854784645*Fhat_nodal[2]+0.3325754854784645*Fhat_nodal[1]+0.449325657467681*Fhat_nodal[0])*dv10; 
  out[4] += (0.2381448361039205*Fhat_nodal[3]-0.2381448361039205*(Fhat_nodal[2]+Fhat_nodal[1])+0.2381448361039205*Fhat_nodal[0])*dv10; 
  out[5] += (0.3889135400150045*Fhat_nodal[3]+0.7291204487348905*(Fhat_nodal[2]+Fhat_nodal[1])+0.3889135400150045*Fhat_nodal[0])*dv10; 
  out[6] += ((-0.4124789556921534*Fhat_nodal[3])+0.4124789556921534*(Fhat_nodal[2]+Fhat_nodal[1])-0.4124789556921534*Fhat_nodal[0])*dv10; 
  out[7] += (0.5800769294653423*Fhat_nodal[3]+0.4293531055383132*Fhat_nodal[2]-0.4293531055383132*Fhat_nodal[1]-0.5800769294653423*Fhat_nodal[0])*dv10; 
  out[8] += (0.5325080420189123*Fhat_nodal[3]-0.5325080420189123*(Fhat_nodal[2]+Fhat_nodal[1])+0.5325080420189123*Fhat_nodal[0])*dv10; 

  } 
  return 0.0;

} 
