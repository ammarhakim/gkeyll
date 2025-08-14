#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_boundary_surfvx_1x1v_ser_p3(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[1]; 

  const double *Ghat = &flux[0]; 
  if (edge == -1) { 

  out[0] += -(0.7071067811865475*Ghat[0]*dv10); 
  out[1] += -(0.7071067811865475*Ghat[1]*dv10); 
  out[2] += -(1.224744871391589*Ghat[0]*dv10); 
  out[3] += -(1.224744871391589*Ghat[1]*dv10); 
  out[4] += -(0.7071067811865475*Ghat[2]*dv10); 
  out[5] += -(1.5811388300841895*Ghat[0]*dv10); 
  out[6] += -(1.224744871391589*Ghat[2]*dv10); 
  out[7] += -(1.5811388300841898*Ghat[1]*dv10); 
  out[8] += -(0.7071067811865475*Ghat[3]*dv10); 
  out[9] += -(1.8708286933869707*Ghat[0]*dv10); 
  out[10] += -(1.2247448713915887*Ghat[3]*dv10); 
  out[11] += -(1.8708286933869707*Ghat[1]*dv10); 

  } else { 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -(1.224744871391589*Ghat[0]*dv10); 
  out[3] += -(1.224744871391589*Ghat[1]*dv10); 
  out[4] += 0.7071067811865475*Ghat[2]*dv10; 
  out[5] += 1.5811388300841895*Ghat[0]*dv10; 
  out[6] += -(1.224744871391589*Ghat[2]*dv10); 
  out[7] += 1.5811388300841898*Ghat[1]*dv10; 
  out[8] += 0.7071067811865475*Ghat[3]*dv10; 
  out[9] += -(1.8708286933869707*Ghat[0]*dv10); 
  out[10] += -(1.2247448713915887*Ghat[3]*dv10); 
  out[11] += -(1.8708286933869707*Ghat[1]*dv10); 

  } 
  return 0.0;

} 
