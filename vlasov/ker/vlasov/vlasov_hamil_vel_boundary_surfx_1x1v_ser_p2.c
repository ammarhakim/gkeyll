#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_boundary_surfx_1x1v_ser_p2(const double *w, const double *dxv, 
  const double *hamil, 
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 
  const double wv = w[1]; 

  double Ghat[3]; 

  if (edge == -1) { 

  if (wv>0) { 

  Ghat[0] = 4.330127018922194*hamil[2]*fskin[6]+1.9364916731037085*hamil[1]*fskin[4]+3.3541019662496847*hamil[2]*fskin[3]+1.9364916731037085*fskin[2]*hamil[2]+1.5*fskin[1]*hamil[1]+0.8660254037844386*fskin[0]*hamil[1]; 
  Ghat[1] = 3.0*hamil[2]*fskin[7]+1.9364916731037085*hamil[1]*fskin[6]+1.7320508075688772*hamil[2]*fskin[5]+4.330127018922193*hamil[2]*fskin[4]+1.5*hamil[1]*fskin[3]+3.3541019662496847*fskin[1]*hamil[2]+1.9364916731037085*fskin[0]*hamil[2]+0.8660254037844386*hamil[1]*fskin[2]; 
  Ghat[2] = 1.5*hamil[1]*fskin[7]+3.872983346207417*hamil[2]*fskin[6]+0.8660254037844386*hamil[1]*fskin[5]+3.0*hamil[2]*fskin[3]+1.7320508075688772*fskin[2]*hamil[2]; 

  } else { 

  Ghat[0] = 4.330127018922194*hamil[2]*fedge[6]+1.9364916731037085*hamil[1]*fedge[4]-3.3541019662496847*hamil[2]*fedge[3]+1.9364916731037085*fedge[2]*hamil[2]-1.5*fedge[1]*hamil[1]+0.8660254037844386*fedge[0]*hamil[1]; 
  Ghat[1] = -(3.0*hamil[2]*fedge[7])+1.9364916731037085*hamil[1]*fedge[6]+1.7320508075688772*hamil[2]*fedge[5]+4.330127018922193*hamil[2]*fedge[4]-1.5*hamil[1]*fedge[3]-3.3541019662496847*fedge[1]*hamil[2]+1.9364916731037085*fedge[0]*hamil[2]+0.8660254037844386*hamil[1]*fedge[2]; 
  Ghat[2] = -(1.5*hamil[1]*fedge[7])+3.872983346207417*hamil[2]*fedge[6]+0.8660254037844386*hamil[1]*fedge[5]-3.0*hamil[2]*fedge[3]+1.7320508075688772*fedge[2]*hamil[2]; 

  } 

  out[0] += -(0.7071067811865475*Ghat[0]*dv10*dx10); 
  out[1] += -(1.224744871391589*Ghat[0]*dv10*dx10); 
  out[2] += -(0.7071067811865475*Ghat[1]*dv10*dx10); 
  out[3] += -(1.224744871391589*Ghat[1]*dv10*dx10); 
  out[4] += -(1.5811388300841895*Ghat[0]*dv10*dx10); 
  out[5] += -(0.7071067811865475*Ghat[2]*dv10*dx10); 
  out[6] += -(1.5811388300841898*Ghat[1]*dv10*dx10); 
  out[7] += -(1.224744871391589*Ghat[2]*dv10*dx10); 

  } else { 

  if (wv>0) { 

  Ghat[0] = 4.330127018922194*hamil[2]*fedge[6]+1.9364916731037085*hamil[1]*fedge[4]+3.3541019662496847*hamil[2]*fedge[3]+1.9364916731037085*fedge[2]*hamil[2]+1.5*fedge[1]*hamil[1]+0.8660254037844386*fedge[0]*hamil[1]; 
  Ghat[1] = 3.0*hamil[2]*fedge[7]+1.9364916731037085*hamil[1]*fedge[6]+1.7320508075688772*hamil[2]*fedge[5]+4.330127018922193*hamil[2]*fedge[4]+1.5*hamil[1]*fedge[3]+3.3541019662496847*fedge[1]*hamil[2]+1.9364916731037085*fedge[0]*hamil[2]+0.8660254037844386*hamil[1]*fedge[2]; 
  Ghat[2] = 1.5*hamil[1]*fedge[7]+3.872983346207417*hamil[2]*fedge[6]+0.8660254037844386*hamil[1]*fedge[5]+3.0*hamil[2]*fedge[3]+1.7320508075688772*fedge[2]*hamil[2]; 

  } else { 

  Ghat[0] = 4.330127018922194*hamil[2]*fskin[6]+1.9364916731037085*hamil[1]*fskin[4]-3.3541019662496847*hamil[2]*fskin[3]+1.9364916731037085*fskin[2]*hamil[2]-1.5*fskin[1]*hamil[1]+0.8660254037844386*fskin[0]*hamil[1]; 
  Ghat[1] = -(3.0*hamil[2]*fskin[7])+1.9364916731037085*hamil[1]*fskin[6]+1.7320508075688772*hamil[2]*fskin[5]+4.330127018922193*hamil[2]*fskin[4]-1.5*hamil[1]*fskin[3]-3.3541019662496847*fskin[1]*hamil[2]+1.9364916731037085*fskin[0]*hamil[2]+0.8660254037844386*hamil[1]*fskin[2]; 
  Ghat[2] = -(1.5*hamil[1]*fskin[7])+3.872983346207417*hamil[2]*fskin[6]+0.8660254037844386*hamil[1]*fskin[5]-3.0*hamil[2]*fskin[3]+1.7320508075688772*fskin[2]*hamil[2]; 

  } 

  out[0] += 0.7071067811865475*Ghat[0]*dv10*dx10; 
  out[1] += -(1.224744871391589*Ghat[0]*dv10*dx10); 
  out[2] += 0.7071067811865475*Ghat[1]*dv10*dx10; 
  out[3] += -(1.224744871391589*Ghat[1]*dv10*dx10); 
  out[4] += 1.5811388300841895*Ghat[0]*dv10*dx10; 
  out[5] += 0.7071067811865475*Ghat[2]*dv10*dx10; 
  out[6] += 1.5811388300841898*Ghat[1]*dv10*dx10; 
  out[7] += -(1.224744871391589*Ghat[2]*dv10*dx10); 

  } 
  return 0.0;

} 
