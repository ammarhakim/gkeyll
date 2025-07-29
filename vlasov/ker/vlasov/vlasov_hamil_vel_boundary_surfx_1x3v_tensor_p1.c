#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_hamil_vel_boundary_surfx_1x3v_tensor_p1(const double *w, const double *dxv, 
  const double *hamil, 
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = 2.0/dxv[0]; 
  const double dv10 = 2.0/dxv[1]; 
  const double wv = w[1]; 

  double Ghat[8]; 

  if (edge == -1) { 

  if (wv>0) { 

  Ghat[0] = 0.75*hamil[7]*fskin[13]+0.4330127018922193*hamil[7]*fskin[10]+0.75*hamil[5]*fskin[8]+0.75*hamil[4]*fskin[6]+0.4330127018922193*fskin[4]*hamil[5]+0.4330127018922193*fskin[3]*hamil[4]+0.75*fskin[1]*hamil[1]+0.4330127018922193*fskin[0]*hamil[1]; 
  Ghat[1] = 0.75*hamil[7]*fskin[15]+0.4330127018922193*hamil[7]*fskin[14]+0.75*hamil[5]*fskin[12]+0.75*hamil[4]*fskin[11]+0.4330127018922193*hamil[5]*fskin[9]+0.4330127018922193*hamil[4]*fskin[7]+0.75*hamil[1]*fskin[5]+0.4330127018922193*hamil[1]*fskin[2]; 
  Ghat[2] = 0.75*hamil[5]*fskin[13]+0.4330127018922193*hamil[5]*fskin[10]+0.75*hamil[7]*fskin[8]+0.4330127018922193*fskin[4]*hamil[7]+0.75*hamil[1]*fskin[6]+0.75*fskin[1]*hamil[4]+0.4330127018922193*fskin[0]*hamil[4]+0.4330127018922193*hamil[1]*fskin[3]; 
  Ghat[3] = 0.75*hamil[4]*fskin[13]+0.4330127018922193*hamil[4]*fskin[10]+0.75*hamil[1]*fskin[8]+0.75*fskin[6]*hamil[7]+0.4330127018922193*fskin[3]*hamil[7]+0.75*fskin[1]*hamil[5]+0.4330127018922193*fskin[0]*hamil[5]+0.4330127018922193*hamil[1]*fskin[4]; 
  Ghat[4] = 0.75*hamil[5]*fskin[15]+0.4330127018922193*hamil[5]*fskin[14]+0.75*hamil[7]*fskin[12]+0.75*hamil[1]*fskin[11]+0.4330127018922193*hamil[7]*fskin[9]+0.4330127018922193*hamil[1]*fskin[7]+0.75*hamil[4]*fskin[5]+0.4330127018922193*fskin[2]*hamil[4]; 
  Ghat[5] = 0.75*hamil[4]*fskin[15]+0.4330127018922193*hamil[4]*fskin[14]+0.75*hamil[1]*fskin[12]+0.75*hamil[7]*fskin[11]+0.4330127018922193*hamil[1]*fskin[9]+0.4330127018922193*fskin[7]*hamil[7]+0.75*fskin[5]*hamil[5]+0.4330127018922193*fskin[2]*hamil[5]; 
  Ghat[6] = 0.75*hamil[1]*fskin[13]+0.4330127018922193*hamil[1]*fskin[10]+0.75*hamil[4]*fskin[8]+0.75*fskin[1]*hamil[7]+0.4330127018922193*fskin[0]*hamil[7]+0.75*hamil[5]*fskin[6]+0.4330127018922193*fskin[3]*hamil[5]+0.4330127018922193*fskin[4]*hamil[4]; 
  Ghat[7] = 0.75*hamil[1]*fskin[15]+0.4330127018922193*hamil[1]*fskin[14]+0.75*hamil[4]*fskin[12]+0.75*hamil[5]*fskin[11]+0.4330127018922193*hamil[4]*fskin[9]+0.75*fskin[5]*hamil[7]+0.4330127018922193*fskin[2]*hamil[7]+0.4330127018922193*hamil[5]*fskin[7]; 

  } else { 

  Ghat[0] = -(0.75*hamil[7]*fedge[13])+0.4330127018922193*hamil[7]*fedge[10]-0.75*hamil[5]*fedge[8]-0.75*hamil[4]*fedge[6]+0.4330127018922193*fedge[4]*hamil[5]+0.4330127018922193*fedge[3]*hamil[4]-0.75*fedge[1]*hamil[1]+0.4330127018922193*fedge[0]*hamil[1]; 
  Ghat[1] = -(0.75*hamil[7]*fedge[15])+0.4330127018922193*hamil[7]*fedge[14]-0.75*hamil[5]*fedge[12]-0.75*hamil[4]*fedge[11]+0.4330127018922193*hamil[5]*fedge[9]+0.4330127018922193*hamil[4]*fedge[7]-0.75*hamil[1]*fedge[5]+0.4330127018922193*hamil[1]*fedge[2]; 
  Ghat[2] = -(0.75*hamil[5]*fedge[13])+0.4330127018922193*hamil[5]*fedge[10]-0.75*hamil[7]*fedge[8]+0.4330127018922193*fedge[4]*hamil[7]-0.75*hamil[1]*fedge[6]-0.75*fedge[1]*hamil[4]+0.4330127018922193*fedge[0]*hamil[4]+0.4330127018922193*hamil[1]*fedge[3]; 
  Ghat[3] = -(0.75*hamil[4]*fedge[13])+0.4330127018922193*hamil[4]*fedge[10]-0.75*hamil[1]*fedge[8]-0.75*fedge[6]*hamil[7]+0.4330127018922193*fedge[3]*hamil[7]-0.75*fedge[1]*hamil[5]+0.4330127018922193*fedge[0]*hamil[5]+0.4330127018922193*hamil[1]*fedge[4]; 
  Ghat[4] = -(0.75*hamil[5]*fedge[15])+0.4330127018922193*hamil[5]*fedge[14]-0.75*hamil[7]*fedge[12]-0.75*hamil[1]*fedge[11]+0.4330127018922193*hamil[7]*fedge[9]+0.4330127018922193*hamil[1]*fedge[7]-0.75*hamil[4]*fedge[5]+0.4330127018922193*fedge[2]*hamil[4]; 
  Ghat[5] = -(0.75*hamil[4]*fedge[15])+0.4330127018922193*hamil[4]*fedge[14]-0.75*hamil[1]*fedge[12]-0.75*hamil[7]*fedge[11]+0.4330127018922193*hamil[1]*fedge[9]+0.4330127018922193*fedge[7]*hamil[7]-0.75*fedge[5]*hamil[5]+0.4330127018922193*fedge[2]*hamil[5]; 
  Ghat[6] = -(0.75*hamil[1]*fedge[13])+0.4330127018922193*hamil[1]*fedge[10]-0.75*hamil[4]*fedge[8]-0.75*fedge[1]*hamil[7]+0.4330127018922193*fedge[0]*hamil[7]-0.75*hamil[5]*fedge[6]+0.4330127018922193*fedge[3]*hamil[5]+0.4330127018922193*fedge[4]*hamil[4]; 
  Ghat[7] = -(0.75*hamil[1]*fedge[15])+0.4330127018922193*hamil[1]*fedge[14]-0.75*hamil[4]*fedge[12]-0.75*hamil[5]*fedge[11]+0.4330127018922193*hamil[4]*fedge[9]-0.75*fedge[5]*hamil[7]+0.4330127018922193*fedge[2]*hamil[7]+0.4330127018922193*hamil[5]*fedge[7]; 

  } 

  out[0] += -(0.7071067811865475*Ghat[0]*dv10*dx10); 
  out[1] += -(1.224744871391589*Ghat[0]*dv10*dx10); 
  out[2] += -(0.7071067811865475*Ghat[1]*dv10*dx10); 
  out[3] += -(0.7071067811865475*Ghat[2]*dv10*dx10); 
  out[4] += -(0.7071067811865475*Ghat[3]*dv10*dx10); 
  out[5] += -(1.224744871391589*Ghat[1]*dv10*dx10); 
  out[6] += -(1.224744871391589*Ghat[2]*dv10*dx10); 
  out[7] += -(0.7071067811865475*Ghat[4]*dv10*dx10); 
  out[8] += -(1.224744871391589*Ghat[3]*dv10*dx10); 
  out[9] += -(0.7071067811865475*Ghat[5]*dv10*dx10); 
  out[10] += -(0.7071067811865475*Ghat[6]*dv10*dx10); 
  out[11] += -(1.224744871391589*Ghat[4]*dv10*dx10); 
  out[12] += -(1.224744871391589*Ghat[5]*dv10*dx10); 
  out[13] += -(1.224744871391589*Ghat[6]*dv10*dx10); 
  out[14] += -(0.7071067811865475*Ghat[7]*dv10*dx10); 
  out[15] += -(1.224744871391589*Ghat[7]*dv10*dx10); 

  } else { 

  if (wv>0) { 

  Ghat[0] = 0.75*hamil[7]*fedge[13]+0.4330127018922193*hamil[7]*fedge[10]+0.75*hamil[5]*fedge[8]+0.75*hamil[4]*fedge[6]+0.4330127018922193*fedge[4]*hamil[5]+0.4330127018922193*fedge[3]*hamil[4]+0.75*fedge[1]*hamil[1]+0.4330127018922193*fedge[0]*hamil[1]; 
  Ghat[1] = 0.75*hamil[7]*fedge[15]+0.4330127018922193*hamil[7]*fedge[14]+0.75*hamil[5]*fedge[12]+0.75*hamil[4]*fedge[11]+0.4330127018922193*hamil[5]*fedge[9]+0.4330127018922193*hamil[4]*fedge[7]+0.75*hamil[1]*fedge[5]+0.4330127018922193*hamil[1]*fedge[2]; 
  Ghat[2] = 0.75*hamil[5]*fedge[13]+0.4330127018922193*hamil[5]*fedge[10]+0.75*hamil[7]*fedge[8]+0.4330127018922193*fedge[4]*hamil[7]+0.75*hamil[1]*fedge[6]+0.75*fedge[1]*hamil[4]+0.4330127018922193*fedge[0]*hamil[4]+0.4330127018922193*hamil[1]*fedge[3]; 
  Ghat[3] = 0.75*hamil[4]*fedge[13]+0.4330127018922193*hamil[4]*fedge[10]+0.75*hamil[1]*fedge[8]+0.75*fedge[6]*hamil[7]+0.4330127018922193*fedge[3]*hamil[7]+0.75*fedge[1]*hamil[5]+0.4330127018922193*fedge[0]*hamil[5]+0.4330127018922193*hamil[1]*fedge[4]; 
  Ghat[4] = 0.75*hamil[5]*fedge[15]+0.4330127018922193*hamil[5]*fedge[14]+0.75*hamil[7]*fedge[12]+0.75*hamil[1]*fedge[11]+0.4330127018922193*hamil[7]*fedge[9]+0.4330127018922193*hamil[1]*fedge[7]+0.75*hamil[4]*fedge[5]+0.4330127018922193*fedge[2]*hamil[4]; 
  Ghat[5] = 0.75*hamil[4]*fedge[15]+0.4330127018922193*hamil[4]*fedge[14]+0.75*hamil[1]*fedge[12]+0.75*hamil[7]*fedge[11]+0.4330127018922193*hamil[1]*fedge[9]+0.4330127018922193*fedge[7]*hamil[7]+0.75*fedge[5]*hamil[5]+0.4330127018922193*fedge[2]*hamil[5]; 
  Ghat[6] = 0.75*hamil[1]*fedge[13]+0.4330127018922193*hamil[1]*fedge[10]+0.75*hamil[4]*fedge[8]+0.75*fedge[1]*hamil[7]+0.4330127018922193*fedge[0]*hamil[7]+0.75*hamil[5]*fedge[6]+0.4330127018922193*fedge[3]*hamil[5]+0.4330127018922193*fedge[4]*hamil[4]; 
  Ghat[7] = 0.75*hamil[1]*fedge[15]+0.4330127018922193*hamil[1]*fedge[14]+0.75*hamil[4]*fedge[12]+0.75*hamil[5]*fedge[11]+0.4330127018922193*hamil[4]*fedge[9]+0.75*fedge[5]*hamil[7]+0.4330127018922193*fedge[2]*hamil[7]+0.4330127018922193*hamil[5]*fedge[7]; 

  } else { 

  Ghat[0] = -(0.75*hamil[7]*fskin[13])+0.4330127018922193*hamil[7]*fskin[10]-0.75*hamil[5]*fskin[8]-0.75*hamil[4]*fskin[6]+0.4330127018922193*fskin[4]*hamil[5]+0.4330127018922193*fskin[3]*hamil[4]-0.75*fskin[1]*hamil[1]+0.4330127018922193*fskin[0]*hamil[1]; 
  Ghat[1] = -(0.75*hamil[7]*fskin[15])+0.4330127018922193*hamil[7]*fskin[14]-0.75*hamil[5]*fskin[12]-0.75*hamil[4]*fskin[11]+0.4330127018922193*hamil[5]*fskin[9]+0.4330127018922193*hamil[4]*fskin[7]-0.75*hamil[1]*fskin[5]+0.4330127018922193*hamil[1]*fskin[2]; 
  Ghat[2] = -(0.75*hamil[5]*fskin[13])+0.4330127018922193*hamil[5]*fskin[10]-0.75*hamil[7]*fskin[8]+0.4330127018922193*fskin[4]*hamil[7]-0.75*hamil[1]*fskin[6]-0.75*fskin[1]*hamil[4]+0.4330127018922193*fskin[0]*hamil[4]+0.4330127018922193*hamil[1]*fskin[3]; 
  Ghat[3] = -(0.75*hamil[4]*fskin[13])+0.4330127018922193*hamil[4]*fskin[10]-0.75*hamil[1]*fskin[8]-0.75*fskin[6]*hamil[7]+0.4330127018922193*fskin[3]*hamil[7]-0.75*fskin[1]*hamil[5]+0.4330127018922193*fskin[0]*hamil[5]+0.4330127018922193*hamil[1]*fskin[4]; 
  Ghat[4] = -(0.75*hamil[5]*fskin[15])+0.4330127018922193*hamil[5]*fskin[14]-0.75*hamil[7]*fskin[12]-0.75*hamil[1]*fskin[11]+0.4330127018922193*hamil[7]*fskin[9]+0.4330127018922193*hamil[1]*fskin[7]-0.75*hamil[4]*fskin[5]+0.4330127018922193*fskin[2]*hamil[4]; 
  Ghat[5] = -(0.75*hamil[4]*fskin[15])+0.4330127018922193*hamil[4]*fskin[14]-0.75*hamil[1]*fskin[12]-0.75*hamil[7]*fskin[11]+0.4330127018922193*hamil[1]*fskin[9]+0.4330127018922193*fskin[7]*hamil[7]-0.75*fskin[5]*hamil[5]+0.4330127018922193*fskin[2]*hamil[5]; 
  Ghat[6] = -(0.75*hamil[1]*fskin[13])+0.4330127018922193*hamil[1]*fskin[10]-0.75*hamil[4]*fskin[8]-0.75*fskin[1]*hamil[7]+0.4330127018922193*fskin[0]*hamil[7]-0.75*hamil[5]*fskin[6]+0.4330127018922193*fskin[3]*hamil[5]+0.4330127018922193*fskin[4]*hamil[4]; 
  Ghat[7] = -(0.75*hamil[1]*fskin[15])+0.4330127018922193*hamil[1]*fskin[14]-0.75*hamil[4]*fskin[12]-0.75*hamil[5]*fskin[11]+0.4330127018922193*hamil[4]*fskin[9]-0.75*fskin[5]*hamil[7]+0.4330127018922193*fskin[2]*hamil[7]+0.4330127018922193*hamil[5]*fskin[7]; 

  } 

  out[0] += 0.7071067811865475*Ghat[0]*dv10*dx10; 
  out[1] += -(1.224744871391589*Ghat[0]*dv10*dx10); 
  out[2] += 0.7071067811865475*Ghat[1]*dv10*dx10; 
  out[3] += 0.7071067811865475*Ghat[2]*dv10*dx10; 
  out[4] += 0.7071067811865475*Ghat[3]*dv10*dx10; 
  out[5] += -(1.224744871391589*Ghat[1]*dv10*dx10); 
  out[6] += -(1.224744871391589*Ghat[2]*dv10*dx10); 
  out[7] += 0.7071067811865475*Ghat[4]*dv10*dx10; 
  out[8] += -(1.224744871391589*Ghat[3]*dv10*dx10); 
  out[9] += 0.7071067811865475*Ghat[5]*dv10*dx10; 
  out[10] += 0.7071067811865475*Ghat[6]*dv10*dx10; 
  out[11] += -(1.224744871391589*Ghat[4]*dv10*dx10); 
  out[12] += -(1.224744871391589*Ghat[5]*dv10*dx10); 
  out[13] += -(1.224744871391589*Ghat[6]*dv10*dx10); 
  out[14] += 0.7071067811865475*Ghat[7]*dv10*dx10; 
  out[15] += -(1.224744871391589*Ghat[7]*dv10*dx10); 

  } 
  return 0.0;

} 
