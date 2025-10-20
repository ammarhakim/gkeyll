#include <gkyl_maxwell_kernels.h> 
GKYL_CU_DH double maxwell_vol_2x_ser_p1(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *q, double* GKYL_RESTRICT out) 
{ 
  const double c2 = meq->c*meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double *ex = &q[0]; 
  const double *ey = &q[4]; 
  const double *ez = &q[8]; 
  const double *bx = &q[12]; 
  const double *by = &q[16]; 
  const double *bz = &q[20]; 
  const double *ph = &q[24]; 
  const double *ps = &q[28]; 
 
  double *outEx = &out[0]; 
  double *outEy = &out[4]; 
  double *outEz = &out[8]; 
  double *outBx = &out[12]; 
  double *outBy = &out[16]; 
  double *outBz = &out[20]; 
  double *outPh = &out[24]; 
  double *outPs = &out[28]; 
 
  double dx0 = 2.0/dx[0]; 
  double dx1 = 2.0/dx[1]; 

  outEx[1] += 1.732050807568877*ph[0]*c2chi*dx0; 
  outEx[2] += -1.732050807568877*bz[0]*c2*dx1; 
  outEx[3] += 1.732050807568877*ph[2]*c2chi*dx0-1.732050807568877*bz[1]*c2*dx1; 

  outEy[1] += 1.732050807568877*bz[0]*c2*dx0; 
  outEy[2] += 1.732050807568877*ph[0]*c2chi*dx1; 
  outEy[3] += 1.732050807568877*ph[1]*c2chi*dx1+1.732050807568877*bz[2]*c2*dx0; 

  outEz[1] += -1.732050807568877*by[0]*c2*dx0; 
  outEz[2] += 1.732050807568877*bx[0]*c2*dx1; 
  outEz[3] += 1.732050807568877*bx[1]*c2*dx1-1.732050807568877*by[2]*c2*dx0; 

  outBx[1] += 1.732050807568877*ps[0]*dx0*gamma; 
  outBx[2] += 1.732050807568877*ez[0]*dx1; 
  outBx[3] += 1.732050807568877*ps[2]*dx0*gamma+1.732050807568877*ez[1]*dx1; 

  outBy[1] += -1.732050807568877*ez[0]*dx0; 
  outBy[2] += 1.732050807568877*ps[0]*dx1*gamma; 
  outBy[3] += 1.732050807568877*ps[1]*dx1*gamma-1.732050807568877*ez[2]*dx0; 

  outBz[1] += 1.732050807568877*ey[0]*dx0; 
  outBz[2] += -1.732050807568877*ex[0]*dx1; 
  outBz[3] += 1.732050807568877*ey[2]*dx0-1.732050807568877*ex[1]*dx1; 

  outPh[1] += 1.732050807568877*ex[0]*chi*dx0; 
  outPh[2] += 1.732050807568877*ey[0]*chi*dx1; 
  outPh[3] += 1.732050807568877*ey[1]*chi*dx1+1.732050807568877*ex[2]*chi*dx0; 

  outPs[1] += 1.732050807568877*bx[0]*c2gamma*dx0; 
  outPs[2] += 1.732050807568877*by[0]*c2gamma*dx1; 
  outPs[3] += 1.732050807568877*by[1]*c2gamma*dx1+1.732050807568877*bx[2]*c2gamma*dx0; 

  double cflFreq = 0.0; 
  cflFreq += meq->c/dx[0]; 
  cflFreq += meq->c/dx[1]; 
  return 3.0*cflFreq; 
} 
