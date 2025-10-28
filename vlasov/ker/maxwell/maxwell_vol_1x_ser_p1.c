#include <gkyl_maxwell_kernels.h> 
GKYL_CU_DH double maxwell_vol_1x_ser_p1(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *q, double* GKYL_RESTRICT out) 
{ 
  const double c2 = meq->c*meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double *ex = &q[0]; 
  const double *ey = &q[2]; 
  const double *ez = &q[4]; 
  const double *bx = &q[6]; 
  const double *by = &q[8]; 
  const double *bz = &q[10]; 
  const double *ph = &q[12]; 
  const double *ps = &q[14]; 
 
  double *outEx = &out[0]; 
  double *outEy = &out[2]; 
  double *outEz = &out[4]; 
  double *outBx = &out[6]; 
  double *outBy = &out[8]; 
  double *outBz = &out[10]; 
  double *outPh = &out[12]; 
  double *outPs = &out[14]; 
 
  double dx0 = 2.0/dx[0]; 

  outEx[1] += 1.732050807568877*ph[0]*c2chi*dx0; 

  outEy[1] += 1.732050807568877*bz[0]*c2*dx0; 

  outEz[1] += -1.732050807568877*by[0]*c2*dx0; 

  outBx[1] += 1.732050807568877*ps[0]*dx0*gamma; 

  outBy[1] += -1.732050807568877*ez[0]*dx0; 

  outBz[1] += 1.732050807568877*ey[0]*dx0; 

  outPh[1] += 1.732050807568877*ex[0]*chi*dx0; 

  outPs[1] += 1.732050807568877*bx[0]*c2gamma*dx0; 

  double cflFreq = 0.0; 
  cflFreq += meq->c/dx[0]; 
  return 3.0*cflFreq; 
} 
