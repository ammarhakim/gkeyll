#include <gkyl_maxwell_kernels.h> 
GKYL_CU_DH double maxwell_surfx_1x_ser_p1(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double dx1 = 2.0/dx[0]; 
  const double *exl = &ql[0]; 
  const double *eyl = &ql[2]; 
  const double *ezl = &ql[4]; 
  const double *bxl = &ql[6]; 
  const double *byl = &ql[8]; 
  const double *bzl = &ql[10]; 
  const double *phl = &ql[12]; 
  const double *psl = &ql[14]; 
 
  const double *exc = &qc[0]; 
  const double *eyc = &qc[2]; 
  const double *ezc = &qc[4]; 
  const double *bxc = &qc[6]; 
  const double *byc = &qc[8]; 
  const double *bzc = &qc[10]; 
  const double *phc = &qc[12]; 
  const double *psc = &qc[14]; 
 
  const double *exr = &qr[0]; 
  const double *eyr = &qr[2]; 
  const double *ezr = &qr[4]; 
  const double *bxr = &qr[6]; 
  const double *byr = &qr[8]; 
  const double *bzr = &qr[10]; 
  const double *phr = &qr[12]; 
  const double *psr = &qr[14]; 
 
  double *outEx = &out[0]; 
  double *outEy = &out[2]; 
  double *outEz = &out[4]; 
  double *outBx = &out[6]; 
  double *outBy = &out[8]; 
  double *outBz = &out[10]; 
  double *outPh = &out[12]; 
  double *outPs = &out[14]; 
 
  double incr_l[2]; 
 
  double incr_r[2]; 
 
  incr_l[0] = (0.4330127018922193*phl[1]-0.4330127018922193*phc[1]+0.25*(phl[0]+phc[0]))*c2chi; 
  incr_l[1] = ((-0.75*phl[1])+0.75*phc[1]-0.4330127018922193*(phl[0]+phc[0]))*c2chi; 

  incr_r[0] = (0.4330127018922193*phr[1]-0.4330127018922193*phc[1]-0.25*(phr[0]+phc[0]))*c2chi; 
  incr_r[1] = (0.75*phr[1]-0.75*phc[1]-0.4330127018922193*(phr[0]+phc[0]))*c2chi; 

  outEx[0] += (incr_r[0]+incr_l[0])*dx1; 
  outEx[1] += (incr_r[1]+incr_l[1])*dx1; 

  incr_l[0] = (0.4330127018922193*bzl[1]-0.4330127018922193*bzc[1]+0.25*(bzl[0]+bzc[0]))*c2; 
  incr_l[1] = ((-0.75*bzl[1])+0.75*bzc[1]-0.4330127018922193*(bzl[0]+bzc[0]))*c2; 

  incr_r[0] = (0.4330127018922193*bzr[1]-0.4330127018922193*bzc[1]-0.25*(bzr[0]+bzc[0]))*c2; 
  incr_r[1] = (0.75*bzr[1]-0.75*bzc[1]-0.4330127018922193*(bzr[0]+bzc[0]))*c2; 

  outEy[0] += (incr_r[0]+incr_l[0])*dx1; 
  outEy[1] += (incr_r[1]+incr_l[1])*dx1; 

  incr_l[0] = ((-0.4330127018922193*byl[1])+0.4330127018922193*byc[1]-0.25*(byl[0]+byc[0]))*c2; 
  incr_l[1] = (0.75*byl[1]-0.75*byc[1]+0.4330127018922193*(byl[0]+byc[0]))*c2; 

  incr_r[0] = ((-0.4330127018922193*byr[1])+0.4330127018922193*byc[1]+0.25*(byr[0]+byc[0]))*c2; 
  incr_r[1] = ((-0.75*byr[1])+0.75*byc[1]+0.4330127018922193*(byr[0]+byc[0]))*c2; 

  outEz[0] += (incr_r[0]+incr_l[0])*dx1; 
  outEz[1] += (incr_r[1]+incr_l[1])*dx1; 

  incr_l[0] = (0.4330127018922193*psl[1]-0.4330127018922193*psc[1]+0.25*(psl[0]+psc[0]))*gamma; 
  incr_l[1] = ((-0.75*psl[1])+0.75*psc[1]-0.4330127018922193*(psl[0]+psc[0]))*gamma; 

  incr_r[0] = (0.4330127018922193*psr[1]-0.4330127018922193*psc[1]-0.25*(psr[0]+psc[0]))*gamma; 
  incr_r[1] = (0.75*psr[1]-0.75*psc[1]-0.4330127018922193*(psr[0]+psc[0]))*gamma; 

  outBx[0] += (incr_r[0]+incr_l[0])*dx1; 
  outBx[1] += (incr_r[1]+incr_l[1])*dx1; 

  incr_l[0] = (-0.4330127018922193*ezl[1])+0.4330127018922193*ezc[1]-0.25*(ezl[0]+ezc[0]); 
  incr_l[1] = 0.75*ezl[1]-0.75*ezc[1]+0.4330127018922193*(ezl[0]+ezc[0]); 

  incr_r[0] = (-0.4330127018922193*ezr[1])+0.4330127018922193*ezc[1]+0.25*(ezr[0]+ezc[0]); 
  incr_r[1] = (-0.75*ezr[1])+0.75*ezc[1]+0.4330127018922193*(ezr[0]+ezc[0]); 

  outBy[0] += (incr_r[0]+incr_l[0])*dx1; 
  outBy[1] += (incr_r[1]+incr_l[1])*dx1; 

  incr_l[0] = 0.4330127018922193*eyl[1]-0.4330127018922193*eyc[1]+0.25*(eyl[0]+eyc[0]); 
  incr_l[1] = (-0.75*eyl[1])+0.75*eyc[1]-0.4330127018922193*(eyl[0]+eyc[0]); 

  incr_r[0] = 0.4330127018922193*eyr[1]-0.4330127018922193*eyc[1]-0.25*(eyr[0]+eyc[0]); 
  incr_r[1] = 0.75*eyr[1]-0.75*eyc[1]-0.4330127018922193*(eyr[0]+eyc[0]); 

  outBz[0] += (incr_r[0]+incr_l[0])*dx1; 
  outBz[1] += (incr_r[1]+incr_l[1])*dx1; 

  incr_l[0] = (0.4330127018922193*exl[1]-0.4330127018922193*exc[1]+0.25*(exl[0]+exc[0]))*chi; 
  incr_l[1] = ((-0.75*exl[1])+0.75*exc[1]-0.4330127018922193*(exl[0]+exc[0]))*chi; 

  incr_r[0] = (0.4330127018922193*exr[1]-0.4330127018922193*exc[1]-0.25*(exr[0]+exc[0]))*chi; 
  incr_r[1] = (0.75*exr[1]-0.75*exc[1]-0.4330127018922193*(exr[0]+exc[0]))*chi; 

  outPh[0] += (incr_r[0]+incr_l[0])*dx1; 
  outPh[1] += (incr_r[1]+incr_l[1])*dx1; 

  incr_l[0] = (0.4330127018922193*bxl[1]-0.4330127018922193*bxc[1]+0.25*(bxl[0]+bxc[0]))*c2gamma; 
  incr_l[1] = ((-0.75*bxl[1])+0.75*bxc[1]-0.4330127018922193*(bxl[0]+bxc[0]))*c2gamma; 

  incr_r[0] = (0.4330127018922193*bxr[1]-0.4330127018922193*bxc[1]-0.25*(bxr[0]+bxc[0]))*c2gamma; 
  incr_r[1] = (0.75*bxr[1]-0.75*bxc[1]-0.4330127018922193*(bxr[0]+bxc[0]))*c2gamma; 

  outPs[0] += (incr_r[0]+incr_l[0])*dx1; 
  outPs[1] += (incr_r[1]+incr_l[1])*dx1; 

  return 0.;

} 
