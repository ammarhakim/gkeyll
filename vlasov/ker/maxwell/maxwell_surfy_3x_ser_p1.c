#include <gkyl_maxwell_kernels.h> 
GKYL_CU_DH double maxwell_surfy_3x_ser_p1(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out) 
{ 
  const double c = meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2 = c*c; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double dx1 = 2.0/dx[1]; 
  const double *exl = &ql[0]; 
  const double *eyl = &ql[8]; 
  const double *ezl = &ql[16]; 
  const double *bxl = &ql[24]; 
  const double *byl = &ql[32]; 
  const double *bzl = &ql[40]; 
  const double *phl = &ql[48]; 
  const double *psl = &ql[56]; 
 
  const double *exc = &qc[0]; 
  const double *eyc = &qc[8]; 
  const double *ezc = &qc[16]; 
  const double *bxc = &qc[24]; 
  const double *byc = &qc[32]; 
  const double *bzc = &qc[40]; 
  const double *phc = &qc[48]; 
  const double *psc = &qc[56]; 
 
  const double *exr = &qr[0]; 
  const double *eyr = &qr[8]; 
  const double *ezr = &qr[16]; 
  const double *bxr = &qr[24]; 
  const double *byr = &qr[32]; 
  const double *bzr = &qr[40]; 
  const double *phr = &qr[48]; 
  const double *psr = &qr[56]; 
 
  double *outEx = &out[0]; 
  double *outEy = &out[8]; 
  double *outEz = &out[16]; 
  double *outBx = &out[24]; 
  double *outBy = &out[32]; 
  double *outBz = &out[40]; 
  double *outPh = &out[48]; 
  double *outPs = &out[56]; 
 
  double incr_l[8]; 
 
  double incr_r[8]; 
 
  incr_l[0] = ((-0.4330127018922193*bzl[2])+0.4330127018922193*bzc[2]-0.25*(bzl[0]+bzc[0]))*c2; 
  incr_l[1] = ((-0.4330127018922193*bzl[4])+0.4330127018922193*bzc[4]-0.25*(bzl[1]+bzc[1]))*c2; 
  incr_l[2] = (0.75*bzl[2]-0.75*bzc[2]+0.4330127018922193*(bzl[0]+bzc[0]))*c2; 
  incr_l[3] = ((-0.4330127018922193*bzl[6])+0.4330127018922193*bzc[6]-0.25*(bzl[3]+bzc[3]))*c2; 
  incr_l[4] = (0.75*bzl[4]-0.75*bzc[4]+0.4330127018922193*(bzl[1]+bzc[1]))*c2; 
  incr_l[5] = ((-0.4330127018922193*bzl[7])+0.4330127018922193*bzc[7]-0.25*(bzl[5]+bzc[5]))*c2; 
  incr_l[6] = (0.75*bzl[6]-0.75*bzc[6]+0.4330127018922193*(bzl[3]+bzc[3]))*c2; 
  incr_l[7] = (0.75*bzl[7]-0.75*bzc[7]+0.4330127018922193*(bzl[5]+bzc[5]))*c2; 

  incr_r[0] = ((-0.4330127018922193*bzr[2])+0.4330127018922193*bzc[2]+0.25*(bzr[0]+bzc[0]))*c2; 
  incr_r[1] = ((-0.4330127018922193*bzr[4])+0.4330127018922193*bzc[4]+0.25*(bzr[1]+bzc[1]))*c2; 
  incr_r[2] = ((-0.75*bzr[2])+0.75*bzc[2]+0.4330127018922193*(bzr[0]+bzc[0]))*c2; 
  incr_r[3] = ((-0.4330127018922193*bzr[6])+0.4330127018922193*bzc[6]+0.25*(bzr[3]+bzc[3]))*c2; 
  incr_r[4] = ((-0.75*bzr[4])+0.75*bzc[4]+0.4330127018922193*(bzr[1]+bzc[1]))*c2; 
  incr_r[5] = ((-0.4330127018922193*bzr[7])+0.4330127018922193*bzc[7]+0.25*(bzr[5]+bzc[5]))*c2; 
  incr_r[6] = ((-0.75*bzr[6])+0.75*bzc[6]+0.4330127018922193*(bzr[3]+bzc[3]))*c2; 
  incr_r[7] = ((-0.75*bzr[7])+0.75*bzc[7]+0.4330127018922193*(bzr[5]+bzc[5]))*c2; 

  outEx[0] += (incr_r[0]+incr_l[0])*dx1; 
  outEx[1] += (incr_r[1]+incr_l[1])*dx1; 
  outEx[2] += (incr_r[2]+incr_l[2])*dx1; 
  outEx[3] += (incr_r[3]+incr_l[3])*dx1; 
  outEx[4] += (incr_r[4]+incr_l[4])*dx1; 
  outEx[5] += (incr_r[5]+incr_l[5])*dx1; 
  outEx[6] += (incr_r[6]+incr_l[6])*dx1; 
  outEx[7] += (incr_r[7]+incr_l[7])*dx1; 

  incr_l[0] = (0.4330127018922193*phl[2]-0.4330127018922193*phc[2]+0.25*(phl[0]+phc[0]))*c2chi; 
  incr_l[1] = (0.4330127018922193*phl[4]-0.4330127018922193*phc[4]+0.25*(phl[1]+phc[1]))*c2chi; 
  incr_l[2] = ((-0.75*phl[2])+0.75*phc[2]-0.4330127018922193*(phl[0]+phc[0]))*c2chi; 
  incr_l[3] = (0.4330127018922193*phl[6]-0.4330127018922193*phc[6]+0.25*(phl[3]+phc[3]))*c2chi; 
  incr_l[4] = ((-0.75*phl[4])+0.75*phc[4]-0.4330127018922193*(phl[1]+phc[1]))*c2chi; 
  incr_l[5] = (0.4330127018922193*phl[7]-0.4330127018922193*phc[7]+0.25*(phl[5]+phc[5]))*c2chi; 
  incr_l[6] = ((-0.75*phl[6])+0.75*phc[6]-0.4330127018922193*(phl[3]+phc[3]))*c2chi; 
  incr_l[7] = ((-0.75*phl[7])+0.75*phc[7]-0.4330127018922193*(phl[5]+phc[5]))*c2chi; 

  incr_r[0] = (0.4330127018922193*phr[2]-0.4330127018922193*phc[2]-0.25*(phr[0]+phc[0]))*c2chi; 
  incr_r[1] = (0.4330127018922193*phr[4]-0.4330127018922193*phc[4]-0.25*(phr[1]+phc[1]))*c2chi; 
  incr_r[2] = (0.75*phr[2]-0.75*phc[2]-0.4330127018922193*(phr[0]+phc[0]))*c2chi; 
  incr_r[3] = (0.4330127018922193*phr[6]-0.4330127018922193*phc[6]-0.25*(phr[3]+phc[3]))*c2chi; 
  incr_r[4] = (0.75*phr[4]-0.75*phc[4]-0.4330127018922193*(phr[1]+phc[1]))*c2chi; 
  incr_r[5] = (0.4330127018922193*phr[7]-0.4330127018922193*phc[7]-0.25*(phr[5]+phc[5]))*c2chi; 
  incr_r[6] = (0.75*phr[6]-0.75*phc[6]-0.4330127018922193*(phr[3]+phc[3]))*c2chi; 
  incr_r[7] = (0.75*phr[7]-0.75*phc[7]-0.4330127018922193*(phr[5]+phc[5]))*c2chi; 

  outEy[0] += (incr_r[0]+incr_l[0])*dx1; 
  outEy[1] += (incr_r[1]+incr_l[1])*dx1; 
  outEy[2] += (incr_r[2]+incr_l[2])*dx1; 
  outEy[3] += (incr_r[3]+incr_l[3])*dx1; 
  outEy[4] += (incr_r[4]+incr_l[4])*dx1; 
  outEy[5] += (incr_r[5]+incr_l[5])*dx1; 
  outEy[6] += (incr_r[6]+incr_l[6])*dx1; 
  outEy[7] += (incr_r[7]+incr_l[7])*dx1; 

  incr_l[0] = (0.4330127018922193*bxl[2]-0.4330127018922193*bxc[2]+0.25*(bxl[0]+bxc[0]))*c2; 
  incr_l[1] = (0.4330127018922193*bxl[4]-0.4330127018922193*bxc[4]+0.25*(bxl[1]+bxc[1]))*c2; 
  incr_l[2] = ((-0.75*bxl[2])+0.75*bxc[2]-0.4330127018922193*(bxl[0]+bxc[0]))*c2; 
  incr_l[3] = (0.4330127018922193*bxl[6]-0.4330127018922193*bxc[6]+0.25*(bxl[3]+bxc[3]))*c2; 
  incr_l[4] = ((-0.75*bxl[4])+0.75*bxc[4]-0.4330127018922193*(bxl[1]+bxc[1]))*c2; 
  incr_l[5] = (0.4330127018922193*bxl[7]-0.4330127018922193*bxc[7]+0.25*(bxl[5]+bxc[5]))*c2; 
  incr_l[6] = ((-0.75*bxl[6])+0.75*bxc[6]-0.4330127018922193*(bxl[3]+bxc[3]))*c2; 
  incr_l[7] = ((-0.75*bxl[7])+0.75*bxc[7]-0.4330127018922193*(bxl[5]+bxc[5]))*c2; 

  incr_r[0] = (0.4330127018922193*bxr[2]-0.4330127018922193*bxc[2]-0.25*(bxr[0]+bxc[0]))*c2; 
  incr_r[1] = (0.4330127018922193*bxr[4]-0.4330127018922193*bxc[4]-0.25*(bxr[1]+bxc[1]))*c2; 
  incr_r[2] = (0.75*bxr[2]-0.75*bxc[2]-0.4330127018922193*(bxr[0]+bxc[0]))*c2; 
  incr_r[3] = (0.4330127018922193*bxr[6]-0.4330127018922193*bxc[6]-0.25*(bxr[3]+bxc[3]))*c2; 
  incr_r[4] = (0.75*bxr[4]-0.75*bxc[4]-0.4330127018922193*(bxr[1]+bxc[1]))*c2; 
  incr_r[5] = (0.4330127018922193*bxr[7]-0.4330127018922193*bxc[7]-0.25*(bxr[5]+bxc[5]))*c2; 
  incr_r[6] = (0.75*bxr[6]-0.75*bxc[6]-0.4330127018922193*(bxr[3]+bxc[3]))*c2; 
  incr_r[7] = (0.75*bxr[7]-0.75*bxc[7]-0.4330127018922193*(bxr[5]+bxc[5]))*c2; 

  outEz[0] += (incr_r[0]+incr_l[0])*dx1; 
  outEz[1] += (incr_r[1]+incr_l[1])*dx1; 
  outEz[2] += (incr_r[2]+incr_l[2])*dx1; 
  outEz[3] += (incr_r[3]+incr_l[3])*dx1; 
  outEz[4] += (incr_r[4]+incr_l[4])*dx1; 
  outEz[5] += (incr_r[5]+incr_l[5])*dx1; 
  outEz[6] += (incr_r[6]+incr_l[6])*dx1; 
  outEz[7] += (incr_r[7]+incr_l[7])*dx1; 

  incr_l[0] = 0.4330127018922193*ezl[2]-0.4330127018922193*ezc[2]+0.25*(ezl[0]+ezc[0]); 
  incr_l[1] = 0.4330127018922193*ezl[4]-0.4330127018922193*ezc[4]+0.25*(ezl[1]+ezc[1]); 
  incr_l[2] = (-0.75*ezl[2])+0.75*ezc[2]-0.4330127018922193*(ezl[0]+ezc[0]); 
  incr_l[3] = 0.4330127018922193*ezl[6]-0.4330127018922193*ezc[6]+0.25*(ezl[3]+ezc[3]); 
  incr_l[4] = (-0.75*ezl[4])+0.75*ezc[4]-0.4330127018922193*(ezl[1]+ezc[1]); 
  incr_l[5] = 0.4330127018922193*ezl[7]-0.4330127018922193*ezc[7]+0.25*(ezl[5]+ezc[5]); 
  incr_l[6] = (-0.75*ezl[6])+0.75*ezc[6]-0.4330127018922193*(ezl[3]+ezc[3]); 
  incr_l[7] = (-0.75*ezl[7])+0.75*ezc[7]-0.4330127018922193*(ezl[5]+ezc[5]); 

  incr_r[0] = 0.4330127018922193*ezr[2]-0.4330127018922193*ezc[2]-0.25*(ezr[0]+ezc[0]); 
  incr_r[1] = 0.4330127018922193*ezr[4]-0.4330127018922193*ezc[4]-0.25*(ezr[1]+ezc[1]); 
  incr_r[2] = 0.75*ezr[2]-0.75*ezc[2]-0.4330127018922193*(ezr[0]+ezc[0]); 
  incr_r[3] = 0.4330127018922193*ezr[6]-0.4330127018922193*ezc[6]-0.25*(ezr[3]+ezc[3]); 
  incr_r[4] = 0.75*ezr[4]-0.75*ezc[4]-0.4330127018922193*(ezr[1]+ezc[1]); 
  incr_r[5] = 0.4330127018922193*ezr[7]-0.4330127018922193*ezc[7]-0.25*(ezr[5]+ezc[5]); 
  incr_r[6] = 0.75*ezr[6]-0.75*ezc[6]-0.4330127018922193*(ezr[3]+ezc[3]); 
  incr_r[7] = 0.75*ezr[7]-0.75*ezc[7]-0.4330127018922193*(ezr[5]+ezc[5]); 

  outBx[0] += (incr_r[0]+incr_l[0])*dx1; 
  outBx[1] += (incr_r[1]+incr_l[1])*dx1; 
  outBx[2] += (incr_r[2]+incr_l[2])*dx1; 
  outBx[3] += (incr_r[3]+incr_l[3])*dx1; 
  outBx[4] += (incr_r[4]+incr_l[4])*dx1; 
  outBx[5] += (incr_r[5]+incr_l[5])*dx1; 
  outBx[6] += (incr_r[6]+incr_l[6])*dx1; 
  outBx[7] += (incr_r[7]+incr_l[7])*dx1; 

  incr_l[0] = (0.4330127018922193*psl[2]-0.4330127018922193*psc[2]+0.25*(psl[0]+psc[0]))*gamma; 
  incr_l[1] = (0.4330127018922193*psl[4]-0.4330127018922193*psc[4]+0.25*(psl[1]+psc[1]))*gamma; 
  incr_l[2] = ((-0.75*psl[2])+0.75*psc[2]-0.4330127018922193*(psl[0]+psc[0]))*gamma; 
  incr_l[3] = (0.4330127018922193*psl[6]-0.4330127018922193*psc[6]+0.25*(psl[3]+psc[3]))*gamma; 
  incr_l[4] = ((-0.75*psl[4])+0.75*psc[4]-0.4330127018922193*(psl[1]+psc[1]))*gamma; 
  incr_l[5] = (0.4330127018922193*psl[7]-0.4330127018922193*psc[7]+0.25*(psl[5]+psc[5]))*gamma; 
  incr_l[6] = ((-0.75*psl[6])+0.75*psc[6]-0.4330127018922193*(psl[3]+psc[3]))*gamma; 
  incr_l[7] = ((-0.75*psl[7])+0.75*psc[7]-0.4330127018922193*(psl[5]+psc[5]))*gamma; 

  incr_r[0] = (0.4330127018922193*psr[2]-0.4330127018922193*psc[2]-0.25*(psr[0]+psc[0]))*gamma; 
  incr_r[1] = (0.4330127018922193*psr[4]-0.4330127018922193*psc[4]-0.25*(psr[1]+psc[1]))*gamma; 
  incr_r[2] = (0.75*psr[2]-0.75*psc[2]-0.4330127018922193*(psr[0]+psc[0]))*gamma; 
  incr_r[3] = (0.4330127018922193*psr[6]-0.4330127018922193*psc[6]-0.25*(psr[3]+psc[3]))*gamma; 
  incr_r[4] = (0.75*psr[4]-0.75*psc[4]-0.4330127018922193*(psr[1]+psc[1]))*gamma; 
  incr_r[5] = (0.4330127018922193*psr[7]-0.4330127018922193*psc[7]-0.25*(psr[5]+psc[5]))*gamma; 
  incr_r[6] = (0.75*psr[6]-0.75*psc[6]-0.4330127018922193*(psr[3]+psc[3]))*gamma; 
  incr_r[7] = (0.75*psr[7]-0.75*psc[7]-0.4330127018922193*(psr[5]+psc[5]))*gamma; 

  outBy[0] += (incr_r[0]+incr_l[0])*dx1; 
  outBy[1] += (incr_r[1]+incr_l[1])*dx1; 
  outBy[2] += (incr_r[2]+incr_l[2])*dx1; 
  outBy[3] += (incr_r[3]+incr_l[3])*dx1; 
  outBy[4] += (incr_r[4]+incr_l[4])*dx1; 
  outBy[5] += (incr_r[5]+incr_l[5])*dx1; 
  outBy[6] += (incr_r[6]+incr_l[6])*dx1; 
  outBy[7] += (incr_r[7]+incr_l[7])*dx1; 

  incr_l[0] = (-0.4330127018922193*exl[2])+0.4330127018922193*exc[2]-0.25*(exl[0]+exc[0]); 
  incr_l[1] = (-0.4330127018922193*exl[4])+0.4330127018922193*exc[4]-0.25*(exl[1]+exc[1]); 
  incr_l[2] = 0.75*exl[2]-0.75*exc[2]+0.4330127018922193*(exl[0]+exc[0]); 
  incr_l[3] = (-0.4330127018922193*exl[6])+0.4330127018922193*exc[6]-0.25*(exl[3]+exc[3]); 
  incr_l[4] = 0.75*exl[4]-0.75*exc[4]+0.4330127018922193*(exl[1]+exc[1]); 
  incr_l[5] = (-0.4330127018922193*exl[7])+0.4330127018922193*exc[7]-0.25*(exl[5]+exc[5]); 
  incr_l[6] = 0.75*exl[6]-0.75*exc[6]+0.4330127018922193*(exl[3]+exc[3]); 
  incr_l[7] = 0.75*exl[7]-0.75*exc[7]+0.4330127018922193*(exl[5]+exc[5]); 

  incr_r[0] = (-0.4330127018922193*exr[2])+0.4330127018922193*exc[2]+0.25*(exr[0]+exc[0]); 
  incr_r[1] = (-0.4330127018922193*exr[4])+0.4330127018922193*exc[4]+0.25*(exr[1]+exc[1]); 
  incr_r[2] = (-0.75*exr[2])+0.75*exc[2]+0.4330127018922193*(exr[0]+exc[0]); 
  incr_r[3] = (-0.4330127018922193*exr[6])+0.4330127018922193*exc[6]+0.25*(exr[3]+exc[3]); 
  incr_r[4] = (-0.75*exr[4])+0.75*exc[4]+0.4330127018922193*(exr[1]+exc[1]); 
  incr_r[5] = (-0.4330127018922193*exr[7])+0.4330127018922193*exc[7]+0.25*(exr[5]+exc[5]); 
  incr_r[6] = (-0.75*exr[6])+0.75*exc[6]+0.4330127018922193*(exr[3]+exc[3]); 
  incr_r[7] = (-0.75*exr[7])+0.75*exc[7]+0.4330127018922193*(exr[5]+exc[5]); 

  outBz[0] += (incr_r[0]+incr_l[0])*dx1; 
  outBz[1] += (incr_r[1]+incr_l[1])*dx1; 
  outBz[2] += (incr_r[2]+incr_l[2])*dx1; 
  outBz[3] += (incr_r[3]+incr_l[3])*dx1; 
  outBz[4] += (incr_r[4]+incr_l[4])*dx1; 
  outBz[5] += (incr_r[5]+incr_l[5])*dx1; 
  outBz[6] += (incr_r[6]+incr_l[6])*dx1; 
  outBz[7] += (incr_r[7]+incr_l[7])*dx1; 

  incr_l[0] = (0.4330127018922193*eyl[2]-0.4330127018922193*eyc[2]+0.25*(eyl[0]+eyc[0]))*chi; 
  incr_l[1] = (0.4330127018922193*eyl[4]-0.4330127018922193*eyc[4]+0.25*(eyl[1]+eyc[1]))*chi; 
  incr_l[2] = ((-0.75*eyl[2])+0.75*eyc[2]-0.4330127018922193*(eyl[0]+eyc[0]))*chi; 
  incr_l[3] = (0.4330127018922193*eyl[6]-0.4330127018922193*eyc[6]+0.25*(eyl[3]+eyc[3]))*chi; 
  incr_l[4] = ((-0.75*eyl[4])+0.75*eyc[4]-0.4330127018922193*(eyl[1]+eyc[1]))*chi; 
  incr_l[5] = (0.4330127018922193*eyl[7]-0.4330127018922193*eyc[7]+0.25*(eyl[5]+eyc[5]))*chi; 
  incr_l[6] = ((-0.75*eyl[6])+0.75*eyc[6]-0.4330127018922193*(eyl[3]+eyc[3]))*chi; 
  incr_l[7] = ((-0.75*eyl[7])+0.75*eyc[7]-0.4330127018922193*(eyl[5]+eyc[5]))*chi; 

  incr_r[0] = (0.4330127018922193*eyr[2]-0.4330127018922193*eyc[2]-0.25*(eyr[0]+eyc[0]))*chi; 
  incr_r[1] = (0.4330127018922193*eyr[4]-0.4330127018922193*eyc[4]-0.25*(eyr[1]+eyc[1]))*chi; 
  incr_r[2] = (0.75*eyr[2]-0.75*eyc[2]-0.4330127018922193*(eyr[0]+eyc[0]))*chi; 
  incr_r[3] = (0.4330127018922193*eyr[6]-0.4330127018922193*eyc[6]-0.25*(eyr[3]+eyc[3]))*chi; 
  incr_r[4] = (0.75*eyr[4]-0.75*eyc[4]-0.4330127018922193*(eyr[1]+eyc[1]))*chi; 
  incr_r[5] = (0.4330127018922193*eyr[7]-0.4330127018922193*eyc[7]-0.25*(eyr[5]+eyc[5]))*chi; 
  incr_r[6] = (0.75*eyr[6]-0.75*eyc[6]-0.4330127018922193*(eyr[3]+eyc[3]))*chi; 
  incr_r[7] = (0.75*eyr[7]-0.75*eyc[7]-0.4330127018922193*(eyr[5]+eyc[5]))*chi; 

  outPh[0] += (incr_r[0]+incr_l[0])*dx1; 
  outPh[1] += (incr_r[1]+incr_l[1])*dx1; 
  outPh[2] += (incr_r[2]+incr_l[2])*dx1; 
  outPh[3] += (incr_r[3]+incr_l[3])*dx1; 
  outPh[4] += (incr_r[4]+incr_l[4])*dx1; 
  outPh[5] += (incr_r[5]+incr_l[5])*dx1; 
  outPh[6] += (incr_r[6]+incr_l[6])*dx1; 
  outPh[7] += (incr_r[7]+incr_l[7])*dx1; 

  incr_l[0] = (0.4330127018922193*byl[2]-0.4330127018922193*byc[2]+0.25*(byl[0]+byc[0]))*c2gamma; 
  incr_l[1] = (0.4330127018922193*byl[4]-0.4330127018922193*byc[4]+0.25*(byl[1]+byc[1]))*c2gamma; 
  incr_l[2] = ((-0.75*byl[2])+0.75*byc[2]-0.4330127018922193*(byl[0]+byc[0]))*c2gamma; 
  incr_l[3] = (0.4330127018922193*byl[6]-0.4330127018922193*byc[6]+0.25*(byl[3]+byc[3]))*c2gamma; 
  incr_l[4] = ((-0.75*byl[4])+0.75*byc[4]-0.4330127018922193*(byl[1]+byc[1]))*c2gamma; 
  incr_l[5] = (0.4330127018922193*byl[7]-0.4330127018922193*byc[7]+0.25*(byl[5]+byc[5]))*c2gamma; 
  incr_l[6] = ((-0.75*byl[6])+0.75*byc[6]-0.4330127018922193*(byl[3]+byc[3]))*c2gamma; 
  incr_l[7] = ((-0.75*byl[7])+0.75*byc[7]-0.4330127018922193*(byl[5]+byc[5]))*c2gamma; 

  incr_r[0] = (0.4330127018922193*byr[2]-0.4330127018922193*byc[2]-0.25*(byr[0]+byc[0]))*c2gamma; 
  incr_r[1] = (0.4330127018922193*byr[4]-0.4330127018922193*byc[4]-0.25*(byr[1]+byc[1]))*c2gamma; 
  incr_r[2] = (0.75*byr[2]-0.75*byc[2]-0.4330127018922193*(byr[0]+byc[0]))*c2gamma; 
  incr_r[3] = (0.4330127018922193*byr[6]-0.4330127018922193*byc[6]-0.25*(byr[3]+byc[3]))*c2gamma; 
  incr_r[4] = (0.75*byr[4]-0.75*byc[4]-0.4330127018922193*(byr[1]+byc[1]))*c2gamma; 
  incr_r[5] = (0.4330127018922193*byr[7]-0.4330127018922193*byc[7]-0.25*(byr[5]+byc[5]))*c2gamma; 
  incr_r[6] = (0.75*byr[6]-0.75*byc[6]-0.4330127018922193*(byr[3]+byc[3]))*c2gamma; 
  incr_r[7] = (0.75*byr[7]-0.75*byc[7]-0.4330127018922193*(byr[5]+byc[5]))*c2gamma; 

  outPs[0] += (incr_r[0]+incr_l[0])*dx1; 
  outPs[1] += (incr_r[1]+incr_l[1])*dx1; 
  outPs[2] += (incr_r[2]+incr_l[2])*dx1; 
  outPs[3] += (incr_r[3]+incr_l[3])*dx1; 
  outPs[4] += (incr_r[4]+incr_l[4])*dx1; 
  outPs[5] += (incr_r[5]+incr_l[5])*dx1; 
  outPs[6] += (incr_r[6]+incr_l[6])*dx1; 
  outPs[7] += (incr_r[7]+incr_l[7])*dx1; 

  return 0.;

} 
