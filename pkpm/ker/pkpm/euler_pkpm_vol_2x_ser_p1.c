#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH double euler_pkpm_vol_2x_ser_p1(const double *w, const double *dxv, const double *prim, const double *p_ij, const double *euler_pkpm, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:    Cell-center coordinates.
  // dxv[NDIM]:  Cell spacing.
  // prim:       Input primitive variables [ux, uy, uz, 3*Txx/m, 3*Tyy/m, 3*Tzz/m, 1/rho div(p_par b), T_perp/m, m/T_perp].
  // p_ij:       Input pressure tensor.
  // euler_pkpm: [rho ux, rho uy, rho uz], Fluid input state vector.
  // out:        Incremented output.

  double dx10 = 2./dxv[0]; 
  double dx11 = 2./dxv[1]; 

  const double *rhoux = &euler_pkpm[0]; 
  const double *rhouy = &euler_pkpm[4]; 
  const double *rhouz = &euler_pkpm[8]; 

  const double *ux = &prim[0]; 
  const double *uy = &prim[4]; 
  const double *uz = &prim[8]; 

  const double *Pxx = &p_ij[0]; 
  const double *Pxy = &p_ij[4]; 
  const double *Pxz = &p_ij[8]; 
  const double *Pyy = &p_ij[12]; 
  const double *Pyz = &p_ij[16]; 
  const double *Pzz = &p_ij[20]; 

  double *outrhoux = &out[0]; 
  double *outrhouy = &out[4]; 
  double *outrhouz = &out[8]; 

  double cflFreq_mid = 0.0; 
  cflFreq_mid += 0.5*3.0*dx10*(fabs(0.5*ux[0])); 
  cflFreq_mid += 0.5*3.0*dx11*(fabs(0.5*uy[0])); 

  outrhoux[1] += 0.8660254037844386*rhoux[3]*ux[3]*dx10+0.8660254037844386*rhoux[2]*ux[2]*dx10+0.8660254037844386*rhoux[1]*ux[1]*dx10+0.8660254037844386*rhoux[0]*ux[0]*dx10+1.732050807568877*Pxx[0]*dx10; 
  outrhoux[2] += 0.8660254037844386*rhoux[3]*uy[3]*dx11+0.8660254037844386*rhoux[2]*uy[2]*dx11+0.8660254037844386*rhoux[1]*uy[1]*dx11+0.8660254037844386*rhoux[0]*uy[0]*dx11+1.732050807568877*Pxy[0]*dx11; 
  outrhoux[3] += 0.8660254037844386*rhoux[2]*uy[3]*dx11+0.8660254037844386*uy[2]*rhoux[3]*dx11+0.8660254037844386*rhoux[0]*uy[1]*dx11+0.8660254037844386*uy[0]*rhoux[1]*dx11+1.732050807568877*Pxy[1]*dx11+0.8660254037844386*rhoux[1]*ux[3]*dx10+0.8660254037844386*ux[1]*rhoux[3]*dx10+0.8660254037844386*rhoux[0]*ux[2]*dx10+0.8660254037844386*ux[0]*rhoux[2]*dx10+1.732050807568877*Pxx[2]*dx10; 

  outrhouy[1] += 0.8660254037844386*rhouy[3]*ux[3]*dx10+0.8660254037844386*rhouy[2]*ux[2]*dx10+0.8660254037844386*rhouy[1]*ux[1]*dx10+0.8660254037844386*rhouy[0]*ux[0]*dx10+1.732050807568877*Pxy[0]*dx10; 
  outrhouy[2] += 0.8660254037844386*rhouy[3]*uy[3]*dx11+0.8660254037844386*rhouy[2]*uy[2]*dx11+0.8660254037844386*rhouy[1]*uy[1]*dx11+0.8660254037844386*rhouy[0]*uy[0]*dx11+1.732050807568877*Pyy[0]*dx11; 
  outrhouy[3] += 0.8660254037844386*rhouy[2]*uy[3]*dx11+0.8660254037844386*uy[2]*rhouy[3]*dx11+0.8660254037844386*rhouy[0]*uy[1]*dx11+0.8660254037844386*uy[0]*rhouy[1]*dx11+1.732050807568877*Pyy[1]*dx11+0.8660254037844386*rhouy[1]*ux[3]*dx10+0.8660254037844386*ux[1]*rhouy[3]*dx10+0.8660254037844386*rhouy[0]*ux[2]*dx10+0.8660254037844386*ux[0]*rhouy[2]*dx10+1.732050807568877*Pxy[2]*dx10; 

  outrhouz[1] += 0.8660254037844386*rhouz[3]*ux[3]*dx10+0.8660254037844386*rhouz[2]*ux[2]*dx10+0.8660254037844386*rhouz[1]*ux[1]*dx10+0.8660254037844386*rhouz[0]*ux[0]*dx10+1.732050807568877*Pxz[0]*dx10; 
  outrhouz[2] += 0.8660254037844386*rhouz[3]*uy[3]*dx11+0.8660254037844386*rhouz[2]*uy[2]*dx11+0.8660254037844386*rhouz[1]*uy[1]*dx11+0.8660254037844386*rhouz[0]*uy[0]*dx11+1.732050807568877*Pyz[0]*dx11; 
  outrhouz[3] += 0.8660254037844386*rhouz[2]*uy[3]*dx11+0.8660254037844386*uy[2]*rhouz[3]*dx11+0.8660254037844386*rhouz[0]*uy[1]*dx11+0.8660254037844386*uy[0]*rhouz[1]*dx11+1.732050807568877*Pyz[1]*dx11+0.8660254037844386*rhouz[1]*ux[3]*dx10+0.8660254037844386*ux[1]*rhouz[3]*dx10+0.8660254037844386*rhouz[0]*ux[2]*dx10+0.8660254037844386*ux[0]*rhouz[2]*dx10+1.732050807568877*Pxz[2]*dx10; 

  return cflFreq_mid; 
} 
