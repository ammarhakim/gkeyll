#include <gkyl_euler_kernels.h> 
GKYL_CU_DH double euler_vol_1x_ser_p2(const double *w, const double *dxv, double gas_gamma, 
    const double *u, const double *p, const double *fluid, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // gas_gamma: Adiabatic index.
  // u:         Input flow velocity [ux, uy, uz].
  // p:         Input pressure .
  // fluid:     [rho, rho ux, rho uy, rho uz, E], Fluid input state vector.
  // out:       Incremented output.

  double dx10 = 2./dxv[0]; 

  const double *rho = &fluid[0]; 
  const double *rhoux = &fluid[3]; 
  const double *rhouy = &fluid[6]; 
  const double *rhouz = &fluid[9]; 
  const double *energy = &fluid[12]; 

  const double *ux = &u[0]; 
  const double *uy = &u[3]; 
  const double *uz = &u[6]; 

  double *outrho = &out[0]; 
  double *outrhoux = &out[3]; 
  double *outrhouy = &out[6]; 
  double *outrhouz = &out[9]; 
  double *outenergy = &out[12]; 

  double cflFreq_mid = 0.0; 
  cflFreq_mid += 0.5*5.0*dx10*(fabs(0.7071067811865475*ux[0]-0.7905694150420947*ux[2]) + sqrt(gas_gamma*0.7071067811865475*p[0]-0.7905694150420947*p[2]/(0.7071067811865475*rho[0]-0.7905694150420947*rho[2]))); 

  outrho[1] += 1.732050807568877*rhoux[0]*dx10; 
  outrho[2] += 3.872983346207417*rhoux[1]*dx10; 

  outrhoux[1] += 1.224744871391589*rhoux[2]*ux[2]*dx10+1.224744871391589*rhoux[1]*ux[1]*dx10+1.224744871391589*rhoux[0]*ux[0]*dx10+1.732050807568877*p[0]*dx10; 
  outrhoux[2] += 2.449489742783178*rhoux[1]*ux[2]*dx10+2.449489742783178*ux[1]*rhoux[2]*dx10+2.738612787525831*rhoux[0]*ux[1]*dx10+2.738612787525831*ux[0]*rhoux[1]*dx10+3.872983346207417*p[1]*dx10; 

  outrhouy[1] += 1.224744871391589*rhouy[2]*ux[2]*dx10+1.224744871391589*rhouy[1]*ux[1]*dx10+1.224744871391589*rhouy[0]*ux[0]*dx10; 
  outrhouy[2] += 2.449489742783178*rhouy[1]*ux[2]*dx10+2.449489742783178*ux[1]*rhouy[2]*dx10+2.738612787525831*rhouy[0]*ux[1]*dx10+2.738612787525831*ux[0]*rhouy[1]*dx10; 

  outrhouz[1] += 1.224744871391589*rhouz[2]*ux[2]*dx10+1.224744871391589*rhouz[1]*ux[1]*dx10+1.224744871391589*rhouz[0]*ux[0]*dx10; 
  outrhouz[2] += 2.449489742783178*rhouz[1]*ux[2]*dx10+2.449489742783178*ux[1]*rhouz[2]*dx10+2.738612787525831*rhouz[0]*ux[1]*dx10+2.738612787525831*ux[0]*rhouz[1]*dx10; 

  outenergy[1] += 1.224744871391589*p[2]*ux[2]*dx10+1.224744871391589*energy[2]*ux[2]*dx10+1.224744871391589*p[1]*ux[1]*dx10+1.224744871391589*energy[1]*ux[1]*dx10+1.224744871391589*p[0]*ux[0]*dx10+1.224744871391589*energy[0]*ux[0]*dx10; 
  outenergy[2] += 2.449489742783178*p[1]*ux[2]*dx10+2.449489742783178*energy[1]*ux[2]*dx10+2.449489742783178*ux[1]*p[2]*dx10+2.449489742783178*ux[1]*energy[2]*dx10+2.738612787525831*p[0]*ux[1]*dx10+2.738612787525831*energy[0]*ux[1]*dx10+2.738612787525831*ux[0]*p[1]*dx10+2.738612787525831*ux[0]*energy[1]*dx10; 

  return cflFreq_mid; 
} 
