#include <gkyl_mat.h> 
#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void fluid_em_coupling_copy_1x_ser_p1(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, 
  struct gkyl_nmat *x, double* GKYL_RESTRICT fluid[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // x:      Input solution vector. 
  // fluid:       [rho, rho ux, rho uy, rho uz...], Fluid output state vector.
  //              Source solve only updates momentum of fluid system. 
  //              (isothermal Euler, Euler, Ten moment). 
  // em:          [Ex, Ey, Ez, Bx, By, Bz], EM output state vector.
  //              Source solve only updates Ex, Ey, Ez. 

  struct gkyl_mat sol = gkyl_nmat_get(x, count); 
  double rhoux_new[2] = {0.0}; 
  double rhouy_new[2] = {0.0}; 
  double rhouz_new[2] = {0.0}; 
  double Ex_new[2] = {0.0}; 
  double Ey_new[2] = {0.0}; 
  double Ez_new[2] = {0.0}; 

  for (int i = 0; i < num_species; ++i) { 
    double *out_rhoux = &fluid[i][2]; 
    double *out_rhouy = &fluid[i][4]; 
    double *out_rhouz = &fluid[i][6]; 

    rhoux_new[0] = gkyl_mat_get(&sol, 0 + i*(6), 0); 
    rhouy_new[0] = gkyl_mat_get(&sol, 2 + i*(6), 0); 
    rhouz_new[0] = gkyl_mat_get(&sol, 4 + i*(6), 0); 

    out_rhoux[0] = 2.0*rhoux_new[0]/qbym[i] - out_rhoux[0]; 
    out_rhouy[0] = 2.0*rhouy_new[0]/qbym[i] - out_rhouy[0]; 
    out_rhouz[0] = 2.0*rhouz_new[0]/qbym[i] - out_rhouz[0]; 

    rhoux_new[1] = gkyl_mat_get(&sol, 1 + i*(6), 0); 
    rhouy_new[1] = gkyl_mat_get(&sol, 3 + i*(6), 0); 
    rhouz_new[1] = gkyl_mat_get(&sol, 5 + i*(6), 0); 

    out_rhoux[1] = 2.0*rhoux_new[1]/qbym[i] - out_rhoux[1]; 
    out_rhouy[1] = 2.0*rhouy_new[1]/qbym[i] - out_rhouy[1]; 
    out_rhouz[1] = 2.0*rhouz_new[1]/qbym[i] - out_rhouz[1]; 

  } 

  double *out_Ex = &em[0]; 
  double *out_Ey = &em[2]; 
  double *out_Ez = &em[4]; 
  Ex_new[0] = gkyl_mat_get(&sol, 0 + num_species*(6), 0); 
  Ey_new[0] = gkyl_mat_get(&sol, 2 + num_species*(6), 0); 
  Ez_new[0] = gkyl_mat_get(&sol, 4 + num_species*(6), 0); 

  out_Ex[0] = 2.0*Ex_new[0]/epsilon0 - out_Ex[0]; 
  out_Ey[0] = 2.0*Ey_new[0]/epsilon0 - out_Ey[0]; 
  out_Ez[0] = 2.0*Ez_new[0]/epsilon0 - out_Ez[0]; 

  Ex_new[1] = gkyl_mat_get(&sol, 1 + num_species*(6), 0); 
  Ey_new[1] = gkyl_mat_get(&sol, 3 + num_species*(6), 0); 
  Ez_new[1] = gkyl_mat_get(&sol, 5 + num_species*(6), 0); 

  out_Ex[1] = 2.0*Ex_new[1]/epsilon0 - out_Ex[1]; 
  out_Ey[1] = 2.0*Ey_new[1]/epsilon0 - out_Ey[1]; 
  out_Ez[1] = 2.0*Ez_new[1]/epsilon0 - out_Ez[1]; 

} 
