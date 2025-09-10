#include <gkyl_mat.h> 
#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void fluid_em_coupling_copy_1x_ser_p3(int count, 
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
  double rhoux_new[4] = {0.0}; 
  double rhouy_new[4] = {0.0}; 
  double rhouz_new[4] = {0.0}; 
  double Ex_new[4] = {0.0}; 
  double Ey_new[4] = {0.0}; 
  double Ez_new[4] = {0.0}; 

  for (int i = 0; i < num_species; ++i) { 
    double *out_rhoux = &fluid[i][4]; 
    double *out_rhouy = &fluid[i][8]; 
    double *out_rhouz = &fluid[i][12]; 

    rhoux_new[0] = gkyl_mat_get(&sol, 0 + i*(12), 0); 
    rhouy_new[0] = gkyl_mat_get(&sol, 4 + i*(12), 0); 
    rhouz_new[0] = gkyl_mat_get(&sol, 8 + i*(12), 0); 

    out_rhoux[0] = 2.0*rhoux_new[0]/qbym[i] - out_rhoux[0]; 
    out_rhouy[0] = 2.0*rhouy_new[0]/qbym[i] - out_rhouy[0]; 
    out_rhouz[0] = 2.0*rhouz_new[0]/qbym[i] - out_rhouz[0]; 

    rhoux_new[1] = gkyl_mat_get(&sol, 1 + i*(12), 0); 
    rhouy_new[1] = gkyl_mat_get(&sol, 5 + i*(12), 0); 
    rhouz_new[1] = gkyl_mat_get(&sol, 9 + i*(12), 0); 

    out_rhoux[1] = 2.0*rhoux_new[1]/qbym[i] - out_rhoux[1]; 
    out_rhouy[1] = 2.0*rhouy_new[1]/qbym[i] - out_rhouy[1]; 
    out_rhouz[1] = 2.0*rhouz_new[1]/qbym[i] - out_rhouz[1]; 

    rhoux_new[2] = gkyl_mat_get(&sol, 2 + i*(12), 0); 
    rhouy_new[2] = gkyl_mat_get(&sol, 6 + i*(12), 0); 
    rhouz_new[2] = gkyl_mat_get(&sol, 10 + i*(12), 0); 

    out_rhoux[2] = 2.0*rhoux_new[2]/qbym[i] - out_rhoux[2]; 
    out_rhouy[2] = 2.0*rhouy_new[2]/qbym[i] - out_rhouy[2]; 
    out_rhouz[2] = 2.0*rhouz_new[2]/qbym[i] - out_rhouz[2]; 

    rhoux_new[3] = gkyl_mat_get(&sol, 3 + i*(12), 0); 
    rhouy_new[3] = gkyl_mat_get(&sol, 7 + i*(12), 0); 
    rhouz_new[3] = gkyl_mat_get(&sol, 11 + i*(12), 0); 

    out_rhoux[3] = 2.0*rhoux_new[3]/qbym[i] - out_rhoux[3]; 
    out_rhouy[3] = 2.0*rhouy_new[3]/qbym[i] - out_rhouy[3]; 
    out_rhouz[3] = 2.0*rhouz_new[3]/qbym[i] - out_rhouz[3]; 

  } 

  double *out_Ex = &em[0]; 
  double *out_Ey = &em[4]; 
  double *out_Ez = &em[8]; 
  Ex_new[0] = gkyl_mat_get(&sol, 0 + num_species*(12), 0); 
  Ey_new[0] = gkyl_mat_get(&sol, 4 + num_species*(12), 0); 
  Ez_new[0] = gkyl_mat_get(&sol, 8 + num_species*(12), 0); 

  out_Ex[0] = 2.0*Ex_new[0]/epsilon0 - out_Ex[0]; 
  out_Ey[0] = 2.0*Ey_new[0]/epsilon0 - out_Ey[0]; 
  out_Ez[0] = 2.0*Ez_new[0]/epsilon0 - out_Ez[0]; 

  Ex_new[1] = gkyl_mat_get(&sol, 1 + num_species*(12), 0); 
  Ey_new[1] = gkyl_mat_get(&sol, 5 + num_species*(12), 0); 
  Ez_new[1] = gkyl_mat_get(&sol, 9 + num_species*(12), 0); 

  out_Ex[1] = 2.0*Ex_new[1]/epsilon0 - out_Ex[1]; 
  out_Ey[1] = 2.0*Ey_new[1]/epsilon0 - out_Ey[1]; 
  out_Ez[1] = 2.0*Ez_new[1]/epsilon0 - out_Ez[1]; 

  Ex_new[2] = gkyl_mat_get(&sol, 2 + num_species*(12), 0); 
  Ey_new[2] = gkyl_mat_get(&sol, 6 + num_species*(12), 0); 
  Ez_new[2] = gkyl_mat_get(&sol, 10 + num_species*(12), 0); 

  out_Ex[2] = 2.0*Ex_new[2]/epsilon0 - out_Ex[2]; 
  out_Ey[2] = 2.0*Ey_new[2]/epsilon0 - out_Ey[2]; 
  out_Ez[2] = 2.0*Ez_new[2]/epsilon0 - out_Ez[2]; 

  Ex_new[3] = gkyl_mat_get(&sol, 3 + num_species*(12), 0); 
  Ey_new[3] = gkyl_mat_get(&sol, 7 + num_species*(12), 0); 
  Ez_new[3] = gkyl_mat_get(&sol, 11 + num_species*(12), 0); 

  out_Ex[3] = 2.0*Ex_new[3]/epsilon0 - out_Ex[3]; 
  out_Ey[3] = 2.0*Ey_new[3]/epsilon0 - out_Ey[3]; 
  out_Ez[3] = 2.0*Ez_new[3]/epsilon0 - out_Ez[3]; 

} 
