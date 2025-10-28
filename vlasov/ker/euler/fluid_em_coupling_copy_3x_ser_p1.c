#include <gkyl_mat.h> 
#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void fluid_em_coupling_copy_3x_ser_p1(int count, 
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
  double rhoux_new[8] = {0.0}; 
  double rhouy_new[8] = {0.0}; 
  double rhouz_new[8] = {0.0}; 
  double Ex_new[8] = {0.0}; 
  double Ey_new[8] = {0.0}; 
  double Ez_new[8] = {0.0}; 

  for (int i = 0; i < num_species; ++i) { 
    double *out_rhoux = &fluid[i][8]; 
    double *out_rhouy = &fluid[i][16]; 
    double *out_rhouz = &fluid[i][24]; 

    rhoux_new[0] = gkyl_mat_get(&sol, 0 + i*(24), 0); 
    rhouy_new[0] = gkyl_mat_get(&sol, 8 + i*(24), 0); 
    rhouz_new[0] = gkyl_mat_get(&sol, 16 + i*(24), 0); 

    out_rhoux[0] = 2.0*rhoux_new[0]/qbym[i] - out_rhoux[0]; 
    out_rhouy[0] = 2.0*rhouy_new[0]/qbym[i] - out_rhouy[0]; 
    out_rhouz[0] = 2.0*rhouz_new[0]/qbym[i] - out_rhouz[0]; 

    rhoux_new[1] = gkyl_mat_get(&sol, 1 + i*(24), 0); 
    rhouy_new[1] = gkyl_mat_get(&sol, 9 + i*(24), 0); 
    rhouz_new[1] = gkyl_mat_get(&sol, 17 + i*(24), 0); 

    out_rhoux[1] = 2.0*rhoux_new[1]/qbym[i] - out_rhoux[1]; 
    out_rhouy[1] = 2.0*rhouy_new[1]/qbym[i] - out_rhouy[1]; 
    out_rhouz[1] = 2.0*rhouz_new[1]/qbym[i] - out_rhouz[1]; 

    rhoux_new[2] = gkyl_mat_get(&sol, 2 + i*(24), 0); 
    rhouy_new[2] = gkyl_mat_get(&sol, 10 + i*(24), 0); 
    rhouz_new[2] = gkyl_mat_get(&sol, 18 + i*(24), 0); 

    out_rhoux[2] = 2.0*rhoux_new[2]/qbym[i] - out_rhoux[2]; 
    out_rhouy[2] = 2.0*rhouy_new[2]/qbym[i] - out_rhouy[2]; 
    out_rhouz[2] = 2.0*rhouz_new[2]/qbym[i] - out_rhouz[2]; 

    rhoux_new[3] = gkyl_mat_get(&sol, 3 + i*(24), 0); 
    rhouy_new[3] = gkyl_mat_get(&sol, 11 + i*(24), 0); 
    rhouz_new[3] = gkyl_mat_get(&sol, 19 + i*(24), 0); 

    out_rhoux[3] = 2.0*rhoux_new[3]/qbym[i] - out_rhoux[3]; 
    out_rhouy[3] = 2.0*rhouy_new[3]/qbym[i] - out_rhouy[3]; 
    out_rhouz[3] = 2.0*rhouz_new[3]/qbym[i] - out_rhouz[3]; 

    rhoux_new[4] = gkyl_mat_get(&sol, 4 + i*(24), 0); 
    rhouy_new[4] = gkyl_mat_get(&sol, 12 + i*(24), 0); 
    rhouz_new[4] = gkyl_mat_get(&sol, 20 + i*(24), 0); 

    out_rhoux[4] = 2.0*rhoux_new[4]/qbym[i] - out_rhoux[4]; 
    out_rhouy[4] = 2.0*rhouy_new[4]/qbym[i] - out_rhouy[4]; 
    out_rhouz[4] = 2.0*rhouz_new[4]/qbym[i] - out_rhouz[4]; 

    rhoux_new[5] = gkyl_mat_get(&sol, 5 + i*(24), 0); 
    rhouy_new[5] = gkyl_mat_get(&sol, 13 + i*(24), 0); 
    rhouz_new[5] = gkyl_mat_get(&sol, 21 + i*(24), 0); 

    out_rhoux[5] = 2.0*rhoux_new[5]/qbym[i] - out_rhoux[5]; 
    out_rhouy[5] = 2.0*rhouy_new[5]/qbym[i] - out_rhouy[5]; 
    out_rhouz[5] = 2.0*rhouz_new[5]/qbym[i] - out_rhouz[5]; 

    rhoux_new[6] = gkyl_mat_get(&sol, 6 + i*(24), 0); 
    rhouy_new[6] = gkyl_mat_get(&sol, 14 + i*(24), 0); 
    rhouz_new[6] = gkyl_mat_get(&sol, 22 + i*(24), 0); 

    out_rhoux[6] = 2.0*rhoux_new[6]/qbym[i] - out_rhoux[6]; 
    out_rhouy[6] = 2.0*rhouy_new[6]/qbym[i] - out_rhouy[6]; 
    out_rhouz[6] = 2.0*rhouz_new[6]/qbym[i] - out_rhouz[6]; 

    rhoux_new[7] = gkyl_mat_get(&sol, 7 + i*(24), 0); 
    rhouy_new[7] = gkyl_mat_get(&sol, 15 + i*(24), 0); 
    rhouz_new[7] = gkyl_mat_get(&sol, 23 + i*(24), 0); 

    out_rhoux[7] = 2.0*rhoux_new[7]/qbym[i] - out_rhoux[7]; 
    out_rhouy[7] = 2.0*rhouy_new[7]/qbym[i] - out_rhouy[7]; 
    out_rhouz[7] = 2.0*rhouz_new[7]/qbym[i] - out_rhouz[7]; 

  } 

  double *out_Ex = &em[0]; 
  double *out_Ey = &em[8]; 
  double *out_Ez = &em[16]; 
  Ex_new[0] = gkyl_mat_get(&sol, 0 + num_species*(24), 0); 
  Ey_new[0] = gkyl_mat_get(&sol, 8 + num_species*(24), 0); 
  Ez_new[0] = gkyl_mat_get(&sol, 16 + num_species*(24), 0); 

  out_Ex[0] = 2.0*Ex_new[0]/epsilon0 - out_Ex[0]; 
  out_Ey[0] = 2.0*Ey_new[0]/epsilon0 - out_Ey[0]; 
  out_Ez[0] = 2.0*Ez_new[0]/epsilon0 - out_Ez[0]; 

  Ex_new[1] = gkyl_mat_get(&sol, 1 + num_species*(24), 0); 
  Ey_new[1] = gkyl_mat_get(&sol, 9 + num_species*(24), 0); 
  Ez_new[1] = gkyl_mat_get(&sol, 17 + num_species*(24), 0); 

  out_Ex[1] = 2.0*Ex_new[1]/epsilon0 - out_Ex[1]; 
  out_Ey[1] = 2.0*Ey_new[1]/epsilon0 - out_Ey[1]; 
  out_Ez[1] = 2.0*Ez_new[1]/epsilon0 - out_Ez[1]; 

  Ex_new[2] = gkyl_mat_get(&sol, 2 + num_species*(24), 0); 
  Ey_new[2] = gkyl_mat_get(&sol, 10 + num_species*(24), 0); 
  Ez_new[2] = gkyl_mat_get(&sol, 18 + num_species*(24), 0); 

  out_Ex[2] = 2.0*Ex_new[2]/epsilon0 - out_Ex[2]; 
  out_Ey[2] = 2.0*Ey_new[2]/epsilon0 - out_Ey[2]; 
  out_Ez[2] = 2.0*Ez_new[2]/epsilon0 - out_Ez[2]; 

  Ex_new[3] = gkyl_mat_get(&sol, 3 + num_species*(24), 0); 
  Ey_new[3] = gkyl_mat_get(&sol, 11 + num_species*(24), 0); 
  Ez_new[3] = gkyl_mat_get(&sol, 19 + num_species*(24), 0); 

  out_Ex[3] = 2.0*Ex_new[3]/epsilon0 - out_Ex[3]; 
  out_Ey[3] = 2.0*Ey_new[3]/epsilon0 - out_Ey[3]; 
  out_Ez[3] = 2.0*Ez_new[3]/epsilon0 - out_Ez[3]; 

  Ex_new[4] = gkyl_mat_get(&sol, 4 + num_species*(24), 0); 
  Ey_new[4] = gkyl_mat_get(&sol, 12 + num_species*(24), 0); 
  Ez_new[4] = gkyl_mat_get(&sol, 20 + num_species*(24), 0); 

  out_Ex[4] = 2.0*Ex_new[4]/epsilon0 - out_Ex[4]; 
  out_Ey[4] = 2.0*Ey_new[4]/epsilon0 - out_Ey[4]; 
  out_Ez[4] = 2.0*Ez_new[4]/epsilon0 - out_Ez[4]; 

  Ex_new[5] = gkyl_mat_get(&sol, 5 + num_species*(24), 0); 
  Ey_new[5] = gkyl_mat_get(&sol, 13 + num_species*(24), 0); 
  Ez_new[5] = gkyl_mat_get(&sol, 21 + num_species*(24), 0); 

  out_Ex[5] = 2.0*Ex_new[5]/epsilon0 - out_Ex[5]; 
  out_Ey[5] = 2.0*Ey_new[5]/epsilon0 - out_Ey[5]; 
  out_Ez[5] = 2.0*Ez_new[5]/epsilon0 - out_Ez[5]; 

  Ex_new[6] = gkyl_mat_get(&sol, 6 + num_species*(24), 0); 
  Ey_new[6] = gkyl_mat_get(&sol, 14 + num_species*(24), 0); 
  Ez_new[6] = gkyl_mat_get(&sol, 22 + num_species*(24), 0); 

  out_Ex[6] = 2.0*Ex_new[6]/epsilon0 - out_Ex[6]; 
  out_Ey[6] = 2.0*Ey_new[6]/epsilon0 - out_Ey[6]; 
  out_Ez[6] = 2.0*Ez_new[6]/epsilon0 - out_Ez[6]; 

  Ex_new[7] = gkyl_mat_get(&sol, 7 + num_species*(24), 0); 
  Ey_new[7] = gkyl_mat_get(&sol, 15 + num_species*(24), 0); 
  Ez_new[7] = gkyl_mat_get(&sol, 23 + num_species*(24), 0); 

  out_Ex[7] = 2.0*Ex_new[7]/epsilon0 - out_Ex[7]; 
  out_Ey[7] = 2.0*Ey_new[7]/epsilon0 - out_Ey[7]; 
  out_Ez[7] = 2.0*Ez_new[7]/epsilon0 - out_Ez[7]; 

} 
