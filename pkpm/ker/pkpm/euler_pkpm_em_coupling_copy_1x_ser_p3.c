#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void euler_pkpm_em_coupling_copy_1x_ser_p3(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, 
  struct gkyl_nmat *x, 
  const double *vlasov_pkpm_moms[GKYL_MAX_SPECIES], const double *pkpm_u[GKYL_MAX_SPECIES], 
  double* GKYL_RESTRICT euler_pkpm[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em) 
{ 
  // count:       integer to indicate which matrix being fetched. 
  // x:           Input solution vector. 
  // vlasov_pkpm_moms: [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model at old time t^n.
  // pkpm_u:      [ux, uy, uz], Input flow velocity at old time t^n.
  // euler_pkpm:  [rho ux, rho uy, rho uz], Fluid output state vector at time t^{n+1}.
  // em:          [Ex, Ey, Ez, Bx, By, Bz], EM output state vector at time t^{n+1}.
  //              Source solve only updates Ex, Ey, Ez. 

  struct gkyl_mat sol = gkyl_nmat_get(x, count); 
  double ux_new[4] = {0.0}; 
  double uy_new[4] = {0.0}; 
  double uz_new[4] = {0.0}; 
  double Ex_new[4] = {0.0}; 
  double Ey_new[4] = {0.0}; 
  double Ez_new[4] = {0.0}; 

  for (int i = 0; i < num_species; ++i) { 
    const double *rho_old = &vlasov_pkpm_moms[i][0]; 
    const double *ux_old = &pkpm_u[i][0]; 
    const double *uy_old = &pkpm_u[i][4]; 
    const double *uz_old = &pkpm_u[i][8]; 
    ux_new[0] = 2.0*gkyl_mat_get(&sol, 0 + i*(12), 0) - ux_old[0]; 
    uy_new[0] = 2.0*gkyl_mat_get(&sol, 4 + i*(12), 0) - uy_old[0]; 
    uz_new[0] = 2.0*gkyl_mat_get(&sol, 8 + i*(12), 0) - uz_old[0]; 

    ux_new[1] = 2.0*gkyl_mat_get(&sol, 1 + i*(12), 0) - ux_old[1]; 
    uy_new[1] = 2.0*gkyl_mat_get(&sol, 5 + i*(12), 0) - uy_old[1]; 
    uz_new[1] = 2.0*gkyl_mat_get(&sol, 9 + i*(12), 0) - uz_old[1]; 

    ux_new[2] = 2.0*gkyl_mat_get(&sol, 2 + i*(12), 0) - ux_old[2]; 
    uy_new[2] = 2.0*gkyl_mat_get(&sol, 6 + i*(12), 0) - uy_old[2]; 
    uz_new[2] = 2.0*gkyl_mat_get(&sol, 10 + i*(12), 0) - uz_old[2]; 

    ux_new[3] = 2.0*gkyl_mat_get(&sol, 3 + i*(12), 0) - ux_old[3]; 
    uy_new[3] = 2.0*gkyl_mat_get(&sol, 7 + i*(12), 0) - uy_old[3]; 
    uz_new[3] = 2.0*gkyl_mat_get(&sol, 11 + i*(12), 0) - uz_old[3]; 

    double *out_rhoux = &euler_pkpm[i][0]; 
    double *out_rhouy = &euler_pkpm[i][4]; 
    double *out_rhouz = &euler_pkpm[i][8]; 

    binop_mul_1d_ser_p3(rho_old, ux_new, out_rhoux); 
    binop_mul_1d_ser_p3(rho_old, uy_new, out_rhouy); 
    binop_mul_1d_ser_p3(rho_old, uz_new, out_rhouz); 
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
