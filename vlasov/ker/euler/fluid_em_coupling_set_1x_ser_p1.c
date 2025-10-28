#include <gkyl_mat.h> 
#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void fluid_em_coupling_set_1x_ser_p1(int count, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0, double dt, 
  struct gkyl_nmat *A_n, struct gkyl_nmat *rhs_n, 
  const double *app_accel[GKYL_MAX_SPECIES], const double *ext_em, const double *app_current, 
  double* GKYL_RESTRICT fluid[GKYL_MAX_SPECIES], double* GKYL_RESTRICT em) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // A:           preallocated LHS matrix. 
  // rhs:         preallocated RHS vector. 
  // app_accel:   Applied accelerations (external forces).
  // ext_em:      Externally applied EM fields.
  // app_current: Applied external currents.
  // fluid:       [rho, rho ux, rho uy, rho uz...], Fluid input state vector.
  //              only need rho and momentum to get update source terms of fluid system. 
  //              (isothermal Euler, Euler, Ten moment). 
  // em:          [Ex, Ey, Ez, Bx, By, Bz], EM input state vector.

  struct gkyl_mat lhs = gkyl_nmat_get(A_n, count); 
  struct gkyl_mat rhs = gkyl_nmat_get(rhs_n, count); 
  // Clear matrix and rhs source solve. 
  gkyl_mat_clear(&lhs, 0.0); gkyl_mat_clear(&rhs, 0.0); 

  double rho[GKYL_MAX_SPECIES][2]; 
  double rhoux[GKYL_MAX_SPECIES][2]; 
  double rhouy[GKYL_MAX_SPECIES][2]; 
  double rhouz[GKYL_MAX_SPECIES][2]; 

  double app_accel_x[GKYL_MAX_SPECIES][2]; 
  double app_accel_y[GKYL_MAX_SPECIES][2]; 
  double app_accel_z[GKYL_MAX_SPECIES][2]; 

  for (int i = 0; i < num_species; ++i) { 
    double *inp_fluid = fluid[i]; 
    const double *inp_app_accel = app_accel[i]; 

    rho[i][0] = inp_fluid[0]; 
    rhoux[i][0] = inp_fluid[2]; 
    rhouy[i][0] = inp_fluid[4]; 
    rhouz[i][0] = inp_fluid[6]; 

    app_accel_x[i][0] = inp_app_accel[0]; 
    app_accel_y[i][0] = inp_app_accel[2]; 
    app_accel_z[i][0] = inp_app_accel[4]; 

    rho[i][1] = inp_fluid[1]; 
    rhoux[i][1] = inp_fluid[3]; 
    rhouy[i][1] = inp_fluid[5]; 
    rhouz[i][1] = inp_fluid[7]; 

    app_accel_x[i][1] = inp_app_accel[1]; 
    app_accel_y[i][1] = inp_app_accel[3]; 
    app_accel_z[i][1] = inp_app_accel[5]; 

  } 

  double *Ex = &em[0]; 
  double *Ey = &em[2]; 
  double *Ez = &em[4]; 
  double *Bx = &em[6]; 
  double *By = &em[8]; 
  double *Bz = &em[10]; 

  const double *ext_Ex = &ext_em[0]; 
  const double *ext_Ey = &ext_em[2]; 
  const double *ext_Ez = &ext_em[4]; 
  const double *ext_Bx = &ext_em[6]; 
  const double *ext_By = &ext_em[8]; 
  const double *ext_Bz = &ext_em[10]; 

  const double *app_curr_x = &app_current[0]; 
  const double *app_curr_y = &app_current[2]; 
  const double *app_curr_z = &app_current[4]; 

  double tot_Bx[2]; 
  double tot_By[2]; 
  double tot_Bz[2]; 
  tot_Bx[0] = Bx[0] + ext_Bx[0]; 
  tot_By[0] = By[0] + ext_By[0]; 
  tot_Bz[0] = Bz[0] + ext_Bz[0]; 
  tot_Bx[1] = Bx[1] + ext_Bx[1]; 
  tot_By[1] = By[1] + ext_By[1]; 
  tot_Bz[1] = Bz[1] + ext_Bz[1]; 

  // Set RHS for momentum equations, including solution at known time-step and external forces. 
  for (int i = 0; i < num_species; ++i) { 

    gkyl_mat_set(&rhs, 0 + i*(6), 0, qbym[i]*rhoux[i][0] + 0.5*dt*qbym[i]*rho[i][0]*(qbym[i]*ext_Ex[0] + app_accel_x[i][0])); 
    gkyl_mat_set(&rhs, 2 + i*(6), 0, qbym[i]*rhouy[i][0] + 0.5*dt*qbym[i]*rho[i][0]*(qbym[i]*ext_Ey[0] + app_accel_y[i][0])); 
    gkyl_mat_set(&rhs, 4 + i*(6), 0, qbym[i]*rhouz[i][0] + 0.5*dt*qbym[i]*rho[i][0]*(qbym[i]*ext_Ez[0] + app_accel_z[i][0])); 

    gkyl_mat_set(&rhs, 1 + i*(6), 0, qbym[i]*rhoux[i][1] + 0.5*dt*qbym[i]*rho[i][1]*(qbym[i]*ext_Ex[1] + app_accel_x[i][1])); 
    gkyl_mat_set(&rhs, 3 + i*(6), 0, qbym[i]*rhouy[i][1] + 0.5*dt*qbym[i]*rho[i][1]*(qbym[i]*ext_Ey[1] + app_accel_y[i][1])); 
    gkyl_mat_set(&rhs, 5 + i*(6), 0, qbym[i]*rhouz[i][1] + 0.5*dt*qbym[i]*rho[i][1]*(qbym[i]*ext_Ez[1] + app_accel_z[i][1])); 

  } 

  // Set RHS for Ampere's Law, including solution at known time-step and applied currents. 
  gkyl_mat_set(&rhs, 0 + num_species*(6), 0, epsilon0*Ex[0] - 0.5*dt*app_curr_x[0]); 
  gkyl_mat_set(&rhs, 2 + num_species*(6), 0, epsilon0*Ey[0] - 0.5*dt*app_curr_y[0]); 
  gkyl_mat_set(&rhs, 4 + num_species*(6), 0, epsilon0*Ez[0] - 0.5*dt*app_curr_z[0]); 

  gkyl_mat_set(&rhs, 1 + num_species*(6), 0, epsilon0*Ex[1] - 0.5*dt*app_curr_x[1]); 
  gkyl_mat_set(&rhs, 3 + num_species*(6), 0, epsilon0*Ey[1] - 0.5*dt*app_curr_y[1]); 
  gkyl_mat_set(&rhs, 5 + num_species*(6), 0, epsilon0*Ez[1] - 0.5*dt*app_curr_z[1]); 


  // Construct LHS. 
  // For momentum equation: J_s^{n+1} - 0.5*dt*(q_s^2/m_s^2*rho_s^n*E^{n+1} + q_s/m_s*J_s^{n+1} x B^n). 
  // For Ampere's Law: epsilon0*E^{n+1} + 0.5*dt*sum_s J_s^{n+1}. 
  for (int s = 0; s < num_species; ++s) { 
 
    double E_field_fac = -0.5*dt*qbym[s]*qbym[s]/epsilon0; 
    double B_field_fac = -0.5*dt*qbym[s]; 
    gkyl_mat_set(&lhs, 0 + s*(6), 0 + s*(6), 1.0); 
    gkyl_mat_set(&lhs, 2 + s*(6), 2 + s*(6), 1.0); 
    gkyl_mat_set(&lhs, 4 + s*(6), 4 + s*(6), 1.0); 
 
    gkyl_mat_set(&lhs, 0 + s*(6), 0 + num_species*(6), E_field_fac*(0.7071067811865475*rho[s][0])); 
    gkyl_mat_set(&lhs, 2 + s*(6), 2 + num_species*(6), E_field_fac*(0.7071067811865475*rho[s][0])); 
    gkyl_mat_set(&lhs, 4 + s*(6), 4 + num_species*(6), E_field_fac*(0.7071067811865475*rho[s][0])); 
 
    gkyl_mat_set(&lhs, 2 + s*(6), 4 + s*(6), B_field_fac*(0.7071067811865475*tot_Bx[0])); 
    gkyl_mat_set(&lhs, 4 + s*(6), 2 + s*(6), -B_field_fac*(0.7071067811865475*tot_Bx[0])); 
 
    gkyl_mat_set(&lhs, 4 + s*(6), 0 + s*(6), B_field_fac*(0.7071067811865475*tot_By[0])); 
    gkyl_mat_set(&lhs, 0 + s*(6), 4 + s*(6), -B_field_fac*(0.7071067811865475*tot_By[0])); 
 
    gkyl_mat_set(&lhs, 0 + s*(6), 2 + s*(6), B_field_fac*(0.7071067811865475*tot_Bz[0])); 
    gkyl_mat_set(&lhs, 2 + s*(6), 0 + s*(6), -B_field_fac*(0.7071067811865475*tot_Bz[0])); 
 
    gkyl_mat_set(&lhs, 0 + num_species*(6), 0 + s*(6), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 2 + num_species*(6), 2 + s*(6), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 4 + num_species*(6), 4 + s*(6), 0.5*dt*(1.0)); 
 
    gkyl_mat_set(&lhs, 0 + s*(6), 1 + s*(6), 0.0); 
    gkyl_mat_set(&lhs, 2 + s*(6), 3 + s*(6), 0.0); 
    gkyl_mat_set(&lhs, 4 + s*(6), 5 + s*(6), 0.0); 
 
    gkyl_mat_set(&lhs, 0 + s*(6), 1 + num_species*(6), E_field_fac*(0.7071067811865475*rho[s][1])); 
    gkyl_mat_set(&lhs, 2 + s*(6), 3 + num_species*(6), E_field_fac*(0.7071067811865475*rho[s][1])); 
    gkyl_mat_set(&lhs, 4 + s*(6), 5 + num_species*(6), E_field_fac*(0.7071067811865475*rho[s][1])); 
 
    gkyl_mat_set(&lhs, 2 + s*(6), 5 + s*(6), B_field_fac*(0.7071067811865475*tot_Bx[1])); 
    gkyl_mat_set(&lhs, 4 + s*(6), 3 + s*(6), -B_field_fac*(0.7071067811865475*tot_Bx[1])); 
 
    gkyl_mat_set(&lhs, 4 + s*(6), 1 + s*(6), B_field_fac*(0.7071067811865475*tot_By[1])); 
    gkyl_mat_set(&lhs, 0 + s*(6), 5 + s*(6), -B_field_fac*(0.7071067811865475*tot_By[1])); 
 
    gkyl_mat_set(&lhs, 0 + s*(6), 3 + s*(6), B_field_fac*(0.7071067811865475*tot_Bz[1])); 
    gkyl_mat_set(&lhs, 2 + s*(6), 1 + s*(6), -B_field_fac*(0.7071067811865475*tot_Bz[1])); 
 
    gkyl_mat_set(&lhs, 0 + num_species*(6), 1 + s*(6), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 2 + num_species*(6), 3 + s*(6), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 4 + num_species*(6), 5 + s*(6), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 1 + s*(6), 0 + s*(6), 0.0); 
    gkyl_mat_set(&lhs, 3 + s*(6), 2 + s*(6), 0.0); 
    gkyl_mat_set(&lhs, 5 + s*(6), 4 + s*(6), 0.0); 
 
    gkyl_mat_set(&lhs, 1 + s*(6), 0 + num_species*(6), E_field_fac*(0.7071067811865475*rho[s][1])); 
    gkyl_mat_set(&lhs, 3 + s*(6), 2 + num_species*(6), E_field_fac*(0.7071067811865475*rho[s][1])); 
    gkyl_mat_set(&lhs, 5 + s*(6), 4 + num_species*(6), E_field_fac*(0.7071067811865475*rho[s][1])); 
 
    gkyl_mat_set(&lhs, 3 + s*(6), 4 + s*(6), B_field_fac*(0.7071067811865475*tot_Bx[1])); 
    gkyl_mat_set(&lhs, 5 + s*(6), 2 + s*(6), -B_field_fac*(0.7071067811865475*tot_Bx[1])); 
 
    gkyl_mat_set(&lhs, 5 + s*(6), 0 + s*(6), B_field_fac*(0.7071067811865475*tot_By[1])); 
    gkyl_mat_set(&lhs, 1 + s*(6), 4 + s*(6), -B_field_fac*(0.7071067811865475*tot_By[1])); 
 
    gkyl_mat_set(&lhs, 1 + s*(6), 2 + s*(6), B_field_fac*(0.7071067811865475*tot_Bz[1])); 
    gkyl_mat_set(&lhs, 3 + s*(6), 0 + s*(6), -B_field_fac*(0.7071067811865475*tot_Bz[1])); 
 
    gkyl_mat_set(&lhs, 1 + num_species*(6), 0 + s*(6), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 3 + num_species*(6), 2 + s*(6), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 5 + num_species*(6), 4 + s*(6), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 1 + s*(6), 1 + s*(6), 1.0); 
    gkyl_mat_set(&lhs, 3 + s*(6), 3 + s*(6), 1.0); 
    gkyl_mat_set(&lhs, 5 + s*(6), 5 + s*(6), 1.0); 
 
    gkyl_mat_set(&lhs, 1 + s*(6), 1 + num_species*(6), E_field_fac*(0.7071067811865475*rho[s][0])); 
    gkyl_mat_set(&lhs, 3 + s*(6), 3 + num_species*(6), E_field_fac*(0.7071067811865475*rho[s][0])); 
    gkyl_mat_set(&lhs, 5 + s*(6), 5 + num_species*(6), E_field_fac*(0.7071067811865475*rho[s][0])); 
 
    gkyl_mat_set(&lhs, 3 + s*(6), 5 + s*(6), B_field_fac*(0.7071067811865475*tot_Bx[0])); 
    gkyl_mat_set(&lhs, 5 + s*(6), 3 + s*(6), -B_field_fac*(0.7071067811865475*tot_Bx[0])); 
 
    gkyl_mat_set(&lhs, 5 + s*(6), 1 + s*(6), B_field_fac*(0.7071067811865475*tot_By[0])); 
    gkyl_mat_set(&lhs, 1 + s*(6), 5 + s*(6), -B_field_fac*(0.7071067811865475*tot_By[0])); 
 
    gkyl_mat_set(&lhs, 1 + s*(6), 3 + s*(6), B_field_fac*(0.7071067811865475*tot_Bz[0])); 
    gkyl_mat_set(&lhs, 3 + s*(6), 1 + s*(6), -B_field_fac*(0.7071067811865475*tot_Bz[0])); 
 
    gkyl_mat_set(&lhs, 1 + num_species*(6), 1 + s*(6), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 3 + num_species*(6), 3 + s*(6), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 5 + num_species*(6), 5 + s*(6), 0.5*dt*(1.0)); 
 
  } 
  gkyl_mat_set(&lhs, 0 + num_species*(6), 0 + num_species*(6), 1.0); 
  gkyl_mat_set(&lhs, 2 + num_species*(6), 2 + num_species*(6), 1.0); 
  gkyl_mat_set(&lhs, 4 + num_species*(6), 4 + num_species*(6), 1.0); 
 
  gkyl_mat_set(&lhs, 0 + num_species*(6), 1 + num_species*(6), 0.0); 
  gkyl_mat_set(&lhs, 2 + num_species*(6), 3 + num_species*(6), 0.0); 
  gkyl_mat_set(&lhs, 4 + num_species*(6), 5 + num_species*(6), 0.0); 
 
  gkyl_mat_set(&lhs, 1 + num_species*(6), 0 + num_species*(6), 0.0); 
  gkyl_mat_set(&lhs, 3 + num_species*(6), 2 + num_species*(6), 0.0); 
  gkyl_mat_set(&lhs, 5 + num_species*(6), 4 + num_species*(6), 0.0); 
 
  gkyl_mat_set(&lhs, 1 + num_species*(6), 1 + num_species*(6), 1.0); 
  gkyl_mat_set(&lhs, 3 + num_species*(6), 3 + num_species*(6), 1.0); 
  gkyl_mat_set(&lhs, 5 + num_species*(6), 5 + num_species*(6), 1.0); 
 
} 
