#include <gkyl_mat.h> 
#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void fluid_em_coupling_set_3x_ser_p1(int count, 
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

  double rho[GKYL_MAX_SPECIES][8]; 
  double rhoux[GKYL_MAX_SPECIES][8]; 
  double rhouy[GKYL_MAX_SPECIES][8]; 
  double rhouz[GKYL_MAX_SPECIES][8]; 

  double app_accel_x[GKYL_MAX_SPECIES][8]; 
  double app_accel_y[GKYL_MAX_SPECIES][8]; 
  double app_accel_z[GKYL_MAX_SPECIES][8]; 

  for (int i = 0; i < num_species; ++i) { 
    double *inp_fluid = fluid[i]; 
    const double *inp_app_accel = app_accel[i]; 

    rho[i][0] = inp_fluid[0]; 
    rhoux[i][0] = inp_fluid[8]; 
    rhouy[i][0] = inp_fluid[16]; 
    rhouz[i][0] = inp_fluid[24]; 

    app_accel_x[i][0] = inp_app_accel[0]; 
    app_accel_y[i][0] = inp_app_accel[8]; 
    app_accel_z[i][0] = inp_app_accel[16]; 

    rho[i][1] = inp_fluid[1]; 
    rhoux[i][1] = inp_fluid[9]; 
    rhouy[i][1] = inp_fluid[17]; 
    rhouz[i][1] = inp_fluid[25]; 

    app_accel_x[i][1] = inp_app_accel[1]; 
    app_accel_y[i][1] = inp_app_accel[9]; 
    app_accel_z[i][1] = inp_app_accel[17]; 

    rho[i][2] = inp_fluid[2]; 
    rhoux[i][2] = inp_fluid[10]; 
    rhouy[i][2] = inp_fluid[18]; 
    rhouz[i][2] = inp_fluid[26]; 

    app_accel_x[i][2] = inp_app_accel[2]; 
    app_accel_y[i][2] = inp_app_accel[10]; 
    app_accel_z[i][2] = inp_app_accel[18]; 

    rho[i][3] = inp_fluid[3]; 
    rhoux[i][3] = inp_fluid[11]; 
    rhouy[i][3] = inp_fluid[19]; 
    rhouz[i][3] = inp_fluid[27]; 

    app_accel_x[i][3] = inp_app_accel[3]; 
    app_accel_y[i][3] = inp_app_accel[11]; 
    app_accel_z[i][3] = inp_app_accel[19]; 

    rho[i][4] = inp_fluid[4]; 
    rhoux[i][4] = inp_fluid[12]; 
    rhouy[i][4] = inp_fluid[20]; 
    rhouz[i][4] = inp_fluid[28]; 

    app_accel_x[i][4] = inp_app_accel[4]; 
    app_accel_y[i][4] = inp_app_accel[12]; 
    app_accel_z[i][4] = inp_app_accel[20]; 

    rho[i][5] = inp_fluid[5]; 
    rhoux[i][5] = inp_fluid[13]; 
    rhouy[i][5] = inp_fluid[21]; 
    rhouz[i][5] = inp_fluid[29]; 

    app_accel_x[i][5] = inp_app_accel[5]; 
    app_accel_y[i][5] = inp_app_accel[13]; 
    app_accel_z[i][5] = inp_app_accel[21]; 

    rho[i][6] = inp_fluid[6]; 
    rhoux[i][6] = inp_fluid[14]; 
    rhouy[i][6] = inp_fluid[22]; 
    rhouz[i][6] = inp_fluid[30]; 

    app_accel_x[i][6] = inp_app_accel[6]; 
    app_accel_y[i][6] = inp_app_accel[14]; 
    app_accel_z[i][6] = inp_app_accel[22]; 

    rho[i][7] = inp_fluid[7]; 
    rhoux[i][7] = inp_fluid[15]; 
    rhouy[i][7] = inp_fluid[23]; 
    rhouz[i][7] = inp_fluid[31]; 

    app_accel_x[i][7] = inp_app_accel[7]; 
    app_accel_y[i][7] = inp_app_accel[15]; 
    app_accel_z[i][7] = inp_app_accel[23]; 

  } 

  double *Ex = &em[0]; 
  double *Ey = &em[8]; 
  double *Ez = &em[16]; 
  double *Bx = &em[24]; 
  double *By = &em[32]; 
  double *Bz = &em[40]; 

  const double *ext_Ex = &ext_em[0]; 
  const double *ext_Ey = &ext_em[8]; 
  const double *ext_Ez = &ext_em[16]; 
  const double *ext_Bx = &ext_em[24]; 
  const double *ext_By = &ext_em[32]; 
  const double *ext_Bz = &ext_em[40]; 

  const double *app_curr_x = &app_current[0]; 
  const double *app_curr_y = &app_current[8]; 
  const double *app_curr_z = &app_current[16]; 

  double tot_Bx[8]; 
  double tot_By[8]; 
  double tot_Bz[8]; 
  tot_Bx[0] = Bx[0] + ext_Bx[0]; 
  tot_By[0] = By[0] + ext_By[0]; 
  tot_Bz[0] = Bz[0] + ext_Bz[0]; 
  tot_Bx[1] = Bx[1] + ext_Bx[1]; 
  tot_By[1] = By[1] + ext_By[1]; 
  tot_Bz[1] = Bz[1] + ext_Bz[1]; 
  tot_Bx[2] = Bx[2] + ext_Bx[2]; 
  tot_By[2] = By[2] + ext_By[2]; 
  tot_Bz[2] = Bz[2] + ext_Bz[2]; 
  tot_Bx[3] = Bx[3] + ext_Bx[3]; 
  tot_By[3] = By[3] + ext_By[3]; 
  tot_Bz[3] = Bz[3] + ext_Bz[3]; 
  tot_Bx[4] = Bx[4] + ext_Bx[4]; 
  tot_By[4] = By[4] + ext_By[4]; 
  tot_Bz[4] = Bz[4] + ext_Bz[4]; 
  tot_Bx[5] = Bx[5] + ext_Bx[5]; 
  tot_By[5] = By[5] + ext_By[5]; 
  tot_Bz[5] = Bz[5] + ext_Bz[5]; 
  tot_Bx[6] = Bx[6] + ext_Bx[6]; 
  tot_By[6] = By[6] + ext_By[6]; 
  tot_Bz[6] = Bz[6] + ext_Bz[6]; 
  tot_Bx[7] = Bx[7] + ext_Bx[7]; 
  tot_By[7] = By[7] + ext_By[7]; 
  tot_Bz[7] = Bz[7] + ext_Bz[7]; 

  // Set RHS for momentum equations, including solution at known time-step and external forces. 
  for (int i = 0; i < num_species; ++i) { 

    gkyl_mat_set(&rhs, 0 + i*(24), 0, qbym[i]*rhoux[i][0] + 0.5*dt*qbym[i]*rho[i][0]*(qbym[i]*ext_Ex[0] + app_accel_x[i][0])); 
    gkyl_mat_set(&rhs, 8 + i*(24), 0, qbym[i]*rhouy[i][0] + 0.5*dt*qbym[i]*rho[i][0]*(qbym[i]*ext_Ey[0] + app_accel_y[i][0])); 
    gkyl_mat_set(&rhs, 16 + i*(24), 0, qbym[i]*rhouz[i][0] + 0.5*dt*qbym[i]*rho[i][0]*(qbym[i]*ext_Ez[0] + app_accel_z[i][0])); 

    gkyl_mat_set(&rhs, 1 + i*(24), 0, qbym[i]*rhoux[i][1] + 0.5*dt*qbym[i]*rho[i][1]*(qbym[i]*ext_Ex[1] + app_accel_x[i][1])); 
    gkyl_mat_set(&rhs, 9 + i*(24), 0, qbym[i]*rhouy[i][1] + 0.5*dt*qbym[i]*rho[i][1]*(qbym[i]*ext_Ey[1] + app_accel_y[i][1])); 
    gkyl_mat_set(&rhs, 17 + i*(24), 0, qbym[i]*rhouz[i][1] + 0.5*dt*qbym[i]*rho[i][1]*(qbym[i]*ext_Ez[1] + app_accel_z[i][1])); 

    gkyl_mat_set(&rhs, 2 + i*(24), 0, qbym[i]*rhoux[i][2] + 0.5*dt*qbym[i]*rho[i][2]*(qbym[i]*ext_Ex[2] + app_accel_x[i][2])); 
    gkyl_mat_set(&rhs, 10 + i*(24), 0, qbym[i]*rhouy[i][2] + 0.5*dt*qbym[i]*rho[i][2]*(qbym[i]*ext_Ey[2] + app_accel_y[i][2])); 
    gkyl_mat_set(&rhs, 18 + i*(24), 0, qbym[i]*rhouz[i][2] + 0.5*dt*qbym[i]*rho[i][2]*(qbym[i]*ext_Ez[2] + app_accel_z[i][2])); 

    gkyl_mat_set(&rhs, 3 + i*(24), 0, qbym[i]*rhoux[i][3] + 0.5*dt*qbym[i]*rho[i][3]*(qbym[i]*ext_Ex[3] + app_accel_x[i][3])); 
    gkyl_mat_set(&rhs, 11 + i*(24), 0, qbym[i]*rhouy[i][3] + 0.5*dt*qbym[i]*rho[i][3]*(qbym[i]*ext_Ey[3] + app_accel_y[i][3])); 
    gkyl_mat_set(&rhs, 19 + i*(24), 0, qbym[i]*rhouz[i][3] + 0.5*dt*qbym[i]*rho[i][3]*(qbym[i]*ext_Ez[3] + app_accel_z[i][3])); 

    gkyl_mat_set(&rhs, 4 + i*(24), 0, qbym[i]*rhoux[i][4] + 0.5*dt*qbym[i]*rho[i][4]*(qbym[i]*ext_Ex[4] + app_accel_x[i][4])); 
    gkyl_mat_set(&rhs, 12 + i*(24), 0, qbym[i]*rhouy[i][4] + 0.5*dt*qbym[i]*rho[i][4]*(qbym[i]*ext_Ey[4] + app_accel_y[i][4])); 
    gkyl_mat_set(&rhs, 20 + i*(24), 0, qbym[i]*rhouz[i][4] + 0.5*dt*qbym[i]*rho[i][4]*(qbym[i]*ext_Ez[4] + app_accel_z[i][4])); 

    gkyl_mat_set(&rhs, 5 + i*(24), 0, qbym[i]*rhoux[i][5] + 0.5*dt*qbym[i]*rho[i][5]*(qbym[i]*ext_Ex[5] + app_accel_x[i][5])); 
    gkyl_mat_set(&rhs, 13 + i*(24), 0, qbym[i]*rhouy[i][5] + 0.5*dt*qbym[i]*rho[i][5]*(qbym[i]*ext_Ey[5] + app_accel_y[i][5])); 
    gkyl_mat_set(&rhs, 21 + i*(24), 0, qbym[i]*rhouz[i][5] + 0.5*dt*qbym[i]*rho[i][5]*(qbym[i]*ext_Ez[5] + app_accel_z[i][5])); 

    gkyl_mat_set(&rhs, 6 + i*(24), 0, qbym[i]*rhoux[i][6] + 0.5*dt*qbym[i]*rho[i][6]*(qbym[i]*ext_Ex[6] + app_accel_x[i][6])); 
    gkyl_mat_set(&rhs, 14 + i*(24), 0, qbym[i]*rhouy[i][6] + 0.5*dt*qbym[i]*rho[i][6]*(qbym[i]*ext_Ey[6] + app_accel_y[i][6])); 
    gkyl_mat_set(&rhs, 22 + i*(24), 0, qbym[i]*rhouz[i][6] + 0.5*dt*qbym[i]*rho[i][6]*(qbym[i]*ext_Ez[6] + app_accel_z[i][6])); 

    gkyl_mat_set(&rhs, 7 + i*(24), 0, qbym[i]*rhoux[i][7] + 0.5*dt*qbym[i]*rho[i][7]*(qbym[i]*ext_Ex[7] + app_accel_x[i][7])); 
    gkyl_mat_set(&rhs, 15 + i*(24), 0, qbym[i]*rhouy[i][7] + 0.5*dt*qbym[i]*rho[i][7]*(qbym[i]*ext_Ey[7] + app_accel_y[i][7])); 
    gkyl_mat_set(&rhs, 23 + i*(24), 0, qbym[i]*rhouz[i][7] + 0.5*dt*qbym[i]*rho[i][7]*(qbym[i]*ext_Ez[7] + app_accel_z[i][7])); 

  } 

  // Set RHS for Ampere's Law, including solution at known time-step and applied currents. 
  gkyl_mat_set(&rhs, 0 + num_species*(24), 0, epsilon0*Ex[0] - 0.5*dt*app_curr_x[0]); 
  gkyl_mat_set(&rhs, 8 + num_species*(24), 0, epsilon0*Ey[0] - 0.5*dt*app_curr_y[0]); 
  gkyl_mat_set(&rhs, 16 + num_species*(24), 0, epsilon0*Ez[0] - 0.5*dt*app_curr_z[0]); 

  gkyl_mat_set(&rhs, 1 + num_species*(24), 0, epsilon0*Ex[1] - 0.5*dt*app_curr_x[1]); 
  gkyl_mat_set(&rhs, 9 + num_species*(24), 0, epsilon0*Ey[1] - 0.5*dt*app_curr_y[1]); 
  gkyl_mat_set(&rhs, 17 + num_species*(24), 0, epsilon0*Ez[1] - 0.5*dt*app_curr_z[1]); 

  gkyl_mat_set(&rhs, 2 + num_species*(24), 0, epsilon0*Ex[2] - 0.5*dt*app_curr_x[2]); 
  gkyl_mat_set(&rhs, 10 + num_species*(24), 0, epsilon0*Ey[2] - 0.5*dt*app_curr_y[2]); 
  gkyl_mat_set(&rhs, 18 + num_species*(24), 0, epsilon0*Ez[2] - 0.5*dt*app_curr_z[2]); 

  gkyl_mat_set(&rhs, 3 + num_species*(24), 0, epsilon0*Ex[3] - 0.5*dt*app_curr_x[3]); 
  gkyl_mat_set(&rhs, 11 + num_species*(24), 0, epsilon0*Ey[3] - 0.5*dt*app_curr_y[3]); 
  gkyl_mat_set(&rhs, 19 + num_species*(24), 0, epsilon0*Ez[3] - 0.5*dt*app_curr_z[3]); 

  gkyl_mat_set(&rhs, 4 + num_species*(24), 0, epsilon0*Ex[4] - 0.5*dt*app_curr_x[4]); 
  gkyl_mat_set(&rhs, 12 + num_species*(24), 0, epsilon0*Ey[4] - 0.5*dt*app_curr_y[4]); 
  gkyl_mat_set(&rhs, 20 + num_species*(24), 0, epsilon0*Ez[4] - 0.5*dt*app_curr_z[4]); 

  gkyl_mat_set(&rhs, 5 + num_species*(24), 0, epsilon0*Ex[5] - 0.5*dt*app_curr_x[5]); 
  gkyl_mat_set(&rhs, 13 + num_species*(24), 0, epsilon0*Ey[5] - 0.5*dt*app_curr_y[5]); 
  gkyl_mat_set(&rhs, 21 + num_species*(24), 0, epsilon0*Ez[5] - 0.5*dt*app_curr_z[5]); 

  gkyl_mat_set(&rhs, 6 + num_species*(24), 0, epsilon0*Ex[6] - 0.5*dt*app_curr_x[6]); 
  gkyl_mat_set(&rhs, 14 + num_species*(24), 0, epsilon0*Ey[6] - 0.5*dt*app_curr_y[6]); 
  gkyl_mat_set(&rhs, 22 + num_species*(24), 0, epsilon0*Ez[6] - 0.5*dt*app_curr_z[6]); 

  gkyl_mat_set(&rhs, 7 + num_species*(24), 0, epsilon0*Ex[7] - 0.5*dt*app_curr_x[7]); 
  gkyl_mat_set(&rhs, 15 + num_species*(24), 0, epsilon0*Ey[7] - 0.5*dt*app_curr_y[7]); 
  gkyl_mat_set(&rhs, 23 + num_species*(24), 0, epsilon0*Ez[7] - 0.5*dt*app_curr_z[7]); 


  // Construct LHS. 
  // For momentum equation: J_s^{n+1} - 0.5*dt*(q_s^2/m_s^2*rho_s^n*E^{n+1} + q_s/m_s*J_s^{n+1} x B^n). 
  // For Ampere's Law: epsilon0*E^{n+1} + 0.5*dt*sum_s J_s^{n+1}. 
  for (int s = 0; s < num_species; ++s) { 
 
    double E_field_fac = -0.5*dt*qbym[s]*qbym[s]/epsilon0; 
    double B_field_fac = -0.5*dt*qbym[s]; 
    gkyl_mat_set(&lhs, 0 + s*(24), 0 + s*(24), 1.0); 
    gkyl_mat_set(&lhs, 8 + s*(24), 8 + s*(24), 1.0); 
    gkyl_mat_set(&lhs, 16 + s*(24), 16 + s*(24), 1.0); 
 
    gkyl_mat_set(&lhs, 0 + s*(24), 0 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][0])); 
    gkyl_mat_set(&lhs, 8 + s*(24), 8 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][0])); 
    gkyl_mat_set(&lhs, 16 + s*(24), 16 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][0])); 
 
    gkyl_mat_set(&lhs, 8 + s*(24), 16 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[0])); 
    gkyl_mat_set(&lhs, 16 + s*(24), 8 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[0])); 
 
    gkyl_mat_set(&lhs, 16 + s*(24), 0 + s*(24), B_field_fac*(0.3535533905932737*tot_By[0])); 
    gkyl_mat_set(&lhs, 0 + s*(24), 16 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[0])); 
 
    gkyl_mat_set(&lhs, 0 + s*(24), 8 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[0])); 
    gkyl_mat_set(&lhs, 8 + s*(24), 0 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[0])); 
 
    gkyl_mat_set(&lhs, 0 + num_species*(24), 0 + s*(24), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 8 + num_species*(24), 8 + s*(24), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 16 + num_species*(24), 16 + s*(24), 0.5*dt*(1.0)); 
 
    gkyl_mat_set(&lhs, 0 + s*(24), 1 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 8 + s*(24), 9 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 16 + s*(24), 17 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 0 + s*(24), 1 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][1])); 
    gkyl_mat_set(&lhs, 8 + s*(24), 9 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][1])); 
    gkyl_mat_set(&lhs, 16 + s*(24), 17 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][1])); 
 
    gkyl_mat_set(&lhs, 8 + s*(24), 17 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[1])); 
    gkyl_mat_set(&lhs, 16 + s*(24), 9 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[1])); 
 
    gkyl_mat_set(&lhs, 16 + s*(24), 1 + s*(24), B_field_fac*(0.3535533905932737*tot_By[1])); 
    gkyl_mat_set(&lhs, 0 + s*(24), 17 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[1])); 
 
    gkyl_mat_set(&lhs, 0 + s*(24), 9 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[1])); 
    gkyl_mat_set(&lhs, 8 + s*(24), 1 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[1])); 
 
    gkyl_mat_set(&lhs, 0 + num_species*(24), 1 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 8 + num_species*(24), 9 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 16 + num_species*(24), 17 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 0 + s*(24), 2 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 8 + s*(24), 10 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 16 + s*(24), 18 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 0 + s*(24), 2 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][2])); 
    gkyl_mat_set(&lhs, 8 + s*(24), 10 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][2])); 
    gkyl_mat_set(&lhs, 16 + s*(24), 18 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][2])); 
 
    gkyl_mat_set(&lhs, 8 + s*(24), 18 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[2])); 
    gkyl_mat_set(&lhs, 16 + s*(24), 10 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[2])); 
 
    gkyl_mat_set(&lhs, 16 + s*(24), 2 + s*(24), B_field_fac*(0.3535533905932737*tot_By[2])); 
    gkyl_mat_set(&lhs, 0 + s*(24), 18 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[2])); 
 
    gkyl_mat_set(&lhs, 0 + s*(24), 10 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[2])); 
    gkyl_mat_set(&lhs, 8 + s*(24), 2 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[2])); 
 
    gkyl_mat_set(&lhs, 0 + num_species*(24), 2 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 8 + num_species*(24), 10 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 16 + num_species*(24), 18 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 0 + s*(24), 3 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 8 + s*(24), 11 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 16 + s*(24), 19 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 0 + s*(24), 3 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][3])); 
    gkyl_mat_set(&lhs, 8 + s*(24), 11 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][3])); 
    gkyl_mat_set(&lhs, 16 + s*(24), 19 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][3])); 
 
    gkyl_mat_set(&lhs, 8 + s*(24), 19 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[3])); 
    gkyl_mat_set(&lhs, 16 + s*(24), 11 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[3])); 
 
    gkyl_mat_set(&lhs, 16 + s*(24), 3 + s*(24), B_field_fac*(0.3535533905932737*tot_By[3])); 
    gkyl_mat_set(&lhs, 0 + s*(24), 19 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[3])); 
 
    gkyl_mat_set(&lhs, 0 + s*(24), 11 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[3])); 
    gkyl_mat_set(&lhs, 8 + s*(24), 3 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[3])); 
 
    gkyl_mat_set(&lhs, 0 + num_species*(24), 3 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 8 + num_species*(24), 11 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 16 + num_species*(24), 19 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 0 + s*(24), 4 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 8 + s*(24), 12 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 16 + s*(24), 20 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 0 + s*(24), 4 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][4])); 
    gkyl_mat_set(&lhs, 8 + s*(24), 12 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][4])); 
    gkyl_mat_set(&lhs, 16 + s*(24), 20 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][4])); 
 
    gkyl_mat_set(&lhs, 8 + s*(24), 20 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[4])); 
    gkyl_mat_set(&lhs, 16 + s*(24), 12 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[4])); 
 
    gkyl_mat_set(&lhs, 16 + s*(24), 4 + s*(24), B_field_fac*(0.3535533905932737*tot_By[4])); 
    gkyl_mat_set(&lhs, 0 + s*(24), 20 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[4])); 
 
    gkyl_mat_set(&lhs, 0 + s*(24), 12 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[4])); 
    gkyl_mat_set(&lhs, 8 + s*(24), 4 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[4])); 
 
    gkyl_mat_set(&lhs, 0 + num_species*(24), 4 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 8 + num_species*(24), 12 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 16 + num_species*(24), 20 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 0 + s*(24), 5 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 8 + s*(24), 13 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 16 + s*(24), 21 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 0 + s*(24), 5 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][5])); 
    gkyl_mat_set(&lhs, 8 + s*(24), 13 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][5])); 
    gkyl_mat_set(&lhs, 16 + s*(24), 21 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][5])); 
 
    gkyl_mat_set(&lhs, 8 + s*(24), 21 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[5])); 
    gkyl_mat_set(&lhs, 16 + s*(24), 13 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[5])); 
 
    gkyl_mat_set(&lhs, 16 + s*(24), 5 + s*(24), B_field_fac*(0.3535533905932737*tot_By[5])); 
    gkyl_mat_set(&lhs, 0 + s*(24), 21 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[5])); 
 
    gkyl_mat_set(&lhs, 0 + s*(24), 13 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[5])); 
    gkyl_mat_set(&lhs, 8 + s*(24), 5 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[5])); 
 
    gkyl_mat_set(&lhs, 0 + num_species*(24), 5 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 8 + num_species*(24), 13 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 16 + num_species*(24), 21 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 0 + s*(24), 6 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 8 + s*(24), 14 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 16 + s*(24), 22 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 0 + s*(24), 6 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][6])); 
    gkyl_mat_set(&lhs, 8 + s*(24), 14 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][6])); 
    gkyl_mat_set(&lhs, 16 + s*(24), 22 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][6])); 
 
    gkyl_mat_set(&lhs, 8 + s*(24), 22 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[6])); 
    gkyl_mat_set(&lhs, 16 + s*(24), 14 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[6])); 
 
    gkyl_mat_set(&lhs, 16 + s*(24), 6 + s*(24), B_field_fac*(0.3535533905932737*tot_By[6])); 
    gkyl_mat_set(&lhs, 0 + s*(24), 22 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[6])); 
 
    gkyl_mat_set(&lhs, 0 + s*(24), 14 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[6])); 
    gkyl_mat_set(&lhs, 8 + s*(24), 6 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[6])); 
 
    gkyl_mat_set(&lhs, 0 + num_species*(24), 6 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 8 + num_species*(24), 14 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 16 + num_species*(24), 22 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 0 + s*(24), 7 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 8 + s*(24), 15 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 16 + s*(24), 23 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 0 + s*(24), 7 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][7])); 
    gkyl_mat_set(&lhs, 8 + s*(24), 15 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][7])); 
    gkyl_mat_set(&lhs, 16 + s*(24), 23 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][7])); 
 
    gkyl_mat_set(&lhs, 8 + s*(24), 23 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[7])); 
    gkyl_mat_set(&lhs, 16 + s*(24), 15 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[7])); 
 
    gkyl_mat_set(&lhs, 16 + s*(24), 7 + s*(24), B_field_fac*(0.3535533905932737*tot_By[7])); 
    gkyl_mat_set(&lhs, 0 + s*(24), 23 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[7])); 
 
    gkyl_mat_set(&lhs, 0 + s*(24), 15 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[7])); 
    gkyl_mat_set(&lhs, 8 + s*(24), 7 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[7])); 
 
    gkyl_mat_set(&lhs, 0 + num_species*(24), 7 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 8 + num_species*(24), 15 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 16 + num_species*(24), 23 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 1 + s*(24), 0 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 9 + s*(24), 8 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 17 + s*(24), 16 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 1 + s*(24), 0 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][1])); 
    gkyl_mat_set(&lhs, 9 + s*(24), 8 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][1])); 
    gkyl_mat_set(&lhs, 17 + s*(24), 16 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][1])); 
 
    gkyl_mat_set(&lhs, 9 + s*(24), 16 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[1])); 
    gkyl_mat_set(&lhs, 17 + s*(24), 8 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[1])); 
 
    gkyl_mat_set(&lhs, 17 + s*(24), 0 + s*(24), B_field_fac*(0.3535533905932737*tot_By[1])); 
    gkyl_mat_set(&lhs, 1 + s*(24), 16 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[1])); 
 
    gkyl_mat_set(&lhs, 1 + s*(24), 8 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[1])); 
    gkyl_mat_set(&lhs, 9 + s*(24), 0 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[1])); 
 
    gkyl_mat_set(&lhs, 1 + num_species*(24), 0 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 9 + num_species*(24), 8 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 17 + num_species*(24), 16 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 1 + s*(24), 1 + s*(24), 1.0); 
    gkyl_mat_set(&lhs, 9 + s*(24), 9 + s*(24), 1.0); 
    gkyl_mat_set(&lhs, 17 + s*(24), 17 + s*(24), 1.0); 
 
    gkyl_mat_set(&lhs, 1 + s*(24), 1 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][0])); 
    gkyl_mat_set(&lhs, 9 + s*(24), 9 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][0])); 
    gkyl_mat_set(&lhs, 17 + s*(24), 17 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][0])); 
 
    gkyl_mat_set(&lhs, 9 + s*(24), 17 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[0])); 
    gkyl_mat_set(&lhs, 17 + s*(24), 9 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[0])); 
 
    gkyl_mat_set(&lhs, 17 + s*(24), 1 + s*(24), B_field_fac*(0.3535533905932737*tot_By[0])); 
    gkyl_mat_set(&lhs, 1 + s*(24), 17 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[0])); 
 
    gkyl_mat_set(&lhs, 1 + s*(24), 9 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[0])); 
    gkyl_mat_set(&lhs, 9 + s*(24), 1 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[0])); 
 
    gkyl_mat_set(&lhs, 1 + num_species*(24), 1 + s*(24), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 9 + num_species*(24), 9 + s*(24), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 17 + num_species*(24), 17 + s*(24), 0.5*dt*(1.0)); 
 
    gkyl_mat_set(&lhs, 1 + s*(24), 2 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 9 + s*(24), 10 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 17 + s*(24), 18 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 1 + s*(24), 2 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][4])); 
    gkyl_mat_set(&lhs, 9 + s*(24), 10 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][4])); 
    gkyl_mat_set(&lhs, 17 + s*(24), 18 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][4])); 
 
    gkyl_mat_set(&lhs, 9 + s*(24), 18 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[4])); 
    gkyl_mat_set(&lhs, 17 + s*(24), 10 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[4])); 
 
    gkyl_mat_set(&lhs, 17 + s*(24), 2 + s*(24), B_field_fac*(0.3535533905932737*tot_By[4])); 
    gkyl_mat_set(&lhs, 1 + s*(24), 18 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[4])); 
 
    gkyl_mat_set(&lhs, 1 + s*(24), 10 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[4])); 
    gkyl_mat_set(&lhs, 9 + s*(24), 2 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[4])); 
 
    gkyl_mat_set(&lhs, 1 + num_species*(24), 2 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 9 + num_species*(24), 10 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 17 + num_species*(24), 18 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 1 + s*(24), 3 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 9 + s*(24), 11 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 17 + s*(24), 19 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 1 + s*(24), 3 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][5])); 
    gkyl_mat_set(&lhs, 9 + s*(24), 11 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][5])); 
    gkyl_mat_set(&lhs, 17 + s*(24), 19 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][5])); 
 
    gkyl_mat_set(&lhs, 9 + s*(24), 19 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[5])); 
    gkyl_mat_set(&lhs, 17 + s*(24), 11 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[5])); 
 
    gkyl_mat_set(&lhs, 17 + s*(24), 3 + s*(24), B_field_fac*(0.3535533905932737*tot_By[5])); 
    gkyl_mat_set(&lhs, 1 + s*(24), 19 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[5])); 
 
    gkyl_mat_set(&lhs, 1 + s*(24), 11 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[5])); 
    gkyl_mat_set(&lhs, 9 + s*(24), 3 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[5])); 
 
    gkyl_mat_set(&lhs, 1 + num_species*(24), 3 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 9 + num_species*(24), 11 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 17 + num_species*(24), 19 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 1 + s*(24), 4 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 9 + s*(24), 12 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 17 + s*(24), 20 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 1 + s*(24), 4 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][2])); 
    gkyl_mat_set(&lhs, 9 + s*(24), 12 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][2])); 
    gkyl_mat_set(&lhs, 17 + s*(24), 20 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][2])); 
 
    gkyl_mat_set(&lhs, 9 + s*(24), 20 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[2])); 
    gkyl_mat_set(&lhs, 17 + s*(24), 12 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[2])); 
 
    gkyl_mat_set(&lhs, 17 + s*(24), 4 + s*(24), B_field_fac*(0.3535533905932737*tot_By[2])); 
    gkyl_mat_set(&lhs, 1 + s*(24), 20 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[2])); 
 
    gkyl_mat_set(&lhs, 1 + s*(24), 12 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[2])); 
    gkyl_mat_set(&lhs, 9 + s*(24), 4 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[2])); 
 
    gkyl_mat_set(&lhs, 1 + num_species*(24), 4 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 9 + num_species*(24), 12 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 17 + num_species*(24), 20 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 1 + s*(24), 5 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 9 + s*(24), 13 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 17 + s*(24), 21 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 1 + s*(24), 5 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][3])); 
    gkyl_mat_set(&lhs, 9 + s*(24), 13 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][3])); 
    gkyl_mat_set(&lhs, 17 + s*(24), 21 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][3])); 
 
    gkyl_mat_set(&lhs, 9 + s*(24), 21 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[3])); 
    gkyl_mat_set(&lhs, 17 + s*(24), 13 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[3])); 
 
    gkyl_mat_set(&lhs, 17 + s*(24), 5 + s*(24), B_field_fac*(0.3535533905932737*tot_By[3])); 
    gkyl_mat_set(&lhs, 1 + s*(24), 21 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[3])); 
 
    gkyl_mat_set(&lhs, 1 + s*(24), 13 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[3])); 
    gkyl_mat_set(&lhs, 9 + s*(24), 5 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[3])); 
 
    gkyl_mat_set(&lhs, 1 + num_species*(24), 5 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 9 + num_species*(24), 13 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 17 + num_species*(24), 21 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 1 + s*(24), 6 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 9 + s*(24), 14 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 17 + s*(24), 22 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 1 + s*(24), 6 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][7])); 
    gkyl_mat_set(&lhs, 9 + s*(24), 14 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][7])); 
    gkyl_mat_set(&lhs, 17 + s*(24), 22 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][7])); 
 
    gkyl_mat_set(&lhs, 9 + s*(24), 22 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[7])); 
    gkyl_mat_set(&lhs, 17 + s*(24), 14 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[7])); 
 
    gkyl_mat_set(&lhs, 17 + s*(24), 6 + s*(24), B_field_fac*(0.3535533905932737*tot_By[7])); 
    gkyl_mat_set(&lhs, 1 + s*(24), 22 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[7])); 
 
    gkyl_mat_set(&lhs, 1 + s*(24), 14 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[7])); 
    gkyl_mat_set(&lhs, 9 + s*(24), 6 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[7])); 
 
    gkyl_mat_set(&lhs, 1 + num_species*(24), 6 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 9 + num_species*(24), 14 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 17 + num_species*(24), 22 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 1 + s*(24), 7 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 9 + s*(24), 15 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 17 + s*(24), 23 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 1 + s*(24), 7 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][6])); 
    gkyl_mat_set(&lhs, 9 + s*(24), 15 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][6])); 
    gkyl_mat_set(&lhs, 17 + s*(24), 23 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][6])); 
 
    gkyl_mat_set(&lhs, 9 + s*(24), 23 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[6])); 
    gkyl_mat_set(&lhs, 17 + s*(24), 15 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[6])); 
 
    gkyl_mat_set(&lhs, 17 + s*(24), 7 + s*(24), B_field_fac*(0.3535533905932737*tot_By[6])); 
    gkyl_mat_set(&lhs, 1 + s*(24), 23 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[6])); 
 
    gkyl_mat_set(&lhs, 1 + s*(24), 15 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[6])); 
    gkyl_mat_set(&lhs, 9 + s*(24), 7 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[6])); 
 
    gkyl_mat_set(&lhs, 1 + num_species*(24), 7 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 9 + num_species*(24), 15 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 17 + num_species*(24), 23 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 2 + s*(24), 0 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 10 + s*(24), 8 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 18 + s*(24), 16 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 2 + s*(24), 0 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][2])); 
    gkyl_mat_set(&lhs, 10 + s*(24), 8 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][2])); 
    gkyl_mat_set(&lhs, 18 + s*(24), 16 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][2])); 
 
    gkyl_mat_set(&lhs, 10 + s*(24), 16 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[2])); 
    gkyl_mat_set(&lhs, 18 + s*(24), 8 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[2])); 
 
    gkyl_mat_set(&lhs, 18 + s*(24), 0 + s*(24), B_field_fac*(0.3535533905932737*tot_By[2])); 
    gkyl_mat_set(&lhs, 2 + s*(24), 16 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[2])); 
 
    gkyl_mat_set(&lhs, 2 + s*(24), 8 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[2])); 
    gkyl_mat_set(&lhs, 10 + s*(24), 0 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[2])); 
 
    gkyl_mat_set(&lhs, 2 + num_species*(24), 0 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 10 + num_species*(24), 8 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 18 + num_species*(24), 16 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 2 + s*(24), 1 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 10 + s*(24), 9 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 18 + s*(24), 17 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 2 + s*(24), 1 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][4])); 
    gkyl_mat_set(&lhs, 10 + s*(24), 9 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][4])); 
    gkyl_mat_set(&lhs, 18 + s*(24), 17 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][4])); 
 
    gkyl_mat_set(&lhs, 10 + s*(24), 17 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[4])); 
    gkyl_mat_set(&lhs, 18 + s*(24), 9 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[4])); 
 
    gkyl_mat_set(&lhs, 18 + s*(24), 1 + s*(24), B_field_fac*(0.3535533905932737*tot_By[4])); 
    gkyl_mat_set(&lhs, 2 + s*(24), 17 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[4])); 
 
    gkyl_mat_set(&lhs, 2 + s*(24), 9 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[4])); 
    gkyl_mat_set(&lhs, 10 + s*(24), 1 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[4])); 
 
    gkyl_mat_set(&lhs, 2 + num_species*(24), 1 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 10 + num_species*(24), 9 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 18 + num_species*(24), 17 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 2 + s*(24), 2 + s*(24), 1.0); 
    gkyl_mat_set(&lhs, 10 + s*(24), 10 + s*(24), 1.0); 
    gkyl_mat_set(&lhs, 18 + s*(24), 18 + s*(24), 1.0); 
 
    gkyl_mat_set(&lhs, 2 + s*(24), 2 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][0])); 
    gkyl_mat_set(&lhs, 10 + s*(24), 10 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][0])); 
    gkyl_mat_set(&lhs, 18 + s*(24), 18 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][0])); 
 
    gkyl_mat_set(&lhs, 10 + s*(24), 18 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[0])); 
    gkyl_mat_set(&lhs, 18 + s*(24), 10 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[0])); 
 
    gkyl_mat_set(&lhs, 18 + s*(24), 2 + s*(24), B_field_fac*(0.3535533905932737*tot_By[0])); 
    gkyl_mat_set(&lhs, 2 + s*(24), 18 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[0])); 
 
    gkyl_mat_set(&lhs, 2 + s*(24), 10 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[0])); 
    gkyl_mat_set(&lhs, 10 + s*(24), 2 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[0])); 
 
    gkyl_mat_set(&lhs, 2 + num_species*(24), 2 + s*(24), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 10 + num_species*(24), 10 + s*(24), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 18 + num_species*(24), 18 + s*(24), 0.5*dt*(1.0)); 
 
    gkyl_mat_set(&lhs, 2 + s*(24), 3 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 10 + s*(24), 11 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 18 + s*(24), 19 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 2 + s*(24), 3 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][6])); 
    gkyl_mat_set(&lhs, 10 + s*(24), 11 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][6])); 
    gkyl_mat_set(&lhs, 18 + s*(24), 19 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][6])); 
 
    gkyl_mat_set(&lhs, 10 + s*(24), 19 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[6])); 
    gkyl_mat_set(&lhs, 18 + s*(24), 11 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[6])); 
 
    gkyl_mat_set(&lhs, 18 + s*(24), 3 + s*(24), B_field_fac*(0.3535533905932737*tot_By[6])); 
    gkyl_mat_set(&lhs, 2 + s*(24), 19 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[6])); 
 
    gkyl_mat_set(&lhs, 2 + s*(24), 11 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[6])); 
    gkyl_mat_set(&lhs, 10 + s*(24), 3 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[6])); 
 
    gkyl_mat_set(&lhs, 2 + num_species*(24), 3 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 10 + num_species*(24), 11 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 18 + num_species*(24), 19 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 2 + s*(24), 4 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 10 + s*(24), 12 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 18 + s*(24), 20 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 2 + s*(24), 4 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][1])); 
    gkyl_mat_set(&lhs, 10 + s*(24), 12 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][1])); 
    gkyl_mat_set(&lhs, 18 + s*(24), 20 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][1])); 
 
    gkyl_mat_set(&lhs, 10 + s*(24), 20 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[1])); 
    gkyl_mat_set(&lhs, 18 + s*(24), 12 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[1])); 
 
    gkyl_mat_set(&lhs, 18 + s*(24), 4 + s*(24), B_field_fac*(0.3535533905932737*tot_By[1])); 
    gkyl_mat_set(&lhs, 2 + s*(24), 20 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[1])); 
 
    gkyl_mat_set(&lhs, 2 + s*(24), 12 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[1])); 
    gkyl_mat_set(&lhs, 10 + s*(24), 4 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[1])); 
 
    gkyl_mat_set(&lhs, 2 + num_species*(24), 4 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 10 + num_species*(24), 12 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 18 + num_species*(24), 20 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 2 + s*(24), 5 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 10 + s*(24), 13 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 18 + s*(24), 21 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 2 + s*(24), 5 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][7])); 
    gkyl_mat_set(&lhs, 10 + s*(24), 13 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][7])); 
    gkyl_mat_set(&lhs, 18 + s*(24), 21 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][7])); 
 
    gkyl_mat_set(&lhs, 10 + s*(24), 21 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[7])); 
    gkyl_mat_set(&lhs, 18 + s*(24), 13 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[7])); 
 
    gkyl_mat_set(&lhs, 18 + s*(24), 5 + s*(24), B_field_fac*(0.3535533905932737*tot_By[7])); 
    gkyl_mat_set(&lhs, 2 + s*(24), 21 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[7])); 
 
    gkyl_mat_set(&lhs, 2 + s*(24), 13 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[7])); 
    gkyl_mat_set(&lhs, 10 + s*(24), 5 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[7])); 
 
    gkyl_mat_set(&lhs, 2 + num_species*(24), 5 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 10 + num_species*(24), 13 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 18 + num_species*(24), 21 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 2 + s*(24), 6 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 10 + s*(24), 14 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 18 + s*(24), 22 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 2 + s*(24), 6 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][3])); 
    gkyl_mat_set(&lhs, 10 + s*(24), 14 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][3])); 
    gkyl_mat_set(&lhs, 18 + s*(24), 22 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][3])); 
 
    gkyl_mat_set(&lhs, 10 + s*(24), 22 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[3])); 
    gkyl_mat_set(&lhs, 18 + s*(24), 14 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[3])); 
 
    gkyl_mat_set(&lhs, 18 + s*(24), 6 + s*(24), B_field_fac*(0.3535533905932737*tot_By[3])); 
    gkyl_mat_set(&lhs, 2 + s*(24), 22 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[3])); 
 
    gkyl_mat_set(&lhs, 2 + s*(24), 14 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[3])); 
    gkyl_mat_set(&lhs, 10 + s*(24), 6 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[3])); 
 
    gkyl_mat_set(&lhs, 2 + num_species*(24), 6 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 10 + num_species*(24), 14 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 18 + num_species*(24), 22 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 2 + s*(24), 7 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 10 + s*(24), 15 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 18 + s*(24), 23 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 2 + s*(24), 7 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][5])); 
    gkyl_mat_set(&lhs, 10 + s*(24), 15 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][5])); 
    gkyl_mat_set(&lhs, 18 + s*(24), 23 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][5])); 
 
    gkyl_mat_set(&lhs, 10 + s*(24), 23 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[5])); 
    gkyl_mat_set(&lhs, 18 + s*(24), 15 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[5])); 
 
    gkyl_mat_set(&lhs, 18 + s*(24), 7 + s*(24), B_field_fac*(0.3535533905932737*tot_By[5])); 
    gkyl_mat_set(&lhs, 2 + s*(24), 23 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[5])); 
 
    gkyl_mat_set(&lhs, 2 + s*(24), 15 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[5])); 
    gkyl_mat_set(&lhs, 10 + s*(24), 7 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[5])); 
 
    gkyl_mat_set(&lhs, 2 + num_species*(24), 7 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 10 + num_species*(24), 15 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 18 + num_species*(24), 23 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 3 + s*(24), 0 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 11 + s*(24), 8 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 19 + s*(24), 16 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 3 + s*(24), 0 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][3])); 
    gkyl_mat_set(&lhs, 11 + s*(24), 8 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][3])); 
    gkyl_mat_set(&lhs, 19 + s*(24), 16 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][3])); 
 
    gkyl_mat_set(&lhs, 11 + s*(24), 16 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[3])); 
    gkyl_mat_set(&lhs, 19 + s*(24), 8 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[3])); 
 
    gkyl_mat_set(&lhs, 19 + s*(24), 0 + s*(24), B_field_fac*(0.3535533905932737*tot_By[3])); 
    gkyl_mat_set(&lhs, 3 + s*(24), 16 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[3])); 
 
    gkyl_mat_set(&lhs, 3 + s*(24), 8 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[3])); 
    gkyl_mat_set(&lhs, 11 + s*(24), 0 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[3])); 
 
    gkyl_mat_set(&lhs, 3 + num_species*(24), 0 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 11 + num_species*(24), 8 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 19 + num_species*(24), 16 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 3 + s*(24), 1 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 11 + s*(24), 9 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 19 + s*(24), 17 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 3 + s*(24), 1 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][5])); 
    gkyl_mat_set(&lhs, 11 + s*(24), 9 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][5])); 
    gkyl_mat_set(&lhs, 19 + s*(24), 17 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][5])); 
 
    gkyl_mat_set(&lhs, 11 + s*(24), 17 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[5])); 
    gkyl_mat_set(&lhs, 19 + s*(24), 9 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[5])); 
 
    gkyl_mat_set(&lhs, 19 + s*(24), 1 + s*(24), B_field_fac*(0.3535533905932737*tot_By[5])); 
    gkyl_mat_set(&lhs, 3 + s*(24), 17 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[5])); 
 
    gkyl_mat_set(&lhs, 3 + s*(24), 9 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[5])); 
    gkyl_mat_set(&lhs, 11 + s*(24), 1 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[5])); 
 
    gkyl_mat_set(&lhs, 3 + num_species*(24), 1 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 11 + num_species*(24), 9 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 19 + num_species*(24), 17 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 3 + s*(24), 2 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 11 + s*(24), 10 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 19 + s*(24), 18 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 3 + s*(24), 2 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][6])); 
    gkyl_mat_set(&lhs, 11 + s*(24), 10 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][6])); 
    gkyl_mat_set(&lhs, 19 + s*(24), 18 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][6])); 
 
    gkyl_mat_set(&lhs, 11 + s*(24), 18 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[6])); 
    gkyl_mat_set(&lhs, 19 + s*(24), 10 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[6])); 
 
    gkyl_mat_set(&lhs, 19 + s*(24), 2 + s*(24), B_field_fac*(0.3535533905932737*tot_By[6])); 
    gkyl_mat_set(&lhs, 3 + s*(24), 18 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[6])); 
 
    gkyl_mat_set(&lhs, 3 + s*(24), 10 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[6])); 
    gkyl_mat_set(&lhs, 11 + s*(24), 2 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[6])); 
 
    gkyl_mat_set(&lhs, 3 + num_species*(24), 2 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 11 + num_species*(24), 10 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 19 + num_species*(24), 18 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 3 + s*(24), 3 + s*(24), 1.0); 
    gkyl_mat_set(&lhs, 11 + s*(24), 11 + s*(24), 1.0); 
    gkyl_mat_set(&lhs, 19 + s*(24), 19 + s*(24), 1.0); 
 
    gkyl_mat_set(&lhs, 3 + s*(24), 3 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][0])); 
    gkyl_mat_set(&lhs, 11 + s*(24), 11 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][0])); 
    gkyl_mat_set(&lhs, 19 + s*(24), 19 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][0])); 
 
    gkyl_mat_set(&lhs, 11 + s*(24), 19 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[0])); 
    gkyl_mat_set(&lhs, 19 + s*(24), 11 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[0])); 
 
    gkyl_mat_set(&lhs, 19 + s*(24), 3 + s*(24), B_field_fac*(0.3535533905932737*tot_By[0])); 
    gkyl_mat_set(&lhs, 3 + s*(24), 19 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[0])); 
 
    gkyl_mat_set(&lhs, 3 + s*(24), 11 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[0])); 
    gkyl_mat_set(&lhs, 11 + s*(24), 3 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[0])); 
 
    gkyl_mat_set(&lhs, 3 + num_species*(24), 3 + s*(24), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 11 + num_species*(24), 11 + s*(24), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 19 + num_species*(24), 19 + s*(24), 0.5*dt*(1.0)); 
 
    gkyl_mat_set(&lhs, 3 + s*(24), 4 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 11 + s*(24), 12 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 19 + s*(24), 20 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 3 + s*(24), 4 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][7])); 
    gkyl_mat_set(&lhs, 11 + s*(24), 12 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][7])); 
    gkyl_mat_set(&lhs, 19 + s*(24), 20 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][7])); 
 
    gkyl_mat_set(&lhs, 11 + s*(24), 20 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[7])); 
    gkyl_mat_set(&lhs, 19 + s*(24), 12 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[7])); 
 
    gkyl_mat_set(&lhs, 19 + s*(24), 4 + s*(24), B_field_fac*(0.3535533905932737*tot_By[7])); 
    gkyl_mat_set(&lhs, 3 + s*(24), 20 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[7])); 
 
    gkyl_mat_set(&lhs, 3 + s*(24), 12 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[7])); 
    gkyl_mat_set(&lhs, 11 + s*(24), 4 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[7])); 
 
    gkyl_mat_set(&lhs, 3 + num_species*(24), 4 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 11 + num_species*(24), 12 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 19 + num_species*(24), 20 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 3 + s*(24), 5 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 11 + s*(24), 13 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 19 + s*(24), 21 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 3 + s*(24), 5 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][1])); 
    gkyl_mat_set(&lhs, 11 + s*(24), 13 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][1])); 
    gkyl_mat_set(&lhs, 19 + s*(24), 21 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][1])); 
 
    gkyl_mat_set(&lhs, 11 + s*(24), 21 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[1])); 
    gkyl_mat_set(&lhs, 19 + s*(24), 13 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[1])); 
 
    gkyl_mat_set(&lhs, 19 + s*(24), 5 + s*(24), B_field_fac*(0.3535533905932737*tot_By[1])); 
    gkyl_mat_set(&lhs, 3 + s*(24), 21 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[1])); 
 
    gkyl_mat_set(&lhs, 3 + s*(24), 13 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[1])); 
    gkyl_mat_set(&lhs, 11 + s*(24), 5 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[1])); 
 
    gkyl_mat_set(&lhs, 3 + num_species*(24), 5 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 11 + num_species*(24), 13 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 19 + num_species*(24), 21 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 3 + s*(24), 6 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 11 + s*(24), 14 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 19 + s*(24), 22 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 3 + s*(24), 6 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][2])); 
    gkyl_mat_set(&lhs, 11 + s*(24), 14 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][2])); 
    gkyl_mat_set(&lhs, 19 + s*(24), 22 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][2])); 
 
    gkyl_mat_set(&lhs, 11 + s*(24), 22 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[2])); 
    gkyl_mat_set(&lhs, 19 + s*(24), 14 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[2])); 
 
    gkyl_mat_set(&lhs, 19 + s*(24), 6 + s*(24), B_field_fac*(0.3535533905932737*tot_By[2])); 
    gkyl_mat_set(&lhs, 3 + s*(24), 22 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[2])); 
 
    gkyl_mat_set(&lhs, 3 + s*(24), 14 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[2])); 
    gkyl_mat_set(&lhs, 11 + s*(24), 6 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[2])); 
 
    gkyl_mat_set(&lhs, 3 + num_species*(24), 6 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 11 + num_species*(24), 14 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 19 + num_species*(24), 22 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 3 + s*(24), 7 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 11 + s*(24), 15 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 19 + s*(24), 23 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 3 + s*(24), 7 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][4])); 
    gkyl_mat_set(&lhs, 11 + s*(24), 15 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][4])); 
    gkyl_mat_set(&lhs, 19 + s*(24), 23 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][4])); 
 
    gkyl_mat_set(&lhs, 11 + s*(24), 23 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[4])); 
    gkyl_mat_set(&lhs, 19 + s*(24), 15 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[4])); 
 
    gkyl_mat_set(&lhs, 19 + s*(24), 7 + s*(24), B_field_fac*(0.3535533905932737*tot_By[4])); 
    gkyl_mat_set(&lhs, 3 + s*(24), 23 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[4])); 
 
    gkyl_mat_set(&lhs, 3 + s*(24), 15 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[4])); 
    gkyl_mat_set(&lhs, 11 + s*(24), 7 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[4])); 
 
    gkyl_mat_set(&lhs, 3 + num_species*(24), 7 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 11 + num_species*(24), 15 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 19 + num_species*(24), 23 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 4 + s*(24), 0 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 12 + s*(24), 8 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 20 + s*(24), 16 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 4 + s*(24), 0 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][4])); 
    gkyl_mat_set(&lhs, 12 + s*(24), 8 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][4])); 
    gkyl_mat_set(&lhs, 20 + s*(24), 16 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][4])); 
 
    gkyl_mat_set(&lhs, 12 + s*(24), 16 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[4])); 
    gkyl_mat_set(&lhs, 20 + s*(24), 8 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[4])); 
 
    gkyl_mat_set(&lhs, 20 + s*(24), 0 + s*(24), B_field_fac*(0.3535533905932737*tot_By[4])); 
    gkyl_mat_set(&lhs, 4 + s*(24), 16 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[4])); 
 
    gkyl_mat_set(&lhs, 4 + s*(24), 8 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[4])); 
    gkyl_mat_set(&lhs, 12 + s*(24), 0 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[4])); 
 
    gkyl_mat_set(&lhs, 4 + num_species*(24), 0 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 12 + num_species*(24), 8 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 20 + num_species*(24), 16 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 4 + s*(24), 1 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 12 + s*(24), 9 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 20 + s*(24), 17 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 4 + s*(24), 1 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][2])); 
    gkyl_mat_set(&lhs, 12 + s*(24), 9 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][2])); 
    gkyl_mat_set(&lhs, 20 + s*(24), 17 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][2])); 
 
    gkyl_mat_set(&lhs, 12 + s*(24), 17 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[2])); 
    gkyl_mat_set(&lhs, 20 + s*(24), 9 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[2])); 
 
    gkyl_mat_set(&lhs, 20 + s*(24), 1 + s*(24), B_field_fac*(0.3535533905932737*tot_By[2])); 
    gkyl_mat_set(&lhs, 4 + s*(24), 17 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[2])); 
 
    gkyl_mat_set(&lhs, 4 + s*(24), 9 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[2])); 
    gkyl_mat_set(&lhs, 12 + s*(24), 1 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[2])); 
 
    gkyl_mat_set(&lhs, 4 + num_species*(24), 1 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 12 + num_species*(24), 9 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 20 + num_species*(24), 17 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 4 + s*(24), 2 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 12 + s*(24), 10 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 20 + s*(24), 18 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 4 + s*(24), 2 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][1])); 
    gkyl_mat_set(&lhs, 12 + s*(24), 10 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][1])); 
    gkyl_mat_set(&lhs, 20 + s*(24), 18 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][1])); 
 
    gkyl_mat_set(&lhs, 12 + s*(24), 18 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[1])); 
    gkyl_mat_set(&lhs, 20 + s*(24), 10 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[1])); 
 
    gkyl_mat_set(&lhs, 20 + s*(24), 2 + s*(24), B_field_fac*(0.3535533905932737*tot_By[1])); 
    gkyl_mat_set(&lhs, 4 + s*(24), 18 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[1])); 
 
    gkyl_mat_set(&lhs, 4 + s*(24), 10 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[1])); 
    gkyl_mat_set(&lhs, 12 + s*(24), 2 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[1])); 
 
    gkyl_mat_set(&lhs, 4 + num_species*(24), 2 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 12 + num_species*(24), 10 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 20 + num_species*(24), 18 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 4 + s*(24), 3 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 12 + s*(24), 11 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 20 + s*(24), 19 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 4 + s*(24), 3 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][7])); 
    gkyl_mat_set(&lhs, 12 + s*(24), 11 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][7])); 
    gkyl_mat_set(&lhs, 20 + s*(24), 19 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][7])); 
 
    gkyl_mat_set(&lhs, 12 + s*(24), 19 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[7])); 
    gkyl_mat_set(&lhs, 20 + s*(24), 11 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[7])); 
 
    gkyl_mat_set(&lhs, 20 + s*(24), 3 + s*(24), B_field_fac*(0.3535533905932737*tot_By[7])); 
    gkyl_mat_set(&lhs, 4 + s*(24), 19 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[7])); 
 
    gkyl_mat_set(&lhs, 4 + s*(24), 11 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[7])); 
    gkyl_mat_set(&lhs, 12 + s*(24), 3 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[7])); 
 
    gkyl_mat_set(&lhs, 4 + num_species*(24), 3 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 12 + num_species*(24), 11 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 20 + num_species*(24), 19 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 4 + s*(24), 4 + s*(24), 1.0); 
    gkyl_mat_set(&lhs, 12 + s*(24), 12 + s*(24), 1.0); 
    gkyl_mat_set(&lhs, 20 + s*(24), 20 + s*(24), 1.0); 
 
    gkyl_mat_set(&lhs, 4 + s*(24), 4 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][0])); 
    gkyl_mat_set(&lhs, 12 + s*(24), 12 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][0])); 
    gkyl_mat_set(&lhs, 20 + s*(24), 20 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][0])); 
 
    gkyl_mat_set(&lhs, 12 + s*(24), 20 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[0])); 
    gkyl_mat_set(&lhs, 20 + s*(24), 12 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[0])); 
 
    gkyl_mat_set(&lhs, 20 + s*(24), 4 + s*(24), B_field_fac*(0.3535533905932737*tot_By[0])); 
    gkyl_mat_set(&lhs, 4 + s*(24), 20 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[0])); 
 
    gkyl_mat_set(&lhs, 4 + s*(24), 12 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[0])); 
    gkyl_mat_set(&lhs, 12 + s*(24), 4 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[0])); 
 
    gkyl_mat_set(&lhs, 4 + num_species*(24), 4 + s*(24), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 12 + num_species*(24), 12 + s*(24), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 20 + num_species*(24), 20 + s*(24), 0.5*dt*(1.0)); 
 
    gkyl_mat_set(&lhs, 4 + s*(24), 5 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 12 + s*(24), 13 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 20 + s*(24), 21 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 4 + s*(24), 5 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][6])); 
    gkyl_mat_set(&lhs, 12 + s*(24), 13 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][6])); 
    gkyl_mat_set(&lhs, 20 + s*(24), 21 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][6])); 
 
    gkyl_mat_set(&lhs, 12 + s*(24), 21 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[6])); 
    gkyl_mat_set(&lhs, 20 + s*(24), 13 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[6])); 
 
    gkyl_mat_set(&lhs, 20 + s*(24), 5 + s*(24), B_field_fac*(0.3535533905932737*tot_By[6])); 
    gkyl_mat_set(&lhs, 4 + s*(24), 21 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[6])); 
 
    gkyl_mat_set(&lhs, 4 + s*(24), 13 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[6])); 
    gkyl_mat_set(&lhs, 12 + s*(24), 5 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[6])); 
 
    gkyl_mat_set(&lhs, 4 + num_species*(24), 5 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 12 + num_species*(24), 13 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 20 + num_species*(24), 21 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 4 + s*(24), 6 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 12 + s*(24), 14 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 20 + s*(24), 22 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 4 + s*(24), 6 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][5])); 
    gkyl_mat_set(&lhs, 12 + s*(24), 14 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][5])); 
    gkyl_mat_set(&lhs, 20 + s*(24), 22 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][5])); 
 
    gkyl_mat_set(&lhs, 12 + s*(24), 22 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[5])); 
    gkyl_mat_set(&lhs, 20 + s*(24), 14 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[5])); 
 
    gkyl_mat_set(&lhs, 20 + s*(24), 6 + s*(24), B_field_fac*(0.3535533905932737*tot_By[5])); 
    gkyl_mat_set(&lhs, 4 + s*(24), 22 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[5])); 
 
    gkyl_mat_set(&lhs, 4 + s*(24), 14 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[5])); 
    gkyl_mat_set(&lhs, 12 + s*(24), 6 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[5])); 
 
    gkyl_mat_set(&lhs, 4 + num_species*(24), 6 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 12 + num_species*(24), 14 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 20 + num_species*(24), 22 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 4 + s*(24), 7 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 12 + s*(24), 15 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 20 + s*(24), 23 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 4 + s*(24), 7 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][3])); 
    gkyl_mat_set(&lhs, 12 + s*(24), 15 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][3])); 
    gkyl_mat_set(&lhs, 20 + s*(24), 23 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][3])); 
 
    gkyl_mat_set(&lhs, 12 + s*(24), 23 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[3])); 
    gkyl_mat_set(&lhs, 20 + s*(24), 15 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[3])); 
 
    gkyl_mat_set(&lhs, 20 + s*(24), 7 + s*(24), B_field_fac*(0.3535533905932737*tot_By[3])); 
    gkyl_mat_set(&lhs, 4 + s*(24), 23 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[3])); 
 
    gkyl_mat_set(&lhs, 4 + s*(24), 15 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[3])); 
    gkyl_mat_set(&lhs, 12 + s*(24), 7 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[3])); 
 
    gkyl_mat_set(&lhs, 4 + num_species*(24), 7 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 12 + num_species*(24), 15 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 20 + num_species*(24), 23 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 5 + s*(24), 0 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 13 + s*(24), 8 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 21 + s*(24), 16 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 5 + s*(24), 0 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][5])); 
    gkyl_mat_set(&lhs, 13 + s*(24), 8 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][5])); 
    gkyl_mat_set(&lhs, 21 + s*(24), 16 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][5])); 
 
    gkyl_mat_set(&lhs, 13 + s*(24), 16 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[5])); 
    gkyl_mat_set(&lhs, 21 + s*(24), 8 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[5])); 
 
    gkyl_mat_set(&lhs, 21 + s*(24), 0 + s*(24), B_field_fac*(0.3535533905932737*tot_By[5])); 
    gkyl_mat_set(&lhs, 5 + s*(24), 16 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[5])); 
 
    gkyl_mat_set(&lhs, 5 + s*(24), 8 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[5])); 
    gkyl_mat_set(&lhs, 13 + s*(24), 0 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[5])); 
 
    gkyl_mat_set(&lhs, 5 + num_species*(24), 0 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 13 + num_species*(24), 8 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 21 + num_species*(24), 16 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 5 + s*(24), 1 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 13 + s*(24), 9 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 21 + s*(24), 17 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 5 + s*(24), 1 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][3])); 
    gkyl_mat_set(&lhs, 13 + s*(24), 9 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][3])); 
    gkyl_mat_set(&lhs, 21 + s*(24), 17 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][3])); 
 
    gkyl_mat_set(&lhs, 13 + s*(24), 17 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[3])); 
    gkyl_mat_set(&lhs, 21 + s*(24), 9 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[3])); 
 
    gkyl_mat_set(&lhs, 21 + s*(24), 1 + s*(24), B_field_fac*(0.3535533905932737*tot_By[3])); 
    gkyl_mat_set(&lhs, 5 + s*(24), 17 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[3])); 
 
    gkyl_mat_set(&lhs, 5 + s*(24), 9 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[3])); 
    gkyl_mat_set(&lhs, 13 + s*(24), 1 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[3])); 
 
    gkyl_mat_set(&lhs, 5 + num_species*(24), 1 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 13 + num_species*(24), 9 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 21 + num_species*(24), 17 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 5 + s*(24), 2 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 13 + s*(24), 10 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 21 + s*(24), 18 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 5 + s*(24), 2 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][7])); 
    gkyl_mat_set(&lhs, 13 + s*(24), 10 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][7])); 
    gkyl_mat_set(&lhs, 21 + s*(24), 18 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][7])); 
 
    gkyl_mat_set(&lhs, 13 + s*(24), 18 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[7])); 
    gkyl_mat_set(&lhs, 21 + s*(24), 10 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[7])); 
 
    gkyl_mat_set(&lhs, 21 + s*(24), 2 + s*(24), B_field_fac*(0.3535533905932737*tot_By[7])); 
    gkyl_mat_set(&lhs, 5 + s*(24), 18 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[7])); 
 
    gkyl_mat_set(&lhs, 5 + s*(24), 10 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[7])); 
    gkyl_mat_set(&lhs, 13 + s*(24), 2 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[7])); 
 
    gkyl_mat_set(&lhs, 5 + num_species*(24), 2 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 13 + num_species*(24), 10 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 21 + num_species*(24), 18 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 5 + s*(24), 3 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 13 + s*(24), 11 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 21 + s*(24), 19 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 5 + s*(24), 3 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][1])); 
    gkyl_mat_set(&lhs, 13 + s*(24), 11 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][1])); 
    gkyl_mat_set(&lhs, 21 + s*(24), 19 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][1])); 
 
    gkyl_mat_set(&lhs, 13 + s*(24), 19 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[1])); 
    gkyl_mat_set(&lhs, 21 + s*(24), 11 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[1])); 
 
    gkyl_mat_set(&lhs, 21 + s*(24), 3 + s*(24), B_field_fac*(0.3535533905932737*tot_By[1])); 
    gkyl_mat_set(&lhs, 5 + s*(24), 19 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[1])); 
 
    gkyl_mat_set(&lhs, 5 + s*(24), 11 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[1])); 
    gkyl_mat_set(&lhs, 13 + s*(24), 3 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[1])); 
 
    gkyl_mat_set(&lhs, 5 + num_species*(24), 3 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 13 + num_species*(24), 11 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 21 + num_species*(24), 19 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 5 + s*(24), 4 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 13 + s*(24), 12 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 21 + s*(24), 20 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 5 + s*(24), 4 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][6])); 
    gkyl_mat_set(&lhs, 13 + s*(24), 12 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][6])); 
    gkyl_mat_set(&lhs, 21 + s*(24), 20 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][6])); 
 
    gkyl_mat_set(&lhs, 13 + s*(24), 20 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[6])); 
    gkyl_mat_set(&lhs, 21 + s*(24), 12 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[6])); 
 
    gkyl_mat_set(&lhs, 21 + s*(24), 4 + s*(24), B_field_fac*(0.3535533905932737*tot_By[6])); 
    gkyl_mat_set(&lhs, 5 + s*(24), 20 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[6])); 
 
    gkyl_mat_set(&lhs, 5 + s*(24), 12 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[6])); 
    gkyl_mat_set(&lhs, 13 + s*(24), 4 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[6])); 
 
    gkyl_mat_set(&lhs, 5 + num_species*(24), 4 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 13 + num_species*(24), 12 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 21 + num_species*(24), 20 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 5 + s*(24), 5 + s*(24), 1.0); 
    gkyl_mat_set(&lhs, 13 + s*(24), 13 + s*(24), 1.0); 
    gkyl_mat_set(&lhs, 21 + s*(24), 21 + s*(24), 1.0); 
 
    gkyl_mat_set(&lhs, 5 + s*(24), 5 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][0])); 
    gkyl_mat_set(&lhs, 13 + s*(24), 13 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][0])); 
    gkyl_mat_set(&lhs, 21 + s*(24), 21 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][0])); 
 
    gkyl_mat_set(&lhs, 13 + s*(24), 21 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[0])); 
    gkyl_mat_set(&lhs, 21 + s*(24), 13 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[0])); 
 
    gkyl_mat_set(&lhs, 21 + s*(24), 5 + s*(24), B_field_fac*(0.3535533905932737*tot_By[0])); 
    gkyl_mat_set(&lhs, 5 + s*(24), 21 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[0])); 
 
    gkyl_mat_set(&lhs, 5 + s*(24), 13 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[0])); 
    gkyl_mat_set(&lhs, 13 + s*(24), 5 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[0])); 
 
    gkyl_mat_set(&lhs, 5 + num_species*(24), 5 + s*(24), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 13 + num_species*(24), 13 + s*(24), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 21 + num_species*(24), 21 + s*(24), 0.5*dt*(1.0)); 
 
    gkyl_mat_set(&lhs, 5 + s*(24), 6 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 13 + s*(24), 14 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 21 + s*(24), 22 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 5 + s*(24), 6 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][4])); 
    gkyl_mat_set(&lhs, 13 + s*(24), 14 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][4])); 
    gkyl_mat_set(&lhs, 21 + s*(24), 22 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][4])); 
 
    gkyl_mat_set(&lhs, 13 + s*(24), 22 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[4])); 
    gkyl_mat_set(&lhs, 21 + s*(24), 14 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[4])); 
 
    gkyl_mat_set(&lhs, 21 + s*(24), 6 + s*(24), B_field_fac*(0.3535533905932737*tot_By[4])); 
    gkyl_mat_set(&lhs, 5 + s*(24), 22 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[4])); 
 
    gkyl_mat_set(&lhs, 5 + s*(24), 14 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[4])); 
    gkyl_mat_set(&lhs, 13 + s*(24), 6 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[4])); 
 
    gkyl_mat_set(&lhs, 5 + num_species*(24), 6 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 13 + num_species*(24), 14 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 21 + num_species*(24), 22 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 5 + s*(24), 7 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 13 + s*(24), 15 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 21 + s*(24), 23 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 5 + s*(24), 7 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][2])); 
    gkyl_mat_set(&lhs, 13 + s*(24), 15 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][2])); 
    gkyl_mat_set(&lhs, 21 + s*(24), 23 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][2])); 
 
    gkyl_mat_set(&lhs, 13 + s*(24), 23 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[2])); 
    gkyl_mat_set(&lhs, 21 + s*(24), 15 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[2])); 
 
    gkyl_mat_set(&lhs, 21 + s*(24), 7 + s*(24), B_field_fac*(0.3535533905932737*tot_By[2])); 
    gkyl_mat_set(&lhs, 5 + s*(24), 23 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[2])); 
 
    gkyl_mat_set(&lhs, 5 + s*(24), 15 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[2])); 
    gkyl_mat_set(&lhs, 13 + s*(24), 7 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[2])); 
 
    gkyl_mat_set(&lhs, 5 + num_species*(24), 7 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 13 + num_species*(24), 15 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 21 + num_species*(24), 23 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 6 + s*(24), 0 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 14 + s*(24), 8 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 22 + s*(24), 16 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 6 + s*(24), 0 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][6])); 
    gkyl_mat_set(&lhs, 14 + s*(24), 8 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][6])); 
    gkyl_mat_set(&lhs, 22 + s*(24), 16 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][6])); 
 
    gkyl_mat_set(&lhs, 14 + s*(24), 16 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[6])); 
    gkyl_mat_set(&lhs, 22 + s*(24), 8 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[6])); 
 
    gkyl_mat_set(&lhs, 22 + s*(24), 0 + s*(24), B_field_fac*(0.3535533905932737*tot_By[6])); 
    gkyl_mat_set(&lhs, 6 + s*(24), 16 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[6])); 
 
    gkyl_mat_set(&lhs, 6 + s*(24), 8 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[6])); 
    gkyl_mat_set(&lhs, 14 + s*(24), 0 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[6])); 
 
    gkyl_mat_set(&lhs, 6 + num_species*(24), 0 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 14 + num_species*(24), 8 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 22 + num_species*(24), 16 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 6 + s*(24), 1 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 14 + s*(24), 9 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 22 + s*(24), 17 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 6 + s*(24), 1 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][7])); 
    gkyl_mat_set(&lhs, 14 + s*(24), 9 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][7])); 
    gkyl_mat_set(&lhs, 22 + s*(24), 17 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][7])); 
 
    gkyl_mat_set(&lhs, 14 + s*(24), 17 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[7])); 
    gkyl_mat_set(&lhs, 22 + s*(24), 9 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[7])); 
 
    gkyl_mat_set(&lhs, 22 + s*(24), 1 + s*(24), B_field_fac*(0.3535533905932737*tot_By[7])); 
    gkyl_mat_set(&lhs, 6 + s*(24), 17 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[7])); 
 
    gkyl_mat_set(&lhs, 6 + s*(24), 9 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[7])); 
    gkyl_mat_set(&lhs, 14 + s*(24), 1 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[7])); 
 
    gkyl_mat_set(&lhs, 6 + num_species*(24), 1 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 14 + num_species*(24), 9 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 22 + num_species*(24), 17 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 6 + s*(24), 2 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 14 + s*(24), 10 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 22 + s*(24), 18 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 6 + s*(24), 2 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][3])); 
    gkyl_mat_set(&lhs, 14 + s*(24), 10 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][3])); 
    gkyl_mat_set(&lhs, 22 + s*(24), 18 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][3])); 
 
    gkyl_mat_set(&lhs, 14 + s*(24), 18 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[3])); 
    gkyl_mat_set(&lhs, 22 + s*(24), 10 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[3])); 
 
    gkyl_mat_set(&lhs, 22 + s*(24), 2 + s*(24), B_field_fac*(0.3535533905932737*tot_By[3])); 
    gkyl_mat_set(&lhs, 6 + s*(24), 18 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[3])); 
 
    gkyl_mat_set(&lhs, 6 + s*(24), 10 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[3])); 
    gkyl_mat_set(&lhs, 14 + s*(24), 2 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[3])); 
 
    gkyl_mat_set(&lhs, 6 + num_species*(24), 2 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 14 + num_species*(24), 10 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 22 + num_species*(24), 18 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 6 + s*(24), 3 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 14 + s*(24), 11 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 22 + s*(24), 19 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 6 + s*(24), 3 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][2])); 
    gkyl_mat_set(&lhs, 14 + s*(24), 11 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][2])); 
    gkyl_mat_set(&lhs, 22 + s*(24), 19 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][2])); 
 
    gkyl_mat_set(&lhs, 14 + s*(24), 19 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[2])); 
    gkyl_mat_set(&lhs, 22 + s*(24), 11 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[2])); 
 
    gkyl_mat_set(&lhs, 22 + s*(24), 3 + s*(24), B_field_fac*(0.3535533905932737*tot_By[2])); 
    gkyl_mat_set(&lhs, 6 + s*(24), 19 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[2])); 
 
    gkyl_mat_set(&lhs, 6 + s*(24), 11 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[2])); 
    gkyl_mat_set(&lhs, 14 + s*(24), 3 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[2])); 
 
    gkyl_mat_set(&lhs, 6 + num_species*(24), 3 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 14 + num_species*(24), 11 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 22 + num_species*(24), 19 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 6 + s*(24), 4 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 14 + s*(24), 12 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 22 + s*(24), 20 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 6 + s*(24), 4 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][5])); 
    gkyl_mat_set(&lhs, 14 + s*(24), 12 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][5])); 
    gkyl_mat_set(&lhs, 22 + s*(24), 20 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][5])); 
 
    gkyl_mat_set(&lhs, 14 + s*(24), 20 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[5])); 
    gkyl_mat_set(&lhs, 22 + s*(24), 12 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[5])); 
 
    gkyl_mat_set(&lhs, 22 + s*(24), 4 + s*(24), B_field_fac*(0.3535533905932737*tot_By[5])); 
    gkyl_mat_set(&lhs, 6 + s*(24), 20 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[5])); 
 
    gkyl_mat_set(&lhs, 6 + s*(24), 12 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[5])); 
    gkyl_mat_set(&lhs, 14 + s*(24), 4 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[5])); 
 
    gkyl_mat_set(&lhs, 6 + num_species*(24), 4 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 14 + num_species*(24), 12 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 22 + num_species*(24), 20 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 6 + s*(24), 5 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 14 + s*(24), 13 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 22 + s*(24), 21 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 6 + s*(24), 5 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][4])); 
    gkyl_mat_set(&lhs, 14 + s*(24), 13 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][4])); 
    gkyl_mat_set(&lhs, 22 + s*(24), 21 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][4])); 
 
    gkyl_mat_set(&lhs, 14 + s*(24), 21 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[4])); 
    gkyl_mat_set(&lhs, 22 + s*(24), 13 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[4])); 
 
    gkyl_mat_set(&lhs, 22 + s*(24), 5 + s*(24), B_field_fac*(0.3535533905932737*tot_By[4])); 
    gkyl_mat_set(&lhs, 6 + s*(24), 21 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[4])); 
 
    gkyl_mat_set(&lhs, 6 + s*(24), 13 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[4])); 
    gkyl_mat_set(&lhs, 14 + s*(24), 5 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[4])); 
 
    gkyl_mat_set(&lhs, 6 + num_species*(24), 5 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 14 + num_species*(24), 13 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 22 + num_species*(24), 21 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 6 + s*(24), 6 + s*(24), 1.0); 
    gkyl_mat_set(&lhs, 14 + s*(24), 14 + s*(24), 1.0); 
    gkyl_mat_set(&lhs, 22 + s*(24), 22 + s*(24), 1.0); 
 
    gkyl_mat_set(&lhs, 6 + s*(24), 6 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][0])); 
    gkyl_mat_set(&lhs, 14 + s*(24), 14 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][0])); 
    gkyl_mat_set(&lhs, 22 + s*(24), 22 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][0])); 
 
    gkyl_mat_set(&lhs, 14 + s*(24), 22 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[0])); 
    gkyl_mat_set(&lhs, 22 + s*(24), 14 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[0])); 
 
    gkyl_mat_set(&lhs, 22 + s*(24), 6 + s*(24), B_field_fac*(0.3535533905932737*tot_By[0])); 
    gkyl_mat_set(&lhs, 6 + s*(24), 22 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[0])); 
 
    gkyl_mat_set(&lhs, 6 + s*(24), 14 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[0])); 
    gkyl_mat_set(&lhs, 14 + s*(24), 6 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[0])); 
 
    gkyl_mat_set(&lhs, 6 + num_species*(24), 6 + s*(24), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 14 + num_species*(24), 14 + s*(24), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 22 + num_species*(24), 22 + s*(24), 0.5*dt*(1.0)); 
 
    gkyl_mat_set(&lhs, 6 + s*(24), 7 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 14 + s*(24), 15 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 22 + s*(24), 23 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 6 + s*(24), 7 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][1])); 
    gkyl_mat_set(&lhs, 14 + s*(24), 15 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][1])); 
    gkyl_mat_set(&lhs, 22 + s*(24), 23 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][1])); 
 
    gkyl_mat_set(&lhs, 14 + s*(24), 23 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[1])); 
    gkyl_mat_set(&lhs, 22 + s*(24), 15 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[1])); 
 
    gkyl_mat_set(&lhs, 22 + s*(24), 7 + s*(24), B_field_fac*(0.3535533905932737*tot_By[1])); 
    gkyl_mat_set(&lhs, 6 + s*(24), 23 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[1])); 
 
    gkyl_mat_set(&lhs, 6 + s*(24), 15 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[1])); 
    gkyl_mat_set(&lhs, 14 + s*(24), 7 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[1])); 
 
    gkyl_mat_set(&lhs, 6 + num_species*(24), 7 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 14 + num_species*(24), 15 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 22 + num_species*(24), 23 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 7 + s*(24), 0 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 15 + s*(24), 8 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 23 + s*(24), 16 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 7 + s*(24), 0 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][7])); 
    gkyl_mat_set(&lhs, 15 + s*(24), 8 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][7])); 
    gkyl_mat_set(&lhs, 23 + s*(24), 16 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][7])); 
 
    gkyl_mat_set(&lhs, 15 + s*(24), 16 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[7])); 
    gkyl_mat_set(&lhs, 23 + s*(24), 8 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[7])); 
 
    gkyl_mat_set(&lhs, 23 + s*(24), 0 + s*(24), B_field_fac*(0.3535533905932737*tot_By[7])); 
    gkyl_mat_set(&lhs, 7 + s*(24), 16 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[7])); 
 
    gkyl_mat_set(&lhs, 7 + s*(24), 8 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[7])); 
    gkyl_mat_set(&lhs, 15 + s*(24), 0 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[7])); 
 
    gkyl_mat_set(&lhs, 7 + num_species*(24), 0 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 15 + num_species*(24), 8 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 23 + num_species*(24), 16 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 7 + s*(24), 1 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 15 + s*(24), 9 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 23 + s*(24), 17 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 7 + s*(24), 1 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][6])); 
    gkyl_mat_set(&lhs, 15 + s*(24), 9 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][6])); 
    gkyl_mat_set(&lhs, 23 + s*(24), 17 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][6])); 
 
    gkyl_mat_set(&lhs, 15 + s*(24), 17 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[6])); 
    gkyl_mat_set(&lhs, 23 + s*(24), 9 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[6])); 
 
    gkyl_mat_set(&lhs, 23 + s*(24), 1 + s*(24), B_field_fac*(0.3535533905932737*tot_By[6])); 
    gkyl_mat_set(&lhs, 7 + s*(24), 17 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[6])); 
 
    gkyl_mat_set(&lhs, 7 + s*(24), 9 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[6])); 
    gkyl_mat_set(&lhs, 15 + s*(24), 1 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[6])); 
 
    gkyl_mat_set(&lhs, 7 + num_species*(24), 1 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 15 + num_species*(24), 9 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 23 + num_species*(24), 17 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 7 + s*(24), 2 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 15 + s*(24), 10 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 23 + s*(24), 18 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 7 + s*(24), 2 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][5])); 
    gkyl_mat_set(&lhs, 15 + s*(24), 10 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][5])); 
    gkyl_mat_set(&lhs, 23 + s*(24), 18 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][5])); 
 
    gkyl_mat_set(&lhs, 15 + s*(24), 18 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[5])); 
    gkyl_mat_set(&lhs, 23 + s*(24), 10 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[5])); 
 
    gkyl_mat_set(&lhs, 23 + s*(24), 2 + s*(24), B_field_fac*(0.3535533905932737*tot_By[5])); 
    gkyl_mat_set(&lhs, 7 + s*(24), 18 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[5])); 
 
    gkyl_mat_set(&lhs, 7 + s*(24), 10 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[5])); 
    gkyl_mat_set(&lhs, 15 + s*(24), 2 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[5])); 
 
    gkyl_mat_set(&lhs, 7 + num_species*(24), 2 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 15 + num_species*(24), 10 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 23 + num_species*(24), 18 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 7 + s*(24), 3 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 15 + s*(24), 11 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 23 + s*(24), 19 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 7 + s*(24), 3 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][4])); 
    gkyl_mat_set(&lhs, 15 + s*(24), 11 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][4])); 
    gkyl_mat_set(&lhs, 23 + s*(24), 19 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][4])); 
 
    gkyl_mat_set(&lhs, 15 + s*(24), 19 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[4])); 
    gkyl_mat_set(&lhs, 23 + s*(24), 11 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[4])); 
 
    gkyl_mat_set(&lhs, 23 + s*(24), 3 + s*(24), B_field_fac*(0.3535533905932737*tot_By[4])); 
    gkyl_mat_set(&lhs, 7 + s*(24), 19 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[4])); 
 
    gkyl_mat_set(&lhs, 7 + s*(24), 11 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[4])); 
    gkyl_mat_set(&lhs, 15 + s*(24), 3 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[4])); 
 
    gkyl_mat_set(&lhs, 7 + num_species*(24), 3 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 15 + num_species*(24), 11 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 23 + num_species*(24), 19 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 7 + s*(24), 4 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 15 + s*(24), 12 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 23 + s*(24), 20 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 7 + s*(24), 4 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][3])); 
    gkyl_mat_set(&lhs, 15 + s*(24), 12 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][3])); 
    gkyl_mat_set(&lhs, 23 + s*(24), 20 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][3])); 
 
    gkyl_mat_set(&lhs, 15 + s*(24), 20 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[3])); 
    gkyl_mat_set(&lhs, 23 + s*(24), 12 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[3])); 
 
    gkyl_mat_set(&lhs, 23 + s*(24), 4 + s*(24), B_field_fac*(0.3535533905932737*tot_By[3])); 
    gkyl_mat_set(&lhs, 7 + s*(24), 20 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[3])); 
 
    gkyl_mat_set(&lhs, 7 + s*(24), 12 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[3])); 
    gkyl_mat_set(&lhs, 15 + s*(24), 4 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[3])); 
 
    gkyl_mat_set(&lhs, 7 + num_species*(24), 4 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 15 + num_species*(24), 12 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 23 + num_species*(24), 20 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 7 + s*(24), 5 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 15 + s*(24), 13 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 23 + s*(24), 21 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 7 + s*(24), 5 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][2])); 
    gkyl_mat_set(&lhs, 15 + s*(24), 13 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][2])); 
    gkyl_mat_set(&lhs, 23 + s*(24), 21 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][2])); 
 
    gkyl_mat_set(&lhs, 15 + s*(24), 21 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[2])); 
    gkyl_mat_set(&lhs, 23 + s*(24), 13 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[2])); 
 
    gkyl_mat_set(&lhs, 23 + s*(24), 5 + s*(24), B_field_fac*(0.3535533905932737*tot_By[2])); 
    gkyl_mat_set(&lhs, 7 + s*(24), 21 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[2])); 
 
    gkyl_mat_set(&lhs, 7 + s*(24), 13 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[2])); 
    gkyl_mat_set(&lhs, 15 + s*(24), 5 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[2])); 
 
    gkyl_mat_set(&lhs, 7 + num_species*(24), 5 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 15 + num_species*(24), 13 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 23 + num_species*(24), 21 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 7 + s*(24), 6 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 15 + s*(24), 14 + s*(24), 0.0); 
    gkyl_mat_set(&lhs, 23 + s*(24), 22 + s*(24), 0.0); 
 
    gkyl_mat_set(&lhs, 7 + s*(24), 6 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][1])); 
    gkyl_mat_set(&lhs, 15 + s*(24), 14 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][1])); 
    gkyl_mat_set(&lhs, 23 + s*(24), 22 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][1])); 
 
    gkyl_mat_set(&lhs, 15 + s*(24), 22 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[1])); 
    gkyl_mat_set(&lhs, 23 + s*(24), 14 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[1])); 
 
    gkyl_mat_set(&lhs, 23 + s*(24), 6 + s*(24), B_field_fac*(0.3535533905932737*tot_By[1])); 
    gkyl_mat_set(&lhs, 7 + s*(24), 22 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[1])); 
 
    gkyl_mat_set(&lhs, 7 + s*(24), 14 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[1])); 
    gkyl_mat_set(&lhs, 15 + s*(24), 6 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[1])); 
 
    gkyl_mat_set(&lhs, 7 + num_species*(24), 6 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 15 + num_species*(24), 14 + s*(24), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 23 + num_species*(24), 22 + s*(24), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 7 + s*(24), 7 + s*(24), 1.0); 
    gkyl_mat_set(&lhs, 15 + s*(24), 15 + s*(24), 1.0); 
    gkyl_mat_set(&lhs, 23 + s*(24), 23 + s*(24), 1.0); 
 
    gkyl_mat_set(&lhs, 7 + s*(24), 7 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][0])); 
    gkyl_mat_set(&lhs, 15 + s*(24), 15 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][0])); 
    gkyl_mat_set(&lhs, 23 + s*(24), 23 + num_species*(24), E_field_fac*(0.3535533905932737*rho[s][0])); 
 
    gkyl_mat_set(&lhs, 15 + s*(24), 23 + s*(24), B_field_fac*(0.3535533905932737*tot_Bx[0])); 
    gkyl_mat_set(&lhs, 23 + s*(24), 15 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bx[0])); 
 
    gkyl_mat_set(&lhs, 23 + s*(24), 7 + s*(24), B_field_fac*(0.3535533905932737*tot_By[0])); 
    gkyl_mat_set(&lhs, 7 + s*(24), 23 + s*(24), -B_field_fac*(0.3535533905932737*tot_By[0])); 
 
    gkyl_mat_set(&lhs, 7 + s*(24), 15 + s*(24), B_field_fac*(0.3535533905932737*tot_Bz[0])); 
    gkyl_mat_set(&lhs, 15 + s*(24), 7 + s*(24), -B_field_fac*(0.3535533905932737*tot_Bz[0])); 
 
    gkyl_mat_set(&lhs, 7 + num_species*(24), 7 + s*(24), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 15 + num_species*(24), 15 + s*(24), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 23 + num_species*(24), 23 + s*(24), 0.5*dt*(1.0)); 
 
  } 
  gkyl_mat_set(&lhs, 0 + num_species*(24), 0 + num_species*(24), 1.0); 
  gkyl_mat_set(&lhs, 8 + num_species*(24), 8 + num_species*(24), 1.0); 
  gkyl_mat_set(&lhs, 16 + num_species*(24), 16 + num_species*(24), 1.0); 
 
  gkyl_mat_set(&lhs, 0 + num_species*(24), 1 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 8 + num_species*(24), 9 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 16 + num_species*(24), 17 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 0 + num_species*(24), 2 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 8 + num_species*(24), 10 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 16 + num_species*(24), 18 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 0 + num_species*(24), 3 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 8 + num_species*(24), 11 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 16 + num_species*(24), 19 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 0 + num_species*(24), 4 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 8 + num_species*(24), 12 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 16 + num_species*(24), 20 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 0 + num_species*(24), 5 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 8 + num_species*(24), 13 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 16 + num_species*(24), 21 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 0 + num_species*(24), 6 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 8 + num_species*(24), 14 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 16 + num_species*(24), 22 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 0 + num_species*(24), 7 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 8 + num_species*(24), 15 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 16 + num_species*(24), 23 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 1 + num_species*(24), 0 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 9 + num_species*(24), 8 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 17 + num_species*(24), 16 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 1 + num_species*(24), 1 + num_species*(24), 1.0); 
  gkyl_mat_set(&lhs, 9 + num_species*(24), 9 + num_species*(24), 1.0); 
  gkyl_mat_set(&lhs, 17 + num_species*(24), 17 + num_species*(24), 1.0); 
 
  gkyl_mat_set(&lhs, 1 + num_species*(24), 2 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 9 + num_species*(24), 10 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 17 + num_species*(24), 18 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 1 + num_species*(24), 3 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 9 + num_species*(24), 11 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 17 + num_species*(24), 19 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 1 + num_species*(24), 4 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 9 + num_species*(24), 12 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 17 + num_species*(24), 20 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 1 + num_species*(24), 5 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 9 + num_species*(24), 13 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 17 + num_species*(24), 21 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 1 + num_species*(24), 6 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 9 + num_species*(24), 14 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 17 + num_species*(24), 22 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 1 + num_species*(24), 7 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 9 + num_species*(24), 15 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 17 + num_species*(24), 23 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 2 + num_species*(24), 0 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 10 + num_species*(24), 8 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 18 + num_species*(24), 16 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 2 + num_species*(24), 1 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 10 + num_species*(24), 9 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 18 + num_species*(24), 17 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 2 + num_species*(24), 2 + num_species*(24), 1.0); 
  gkyl_mat_set(&lhs, 10 + num_species*(24), 10 + num_species*(24), 1.0); 
  gkyl_mat_set(&lhs, 18 + num_species*(24), 18 + num_species*(24), 1.0); 
 
  gkyl_mat_set(&lhs, 2 + num_species*(24), 3 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 10 + num_species*(24), 11 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 18 + num_species*(24), 19 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 2 + num_species*(24), 4 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 10 + num_species*(24), 12 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 18 + num_species*(24), 20 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 2 + num_species*(24), 5 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 10 + num_species*(24), 13 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 18 + num_species*(24), 21 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 2 + num_species*(24), 6 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 10 + num_species*(24), 14 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 18 + num_species*(24), 22 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 2 + num_species*(24), 7 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 10 + num_species*(24), 15 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 18 + num_species*(24), 23 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 3 + num_species*(24), 0 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 11 + num_species*(24), 8 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 19 + num_species*(24), 16 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 3 + num_species*(24), 1 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 11 + num_species*(24), 9 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 19 + num_species*(24), 17 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 3 + num_species*(24), 2 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 11 + num_species*(24), 10 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 19 + num_species*(24), 18 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 3 + num_species*(24), 3 + num_species*(24), 1.0); 
  gkyl_mat_set(&lhs, 11 + num_species*(24), 11 + num_species*(24), 1.0); 
  gkyl_mat_set(&lhs, 19 + num_species*(24), 19 + num_species*(24), 1.0); 
 
  gkyl_mat_set(&lhs, 3 + num_species*(24), 4 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 11 + num_species*(24), 12 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 19 + num_species*(24), 20 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 3 + num_species*(24), 5 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 11 + num_species*(24), 13 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 19 + num_species*(24), 21 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 3 + num_species*(24), 6 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 11 + num_species*(24), 14 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 19 + num_species*(24), 22 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 3 + num_species*(24), 7 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 11 + num_species*(24), 15 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 19 + num_species*(24), 23 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 4 + num_species*(24), 0 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 12 + num_species*(24), 8 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 20 + num_species*(24), 16 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 4 + num_species*(24), 1 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 12 + num_species*(24), 9 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 20 + num_species*(24), 17 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 4 + num_species*(24), 2 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 12 + num_species*(24), 10 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 20 + num_species*(24), 18 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 4 + num_species*(24), 3 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 12 + num_species*(24), 11 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 20 + num_species*(24), 19 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 4 + num_species*(24), 4 + num_species*(24), 1.0); 
  gkyl_mat_set(&lhs, 12 + num_species*(24), 12 + num_species*(24), 1.0); 
  gkyl_mat_set(&lhs, 20 + num_species*(24), 20 + num_species*(24), 1.0); 
 
  gkyl_mat_set(&lhs, 4 + num_species*(24), 5 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 12 + num_species*(24), 13 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 20 + num_species*(24), 21 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 4 + num_species*(24), 6 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 12 + num_species*(24), 14 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 20 + num_species*(24), 22 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 4 + num_species*(24), 7 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 12 + num_species*(24), 15 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 20 + num_species*(24), 23 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 5 + num_species*(24), 0 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 13 + num_species*(24), 8 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 21 + num_species*(24), 16 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 5 + num_species*(24), 1 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 13 + num_species*(24), 9 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 21 + num_species*(24), 17 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 5 + num_species*(24), 2 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 13 + num_species*(24), 10 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 21 + num_species*(24), 18 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 5 + num_species*(24), 3 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 13 + num_species*(24), 11 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 21 + num_species*(24), 19 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 5 + num_species*(24), 4 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 13 + num_species*(24), 12 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 21 + num_species*(24), 20 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 5 + num_species*(24), 5 + num_species*(24), 1.0); 
  gkyl_mat_set(&lhs, 13 + num_species*(24), 13 + num_species*(24), 1.0); 
  gkyl_mat_set(&lhs, 21 + num_species*(24), 21 + num_species*(24), 1.0); 
 
  gkyl_mat_set(&lhs, 5 + num_species*(24), 6 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 13 + num_species*(24), 14 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 21 + num_species*(24), 22 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 5 + num_species*(24), 7 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 13 + num_species*(24), 15 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 21 + num_species*(24), 23 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 6 + num_species*(24), 0 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 14 + num_species*(24), 8 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 22 + num_species*(24), 16 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 6 + num_species*(24), 1 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 14 + num_species*(24), 9 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 22 + num_species*(24), 17 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 6 + num_species*(24), 2 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 14 + num_species*(24), 10 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 22 + num_species*(24), 18 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 6 + num_species*(24), 3 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 14 + num_species*(24), 11 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 22 + num_species*(24), 19 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 6 + num_species*(24), 4 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 14 + num_species*(24), 12 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 22 + num_species*(24), 20 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 6 + num_species*(24), 5 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 14 + num_species*(24), 13 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 22 + num_species*(24), 21 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 6 + num_species*(24), 6 + num_species*(24), 1.0); 
  gkyl_mat_set(&lhs, 14 + num_species*(24), 14 + num_species*(24), 1.0); 
  gkyl_mat_set(&lhs, 22 + num_species*(24), 22 + num_species*(24), 1.0); 
 
  gkyl_mat_set(&lhs, 6 + num_species*(24), 7 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 14 + num_species*(24), 15 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 22 + num_species*(24), 23 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 7 + num_species*(24), 0 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 15 + num_species*(24), 8 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 23 + num_species*(24), 16 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 7 + num_species*(24), 1 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 15 + num_species*(24), 9 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 23 + num_species*(24), 17 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 7 + num_species*(24), 2 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 15 + num_species*(24), 10 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 23 + num_species*(24), 18 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 7 + num_species*(24), 3 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 15 + num_species*(24), 11 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 23 + num_species*(24), 19 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 7 + num_species*(24), 4 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 15 + num_species*(24), 12 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 23 + num_species*(24), 20 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 7 + num_species*(24), 5 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 15 + num_species*(24), 13 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 23 + num_species*(24), 21 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 7 + num_species*(24), 6 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 15 + num_species*(24), 14 + num_species*(24), 0.0); 
  gkyl_mat_set(&lhs, 23 + num_species*(24), 22 + num_species*(24), 0.0); 
 
  gkyl_mat_set(&lhs, 7 + num_species*(24), 7 + num_species*(24), 1.0); 
  gkyl_mat_set(&lhs, 15 + num_species*(24), 15 + num_species*(24), 1.0); 
  gkyl_mat_set(&lhs, 23 + num_species*(24), 23 + num_species*(24), 1.0); 
 
} 
