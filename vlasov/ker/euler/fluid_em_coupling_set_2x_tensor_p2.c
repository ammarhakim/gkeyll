#include <gkyl_mat.h> 
#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void fluid_em_coupling_set_2x_tensor_p2(int count, 
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

  double rho[GKYL_MAX_SPECIES][9]; 
  double rhoux[GKYL_MAX_SPECIES][9]; 
  double rhouy[GKYL_MAX_SPECIES][9]; 
  double rhouz[GKYL_MAX_SPECIES][9]; 

  double app_accel_x[GKYL_MAX_SPECIES][9]; 
  double app_accel_y[GKYL_MAX_SPECIES][9]; 
  double app_accel_z[GKYL_MAX_SPECIES][9]; 

  for (int i = 0; i < num_species; ++i) { 
    double *inp_fluid = fluid[i]; 
    const double *inp_app_accel = app_accel[i]; 

    rho[i][0] = inp_fluid[0]; 
    rhoux[i][0] = inp_fluid[9]; 
    rhouy[i][0] = inp_fluid[18]; 
    rhouz[i][0] = inp_fluid[27]; 

    app_accel_x[i][0] = inp_app_accel[0]; 
    app_accel_y[i][0] = inp_app_accel[9]; 
    app_accel_z[i][0] = inp_app_accel[18]; 

    rho[i][1] = inp_fluid[1]; 
    rhoux[i][1] = inp_fluid[10]; 
    rhouy[i][1] = inp_fluid[19]; 
    rhouz[i][1] = inp_fluid[28]; 

    app_accel_x[i][1] = inp_app_accel[1]; 
    app_accel_y[i][1] = inp_app_accel[10]; 
    app_accel_z[i][1] = inp_app_accel[19]; 

    rho[i][2] = inp_fluid[2]; 
    rhoux[i][2] = inp_fluid[11]; 
    rhouy[i][2] = inp_fluid[20]; 
    rhouz[i][2] = inp_fluid[29]; 

    app_accel_x[i][2] = inp_app_accel[2]; 
    app_accel_y[i][2] = inp_app_accel[11]; 
    app_accel_z[i][2] = inp_app_accel[20]; 

    rho[i][3] = inp_fluid[3]; 
    rhoux[i][3] = inp_fluid[12]; 
    rhouy[i][3] = inp_fluid[21]; 
    rhouz[i][3] = inp_fluid[30]; 

    app_accel_x[i][3] = inp_app_accel[3]; 
    app_accel_y[i][3] = inp_app_accel[12]; 
    app_accel_z[i][3] = inp_app_accel[21]; 

    rho[i][4] = inp_fluid[4]; 
    rhoux[i][4] = inp_fluid[13]; 
    rhouy[i][4] = inp_fluid[22]; 
    rhouz[i][4] = inp_fluid[31]; 

    app_accel_x[i][4] = inp_app_accel[4]; 
    app_accel_y[i][4] = inp_app_accel[13]; 
    app_accel_z[i][4] = inp_app_accel[22]; 

    rho[i][5] = inp_fluid[5]; 
    rhoux[i][5] = inp_fluid[14]; 
    rhouy[i][5] = inp_fluid[23]; 
    rhouz[i][5] = inp_fluid[32]; 

    app_accel_x[i][5] = inp_app_accel[5]; 
    app_accel_y[i][5] = inp_app_accel[14]; 
    app_accel_z[i][5] = inp_app_accel[23]; 

    rho[i][6] = inp_fluid[6]; 
    rhoux[i][6] = inp_fluid[15]; 
    rhouy[i][6] = inp_fluid[24]; 
    rhouz[i][6] = inp_fluid[33]; 

    app_accel_x[i][6] = inp_app_accel[6]; 
    app_accel_y[i][6] = inp_app_accel[15]; 
    app_accel_z[i][6] = inp_app_accel[24]; 

    rho[i][7] = inp_fluid[7]; 
    rhoux[i][7] = inp_fluid[16]; 
    rhouy[i][7] = inp_fluid[25]; 
    rhouz[i][7] = inp_fluid[34]; 

    app_accel_x[i][7] = inp_app_accel[7]; 
    app_accel_y[i][7] = inp_app_accel[16]; 
    app_accel_z[i][7] = inp_app_accel[25]; 

    rho[i][8] = inp_fluid[8]; 
    rhoux[i][8] = inp_fluid[17]; 
    rhouy[i][8] = inp_fluid[26]; 
    rhouz[i][8] = inp_fluid[35]; 

    app_accel_x[i][8] = inp_app_accel[8]; 
    app_accel_y[i][8] = inp_app_accel[17]; 
    app_accel_z[i][8] = inp_app_accel[26]; 

  } 

  double *Ex = &em[0]; 
  double *Ey = &em[9]; 
  double *Ez = &em[18]; 
  double *Bx = &em[27]; 
  double *By = &em[36]; 
  double *Bz = &em[45]; 

  const double *ext_Ex = &ext_em[0]; 
  const double *ext_Ey = &ext_em[9]; 
  const double *ext_Ez = &ext_em[18]; 
  const double *ext_Bx = &ext_em[27]; 
  const double *ext_By = &ext_em[36]; 
  const double *ext_Bz = &ext_em[45]; 

  const double *app_curr_x = &app_current[0]; 
  const double *app_curr_y = &app_current[9]; 
  const double *app_curr_z = &app_current[18]; 

  double tot_Bx[9]; 
  double tot_By[9]; 
  double tot_Bz[9]; 
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
  tot_Bx[8] = Bx[8] + ext_Bx[8]; 
  tot_By[8] = By[8] + ext_By[8]; 
  tot_Bz[8] = Bz[8] + ext_Bz[8]; 

  // Set RHS for momentum equations, including solution at known time-step and external forces. 
  for (int i = 0; i < num_species; ++i) { 

    gkyl_mat_set(&rhs, 0 + i*(27), 0, qbym[i]*rhoux[i][0] + 0.5*dt*qbym[i]*rho[i][0]*(qbym[i]*ext_Ex[0] + app_accel_x[i][0])); 
    gkyl_mat_set(&rhs, 9 + i*(27), 0, qbym[i]*rhouy[i][0] + 0.5*dt*qbym[i]*rho[i][0]*(qbym[i]*ext_Ey[0] + app_accel_y[i][0])); 
    gkyl_mat_set(&rhs, 18 + i*(27), 0, qbym[i]*rhouz[i][0] + 0.5*dt*qbym[i]*rho[i][0]*(qbym[i]*ext_Ez[0] + app_accel_z[i][0])); 

    gkyl_mat_set(&rhs, 1 + i*(27), 0, qbym[i]*rhoux[i][1] + 0.5*dt*qbym[i]*rho[i][1]*(qbym[i]*ext_Ex[1] + app_accel_x[i][1])); 
    gkyl_mat_set(&rhs, 10 + i*(27), 0, qbym[i]*rhouy[i][1] + 0.5*dt*qbym[i]*rho[i][1]*(qbym[i]*ext_Ey[1] + app_accel_y[i][1])); 
    gkyl_mat_set(&rhs, 19 + i*(27), 0, qbym[i]*rhouz[i][1] + 0.5*dt*qbym[i]*rho[i][1]*(qbym[i]*ext_Ez[1] + app_accel_z[i][1])); 

    gkyl_mat_set(&rhs, 2 + i*(27), 0, qbym[i]*rhoux[i][2] + 0.5*dt*qbym[i]*rho[i][2]*(qbym[i]*ext_Ex[2] + app_accel_x[i][2])); 
    gkyl_mat_set(&rhs, 11 + i*(27), 0, qbym[i]*rhouy[i][2] + 0.5*dt*qbym[i]*rho[i][2]*(qbym[i]*ext_Ey[2] + app_accel_y[i][2])); 
    gkyl_mat_set(&rhs, 20 + i*(27), 0, qbym[i]*rhouz[i][2] + 0.5*dt*qbym[i]*rho[i][2]*(qbym[i]*ext_Ez[2] + app_accel_z[i][2])); 

    gkyl_mat_set(&rhs, 3 + i*(27), 0, qbym[i]*rhoux[i][3] + 0.5*dt*qbym[i]*rho[i][3]*(qbym[i]*ext_Ex[3] + app_accel_x[i][3])); 
    gkyl_mat_set(&rhs, 12 + i*(27), 0, qbym[i]*rhouy[i][3] + 0.5*dt*qbym[i]*rho[i][3]*(qbym[i]*ext_Ey[3] + app_accel_y[i][3])); 
    gkyl_mat_set(&rhs, 21 + i*(27), 0, qbym[i]*rhouz[i][3] + 0.5*dt*qbym[i]*rho[i][3]*(qbym[i]*ext_Ez[3] + app_accel_z[i][3])); 

    gkyl_mat_set(&rhs, 4 + i*(27), 0, qbym[i]*rhoux[i][4] + 0.5*dt*qbym[i]*rho[i][4]*(qbym[i]*ext_Ex[4] + app_accel_x[i][4])); 
    gkyl_mat_set(&rhs, 13 + i*(27), 0, qbym[i]*rhouy[i][4] + 0.5*dt*qbym[i]*rho[i][4]*(qbym[i]*ext_Ey[4] + app_accel_y[i][4])); 
    gkyl_mat_set(&rhs, 22 + i*(27), 0, qbym[i]*rhouz[i][4] + 0.5*dt*qbym[i]*rho[i][4]*(qbym[i]*ext_Ez[4] + app_accel_z[i][4])); 

    gkyl_mat_set(&rhs, 5 + i*(27), 0, qbym[i]*rhoux[i][5] + 0.5*dt*qbym[i]*rho[i][5]*(qbym[i]*ext_Ex[5] + app_accel_x[i][5])); 
    gkyl_mat_set(&rhs, 14 + i*(27), 0, qbym[i]*rhouy[i][5] + 0.5*dt*qbym[i]*rho[i][5]*(qbym[i]*ext_Ey[5] + app_accel_y[i][5])); 
    gkyl_mat_set(&rhs, 23 + i*(27), 0, qbym[i]*rhouz[i][5] + 0.5*dt*qbym[i]*rho[i][5]*(qbym[i]*ext_Ez[5] + app_accel_z[i][5])); 

    gkyl_mat_set(&rhs, 6 + i*(27), 0, qbym[i]*rhoux[i][6] + 0.5*dt*qbym[i]*rho[i][6]*(qbym[i]*ext_Ex[6] + app_accel_x[i][6])); 
    gkyl_mat_set(&rhs, 15 + i*(27), 0, qbym[i]*rhouy[i][6] + 0.5*dt*qbym[i]*rho[i][6]*(qbym[i]*ext_Ey[6] + app_accel_y[i][6])); 
    gkyl_mat_set(&rhs, 24 + i*(27), 0, qbym[i]*rhouz[i][6] + 0.5*dt*qbym[i]*rho[i][6]*(qbym[i]*ext_Ez[6] + app_accel_z[i][6])); 

    gkyl_mat_set(&rhs, 7 + i*(27), 0, qbym[i]*rhoux[i][7] + 0.5*dt*qbym[i]*rho[i][7]*(qbym[i]*ext_Ex[7] + app_accel_x[i][7])); 
    gkyl_mat_set(&rhs, 16 + i*(27), 0, qbym[i]*rhouy[i][7] + 0.5*dt*qbym[i]*rho[i][7]*(qbym[i]*ext_Ey[7] + app_accel_y[i][7])); 
    gkyl_mat_set(&rhs, 25 + i*(27), 0, qbym[i]*rhouz[i][7] + 0.5*dt*qbym[i]*rho[i][7]*(qbym[i]*ext_Ez[7] + app_accel_z[i][7])); 

    gkyl_mat_set(&rhs, 8 + i*(27), 0, qbym[i]*rhoux[i][8] + 0.5*dt*qbym[i]*rho[i][8]*(qbym[i]*ext_Ex[8] + app_accel_x[i][8])); 
    gkyl_mat_set(&rhs, 17 + i*(27), 0, qbym[i]*rhouy[i][8] + 0.5*dt*qbym[i]*rho[i][8]*(qbym[i]*ext_Ey[8] + app_accel_y[i][8])); 
    gkyl_mat_set(&rhs, 26 + i*(27), 0, qbym[i]*rhouz[i][8] + 0.5*dt*qbym[i]*rho[i][8]*(qbym[i]*ext_Ez[8] + app_accel_z[i][8])); 

  } 

  // Set RHS for Ampere's Law, including solution at known time-step and applied currents. 
  gkyl_mat_set(&rhs, 0 + num_species*(27), 0, epsilon0*Ex[0] - 0.5*dt*app_curr_x[0]); 
  gkyl_mat_set(&rhs, 9 + num_species*(27), 0, epsilon0*Ey[0] - 0.5*dt*app_curr_y[0]); 
  gkyl_mat_set(&rhs, 18 + num_species*(27), 0, epsilon0*Ez[0] - 0.5*dt*app_curr_z[0]); 

  gkyl_mat_set(&rhs, 1 + num_species*(27), 0, epsilon0*Ex[1] - 0.5*dt*app_curr_x[1]); 
  gkyl_mat_set(&rhs, 10 + num_species*(27), 0, epsilon0*Ey[1] - 0.5*dt*app_curr_y[1]); 
  gkyl_mat_set(&rhs, 19 + num_species*(27), 0, epsilon0*Ez[1] - 0.5*dt*app_curr_z[1]); 

  gkyl_mat_set(&rhs, 2 + num_species*(27), 0, epsilon0*Ex[2] - 0.5*dt*app_curr_x[2]); 
  gkyl_mat_set(&rhs, 11 + num_species*(27), 0, epsilon0*Ey[2] - 0.5*dt*app_curr_y[2]); 
  gkyl_mat_set(&rhs, 20 + num_species*(27), 0, epsilon0*Ez[2] - 0.5*dt*app_curr_z[2]); 

  gkyl_mat_set(&rhs, 3 + num_species*(27), 0, epsilon0*Ex[3] - 0.5*dt*app_curr_x[3]); 
  gkyl_mat_set(&rhs, 12 + num_species*(27), 0, epsilon0*Ey[3] - 0.5*dt*app_curr_y[3]); 
  gkyl_mat_set(&rhs, 21 + num_species*(27), 0, epsilon0*Ez[3] - 0.5*dt*app_curr_z[3]); 

  gkyl_mat_set(&rhs, 4 + num_species*(27), 0, epsilon0*Ex[4] - 0.5*dt*app_curr_x[4]); 
  gkyl_mat_set(&rhs, 13 + num_species*(27), 0, epsilon0*Ey[4] - 0.5*dt*app_curr_y[4]); 
  gkyl_mat_set(&rhs, 22 + num_species*(27), 0, epsilon0*Ez[4] - 0.5*dt*app_curr_z[4]); 

  gkyl_mat_set(&rhs, 5 + num_species*(27), 0, epsilon0*Ex[5] - 0.5*dt*app_curr_x[5]); 
  gkyl_mat_set(&rhs, 14 + num_species*(27), 0, epsilon0*Ey[5] - 0.5*dt*app_curr_y[5]); 
  gkyl_mat_set(&rhs, 23 + num_species*(27), 0, epsilon0*Ez[5] - 0.5*dt*app_curr_z[5]); 

  gkyl_mat_set(&rhs, 6 + num_species*(27), 0, epsilon0*Ex[6] - 0.5*dt*app_curr_x[6]); 
  gkyl_mat_set(&rhs, 15 + num_species*(27), 0, epsilon0*Ey[6] - 0.5*dt*app_curr_y[6]); 
  gkyl_mat_set(&rhs, 24 + num_species*(27), 0, epsilon0*Ez[6] - 0.5*dt*app_curr_z[6]); 

  gkyl_mat_set(&rhs, 7 + num_species*(27), 0, epsilon0*Ex[7] - 0.5*dt*app_curr_x[7]); 
  gkyl_mat_set(&rhs, 16 + num_species*(27), 0, epsilon0*Ey[7] - 0.5*dt*app_curr_y[7]); 
  gkyl_mat_set(&rhs, 25 + num_species*(27), 0, epsilon0*Ez[7] - 0.5*dt*app_curr_z[7]); 

  gkyl_mat_set(&rhs, 8 + num_species*(27), 0, epsilon0*Ex[8] - 0.5*dt*app_curr_x[8]); 
  gkyl_mat_set(&rhs, 17 + num_species*(27), 0, epsilon0*Ey[8] - 0.5*dt*app_curr_y[8]); 
  gkyl_mat_set(&rhs, 26 + num_species*(27), 0, epsilon0*Ez[8] - 0.5*dt*app_curr_z[8]); 


  // Construct LHS. 
  // For momentum equation: J_s^{n+1} - 0.5*dt*(q_s^2/m_s^2*rho_s^n*E^{n+1} + q_s/m_s*J_s^{n+1} x B^n). 
  // For Ampere's Law: epsilon0*E^{n+1} + 0.5*dt*sum_s J_s^{n+1}. 
  for (int s = 0; s < num_species; ++s) { 
 
    double E_field_fac = -0.5*dt*qbym[s]*qbym[s]/epsilon0; 
    double B_field_fac = -0.5*dt*qbym[s]; 
    gkyl_mat_set(&lhs, 0 + s*(27), 0 + s*(27), 1.0); 
    gkyl_mat_set(&lhs, 9 + s*(27), 9 + s*(27), 1.0); 
    gkyl_mat_set(&lhs, 18 + s*(27), 18 + s*(27), 1.0); 
 
    gkyl_mat_set(&lhs, 0 + s*(27), 0 + num_species*(27), E_field_fac*(0.5*rho[s][0])); 
    gkyl_mat_set(&lhs, 9 + s*(27), 9 + num_species*(27), E_field_fac*(0.5*rho[s][0])); 
    gkyl_mat_set(&lhs, 18 + s*(27), 18 + num_species*(27), E_field_fac*(0.5*rho[s][0])); 
 
    gkyl_mat_set(&lhs, 9 + s*(27), 18 + s*(27), B_field_fac*(0.5*tot_Bx[0])); 
    gkyl_mat_set(&lhs, 18 + s*(27), 9 + s*(27), -B_field_fac*(0.5*tot_Bx[0])); 
 
    gkyl_mat_set(&lhs, 18 + s*(27), 0 + s*(27), B_field_fac*(0.5*tot_By[0])); 
    gkyl_mat_set(&lhs, 0 + s*(27), 18 + s*(27), -B_field_fac*(0.5*tot_By[0])); 
 
    gkyl_mat_set(&lhs, 0 + s*(27), 9 + s*(27), B_field_fac*(0.5*tot_Bz[0])); 
    gkyl_mat_set(&lhs, 9 + s*(27), 0 + s*(27), -B_field_fac*(0.5*tot_Bz[0])); 
 
    gkyl_mat_set(&lhs, 0 + num_species*(27), 0 + s*(27), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 9 + num_species*(27), 9 + s*(27), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 18 + num_species*(27), 18 + s*(27), 0.5*dt*(1.0)); 
 
    gkyl_mat_set(&lhs, 0 + s*(27), 1 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 9 + s*(27), 10 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 18 + s*(27), 19 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 0 + s*(27), 1 + num_species*(27), E_field_fac*(0.5*rho[s][1])); 
    gkyl_mat_set(&lhs, 9 + s*(27), 10 + num_species*(27), E_field_fac*(0.5*rho[s][1])); 
    gkyl_mat_set(&lhs, 18 + s*(27), 19 + num_species*(27), E_field_fac*(0.5*rho[s][1])); 
 
    gkyl_mat_set(&lhs, 9 + s*(27), 19 + s*(27), B_field_fac*(0.5*tot_Bx[1])); 
    gkyl_mat_set(&lhs, 18 + s*(27), 10 + s*(27), -B_field_fac*(0.5*tot_Bx[1])); 
 
    gkyl_mat_set(&lhs, 18 + s*(27), 1 + s*(27), B_field_fac*(0.5*tot_By[1])); 
    gkyl_mat_set(&lhs, 0 + s*(27), 19 + s*(27), -B_field_fac*(0.5*tot_By[1])); 
 
    gkyl_mat_set(&lhs, 0 + s*(27), 10 + s*(27), B_field_fac*(0.5*tot_Bz[1])); 
    gkyl_mat_set(&lhs, 9 + s*(27), 1 + s*(27), -B_field_fac*(0.5*tot_Bz[1])); 
 
    gkyl_mat_set(&lhs, 0 + num_species*(27), 1 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 9 + num_species*(27), 10 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 18 + num_species*(27), 19 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 0 + s*(27), 2 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 9 + s*(27), 11 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 18 + s*(27), 20 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 0 + s*(27), 2 + num_species*(27), E_field_fac*(0.5*rho[s][2])); 
    gkyl_mat_set(&lhs, 9 + s*(27), 11 + num_species*(27), E_field_fac*(0.5*rho[s][2])); 
    gkyl_mat_set(&lhs, 18 + s*(27), 20 + num_species*(27), E_field_fac*(0.5*rho[s][2])); 
 
    gkyl_mat_set(&lhs, 9 + s*(27), 20 + s*(27), B_field_fac*(0.5*tot_Bx[2])); 
    gkyl_mat_set(&lhs, 18 + s*(27), 11 + s*(27), -B_field_fac*(0.5*tot_Bx[2])); 
 
    gkyl_mat_set(&lhs, 18 + s*(27), 2 + s*(27), B_field_fac*(0.5*tot_By[2])); 
    gkyl_mat_set(&lhs, 0 + s*(27), 20 + s*(27), -B_field_fac*(0.5*tot_By[2])); 
 
    gkyl_mat_set(&lhs, 0 + s*(27), 11 + s*(27), B_field_fac*(0.5*tot_Bz[2])); 
    gkyl_mat_set(&lhs, 9 + s*(27), 2 + s*(27), -B_field_fac*(0.5*tot_Bz[2])); 
 
    gkyl_mat_set(&lhs, 0 + num_species*(27), 2 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 9 + num_species*(27), 11 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 18 + num_species*(27), 20 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 0 + s*(27), 3 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 9 + s*(27), 12 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 18 + s*(27), 21 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 0 + s*(27), 3 + num_species*(27), E_field_fac*(0.5*rho[s][3])); 
    gkyl_mat_set(&lhs, 9 + s*(27), 12 + num_species*(27), E_field_fac*(0.5*rho[s][3])); 
    gkyl_mat_set(&lhs, 18 + s*(27), 21 + num_species*(27), E_field_fac*(0.5*rho[s][3])); 
 
    gkyl_mat_set(&lhs, 9 + s*(27), 21 + s*(27), B_field_fac*(0.5*tot_Bx[3])); 
    gkyl_mat_set(&lhs, 18 + s*(27), 12 + s*(27), -B_field_fac*(0.5*tot_Bx[3])); 
 
    gkyl_mat_set(&lhs, 18 + s*(27), 3 + s*(27), B_field_fac*(0.5*tot_By[3])); 
    gkyl_mat_set(&lhs, 0 + s*(27), 21 + s*(27), -B_field_fac*(0.5*tot_By[3])); 
 
    gkyl_mat_set(&lhs, 0 + s*(27), 12 + s*(27), B_field_fac*(0.5*tot_Bz[3])); 
    gkyl_mat_set(&lhs, 9 + s*(27), 3 + s*(27), -B_field_fac*(0.5*tot_Bz[3])); 
 
    gkyl_mat_set(&lhs, 0 + num_species*(27), 3 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 9 + num_species*(27), 12 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 18 + num_species*(27), 21 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 0 + s*(27), 4 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 9 + s*(27), 13 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 18 + s*(27), 22 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 0 + s*(27), 4 + num_species*(27), E_field_fac*(0.5*rho[s][4])); 
    gkyl_mat_set(&lhs, 9 + s*(27), 13 + num_species*(27), E_field_fac*(0.5*rho[s][4])); 
    gkyl_mat_set(&lhs, 18 + s*(27), 22 + num_species*(27), E_field_fac*(0.5*rho[s][4])); 
 
    gkyl_mat_set(&lhs, 9 + s*(27), 22 + s*(27), B_field_fac*(0.5*tot_Bx[4])); 
    gkyl_mat_set(&lhs, 18 + s*(27), 13 + s*(27), -B_field_fac*(0.5*tot_Bx[4])); 
 
    gkyl_mat_set(&lhs, 18 + s*(27), 4 + s*(27), B_field_fac*(0.5*tot_By[4])); 
    gkyl_mat_set(&lhs, 0 + s*(27), 22 + s*(27), -B_field_fac*(0.5*tot_By[4])); 
 
    gkyl_mat_set(&lhs, 0 + s*(27), 13 + s*(27), B_field_fac*(0.5*tot_Bz[4])); 
    gkyl_mat_set(&lhs, 9 + s*(27), 4 + s*(27), -B_field_fac*(0.5*tot_Bz[4])); 
 
    gkyl_mat_set(&lhs, 0 + num_species*(27), 4 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 9 + num_species*(27), 13 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 18 + num_species*(27), 22 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 0 + s*(27), 5 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 9 + s*(27), 14 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 18 + s*(27), 23 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 0 + s*(27), 5 + num_species*(27), E_field_fac*(0.5*rho[s][5])); 
    gkyl_mat_set(&lhs, 9 + s*(27), 14 + num_species*(27), E_field_fac*(0.5*rho[s][5])); 
    gkyl_mat_set(&lhs, 18 + s*(27), 23 + num_species*(27), E_field_fac*(0.5*rho[s][5])); 
 
    gkyl_mat_set(&lhs, 9 + s*(27), 23 + s*(27), B_field_fac*(0.5*tot_Bx[5])); 
    gkyl_mat_set(&lhs, 18 + s*(27), 14 + s*(27), -B_field_fac*(0.5*tot_Bx[5])); 
 
    gkyl_mat_set(&lhs, 18 + s*(27), 5 + s*(27), B_field_fac*(0.5*tot_By[5])); 
    gkyl_mat_set(&lhs, 0 + s*(27), 23 + s*(27), -B_field_fac*(0.5*tot_By[5])); 
 
    gkyl_mat_set(&lhs, 0 + s*(27), 14 + s*(27), B_field_fac*(0.5*tot_Bz[5])); 
    gkyl_mat_set(&lhs, 9 + s*(27), 5 + s*(27), -B_field_fac*(0.5*tot_Bz[5])); 
 
    gkyl_mat_set(&lhs, 0 + num_species*(27), 5 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 9 + num_species*(27), 14 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 18 + num_species*(27), 23 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 0 + s*(27), 6 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 9 + s*(27), 15 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 18 + s*(27), 24 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 0 + s*(27), 6 + num_species*(27), E_field_fac*(0.5*rho[s][6])); 
    gkyl_mat_set(&lhs, 9 + s*(27), 15 + num_species*(27), E_field_fac*(0.5*rho[s][6])); 
    gkyl_mat_set(&lhs, 18 + s*(27), 24 + num_species*(27), E_field_fac*(0.5*rho[s][6])); 
 
    gkyl_mat_set(&lhs, 9 + s*(27), 24 + s*(27), B_field_fac*(0.5*tot_Bx[6])); 
    gkyl_mat_set(&lhs, 18 + s*(27), 15 + s*(27), -B_field_fac*(0.5*tot_Bx[6])); 
 
    gkyl_mat_set(&lhs, 18 + s*(27), 6 + s*(27), B_field_fac*(0.5*tot_By[6])); 
    gkyl_mat_set(&lhs, 0 + s*(27), 24 + s*(27), -B_field_fac*(0.5*tot_By[6])); 
 
    gkyl_mat_set(&lhs, 0 + s*(27), 15 + s*(27), B_field_fac*(0.5*tot_Bz[6])); 
    gkyl_mat_set(&lhs, 9 + s*(27), 6 + s*(27), -B_field_fac*(0.5*tot_Bz[6])); 
 
    gkyl_mat_set(&lhs, 0 + num_species*(27), 6 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 9 + num_species*(27), 15 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 18 + num_species*(27), 24 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 0 + s*(27), 7 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 9 + s*(27), 16 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 18 + s*(27), 25 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 0 + s*(27), 7 + num_species*(27), E_field_fac*(0.5*rho[s][7])); 
    gkyl_mat_set(&lhs, 9 + s*(27), 16 + num_species*(27), E_field_fac*(0.5*rho[s][7])); 
    gkyl_mat_set(&lhs, 18 + s*(27), 25 + num_species*(27), E_field_fac*(0.5*rho[s][7])); 
 
    gkyl_mat_set(&lhs, 9 + s*(27), 25 + s*(27), B_field_fac*(0.5*tot_Bx[7])); 
    gkyl_mat_set(&lhs, 18 + s*(27), 16 + s*(27), -B_field_fac*(0.5*tot_Bx[7])); 
 
    gkyl_mat_set(&lhs, 18 + s*(27), 7 + s*(27), B_field_fac*(0.5*tot_By[7])); 
    gkyl_mat_set(&lhs, 0 + s*(27), 25 + s*(27), -B_field_fac*(0.5*tot_By[7])); 
 
    gkyl_mat_set(&lhs, 0 + s*(27), 16 + s*(27), B_field_fac*(0.5*tot_Bz[7])); 
    gkyl_mat_set(&lhs, 9 + s*(27), 7 + s*(27), -B_field_fac*(0.5*tot_Bz[7])); 
 
    gkyl_mat_set(&lhs, 0 + num_species*(27), 7 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 9 + num_species*(27), 16 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 18 + num_species*(27), 25 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 0 + s*(27), 8 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 9 + s*(27), 17 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 18 + s*(27), 26 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 0 + s*(27), 8 + num_species*(27), E_field_fac*(0.5*rho[s][8])); 
    gkyl_mat_set(&lhs, 9 + s*(27), 17 + num_species*(27), E_field_fac*(0.5*rho[s][8])); 
    gkyl_mat_set(&lhs, 18 + s*(27), 26 + num_species*(27), E_field_fac*(0.5*rho[s][8])); 
 
    gkyl_mat_set(&lhs, 9 + s*(27), 26 + s*(27), B_field_fac*(0.5*tot_Bx[8])); 
    gkyl_mat_set(&lhs, 18 + s*(27), 17 + s*(27), -B_field_fac*(0.5*tot_Bx[8])); 
 
    gkyl_mat_set(&lhs, 18 + s*(27), 8 + s*(27), B_field_fac*(0.5*tot_By[8])); 
    gkyl_mat_set(&lhs, 0 + s*(27), 26 + s*(27), -B_field_fac*(0.5*tot_By[8])); 
 
    gkyl_mat_set(&lhs, 0 + s*(27), 17 + s*(27), B_field_fac*(0.5*tot_Bz[8])); 
    gkyl_mat_set(&lhs, 9 + s*(27), 8 + s*(27), -B_field_fac*(0.5*tot_Bz[8])); 
 
    gkyl_mat_set(&lhs, 0 + num_species*(27), 8 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 9 + num_species*(27), 17 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 18 + num_species*(27), 26 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 1 + s*(27), 0 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 10 + s*(27), 9 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 19 + s*(27), 18 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 1 + s*(27), 0 + num_species*(27), E_field_fac*(0.5*rho[s][1])); 
    gkyl_mat_set(&lhs, 10 + s*(27), 9 + num_species*(27), E_field_fac*(0.5*rho[s][1])); 
    gkyl_mat_set(&lhs, 19 + s*(27), 18 + num_species*(27), E_field_fac*(0.5*rho[s][1])); 
 
    gkyl_mat_set(&lhs, 10 + s*(27), 18 + s*(27), B_field_fac*(0.5*tot_Bx[1])); 
    gkyl_mat_set(&lhs, 19 + s*(27), 9 + s*(27), -B_field_fac*(0.5*tot_Bx[1])); 
 
    gkyl_mat_set(&lhs, 19 + s*(27), 0 + s*(27), B_field_fac*(0.5*tot_By[1])); 
    gkyl_mat_set(&lhs, 1 + s*(27), 18 + s*(27), -B_field_fac*(0.5*tot_By[1])); 
 
    gkyl_mat_set(&lhs, 1 + s*(27), 9 + s*(27), B_field_fac*(0.5*tot_Bz[1])); 
    gkyl_mat_set(&lhs, 10 + s*(27), 0 + s*(27), -B_field_fac*(0.5*tot_Bz[1])); 
 
    gkyl_mat_set(&lhs, 1 + num_species*(27), 0 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 10 + num_species*(27), 9 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 19 + num_species*(27), 18 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 1 + s*(27), 1 + s*(27), 1.0); 
    gkyl_mat_set(&lhs, 10 + s*(27), 10 + s*(27), 1.0); 
    gkyl_mat_set(&lhs, 19 + s*(27), 19 + s*(27), 1.0); 
 
    gkyl_mat_set(&lhs, 1 + s*(27), 1 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][4]+0.5*rho[s][0])); 
    gkyl_mat_set(&lhs, 10 + s*(27), 10 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][4]+0.5*rho[s][0])); 
    gkyl_mat_set(&lhs, 19 + s*(27), 19 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][4]+0.5*rho[s][0])); 
 
    gkyl_mat_set(&lhs, 10 + s*(27), 19 + s*(27), B_field_fac*(0.4472135954999579*tot_Bx[4]+0.5*tot_Bx[0])); 
    gkyl_mat_set(&lhs, 19 + s*(27), 10 + s*(27), -B_field_fac*(0.4472135954999579*tot_Bx[4]+0.5*tot_Bx[0])); 
 
    gkyl_mat_set(&lhs, 19 + s*(27), 1 + s*(27), B_field_fac*(0.4472135954999579*tot_By[4]+0.5*tot_By[0])); 
    gkyl_mat_set(&lhs, 1 + s*(27), 19 + s*(27), -B_field_fac*(0.4472135954999579*tot_By[4]+0.5*tot_By[0])); 
 
    gkyl_mat_set(&lhs, 1 + s*(27), 10 + s*(27), B_field_fac*(0.4472135954999579*tot_Bz[4]+0.5*tot_Bz[0])); 
    gkyl_mat_set(&lhs, 10 + s*(27), 1 + s*(27), -B_field_fac*(0.4472135954999579*tot_Bz[4]+0.5*tot_Bz[0])); 
 
    gkyl_mat_set(&lhs, 1 + num_species*(27), 1 + s*(27), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 10 + num_species*(27), 10 + s*(27), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 19 + num_species*(27), 19 + s*(27), 0.5*dt*(1.0)); 
 
    gkyl_mat_set(&lhs, 1 + s*(27), 2 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 10 + s*(27), 11 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 19 + s*(27), 20 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 1 + s*(27), 2 + num_species*(27), E_field_fac*(0.5*rho[s][3])); 
    gkyl_mat_set(&lhs, 10 + s*(27), 11 + num_species*(27), E_field_fac*(0.5*rho[s][3])); 
    gkyl_mat_set(&lhs, 19 + s*(27), 20 + num_species*(27), E_field_fac*(0.5*rho[s][3])); 
 
    gkyl_mat_set(&lhs, 10 + s*(27), 20 + s*(27), B_field_fac*(0.5*tot_Bx[3])); 
    gkyl_mat_set(&lhs, 19 + s*(27), 11 + s*(27), -B_field_fac*(0.5*tot_Bx[3])); 
 
    gkyl_mat_set(&lhs, 19 + s*(27), 2 + s*(27), B_field_fac*(0.5*tot_By[3])); 
    gkyl_mat_set(&lhs, 1 + s*(27), 20 + s*(27), -B_field_fac*(0.5*tot_By[3])); 
 
    gkyl_mat_set(&lhs, 1 + s*(27), 11 + s*(27), B_field_fac*(0.5*tot_Bz[3])); 
    gkyl_mat_set(&lhs, 10 + s*(27), 2 + s*(27), -B_field_fac*(0.5*tot_Bz[3])); 
 
    gkyl_mat_set(&lhs, 1 + num_species*(27), 2 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 10 + num_species*(27), 11 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 19 + num_species*(27), 20 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 1 + s*(27), 3 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 10 + s*(27), 12 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 19 + s*(27), 21 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 1 + s*(27), 3 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][6]+0.5*rho[s][2])); 
    gkyl_mat_set(&lhs, 10 + s*(27), 12 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][6]+0.5*rho[s][2])); 
    gkyl_mat_set(&lhs, 19 + s*(27), 21 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][6]+0.5*rho[s][2])); 
 
    gkyl_mat_set(&lhs, 10 + s*(27), 21 + s*(27), B_field_fac*(0.447213595499958*tot_Bx[6]+0.5*tot_Bx[2])); 
    gkyl_mat_set(&lhs, 19 + s*(27), 12 + s*(27), -B_field_fac*(0.447213595499958*tot_Bx[6]+0.5*tot_Bx[2])); 
 
    gkyl_mat_set(&lhs, 19 + s*(27), 3 + s*(27), B_field_fac*(0.447213595499958*tot_By[6]+0.5*tot_By[2])); 
    gkyl_mat_set(&lhs, 1 + s*(27), 21 + s*(27), -B_field_fac*(0.447213595499958*tot_By[6]+0.5*tot_By[2])); 
 
    gkyl_mat_set(&lhs, 1 + s*(27), 12 + s*(27), B_field_fac*(0.447213595499958*tot_Bz[6]+0.5*tot_Bz[2])); 
    gkyl_mat_set(&lhs, 10 + s*(27), 3 + s*(27), -B_field_fac*(0.447213595499958*tot_Bz[6]+0.5*tot_Bz[2])); 
 
    gkyl_mat_set(&lhs, 1 + num_species*(27), 3 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 10 + num_species*(27), 12 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 19 + num_species*(27), 21 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 1 + s*(27), 4 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 10 + s*(27), 13 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 19 + s*(27), 22 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 1 + s*(27), 4 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][1])); 
    gkyl_mat_set(&lhs, 10 + s*(27), 13 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][1])); 
    gkyl_mat_set(&lhs, 19 + s*(27), 22 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][1])); 
 
    gkyl_mat_set(&lhs, 10 + s*(27), 22 + s*(27), B_field_fac*(0.4472135954999579*tot_Bx[1])); 
    gkyl_mat_set(&lhs, 19 + s*(27), 13 + s*(27), -B_field_fac*(0.4472135954999579*tot_Bx[1])); 
 
    gkyl_mat_set(&lhs, 19 + s*(27), 4 + s*(27), B_field_fac*(0.4472135954999579*tot_By[1])); 
    gkyl_mat_set(&lhs, 1 + s*(27), 22 + s*(27), -B_field_fac*(0.4472135954999579*tot_By[1])); 
 
    gkyl_mat_set(&lhs, 1 + s*(27), 13 + s*(27), B_field_fac*(0.4472135954999579*tot_Bz[1])); 
    gkyl_mat_set(&lhs, 10 + s*(27), 4 + s*(27), -B_field_fac*(0.4472135954999579*tot_Bz[1])); 
 
    gkyl_mat_set(&lhs, 1 + num_species*(27), 4 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 10 + num_species*(27), 13 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 19 + num_species*(27), 22 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 1 + s*(27), 5 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 10 + s*(27), 14 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 19 + s*(27), 23 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 1 + s*(27), 5 + num_species*(27), E_field_fac*(0.5000000000000001*rho[s][7])); 
    gkyl_mat_set(&lhs, 10 + s*(27), 14 + num_species*(27), E_field_fac*(0.5000000000000001*rho[s][7])); 
    gkyl_mat_set(&lhs, 19 + s*(27), 23 + num_species*(27), E_field_fac*(0.5000000000000001*rho[s][7])); 
 
    gkyl_mat_set(&lhs, 10 + s*(27), 23 + s*(27), B_field_fac*(0.5000000000000001*tot_Bx[7])); 
    gkyl_mat_set(&lhs, 19 + s*(27), 14 + s*(27), -B_field_fac*(0.5000000000000001*tot_Bx[7])); 
 
    gkyl_mat_set(&lhs, 19 + s*(27), 5 + s*(27), B_field_fac*(0.5000000000000001*tot_By[7])); 
    gkyl_mat_set(&lhs, 1 + s*(27), 23 + s*(27), -B_field_fac*(0.5000000000000001*tot_By[7])); 
 
    gkyl_mat_set(&lhs, 1 + s*(27), 14 + s*(27), B_field_fac*(0.5000000000000001*tot_Bz[7])); 
    gkyl_mat_set(&lhs, 10 + s*(27), 5 + s*(27), -B_field_fac*(0.5000000000000001*tot_Bz[7])); 
 
    gkyl_mat_set(&lhs, 1 + num_species*(27), 5 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 10 + num_species*(27), 14 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 19 + num_species*(27), 23 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 1 + s*(27), 6 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 10 + s*(27), 15 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 19 + s*(27), 24 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 1 + s*(27), 6 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][3])); 
    gkyl_mat_set(&lhs, 10 + s*(27), 15 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][3])); 
    gkyl_mat_set(&lhs, 19 + s*(27), 24 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][3])); 
 
    gkyl_mat_set(&lhs, 10 + s*(27), 24 + s*(27), B_field_fac*(0.447213595499958*tot_Bx[3])); 
    gkyl_mat_set(&lhs, 19 + s*(27), 15 + s*(27), -B_field_fac*(0.447213595499958*tot_Bx[3])); 
 
    gkyl_mat_set(&lhs, 19 + s*(27), 6 + s*(27), B_field_fac*(0.447213595499958*tot_By[3])); 
    gkyl_mat_set(&lhs, 1 + s*(27), 24 + s*(27), -B_field_fac*(0.447213595499958*tot_By[3])); 
 
    gkyl_mat_set(&lhs, 1 + s*(27), 15 + s*(27), B_field_fac*(0.447213595499958*tot_Bz[3])); 
    gkyl_mat_set(&lhs, 10 + s*(27), 6 + s*(27), -B_field_fac*(0.447213595499958*tot_Bz[3])); 
 
    gkyl_mat_set(&lhs, 1 + num_species*(27), 6 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 10 + num_species*(27), 15 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 19 + num_species*(27), 24 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 1 + s*(27), 7 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 10 + s*(27), 16 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 19 + s*(27), 25 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 1 + s*(27), 7 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][8]+0.5000000000000001*rho[s][5])); 
    gkyl_mat_set(&lhs, 10 + s*(27), 16 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][8]+0.5000000000000001*rho[s][5])); 
    gkyl_mat_set(&lhs, 19 + s*(27), 25 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][8]+0.5000000000000001*rho[s][5])); 
 
    gkyl_mat_set(&lhs, 10 + s*(27), 25 + s*(27), B_field_fac*(0.447213595499958*tot_Bx[8]+0.5000000000000001*tot_Bx[5])); 
    gkyl_mat_set(&lhs, 19 + s*(27), 16 + s*(27), -B_field_fac*(0.447213595499958*tot_Bx[8]+0.5000000000000001*tot_Bx[5])); 
 
    gkyl_mat_set(&lhs, 19 + s*(27), 7 + s*(27), B_field_fac*(0.447213595499958*tot_By[8]+0.5000000000000001*tot_By[5])); 
    gkyl_mat_set(&lhs, 1 + s*(27), 25 + s*(27), -B_field_fac*(0.447213595499958*tot_By[8]+0.5000000000000001*tot_By[5])); 
 
    gkyl_mat_set(&lhs, 1 + s*(27), 16 + s*(27), B_field_fac*(0.447213595499958*tot_Bz[8]+0.5000000000000001*tot_Bz[5])); 
    gkyl_mat_set(&lhs, 10 + s*(27), 7 + s*(27), -B_field_fac*(0.447213595499958*tot_Bz[8]+0.5000000000000001*tot_Bz[5])); 
 
    gkyl_mat_set(&lhs, 1 + num_species*(27), 7 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 10 + num_species*(27), 16 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 19 + num_species*(27), 25 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 1 + s*(27), 8 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 10 + s*(27), 17 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 19 + s*(27), 26 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 1 + s*(27), 8 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][7])); 
    gkyl_mat_set(&lhs, 10 + s*(27), 17 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][7])); 
    gkyl_mat_set(&lhs, 19 + s*(27), 26 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][7])); 
 
    gkyl_mat_set(&lhs, 10 + s*(27), 26 + s*(27), B_field_fac*(0.447213595499958*tot_Bx[7])); 
    gkyl_mat_set(&lhs, 19 + s*(27), 17 + s*(27), -B_field_fac*(0.447213595499958*tot_Bx[7])); 
 
    gkyl_mat_set(&lhs, 19 + s*(27), 8 + s*(27), B_field_fac*(0.447213595499958*tot_By[7])); 
    gkyl_mat_set(&lhs, 1 + s*(27), 26 + s*(27), -B_field_fac*(0.447213595499958*tot_By[7])); 
 
    gkyl_mat_set(&lhs, 1 + s*(27), 17 + s*(27), B_field_fac*(0.447213595499958*tot_Bz[7])); 
    gkyl_mat_set(&lhs, 10 + s*(27), 8 + s*(27), -B_field_fac*(0.447213595499958*tot_Bz[7])); 
 
    gkyl_mat_set(&lhs, 1 + num_species*(27), 8 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 10 + num_species*(27), 17 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 19 + num_species*(27), 26 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 2 + s*(27), 0 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 11 + s*(27), 9 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 20 + s*(27), 18 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 2 + s*(27), 0 + num_species*(27), E_field_fac*(0.5*rho[s][2])); 
    gkyl_mat_set(&lhs, 11 + s*(27), 9 + num_species*(27), E_field_fac*(0.5*rho[s][2])); 
    gkyl_mat_set(&lhs, 20 + s*(27), 18 + num_species*(27), E_field_fac*(0.5*rho[s][2])); 
 
    gkyl_mat_set(&lhs, 11 + s*(27), 18 + s*(27), B_field_fac*(0.5*tot_Bx[2])); 
    gkyl_mat_set(&lhs, 20 + s*(27), 9 + s*(27), -B_field_fac*(0.5*tot_Bx[2])); 
 
    gkyl_mat_set(&lhs, 20 + s*(27), 0 + s*(27), B_field_fac*(0.5*tot_By[2])); 
    gkyl_mat_set(&lhs, 2 + s*(27), 18 + s*(27), -B_field_fac*(0.5*tot_By[2])); 
 
    gkyl_mat_set(&lhs, 2 + s*(27), 9 + s*(27), B_field_fac*(0.5*tot_Bz[2])); 
    gkyl_mat_set(&lhs, 11 + s*(27), 0 + s*(27), -B_field_fac*(0.5*tot_Bz[2])); 
 
    gkyl_mat_set(&lhs, 2 + num_species*(27), 0 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 11 + num_species*(27), 9 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 20 + num_species*(27), 18 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 2 + s*(27), 1 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 11 + s*(27), 10 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 20 + s*(27), 19 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 2 + s*(27), 1 + num_species*(27), E_field_fac*(0.5*rho[s][3])); 
    gkyl_mat_set(&lhs, 11 + s*(27), 10 + num_species*(27), E_field_fac*(0.5*rho[s][3])); 
    gkyl_mat_set(&lhs, 20 + s*(27), 19 + num_species*(27), E_field_fac*(0.5*rho[s][3])); 
 
    gkyl_mat_set(&lhs, 11 + s*(27), 19 + s*(27), B_field_fac*(0.5*tot_Bx[3])); 
    gkyl_mat_set(&lhs, 20 + s*(27), 10 + s*(27), -B_field_fac*(0.5*tot_Bx[3])); 
 
    gkyl_mat_set(&lhs, 20 + s*(27), 1 + s*(27), B_field_fac*(0.5*tot_By[3])); 
    gkyl_mat_set(&lhs, 2 + s*(27), 19 + s*(27), -B_field_fac*(0.5*tot_By[3])); 
 
    gkyl_mat_set(&lhs, 2 + s*(27), 10 + s*(27), B_field_fac*(0.5*tot_Bz[3])); 
    gkyl_mat_set(&lhs, 11 + s*(27), 1 + s*(27), -B_field_fac*(0.5*tot_Bz[3])); 
 
    gkyl_mat_set(&lhs, 2 + num_species*(27), 1 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 11 + num_species*(27), 10 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 20 + num_species*(27), 19 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 2 + s*(27), 2 + s*(27), 1.0); 
    gkyl_mat_set(&lhs, 11 + s*(27), 11 + s*(27), 1.0); 
    gkyl_mat_set(&lhs, 20 + s*(27), 20 + s*(27), 1.0); 
 
    gkyl_mat_set(&lhs, 2 + s*(27), 2 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][5]+0.5*rho[s][0])); 
    gkyl_mat_set(&lhs, 11 + s*(27), 11 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][5]+0.5*rho[s][0])); 
    gkyl_mat_set(&lhs, 20 + s*(27), 20 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][5]+0.5*rho[s][0])); 
 
    gkyl_mat_set(&lhs, 11 + s*(27), 20 + s*(27), B_field_fac*(0.4472135954999579*tot_Bx[5]+0.5*tot_Bx[0])); 
    gkyl_mat_set(&lhs, 20 + s*(27), 11 + s*(27), -B_field_fac*(0.4472135954999579*tot_Bx[5]+0.5*tot_Bx[0])); 
 
    gkyl_mat_set(&lhs, 20 + s*(27), 2 + s*(27), B_field_fac*(0.4472135954999579*tot_By[5]+0.5*tot_By[0])); 
    gkyl_mat_set(&lhs, 2 + s*(27), 20 + s*(27), -B_field_fac*(0.4472135954999579*tot_By[5]+0.5*tot_By[0])); 
 
    gkyl_mat_set(&lhs, 2 + s*(27), 11 + s*(27), B_field_fac*(0.4472135954999579*tot_Bz[5]+0.5*tot_Bz[0])); 
    gkyl_mat_set(&lhs, 11 + s*(27), 2 + s*(27), -B_field_fac*(0.4472135954999579*tot_Bz[5]+0.5*tot_Bz[0])); 
 
    gkyl_mat_set(&lhs, 2 + num_species*(27), 2 + s*(27), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 11 + num_species*(27), 11 + s*(27), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 20 + num_species*(27), 20 + s*(27), 0.5*dt*(1.0)); 
 
    gkyl_mat_set(&lhs, 2 + s*(27), 3 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 11 + s*(27), 12 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 20 + s*(27), 21 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 2 + s*(27), 3 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][7]+0.5*rho[s][1])); 
    gkyl_mat_set(&lhs, 11 + s*(27), 12 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][7]+0.5*rho[s][1])); 
    gkyl_mat_set(&lhs, 20 + s*(27), 21 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][7]+0.5*rho[s][1])); 
 
    gkyl_mat_set(&lhs, 11 + s*(27), 21 + s*(27), B_field_fac*(0.447213595499958*tot_Bx[7]+0.5*tot_Bx[1])); 
    gkyl_mat_set(&lhs, 20 + s*(27), 12 + s*(27), -B_field_fac*(0.447213595499958*tot_Bx[7]+0.5*tot_Bx[1])); 
 
    gkyl_mat_set(&lhs, 20 + s*(27), 3 + s*(27), B_field_fac*(0.447213595499958*tot_By[7]+0.5*tot_By[1])); 
    gkyl_mat_set(&lhs, 2 + s*(27), 21 + s*(27), -B_field_fac*(0.447213595499958*tot_By[7]+0.5*tot_By[1])); 
 
    gkyl_mat_set(&lhs, 2 + s*(27), 12 + s*(27), B_field_fac*(0.447213595499958*tot_Bz[7]+0.5*tot_Bz[1])); 
    gkyl_mat_set(&lhs, 11 + s*(27), 3 + s*(27), -B_field_fac*(0.447213595499958*tot_Bz[7]+0.5*tot_Bz[1])); 
 
    gkyl_mat_set(&lhs, 2 + num_species*(27), 3 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 11 + num_species*(27), 12 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 20 + num_species*(27), 21 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 2 + s*(27), 4 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 11 + s*(27), 13 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 20 + s*(27), 22 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 2 + s*(27), 4 + num_species*(27), E_field_fac*(0.5000000000000001*rho[s][6])); 
    gkyl_mat_set(&lhs, 11 + s*(27), 13 + num_species*(27), E_field_fac*(0.5000000000000001*rho[s][6])); 
    gkyl_mat_set(&lhs, 20 + s*(27), 22 + num_species*(27), E_field_fac*(0.5000000000000001*rho[s][6])); 
 
    gkyl_mat_set(&lhs, 11 + s*(27), 22 + s*(27), B_field_fac*(0.5000000000000001*tot_Bx[6])); 
    gkyl_mat_set(&lhs, 20 + s*(27), 13 + s*(27), -B_field_fac*(0.5000000000000001*tot_Bx[6])); 
 
    gkyl_mat_set(&lhs, 20 + s*(27), 4 + s*(27), B_field_fac*(0.5000000000000001*tot_By[6])); 
    gkyl_mat_set(&lhs, 2 + s*(27), 22 + s*(27), -B_field_fac*(0.5000000000000001*tot_By[6])); 
 
    gkyl_mat_set(&lhs, 2 + s*(27), 13 + s*(27), B_field_fac*(0.5000000000000001*tot_Bz[6])); 
    gkyl_mat_set(&lhs, 11 + s*(27), 4 + s*(27), -B_field_fac*(0.5000000000000001*tot_Bz[6])); 
 
    gkyl_mat_set(&lhs, 2 + num_species*(27), 4 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 11 + num_species*(27), 13 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 20 + num_species*(27), 22 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 2 + s*(27), 5 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 11 + s*(27), 14 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 20 + s*(27), 23 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 2 + s*(27), 5 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][2])); 
    gkyl_mat_set(&lhs, 11 + s*(27), 14 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][2])); 
    gkyl_mat_set(&lhs, 20 + s*(27), 23 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][2])); 
 
    gkyl_mat_set(&lhs, 11 + s*(27), 23 + s*(27), B_field_fac*(0.4472135954999579*tot_Bx[2])); 
    gkyl_mat_set(&lhs, 20 + s*(27), 14 + s*(27), -B_field_fac*(0.4472135954999579*tot_Bx[2])); 
 
    gkyl_mat_set(&lhs, 20 + s*(27), 5 + s*(27), B_field_fac*(0.4472135954999579*tot_By[2])); 
    gkyl_mat_set(&lhs, 2 + s*(27), 23 + s*(27), -B_field_fac*(0.4472135954999579*tot_By[2])); 
 
    gkyl_mat_set(&lhs, 2 + s*(27), 14 + s*(27), B_field_fac*(0.4472135954999579*tot_Bz[2])); 
    gkyl_mat_set(&lhs, 11 + s*(27), 5 + s*(27), -B_field_fac*(0.4472135954999579*tot_Bz[2])); 
 
    gkyl_mat_set(&lhs, 2 + num_species*(27), 5 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 11 + num_species*(27), 14 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 20 + num_species*(27), 23 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 2 + s*(27), 6 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 11 + s*(27), 15 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 20 + s*(27), 24 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 2 + s*(27), 6 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][8]+0.5000000000000001*rho[s][4])); 
    gkyl_mat_set(&lhs, 11 + s*(27), 15 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][8]+0.5000000000000001*rho[s][4])); 
    gkyl_mat_set(&lhs, 20 + s*(27), 24 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][8]+0.5000000000000001*rho[s][4])); 
 
    gkyl_mat_set(&lhs, 11 + s*(27), 24 + s*(27), B_field_fac*(0.447213595499958*tot_Bx[8]+0.5000000000000001*tot_Bx[4])); 
    gkyl_mat_set(&lhs, 20 + s*(27), 15 + s*(27), -B_field_fac*(0.447213595499958*tot_Bx[8]+0.5000000000000001*tot_Bx[4])); 
 
    gkyl_mat_set(&lhs, 20 + s*(27), 6 + s*(27), B_field_fac*(0.447213595499958*tot_By[8]+0.5000000000000001*tot_By[4])); 
    gkyl_mat_set(&lhs, 2 + s*(27), 24 + s*(27), -B_field_fac*(0.447213595499958*tot_By[8]+0.5000000000000001*tot_By[4])); 
 
    gkyl_mat_set(&lhs, 2 + s*(27), 15 + s*(27), B_field_fac*(0.447213595499958*tot_Bz[8]+0.5000000000000001*tot_Bz[4])); 
    gkyl_mat_set(&lhs, 11 + s*(27), 6 + s*(27), -B_field_fac*(0.447213595499958*tot_Bz[8]+0.5000000000000001*tot_Bz[4])); 
 
    gkyl_mat_set(&lhs, 2 + num_species*(27), 6 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 11 + num_species*(27), 15 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 20 + num_species*(27), 24 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 2 + s*(27), 7 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 11 + s*(27), 16 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 20 + s*(27), 25 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 2 + s*(27), 7 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][3])); 
    gkyl_mat_set(&lhs, 11 + s*(27), 16 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][3])); 
    gkyl_mat_set(&lhs, 20 + s*(27), 25 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][3])); 
 
    gkyl_mat_set(&lhs, 11 + s*(27), 25 + s*(27), B_field_fac*(0.447213595499958*tot_Bx[3])); 
    gkyl_mat_set(&lhs, 20 + s*(27), 16 + s*(27), -B_field_fac*(0.447213595499958*tot_Bx[3])); 
 
    gkyl_mat_set(&lhs, 20 + s*(27), 7 + s*(27), B_field_fac*(0.447213595499958*tot_By[3])); 
    gkyl_mat_set(&lhs, 2 + s*(27), 25 + s*(27), -B_field_fac*(0.447213595499958*tot_By[3])); 
 
    gkyl_mat_set(&lhs, 2 + s*(27), 16 + s*(27), B_field_fac*(0.447213595499958*tot_Bz[3])); 
    gkyl_mat_set(&lhs, 11 + s*(27), 7 + s*(27), -B_field_fac*(0.447213595499958*tot_Bz[3])); 
 
    gkyl_mat_set(&lhs, 2 + num_species*(27), 7 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 11 + num_species*(27), 16 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 20 + num_species*(27), 25 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 2 + s*(27), 8 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 11 + s*(27), 17 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 20 + s*(27), 26 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 2 + s*(27), 8 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][6])); 
    gkyl_mat_set(&lhs, 11 + s*(27), 17 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][6])); 
    gkyl_mat_set(&lhs, 20 + s*(27), 26 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][6])); 
 
    gkyl_mat_set(&lhs, 11 + s*(27), 26 + s*(27), B_field_fac*(0.447213595499958*tot_Bx[6])); 
    gkyl_mat_set(&lhs, 20 + s*(27), 17 + s*(27), -B_field_fac*(0.447213595499958*tot_Bx[6])); 
 
    gkyl_mat_set(&lhs, 20 + s*(27), 8 + s*(27), B_field_fac*(0.447213595499958*tot_By[6])); 
    gkyl_mat_set(&lhs, 2 + s*(27), 26 + s*(27), -B_field_fac*(0.447213595499958*tot_By[6])); 
 
    gkyl_mat_set(&lhs, 2 + s*(27), 17 + s*(27), B_field_fac*(0.447213595499958*tot_Bz[6])); 
    gkyl_mat_set(&lhs, 11 + s*(27), 8 + s*(27), -B_field_fac*(0.447213595499958*tot_Bz[6])); 
 
    gkyl_mat_set(&lhs, 2 + num_species*(27), 8 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 11 + num_species*(27), 17 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 20 + num_species*(27), 26 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 3 + s*(27), 0 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 12 + s*(27), 9 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 21 + s*(27), 18 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 3 + s*(27), 0 + num_species*(27), E_field_fac*(0.5*rho[s][3])); 
    gkyl_mat_set(&lhs, 12 + s*(27), 9 + num_species*(27), E_field_fac*(0.5*rho[s][3])); 
    gkyl_mat_set(&lhs, 21 + s*(27), 18 + num_species*(27), E_field_fac*(0.5*rho[s][3])); 
 
    gkyl_mat_set(&lhs, 12 + s*(27), 18 + s*(27), B_field_fac*(0.5*tot_Bx[3])); 
    gkyl_mat_set(&lhs, 21 + s*(27), 9 + s*(27), -B_field_fac*(0.5*tot_Bx[3])); 
 
    gkyl_mat_set(&lhs, 21 + s*(27), 0 + s*(27), B_field_fac*(0.5*tot_By[3])); 
    gkyl_mat_set(&lhs, 3 + s*(27), 18 + s*(27), -B_field_fac*(0.5*tot_By[3])); 
 
    gkyl_mat_set(&lhs, 3 + s*(27), 9 + s*(27), B_field_fac*(0.5*tot_Bz[3])); 
    gkyl_mat_set(&lhs, 12 + s*(27), 0 + s*(27), -B_field_fac*(0.5*tot_Bz[3])); 
 
    gkyl_mat_set(&lhs, 3 + num_species*(27), 0 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 12 + num_species*(27), 9 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 21 + num_species*(27), 18 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 3 + s*(27), 1 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 12 + s*(27), 10 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 21 + s*(27), 19 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 3 + s*(27), 1 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][6]+0.5*rho[s][2])); 
    gkyl_mat_set(&lhs, 12 + s*(27), 10 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][6]+0.5*rho[s][2])); 
    gkyl_mat_set(&lhs, 21 + s*(27), 19 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][6]+0.5*rho[s][2])); 
 
    gkyl_mat_set(&lhs, 12 + s*(27), 19 + s*(27), B_field_fac*(0.447213595499958*tot_Bx[6]+0.5*tot_Bx[2])); 
    gkyl_mat_set(&lhs, 21 + s*(27), 10 + s*(27), -B_field_fac*(0.447213595499958*tot_Bx[6]+0.5*tot_Bx[2])); 
 
    gkyl_mat_set(&lhs, 21 + s*(27), 1 + s*(27), B_field_fac*(0.447213595499958*tot_By[6]+0.5*tot_By[2])); 
    gkyl_mat_set(&lhs, 3 + s*(27), 19 + s*(27), -B_field_fac*(0.447213595499958*tot_By[6]+0.5*tot_By[2])); 
 
    gkyl_mat_set(&lhs, 3 + s*(27), 10 + s*(27), B_field_fac*(0.447213595499958*tot_Bz[6]+0.5*tot_Bz[2])); 
    gkyl_mat_set(&lhs, 12 + s*(27), 1 + s*(27), -B_field_fac*(0.447213595499958*tot_Bz[6]+0.5*tot_Bz[2])); 
 
    gkyl_mat_set(&lhs, 3 + num_species*(27), 1 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 12 + num_species*(27), 10 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 21 + num_species*(27), 19 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 3 + s*(27), 2 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 12 + s*(27), 11 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 21 + s*(27), 20 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 3 + s*(27), 2 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][7]+0.5*rho[s][1])); 
    gkyl_mat_set(&lhs, 12 + s*(27), 11 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][7]+0.5*rho[s][1])); 
    gkyl_mat_set(&lhs, 21 + s*(27), 20 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][7]+0.5*rho[s][1])); 
 
    gkyl_mat_set(&lhs, 12 + s*(27), 20 + s*(27), B_field_fac*(0.447213595499958*tot_Bx[7]+0.5*tot_Bx[1])); 
    gkyl_mat_set(&lhs, 21 + s*(27), 11 + s*(27), -B_field_fac*(0.447213595499958*tot_Bx[7]+0.5*tot_Bx[1])); 
 
    gkyl_mat_set(&lhs, 21 + s*(27), 2 + s*(27), B_field_fac*(0.447213595499958*tot_By[7]+0.5*tot_By[1])); 
    gkyl_mat_set(&lhs, 3 + s*(27), 20 + s*(27), -B_field_fac*(0.447213595499958*tot_By[7]+0.5*tot_By[1])); 
 
    gkyl_mat_set(&lhs, 3 + s*(27), 11 + s*(27), B_field_fac*(0.447213595499958*tot_Bz[7]+0.5*tot_Bz[1])); 
    gkyl_mat_set(&lhs, 12 + s*(27), 2 + s*(27), -B_field_fac*(0.447213595499958*tot_Bz[7]+0.5*tot_Bz[1])); 
 
    gkyl_mat_set(&lhs, 3 + num_species*(27), 2 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 12 + num_species*(27), 11 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 21 + num_species*(27), 20 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 3 + s*(27), 3 + s*(27), 1.0); 
    gkyl_mat_set(&lhs, 12 + s*(27), 12 + s*(27), 1.0); 
    gkyl_mat_set(&lhs, 21 + s*(27), 21 + s*(27), 1.0); 
 
    gkyl_mat_set(&lhs, 3 + s*(27), 3 + num_species*(27), E_field_fac*(0.4*rho[s][8]+0.4472135954999579*rho[s][5]+0.4472135954999579*rho[s][4]+0.5*rho[s][0])); 
    gkyl_mat_set(&lhs, 12 + s*(27), 12 + num_species*(27), E_field_fac*(0.4*rho[s][8]+0.4472135954999579*rho[s][5]+0.4472135954999579*rho[s][4]+0.5*rho[s][0])); 
    gkyl_mat_set(&lhs, 21 + s*(27), 21 + num_species*(27), E_field_fac*(0.4*rho[s][8]+0.4472135954999579*rho[s][5]+0.4472135954999579*rho[s][4]+0.5*rho[s][0])); 
 
    gkyl_mat_set(&lhs, 12 + s*(27), 21 + s*(27), B_field_fac*(0.4*tot_Bx[8]+0.4472135954999579*tot_Bx[5]+0.4472135954999579*tot_Bx[4]+0.5*tot_Bx[0])); 
    gkyl_mat_set(&lhs, 21 + s*(27), 12 + s*(27), -B_field_fac*(0.4*tot_Bx[8]+0.4472135954999579*tot_Bx[5]+0.4472135954999579*tot_Bx[4]+0.5*tot_Bx[0])); 
 
    gkyl_mat_set(&lhs, 21 + s*(27), 3 + s*(27), B_field_fac*(0.4*tot_By[8]+0.4472135954999579*tot_By[5]+0.4472135954999579*tot_By[4]+0.5*tot_By[0])); 
    gkyl_mat_set(&lhs, 3 + s*(27), 21 + s*(27), -B_field_fac*(0.4*tot_By[8]+0.4472135954999579*tot_By[5]+0.4472135954999579*tot_By[4]+0.5*tot_By[0])); 
 
    gkyl_mat_set(&lhs, 3 + s*(27), 12 + s*(27), B_field_fac*(0.4*tot_Bz[8]+0.4472135954999579*tot_Bz[5]+0.4472135954999579*tot_Bz[4]+0.5*tot_Bz[0])); 
    gkyl_mat_set(&lhs, 12 + s*(27), 3 + s*(27), -B_field_fac*(0.4*tot_Bz[8]+0.4472135954999579*tot_Bz[5]+0.4472135954999579*tot_Bz[4]+0.5*tot_Bz[0])); 
 
    gkyl_mat_set(&lhs, 3 + num_species*(27), 3 + s*(27), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 12 + num_species*(27), 12 + s*(27), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 21 + num_species*(27), 21 + s*(27), 0.5*dt*(1.0)); 
 
    gkyl_mat_set(&lhs, 3 + s*(27), 4 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 12 + s*(27), 13 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 21 + s*(27), 22 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 3 + s*(27), 4 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][3])); 
    gkyl_mat_set(&lhs, 12 + s*(27), 13 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][3])); 
    gkyl_mat_set(&lhs, 21 + s*(27), 22 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][3])); 
 
    gkyl_mat_set(&lhs, 12 + s*(27), 22 + s*(27), B_field_fac*(0.4472135954999579*tot_Bx[3])); 
    gkyl_mat_set(&lhs, 21 + s*(27), 13 + s*(27), -B_field_fac*(0.4472135954999579*tot_Bx[3])); 
 
    gkyl_mat_set(&lhs, 21 + s*(27), 4 + s*(27), B_field_fac*(0.4472135954999579*tot_By[3])); 
    gkyl_mat_set(&lhs, 3 + s*(27), 22 + s*(27), -B_field_fac*(0.4472135954999579*tot_By[3])); 
 
    gkyl_mat_set(&lhs, 3 + s*(27), 13 + s*(27), B_field_fac*(0.4472135954999579*tot_Bz[3])); 
    gkyl_mat_set(&lhs, 12 + s*(27), 4 + s*(27), -B_field_fac*(0.4472135954999579*tot_Bz[3])); 
 
    gkyl_mat_set(&lhs, 3 + num_species*(27), 4 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 12 + num_species*(27), 13 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 21 + num_species*(27), 22 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 3 + s*(27), 5 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 12 + s*(27), 14 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 21 + s*(27), 23 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 3 + s*(27), 5 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][3])); 
    gkyl_mat_set(&lhs, 12 + s*(27), 14 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][3])); 
    gkyl_mat_set(&lhs, 21 + s*(27), 23 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][3])); 
 
    gkyl_mat_set(&lhs, 12 + s*(27), 23 + s*(27), B_field_fac*(0.4472135954999579*tot_Bx[3])); 
    gkyl_mat_set(&lhs, 21 + s*(27), 14 + s*(27), -B_field_fac*(0.4472135954999579*tot_Bx[3])); 
 
    gkyl_mat_set(&lhs, 21 + s*(27), 5 + s*(27), B_field_fac*(0.4472135954999579*tot_By[3])); 
    gkyl_mat_set(&lhs, 3 + s*(27), 23 + s*(27), -B_field_fac*(0.4472135954999579*tot_By[3])); 
 
    gkyl_mat_set(&lhs, 3 + s*(27), 14 + s*(27), B_field_fac*(0.4472135954999579*tot_Bz[3])); 
    gkyl_mat_set(&lhs, 12 + s*(27), 5 + s*(27), -B_field_fac*(0.4472135954999579*tot_Bz[3])); 
 
    gkyl_mat_set(&lhs, 3 + num_species*(27), 5 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 12 + num_species*(27), 14 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 21 + num_species*(27), 23 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 3 + s*(27), 6 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 12 + s*(27), 15 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 21 + s*(27), 24 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 3 + s*(27), 6 + num_species*(27), E_field_fac*(0.4*rho[s][7]+0.447213595499958*rho[s][1])); 
    gkyl_mat_set(&lhs, 12 + s*(27), 15 + num_species*(27), E_field_fac*(0.4*rho[s][7]+0.447213595499958*rho[s][1])); 
    gkyl_mat_set(&lhs, 21 + s*(27), 24 + num_species*(27), E_field_fac*(0.4*rho[s][7]+0.447213595499958*rho[s][1])); 
 
    gkyl_mat_set(&lhs, 12 + s*(27), 24 + s*(27), B_field_fac*(0.4*tot_Bx[7]+0.447213595499958*tot_Bx[1])); 
    gkyl_mat_set(&lhs, 21 + s*(27), 15 + s*(27), -B_field_fac*(0.4*tot_Bx[7]+0.447213595499958*tot_Bx[1])); 
 
    gkyl_mat_set(&lhs, 21 + s*(27), 6 + s*(27), B_field_fac*(0.4*tot_By[7]+0.447213595499958*tot_By[1])); 
    gkyl_mat_set(&lhs, 3 + s*(27), 24 + s*(27), -B_field_fac*(0.4*tot_By[7]+0.447213595499958*tot_By[1])); 
 
    gkyl_mat_set(&lhs, 3 + s*(27), 15 + s*(27), B_field_fac*(0.4*tot_Bz[7]+0.447213595499958*tot_Bz[1])); 
    gkyl_mat_set(&lhs, 12 + s*(27), 6 + s*(27), -B_field_fac*(0.4*tot_Bz[7]+0.447213595499958*tot_Bz[1])); 
 
    gkyl_mat_set(&lhs, 3 + num_species*(27), 6 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 12 + num_species*(27), 15 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 21 + num_species*(27), 24 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 3 + s*(27), 7 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 12 + s*(27), 16 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 21 + s*(27), 25 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 3 + s*(27), 7 + num_species*(27), E_field_fac*(0.4*rho[s][6]+0.447213595499958*rho[s][2])); 
    gkyl_mat_set(&lhs, 12 + s*(27), 16 + num_species*(27), E_field_fac*(0.4*rho[s][6]+0.447213595499958*rho[s][2])); 
    gkyl_mat_set(&lhs, 21 + s*(27), 25 + num_species*(27), E_field_fac*(0.4*rho[s][6]+0.447213595499958*rho[s][2])); 
 
    gkyl_mat_set(&lhs, 12 + s*(27), 25 + s*(27), B_field_fac*(0.4*tot_Bx[6]+0.447213595499958*tot_Bx[2])); 
    gkyl_mat_set(&lhs, 21 + s*(27), 16 + s*(27), -B_field_fac*(0.4*tot_Bx[6]+0.447213595499958*tot_Bx[2])); 
 
    gkyl_mat_set(&lhs, 21 + s*(27), 7 + s*(27), B_field_fac*(0.4*tot_By[6]+0.447213595499958*tot_By[2])); 
    gkyl_mat_set(&lhs, 3 + s*(27), 25 + s*(27), -B_field_fac*(0.4*tot_By[6]+0.447213595499958*tot_By[2])); 
 
    gkyl_mat_set(&lhs, 3 + s*(27), 16 + s*(27), B_field_fac*(0.4*tot_Bz[6]+0.447213595499958*tot_Bz[2])); 
    gkyl_mat_set(&lhs, 12 + s*(27), 7 + s*(27), -B_field_fac*(0.4*tot_Bz[6]+0.447213595499958*tot_Bz[2])); 
 
    gkyl_mat_set(&lhs, 3 + num_species*(27), 7 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 12 + num_species*(27), 16 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 21 + num_species*(27), 25 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 3 + s*(27), 8 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 12 + s*(27), 17 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 21 + s*(27), 26 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 3 + s*(27), 8 + num_species*(27), E_field_fac*(0.4*rho[s][3])); 
    gkyl_mat_set(&lhs, 12 + s*(27), 17 + num_species*(27), E_field_fac*(0.4*rho[s][3])); 
    gkyl_mat_set(&lhs, 21 + s*(27), 26 + num_species*(27), E_field_fac*(0.4*rho[s][3])); 
 
    gkyl_mat_set(&lhs, 12 + s*(27), 26 + s*(27), B_field_fac*(0.4*tot_Bx[3])); 
    gkyl_mat_set(&lhs, 21 + s*(27), 17 + s*(27), -B_field_fac*(0.4*tot_Bx[3])); 
 
    gkyl_mat_set(&lhs, 21 + s*(27), 8 + s*(27), B_field_fac*(0.4*tot_By[3])); 
    gkyl_mat_set(&lhs, 3 + s*(27), 26 + s*(27), -B_field_fac*(0.4*tot_By[3])); 
 
    gkyl_mat_set(&lhs, 3 + s*(27), 17 + s*(27), B_field_fac*(0.4*tot_Bz[3])); 
    gkyl_mat_set(&lhs, 12 + s*(27), 8 + s*(27), -B_field_fac*(0.4*tot_Bz[3])); 
 
    gkyl_mat_set(&lhs, 3 + num_species*(27), 8 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 12 + num_species*(27), 17 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 21 + num_species*(27), 26 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 4 + s*(27), 0 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 13 + s*(27), 9 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 22 + s*(27), 18 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 4 + s*(27), 0 + num_species*(27), E_field_fac*(0.5*rho[s][4])); 
    gkyl_mat_set(&lhs, 13 + s*(27), 9 + num_species*(27), E_field_fac*(0.5*rho[s][4])); 
    gkyl_mat_set(&lhs, 22 + s*(27), 18 + num_species*(27), E_field_fac*(0.5*rho[s][4])); 
 
    gkyl_mat_set(&lhs, 13 + s*(27), 18 + s*(27), B_field_fac*(0.5*tot_Bx[4])); 
    gkyl_mat_set(&lhs, 22 + s*(27), 9 + s*(27), -B_field_fac*(0.5*tot_Bx[4])); 
 
    gkyl_mat_set(&lhs, 22 + s*(27), 0 + s*(27), B_field_fac*(0.5*tot_By[4])); 
    gkyl_mat_set(&lhs, 4 + s*(27), 18 + s*(27), -B_field_fac*(0.5*tot_By[4])); 
 
    gkyl_mat_set(&lhs, 4 + s*(27), 9 + s*(27), B_field_fac*(0.5*tot_Bz[4])); 
    gkyl_mat_set(&lhs, 13 + s*(27), 0 + s*(27), -B_field_fac*(0.5*tot_Bz[4])); 
 
    gkyl_mat_set(&lhs, 4 + num_species*(27), 0 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 13 + num_species*(27), 9 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 22 + num_species*(27), 18 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 4 + s*(27), 1 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 13 + s*(27), 10 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 22 + s*(27), 19 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 4 + s*(27), 1 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][1])); 
    gkyl_mat_set(&lhs, 13 + s*(27), 10 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][1])); 
    gkyl_mat_set(&lhs, 22 + s*(27), 19 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][1])); 
 
    gkyl_mat_set(&lhs, 13 + s*(27), 19 + s*(27), B_field_fac*(0.4472135954999579*tot_Bx[1])); 
    gkyl_mat_set(&lhs, 22 + s*(27), 10 + s*(27), -B_field_fac*(0.4472135954999579*tot_Bx[1])); 
 
    gkyl_mat_set(&lhs, 22 + s*(27), 1 + s*(27), B_field_fac*(0.4472135954999579*tot_By[1])); 
    gkyl_mat_set(&lhs, 4 + s*(27), 19 + s*(27), -B_field_fac*(0.4472135954999579*tot_By[1])); 
 
    gkyl_mat_set(&lhs, 4 + s*(27), 10 + s*(27), B_field_fac*(0.4472135954999579*tot_Bz[1])); 
    gkyl_mat_set(&lhs, 13 + s*(27), 1 + s*(27), -B_field_fac*(0.4472135954999579*tot_Bz[1])); 
 
    gkyl_mat_set(&lhs, 4 + num_species*(27), 1 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 13 + num_species*(27), 10 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 22 + num_species*(27), 19 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 4 + s*(27), 2 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 13 + s*(27), 11 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 22 + s*(27), 20 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 4 + s*(27), 2 + num_species*(27), E_field_fac*(0.5000000000000001*rho[s][6])); 
    gkyl_mat_set(&lhs, 13 + s*(27), 11 + num_species*(27), E_field_fac*(0.5000000000000001*rho[s][6])); 
    gkyl_mat_set(&lhs, 22 + s*(27), 20 + num_species*(27), E_field_fac*(0.5000000000000001*rho[s][6])); 
 
    gkyl_mat_set(&lhs, 13 + s*(27), 20 + s*(27), B_field_fac*(0.5000000000000001*tot_Bx[6])); 
    gkyl_mat_set(&lhs, 22 + s*(27), 11 + s*(27), -B_field_fac*(0.5000000000000001*tot_Bx[6])); 
 
    gkyl_mat_set(&lhs, 22 + s*(27), 2 + s*(27), B_field_fac*(0.5000000000000001*tot_By[6])); 
    gkyl_mat_set(&lhs, 4 + s*(27), 20 + s*(27), -B_field_fac*(0.5000000000000001*tot_By[6])); 
 
    gkyl_mat_set(&lhs, 4 + s*(27), 11 + s*(27), B_field_fac*(0.5000000000000001*tot_Bz[6])); 
    gkyl_mat_set(&lhs, 13 + s*(27), 2 + s*(27), -B_field_fac*(0.5000000000000001*tot_Bz[6])); 
 
    gkyl_mat_set(&lhs, 4 + num_species*(27), 2 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 13 + num_species*(27), 11 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 22 + num_species*(27), 20 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 4 + s*(27), 3 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 13 + s*(27), 12 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 22 + s*(27), 21 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 4 + s*(27), 3 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][3])); 
    gkyl_mat_set(&lhs, 13 + s*(27), 12 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][3])); 
    gkyl_mat_set(&lhs, 22 + s*(27), 21 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][3])); 
 
    gkyl_mat_set(&lhs, 13 + s*(27), 21 + s*(27), B_field_fac*(0.4472135954999579*tot_Bx[3])); 
    gkyl_mat_set(&lhs, 22 + s*(27), 12 + s*(27), -B_field_fac*(0.4472135954999579*tot_Bx[3])); 
 
    gkyl_mat_set(&lhs, 22 + s*(27), 3 + s*(27), B_field_fac*(0.4472135954999579*tot_By[3])); 
    gkyl_mat_set(&lhs, 4 + s*(27), 21 + s*(27), -B_field_fac*(0.4472135954999579*tot_By[3])); 
 
    gkyl_mat_set(&lhs, 4 + s*(27), 12 + s*(27), B_field_fac*(0.4472135954999579*tot_Bz[3])); 
    gkyl_mat_set(&lhs, 13 + s*(27), 3 + s*(27), -B_field_fac*(0.4472135954999579*tot_Bz[3])); 
 
    gkyl_mat_set(&lhs, 4 + num_species*(27), 3 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 13 + num_species*(27), 12 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 22 + num_species*(27), 21 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 4 + s*(27), 4 + s*(27), 1.0); 
    gkyl_mat_set(&lhs, 13 + s*(27), 13 + s*(27), 1.0); 
    gkyl_mat_set(&lhs, 22 + s*(27), 22 + s*(27), 1.0); 
 
    gkyl_mat_set(&lhs, 4 + s*(27), 4 + num_species*(27), E_field_fac*(0.31943828249997*rho[s][4]+0.5*rho[s][0])); 
    gkyl_mat_set(&lhs, 13 + s*(27), 13 + num_species*(27), E_field_fac*(0.31943828249997*rho[s][4]+0.5*rho[s][0])); 
    gkyl_mat_set(&lhs, 22 + s*(27), 22 + num_species*(27), E_field_fac*(0.31943828249997*rho[s][4]+0.5*rho[s][0])); 
 
    gkyl_mat_set(&lhs, 13 + s*(27), 22 + s*(27), B_field_fac*(0.31943828249997*tot_Bx[4]+0.5*tot_Bx[0])); 
    gkyl_mat_set(&lhs, 22 + s*(27), 13 + s*(27), -B_field_fac*(0.31943828249997*tot_Bx[4]+0.5*tot_Bx[0])); 
 
    gkyl_mat_set(&lhs, 22 + s*(27), 4 + s*(27), B_field_fac*(0.31943828249997*tot_By[4]+0.5*tot_By[0])); 
    gkyl_mat_set(&lhs, 4 + s*(27), 22 + s*(27), -B_field_fac*(0.31943828249997*tot_By[4]+0.5*tot_By[0])); 
 
    gkyl_mat_set(&lhs, 4 + s*(27), 13 + s*(27), B_field_fac*(0.31943828249997*tot_Bz[4]+0.5*tot_Bz[0])); 
    gkyl_mat_set(&lhs, 13 + s*(27), 4 + s*(27), -B_field_fac*(0.31943828249997*tot_Bz[4]+0.5*tot_Bz[0])); 
 
    gkyl_mat_set(&lhs, 4 + num_species*(27), 4 + s*(27), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 13 + num_species*(27), 13 + s*(27), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 22 + num_species*(27), 22 + s*(27), 0.5*dt*(1.0)); 
 
    gkyl_mat_set(&lhs, 4 + s*(27), 5 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 13 + s*(27), 14 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 22 + s*(27), 23 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 4 + s*(27), 5 + num_species*(27), E_field_fac*(0.5*rho[s][8])); 
    gkyl_mat_set(&lhs, 13 + s*(27), 14 + num_species*(27), E_field_fac*(0.5*rho[s][8])); 
    gkyl_mat_set(&lhs, 22 + s*(27), 23 + num_species*(27), E_field_fac*(0.5*rho[s][8])); 
 
    gkyl_mat_set(&lhs, 13 + s*(27), 23 + s*(27), B_field_fac*(0.5*tot_Bx[8])); 
    gkyl_mat_set(&lhs, 22 + s*(27), 14 + s*(27), -B_field_fac*(0.5*tot_Bx[8])); 
 
    gkyl_mat_set(&lhs, 22 + s*(27), 5 + s*(27), B_field_fac*(0.5*tot_By[8])); 
    gkyl_mat_set(&lhs, 4 + s*(27), 23 + s*(27), -B_field_fac*(0.5*tot_By[8])); 
 
    gkyl_mat_set(&lhs, 4 + s*(27), 14 + s*(27), B_field_fac*(0.5*tot_Bz[8])); 
    gkyl_mat_set(&lhs, 13 + s*(27), 5 + s*(27), -B_field_fac*(0.5*tot_Bz[8])); 
 
    gkyl_mat_set(&lhs, 4 + num_species*(27), 5 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 13 + num_species*(27), 14 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 22 + num_species*(27), 23 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 4 + s*(27), 6 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 13 + s*(27), 15 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 22 + s*(27), 24 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 4 + s*(27), 6 + num_species*(27), E_field_fac*(0.31943828249997*rho[s][6]+0.5000000000000001*rho[s][2])); 
    gkyl_mat_set(&lhs, 13 + s*(27), 15 + num_species*(27), E_field_fac*(0.31943828249997*rho[s][6]+0.5000000000000001*rho[s][2])); 
    gkyl_mat_set(&lhs, 22 + s*(27), 24 + num_species*(27), E_field_fac*(0.31943828249997*rho[s][6]+0.5000000000000001*rho[s][2])); 
 
    gkyl_mat_set(&lhs, 13 + s*(27), 24 + s*(27), B_field_fac*(0.31943828249997*tot_Bx[6]+0.5000000000000001*tot_Bx[2])); 
    gkyl_mat_set(&lhs, 22 + s*(27), 15 + s*(27), -B_field_fac*(0.31943828249997*tot_Bx[6]+0.5000000000000001*tot_Bx[2])); 
 
    gkyl_mat_set(&lhs, 22 + s*(27), 6 + s*(27), B_field_fac*(0.31943828249997*tot_By[6]+0.5000000000000001*tot_By[2])); 
    gkyl_mat_set(&lhs, 4 + s*(27), 24 + s*(27), -B_field_fac*(0.31943828249997*tot_By[6]+0.5000000000000001*tot_By[2])); 
 
    gkyl_mat_set(&lhs, 4 + s*(27), 15 + s*(27), B_field_fac*(0.31943828249997*tot_Bz[6]+0.5000000000000001*tot_Bz[2])); 
    gkyl_mat_set(&lhs, 13 + s*(27), 6 + s*(27), -B_field_fac*(0.31943828249997*tot_Bz[6]+0.5000000000000001*tot_Bz[2])); 
 
    gkyl_mat_set(&lhs, 4 + num_species*(27), 6 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 13 + num_species*(27), 15 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 22 + num_species*(27), 24 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 4 + s*(27), 7 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 13 + s*(27), 16 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 22 + s*(27), 25 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 4 + s*(27), 7 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][7])); 
    gkyl_mat_set(&lhs, 13 + s*(27), 16 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][7])); 
    gkyl_mat_set(&lhs, 22 + s*(27), 25 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][7])); 
 
    gkyl_mat_set(&lhs, 13 + s*(27), 25 + s*(27), B_field_fac*(0.4472135954999579*tot_Bx[7])); 
    gkyl_mat_set(&lhs, 22 + s*(27), 16 + s*(27), -B_field_fac*(0.4472135954999579*tot_Bx[7])); 
 
    gkyl_mat_set(&lhs, 22 + s*(27), 7 + s*(27), B_field_fac*(0.4472135954999579*tot_By[7])); 
    gkyl_mat_set(&lhs, 4 + s*(27), 25 + s*(27), -B_field_fac*(0.4472135954999579*tot_By[7])); 
 
    gkyl_mat_set(&lhs, 4 + s*(27), 16 + s*(27), B_field_fac*(0.4472135954999579*tot_Bz[7])); 
    gkyl_mat_set(&lhs, 13 + s*(27), 7 + s*(27), -B_field_fac*(0.4472135954999579*tot_Bz[7])); 
 
    gkyl_mat_set(&lhs, 4 + num_species*(27), 7 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 13 + num_species*(27), 16 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 22 + num_species*(27), 25 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 4 + s*(27), 8 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 13 + s*(27), 17 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 22 + s*(27), 26 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 4 + s*(27), 8 + num_species*(27), E_field_fac*(0.31943828249997*rho[s][8]+0.5*rho[s][5])); 
    gkyl_mat_set(&lhs, 13 + s*(27), 17 + num_species*(27), E_field_fac*(0.31943828249997*rho[s][8]+0.5*rho[s][5])); 
    gkyl_mat_set(&lhs, 22 + s*(27), 26 + num_species*(27), E_field_fac*(0.31943828249997*rho[s][8]+0.5*rho[s][5])); 
 
    gkyl_mat_set(&lhs, 13 + s*(27), 26 + s*(27), B_field_fac*(0.31943828249997*tot_Bx[8]+0.5*tot_Bx[5])); 
    gkyl_mat_set(&lhs, 22 + s*(27), 17 + s*(27), -B_field_fac*(0.31943828249997*tot_Bx[8]+0.5*tot_Bx[5])); 
 
    gkyl_mat_set(&lhs, 22 + s*(27), 8 + s*(27), B_field_fac*(0.31943828249997*tot_By[8]+0.5*tot_By[5])); 
    gkyl_mat_set(&lhs, 4 + s*(27), 26 + s*(27), -B_field_fac*(0.31943828249997*tot_By[8]+0.5*tot_By[5])); 
 
    gkyl_mat_set(&lhs, 4 + s*(27), 17 + s*(27), B_field_fac*(0.31943828249997*tot_Bz[8]+0.5*tot_Bz[5])); 
    gkyl_mat_set(&lhs, 13 + s*(27), 8 + s*(27), -B_field_fac*(0.31943828249997*tot_Bz[8]+0.5*tot_Bz[5])); 
 
    gkyl_mat_set(&lhs, 4 + num_species*(27), 8 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 13 + num_species*(27), 17 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 22 + num_species*(27), 26 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 5 + s*(27), 0 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 14 + s*(27), 9 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 23 + s*(27), 18 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 5 + s*(27), 0 + num_species*(27), E_field_fac*(0.5*rho[s][5])); 
    gkyl_mat_set(&lhs, 14 + s*(27), 9 + num_species*(27), E_field_fac*(0.5*rho[s][5])); 
    gkyl_mat_set(&lhs, 23 + s*(27), 18 + num_species*(27), E_field_fac*(0.5*rho[s][5])); 
 
    gkyl_mat_set(&lhs, 14 + s*(27), 18 + s*(27), B_field_fac*(0.5*tot_Bx[5])); 
    gkyl_mat_set(&lhs, 23 + s*(27), 9 + s*(27), -B_field_fac*(0.5*tot_Bx[5])); 
 
    gkyl_mat_set(&lhs, 23 + s*(27), 0 + s*(27), B_field_fac*(0.5*tot_By[5])); 
    gkyl_mat_set(&lhs, 5 + s*(27), 18 + s*(27), -B_field_fac*(0.5*tot_By[5])); 
 
    gkyl_mat_set(&lhs, 5 + s*(27), 9 + s*(27), B_field_fac*(0.5*tot_Bz[5])); 
    gkyl_mat_set(&lhs, 14 + s*(27), 0 + s*(27), -B_field_fac*(0.5*tot_Bz[5])); 
 
    gkyl_mat_set(&lhs, 5 + num_species*(27), 0 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 14 + num_species*(27), 9 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 23 + num_species*(27), 18 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 5 + s*(27), 1 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 14 + s*(27), 10 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 23 + s*(27), 19 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 5 + s*(27), 1 + num_species*(27), E_field_fac*(0.5000000000000001*rho[s][7])); 
    gkyl_mat_set(&lhs, 14 + s*(27), 10 + num_species*(27), E_field_fac*(0.5000000000000001*rho[s][7])); 
    gkyl_mat_set(&lhs, 23 + s*(27), 19 + num_species*(27), E_field_fac*(0.5000000000000001*rho[s][7])); 
 
    gkyl_mat_set(&lhs, 14 + s*(27), 19 + s*(27), B_field_fac*(0.5000000000000001*tot_Bx[7])); 
    gkyl_mat_set(&lhs, 23 + s*(27), 10 + s*(27), -B_field_fac*(0.5000000000000001*tot_Bx[7])); 
 
    gkyl_mat_set(&lhs, 23 + s*(27), 1 + s*(27), B_field_fac*(0.5000000000000001*tot_By[7])); 
    gkyl_mat_set(&lhs, 5 + s*(27), 19 + s*(27), -B_field_fac*(0.5000000000000001*tot_By[7])); 
 
    gkyl_mat_set(&lhs, 5 + s*(27), 10 + s*(27), B_field_fac*(0.5000000000000001*tot_Bz[7])); 
    gkyl_mat_set(&lhs, 14 + s*(27), 1 + s*(27), -B_field_fac*(0.5000000000000001*tot_Bz[7])); 
 
    gkyl_mat_set(&lhs, 5 + num_species*(27), 1 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 14 + num_species*(27), 10 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 23 + num_species*(27), 19 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 5 + s*(27), 2 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 14 + s*(27), 11 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 23 + s*(27), 20 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 5 + s*(27), 2 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][2])); 
    gkyl_mat_set(&lhs, 14 + s*(27), 11 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][2])); 
    gkyl_mat_set(&lhs, 23 + s*(27), 20 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][2])); 
 
    gkyl_mat_set(&lhs, 14 + s*(27), 20 + s*(27), B_field_fac*(0.4472135954999579*tot_Bx[2])); 
    gkyl_mat_set(&lhs, 23 + s*(27), 11 + s*(27), -B_field_fac*(0.4472135954999579*tot_Bx[2])); 
 
    gkyl_mat_set(&lhs, 23 + s*(27), 2 + s*(27), B_field_fac*(0.4472135954999579*tot_By[2])); 
    gkyl_mat_set(&lhs, 5 + s*(27), 20 + s*(27), -B_field_fac*(0.4472135954999579*tot_By[2])); 
 
    gkyl_mat_set(&lhs, 5 + s*(27), 11 + s*(27), B_field_fac*(0.4472135954999579*tot_Bz[2])); 
    gkyl_mat_set(&lhs, 14 + s*(27), 2 + s*(27), -B_field_fac*(0.4472135954999579*tot_Bz[2])); 
 
    gkyl_mat_set(&lhs, 5 + num_species*(27), 2 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 14 + num_species*(27), 11 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 23 + num_species*(27), 20 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 5 + s*(27), 3 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 14 + s*(27), 12 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 23 + s*(27), 21 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 5 + s*(27), 3 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][3])); 
    gkyl_mat_set(&lhs, 14 + s*(27), 12 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][3])); 
    gkyl_mat_set(&lhs, 23 + s*(27), 21 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][3])); 
 
    gkyl_mat_set(&lhs, 14 + s*(27), 21 + s*(27), B_field_fac*(0.4472135954999579*tot_Bx[3])); 
    gkyl_mat_set(&lhs, 23 + s*(27), 12 + s*(27), -B_field_fac*(0.4472135954999579*tot_Bx[3])); 
 
    gkyl_mat_set(&lhs, 23 + s*(27), 3 + s*(27), B_field_fac*(0.4472135954999579*tot_By[3])); 
    gkyl_mat_set(&lhs, 5 + s*(27), 21 + s*(27), -B_field_fac*(0.4472135954999579*tot_By[3])); 
 
    gkyl_mat_set(&lhs, 5 + s*(27), 12 + s*(27), B_field_fac*(0.4472135954999579*tot_Bz[3])); 
    gkyl_mat_set(&lhs, 14 + s*(27), 3 + s*(27), -B_field_fac*(0.4472135954999579*tot_Bz[3])); 
 
    gkyl_mat_set(&lhs, 5 + num_species*(27), 3 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 14 + num_species*(27), 12 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 23 + num_species*(27), 21 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 5 + s*(27), 4 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 14 + s*(27), 13 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 23 + s*(27), 22 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 5 + s*(27), 4 + num_species*(27), E_field_fac*(0.5*rho[s][8])); 
    gkyl_mat_set(&lhs, 14 + s*(27), 13 + num_species*(27), E_field_fac*(0.5*rho[s][8])); 
    gkyl_mat_set(&lhs, 23 + s*(27), 22 + num_species*(27), E_field_fac*(0.5*rho[s][8])); 
 
    gkyl_mat_set(&lhs, 14 + s*(27), 22 + s*(27), B_field_fac*(0.5*tot_Bx[8])); 
    gkyl_mat_set(&lhs, 23 + s*(27), 13 + s*(27), -B_field_fac*(0.5*tot_Bx[8])); 
 
    gkyl_mat_set(&lhs, 23 + s*(27), 4 + s*(27), B_field_fac*(0.5*tot_By[8])); 
    gkyl_mat_set(&lhs, 5 + s*(27), 22 + s*(27), -B_field_fac*(0.5*tot_By[8])); 
 
    gkyl_mat_set(&lhs, 5 + s*(27), 13 + s*(27), B_field_fac*(0.5*tot_Bz[8])); 
    gkyl_mat_set(&lhs, 14 + s*(27), 4 + s*(27), -B_field_fac*(0.5*tot_Bz[8])); 
 
    gkyl_mat_set(&lhs, 5 + num_species*(27), 4 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 14 + num_species*(27), 13 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 23 + num_species*(27), 22 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 5 + s*(27), 5 + s*(27), 1.0); 
    gkyl_mat_set(&lhs, 14 + s*(27), 14 + s*(27), 1.0); 
    gkyl_mat_set(&lhs, 23 + s*(27), 23 + s*(27), 1.0); 
 
    gkyl_mat_set(&lhs, 5 + s*(27), 5 + num_species*(27), E_field_fac*(0.31943828249997*rho[s][5]+0.5*rho[s][0])); 
    gkyl_mat_set(&lhs, 14 + s*(27), 14 + num_species*(27), E_field_fac*(0.31943828249997*rho[s][5]+0.5*rho[s][0])); 
    gkyl_mat_set(&lhs, 23 + s*(27), 23 + num_species*(27), E_field_fac*(0.31943828249997*rho[s][5]+0.5*rho[s][0])); 
 
    gkyl_mat_set(&lhs, 14 + s*(27), 23 + s*(27), B_field_fac*(0.31943828249997*tot_Bx[5]+0.5*tot_Bx[0])); 
    gkyl_mat_set(&lhs, 23 + s*(27), 14 + s*(27), -B_field_fac*(0.31943828249997*tot_Bx[5]+0.5*tot_Bx[0])); 
 
    gkyl_mat_set(&lhs, 23 + s*(27), 5 + s*(27), B_field_fac*(0.31943828249997*tot_By[5]+0.5*tot_By[0])); 
    gkyl_mat_set(&lhs, 5 + s*(27), 23 + s*(27), -B_field_fac*(0.31943828249997*tot_By[5]+0.5*tot_By[0])); 
 
    gkyl_mat_set(&lhs, 5 + s*(27), 14 + s*(27), B_field_fac*(0.31943828249997*tot_Bz[5]+0.5*tot_Bz[0])); 
    gkyl_mat_set(&lhs, 14 + s*(27), 5 + s*(27), -B_field_fac*(0.31943828249997*tot_Bz[5]+0.5*tot_Bz[0])); 
 
    gkyl_mat_set(&lhs, 5 + num_species*(27), 5 + s*(27), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 14 + num_species*(27), 14 + s*(27), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 23 + num_species*(27), 23 + s*(27), 0.5*dt*(1.0)); 
 
    gkyl_mat_set(&lhs, 5 + s*(27), 6 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 14 + s*(27), 15 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 23 + s*(27), 24 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 5 + s*(27), 6 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][6])); 
    gkyl_mat_set(&lhs, 14 + s*(27), 15 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][6])); 
    gkyl_mat_set(&lhs, 23 + s*(27), 24 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][6])); 
 
    gkyl_mat_set(&lhs, 14 + s*(27), 24 + s*(27), B_field_fac*(0.4472135954999579*tot_Bx[6])); 
    gkyl_mat_set(&lhs, 23 + s*(27), 15 + s*(27), -B_field_fac*(0.4472135954999579*tot_Bx[6])); 
 
    gkyl_mat_set(&lhs, 23 + s*(27), 6 + s*(27), B_field_fac*(0.4472135954999579*tot_By[6])); 
    gkyl_mat_set(&lhs, 5 + s*(27), 24 + s*(27), -B_field_fac*(0.4472135954999579*tot_By[6])); 
 
    gkyl_mat_set(&lhs, 5 + s*(27), 15 + s*(27), B_field_fac*(0.4472135954999579*tot_Bz[6])); 
    gkyl_mat_set(&lhs, 14 + s*(27), 6 + s*(27), -B_field_fac*(0.4472135954999579*tot_Bz[6])); 
 
    gkyl_mat_set(&lhs, 5 + num_species*(27), 6 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 14 + num_species*(27), 15 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 23 + num_species*(27), 24 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 5 + s*(27), 7 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 14 + s*(27), 16 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 23 + s*(27), 25 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 5 + s*(27), 7 + num_species*(27), E_field_fac*(0.31943828249997*rho[s][7]+0.5000000000000001*rho[s][1])); 
    gkyl_mat_set(&lhs, 14 + s*(27), 16 + num_species*(27), E_field_fac*(0.31943828249997*rho[s][7]+0.5000000000000001*rho[s][1])); 
    gkyl_mat_set(&lhs, 23 + s*(27), 25 + num_species*(27), E_field_fac*(0.31943828249997*rho[s][7]+0.5000000000000001*rho[s][1])); 
 
    gkyl_mat_set(&lhs, 14 + s*(27), 25 + s*(27), B_field_fac*(0.31943828249997*tot_Bx[7]+0.5000000000000001*tot_Bx[1])); 
    gkyl_mat_set(&lhs, 23 + s*(27), 16 + s*(27), -B_field_fac*(0.31943828249997*tot_Bx[7]+0.5000000000000001*tot_Bx[1])); 
 
    gkyl_mat_set(&lhs, 23 + s*(27), 7 + s*(27), B_field_fac*(0.31943828249997*tot_By[7]+0.5000000000000001*tot_By[1])); 
    gkyl_mat_set(&lhs, 5 + s*(27), 25 + s*(27), -B_field_fac*(0.31943828249997*tot_By[7]+0.5000000000000001*tot_By[1])); 
 
    gkyl_mat_set(&lhs, 5 + s*(27), 16 + s*(27), B_field_fac*(0.31943828249997*tot_Bz[7]+0.5000000000000001*tot_Bz[1])); 
    gkyl_mat_set(&lhs, 14 + s*(27), 7 + s*(27), -B_field_fac*(0.31943828249997*tot_Bz[7]+0.5000000000000001*tot_Bz[1])); 
 
    gkyl_mat_set(&lhs, 5 + num_species*(27), 7 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 14 + num_species*(27), 16 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 23 + num_species*(27), 25 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 5 + s*(27), 8 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 14 + s*(27), 17 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 23 + s*(27), 26 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 5 + s*(27), 8 + num_species*(27), E_field_fac*(0.31943828249997*rho[s][8]+0.5*rho[s][4])); 
    gkyl_mat_set(&lhs, 14 + s*(27), 17 + num_species*(27), E_field_fac*(0.31943828249997*rho[s][8]+0.5*rho[s][4])); 
    gkyl_mat_set(&lhs, 23 + s*(27), 26 + num_species*(27), E_field_fac*(0.31943828249997*rho[s][8]+0.5*rho[s][4])); 
 
    gkyl_mat_set(&lhs, 14 + s*(27), 26 + s*(27), B_field_fac*(0.31943828249997*tot_Bx[8]+0.5*tot_Bx[4])); 
    gkyl_mat_set(&lhs, 23 + s*(27), 17 + s*(27), -B_field_fac*(0.31943828249997*tot_Bx[8]+0.5*tot_Bx[4])); 
 
    gkyl_mat_set(&lhs, 23 + s*(27), 8 + s*(27), B_field_fac*(0.31943828249997*tot_By[8]+0.5*tot_By[4])); 
    gkyl_mat_set(&lhs, 5 + s*(27), 26 + s*(27), -B_field_fac*(0.31943828249997*tot_By[8]+0.5*tot_By[4])); 
 
    gkyl_mat_set(&lhs, 5 + s*(27), 17 + s*(27), B_field_fac*(0.31943828249997*tot_Bz[8]+0.5*tot_Bz[4])); 
    gkyl_mat_set(&lhs, 14 + s*(27), 8 + s*(27), -B_field_fac*(0.31943828249997*tot_Bz[8]+0.5*tot_Bz[4])); 
 
    gkyl_mat_set(&lhs, 5 + num_species*(27), 8 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 14 + num_species*(27), 17 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 23 + num_species*(27), 26 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 6 + s*(27), 0 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 15 + s*(27), 9 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 24 + s*(27), 18 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 6 + s*(27), 0 + num_species*(27), E_field_fac*(0.5*rho[s][6])); 
    gkyl_mat_set(&lhs, 15 + s*(27), 9 + num_species*(27), E_field_fac*(0.5*rho[s][6])); 
    gkyl_mat_set(&lhs, 24 + s*(27), 18 + num_species*(27), E_field_fac*(0.5*rho[s][6])); 
 
    gkyl_mat_set(&lhs, 15 + s*(27), 18 + s*(27), B_field_fac*(0.5*tot_Bx[6])); 
    gkyl_mat_set(&lhs, 24 + s*(27), 9 + s*(27), -B_field_fac*(0.5*tot_Bx[6])); 
 
    gkyl_mat_set(&lhs, 24 + s*(27), 0 + s*(27), B_field_fac*(0.5*tot_By[6])); 
    gkyl_mat_set(&lhs, 6 + s*(27), 18 + s*(27), -B_field_fac*(0.5*tot_By[6])); 
 
    gkyl_mat_set(&lhs, 6 + s*(27), 9 + s*(27), B_field_fac*(0.5*tot_Bz[6])); 
    gkyl_mat_set(&lhs, 15 + s*(27), 0 + s*(27), -B_field_fac*(0.5*tot_Bz[6])); 
 
    gkyl_mat_set(&lhs, 6 + num_species*(27), 0 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 15 + num_species*(27), 9 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 24 + num_species*(27), 18 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 6 + s*(27), 1 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 15 + s*(27), 10 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 24 + s*(27), 19 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 6 + s*(27), 1 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][3])); 
    gkyl_mat_set(&lhs, 15 + s*(27), 10 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][3])); 
    gkyl_mat_set(&lhs, 24 + s*(27), 19 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][3])); 
 
    gkyl_mat_set(&lhs, 15 + s*(27), 19 + s*(27), B_field_fac*(0.447213595499958*tot_Bx[3])); 
    gkyl_mat_set(&lhs, 24 + s*(27), 10 + s*(27), -B_field_fac*(0.447213595499958*tot_Bx[3])); 
 
    gkyl_mat_set(&lhs, 24 + s*(27), 1 + s*(27), B_field_fac*(0.447213595499958*tot_By[3])); 
    gkyl_mat_set(&lhs, 6 + s*(27), 19 + s*(27), -B_field_fac*(0.447213595499958*tot_By[3])); 
 
    gkyl_mat_set(&lhs, 6 + s*(27), 10 + s*(27), B_field_fac*(0.447213595499958*tot_Bz[3])); 
    gkyl_mat_set(&lhs, 15 + s*(27), 1 + s*(27), -B_field_fac*(0.447213595499958*tot_Bz[3])); 
 
    gkyl_mat_set(&lhs, 6 + num_species*(27), 1 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 15 + num_species*(27), 10 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 24 + num_species*(27), 19 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 6 + s*(27), 2 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 15 + s*(27), 11 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 24 + s*(27), 20 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 6 + s*(27), 2 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][8]+0.5000000000000001*rho[s][4])); 
    gkyl_mat_set(&lhs, 15 + s*(27), 11 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][8]+0.5000000000000001*rho[s][4])); 
    gkyl_mat_set(&lhs, 24 + s*(27), 20 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][8]+0.5000000000000001*rho[s][4])); 
 
    gkyl_mat_set(&lhs, 15 + s*(27), 20 + s*(27), B_field_fac*(0.447213595499958*tot_Bx[8]+0.5000000000000001*tot_Bx[4])); 
    gkyl_mat_set(&lhs, 24 + s*(27), 11 + s*(27), -B_field_fac*(0.447213595499958*tot_Bx[8]+0.5000000000000001*tot_Bx[4])); 
 
    gkyl_mat_set(&lhs, 24 + s*(27), 2 + s*(27), B_field_fac*(0.447213595499958*tot_By[8]+0.5000000000000001*tot_By[4])); 
    gkyl_mat_set(&lhs, 6 + s*(27), 20 + s*(27), -B_field_fac*(0.447213595499958*tot_By[8]+0.5000000000000001*tot_By[4])); 
 
    gkyl_mat_set(&lhs, 6 + s*(27), 11 + s*(27), B_field_fac*(0.447213595499958*tot_Bz[8]+0.5000000000000001*tot_Bz[4])); 
    gkyl_mat_set(&lhs, 15 + s*(27), 2 + s*(27), -B_field_fac*(0.447213595499958*tot_Bz[8]+0.5000000000000001*tot_Bz[4])); 
 
    gkyl_mat_set(&lhs, 6 + num_species*(27), 2 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 15 + num_species*(27), 11 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 24 + num_species*(27), 20 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 6 + s*(27), 3 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 15 + s*(27), 12 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 24 + s*(27), 21 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 6 + s*(27), 3 + num_species*(27), E_field_fac*(0.4*rho[s][7]+0.447213595499958*rho[s][1])); 
    gkyl_mat_set(&lhs, 15 + s*(27), 12 + num_species*(27), E_field_fac*(0.4*rho[s][7]+0.447213595499958*rho[s][1])); 
    gkyl_mat_set(&lhs, 24 + s*(27), 21 + num_species*(27), E_field_fac*(0.4*rho[s][7]+0.447213595499958*rho[s][1])); 
 
    gkyl_mat_set(&lhs, 15 + s*(27), 21 + s*(27), B_field_fac*(0.4*tot_Bx[7]+0.447213595499958*tot_Bx[1])); 
    gkyl_mat_set(&lhs, 24 + s*(27), 12 + s*(27), -B_field_fac*(0.4*tot_Bx[7]+0.447213595499958*tot_Bx[1])); 
 
    gkyl_mat_set(&lhs, 24 + s*(27), 3 + s*(27), B_field_fac*(0.4*tot_By[7]+0.447213595499958*tot_By[1])); 
    gkyl_mat_set(&lhs, 6 + s*(27), 21 + s*(27), -B_field_fac*(0.4*tot_By[7]+0.447213595499958*tot_By[1])); 
 
    gkyl_mat_set(&lhs, 6 + s*(27), 12 + s*(27), B_field_fac*(0.4*tot_Bz[7]+0.447213595499958*tot_Bz[1])); 
    gkyl_mat_set(&lhs, 15 + s*(27), 3 + s*(27), -B_field_fac*(0.4*tot_Bz[7]+0.447213595499958*tot_Bz[1])); 
 
    gkyl_mat_set(&lhs, 6 + num_species*(27), 3 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 15 + num_species*(27), 12 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 24 + num_species*(27), 21 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 6 + s*(27), 4 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 15 + s*(27), 13 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 24 + s*(27), 22 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 6 + s*(27), 4 + num_species*(27), E_field_fac*(0.31943828249997*rho[s][6]+0.5000000000000001*rho[s][2])); 
    gkyl_mat_set(&lhs, 15 + s*(27), 13 + num_species*(27), E_field_fac*(0.31943828249997*rho[s][6]+0.5000000000000001*rho[s][2])); 
    gkyl_mat_set(&lhs, 24 + s*(27), 22 + num_species*(27), E_field_fac*(0.31943828249997*rho[s][6]+0.5000000000000001*rho[s][2])); 
 
    gkyl_mat_set(&lhs, 15 + s*(27), 22 + s*(27), B_field_fac*(0.31943828249997*tot_Bx[6]+0.5000000000000001*tot_Bx[2])); 
    gkyl_mat_set(&lhs, 24 + s*(27), 13 + s*(27), -B_field_fac*(0.31943828249997*tot_Bx[6]+0.5000000000000001*tot_Bx[2])); 
 
    gkyl_mat_set(&lhs, 24 + s*(27), 4 + s*(27), B_field_fac*(0.31943828249997*tot_By[6]+0.5000000000000001*tot_By[2])); 
    gkyl_mat_set(&lhs, 6 + s*(27), 22 + s*(27), -B_field_fac*(0.31943828249997*tot_By[6]+0.5000000000000001*tot_By[2])); 
 
    gkyl_mat_set(&lhs, 6 + s*(27), 13 + s*(27), B_field_fac*(0.31943828249997*tot_Bz[6]+0.5000000000000001*tot_Bz[2])); 
    gkyl_mat_set(&lhs, 15 + s*(27), 4 + s*(27), -B_field_fac*(0.31943828249997*tot_Bz[6]+0.5000000000000001*tot_Bz[2])); 
 
    gkyl_mat_set(&lhs, 6 + num_species*(27), 4 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 15 + num_species*(27), 13 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 24 + num_species*(27), 22 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 6 + s*(27), 5 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 15 + s*(27), 14 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 24 + s*(27), 23 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 6 + s*(27), 5 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][6])); 
    gkyl_mat_set(&lhs, 15 + s*(27), 14 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][6])); 
    gkyl_mat_set(&lhs, 24 + s*(27), 23 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][6])); 
 
    gkyl_mat_set(&lhs, 15 + s*(27), 23 + s*(27), B_field_fac*(0.4472135954999579*tot_Bx[6])); 
    gkyl_mat_set(&lhs, 24 + s*(27), 14 + s*(27), -B_field_fac*(0.4472135954999579*tot_Bx[6])); 
 
    gkyl_mat_set(&lhs, 24 + s*(27), 5 + s*(27), B_field_fac*(0.4472135954999579*tot_By[6])); 
    gkyl_mat_set(&lhs, 6 + s*(27), 23 + s*(27), -B_field_fac*(0.4472135954999579*tot_By[6])); 
 
    gkyl_mat_set(&lhs, 6 + s*(27), 14 + s*(27), B_field_fac*(0.4472135954999579*tot_Bz[6])); 
    gkyl_mat_set(&lhs, 15 + s*(27), 5 + s*(27), -B_field_fac*(0.4472135954999579*tot_Bz[6])); 
 
    gkyl_mat_set(&lhs, 6 + num_species*(27), 5 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 15 + num_species*(27), 14 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 24 + num_species*(27), 23 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 6 + s*(27), 6 + s*(27), 1.0); 
    gkyl_mat_set(&lhs, 15 + s*(27), 15 + s*(27), 1.0); 
    gkyl_mat_set(&lhs, 24 + s*(27), 24 + s*(27), 1.0); 
 
    gkyl_mat_set(&lhs, 6 + s*(27), 6 + num_species*(27), E_field_fac*(0.2857142857142857*rho[s][8]+0.4472135954999579*rho[s][5]+0.31943828249997*rho[s][4]+0.5*rho[s][0])); 
    gkyl_mat_set(&lhs, 15 + s*(27), 15 + num_species*(27), E_field_fac*(0.2857142857142857*rho[s][8]+0.4472135954999579*rho[s][5]+0.31943828249997*rho[s][4]+0.5*rho[s][0])); 
    gkyl_mat_set(&lhs, 24 + s*(27), 24 + num_species*(27), E_field_fac*(0.2857142857142857*rho[s][8]+0.4472135954999579*rho[s][5]+0.31943828249997*rho[s][4]+0.5*rho[s][0])); 
 
    gkyl_mat_set(&lhs, 15 + s*(27), 24 + s*(27), B_field_fac*(0.2857142857142857*tot_Bx[8]+0.4472135954999579*tot_Bx[5]+0.31943828249997*tot_Bx[4]+0.5*tot_Bx[0])); 
    gkyl_mat_set(&lhs, 24 + s*(27), 15 + s*(27), -B_field_fac*(0.2857142857142857*tot_Bx[8]+0.4472135954999579*tot_Bx[5]+0.31943828249997*tot_Bx[4]+0.5*tot_Bx[0])); 
 
    gkyl_mat_set(&lhs, 24 + s*(27), 6 + s*(27), B_field_fac*(0.2857142857142857*tot_By[8]+0.4472135954999579*tot_By[5]+0.31943828249997*tot_By[4]+0.5*tot_By[0])); 
    gkyl_mat_set(&lhs, 6 + s*(27), 24 + s*(27), -B_field_fac*(0.2857142857142857*tot_By[8]+0.4472135954999579*tot_By[5]+0.31943828249997*tot_By[4]+0.5*tot_By[0])); 
 
    gkyl_mat_set(&lhs, 6 + s*(27), 15 + s*(27), B_field_fac*(0.2857142857142857*tot_Bz[8]+0.4472135954999579*tot_Bz[5]+0.31943828249997*tot_Bz[4]+0.5*tot_Bz[0])); 
    gkyl_mat_set(&lhs, 15 + s*(27), 6 + s*(27), -B_field_fac*(0.2857142857142857*tot_Bz[8]+0.4472135954999579*tot_Bz[5]+0.31943828249997*tot_Bz[4]+0.5*tot_Bz[0])); 
 
    gkyl_mat_set(&lhs, 6 + num_species*(27), 6 + s*(27), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 15 + num_species*(27), 15 + s*(27), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 24 + num_species*(27), 24 + s*(27), 0.5*dt*(1.0)); 
 
    gkyl_mat_set(&lhs, 6 + s*(27), 7 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 15 + s*(27), 16 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 24 + s*(27), 25 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 6 + s*(27), 7 + num_species*(27), E_field_fac*(0.4*rho[s][3])); 
    gkyl_mat_set(&lhs, 15 + s*(27), 16 + num_species*(27), E_field_fac*(0.4*rho[s][3])); 
    gkyl_mat_set(&lhs, 24 + s*(27), 25 + num_species*(27), E_field_fac*(0.4*rho[s][3])); 
 
    gkyl_mat_set(&lhs, 15 + s*(27), 25 + s*(27), B_field_fac*(0.4*tot_Bx[3])); 
    gkyl_mat_set(&lhs, 24 + s*(27), 16 + s*(27), -B_field_fac*(0.4*tot_Bx[3])); 
 
    gkyl_mat_set(&lhs, 24 + s*(27), 7 + s*(27), B_field_fac*(0.4*tot_By[3])); 
    gkyl_mat_set(&lhs, 6 + s*(27), 25 + s*(27), -B_field_fac*(0.4*tot_By[3])); 
 
    gkyl_mat_set(&lhs, 6 + s*(27), 16 + s*(27), B_field_fac*(0.4*tot_Bz[3])); 
    gkyl_mat_set(&lhs, 15 + s*(27), 7 + s*(27), -B_field_fac*(0.4*tot_Bz[3])); 
 
    gkyl_mat_set(&lhs, 6 + num_species*(27), 7 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 15 + num_species*(27), 16 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 24 + num_species*(27), 25 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 6 + s*(27), 8 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 15 + s*(27), 17 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 24 + s*(27), 26 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 6 + s*(27), 8 + num_species*(27), E_field_fac*(0.2857142857142857*rho[s][6]+0.447213595499958*rho[s][2])); 
    gkyl_mat_set(&lhs, 15 + s*(27), 17 + num_species*(27), E_field_fac*(0.2857142857142857*rho[s][6]+0.447213595499958*rho[s][2])); 
    gkyl_mat_set(&lhs, 24 + s*(27), 26 + num_species*(27), E_field_fac*(0.2857142857142857*rho[s][6]+0.447213595499958*rho[s][2])); 
 
    gkyl_mat_set(&lhs, 15 + s*(27), 26 + s*(27), B_field_fac*(0.2857142857142857*tot_Bx[6]+0.447213595499958*tot_Bx[2])); 
    gkyl_mat_set(&lhs, 24 + s*(27), 17 + s*(27), -B_field_fac*(0.2857142857142857*tot_Bx[6]+0.447213595499958*tot_Bx[2])); 
 
    gkyl_mat_set(&lhs, 24 + s*(27), 8 + s*(27), B_field_fac*(0.2857142857142857*tot_By[6]+0.447213595499958*tot_By[2])); 
    gkyl_mat_set(&lhs, 6 + s*(27), 26 + s*(27), -B_field_fac*(0.2857142857142857*tot_By[6]+0.447213595499958*tot_By[2])); 
 
    gkyl_mat_set(&lhs, 6 + s*(27), 17 + s*(27), B_field_fac*(0.2857142857142857*tot_Bz[6]+0.447213595499958*tot_Bz[2])); 
    gkyl_mat_set(&lhs, 15 + s*(27), 8 + s*(27), -B_field_fac*(0.2857142857142857*tot_Bz[6]+0.447213595499958*tot_Bz[2])); 
 
    gkyl_mat_set(&lhs, 6 + num_species*(27), 8 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 15 + num_species*(27), 17 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 24 + num_species*(27), 26 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 7 + s*(27), 0 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 16 + s*(27), 9 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 25 + s*(27), 18 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 7 + s*(27), 0 + num_species*(27), E_field_fac*(0.5*rho[s][7])); 
    gkyl_mat_set(&lhs, 16 + s*(27), 9 + num_species*(27), E_field_fac*(0.5*rho[s][7])); 
    gkyl_mat_set(&lhs, 25 + s*(27), 18 + num_species*(27), E_field_fac*(0.5*rho[s][7])); 
 
    gkyl_mat_set(&lhs, 16 + s*(27), 18 + s*(27), B_field_fac*(0.5*tot_Bx[7])); 
    gkyl_mat_set(&lhs, 25 + s*(27), 9 + s*(27), -B_field_fac*(0.5*tot_Bx[7])); 
 
    gkyl_mat_set(&lhs, 25 + s*(27), 0 + s*(27), B_field_fac*(0.5*tot_By[7])); 
    gkyl_mat_set(&lhs, 7 + s*(27), 18 + s*(27), -B_field_fac*(0.5*tot_By[7])); 
 
    gkyl_mat_set(&lhs, 7 + s*(27), 9 + s*(27), B_field_fac*(0.5*tot_Bz[7])); 
    gkyl_mat_set(&lhs, 16 + s*(27), 0 + s*(27), -B_field_fac*(0.5*tot_Bz[7])); 
 
    gkyl_mat_set(&lhs, 7 + num_species*(27), 0 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 16 + num_species*(27), 9 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 25 + num_species*(27), 18 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 7 + s*(27), 1 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 16 + s*(27), 10 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 25 + s*(27), 19 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 7 + s*(27), 1 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][8]+0.5000000000000001*rho[s][5])); 
    gkyl_mat_set(&lhs, 16 + s*(27), 10 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][8]+0.5000000000000001*rho[s][5])); 
    gkyl_mat_set(&lhs, 25 + s*(27), 19 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][8]+0.5000000000000001*rho[s][5])); 
 
    gkyl_mat_set(&lhs, 16 + s*(27), 19 + s*(27), B_field_fac*(0.447213595499958*tot_Bx[8]+0.5000000000000001*tot_Bx[5])); 
    gkyl_mat_set(&lhs, 25 + s*(27), 10 + s*(27), -B_field_fac*(0.447213595499958*tot_Bx[8]+0.5000000000000001*tot_Bx[5])); 
 
    gkyl_mat_set(&lhs, 25 + s*(27), 1 + s*(27), B_field_fac*(0.447213595499958*tot_By[8]+0.5000000000000001*tot_By[5])); 
    gkyl_mat_set(&lhs, 7 + s*(27), 19 + s*(27), -B_field_fac*(0.447213595499958*tot_By[8]+0.5000000000000001*tot_By[5])); 
 
    gkyl_mat_set(&lhs, 7 + s*(27), 10 + s*(27), B_field_fac*(0.447213595499958*tot_Bz[8]+0.5000000000000001*tot_Bz[5])); 
    gkyl_mat_set(&lhs, 16 + s*(27), 1 + s*(27), -B_field_fac*(0.447213595499958*tot_Bz[8]+0.5000000000000001*tot_Bz[5])); 
 
    gkyl_mat_set(&lhs, 7 + num_species*(27), 1 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 16 + num_species*(27), 10 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 25 + num_species*(27), 19 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 7 + s*(27), 2 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 16 + s*(27), 11 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 25 + s*(27), 20 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 7 + s*(27), 2 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][3])); 
    gkyl_mat_set(&lhs, 16 + s*(27), 11 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][3])); 
    gkyl_mat_set(&lhs, 25 + s*(27), 20 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][3])); 
 
    gkyl_mat_set(&lhs, 16 + s*(27), 20 + s*(27), B_field_fac*(0.447213595499958*tot_Bx[3])); 
    gkyl_mat_set(&lhs, 25 + s*(27), 11 + s*(27), -B_field_fac*(0.447213595499958*tot_Bx[3])); 
 
    gkyl_mat_set(&lhs, 25 + s*(27), 2 + s*(27), B_field_fac*(0.447213595499958*tot_By[3])); 
    gkyl_mat_set(&lhs, 7 + s*(27), 20 + s*(27), -B_field_fac*(0.447213595499958*tot_By[3])); 
 
    gkyl_mat_set(&lhs, 7 + s*(27), 11 + s*(27), B_field_fac*(0.447213595499958*tot_Bz[3])); 
    gkyl_mat_set(&lhs, 16 + s*(27), 2 + s*(27), -B_field_fac*(0.447213595499958*tot_Bz[3])); 
 
    gkyl_mat_set(&lhs, 7 + num_species*(27), 2 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 16 + num_species*(27), 11 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 25 + num_species*(27), 20 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 7 + s*(27), 3 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 16 + s*(27), 12 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 25 + s*(27), 21 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 7 + s*(27), 3 + num_species*(27), E_field_fac*(0.4*rho[s][6]+0.447213595499958*rho[s][2])); 
    gkyl_mat_set(&lhs, 16 + s*(27), 12 + num_species*(27), E_field_fac*(0.4*rho[s][6]+0.447213595499958*rho[s][2])); 
    gkyl_mat_set(&lhs, 25 + s*(27), 21 + num_species*(27), E_field_fac*(0.4*rho[s][6]+0.447213595499958*rho[s][2])); 
 
    gkyl_mat_set(&lhs, 16 + s*(27), 21 + s*(27), B_field_fac*(0.4*tot_Bx[6]+0.447213595499958*tot_Bx[2])); 
    gkyl_mat_set(&lhs, 25 + s*(27), 12 + s*(27), -B_field_fac*(0.4*tot_Bx[6]+0.447213595499958*tot_Bx[2])); 
 
    gkyl_mat_set(&lhs, 25 + s*(27), 3 + s*(27), B_field_fac*(0.4*tot_By[6]+0.447213595499958*tot_By[2])); 
    gkyl_mat_set(&lhs, 7 + s*(27), 21 + s*(27), -B_field_fac*(0.4*tot_By[6]+0.447213595499958*tot_By[2])); 
 
    gkyl_mat_set(&lhs, 7 + s*(27), 12 + s*(27), B_field_fac*(0.4*tot_Bz[6]+0.447213595499958*tot_Bz[2])); 
    gkyl_mat_set(&lhs, 16 + s*(27), 3 + s*(27), -B_field_fac*(0.4*tot_Bz[6]+0.447213595499958*tot_Bz[2])); 
 
    gkyl_mat_set(&lhs, 7 + num_species*(27), 3 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 16 + num_species*(27), 12 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 25 + num_species*(27), 21 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 7 + s*(27), 4 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 16 + s*(27), 13 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 25 + s*(27), 22 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 7 + s*(27), 4 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][7])); 
    gkyl_mat_set(&lhs, 16 + s*(27), 13 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][7])); 
    gkyl_mat_set(&lhs, 25 + s*(27), 22 + num_species*(27), E_field_fac*(0.4472135954999579*rho[s][7])); 
 
    gkyl_mat_set(&lhs, 16 + s*(27), 22 + s*(27), B_field_fac*(0.4472135954999579*tot_Bx[7])); 
    gkyl_mat_set(&lhs, 25 + s*(27), 13 + s*(27), -B_field_fac*(0.4472135954999579*tot_Bx[7])); 
 
    gkyl_mat_set(&lhs, 25 + s*(27), 4 + s*(27), B_field_fac*(0.4472135954999579*tot_By[7])); 
    gkyl_mat_set(&lhs, 7 + s*(27), 22 + s*(27), -B_field_fac*(0.4472135954999579*tot_By[7])); 
 
    gkyl_mat_set(&lhs, 7 + s*(27), 13 + s*(27), B_field_fac*(0.4472135954999579*tot_Bz[7])); 
    gkyl_mat_set(&lhs, 16 + s*(27), 4 + s*(27), -B_field_fac*(0.4472135954999579*tot_Bz[7])); 
 
    gkyl_mat_set(&lhs, 7 + num_species*(27), 4 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 16 + num_species*(27), 13 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 25 + num_species*(27), 22 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 7 + s*(27), 5 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 16 + s*(27), 14 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 25 + s*(27), 23 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 7 + s*(27), 5 + num_species*(27), E_field_fac*(0.31943828249997*rho[s][7]+0.5000000000000001*rho[s][1])); 
    gkyl_mat_set(&lhs, 16 + s*(27), 14 + num_species*(27), E_field_fac*(0.31943828249997*rho[s][7]+0.5000000000000001*rho[s][1])); 
    gkyl_mat_set(&lhs, 25 + s*(27), 23 + num_species*(27), E_field_fac*(0.31943828249997*rho[s][7]+0.5000000000000001*rho[s][1])); 
 
    gkyl_mat_set(&lhs, 16 + s*(27), 23 + s*(27), B_field_fac*(0.31943828249997*tot_Bx[7]+0.5000000000000001*tot_Bx[1])); 
    gkyl_mat_set(&lhs, 25 + s*(27), 14 + s*(27), -B_field_fac*(0.31943828249997*tot_Bx[7]+0.5000000000000001*tot_Bx[1])); 
 
    gkyl_mat_set(&lhs, 25 + s*(27), 5 + s*(27), B_field_fac*(0.31943828249997*tot_By[7]+0.5000000000000001*tot_By[1])); 
    gkyl_mat_set(&lhs, 7 + s*(27), 23 + s*(27), -B_field_fac*(0.31943828249997*tot_By[7]+0.5000000000000001*tot_By[1])); 
 
    gkyl_mat_set(&lhs, 7 + s*(27), 14 + s*(27), B_field_fac*(0.31943828249997*tot_Bz[7]+0.5000000000000001*tot_Bz[1])); 
    gkyl_mat_set(&lhs, 16 + s*(27), 5 + s*(27), -B_field_fac*(0.31943828249997*tot_Bz[7]+0.5000000000000001*tot_Bz[1])); 
 
    gkyl_mat_set(&lhs, 7 + num_species*(27), 5 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 16 + num_species*(27), 14 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 25 + num_species*(27), 23 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 7 + s*(27), 6 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 16 + s*(27), 15 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 25 + s*(27), 24 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 7 + s*(27), 6 + num_species*(27), E_field_fac*(0.4*rho[s][3])); 
    gkyl_mat_set(&lhs, 16 + s*(27), 15 + num_species*(27), E_field_fac*(0.4*rho[s][3])); 
    gkyl_mat_set(&lhs, 25 + s*(27), 24 + num_species*(27), E_field_fac*(0.4*rho[s][3])); 
 
    gkyl_mat_set(&lhs, 16 + s*(27), 24 + s*(27), B_field_fac*(0.4*tot_Bx[3])); 
    gkyl_mat_set(&lhs, 25 + s*(27), 15 + s*(27), -B_field_fac*(0.4*tot_Bx[3])); 
 
    gkyl_mat_set(&lhs, 25 + s*(27), 6 + s*(27), B_field_fac*(0.4*tot_By[3])); 
    gkyl_mat_set(&lhs, 7 + s*(27), 24 + s*(27), -B_field_fac*(0.4*tot_By[3])); 
 
    gkyl_mat_set(&lhs, 7 + s*(27), 15 + s*(27), B_field_fac*(0.4*tot_Bz[3])); 
    gkyl_mat_set(&lhs, 16 + s*(27), 6 + s*(27), -B_field_fac*(0.4*tot_Bz[3])); 
 
    gkyl_mat_set(&lhs, 7 + num_species*(27), 6 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 16 + num_species*(27), 15 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 25 + num_species*(27), 24 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 7 + s*(27), 7 + s*(27), 1.0); 
    gkyl_mat_set(&lhs, 16 + s*(27), 16 + s*(27), 1.0); 
    gkyl_mat_set(&lhs, 25 + s*(27), 25 + s*(27), 1.0); 
 
    gkyl_mat_set(&lhs, 7 + s*(27), 7 + num_species*(27), E_field_fac*(0.2857142857142857*rho[s][8]+0.31943828249997*rho[s][5]+0.4472135954999579*rho[s][4]+0.5*rho[s][0])); 
    gkyl_mat_set(&lhs, 16 + s*(27), 16 + num_species*(27), E_field_fac*(0.2857142857142857*rho[s][8]+0.31943828249997*rho[s][5]+0.4472135954999579*rho[s][4]+0.5*rho[s][0])); 
    gkyl_mat_set(&lhs, 25 + s*(27), 25 + num_species*(27), E_field_fac*(0.2857142857142857*rho[s][8]+0.31943828249997*rho[s][5]+0.4472135954999579*rho[s][4]+0.5*rho[s][0])); 
 
    gkyl_mat_set(&lhs, 16 + s*(27), 25 + s*(27), B_field_fac*(0.2857142857142857*tot_Bx[8]+0.31943828249997*tot_Bx[5]+0.4472135954999579*tot_Bx[4]+0.5*tot_Bx[0])); 
    gkyl_mat_set(&lhs, 25 + s*(27), 16 + s*(27), -B_field_fac*(0.2857142857142857*tot_Bx[8]+0.31943828249997*tot_Bx[5]+0.4472135954999579*tot_Bx[4]+0.5*tot_Bx[0])); 
 
    gkyl_mat_set(&lhs, 25 + s*(27), 7 + s*(27), B_field_fac*(0.2857142857142857*tot_By[8]+0.31943828249997*tot_By[5]+0.4472135954999579*tot_By[4]+0.5*tot_By[0])); 
    gkyl_mat_set(&lhs, 7 + s*(27), 25 + s*(27), -B_field_fac*(0.2857142857142857*tot_By[8]+0.31943828249997*tot_By[5]+0.4472135954999579*tot_By[4]+0.5*tot_By[0])); 
 
    gkyl_mat_set(&lhs, 7 + s*(27), 16 + s*(27), B_field_fac*(0.2857142857142857*tot_Bz[8]+0.31943828249997*tot_Bz[5]+0.4472135954999579*tot_Bz[4]+0.5*tot_Bz[0])); 
    gkyl_mat_set(&lhs, 16 + s*(27), 7 + s*(27), -B_field_fac*(0.2857142857142857*tot_Bz[8]+0.31943828249997*tot_Bz[5]+0.4472135954999579*tot_Bz[4]+0.5*tot_Bz[0])); 
 
    gkyl_mat_set(&lhs, 7 + num_species*(27), 7 + s*(27), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 16 + num_species*(27), 16 + s*(27), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 25 + num_species*(27), 25 + s*(27), 0.5*dt*(1.0)); 
 
    gkyl_mat_set(&lhs, 7 + s*(27), 8 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 16 + s*(27), 17 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 25 + s*(27), 26 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 7 + s*(27), 8 + num_species*(27), E_field_fac*(0.2857142857142857*rho[s][7]+0.447213595499958*rho[s][1])); 
    gkyl_mat_set(&lhs, 16 + s*(27), 17 + num_species*(27), E_field_fac*(0.2857142857142857*rho[s][7]+0.447213595499958*rho[s][1])); 
    gkyl_mat_set(&lhs, 25 + s*(27), 26 + num_species*(27), E_field_fac*(0.2857142857142857*rho[s][7]+0.447213595499958*rho[s][1])); 
 
    gkyl_mat_set(&lhs, 16 + s*(27), 26 + s*(27), B_field_fac*(0.2857142857142857*tot_Bx[7]+0.447213595499958*tot_Bx[1])); 
    gkyl_mat_set(&lhs, 25 + s*(27), 17 + s*(27), -B_field_fac*(0.2857142857142857*tot_Bx[7]+0.447213595499958*tot_Bx[1])); 
 
    gkyl_mat_set(&lhs, 25 + s*(27), 8 + s*(27), B_field_fac*(0.2857142857142857*tot_By[7]+0.447213595499958*tot_By[1])); 
    gkyl_mat_set(&lhs, 7 + s*(27), 26 + s*(27), -B_field_fac*(0.2857142857142857*tot_By[7]+0.447213595499958*tot_By[1])); 
 
    gkyl_mat_set(&lhs, 7 + s*(27), 17 + s*(27), B_field_fac*(0.2857142857142857*tot_Bz[7]+0.447213595499958*tot_Bz[1])); 
    gkyl_mat_set(&lhs, 16 + s*(27), 8 + s*(27), -B_field_fac*(0.2857142857142857*tot_Bz[7]+0.447213595499958*tot_Bz[1])); 
 
    gkyl_mat_set(&lhs, 7 + num_species*(27), 8 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 16 + num_species*(27), 17 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 25 + num_species*(27), 26 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 8 + s*(27), 0 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 17 + s*(27), 9 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 26 + s*(27), 18 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 8 + s*(27), 0 + num_species*(27), E_field_fac*(0.5*rho[s][8])); 
    gkyl_mat_set(&lhs, 17 + s*(27), 9 + num_species*(27), E_field_fac*(0.5*rho[s][8])); 
    gkyl_mat_set(&lhs, 26 + s*(27), 18 + num_species*(27), E_field_fac*(0.5*rho[s][8])); 
 
    gkyl_mat_set(&lhs, 17 + s*(27), 18 + s*(27), B_field_fac*(0.5*tot_Bx[8])); 
    gkyl_mat_set(&lhs, 26 + s*(27), 9 + s*(27), -B_field_fac*(0.5*tot_Bx[8])); 
 
    gkyl_mat_set(&lhs, 26 + s*(27), 0 + s*(27), B_field_fac*(0.5*tot_By[8])); 
    gkyl_mat_set(&lhs, 8 + s*(27), 18 + s*(27), -B_field_fac*(0.5*tot_By[8])); 
 
    gkyl_mat_set(&lhs, 8 + s*(27), 9 + s*(27), B_field_fac*(0.5*tot_Bz[8])); 
    gkyl_mat_set(&lhs, 17 + s*(27), 0 + s*(27), -B_field_fac*(0.5*tot_Bz[8])); 
 
    gkyl_mat_set(&lhs, 8 + num_species*(27), 0 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 17 + num_species*(27), 9 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 26 + num_species*(27), 18 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 8 + s*(27), 1 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 17 + s*(27), 10 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 26 + s*(27), 19 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 8 + s*(27), 1 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][7])); 
    gkyl_mat_set(&lhs, 17 + s*(27), 10 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][7])); 
    gkyl_mat_set(&lhs, 26 + s*(27), 19 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][7])); 
 
    gkyl_mat_set(&lhs, 17 + s*(27), 19 + s*(27), B_field_fac*(0.447213595499958*tot_Bx[7])); 
    gkyl_mat_set(&lhs, 26 + s*(27), 10 + s*(27), -B_field_fac*(0.447213595499958*tot_Bx[7])); 
 
    gkyl_mat_set(&lhs, 26 + s*(27), 1 + s*(27), B_field_fac*(0.447213595499958*tot_By[7])); 
    gkyl_mat_set(&lhs, 8 + s*(27), 19 + s*(27), -B_field_fac*(0.447213595499958*tot_By[7])); 
 
    gkyl_mat_set(&lhs, 8 + s*(27), 10 + s*(27), B_field_fac*(0.447213595499958*tot_Bz[7])); 
    gkyl_mat_set(&lhs, 17 + s*(27), 1 + s*(27), -B_field_fac*(0.447213595499958*tot_Bz[7])); 
 
    gkyl_mat_set(&lhs, 8 + num_species*(27), 1 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 17 + num_species*(27), 10 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 26 + num_species*(27), 19 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 8 + s*(27), 2 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 17 + s*(27), 11 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 26 + s*(27), 20 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 8 + s*(27), 2 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][6])); 
    gkyl_mat_set(&lhs, 17 + s*(27), 11 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][6])); 
    gkyl_mat_set(&lhs, 26 + s*(27), 20 + num_species*(27), E_field_fac*(0.447213595499958*rho[s][6])); 
 
    gkyl_mat_set(&lhs, 17 + s*(27), 20 + s*(27), B_field_fac*(0.447213595499958*tot_Bx[6])); 
    gkyl_mat_set(&lhs, 26 + s*(27), 11 + s*(27), -B_field_fac*(0.447213595499958*tot_Bx[6])); 
 
    gkyl_mat_set(&lhs, 26 + s*(27), 2 + s*(27), B_field_fac*(0.447213595499958*tot_By[6])); 
    gkyl_mat_set(&lhs, 8 + s*(27), 20 + s*(27), -B_field_fac*(0.447213595499958*tot_By[6])); 
 
    gkyl_mat_set(&lhs, 8 + s*(27), 11 + s*(27), B_field_fac*(0.447213595499958*tot_Bz[6])); 
    gkyl_mat_set(&lhs, 17 + s*(27), 2 + s*(27), -B_field_fac*(0.447213595499958*tot_Bz[6])); 
 
    gkyl_mat_set(&lhs, 8 + num_species*(27), 2 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 17 + num_species*(27), 11 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 26 + num_species*(27), 20 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 8 + s*(27), 3 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 17 + s*(27), 12 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 26 + s*(27), 21 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 8 + s*(27), 3 + num_species*(27), E_field_fac*(0.4*rho[s][3])); 
    gkyl_mat_set(&lhs, 17 + s*(27), 12 + num_species*(27), E_field_fac*(0.4*rho[s][3])); 
    gkyl_mat_set(&lhs, 26 + s*(27), 21 + num_species*(27), E_field_fac*(0.4*rho[s][3])); 
 
    gkyl_mat_set(&lhs, 17 + s*(27), 21 + s*(27), B_field_fac*(0.4*tot_Bx[3])); 
    gkyl_mat_set(&lhs, 26 + s*(27), 12 + s*(27), -B_field_fac*(0.4*tot_Bx[3])); 
 
    gkyl_mat_set(&lhs, 26 + s*(27), 3 + s*(27), B_field_fac*(0.4*tot_By[3])); 
    gkyl_mat_set(&lhs, 8 + s*(27), 21 + s*(27), -B_field_fac*(0.4*tot_By[3])); 
 
    gkyl_mat_set(&lhs, 8 + s*(27), 12 + s*(27), B_field_fac*(0.4*tot_Bz[3])); 
    gkyl_mat_set(&lhs, 17 + s*(27), 3 + s*(27), -B_field_fac*(0.4*tot_Bz[3])); 
 
    gkyl_mat_set(&lhs, 8 + num_species*(27), 3 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 17 + num_species*(27), 12 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 26 + num_species*(27), 21 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 8 + s*(27), 4 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 17 + s*(27), 13 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 26 + s*(27), 22 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 8 + s*(27), 4 + num_species*(27), E_field_fac*(0.31943828249997*rho[s][8]+0.5*rho[s][5])); 
    gkyl_mat_set(&lhs, 17 + s*(27), 13 + num_species*(27), E_field_fac*(0.31943828249997*rho[s][8]+0.5*rho[s][5])); 
    gkyl_mat_set(&lhs, 26 + s*(27), 22 + num_species*(27), E_field_fac*(0.31943828249997*rho[s][8]+0.5*rho[s][5])); 
 
    gkyl_mat_set(&lhs, 17 + s*(27), 22 + s*(27), B_field_fac*(0.31943828249997*tot_Bx[8]+0.5*tot_Bx[5])); 
    gkyl_mat_set(&lhs, 26 + s*(27), 13 + s*(27), -B_field_fac*(0.31943828249997*tot_Bx[8]+0.5*tot_Bx[5])); 
 
    gkyl_mat_set(&lhs, 26 + s*(27), 4 + s*(27), B_field_fac*(0.31943828249997*tot_By[8]+0.5*tot_By[5])); 
    gkyl_mat_set(&lhs, 8 + s*(27), 22 + s*(27), -B_field_fac*(0.31943828249997*tot_By[8]+0.5*tot_By[5])); 
 
    gkyl_mat_set(&lhs, 8 + s*(27), 13 + s*(27), B_field_fac*(0.31943828249997*tot_Bz[8]+0.5*tot_Bz[5])); 
    gkyl_mat_set(&lhs, 17 + s*(27), 4 + s*(27), -B_field_fac*(0.31943828249997*tot_Bz[8]+0.5*tot_Bz[5])); 
 
    gkyl_mat_set(&lhs, 8 + num_species*(27), 4 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 17 + num_species*(27), 13 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 26 + num_species*(27), 22 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 8 + s*(27), 5 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 17 + s*(27), 14 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 26 + s*(27), 23 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 8 + s*(27), 5 + num_species*(27), E_field_fac*(0.31943828249997*rho[s][8]+0.5*rho[s][4])); 
    gkyl_mat_set(&lhs, 17 + s*(27), 14 + num_species*(27), E_field_fac*(0.31943828249997*rho[s][8]+0.5*rho[s][4])); 
    gkyl_mat_set(&lhs, 26 + s*(27), 23 + num_species*(27), E_field_fac*(0.31943828249997*rho[s][8]+0.5*rho[s][4])); 
 
    gkyl_mat_set(&lhs, 17 + s*(27), 23 + s*(27), B_field_fac*(0.31943828249997*tot_Bx[8]+0.5*tot_Bx[4])); 
    gkyl_mat_set(&lhs, 26 + s*(27), 14 + s*(27), -B_field_fac*(0.31943828249997*tot_Bx[8]+0.5*tot_Bx[4])); 
 
    gkyl_mat_set(&lhs, 26 + s*(27), 5 + s*(27), B_field_fac*(0.31943828249997*tot_By[8]+0.5*tot_By[4])); 
    gkyl_mat_set(&lhs, 8 + s*(27), 23 + s*(27), -B_field_fac*(0.31943828249997*tot_By[8]+0.5*tot_By[4])); 
 
    gkyl_mat_set(&lhs, 8 + s*(27), 14 + s*(27), B_field_fac*(0.31943828249997*tot_Bz[8]+0.5*tot_Bz[4])); 
    gkyl_mat_set(&lhs, 17 + s*(27), 5 + s*(27), -B_field_fac*(0.31943828249997*tot_Bz[8]+0.5*tot_Bz[4])); 
 
    gkyl_mat_set(&lhs, 8 + num_species*(27), 5 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 17 + num_species*(27), 14 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 26 + num_species*(27), 23 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 8 + s*(27), 6 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 17 + s*(27), 15 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 26 + s*(27), 24 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 8 + s*(27), 6 + num_species*(27), E_field_fac*(0.2857142857142857*rho[s][6]+0.447213595499958*rho[s][2])); 
    gkyl_mat_set(&lhs, 17 + s*(27), 15 + num_species*(27), E_field_fac*(0.2857142857142857*rho[s][6]+0.447213595499958*rho[s][2])); 
    gkyl_mat_set(&lhs, 26 + s*(27), 24 + num_species*(27), E_field_fac*(0.2857142857142857*rho[s][6]+0.447213595499958*rho[s][2])); 
 
    gkyl_mat_set(&lhs, 17 + s*(27), 24 + s*(27), B_field_fac*(0.2857142857142857*tot_Bx[6]+0.447213595499958*tot_Bx[2])); 
    gkyl_mat_set(&lhs, 26 + s*(27), 15 + s*(27), -B_field_fac*(0.2857142857142857*tot_Bx[6]+0.447213595499958*tot_Bx[2])); 
 
    gkyl_mat_set(&lhs, 26 + s*(27), 6 + s*(27), B_field_fac*(0.2857142857142857*tot_By[6]+0.447213595499958*tot_By[2])); 
    gkyl_mat_set(&lhs, 8 + s*(27), 24 + s*(27), -B_field_fac*(0.2857142857142857*tot_By[6]+0.447213595499958*tot_By[2])); 
 
    gkyl_mat_set(&lhs, 8 + s*(27), 15 + s*(27), B_field_fac*(0.2857142857142857*tot_Bz[6]+0.447213595499958*tot_Bz[2])); 
    gkyl_mat_set(&lhs, 17 + s*(27), 6 + s*(27), -B_field_fac*(0.2857142857142857*tot_Bz[6]+0.447213595499958*tot_Bz[2])); 
 
    gkyl_mat_set(&lhs, 8 + num_species*(27), 6 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 17 + num_species*(27), 15 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 26 + num_species*(27), 24 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 8 + s*(27), 7 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 17 + s*(27), 16 + s*(27), 0.0); 
    gkyl_mat_set(&lhs, 26 + s*(27), 25 + s*(27), 0.0); 
 
    gkyl_mat_set(&lhs, 8 + s*(27), 7 + num_species*(27), E_field_fac*(0.2857142857142857*rho[s][7]+0.447213595499958*rho[s][1])); 
    gkyl_mat_set(&lhs, 17 + s*(27), 16 + num_species*(27), E_field_fac*(0.2857142857142857*rho[s][7]+0.447213595499958*rho[s][1])); 
    gkyl_mat_set(&lhs, 26 + s*(27), 25 + num_species*(27), E_field_fac*(0.2857142857142857*rho[s][7]+0.447213595499958*rho[s][1])); 
 
    gkyl_mat_set(&lhs, 17 + s*(27), 25 + s*(27), B_field_fac*(0.2857142857142857*tot_Bx[7]+0.447213595499958*tot_Bx[1])); 
    gkyl_mat_set(&lhs, 26 + s*(27), 16 + s*(27), -B_field_fac*(0.2857142857142857*tot_Bx[7]+0.447213595499958*tot_Bx[1])); 
 
    gkyl_mat_set(&lhs, 26 + s*(27), 7 + s*(27), B_field_fac*(0.2857142857142857*tot_By[7]+0.447213595499958*tot_By[1])); 
    gkyl_mat_set(&lhs, 8 + s*(27), 25 + s*(27), -B_field_fac*(0.2857142857142857*tot_By[7]+0.447213595499958*tot_By[1])); 
 
    gkyl_mat_set(&lhs, 8 + s*(27), 16 + s*(27), B_field_fac*(0.2857142857142857*tot_Bz[7]+0.447213595499958*tot_Bz[1])); 
    gkyl_mat_set(&lhs, 17 + s*(27), 7 + s*(27), -B_field_fac*(0.2857142857142857*tot_Bz[7]+0.447213595499958*tot_Bz[1])); 
 
    gkyl_mat_set(&lhs, 8 + num_species*(27), 7 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 17 + num_species*(27), 16 + s*(27), 0.5*dt*(0.0)); 
    gkyl_mat_set(&lhs, 26 + num_species*(27), 25 + s*(27), 0.5*dt*(0.0)); 
 
    gkyl_mat_set(&lhs, 8 + s*(27), 8 + s*(27), 1.0); 
    gkyl_mat_set(&lhs, 17 + s*(27), 17 + s*(27), 1.0); 
    gkyl_mat_set(&lhs, 26 + s*(27), 26 + s*(27), 1.0); 
 
    gkyl_mat_set(&lhs, 8 + s*(27), 8 + num_species*(27), E_field_fac*(0.2040816326530612*rho[s][8]+0.31943828249997*rho[s][5]+0.31943828249997*rho[s][4]+0.5*rho[s][0])); 
    gkyl_mat_set(&lhs, 17 + s*(27), 17 + num_species*(27), E_field_fac*(0.2040816326530612*rho[s][8]+0.31943828249997*rho[s][5]+0.31943828249997*rho[s][4]+0.5*rho[s][0])); 
    gkyl_mat_set(&lhs, 26 + s*(27), 26 + num_species*(27), E_field_fac*(0.2040816326530612*rho[s][8]+0.31943828249997*rho[s][5]+0.31943828249997*rho[s][4]+0.5*rho[s][0])); 
 
    gkyl_mat_set(&lhs, 17 + s*(27), 26 + s*(27), B_field_fac*(0.2040816326530612*tot_Bx[8]+0.31943828249997*tot_Bx[5]+0.31943828249997*tot_Bx[4]+0.5*tot_Bx[0])); 
    gkyl_mat_set(&lhs, 26 + s*(27), 17 + s*(27), -B_field_fac*(0.2040816326530612*tot_Bx[8]+0.31943828249997*tot_Bx[5]+0.31943828249997*tot_Bx[4]+0.5*tot_Bx[0])); 
 
    gkyl_mat_set(&lhs, 26 + s*(27), 8 + s*(27), B_field_fac*(0.2040816326530612*tot_By[8]+0.31943828249997*tot_By[5]+0.31943828249997*tot_By[4]+0.5*tot_By[0])); 
    gkyl_mat_set(&lhs, 8 + s*(27), 26 + s*(27), -B_field_fac*(0.2040816326530612*tot_By[8]+0.31943828249997*tot_By[5]+0.31943828249997*tot_By[4]+0.5*tot_By[0])); 
 
    gkyl_mat_set(&lhs, 8 + s*(27), 17 + s*(27), B_field_fac*(0.2040816326530612*tot_Bz[8]+0.31943828249997*tot_Bz[5]+0.31943828249997*tot_Bz[4]+0.5*tot_Bz[0])); 
    gkyl_mat_set(&lhs, 17 + s*(27), 8 + s*(27), -B_field_fac*(0.2040816326530612*tot_Bz[8]+0.31943828249997*tot_Bz[5]+0.31943828249997*tot_Bz[4]+0.5*tot_Bz[0])); 
 
    gkyl_mat_set(&lhs, 8 + num_species*(27), 8 + s*(27), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 17 + num_species*(27), 17 + s*(27), 0.5*dt*(1.0)); 
    gkyl_mat_set(&lhs, 26 + num_species*(27), 26 + s*(27), 0.5*dt*(1.0)); 
 
  } 
  gkyl_mat_set(&lhs, 0 + num_species*(27), 0 + num_species*(27), 1.0); 
  gkyl_mat_set(&lhs, 9 + num_species*(27), 9 + num_species*(27), 1.0); 
  gkyl_mat_set(&lhs, 18 + num_species*(27), 18 + num_species*(27), 1.0); 
 
  gkyl_mat_set(&lhs, 0 + num_species*(27), 1 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 9 + num_species*(27), 10 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 18 + num_species*(27), 19 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 0 + num_species*(27), 2 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 9 + num_species*(27), 11 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 18 + num_species*(27), 20 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 0 + num_species*(27), 3 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 9 + num_species*(27), 12 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 18 + num_species*(27), 21 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 0 + num_species*(27), 4 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 9 + num_species*(27), 13 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 18 + num_species*(27), 22 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 0 + num_species*(27), 5 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 9 + num_species*(27), 14 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 18 + num_species*(27), 23 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 0 + num_species*(27), 6 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 9 + num_species*(27), 15 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 18 + num_species*(27), 24 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 0 + num_species*(27), 7 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 9 + num_species*(27), 16 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 18 + num_species*(27), 25 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 0 + num_species*(27), 8 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 9 + num_species*(27), 17 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 18 + num_species*(27), 26 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 1 + num_species*(27), 0 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 10 + num_species*(27), 9 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 19 + num_species*(27), 18 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 1 + num_species*(27), 1 + num_species*(27), 1.0); 
  gkyl_mat_set(&lhs, 10 + num_species*(27), 10 + num_species*(27), 1.0); 
  gkyl_mat_set(&lhs, 19 + num_species*(27), 19 + num_species*(27), 1.0); 
 
  gkyl_mat_set(&lhs, 1 + num_species*(27), 2 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 10 + num_species*(27), 11 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 19 + num_species*(27), 20 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 1 + num_species*(27), 3 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 10 + num_species*(27), 12 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 19 + num_species*(27), 21 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 1 + num_species*(27), 4 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 10 + num_species*(27), 13 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 19 + num_species*(27), 22 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 1 + num_species*(27), 5 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 10 + num_species*(27), 14 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 19 + num_species*(27), 23 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 1 + num_species*(27), 6 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 10 + num_species*(27), 15 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 19 + num_species*(27), 24 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 1 + num_species*(27), 7 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 10 + num_species*(27), 16 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 19 + num_species*(27), 25 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 1 + num_species*(27), 8 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 10 + num_species*(27), 17 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 19 + num_species*(27), 26 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 2 + num_species*(27), 0 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 11 + num_species*(27), 9 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 20 + num_species*(27), 18 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 2 + num_species*(27), 1 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 11 + num_species*(27), 10 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 20 + num_species*(27), 19 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 2 + num_species*(27), 2 + num_species*(27), 1.0); 
  gkyl_mat_set(&lhs, 11 + num_species*(27), 11 + num_species*(27), 1.0); 
  gkyl_mat_set(&lhs, 20 + num_species*(27), 20 + num_species*(27), 1.0); 
 
  gkyl_mat_set(&lhs, 2 + num_species*(27), 3 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 11 + num_species*(27), 12 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 20 + num_species*(27), 21 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 2 + num_species*(27), 4 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 11 + num_species*(27), 13 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 20 + num_species*(27), 22 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 2 + num_species*(27), 5 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 11 + num_species*(27), 14 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 20 + num_species*(27), 23 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 2 + num_species*(27), 6 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 11 + num_species*(27), 15 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 20 + num_species*(27), 24 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 2 + num_species*(27), 7 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 11 + num_species*(27), 16 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 20 + num_species*(27), 25 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 2 + num_species*(27), 8 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 11 + num_species*(27), 17 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 20 + num_species*(27), 26 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 3 + num_species*(27), 0 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 12 + num_species*(27), 9 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 21 + num_species*(27), 18 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 3 + num_species*(27), 1 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 12 + num_species*(27), 10 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 21 + num_species*(27), 19 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 3 + num_species*(27), 2 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 12 + num_species*(27), 11 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 21 + num_species*(27), 20 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 3 + num_species*(27), 3 + num_species*(27), 1.0); 
  gkyl_mat_set(&lhs, 12 + num_species*(27), 12 + num_species*(27), 1.0); 
  gkyl_mat_set(&lhs, 21 + num_species*(27), 21 + num_species*(27), 1.0); 
 
  gkyl_mat_set(&lhs, 3 + num_species*(27), 4 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 12 + num_species*(27), 13 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 21 + num_species*(27), 22 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 3 + num_species*(27), 5 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 12 + num_species*(27), 14 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 21 + num_species*(27), 23 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 3 + num_species*(27), 6 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 12 + num_species*(27), 15 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 21 + num_species*(27), 24 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 3 + num_species*(27), 7 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 12 + num_species*(27), 16 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 21 + num_species*(27), 25 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 3 + num_species*(27), 8 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 12 + num_species*(27), 17 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 21 + num_species*(27), 26 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 4 + num_species*(27), 0 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 13 + num_species*(27), 9 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 22 + num_species*(27), 18 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 4 + num_species*(27), 1 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 13 + num_species*(27), 10 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 22 + num_species*(27), 19 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 4 + num_species*(27), 2 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 13 + num_species*(27), 11 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 22 + num_species*(27), 20 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 4 + num_species*(27), 3 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 13 + num_species*(27), 12 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 22 + num_species*(27), 21 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 4 + num_species*(27), 4 + num_species*(27), 1.0); 
  gkyl_mat_set(&lhs, 13 + num_species*(27), 13 + num_species*(27), 1.0); 
  gkyl_mat_set(&lhs, 22 + num_species*(27), 22 + num_species*(27), 1.0); 
 
  gkyl_mat_set(&lhs, 4 + num_species*(27), 5 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 13 + num_species*(27), 14 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 22 + num_species*(27), 23 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 4 + num_species*(27), 6 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 13 + num_species*(27), 15 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 22 + num_species*(27), 24 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 4 + num_species*(27), 7 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 13 + num_species*(27), 16 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 22 + num_species*(27), 25 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 4 + num_species*(27), 8 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 13 + num_species*(27), 17 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 22 + num_species*(27), 26 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 5 + num_species*(27), 0 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 14 + num_species*(27), 9 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 23 + num_species*(27), 18 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 5 + num_species*(27), 1 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 14 + num_species*(27), 10 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 23 + num_species*(27), 19 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 5 + num_species*(27), 2 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 14 + num_species*(27), 11 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 23 + num_species*(27), 20 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 5 + num_species*(27), 3 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 14 + num_species*(27), 12 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 23 + num_species*(27), 21 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 5 + num_species*(27), 4 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 14 + num_species*(27), 13 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 23 + num_species*(27), 22 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 5 + num_species*(27), 5 + num_species*(27), 1.0); 
  gkyl_mat_set(&lhs, 14 + num_species*(27), 14 + num_species*(27), 1.0); 
  gkyl_mat_set(&lhs, 23 + num_species*(27), 23 + num_species*(27), 1.0); 
 
  gkyl_mat_set(&lhs, 5 + num_species*(27), 6 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 14 + num_species*(27), 15 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 23 + num_species*(27), 24 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 5 + num_species*(27), 7 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 14 + num_species*(27), 16 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 23 + num_species*(27), 25 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 5 + num_species*(27), 8 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 14 + num_species*(27), 17 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 23 + num_species*(27), 26 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 6 + num_species*(27), 0 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 15 + num_species*(27), 9 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 24 + num_species*(27), 18 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 6 + num_species*(27), 1 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 15 + num_species*(27), 10 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 24 + num_species*(27), 19 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 6 + num_species*(27), 2 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 15 + num_species*(27), 11 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 24 + num_species*(27), 20 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 6 + num_species*(27), 3 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 15 + num_species*(27), 12 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 24 + num_species*(27), 21 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 6 + num_species*(27), 4 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 15 + num_species*(27), 13 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 24 + num_species*(27), 22 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 6 + num_species*(27), 5 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 15 + num_species*(27), 14 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 24 + num_species*(27), 23 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 6 + num_species*(27), 6 + num_species*(27), 1.0); 
  gkyl_mat_set(&lhs, 15 + num_species*(27), 15 + num_species*(27), 1.0); 
  gkyl_mat_set(&lhs, 24 + num_species*(27), 24 + num_species*(27), 1.0); 
 
  gkyl_mat_set(&lhs, 6 + num_species*(27), 7 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 15 + num_species*(27), 16 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 24 + num_species*(27), 25 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 6 + num_species*(27), 8 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 15 + num_species*(27), 17 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 24 + num_species*(27), 26 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 7 + num_species*(27), 0 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 16 + num_species*(27), 9 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 25 + num_species*(27), 18 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 7 + num_species*(27), 1 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 16 + num_species*(27), 10 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 25 + num_species*(27), 19 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 7 + num_species*(27), 2 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 16 + num_species*(27), 11 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 25 + num_species*(27), 20 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 7 + num_species*(27), 3 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 16 + num_species*(27), 12 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 25 + num_species*(27), 21 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 7 + num_species*(27), 4 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 16 + num_species*(27), 13 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 25 + num_species*(27), 22 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 7 + num_species*(27), 5 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 16 + num_species*(27), 14 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 25 + num_species*(27), 23 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 7 + num_species*(27), 6 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 16 + num_species*(27), 15 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 25 + num_species*(27), 24 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 7 + num_species*(27), 7 + num_species*(27), 1.0); 
  gkyl_mat_set(&lhs, 16 + num_species*(27), 16 + num_species*(27), 1.0); 
  gkyl_mat_set(&lhs, 25 + num_species*(27), 25 + num_species*(27), 1.0); 
 
  gkyl_mat_set(&lhs, 7 + num_species*(27), 8 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 16 + num_species*(27), 17 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 25 + num_species*(27), 26 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 8 + num_species*(27), 0 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 17 + num_species*(27), 9 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 26 + num_species*(27), 18 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 8 + num_species*(27), 1 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 17 + num_species*(27), 10 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 26 + num_species*(27), 19 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 8 + num_species*(27), 2 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 17 + num_species*(27), 11 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 26 + num_species*(27), 20 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 8 + num_species*(27), 3 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 17 + num_species*(27), 12 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 26 + num_species*(27), 21 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 8 + num_species*(27), 4 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 17 + num_species*(27), 13 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 26 + num_species*(27), 22 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 8 + num_species*(27), 5 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 17 + num_species*(27), 14 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 26 + num_species*(27), 23 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 8 + num_species*(27), 6 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 17 + num_species*(27), 15 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 26 + num_species*(27), 24 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 8 + num_species*(27), 7 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 17 + num_species*(27), 16 + num_species*(27), 0.0); 
  gkyl_mat_set(&lhs, 26 + num_species*(27), 25 + num_species*(27), 0.0); 
 
  gkyl_mat_set(&lhs, 8 + num_species*(27), 8 + num_species*(27), 1.0); 
  gkyl_mat_set(&lhs, 17 + num_species*(27), 17 + num_species*(27), 1.0); 
  gkyl_mat_set(&lhs, 26 + num_species*(27), 26 + num_species*(27), 1.0); 
 
} 
