#include <acutest.h>

#include <gkyl_util.h>
#include <gkyl_wv_mhd_rgfm.h>
#include <gkyl_wv_mhd_rgfm_priv.h>

void
test_mhd_rgfm_twospecies_basic()
{
  double gas_gamma1 = 1.4;
  double gas_gamma2 = 1.67;
  double light_speed = 1.0;
  double b_fact = 0.8;

  double *gas_gamma_s = gkyl_malloc(sizeof(double[2]));
  gas_gamma_s[0] = gas_gamma1;
  gas_gamma_s[1] = gas_gamma2;

  struct gkyl_wv_eqn *mhd_rgfm = gkyl_wv_mhd_rgfm_new(2, gas_gamma_s, light_speed, b_fact, 0, 0.0, false);

  TEST_CHECK( mhd_rgfm->num_equations == 13 );
  TEST_CHECK( mhd_rgfm->num_waves == 2 );

  double phi1 = 0.75, rho1 = 1.0, rho2 = 2.0, vx_total = 0.1, vy_total = 0.2, vz_total = 0.3, p_total = 1.5;
  double Bx_total = 0.3, By_total = 0.2, Bz_total = 0.1, psi_total = 0.05;
  
  double rho_total = (phi1 * rho1) + ((1.0 - phi1) * rho2);
  double E1 = (p_total / (gas_gamma1 - 1.0)) + (0.5 * rho1 * ((vx_total * vx_total) + (vy_total * vy_total) + (vz_total * vz_total))) +
    (0.5 * ((Bx_total * Bx_total) + (By_total * By_total) + (Bz_total * Bz_total)));
  double E2 = (p_total / (gas_gamma2 - 1.0)) + (0.5 * rho2 * ((vx_total * vx_total) + (vy_total * vy_total) + (vz_total * vz_total))) +
    (0.5 * ((Bx_total * Bx_total) + (By_total * By_total) + (Bz_total * Bz_total)));
  double E_total = (phi1 * E1) + ((1.0 - phi1) * E2);

  double q[13] = { rho_total, rho_total * vx_total, rho_total * vy_total, rho_total * vz_total, E_total,
    Bx_total, By_total, Bz_total, psi_total, rho_total * phi1, phi1 * rho1, (1.0 - phi1) * rho2, 0.0 };

  double prims[13];
  gkyl_mhd_rgfm_prim_vars(2, gas_gamma_s, q, prims);

  TEST_CHECK( gkyl_compare(prims[0], rho_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[1], vx_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[2], vy_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[3], vz_total, 1e-16) );

  // For now, we check only that the reconstructed interface pressure is of the correct order of magnitude.
  // This error tolerance can be reduced once we have introduced more physical boundary conditions into the system.
  TEST_CHECK( gkyl_compare(prims[4], p_total, 1e-1) );

  TEST_CHECK( gkyl_compare(prims[5], Bx_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[6], By_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[7], Bz_total, 1e-16) );

  TEST_CHECK( gkyl_compare(prims[9], phi1, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[10], rho1, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[11], rho2, 1e-16) );
  
  p_total = prims[4]; // Use the reconstructed interface pressure.

  double fluxes[3][13] = {
   { rho_total * vx_total,
    (rho_total * (vx_total * vx_total)) + (p_total + (0.5 * ((Bx_total * Bx_total) + (By_total * By_total) + (Bz_total * Bz_total)))) - (Bx_total * Bx_total),
    (rho_total * (vx_total * vy_total)) - (Bx_total * By_total),
    (rho_total * (vx_total * vz_total)) - (Bx_total * Bz_total),
    (E_total * vx_total) + (vx_total * (p_total + (0.5 * ((Bx_total * Bx_total) + (By_total * By_total) + (Bz_total * Bz_total))))) -
      (Bx_total * ((vx_total * Bx_total) + (vy_total * By_total) + (vz_total * Bz_total))),
    b_fact * psi_total,
    (By_total * vx_total) - (Bx_total * vy_total),
    (Bz_total * vx_total) - (Bx_total * vz_total),
    b_fact * (light_speed * light_speed) * Bx_total,
    rho_total * vx_total * phi1,
    phi1 * (vx_total * rho1),
    (1.0 - phi1) * (vx_total * rho2),
    0.0 },
   { rho_total * vy_total,
    (rho_total * (vy_total * vx_total)) - (By_total * Bx_total),
    (rho_total * (vy_total * vy_total)) + (p_total + (0.5 * ((Bx_total * Bx_total) + (By_total * By_total) + (Bz_total * Bz_total)))) - (By_total * By_total),
    (rho_total * (vy_total * vz_total)) - (By_total * Bz_total),
    (E_total * vy_total) + (vy_total * (p_total + (0.5 * ((Bx_total * Bx_total) + (By_total * By_total) + (Bz_total * Bz_total))))) -
      (By_total * ((vx_total * Bx_total) + (vy_total * By_total) + (vz_total * Bz_total))),
    (Bx_total * vy_total) - (By_total * vx_total),
    b_fact * psi_total,
    (Bz_total * vy_total) - (By_total * vz_total),
    b_fact * (light_speed * light_speed) * By_total,
    rho_total * vy_total * phi1,
    phi1 * (vy_total * rho1),
    (1.0 - phi1) * (vy_total * rho2),
    0.0 },
   { rho_total * vz_total,
    (rho_total * (vz_total * vx_total)) - (Bz_total * Bx_total),
    (rho_total * (vz_total * vy_total)) - (Bz_total * By_total),
    (rho_total * (vz_total * vz_total)) + (p_total + (0.5 * ((Bx_total * Bx_total) + (By_total * By_total) + (Bz_total * Bz_total)))) - (Bz_total * Bz_total),
    (E_total * vz_total) + (vz_total * (p_total + (0.5 * ((Bx_total * Bx_total) + (By_total * By_total) + (Bz_total * Bz_total))))) -
      (Bz_total * ((vx_total * Bx_total) + (vy_total * By_total) + (vz_total * Bz_total))),
    (Bx_total * vz_total) - (Bz_total * vx_total),
    (By_total * vz_total) - (Bz_total * vy_total),
    b_fact * psi_total,
    b_fact * (light_speed * light_speed) * Bz_total,
    rho_total * vz_total * phi1,
    phi1 * (vz_total * rho1),
    (1.0 - phi1) * (vz_total * rho2),
    0.0 },
  };

  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, 1.0, 0.0 },
    { 0.0, 0.0, 1.0 },
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 },
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, -1.0 },
    { 0.0, 1.0, 0.0 },
  };

  double q_local[13], flux_local[13], flux[13];
  for (int d = 0; d < 3; d++) {
    mhd_rgfm->rotate_to_local_func(mhd_rgfm, tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_mhd_rgfm_flux(2, gas_gamma_s, light_speed, b_fact, q_local, flux_local);
    mhd_rgfm->rotate_to_global_func(mhd_rgfm, tau1[d], tau2[d], norm[d], flux_local, flux);

    for (int i = 0; i < 13; i++) {
      TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-15) );
    }
  }

  double q_l[13], q_g[13];
  for (int d = 0; d < 3; d++) {
    gkyl_wv_eqn_rotate_to_local(mhd_rgfm, tau1[d], tau2[d], norm[d], q, q_l);
    gkyl_wv_eqn_rotate_to_global(mhd_rgfm, tau1[d], tau2[d], norm[d], q_l, q_g);

    for (int i = 0; i < 13; i++) {
      TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-16) );
    }

    double w1[13], q1[13];
    mhd_rgfm->cons_to_riem(mhd_rgfm, q_local, q_local, w1);
    mhd_rgfm->riem_to_cons(mhd_rgfm, q_local, w1, q1);

    for (int i = 0; i < 13; i++) {
      TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-16) );
    }
  }

  gkyl_wv_eqn_release(mhd_rgfm);
  gkyl_free(gas_gamma_s);
}

void
test_mhd_rgfm_threespecies_basic()
{
  double gas_gamma1 = 1.4;
  double gas_gamma2 = 1.67;
  double gas_gamma3 = 1.9;
  double light_speed = 1.0;
  double b_fact = 0.8;

  double *gas_gamma_s = gkyl_malloc(sizeof(double[3]));
  gas_gamma_s[0] = gas_gamma1;
  gas_gamma_s[1] = gas_gamma2;
  gas_gamma_s[2] = gas_gamma3;

  struct gkyl_wv_eqn *mhd_rgfm = gkyl_wv_mhd_rgfm_new(3, gas_gamma_s, light_speed, b_fact, 0, 0.0, false);

  TEST_CHECK( mhd_rgfm->num_equations == 15 );
  TEST_CHECK( mhd_rgfm->num_waves == 2 );

  double phi1 = 0.5, phi2 = 0.3, rho1 = 1.0, rho2 = 2.0, rho3 = 3.0, vx_total = 0.1, vy_total = 0.2, vz_total = 0.3, p_total = 1.5;
  double Bx_total = 0.3, By_total = 0.2, Bz_total = 0.1, psi_total = 0.05;

  double rho_total = (phi1 * rho1) + (phi2 * rho2) + ((1.0 - (phi1 + phi2)) * rho3);
  double E1 = (p_total / (gas_gamma1 - 1.0)) + (0.5 * rho1 * ((vx_total * vx_total) + (vy_total * vy_total) + (vz_total * vz_total))) +
    (0.5 * ((Bx_total * Bx_total) + (By_total * By_total) + (Bz_total * Bz_total)));
  double E2 = (p_total / (gas_gamma2 - 1.0)) + (0.5 * rho2 * ((vx_total * vx_total) + (vy_total * vy_total) + (vz_total * vz_total))) +
    (0.5 * ((Bx_total * Bx_total) + (By_total * By_total) + (Bz_total * Bz_total)));
  double E3 = (p_total / (gas_gamma3 - 1.0)) + (0.5 * rho3 * ((vx_total * vx_total) + (vy_total * vy_total) + (vz_total * vz_total))) +
    (0.5 * ((Bx_total * Bx_total) + (By_total * By_total) + (Bz_total * Bz_total)));
  double E_total = (phi1 * E1) + (phi2 * E2) + ((1.0 - (phi1 + phi2)) * rho3);

  double q[15] = { rho_total, rho_total * vx_total, rho_total * vy_total, rho_total * vz_total, E_total,
    Bx_total, By_total, Bz_total, psi_total, rho_total * phi1, rho_total * phi2, phi1 * rho1, phi2 * rho2, (1.0 - (phi1 + phi2)) * rho3, 0.0 };

  double prims[15];
  gkyl_mhd_rgfm_prim_vars(3, gas_gamma_s, q, prims);

  TEST_CHECK( gkyl_compare(prims[0], rho_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[1], vx_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[2], vy_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[3], vz_total, 1e-16) );

  // For now, we check only that the reconstructed rgfm pressure is of the correct order of magnitude.
  // This error tolerance can be reduced once we have introduced more physical rgfm rules into the system.
  TEST_CHECK( gkyl_compare(prims[4], p_total, 1e-1) ); 

  TEST_CHECK( gkyl_compare(prims[5], Bx_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[6], By_total, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[7], Bz_total, 1e-16) );

  TEST_CHECK( gkyl_compare(prims[9], phi1, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[10], phi2, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[11], rho1, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[12], rho2, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[13], rho3, 1e-16) );
  
  p_total = prims[4]; // Use the reconstructed interface pressure.

  double fluxes[3][15] = {
   { rho_total * vx_total,
    (rho_total * (vx_total * vx_total)) + (p_total + (0.5 * ((Bx_total * Bx_total) + (By_total * By_total) + (Bz_total * Bz_total)))) - (Bx_total * Bx_total),
    (rho_total * (vx_total * vy_total)) - (Bx_total * By_total),
    (rho_total * (vx_total * vz_total)) - (Bx_total * Bz_total),
    (E_total * vx_total) + (vx_total * (p_total + (0.5 * ((Bx_total * Bx_total) + (By_total * By_total) + (Bz_total * Bz_total))))) -
      (Bx_total * ((vx_total * Bx_total) + (vy_total * By_total) + (vz_total * Bz_total))),
    b_fact * psi_total,
    (By_total * vx_total) - (Bx_total * vy_total),
    (Bz_total * vx_total) - (Bx_total * vz_total),
    b_fact * (light_speed * light_speed) * Bx_total,
    rho_total * vx_total * phi1,
    rho_total * vx_total * phi2,
    phi1 * (vx_total * rho1),
    phi2 * (vx_total * rho2),
    (1.0 - (phi1 + phi2)) * (vx_total * rho3),
    0.0 },
   { rho_total * vy_total,
    (rho_total * (vy_total * vx_total)) - (By_total * Bx_total),
    (rho_total * (vy_total * vy_total)) + (p_total + (0.5 * ((Bx_total * Bx_total) + (By_total * By_total) + (Bz_total * Bz_total)))) - (By_total * By_total),
    (rho_total * (vy_total * vz_total)) - (By_total * Bz_total),
    (E_total * vy_total) + (vy_total * (p_total + (0.5 * ((Bx_total * Bx_total) + (By_total * By_total) + (Bz_total * Bz_total))))) -
      (By_total * ((vx_total * Bx_total) + (vy_total * By_total) + (vz_total * Bz_total))),
    (Bx_total * vy_total) - (By_total * vx_total),
    b_fact * psi_total,
    (Bz_total * vy_total) - (By_total * vz_total),
    b_fact * (light_speed * light_speed) * By_total,
    rho_total * vy_total * phi1,
    rho_total * vy_total * phi2,
    phi1 * (vy_total * rho1),
    phi2 * (vy_total * rho2),
    (1.0 - (phi1 + phi2)) * (vy_total * rho3),
    0.0 },
   { rho_total * vz_total,
    (rho_total * (vz_total * vx_total)) - (Bz_total * Bx_total),
    (rho_total * (vz_total * vy_total)) - (Bz_total * By_total),
    (rho_total * (vz_total * vz_total)) + (p_total + (0.5 * ((Bx_total * Bx_total) + (By_total * By_total) + (Bz_total * Bz_total)))) - (Bz_total * Bz_total),
    (E_total * vz_total) + (vz_total * (p_total + (0.5 * ((Bx_total * Bx_total) + (By_total * By_total) + (Bz_total * Bz_total))))) -
      (Bz_total * ((vx_total * Bx_total) + (vy_total * By_total) + (vz_total * Bz_total))),
    (Bx_total * vz_total) - (Bz_total * vx_total),
    (By_total * vz_total) - (Bz_total * vy_total),
    b_fact * psi_total,
    b_fact * (light_speed * light_speed) * Bz_total,
    rho_total * vz_total * phi1,
    rho_total * vz_total * phi2,
    phi1 * (vz_total * rho1),
    phi2 * (vz_total * rho2),
    (1.0 - (phi1 + phi2)) * (vz_total * rho3),
    0.0 },
  };

  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, 1.0, 0.0 },
    { 0.0, 0.0, 1.0 },
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 },
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, -1.0 },
    { 0.0, 1.0, 0.0 },
  };

  double q_local[15], flux_local[15], flux[15];
  for (int d = 0; d < 3; d++) {
    mhd_rgfm->rotate_to_local_func(mhd_rgfm, tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_mhd_rgfm_flux(3, gas_gamma_s, light_speed, b_fact, q_local, flux_local);
    mhd_rgfm->rotate_to_global_func(mhd_rgfm, tau1[d], tau2[d], norm[d], flux_local, flux);

    for (int i = 0; i < 15; i++) {
      TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-15) );
    }
  }

  double q_l[15], q_g[15];
  for (int d = 0; d < 3; d++) {
    gkyl_wv_eqn_rotate_to_local(mhd_rgfm, tau1[d], tau2[d], norm[d], q, q_l);
    gkyl_wv_eqn_rotate_to_global(mhd_rgfm, tau1[d], tau2[d], norm[d], q_l, q_g);

    for (int i = 0; i < 15; i++) {
      TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-16) );
    }

    double w1[15], q1[15];
    mhd_rgfm->cons_to_riem(mhd_rgfm, q_local, q_local, w1);
    mhd_rgfm->riem_to_cons(mhd_rgfm, q_local, w1, q1);

    for (int i = 0; i < 15; i++) {
      TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-16) );
    }
  }

  gkyl_wv_eqn_release(mhd_rgfm);
  gkyl_free(gas_gamma_s);
}

void
test_mhd_rgfm_twospecies_waves()
{
  double gas_gamma1 = 1.4;
  double gas_gamma2 = 1.67;
  double light_speed = 1.0;
  double b_fact = 0.8;
  
  double *gas_gamma_s = gkyl_malloc(sizeof(double[2]));
  gas_gamma_s[0] = gas_gamma1;
  gas_gamma_s[1] = gas_gamma2;

  struct gkyl_wv_eqn *mhd_rgfm = gkyl_wv_mhd_rgfm_new(2, gas_gamma_s, light_speed, b_fact, 0, 0.0, false);
  
  double phi1_l = 0.75, rho1_l = 1.0, rho2_l = 2.0, vx_total_l = 0.1, vy_total_l = 0.2, vz_total_l = 0.3, p_total_l = 1.5;
  double Bx_total_l = 0.3, By_total_l = 0.2, Bz_total_l = 0.1, psi_total_l = 0.05;

  double rho_total_l = (phi1_l * rho1_l) + ((1.0 - phi1_l) * rho2_l);
  double E1_l = (p_total_l / (gas_gamma1 - 1.0)) + (0.5 * rho1_l * ((vx_total_l * vx_total_l) + (vy_total_l * vy_total_l) + (vz_total_l * vz_total_l))) +
    (0.5 * ((Bx_total_l * Bx_total_l) + (By_total_l * By_total_l) + (Bz_total_l * Bz_total_l)));
  double E2_l = (p_total_l / (gas_gamma2 - 1.0)) + (0.5 * rho2_l * ((vx_total_l * vx_total_l) + (vy_total_l * vy_total_l) + (vz_total_l * vz_total_l))) +
    (0.5 * ((Bx_total_l * Bx_total_l) + (By_total_l * By_total_l) + (Bz_total_l * Bz_total_l)));
  double E_total_l = (phi1_l * E1_l) + ((1.0 - phi1_l) * E2_l);

  double phi1_r = 0.25, rho1_r = 0.1, rho2_r = 0.2, vx_total_r = 1.0, vy_total_r = 2.0, vz_total_r = 3.0, p_total_r = 0.15;
  double Bx_total_r = 0.6, By_total_r = 0.4, Bz_total_r = 0.2, psi_total_r = 0.1;

  double rho_total_r = (phi1_r * rho1_r) + ((1.0 - phi1_r) * rho2_r);
  double E1_r = (p_total_r / (gas_gamma1 - 1.0)) + (0.5 * rho1_r * ((vx_total_r * vx_total_r) + (vy_total_r * vy_total_r) + (vz_total_r * vz_total_r))) +
    (0.5 * ((Bx_total_r * Bx_total_r) + (By_total_r * By_total_r) + (Bz_total_r * Bz_total_r)));
  double E2_r = (p_total_r / (gas_gamma2 - 1.0)) + (0.5 * rho2_r * ((vx_total_r * vx_total_r) + (vy_total_r * vy_total_r) + (vz_total_r * vz_total_r))) +
    (0.5 * ((Bx_total_r * Bx_total_r) + (By_total_r * By_total_r) + (Bz_total_r * Bz_total_r)));
  double E_total_r = (phi1_r * E1_r) + ((1.0 - phi1_r) * E2_r);

  double ql[13] = { rho_total_l, rho_total_l * vx_total_l, rho_total_l * vy_total_l, rho_total_l * vz_total_l, E_total_l,
    Bx_total_l, By_total_l, Bz_total_l, psi_total_l, rho_total_l * phi1_l, phi1_l * rho1_l, (1.0 - phi1_l) * rho2_l, 0.0 };
  double qr[13] = { rho_total_r, rho_total_r * vx_total_r, rho_total_r * vy_total_r, rho_total_r * vz_total_r, E_total_r,
    Bx_total_r, By_total_r, Bz_total_r, psi_total_r, rho_total_r * phi1_r, phi1_r * rho1_r, (1.0 - phi1_r) * rho2_r, 0.0 };
  
  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, -1.0, 0.0 },
    { 0.0, 0.0, 1.0 },
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 },
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, 1.0 },
    { 0.0, 1.0, 0.0 }
  };
  
  for (int d = 0; d < 3; d++) {
    double speeds[2], waves[2 * 13], waves_local[2 * 13];

    double ql_local[13], qr_local[13];
    gkyl_wv_eqn_rotate_to_local(mhd_rgfm, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(mhd_rgfm, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[13];
    for (int i = 0; i < 13; i++) {
      delta[i] = qr_local[i] - ql_local[i];
    }

    gkyl_wv_eqn_waves(mhd_rgfm, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    double apdq_local[13], amdq_local[13];
    gkyl_wv_eqn_qfluct(mhd_rgfm, GKYL_WV_HIGH_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

    for (int i = 0; i < 2; i++) {
      gkyl_wv_eqn_rotate_to_global(mhd_rgfm, tau1[d], tau2[d], norm[d], &waves_local[i * 13], &waves[i * 13]);
    }

    double apdq[13], amdq[13];
    gkyl_wv_eqn_rotate_to_global(mhd_rgfm, tau1[d], tau2[d], norm[d], apdq_local, apdq);
    gkyl_wv_eqn_rotate_to_global(mhd_rgfm, tau1[d], tau2[d], norm[d], amdq_local, amdq);

    double fl_local[13], fr_local[13];
    gkyl_mhd_rgfm_flux(2, gas_gamma_s, light_speed, b_fact, ql_local, fl_local);
    gkyl_mhd_rgfm_flux(2, gas_gamma_s, light_speed, b_fact, qr_local, fr_local);

    double fl[13], fr[13];
    gkyl_wv_eqn_rotate_to_global(mhd_rgfm, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(mhd_rgfm, tau1[d], tau2[d], norm[d], fr_local, fr);

    for (int i = 0; i < 13; i++) {
      TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-14) );
    }
  }

  gkyl_wv_eqn_release(mhd_rgfm);
  gkyl_free(gas_gamma_s);
}

void
test_mhd_rgfm_twospecies_waves_2()
{
  double gas_gamma1 = 1.2;
  double gas_gamma2 = 1.7;
  double light_speed = 1.0;
  double b_fact = 0.8;
  
  double *gas_gamma_s = gkyl_malloc(sizeof(double[2]));
  gas_gamma_s[0] = gas_gamma1;
  gas_gamma_s[1] = gas_gamma2;

  struct gkyl_wv_eqn *mhd_rgfm = gkyl_wv_mhd_rgfm_new(2, gas_gamma_s, light_speed, b_fact, 0, 0.0, false);
  
  double phi1_l = 0.9, rho1_l = 1.0, rho2_l = 2.0, vx_total_l = 0.1, vy_total_l = 0.2, vz_total_l = 0.3, p_total_l = 1.5;
  double Bx_total_l = 3.0, By_total_l = 2.0, Bz_total_l = 1.0, psi_total_l = 0.05;

  double rho_total_l = (phi1_l * rho1_l) + ((1.0 - phi1_l) * rho2_l);
  double E1_l = (p_total_l / (gas_gamma1 - 1.0)) + (0.5 * rho1_l * ((vx_total_l * vx_total_l) + (vy_total_l * vy_total_l) + (vz_total_l * vz_total_l))) +
    (0.5 * ((Bx_total_l * Bx_total_l) + (By_total_l * By_total_l) + (Bz_total_l * Bz_total_l)));
  double E2_l = (p_total_l / (gas_gamma2 - 1.0)) + (0.5 * rho2_l * ((vx_total_l * vx_total_l) + (vy_total_l * vy_total_l) + (vz_total_l * vz_total_l))) +
    (0.5 * ((Bx_total_l * Bx_total_l) + (By_total_l * By_total_l) + (Bz_total_l * Bz_total_l)));
  double E_total_l = (phi1_l * E1_l) + ((1.0 - phi1_l) * E2_l);

  double phi1_r = 0.1, rho1_r = 0.01, rho2_r = 0.02, vx_total_r = 1.0, vy_total_r = 2.0, vz_total_r = 3.0, p_total_r = 15.0;
  double Bx_total_r = 0.3, By_total_r = 0.2, Bz_total_r = 0.1, psi_total_r = 0.1;

  double rho_total_r = (phi1_r * rho1_r) + ((1.0 - phi1_r) * rho2_r);
  double E1_r = (p_total_r / (gas_gamma1 - 1.0)) + (0.5 * rho1_r * ((vx_total_r * vx_total_r) + (vy_total_r * vy_total_r) + (vz_total_r * vz_total_r))) +
    (0.5 * ((Bx_total_r * Bx_total_r) + (By_total_r * By_total_r) + (Bz_total_r * Bz_total_r)));
  double E2_r = (p_total_r / (gas_gamma2 - 1.0)) + (0.5 * rho2_r * ((vx_total_r * vx_total_r) + (vy_total_r * vy_total_r) + (vz_total_r * vz_total_r))) +
    (0.5 * ((Bx_total_r * Bx_total_r) + (By_total_r * By_total_r) + (Bz_total_r * Bz_total_r)));
  double E_total_r = (phi1_r * E1_r) + ((1.0 - phi1_r) * E2_r);

  double ql[13] = { rho_total_l, rho_total_l * vx_total_l, rho_total_l * vy_total_l, rho_total_l * vz_total_l, E_total_l,
    Bx_total_l, By_total_l, Bz_total_l, psi_total_l, rho_total_l * phi1_l, phi1_l * rho1_l, (1.0 - phi1_l) * rho2_l, 0.0 };
  double qr[13] = { rho_total_r, rho_total_r * vx_total_r, rho_total_r * vy_total_r, rho_total_r * vz_total_r, E_total_r,
    Bx_total_r, By_total_r, Bz_total_r, psi_total_r, rho_total_r * phi1_r, phi1_r * rho1_r, (1.0 - phi1_r) * rho2_r, 0.0 };
  
  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, -1.0, 0.0 },
    { 0.0, 0.0, 1.0 },
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 },
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, 1.0 },
    { 0.0, 1.0, 0.0 }
  };
  
  for (int d = 0; d < 3; d++) {
    double speeds[2], waves[2 * 13], waves_local[2 * 13];

    double ql_local[13], qr_local[13];
    gkyl_wv_eqn_rotate_to_local(mhd_rgfm, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(mhd_rgfm, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[13];
    for (int i = 0; i < 13; i++) {
      delta[i] = qr_local[i] - ql_local[i];
    }

    gkyl_wv_eqn_waves(mhd_rgfm, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    double apdq_local[13], amdq_local[13];
    gkyl_wv_eqn_qfluct(mhd_rgfm, GKYL_WV_HIGH_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

    for (int i = 0; i < 2; i++) {
      gkyl_wv_eqn_rotate_to_global(mhd_rgfm, tau1[d], tau2[d], norm[d], &waves_local[i * 13], &waves[i * 13]);
    }

    double apdq[13], amdq[13];
    gkyl_wv_eqn_rotate_to_global(mhd_rgfm, tau1[d], tau2[d], norm[d], apdq_local, apdq);
    gkyl_wv_eqn_rotate_to_global(mhd_rgfm, tau1[d], tau2[d], norm[d], amdq_local, amdq);

    double fl_local[13], fr_local[13];
    gkyl_mhd_rgfm_flux(2, gas_gamma_s, light_speed, b_fact, ql_local, fl_local);
    gkyl_mhd_rgfm_flux(2, gas_gamma_s, light_speed, b_fact, qr_local, fr_local);

    double fl[13], fr[13];
    gkyl_wv_eqn_rotate_to_global(mhd_rgfm, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(mhd_rgfm, tau1[d], tau2[d], norm[d], fr_local, fr);

    for (int i = 0; i < 13; i++) {
      TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-12) );
    }
  }

  gkyl_wv_eqn_release(mhd_rgfm);
  gkyl_free(gas_gamma_s);
}

void
test_mhd_rgfm_threespecies_waves()
{
  double gas_gamma1 = 1.4;
  double gas_gamma2 = 1.67;
  double gas_gamma3 = 1.9;
  double light_speed = 1.0;
  double b_fact = 0.8;
  
  double *gas_gamma_s = gkyl_malloc(sizeof(double[3]));
  gas_gamma_s[0] = gas_gamma1;
  gas_gamma_s[1] = gas_gamma2;
  gas_gamma_s[2] = gas_gamma3;

  struct gkyl_wv_eqn *mhd_rgfm = gkyl_wv_mhd_rgfm_new(3, gas_gamma_s, light_speed, b_fact, 0, 0.0, false);
  
  double phi1_l = 0.5, phi2_l = 0.3, rho1_l = 1.0, rho2_l = 2.0, rho3_l = 3.0, vx_total_l = 0.1, vy_total_l = 0.2, vz_total_l = 0.3, p_total_l = 1.5;
  double Bx_total_l = 0.3, By_total_l = 0.2, Bz_total_l = 0.1, psi_total_l = 0.05;

  double rho_total_l = (phi1_l * rho1_l) + (phi2_l * rho2_l) + ((1.0 - (phi1_l + phi2_l)) * rho3_l);
  double E1_l = (p_total_l / (gas_gamma1 - 1.0)) + (0.5 * rho1_l * ((vx_total_l * vx_total_l) + (vy_total_l * vy_total_l) + (vz_total_l * vz_total_l))) +
    (0.5 * ((Bx_total_l * Bx_total_l) + (By_total_l * By_total_l) + (Bz_total_l * Bz_total_l)));
  double E2_l = (p_total_l / (gas_gamma2 - 1.0)) + (0.5 * rho2_l * ((vx_total_l * vx_total_l) + (vy_total_l * vy_total_l) + (vz_total_l * vz_total_l))) +
    (0.5 * ((Bx_total_l * Bx_total_l) + (By_total_l * By_total_l) + (Bz_total_l * Bz_total_l)));
  double E3_l = (p_total_l / (gas_gamma3 - 1.0)) + (0.5 * rho3_l * ((vx_total_l * vx_total_l) + (vy_total_l * vy_total_l) + (vz_total_l * vz_total_l))) +
    (0.5 * ((Bx_total_l * Bx_total_l) + (By_total_l * By_total_l) + (Bz_total_l * Bz_total_l)));
  double E_total_l = (phi1_l * E1_l) + (phi2_l * E2_l) + ((1.0 - (phi1_l + phi2_l)) * E3_l);

  double phi1_r = 0.4, phi2_r = 0.2, rho1_r = 0.1, rho2_r = 0.2, rho3_r = 0.3, vx_total_r = 1.0, vy_total_r = 2.0, vz_total_r = 3.0, p_total_r = 0.15;
  double Bx_total_r = 0.6, By_total_r = 0.4, Bz_total_r = 0.2, psi_total_r = 0.1;

  double rho_total_r = (phi1_r * rho1_r) + (phi2_r * rho2_r) + ((1.0 - (phi1_r + phi2_r)) * rho3_r);
  double E1_r = (p_total_r / (gas_gamma1 - 1.0)) + (0.5 * rho1_r * ((vx_total_r * vx_total_r) + (vy_total_r * vy_total_r) + (vz_total_r * vz_total_r))) +
    (0.5 * ((Bx_total_r * Bx_total_r) + (By_total_r * By_total_r) + (Bz_total_r * Bz_total_r)));
  double E2_r = (p_total_r / (gas_gamma2 - 1.0)) + (0.5 * rho2_r * ((vx_total_r * vx_total_r) + (vy_total_r * vy_total_r) + (vz_total_r * vz_total_r))) +
    (0.5 * ((Bx_total_r * Bx_total_r) + (By_total_r * By_total_r) + (Bz_total_r * Bz_total_r)));
  double E3_r = (p_total_r / (gas_gamma3 - 1.0)) + (0.5 * rho3_r * ((vx_total_r * vx_total_r) + (vy_total_r * vy_total_r) + (vz_total_r * vz_total_r))) +
    (0.5 * ((Bx_total_r * Bx_total_r) + (By_total_r * By_total_r) + (Bz_total_r * Bz_total_r)));
  double E_total_r = (phi1_r * E1_r) + (phi2_r * E2_r) + ((1.0 - (phi1_r + phi2_r)) * E3_r);

  double ql[15] = { rho_total_l, rho_total_l * vx_total_l, rho_total_l * vy_total_l, rho_total_l * vz_total_l, E_total_l,
    Bx_total_l, By_total_l, Bz_total_l, psi_total_l, rho_total_l * phi1_l, rho_total_l * phi2_l, phi1_l * rho1_l, phi2_l * rho2_l, (1.0 - (phi1_l + phi2_l)) * rho3_l, 0.0 };
  double qr[15] = { rho_total_r, rho_total_r * vx_total_r, rho_total_r * vy_total_r, rho_total_r * vz_total_r, E_total_r,
    Bx_total_r, By_total_r, Bz_total_r, psi_total_r, rho_total_r * phi1_r, rho_total_r * phi2_r, phi1_r * rho1_r, phi2_r * rho2_r, (1.0 - (phi1_r + phi2_r)) * rho3_r, 0.0 };
  
  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, -1.0, 0.0 },
    { 0.0, 0.0, 1.0 },
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 },
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, 1.0 },
    { 0.0, 1.0, 0.0 }
  };
  
  for (int d = 0; d < 3; d++) {
    double speeds[2], waves[2 * 15], waves_local[2 * 15];

    double ql_local[15], qr_local[15];
    gkyl_wv_eqn_rotate_to_local(mhd_rgfm, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(mhd_rgfm, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[15];
    for (int i = 0; i < 15; i++) {
      delta[i] = qr_local[i] - ql_local[i];
    }

    gkyl_wv_eqn_waves(mhd_rgfm, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    double apdq_local[15], amdq_local[15];
    gkyl_wv_eqn_qfluct(mhd_rgfm, GKYL_WV_HIGH_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

    for (int i = 0; i < 2; i++) {
      gkyl_wv_eqn_rotate_to_global(mhd_rgfm, tau1[d], tau2[d], norm[d], &waves_local[i * 15], &waves[i * 15]);
    }

    double apdq[15], amdq[15];
    gkyl_wv_eqn_rotate_to_global(mhd_rgfm, tau1[d], tau2[d], norm[d], apdq_local, apdq);
    gkyl_wv_eqn_rotate_to_global(mhd_rgfm, tau1[d], tau2[d], norm[d], amdq_local, amdq);

    double fl_local[15], fr_local[15];
    gkyl_mhd_rgfm_flux(3, gas_gamma_s, light_speed, b_fact, ql_local, fl_local);
    gkyl_mhd_rgfm_flux(3, gas_gamma_s, light_speed, b_fact, qr_local, fr_local);

    double fl[15], fr[15];
    gkyl_wv_eqn_rotate_to_global(mhd_rgfm, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(mhd_rgfm, tau1[d], tau2[d], norm[d], fr_local, fr);

    for (int i = 0; i < 15; i++) {
      TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-13) );
    }
  }

  gkyl_wv_eqn_release(mhd_rgfm);
  gkyl_free(gas_gamma_s);
}

void
test_mhd_rgfm_threespecies_waves_2()
{
  double gas_gamma1 = 1.1;
  double gas_gamma2 = 1.5;
  double gas_gamma3 = 1.9;
  double light_speed = 1.0;
  double b_fact = 0.8;
  
  double *gas_gamma_s = gkyl_malloc(sizeof(double[3]));
  gas_gamma_s[0] = gas_gamma1;
  gas_gamma_s[1] = gas_gamma2;
  gas_gamma_s[2] = gas_gamma3;

  struct gkyl_wv_eqn *mhd_rgfm = gkyl_wv_mhd_rgfm_new(3, gas_gamma_s, light_speed, b_fact, 0, 0.0, false);
  
  double phi1_l = 0.8, phi2_l = 0.1, rho1_l = 1.0, rho2_l = 2.0, rho3_l = 3.0, vx_total_l = 0.1, vy_total_l = 0.2, vz_total_l = 0.3, p_total_l = 1.5;
  double Bx_total_l = 3.0, By_total_l = 2.0, Bz_total_l = 1.0, psi_total_l = 0.05;

  double rho_total_l = (phi1_l * rho1_l) + (phi2_l * rho2_l) + ((1.0 - (phi1_l + phi2_l)) * rho3_l);
  double E1_l = (p_total_l / (gas_gamma1 - 1.0)) + (0.5 * rho1_l * ((vx_total_l * vx_total_l) + (vy_total_l * vy_total_l) + (vz_total_l * vz_total_l))) +
    (0.5 * ((Bx_total_l * Bx_total_l) + (By_total_l * By_total_l) + (Bz_total_l * Bz_total_l)));
  double E2_l = (p_total_l / (gas_gamma2 - 1.0)) + (0.5 * rho2_l * ((vx_total_l * vx_total_l) + (vy_total_l * vy_total_l) + (vz_total_l * vz_total_l))) +
    (0.5 * ((Bx_total_l * Bx_total_l) + (By_total_l * By_total_l) + (Bz_total_l * Bz_total_l)));
  double E3_l = (p_total_l / (gas_gamma3 - 1.0)) + (0.5 * rho3_l * ((vx_total_l * vx_total_l) + (vy_total_l * vy_total_l) + (vz_total_l * vz_total_l))) +
    (0.5 * ((Bx_total_l * Bx_total_l) + (By_total_l * By_total_l) + (Bz_total_l * Bz_total_l)));
  double E_total_l = (phi1_l * E1_l) + (phi2_l * E2_l) + ((1.0 - (phi1_l + phi2_l)) * E3_l);

  double phi1_r = 0.75, phi2_r = 0.2, rho1_r = 0.01, rho2_r = 0.02, rho3_r = 0.03, vx_total_r = 1.0, vy_total_r = 2.0, vz_total_r = 3.0, p_total_r = 15.0;
  double Bx_total_r = 0.3, By_total_r = 0.2, Bz_total_r = 0.1, psi_total_r = 0.1;

  double rho_total_r = (phi1_r * rho1_r) + (phi2_r * rho2_r) + ((1.0 - (phi1_r + phi2_r)) * rho3_r);
  double E1_r = (p_total_r / (gas_gamma1 - 1.0)) + (0.5 * rho1_r * ((vx_total_r * vx_total_r) + (vy_total_r * vy_total_r) + (vz_total_r * vz_total_r))) +
    (0.5 * ((Bx_total_r * Bx_total_r) + (By_total_r * By_total_r) + (Bz_total_r * Bz_total_r)));
  double E2_r = (p_total_r / (gas_gamma2 - 1.0)) + (0.5 * rho2_r * ((vx_total_r * vx_total_r) + (vy_total_r * vy_total_r) + (vz_total_r * vz_total_r))) +
    (0.5 * ((Bx_total_r * Bx_total_r) + (By_total_r * By_total_r) + (Bz_total_r * Bz_total_r)));
  double E3_r = (p_total_r / (gas_gamma3 - 1.0)) + (0.5 * rho3_r * ((vx_total_r * vx_total_r) + (vy_total_r * vy_total_r) + (vz_total_r * vz_total_r))) +
    (0.5 * ((Bx_total_r * Bx_total_r) + (By_total_r * By_total_r) + (Bz_total_r * Bz_total_r)));
  double E_total_r = (phi1_r * E1_r) + (phi2_r * E2_r) + ((1.0 - (phi1_r + phi2_r)) * E3_r);

  double ql[15] = { rho_total_l, rho_total_l * vx_total_l, rho_total_l * vy_total_l, rho_total_l * vz_total_l, E_total_l,
    Bx_total_l, By_total_l, Bz_total_l, psi_total_l, rho_total_l * phi1_l, rho_total_l * phi2_l, phi1_l * rho1_l, phi2_l * rho2_l, (1.0 - (phi1_l + phi2_l)) * rho3_l, 0.0 };
  double qr[15] = { rho_total_r, rho_total_r * vx_total_r, rho_total_r * vy_total_r, rho_total_r * vz_total_r, E_total_r,
    Bx_total_r, By_total_r, Bz_total_r, psi_total_r, rho_total_r * phi1_r, rho_total_r * phi2_r, phi1_r * rho1_r, phi2_r * rho2_r, (1.0 - (phi1_r + phi2_r)) * rho3_r, 0.0 };
  
  double norm[3][3] = {
    { 1.0, 0.0, 0.0 },
    { 0.0, -1.0, 0.0 },
    { 0.0, 0.0, 1.0 },
  };

  double tau1[3][3] = {
    { 0.0, 1.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 },
  };

  double tau2[3][3] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, 1.0 },
    { 0.0, 1.0, 0.0 }
  };
  
  for (int d = 0; d < 3; d++) {
    double speeds[2], waves[2 * 15], waves_local[2 * 15];

    double ql_local[15], qr_local[15];
    gkyl_wv_eqn_rotate_to_local(mhd_rgfm, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(mhd_rgfm, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[15];
    for (int i = 0; i < 15; i++) {
      delta[i] = qr_local[i] - ql_local[i];
    }

    gkyl_wv_eqn_waves(mhd_rgfm, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    double apdq_local[15], amdq_local[15];
    gkyl_wv_eqn_qfluct(mhd_rgfm, GKYL_WV_HIGH_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

    for (int i = 0; i < 2; i++) {
      gkyl_wv_eqn_rotate_to_global(mhd_rgfm, tau1[d], tau2[d], norm[d], &waves_local[i * 15], &waves[i * 15]);
    }

    double apdq[15], amdq[15];
    gkyl_wv_eqn_rotate_to_global(mhd_rgfm, tau1[d], tau2[d], norm[d], apdq_local, apdq);
    gkyl_wv_eqn_rotate_to_global(mhd_rgfm, tau1[d], tau2[d], norm[d], amdq_local, amdq);

    double fl_local[15], fr_local[15];
    gkyl_mhd_rgfm_flux(3, gas_gamma_s, light_speed, b_fact, ql_local, fl_local);
    gkyl_mhd_rgfm_flux(3, gas_gamma_s, light_speed, b_fact, qr_local, fr_local);

    double fl[15], fr[15];
    gkyl_wv_eqn_rotate_to_global(mhd_rgfm, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(mhd_rgfm, tau1[d], tau2[d], norm[d], fr_local, fr);

    for (int i = 0; i < 15; i++) {
      TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-12) );
    }
  }

  gkyl_wv_eqn_release(mhd_rgfm);
  gkyl_free(gas_gamma_s);
}

TEST_LIST = {
  { "mhd_rgfm_twospecies_basic", test_mhd_rgfm_twospecies_basic },
  { "mhd_rgfm_threespecies_basic", test_mhd_rgfm_threespecies_basic },
  { "mhd_rgfm_twospecies_waves", test_mhd_rgfm_twospecies_waves },
  { "mhd_rgfm_twospecies_waves_2", test_mhd_rgfm_twospecies_waves_2 },
  { "mhd_rgfm_threespecies_waves", test_mhd_rgfm_threespecies_waves },
  { "mhd_rgfm_threespecies_waves_2", test_mhd_rgfm_threespecies_waves_2 },
  { NULL, NULL },
};