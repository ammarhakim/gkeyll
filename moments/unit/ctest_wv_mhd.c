#include <acutest.h>

#include <gkyl_util.h>
#include <gkyl_wv_mhd.h>
#include <gkyl_wv_mhd_priv.h>

void
test_mhd_basic()
{
  double gas_gamma = 1.4;
  double light_speed = 1.0;
  double b_fact = 0.8;

  struct gkyl_wv_eqn *mhd = gkyl_wv_mhd_new(gas_gamma, light_speed, b_fact, false);

  TEST_CHECK( mhd->num_equations == 9 );
  TEST_CHECK( mhd->num_waves == 2 );

  double rho = 1.0, vx = 0.1, vy = 0.2, vz = 0.3, p = 1.5;
  double Bx = 0.3, By = 0.2, Bz = 0.1, psi = 0.05;
  
  double Etot = (p / (gas_gamma - 1.0)) + (0.5 * rho * ((vx * vx) + (vy * vy) + (vz * vz))) + (0.5 * ((Bx * Bx) + (By * By) + (Bz * Bz)));

  double q[9] = { rho, rho * vx, rho * vy, rho * vz, Etot, Bx, By, Bz, psi };

  double prims[9];
  gkyl_mhd_prim_vars(gas_gamma, q, prims);

  TEST_CHECK( gkyl_compare(prims[0], rho, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[1], vx, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[2], vy, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[3], vz, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[4], p, 1e-16) );

  TEST_CHECK( gkyl_compare(prims[5], Bx, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[6], By, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[7], Bz, 1e-16) );
  TEST_CHECK( gkyl_compare(prims[8], psi, 1e-16) );

  double fluxes[3][9] = {
   { rho * vx,
    (rho * (vx * vx)) + (p + (0.5 * ((Bx * Bx) + (By * By) + (Bz * Bz)))) - (Bx * Bx),
    (rho * (vx * vy)) - (Bx * By),
    (rho * (vx * vz)) - (Bx * Bz),
    (Etot * vx) + (vx * (p + (0.5 * ((Bx * Bx) + (By * By) + (Bz * Bz))))) - (Bx * ((vx * Bx) + (vy * By) + (vz * Bz))),
    b_fact * psi,
    (By * vx) - (Bx * vy),
    (Bz * vx) - (Bx * vz),
    b_fact * (light_speed * light_speed) * Bx },
   { rho * vy,
    (rho * (vy * vx)) - (By * Bx),
    (rho * (vy * vy)) + (p + (0.5 * ((Bx * Bx) + (By * By) + (Bz * Bz)))) - (By * By),
    (rho * (vy * vz)) - (By * Bz),
    (Etot * vy) + (vy * (p + (0.5 * ((Bx * Bx) + (By * By) + (Bz * Bz))))) - (By * ((vx * Bx) + (vy * By) + (vz * Bz))),
    (Bx * vy) - (By * vx),
    b_fact * psi,
    (Bz * vy) - (By * vz),
    b_fact * (light_speed * light_speed) * By },
   { rho * vz,
    (rho * (vz * vx)) - (Bz * Bx),
    (rho * (vz * vy)) - (Bz * By),
    (rho * (vz * vz)) + (p + (0.5 * ((Bx * Bx) + (By * By) + (Bz * Bz)))) - (Bz * Bz),
    (Etot * vz) + (vz * (p + (0.5 * ((Bx * Bx) + (By * By) + (Bz * Bz))))) - (Bz * ((vx * Bx) + (vy * By) + (vz * Bz))),
    (Bx * vz) - (Bz * vx),
    (By * vz) - (Bz * vy),
    b_fact * psi,
    b_fact * (light_speed * light_speed) * Bz },
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

  double q_local[9], flux_local[9], flux[9];
  for (int d = 0; d < 3; d++) {
    mhd->rotate_to_local_func(mhd, tau1[d], tau2[d], norm[d], q, q_local);
    gkyl_mhd_flux(gas_gamma, light_speed, b_fact, q_local, flux_local);
    mhd->rotate_to_global_func(mhd, tau1[d], tau2[d], norm[d], flux_local, flux);

    for (int i = 0; i < 9; i++) {
      TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-15) );
    }
  }

  double q_l[9], q_g[9];
  for (int d = 0; d < 3; d++) {
    gkyl_wv_eqn_rotate_to_local(mhd, tau1[d], tau2[d], norm[d], q, q_l);
    gkyl_wv_eqn_rotate_to_global(mhd, tau1[d], tau2[d], norm[d], q_l, q_g);

    for (int i = 0; i < 9; i++) {
      TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-16) );
    }

    double w1[9], q1[9];
    mhd->cons_to_riem(mhd, q_local, q_local, w1);
    mhd->riem_to_cons(mhd, q_local, w1, q1);

    for (int i = 0; i < 9; i++) {
      TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-16) );
    }
  }

  gkyl_wv_eqn_release(mhd);
}

void
test_mhd_waves()
{
  double gas_gamma = 1.4;
  double light_speed = 1.0;
  double b_fact = 0.8;

  struct gkyl_wv_eqn *mhd = gkyl_wv_mhd_new(gas_gamma, light_speed, b_fact, false);
  
  double rho_l = 1.0, vx_l = 0.1, vy_l = 0.2, vz_l = 0.3, p_l = 1.5;
  double Bx_l = 0.3, By_l = 0.2, Bz_l = 0.1, psi_l = 0.05;

  double Etot_l = (p_l / (gas_gamma - 1.0)) + (0.5 * rho_l * ((vx_l * vx_l) + (vy_l * vy_l) + (vz_l * vz_l))) + (0.5 * ((Bx_l * Bx_l) + (By_l * By_l) + (Bz_l * Bz_l)));

  double rho_r = 0.1, vx_r = 1.0, vy_r = 2.0, vz_r = 3.0, p_r = 0.15;
  double Bx_r = 0.6, By_r = 0.4, Bz_r = 0.2, psi_r = 0.1;

  double Etot_r = (p_r / (gas_gamma - 1.0)) + (0.5 * rho_r * ((vx_r * vx_r) + (vy_r * vy_r) + (vz_r * vz_r))) + (0.5 * ((Bx_r * Bx_r) + (By_r * By_r) + (Bz_r * Bz_r)));

  double ql[9] = { rho_l, rho_l * vx_l, rho_l * vy_l, rho_l * vz_l, Etot_l, Bx_l, By_l, Bz_l, psi_l };
  double qr[9] = { rho_r, rho_r * vx_r, rho_r * vy_r, rho_r * vz_r, Etot_r, Bx_r, By_r, Bz_r, psi_r };
  
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
    double speeds[2], waves[2 * 9], waves_local[2 * 9];

    double ql_local[9], qr_local[9];
    gkyl_wv_eqn_rotate_to_local(mhd, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(mhd, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[9];
    for (int i = 0; i < 9; i++) {
      delta[i] = qr_local[i] - ql_local[i];
    }

    gkyl_wv_eqn_waves(mhd, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    double apdq_local[9], amdq_local[9];
    gkyl_wv_eqn_qfluct(mhd, GKYL_WV_HIGH_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

    for (int i = 0; i < 2; i++) {
      gkyl_wv_eqn_rotate_to_global(mhd, tau1[d], tau2[d], norm[d], &waves_local[i * 13], &waves[i * 13]);
    }

    double apdq[9], amdq[9];
    gkyl_wv_eqn_rotate_to_global(mhd, tau1[d], tau2[d], norm[d], apdq_local, apdq);
    gkyl_wv_eqn_rotate_to_global(mhd, tau1[d], tau2[d], norm[d], amdq_local, amdq);

    double fl_local[9], fr_local[9];
    gkyl_mhd_flux(gas_gamma, light_speed, b_fact, ql_local, fl_local);
    gkyl_mhd_flux(gas_gamma, light_speed, b_fact, qr_local, fr_local);

    double fl[9], fr[9];
    gkyl_wv_eqn_rotate_to_global(mhd, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(mhd, tau1[d], tau2[d], norm[d], fr_local, fr);

    for (int i = 0; i < 9; i++) {
      TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-14) );
    }
  }

  gkyl_wv_eqn_release(mhd);
}

void
test_mhd_waves_2()
{
  double gas_gamma = 1.4;
  double light_speed = 1.0;
  double b_fact = 0.8;

  struct gkyl_wv_eqn *mhd = gkyl_wv_mhd_new(gas_gamma, light_speed, b_fact, false);
  
  double rho_l = 1.0, vx_l = 0.1, vy_l = 0.2, vz_l = 0.3, p_l = 1.5;
  double Bx_l = 3.0, By_l = 2.0, Bz_l = 1.0, psi_l = 0.05;

  double Etot_l = (p_l / (gas_gamma - 1.0)) + (0.5 * rho_l * ((vx_l * vx_l) + (vy_l * vy_l) + (vz_l * vz_l))) + (0.5 * ((Bx_l * Bx_l) + (By_l * By_l) + (Bz_l * Bz_l)));

  double rho_r = 0.01, vx_r = 1.0, vy_r = 2.0, vz_r = 3.0, p_r = 15.0;
  double Bx_r = 0.3, By_r = 0.2, Bz_r = 0.1, psi_r = 0.1;

  double Etot_r = (p_r / (gas_gamma - 1.0)) + (0.5 * rho_r * ((vx_r * vx_r) + (vy_r * vy_r) + (vz_r * vz_r))) + (0.5 * ((Bx_r * Bx_r) + (By_r * By_r) + (Bz_r * Bz_r)));

  double ql[9] = { rho_l, rho_l * vx_l, rho_l * vy_l, rho_l * vz_l, Etot_l, Bx_l, By_l, Bz_l, psi_l };
  double qr[9] = { rho_r, rho_r * vx_r, rho_r * vy_r, rho_r * vz_r, Etot_r, Bx_r, By_r, Bz_r, psi_r };
  
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
    double speeds[2], waves[2 * 9], waves_local[2 * 9];

    double ql_local[9], qr_local[9];
    gkyl_wv_eqn_rotate_to_local(mhd, tau1[d], tau2[d], norm[d], ql, ql_local);
    gkyl_wv_eqn_rotate_to_local(mhd, tau1[d], tau2[d], norm[d], qr, qr_local);

    double delta[9];
    for (int i = 0; i < 9; i++) {
      delta[i] = qr_local[i] - ql_local[i];
    }

    gkyl_wv_eqn_waves(mhd, GKYL_WV_HIGH_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

    double apdq_local[9], amdq_local[9];
    gkyl_wv_eqn_qfluct(mhd, GKYL_WV_HIGH_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

    for (int i = 0; i < 2; i++) {
      gkyl_wv_eqn_rotate_to_global(mhd, tau1[d], tau2[d], norm[d], &waves_local[i * 13], &waves[i * 13]);
    }

    double apdq[9], amdq[9];
    gkyl_wv_eqn_rotate_to_global(mhd, tau1[d], tau2[d], norm[d], apdq_local, apdq);
    gkyl_wv_eqn_rotate_to_global(mhd, tau1[d], tau2[d], norm[d], amdq_local, amdq);

    double fl_local[9], fr_local[9];
    gkyl_mhd_flux(gas_gamma, light_speed, b_fact, ql_local, fl_local);
    gkyl_mhd_flux(gas_gamma, light_speed, b_fact, qr_local, fr_local);

    double fl[9], fr[9];
    gkyl_wv_eqn_rotate_to_global(mhd, tau1[d], tau2[d], norm[d], fl_local, fl);
    gkyl_wv_eqn_rotate_to_global(mhd, tau1[d], tau2[d], norm[d], fr_local, fr);

    for (int i = 0; i < 9; i++) {
      TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-12) );
    }
  }

  gkyl_wv_eqn_release(mhd);
}

TEST_LIST = {
  { "mhd_basic", test_mhd_basic },
  { "mhd_waves", test_mhd_waves },
  { "mhd_waves_2", test_mhd_waves_2 },
  { NULL, NULL },
};