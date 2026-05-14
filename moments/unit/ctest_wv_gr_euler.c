#include <acutest.h>

#include <gkyl_util.h>
#include <gkyl_wv_gr_euler.h>
#include <gkyl_wv_gr_euler_priv.h>
#include <gkyl_gr_minkowski.h>
#include <gkyl_gr_blackhole.h>

#include "prim_vars_stringent_data.h"

void
test_gr_euler_basic_minkowski()
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  struct gkyl_wv_eqn *gr_euler = gkyl_wv_gr_euler_new(gas_gamma, GKYL_STATIC_GAUGE, 0, spacetime, false);

  TEST_CHECK( gr_euler->num_equations == 71 );
  TEST_CHECK( gr_euler->num_waves == 2 );

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double rho = 1.0, u = 0.1, v = 0.2, w = 0.3, p = 1.5;

      double spatial_det, lapse;
      double *shift = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region;

      double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
      }

      double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
      }

      double *lapse_der = gkyl_malloc(sizeof(double[3]));
      double **shift_der = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        shift_der[i] = gkyl_malloc(sizeof(double[3]));
      }

      double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

        for (int j = 0; j < 3; j++) {
          spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
        }
      }

      spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &spatial_det);
      spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse);
      spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift);
      spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);
      
      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

      spacetime->lapse_function_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
      spacetime->shift_vector_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

      double *vel = gkyl_malloc(sizeof(double[3]));
      vel[0] = u; vel[1] = v; vel[2] = w;

      double v_sq = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          v_sq += spatial_metric[i][j] * vel[i] * vel[j];
        }
      }

      double W = 1.0 / sqrt(1.0 - v_sq);
      double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));

      double q[71];
      q[0] = sqrt(spatial_det) * rho * W;
      q[1] = sqrt(spatial_det) * rho * h * (W * W) * u;
      q[2] = sqrt(spatial_det) * rho * h * (W * W) * v;
      q[3] = sqrt(spatial_det) * rho * h * (W * W) * w;
      q[4] = sqrt(spatial_det) * ((rho * h * (W * W)) - p - (rho * W));

      q[5] = lapse;
      q[6] = shift[0]; q[7] = shift[1]; q[8] = shift[2];

      q[9] = spatial_metric[0][0]; q[10] = spatial_metric[0][1]; q[11] = spatial_metric[0][2];
      q[12] = spatial_metric[1][0]; q[13] = spatial_metric[1][1]; q[14] = spatial_metric[1][2];
      q[15] = spatial_metric[2][0]; q[16] = spatial_metric[2][1]; q[17] = spatial_metric[2][2];

      q[18] = extrinsic_curvature[0][0]; q[19] = extrinsic_curvature[0][1]; q[20] = extrinsic_curvature[0][2];
      q[21] = extrinsic_curvature[1][0]; q[22] = extrinsic_curvature[1][1]; q[23] = extrinsic_curvature[1][2];
      q[24] = extrinsic_curvature[2][0]; q[25] = extrinsic_curvature[2][1]; q[26] = extrinsic_curvature[2][2];

      q[27] = 1.0;

      q[28] = lapse_der[0]; q[29] = lapse_der[1]; q[30] = lapse_der[2];
      q[31] = shift_der[0][0]; q[32] = shift_der[0][1]; q[33] = shift_der[0][2];
      q[34] = shift_der[1][0]; q[35] = shift_der[1][1]; q[36] = shift_der[1][2];
      q[37] = shift_der[2][0]; q[38] = shift_der[2][1]; q[39] = shift_der[2][2];

      q[40] = spatial_metric_der[0][0][0]; q[41] = spatial_metric_der[0][0][1]; q[42] = spatial_metric_der[0][0][2];
      q[43] = spatial_metric_der[0][1][0]; q[44] = spatial_metric_der[0][1][1]; q[45] = spatial_metric_der[0][1][2];
      q[46] = spatial_metric_der[0][2][0]; q[47] = spatial_metric_der[0][2][1]; q[48] = spatial_metric_der[0][2][2];

      q[49] = spatial_metric_der[1][0][0]; q[50] = spatial_metric_der[1][0][1]; q[51] = spatial_metric_der[1][0][2];
      q[52] = spatial_metric_der[1][1][0]; q[53] = spatial_metric_der[1][1][1]; q[54] = spatial_metric_der[1][1][2];
      q[55] = spatial_metric_der[1][2][0]; q[56] = spatial_metric_der[1][2][1]; q[57] = spatial_metric_der[1][2][2];

      q[58] = spatial_metric_der[2][0][0]; q[59] = spatial_metric_der[2][0][1]; q[60] = spatial_metric_der[2][0][2];
      q[61] = spatial_metric_der[2][1][0]; q[62] = spatial_metric_der[2][1][1]; q[63] = spatial_metric_der[2][1][2];
      q[64] = spatial_metric_der[2][2][0]; q[65] = spatial_metric_der[2][2][1]; q[66] = spatial_metric_der[2][2][2];

      q[67] = 0.0;
      q[68] = x; q[69] = y; q[70] = 0.0;

      double prims[71];
      gkyl_gr_euler_prim_vars(gas_gamma, q, prims);
      
      TEST_CHECK( gkyl_compare(prims[0], rho, 1e-13) );
      TEST_CHECK( gkyl_compare(prims[1], u, 1e-13) );
      TEST_CHECK( gkyl_compare(prims[2], v, 1e-13) );
      TEST_CHECK( gkyl_compare(prims[3], w, 1e-13) );
      TEST_CHECK( gkyl_compare(prims[4], p, 1e-13) );

      double fluxes[3][5] = {
        { (lapse * sqrt(spatial_det)) * (rho * W * (vel[0] - (shift[0] / lapse))),
          (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[0] * (vel[0] - (shift[0] / lapse))) + p),
          (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[1] * (vel[0] - (shift[0] / lapse)))),
          (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[2] * (vel[0] - (shift[0] / lapse)))),
          (lapse * sqrt(spatial_det)) * (((rho * h * (W * W)) - p - (rho * W)) * (vel[0] - (shift[0] / lapse)) + (p * vel[0])) },
        { (lapse * sqrt(spatial_det)) * (rho * W * (vel[1] - (shift[1] / lapse))),
          (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[0] * (vel[1] - (shift[1] / lapse)))),
          (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[1] * (vel[1] - (shift[1] / lapse))) + p),
          (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[2] * (vel[1] - (shift[1] / lapse)))),
          (lapse * sqrt(spatial_det)) * (((rho * h * (W * W)) - p - (rho * W)) * (vel[1] - (shift[1] / lapse)) + (p * vel[1])) },
        { (lapse * sqrt(spatial_det)) * (rho * W * (vel[2] - (shift[2] / lapse))),
          (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[0] * (vel[2] - (shift[2] / lapse)))),
          (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[1] * (vel[2] - (shift[2] / lapse)))),
          (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[2] * (vel[2] - (shift[2] / lapse))) + p),
          (lapse * sqrt(spatial_det)) * (((rho * h * (W * W)) - p - (rho * W)) * (vel[2] - (shift[2] / lapse)) + (p * vel[2])) },
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

      double q_local[71], flux_local[71], flux[71];
      for (int d = 0; d < 3; d++) {
        gr_euler->rotate_to_local_func(gr_euler, tau1[d], tau2[d], norm[d], q, q_local);
        gkyl_gr_euler_flux(gas_gamma, q_local, flux_local);
        gr_euler->rotate_to_global_func(gr_euler, tau1[d], tau2[d], norm[d], flux_local, flux);

        for (int i = 0; i < 5; i++) {
          TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-8) );
        }
      }

      double q_l[71], q_g[71];
      for (int d = 0; d < 3; d++) {
        gkyl_wv_eqn_rotate_to_local(gr_euler, tau1[d], tau2[d], norm[d], q, q_l);
        gkyl_wv_eqn_rotate_to_global(gr_euler, tau1[d], tau2[d], norm[d], q_l, q_g);

        for (int i = 0; i < 71; i++) {
          TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-13) );
        }

        double w1[71], q1[71];
        gr_euler->cons_to_riem(gr_euler, q_local, q_local, w1);
        gr_euler->riem_to_cons(gr_euler, q_local, w1, q1);

        for (int i = 0; i < 71; i++) {
          TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-13) );
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric[i]);
        gkyl_free(extrinsic_curvature[i]);
        gkyl_free(shift_der[i]);
    
        for (int j = 0; j < 3; j++) {
          gkyl_free(spatial_metric_der[i][j]);
        }
        gkyl_free(spatial_metric_der[i]);
      }
      gkyl_free(spatial_metric);
      gkyl_free(extrinsic_curvature);
      gkyl_free(shift);
      gkyl_free(vel);
      gkyl_free(lapse_der);
      gkyl_free(shift_der);
      gkyl_free(spatial_metric_der);
    }
  }

  gkyl_wv_eqn_release(gr_euler);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_basic_schwarzschild()
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  struct gkyl_wv_eqn *gr_euler = gkyl_wv_gr_euler_new(gas_gamma, GKYL_STATIC_GAUGE, 0, spacetime, false);

  TEST_CHECK( gr_euler->num_equations == 71 );
  TEST_CHECK( gr_euler->num_waves == 2 );

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double rho = 1.0, u = 0.1, v = 0.2, w = 0.3, p = 1.5;

      double spatial_det, lapse;
      double *shift = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region;

      double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
      }

      double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
      }

      double *lapse_der = gkyl_malloc(sizeof(double[3]));
      double **shift_der = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        shift_der[i] = gkyl_malloc(sizeof(double[3]));
      }

      double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

        for (int j = 0; j < 3; j++) {
          spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
        }
      }

      spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &spatial_det);
      spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse);
      spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift);
      spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);
      
      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

      spacetime->lapse_function_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
      spacetime->shift_vector_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

      double *vel = gkyl_malloc(sizeof(double[3]));
      vel[0] = u; vel[1] = v; vel[2] = w;

      double v_sq = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          v_sq += spatial_metric[i][j] * vel[i] * vel[j];
        }
      }

      double W = 1.0 / sqrt(1.0 - v_sq);
      double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));

      if (!in_excision_region) {
        double q[71];
        q[0] = sqrt(spatial_det) * rho * W;
        q[1] = sqrt(spatial_det) * rho * h * (W * W) * u;
        q[2] = sqrt(spatial_det) * rho * h * (W * W) * v;
        q[3] = sqrt(spatial_det) * rho * h * (W * W) * w;
        q[4] = sqrt(spatial_det) * ((rho * h * (W * W)) - p - (rho * W));

        q[5] = lapse;
        q[6] = shift[0]; q[7] = shift[1]; q[8] = shift[2];

        q[9] = spatial_metric[0][0]; q[10] = spatial_metric[0][1]; q[11] = spatial_metric[0][2];
        q[12] = spatial_metric[1][0]; q[13] = spatial_metric[1][1]; q[14] = spatial_metric[1][2];
        q[15] = spatial_metric[2][0]; q[16] = spatial_metric[2][1]; q[17] = spatial_metric[2][2];

        q[18] = extrinsic_curvature[0][0]; q[19] = extrinsic_curvature[0][1]; q[20] = extrinsic_curvature[0][2];
        q[21] = extrinsic_curvature[1][0]; q[22] = extrinsic_curvature[1][1]; q[23] = extrinsic_curvature[1][2];
        q[24] = extrinsic_curvature[2][0]; q[25] = extrinsic_curvature[2][1]; q[26] = extrinsic_curvature[2][2];

        q[27] = 1.0;

        q[28] = lapse_der[0]; q[29] = lapse_der[1]; q[30] = lapse_der[2];
        q[31] = shift_der[0][0]; q[32] = shift_der[0][1]; q[33] = shift_der[0][2];
        q[34] = shift_der[1][0]; q[35] = shift_der[1][1]; q[36] = shift_der[1][2];
        q[37] = shift_der[2][0]; q[38] = shift_der[2][1]; q[39] = shift_der[2][2];

        q[40] = spatial_metric_der[0][0][0]; q[41] = spatial_metric_der[0][0][1]; q[42] = spatial_metric_der[0][0][2];
        q[43] = spatial_metric_der[0][1][0]; q[44] = spatial_metric_der[0][1][1]; q[45] = spatial_metric_der[0][1][2];
        q[46] = spatial_metric_der[0][2][0]; q[47] = spatial_metric_der[0][2][1]; q[48] = spatial_metric_der[0][2][2];

        q[49] = spatial_metric_der[1][0][0]; q[50] = spatial_metric_der[1][0][1]; q[51] = spatial_metric_der[1][0][2];
        q[52] = spatial_metric_der[1][1][0]; q[53] = spatial_metric_der[1][1][1]; q[54] = spatial_metric_der[1][1][2];
        q[55] = spatial_metric_der[1][2][0]; q[56] = spatial_metric_der[1][2][1]; q[57] = spatial_metric_der[1][2][2];

        q[58] = spatial_metric_der[2][0][0]; q[59] = spatial_metric_der[2][0][1]; q[60] = spatial_metric_der[2][0][2];
        q[61] = spatial_metric_der[2][1][0]; q[62] = spatial_metric_der[2][1][1]; q[63] = spatial_metric_der[2][1][2];
        q[64] = spatial_metric_der[2][2][0]; q[65] = spatial_metric_der[2][2][1]; q[66] = spatial_metric_der[2][2][2];

        q[67] = 0.0;
        q[68] = x; q[69] = y; q[70] = 0.0;

        double prims[71];
        gkyl_gr_euler_prim_vars(gas_gamma, q, prims);
        
        TEST_CHECK( gkyl_compare(prims[0], rho, 1e-12) );
        TEST_CHECK( gkyl_compare(prims[1], u, 1e-12) );
        TEST_CHECK( gkyl_compare(prims[2], v, 1e-12) );
        TEST_CHECK( gkyl_compare(prims[3], w, 1e-12) );
        TEST_CHECK( gkyl_compare(prims[4], p, 1e-12) );

        double fluxes[3][5] = {
          { (lapse * sqrt(spatial_det)) * (rho * W * (vel[0] - (shift[0] / lapse))),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[0] * (vel[0] - (shift[0] / lapse))) + p),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[1] * (vel[0] - (shift[0] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[2] * (vel[0] - (shift[0] / lapse)))),
            (lapse * sqrt(spatial_det)) * (((rho * h * (W * W)) - p - (rho * W)) * (vel[0] - (shift[0] / lapse)) + (p * vel[0])) },
          { (lapse * sqrt(spatial_det)) * (rho * W * (vel[1] - (shift[1] / lapse))),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[0] * (vel[1] - (shift[1] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[1] * (vel[1] - (shift[1] / lapse))) + p),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[2] * (vel[1] - (shift[1] / lapse)))),
            (lapse * sqrt(spatial_det)) * (((rho * h * (W * W)) - p - (rho * W)) * (vel[1] - (shift[1] / lapse)) + (p * vel[1])) },
          { (lapse * sqrt(spatial_det)) * (rho * W * (vel[2] - (shift[2] / lapse))),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[0] * (vel[2] - (shift[2] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[1] * (vel[2] - (shift[2] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[2] * (vel[2] - (shift[2] / lapse))) + p),
            (lapse * sqrt(spatial_det)) * (((rho * h * (W * W)) - p - (rho * W)) * (vel[2] - (shift[2] / lapse)) + (p * vel[2])) },
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

        double q_local[71], flux_local[71], flux[71];
        for (int d = 0; d < 3; d++) {
          gr_euler->rotate_to_local_func(gr_euler, tau1[d], tau2[d], norm[d], q, q_local);
          gkyl_gr_euler_flux(gas_gamma, q_local, flux_local);
          gr_euler->rotate_to_global_func(gr_euler, tau1[d], tau2[d], norm[d], flux_local, flux);

          for (int i = 0; i < 5; i++) {
            TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-1) );
          }
        }

        double q_l[71], q_g[71];
        for (int d = 0; d < 3; d++) {
          gkyl_wv_eqn_rotate_to_local(gr_euler, tau1[d], tau2[d], norm[d], q, q_l);
          gkyl_wv_eqn_rotate_to_global(gr_euler, tau1[d], tau2[d], norm[d], q_l, q_g);

          for (int i = 0; i < 71; i++) {
            TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-13) );
          }

          double w1[71], q1[71];
          gr_euler->cons_to_riem(gr_euler, q_local, q_local, w1);
          gr_euler->riem_to_cons(gr_euler, q_local, w1, q1);

          for (int i = 0; i < 71; i++) {
            TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-13) );
          }
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric[i]);
        gkyl_free(extrinsic_curvature[i]);
        gkyl_free(shift_der[i]);
    
        for (int j = 0; j < 3; j++) {
          gkyl_free(spatial_metric_der[i][j]);
        }
        gkyl_free(spatial_metric_der[i]);
      }
      gkyl_free(spatial_metric);
      gkyl_free(extrinsic_curvature);
      gkyl_free(shift);
      gkyl_free(vel);
      gkyl_free(lapse_der);
      gkyl_free(shift_der);
      gkyl_free(spatial_metric_der);
    }
  }

  gkyl_wv_eqn_release(gr_euler);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_basic_kerr()
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);
  struct gkyl_wv_eqn *gr_euler = gkyl_wv_gr_euler_new(gas_gamma, GKYL_STATIC_GAUGE, 0, spacetime, false);

  TEST_CHECK( gr_euler->num_equations == 71 );
  TEST_CHECK( gr_euler->num_waves == 2 );

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double rho = 1.0, u = 0.1, v = 0.2, w = 0.3, p = 1.5;

      double spatial_det, lapse;
      double *shift = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region;

      double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
      }

      double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
      }

      double *lapse_der = gkyl_malloc(sizeof(double[3]));
      double **shift_der = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        shift_der[i] = gkyl_malloc(sizeof(double[3]));
      }
      
      double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));

        for (int j = 0; j < 3; j++) {
          spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
        }
      }

      spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &spatial_det);
      spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse);
      spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift);
      spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);
      
      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &extrinsic_curvature);

      spacetime->lapse_function_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der);
      spacetime->shift_vector_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der);

      double *vel = gkyl_malloc(sizeof(double[3]));
      vel[0] = u; vel[1] = v; vel[2] = w;

      double v_sq = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          v_sq += spatial_metric[i][j] * vel[i] * vel[j];
        }
      }

      double W = 1.0 / sqrt(1.0 - v_sq);
      double h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)));

      if (!in_excision_region) {
        double q[71];
        q[0] = sqrt(spatial_det) * rho * W;
        q[1] = sqrt(spatial_det) * rho * h * (W * W) * u;
        q[2] = sqrt(spatial_det) * rho * h * (W * W) * v;
        q[3] = sqrt(spatial_det) * rho * h * (W * W) * w;
        q[4] = sqrt(spatial_det) * ((rho * h * (W * W)) - p - (rho * W));

        q[5] = lapse;
        q[6] = shift[0]; q[7] = shift[1]; q[8] = shift[2];

        q[9] = spatial_metric[0][0]; q[10] = spatial_metric[0][1]; q[11] = spatial_metric[0][2];
        q[12] = spatial_metric[1][0]; q[13] = spatial_metric[1][1]; q[14] = spatial_metric[1][2];
        q[15] = spatial_metric[2][0]; q[16] = spatial_metric[2][1]; q[17] = spatial_metric[2][2];

        q[18] = extrinsic_curvature[0][0]; q[19] = extrinsic_curvature[0][1]; q[20] = extrinsic_curvature[0][2];
        q[21] = extrinsic_curvature[1][0]; q[22] = extrinsic_curvature[1][1]; q[23] = extrinsic_curvature[1][2];
        q[24] = extrinsic_curvature[2][0]; q[25] = extrinsic_curvature[2][1]; q[26] = extrinsic_curvature[2][2];

        q[27] = 1.0;

        q[28] = lapse_der[0]; q[29] = lapse_der[1]; q[30] = lapse_der[2];
        q[31] = shift_der[0][0]; q[32] = shift_der[0][1]; q[33] = shift_der[0][2];
        q[34] = shift_der[1][0]; q[35] = shift_der[1][1]; q[36] = shift_der[1][2];
        q[37] = shift_der[2][0]; q[38] = shift_der[2][1]; q[39] = shift_der[2][2];

        q[40] = spatial_metric_der[0][0][0]; q[41] = spatial_metric_der[0][0][1]; q[42] = spatial_metric_der[0][0][2];
        q[43] = spatial_metric_der[0][1][0]; q[44] = spatial_metric_der[0][1][1]; q[45] = spatial_metric_der[0][1][2];
        q[46] = spatial_metric_der[0][2][0]; q[47] = spatial_metric_der[0][2][1]; q[48] = spatial_metric_der[0][2][2];

        q[49] = spatial_metric_der[1][0][0]; q[50] = spatial_metric_der[1][0][1]; q[51] = spatial_metric_der[1][0][2];
        q[52] = spatial_metric_der[1][1][0]; q[53] = spatial_metric_der[1][1][1]; q[54] = spatial_metric_der[1][1][2];
        q[55] = spatial_metric_der[1][2][0]; q[56] = spatial_metric_der[1][2][1]; q[57] = spatial_metric_der[1][2][2];

        q[58] = spatial_metric_der[2][0][0]; q[59] = spatial_metric_der[2][0][1]; q[60] = spatial_metric_der[2][0][2];
        q[61] = spatial_metric_der[2][1][0]; q[62] = spatial_metric_der[2][1][1]; q[63] = spatial_metric_der[2][1][2];
        q[64] = spatial_metric_der[2][2][0]; q[65] = spatial_metric_der[2][2][1]; q[66] = spatial_metric_der[2][2][2];

        q[67] = 0.0;
        q[68] = x; q[69] = y; q[70] = 0.0;

        double prims[71];
        gkyl_gr_euler_prim_vars(gas_gamma, q, prims);
        
        TEST_CHECK( gkyl_compare(prims[0], rho, 1e-12) );
        TEST_CHECK( gkyl_compare(prims[1], u, 1e-12) );
        TEST_CHECK( gkyl_compare(prims[2], v, 1e-12) );
        TEST_CHECK( gkyl_compare(prims[3], w, 1e-12) );
        TEST_CHECK( gkyl_compare(prims[4], p, 1e-12) );

        double fluxes[3][5] = {
          { (lapse * sqrt(spatial_det)) * (rho * W * (vel[0] - (shift[0] / lapse))),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[0] * (vel[0] - (shift[0] / lapse))) + p),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[1] * (vel[0] - (shift[0] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[2] * (vel[0] - (shift[0] / lapse)))),
            (lapse * sqrt(spatial_det)) * (((rho * h * (W * W)) - p - (rho * W)) * (vel[0] - (shift[0] / lapse)) + (p * vel[0])) },
          { (lapse * sqrt(spatial_det)) * (rho * W * (vel[1] - (shift[1] / lapse))),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[0] * (vel[1] - (shift[1] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[1] * (vel[1] - (shift[1] / lapse))) + p),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[2] * (vel[1] - (shift[1] / lapse)))),
            (lapse * sqrt(spatial_det)) * (((rho * h * (W * W)) - p - (rho * W)) * (vel[1] - (shift[1] / lapse)) + (p * vel[1])) },
          { (lapse * sqrt(spatial_det)) * (rho * W * (vel[2] - (shift[2] / lapse))),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[0] * (vel[2] - (shift[2] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[1] * (vel[2] - (shift[2] / lapse)))),
            (lapse * sqrt(spatial_det)) * (rho * h * (W * W) * (vel[2] * (vel[2] - (shift[2] / lapse))) + p),
            (lapse * sqrt(spatial_det)) * (((rho * h * (W * W)) - p - (rho * W)) * (vel[2] - (shift[2] / lapse)) + (p * vel[2])) },
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

        double q_local[71], flux_local[71], flux[71];
        for (int d = 0; d < 3; d++) {
          gr_euler->rotate_to_local_func(gr_euler, tau1[d], tau2[d], norm[d], q, q_local);
          gkyl_gr_euler_flux(gas_gamma, q_local, flux_local);
          gr_euler->rotate_to_global_func(gr_euler, tau1[d], tau2[d], norm[d], flux_local, flux);

          for (int i = 0; i < 5; i++) {
            TEST_CHECK( gkyl_compare(flux[i], fluxes[d][i], 1e-1) );
          }
        }
        
        double q_l[71], q_g[71];
        for (int d = 0; d < 3; d++) {
          gkyl_wv_eqn_rotate_to_local(gr_euler, tau1[d], tau2[d], norm[d], q, q_l);
          gkyl_wv_eqn_rotate_to_global(gr_euler, tau1[d], tau2[d], norm[d], q_l, q_g);

          for (int i = 0; i < 71; i++) {
            TEST_CHECK( gkyl_compare(q[i], q_g[i], 1e-13) );
          }

          double w1[71], q1[71];
          gr_euler->cons_to_riem(gr_euler, q_local, q_local, w1);
          gr_euler->riem_to_cons(gr_euler, q_local, w1, q1);

          for (int i = 0; i < 71; i++) {
            TEST_CHECK( gkyl_compare(q_local[i], q1[i], 1e-13) );
          }
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric[i]);
        gkyl_free(extrinsic_curvature[i]);
        gkyl_free(shift_der[i]);
    
        for (int j = 0; j < 3; j++) {
          gkyl_free(spatial_metric_der[i][j]);
        }
        gkyl_free(spatial_metric_der[i]);
      }
      gkyl_free(spatial_metric);
      gkyl_free(extrinsic_curvature);
      gkyl_free(shift);
      gkyl_free(vel);
      gkyl_free(lapse_der);
      gkyl_free(shift_der);
      gkyl_free(spatial_metric_der);
    }
  }

  gkyl_wv_eqn_release(gr_euler);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_waves_minkowski()
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  struct gkyl_wv_eqn *gr_euler = gkyl_wv_gr_euler_new(gas_gamma, GKYL_STATIC_GAUGE, 0, spacetime, false);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double rho_l = 1.0, u_l = 0.1, v_l = 0.2, w_l = 0.3, p_l = 1.5;
      double rho_r = 0.1, u_r = 0.2, v_r = 0.3, w_r = 0.4, p_r = 0.15;

      double spatial_det_l, spatial_det_r;
      double lapse_l, lapse_r;
      double *shift_l = gkyl_malloc(sizeof(double[3]));
      double *shift_r = gkyl_malloc(sizeof(double[3]));

      double **spatial_metric_l = gkyl_malloc(sizeof(double*[3]));
      double **spatial_metric_r = gkyl_malloc(sizeof(double*[3]));
      double **extrinsic_curvature_l = gkyl_malloc(sizeof(double*[3]));
      double **extrinsic_curvature_r = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_l[i] = gkyl_malloc(sizeof(double[3]));
        spatial_metric_r[i] = gkyl_malloc(sizeof(double[3]));
        extrinsic_curvature_l[i] = gkyl_malloc(sizeof(double[3]));
        extrinsic_curvature_r[i] = gkyl_malloc(sizeof(double[3]));
      }

      double *lapse_der_l = gkyl_malloc(sizeof(double[3]));
      double *lapse_der_r = gkyl_malloc(sizeof(double[3]));
      double **shift_der_l = gkyl_malloc(sizeof(double*[3]));
      double **shift_der_r = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        shift_der_l[i] = gkyl_malloc(sizeof(double[3]));
        shift_der_r[i] = gkyl_malloc(sizeof(double[3]));
      }

      double ***spatial_metric_der_l = gkyl_malloc(sizeof(double**[3]));
      double ***spatial_metric_der_r = gkyl_malloc(sizeof(double**[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_der_l[i] = gkyl_malloc(sizeof(double*[3]));
        spatial_metric_der_r[i] = gkyl_malloc(sizeof(double*[3]));

        for (int j = 0; j < 3; j++) {
          spatial_metric_der_l[i][j] = gkyl_malloc(sizeof(double[3]));
          spatial_metric_der_r[i][j] = gkyl_malloc(sizeof(double[3]));
        }
      }

      spacetime->spatial_metric_det_func(spacetime, 0.0, x - 0.1, y, 0.0, &spatial_det_l);
      spacetime->spatial_metric_det_func(spacetime, 0.0, x + 0.1, y, 0.0, &spatial_det_r);
      spacetime->lapse_function_func(spacetime, 0.0, x - 0.1, y, 0.0, &lapse_l);
      spacetime->lapse_function_func(spacetime, 0.0, x + 0.1, y, 0.0, &lapse_r);
      spacetime->shift_vector_func(spacetime, 0.0, x - 0.1, y, 0.0, &shift_l);
      spacetime->shift_vector_func(spacetime, 0.0, x + 0.1, y, 0.0, &shift_r);

      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x - 0.1, y, 0.0, &spatial_metric_l);
      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x + 0.1, y, 0.0, &spatial_metric_r);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x - 0.1, y, 0.0, 0.1, 0.1, 0.1, &extrinsic_curvature_l);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x + 0.1, y, 0.0, 0.1, 0.1, 0.1, &extrinsic_curvature_r);

      spacetime->lapse_function_der_func(spacetime, 0.0, x - 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der_l);
      spacetime->lapse_function_der_func(spacetime, 0.0, x + 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der_r);
      spacetime->shift_vector_der_func(spacetime, 0.0, x - 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der_l);
      spacetime->shift_vector_der_func(spacetime, 0.0, x + 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der_r);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x - 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der_l);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x + 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der_r);

      double *vel_l = gkyl_malloc(sizeof(double[3]));
      double *vel_r = gkyl_malloc(sizeof(double[3]));
      vel_l[0] = u_l; vel_l[1] = v_l; vel_l[2] = w_l;
      vel_r[0] = u_r; vel_r[1] = v_r; vel_r[2] = w_r;

      double v_sq_l = 0.0, v_sq_r = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          v_sq_l += spatial_metric_l[i][j] * vel_l[i] * vel_l[j];
          v_sq_r += spatial_metric_r[i][j] * vel_r[i] * vel_r[j];
        }
      }

      double W_l = 1.0 / sqrt(1.0 - v_sq_l);
      double W_r = 1.0 / sqrt(1.0 - v_sq_r);
      double h_l = 1.0 + ((p_l / rho_l) * (gas_gamma / (gas_gamma - 1.0)));
      double h_r = 1.0 + ((p_r / rho_r) * (gas_gamma / (gas_gamma - 1.0)));

      double ql[71], qr[71];
      ql[0] = sqrt(spatial_det_l) * rho_l * W_l;
      ql[1] = sqrt(spatial_det_l) * rho_l * h_l * (W_l * W_l) * u_l;
      ql[2] = sqrt(spatial_det_l) * rho_l * h_l * (W_l * W_l) * v_l;
      ql[3] = sqrt(spatial_det_l) * rho_l * h_l * (W_l * W_l) * w_l;
      ql[4] = sqrt(spatial_det_l) * ((rho_l * h_l * (W_l * W_l)) - p_l - (rho_l * W_l));

      ql[5] = lapse_l;
      ql[6] = shift_l[0]; ql[7] = shift_l[1]; ql[8] = shift_l[2];

      ql[9] = spatial_metric_l[0][0]; ql[10] = spatial_metric_l[0][1]; ql[11] = spatial_metric_l[0][2];
      ql[12] = spatial_metric_l[1][0]; ql[13] = spatial_metric_l[1][1]; ql[14] = spatial_metric_l[1][2];
      ql[15] = spatial_metric_l[2][0]; ql[16] = spatial_metric_l[2][1]; ql[17] = spatial_metric_l[2][2];

      ql[18] = extrinsic_curvature_l[0][0]; ql[19] = extrinsic_curvature_l[0][1]; ql[20] = extrinsic_curvature_l[0][2];
      ql[21] = extrinsic_curvature_l[1][0]; ql[22] = extrinsic_curvature_l[1][1]; ql[23] = extrinsic_curvature_l[1][2];
      ql[24] = extrinsic_curvature_l[2][0]; ql[25] = extrinsic_curvature_l[2][1]; ql[26] = extrinsic_curvature_l[2][2];

      ql[27] = 1.0;

      ql[28] = lapse_der_l[0]; ql[29] = lapse_der_l[1]; ql[30] = lapse_der_l[2];
      ql[31] = shift_der_l[0][0]; ql[32] = shift_der_l[0][1]; ql[33] = shift_der_l[0][2];
      ql[34] = shift_der_l[1][0]; ql[35] = shift_der_l[1][1]; ql[36] = shift_der_l[1][2];
      ql[37] = shift_der_l[2][0]; ql[38] = shift_der_l[2][1]; ql[39] = shift_der_l[2][2];

      ql[40] = spatial_metric_der_l[0][0][0]; ql[41] = spatial_metric_der_l[0][0][1]; ql[42] = spatial_metric_der_l[0][0][2];
      ql[43] = spatial_metric_der_l[0][1][0]; ql[44] = spatial_metric_der_l[0][1][1]; ql[45] = spatial_metric_der_l[0][1][2];
      ql[46] = spatial_metric_der_l[0][2][0]; ql[47] = spatial_metric_der_l[0][2][1]; ql[48] = spatial_metric_der_l[0][2][2];

      ql[49] = spatial_metric_der_l[1][0][0]; ql[50] = spatial_metric_der_l[1][0][1]; ql[51] = spatial_metric_der_l[1][0][2];
      ql[52] = spatial_metric_der_l[1][1][0]; ql[53] = spatial_metric_der_l[1][1][1]; ql[54] = spatial_metric_der_l[1][1][2];
      ql[55] = spatial_metric_der_l[1][2][0]; ql[56] = spatial_metric_der_l[1][2][1]; ql[57] = spatial_metric_der_l[1][2][2];

      ql[58] = spatial_metric_der_l[2][0][0]; ql[59] = spatial_metric_der_l[2][0][1]; ql[60] = spatial_metric_der_l[2][0][2];
      ql[61] = spatial_metric_der_l[2][1][0]; ql[62] = spatial_metric_der_l[2][1][1]; ql[63] = spatial_metric_der_l[2][1][2];
      ql[64] = spatial_metric_der_l[2][2][0]; ql[65] = spatial_metric_der_l[2][2][1]; ql[66] = spatial_metric_der_l[2][2][2];

      ql[67] = 0.0;
      ql[68] = x - 0.5; ql[69] = y; ql[70] = 0.0;

      qr[0] = sqrt(spatial_det_r) * rho_r * W_r;
      qr[1] = sqrt(spatial_det_r) * rho_r * h_r * (W_r * W_r) * u_r;
      qr[2] = sqrt(spatial_det_r) * rho_r * h_r * (W_r * W_r) * v_r;
      qr[3] = sqrt(spatial_det_r) * rho_r * h_r * (W_r * W_r) * w_r;
      qr[4] = sqrt(spatial_det_r) * ((rho_r * h_r * (W_r * W_r)) - p_r * (rho_r * W_r));

      qr[5] = lapse_r;
      qr[6] = shift_r[0]; qr[7] = shift_r[1]; qr[8] = shift_r[2];

      qr[9] = spatial_metric_r[0][0]; qr[10] = spatial_metric_r[0][1]; qr[11] = spatial_metric_r[0][2];
      qr[12] = spatial_metric_r[1][0]; qr[13] = spatial_metric_r[1][1]; qr[14] = spatial_metric_r[1][2];
      qr[15] = spatial_metric_r[2][0]; qr[16] = spatial_metric_r[2][1]; qr[17] = spatial_metric_r[2][2];

      qr[18] = extrinsic_curvature_r[0][0]; qr[19] = extrinsic_curvature_r[0][1]; qr[20] = extrinsic_curvature_r[0][2];
      qr[21] = extrinsic_curvature_r[1][0]; qr[22] = extrinsic_curvature_r[1][1]; qr[23] = extrinsic_curvature_r[1][2];
      qr[24] = extrinsic_curvature_r[2][0]; qr[25] = extrinsic_curvature_r[2][1]; qr[26] = extrinsic_curvature_r[2][2];

      qr[27] = 1.0;

      qr[28] = lapse_der_r[0]; qr[29] = lapse_der_r[1]; qr[30] = lapse_der_r[2];
      qr[31] = shift_der_r[0][0]; qr[32] = shift_der_r[0][1]; qr[33] = shift_der_r[0][2];
      qr[34] = shift_der_r[1][0]; qr[35] = shift_der_r[1][1]; qr[36] = shift_der_r[1][2];
      qr[37] = shift_der_r[2][0]; qr[38] = shift_der_r[2][1]; qr[39] = shift_der_r[2][2];

      qr[40] = spatial_metric_der_r[0][0][0]; qr[41] = spatial_metric_der_r[0][0][1]; qr[42] = spatial_metric_der_r[0][0][2];
      qr[43] = spatial_metric_der_r[0][1][0]; qr[44] = spatial_metric_der_r[0][1][1]; qr[45] = spatial_metric_der_r[0][1][2];
      qr[46] = spatial_metric_der_r[0][2][0]; qr[47] = spatial_metric_der_r[0][2][1]; qr[48] = spatial_metric_der_r[0][2][2];

      qr[49] = spatial_metric_der_r[1][0][0]; qr[50] = spatial_metric_der_r[1][0][1]; qr[51] = spatial_metric_der_r[1][0][2];
      qr[52] = spatial_metric_der_r[1][1][0]; qr[53] = spatial_metric_der_r[1][1][1]; qr[54] = spatial_metric_der_r[1][1][2];
      qr[55] = spatial_metric_der_r[1][2][0]; qr[56] = spatial_metric_der_r[1][2][1]; qr[57] = spatial_metric_der_r[1][2][2];

      qr[58] = spatial_metric_der_r[2][0][0]; qr[59] = spatial_metric_der_r[2][0][1]; qr[60] = spatial_metric_der_r[2][0][2];
      qr[61] = spatial_metric_der_r[2][1][0]; qr[62] = spatial_metric_der_r[2][1][1]; qr[63] = spatial_metric_der_r[2][1][2];
      qr[64] = spatial_metric_der_r[2][2][0]; qr[65] = spatial_metric_der_r[2][2][1]; qr[66] = spatial_metric_der_r[2][2][2];

      qr[67] = 0.0;
      qr[68] = x + 0.5; qr[69] = y; qr[70] = 0.0;

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

      for (int d = 0; d < 3; d++) {
        double speeds[3], waves[3 * 71], waves_local[3 * 71];

        double ql_local[71], qr_local[71];
        gkyl_wv_eqn_rotate_to_local(gr_euler, tau1[d], tau2[d], norm[d], ql, ql_local);
        gkyl_wv_eqn_rotate_to_local(gr_euler, tau1[d], tau2[d], norm[d], qr, qr_local);

        double delta[71];
        for (int i = 0; i < 71; i++) {
          delta[i] = qr_local[i] - ql_local[i];
        }

        gkyl_wv_eqn_waves(gr_euler, GKYL_WV_LOW_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

        double apdq_local[71], amdq_local[71];
        gkyl_wv_eqn_qfluct(gr_euler, GKYL_WV_LOW_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

        for (int i = 0; i < 3; i++) {
          gkyl_wv_eqn_rotate_to_global(gr_euler, tau1[d], tau2[d], norm[d], &waves_local[i * 71], &waves[i * 71]);
        }

        double apdq[71], amdq[71];
        gkyl_wv_eqn_rotate_to_global(gr_euler, tau1[d], tau2[d], norm[d], apdq_local, apdq);
        gkyl_wv_eqn_rotate_to_global(gr_euler, tau1[d], tau2[d], norm[d], amdq_local, amdq);

        double fl_local[71], fr_local[71];
        gkyl_gr_euler_flux(gas_gamma, ql_local, fl_local);
        gkyl_gr_euler_flux(gas_gamma, qr_local, fr_local);

        double fl[71], fr[71];
        gkyl_wv_eqn_rotate_to_global(gr_euler, tau1[d], tau2[d], norm[d], fl_local, fl);
        gkyl_wv_eqn_rotate_to_global(gr_euler, tau1[d], tau2[d], norm[d], fr_local, fr);

        for (int i = 0; i < 71; i++) {
          TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-13) );
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric_l[i]);
        gkyl_free(spatial_metric_r[i]);
        gkyl_free(extrinsic_curvature_l[i]);
        gkyl_free(extrinsic_curvature_r[i]);
        gkyl_free(shift_der_l[i]);
        gkyl_free(shift_der_r[i]);
    
        for (int j = 0; j < 3; j++) {
          gkyl_free(spatial_metric_der_l[i][j]);
          gkyl_free(spatial_metric_der_r[i][j]);
        }
        gkyl_free(spatial_metric_der_l[i]);
        gkyl_free(spatial_metric_der_r[i]);
      }
      gkyl_free(spatial_metric_l);
      gkyl_free(spatial_metric_r);
      gkyl_free(extrinsic_curvature_l);
      gkyl_free(extrinsic_curvature_r);
      gkyl_free(shift_l);
      gkyl_free(shift_r);
      gkyl_free(vel_l);
      gkyl_free(vel_r);
      gkyl_free(lapse_der_l);
      gkyl_free(lapse_der_r);
      gkyl_free(shift_der_l);
      gkyl_free(shift_der_r);
      gkyl_free(spatial_metric_der_l);
      gkyl_free(spatial_metric_der_r);
    }
  }

  gkyl_wv_eqn_release(gr_euler);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_waves_schwarzschild()
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  struct gkyl_wv_eqn *gr_euler = gkyl_wv_gr_euler_new(gas_gamma, GKYL_STATIC_GAUGE, 0, spacetime, false);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double rho_l = 1.0, u_l = 0.1, v_l = 0.2, w_l = 0.3, p_l = 1.5;
      double rho_r = 0.1, u_r = 0.2, v_r = 0.3, w_r = 0.4, p_r = 0.15;

      double spatial_det_l, spatial_det_r;
      double lapse_l, lapse_r;
      double *shift_l = gkyl_malloc(sizeof(double[3]));
      double *shift_r = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region_l, in_excision_region_r;

      double **spatial_metric_l = gkyl_malloc(sizeof(double*[3]));
      double **spatial_metric_r = gkyl_malloc(sizeof(double*[3]));
      double **extrinsic_curvature_l = gkyl_malloc(sizeof(double*[3]));
      double **extrinsic_curvature_r = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_l[i] = gkyl_malloc(sizeof(double[3]));
        spatial_metric_r[i] = gkyl_malloc(sizeof(double[3]));
        extrinsic_curvature_l[i] = gkyl_malloc(sizeof(double[3]));
        extrinsic_curvature_r[i] = gkyl_malloc(sizeof(double[3]));
      }

      double *lapse_der_l = gkyl_malloc(sizeof(double[3]));
      double *lapse_der_r = gkyl_malloc(sizeof(double[3]));
      double **shift_der_l = gkyl_malloc(sizeof(double*[3]));
      double **shift_der_r = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        shift_der_l[i] = gkyl_malloc(sizeof(double[3]));
        shift_der_r[i] = gkyl_malloc(sizeof(double[3]));
      }

      double ***spatial_metric_der_l = gkyl_malloc(sizeof(double**[3]));
      double ***spatial_metric_der_r = gkyl_malloc(sizeof(double**[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_der_l[i] = gkyl_malloc(sizeof(double*[3]));
        spatial_metric_der_r[i] = gkyl_malloc(sizeof(double*[3]));

        for (int j = 0; j < 3; j++) {
          spatial_metric_der_l[i][j] = gkyl_malloc(sizeof(double[3]));
          spatial_metric_der_r[i][j] = gkyl_malloc(sizeof(double[3]));
        }
      }

      spacetime->spatial_metric_det_func(spacetime, 0.0, x - 0.1, y, 0.0, &spatial_det_l);
      spacetime->spatial_metric_det_func(spacetime, 0.0, x + 0.1, y, 0.0, &spatial_det_r);
      spacetime->lapse_function_func(spacetime, 0.0, x - 0.1, y, 0.0, &lapse_l);
      spacetime->lapse_function_func(spacetime, 0.0, x + 0.1, y, 0.0, &lapse_r);
      spacetime->shift_vector_func(spacetime, 0.0, x - 0.1, y, 0.0, &shift_l);
      spacetime->shift_vector_func(spacetime, 0.0, x + 0.1, y, 0.0, &shift_r);
      spacetime->excision_region_func(spacetime, 0.0, x - 0.1, y, 0.0, &in_excision_region_l);
      spacetime->excision_region_func(spacetime, 0.0, x + 0.1, y, 0.0, &in_excision_region_r);

      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x - 0.1, y, 0.0, &spatial_metric_l);
      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x + 0.1, y, 0.0, &spatial_metric_r);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x - 0.1, y, 0.0, 0.1, 0.1, 0.1, &extrinsic_curvature_l);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x + 0.1, y, 0.0, 0.1, 0.1, 0.1, &extrinsic_curvature_r);

      spacetime->lapse_function_der_func(spacetime, 0.0, x - 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der_l);
      spacetime->lapse_function_der_func(spacetime, 0.0, x + 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der_r);
      spacetime->shift_vector_der_func(spacetime, 0.0, x - 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der_l);
      spacetime->shift_vector_der_func(spacetime, 0.0, x + 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der_r);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x - 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der_l);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x + 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der_r);

      double *vel_l = gkyl_malloc(sizeof(double[3]));
      double *vel_r = gkyl_malloc(sizeof(double[3]));
      vel_l[0] = u_l; vel_l[1] = v_l; vel_l[2] = w_l;
      vel_r[0] = u_r; vel_r[1] = v_r; vel_r[2] = w_r;

      double v_sq_l = 0.0, v_sq_r = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          v_sq_l += spatial_metric_l[i][j] * vel_l[i] * vel_l[j];
          v_sq_r += spatial_metric_r[i][j] * vel_r[i] * vel_r[j];
        }
      }

      double W_l = 1.0 / sqrt(1.0 - v_sq_l);
      double W_r = 1.0 / sqrt(1.0 - v_sq_r);
      double h_l = 1.0 + ((p_l / rho_l) * (gas_gamma / (gas_gamma - 1.0)));
      double h_r = 1.0 + ((p_r / rho_r) * (gas_gamma / (gas_gamma - 1.0)));

      if (!in_excision_region_l && !in_excision_region_r) {
        double ql[71], qr[71];
        ql[0] = sqrt(spatial_det_l) * rho_l * W_l;
        ql[1] = sqrt(spatial_det_l) * rho_l * h_l * (W_l * W_l) * u_l;
        ql[2] = sqrt(spatial_det_l) * rho_l * h_l * (W_l * W_l) * v_l;
        ql[3] = sqrt(spatial_det_l) * rho_l * h_l * (W_l * W_l) * w_l;
        ql[4] = sqrt(spatial_det_l) * ((rho_l * h_l * (W_l * W_l)) - p_l - (rho_l * W_l));

        ql[5] = lapse_l;
        ql[6] = shift_l[0]; ql[7] = shift_l[1]; ql[8] = shift_l[2];

        ql[9] = spatial_metric_l[0][0]; ql[10] = spatial_metric_l[0][1]; ql[11] = spatial_metric_l[0][2];
        ql[12] = spatial_metric_l[1][0]; ql[13] = spatial_metric_l[1][1]; ql[14] = spatial_metric_l[1][2];
        ql[15] = spatial_metric_l[2][0]; ql[16] = spatial_metric_l[2][1]; ql[17] = spatial_metric_l[2][2];

        ql[18] = extrinsic_curvature_l[0][0]; ql[19] = extrinsic_curvature_l[0][1]; ql[20] = extrinsic_curvature_l[0][2];
        ql[21] = extrinsic_curvature_l[1][0]; ql[22] = extrinsic_curvature_l[1][1]; ql[23] = extrinsic_curvature_l[1][2];
        ql[24] = extrinsic_curvature_l[2][0]; ql[25] = extrinsic_curvature_l[2][1]; ql[26] = extrinsic_curvature_l[2][2];

        ql[27] = 1.0;

        ql[28] = lapse_der_l[0]; ql[29] = lapse_der_l[1]; ql[30] = lapse_der_l[2];
        ql[31] = shift_der_l[0][0]; ql[32] = shift_der_l[0][1]; ql[33] = shift_der_l[0][2];
        ql[34] = shift_der_l[1][0]; ql[35] = shift_der_l[1][1]; ql[36] = shift_der_l[1][2];
        ql[37] = shift_der_l[2][0]; ql[38] = shift_der_l[2][1]; ql[39] = shift_der_l[2][2];

        ql[40] = spatial_metric_der_l[0][0][0]; ql[41] = spatial_metric_der_l[0][0][1]; ql[42] = spatial_metric_der_l[0][0][2];
        ql[43] = spatial_metric_der_l[0][1][0]; ql[44] = spatial_metric_der_l[0][1][1]; ql[45] = spatial_metric_der_l[0][1][2];
        ql[46] = spatial_metric_der_l[0][2][0]; ql[47] = spatial_metric_der_l[0][2][1]; ql[48] = spatial_metric_der_l[0][2][2];

        ql[49] = spatial_metric_der_l[1][0][0]; ql[50] = spatial_metric_der_l[1][0][1]; ql[51] = spatial_metric_der_l[1][0][2];
        ql[52] = spatial_metric_der_l[1][1][0]; ql[53] = spatial_metric_der_l[1][1][1]; ql[54] = spatial_metric_der_l[1][1][2];
        ql[55] = spatial_metric_der_l[1][2][0]; ql[56] = spatial_metric_der_l[1][2][1]; ql[57] = spatial_metric_der_l[1][2][2];

        ql[58] = spatial_metric_der_l[2][0][0]; ql[59] = spatial_metric_der_l[2][0][1]; ql[60] = spatial_metric_der_l[2][0][2];
        ql[61] = spatial_metric_der_l[2][1][0]; ql[62] = spatial_metric_der_l[2][1][1]; ql[63] = spatial_metric_der_l[2][1][2];
        ql[64] = spatial_metric_der_l[2][2][0]; ql[65] = spatial_metric_der_l[2][2][1]; ql[66] = spatial_metric_der_l[2][2][2];

        ql[67] = 0.0;
        ql[68] = x - 0.5; ql[69] = y; ql[70] = 0.0;

        qr[0] = sqrt(spatial_det_r) * rho_r * W_r;
        qr[1] = sqrt(spatial_det_r) * rho_r * h_r * (W_r * W_r) * u_r;
        qr[2] = sqrt(spatial_det_r) * rho_r * h_r * (W_r * W_r) * v_r;
        qr[3] = sqrt(spatial_det_r) * rho_r * h_r * (W_r * W_r) * w_r;
        qr[4] = sqrt(spatial_det_r) * ((rho_r * h_r * (W_r * W_r)) - p_r * (rho_r * W_r));

        qr[5] = lapse_r;
        qr[6] = shift_r[0]; qr[7] = shift_r[1]; qr[8] = shift_r[2];

        qr[9] = spatial_metric_r[0][0]; qr[10] = spatial_metric_r[0][1]; qr[11] = spatial_metric_r[0][2];
        qr[12] = spatial_metric_r[1][0]; qr[13] = spatial_metric_r[1][1]; qr[14] = spatial_metric_r[1][2];
        qr[15] = spatial_metric_r[2][0]; qr[16] = spatial_metric_r[2][1]; qr[17] = spatial_metric_r[2][2];

        qr[18] = extrinsic_curvature_r[0][0]; qr[19] = extrinsic_curvature_r[0][1]; qr[20] = extrinsic_curvature_r[0][2];
        qr[21] = extrinsic_curvature_r[1][0]; qr[22] = extrinsic_curvature_r[1][1]; qr[23] = extrinsic_curvature_r[1][2];
        qr[24] = extrinsic_curvature_r[2][0]; qr[25] = extrinsic_curvature_r[2][1]; qr[26] = extrinsic_curvature_r[2][2];

        qr[27] = 1.0;

        qr[28] = lapse_der_r[0]; qr[29] = lapse_der_r[1]; qr[30] = lapse_der_r[2];
        qr[31] = shift_der_r[0][0]; qr[32] = shift_der_r[0][1]; qr[33] = shift_der_r[0][2];
        qr[34] = shift_der_r[1][0]; qr[35] = shift_der_r[1][1]; qr[36] = shift_der_r[1][2];
        qr[37] = shift_der_r[2][0]; qr[38] = shift_der_r[2][1]; qr[39] = shift_der_r[2][2];

        qr[40] = spatial_metric_der_r[0][0][0]; qr[41] = spatial_metric_der_r[0][0][1]; qr[42] = spatial_metric_der_r[0][0][2];
        qr[43] = spatial_metric_der_r[0][1][0]; qr[44] = spatial_metric_der_r[0][1][1]; qr[45] = spatial_metric_der_r[0][1][2];
        qr[46] = spatial_metric_der_r[0][2][0]; qr[47] = spatial_metric_der_r[0][2][1]; qr[48] = spatial_metric_der_r[0][2][2];

        qr[49] = spatial_metric_der_r[1][0][0]; qr[50] = spatial_metric_der_r[1][0][1]; qr[51] = spatial_metric_der_r[1][0][2];
        qr[52] = spatial_metric_der_r[1][1][0]; qr[53] = spatial_metric_der_r[1][1][1]; qr[54] = spatial_metric_der_r[1][1][2];
        qr[55] = spatial_metric_der_r[1][2][0]; qr[56] = spatial_metric_der_r[1][2][1]; qr[57] = spatial_metric_der_r[1][2][2];

        qr[58] = spatial_metric_der_r[2][0][0]; qr[59] = spatial_metric_der_r[2][0][1]; qr[60] = spatial_metric_der_r[2][0][2];
        qr[61] = spatial_metric_der_r[2][1][0]; qr[62] = spatial_metric_der_r[2][1][1]; qr[63] = spatial_metric_der_r[2][1][2];
        qr[64] = spatial_metric_der_r[2][2][0]; qr[65] = spatial_metric_der_r[2][2][1]; qr[66] = spatial_metric_der_r[2][2][2];

        qr[67] = 0.0;
        qr[68] = x + 0.5; qr[69] = y; qr[70] = 0.0;

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

        for (int d = 0; d < 3; d++) {
          double speeds[3], waves[3 * 71], waves_local[3 * 71];

          double ql_local[71], qr_local[71];
          gkyl_wv_eqn_rotate_to_local(gr_euler, tau1[d], tau2[d], norm[d], ql, ql_local);
          gkyl_wv_eqn_rotate_to_local(gr_euler, tau1[d], tau2[d], norm[d], qr, qr_local);

          double delta[71];
          for (int i = 0; i < 71; i++) {
            delta[i] = qr_local[i] - ql_local[i];
          }

          gkyl_wv_eqn_waves(gr_euler, GKYL_WV_LOW_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

          double apdq_local[71], amdq_local[71];
          gkyl_wv_eqn_qfluct(gr_euler, GKYL_WV_LOW_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

          for (int i = 0; i < 3; i++) {
            gkyl_wv_eqn_rotate_to_global(gr_euler, tau1[d], tau2[d], norm[d], &waves_local[i * 71], &waves[i * 71]);
          }

          double apdq[71], amdq[71];
          gkyl_wv_eqn_rotate_to_global(gr_euler, tau1[d], tau2[d], norm[d], apdq_local, apdq);
          gkyl_wv_eqn_rotate_to_global(gr_euler, tau1[d], tau2[d], norm[d], amdq_local, amdq);

          double fl_local[71], fr_local[71];
          gkyl_gr_euler_flux(gas_gamma, ql_local, fl_local);
          gkyl_gr_euler_flux(gas_gamma, qr_local, fr_local);

          double fl[71], fr[71];
          gkyl_wv_eqn_rotate_to_global(gr_euler, tau1[d], tau2[d], norm[d], fl_local, fl);
          gkyl_wv_eqn_rotate_to_global(gr_euler, tau1[d], tau2[d], norm[d], fr_local, fr);

          for (int i = 0; i < 71; i++) {
            TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-12) );
          }
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric_l[i]);
        gkyl_free(spatial_metric_r[i]);
        gkyl_free(extrinsic_curvature_l[i]);
        gkyl_free(extrinsic_curvature_r[i]);
        gkyl_free(shift_der_l[i]);
        gkyl_free(shift_der_r[i]);
    
        for (int j = 0; j < 3; j++) {
          gkyl_free(spatial_metric_der_l[i][j]);
          gkyl_free(spatial_metric_der_r[i][j]);
        }
        gkyl_free(spatial_metric_der_l[i]);
        gkyl_free(spatial_metric_der_r[i]);
      }
      gkyl_free(spatial_metric_l);
      gkyl_free(spatial_metric_r);
      gkyl_free(extrinsic_curvature_l);
      gkyl_free(extrinsic_curvature_r);
      gkyl_free(shift_l);
      gkyl_free(shift_r);
      gkyl_free(vel_l);
      gkyl_free(vel_r);
      gkyl_free(lapse_der_l);
      gkyl_free(lapse_der_r);
      gkyl_free(shift_der_l);
      gkyl_free(shift_der_r);
      gkyl_free(spatial_metric_der_l);
      gkyl_free(spatial_metric_der_r);
    }
  }

  gkyl_wv_eqn_release(gr_euler);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_euler_waves_kerr()
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);
  struct gkyl_wv_eqn *gr_euler = gkyl_wv_gr_euler_new(gas_gamma, GKYL_STATIC_GAUGE, 0, spacetime, false);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double rho_l = 1.0, u_l = 0.1, v_l = 0.2, w_l = 0.3, p_l = 1.5;
      double rho_r = 0.1, u_r = 0.2, v_r = 0.3, w_r = 0.4, p_r = 0.15;

      double spatial_det_l, spatial_det_r;
      double lapse_l, lapse_r;
      double *shift_l = gkyl_malloc(sizeof(double[3]));
      double *shift_r = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region_l, in_excision_region_r;

      double **spatial_metric_l = gkyl_malloc(sizeof(double*[3]));
      double **spatial_metric_r = gkyl_malloc(sizeof(double*[3]));
      double **extrinsic_curvature_l = gkyl_malloc(sizeof(double*[3]));
      double **extrinsic_curvature_r = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_l[i] = gkyl_malloc(sizeof(double[3]));
        spatial_metric_r[i] = gkyl_malloc(sizeof(double[3]));
        extrinsic_curvature_l[i] = gkyl_malloc(sizeof(double[3]));
        extrinsic_curvature_r[i] = gkyl_malloc(sizeof(double[3]));
      }

      double *lapse_der_l = gkyl_malloc(sizeof(double[3]));
      double *lapse_der_r = gkyl_malloc(sizeof(double[3]));
      double **shift_der_l = gkyl_malloc(sizeof(double*[3]));
      double **shift_der_r = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) {
        shift_der_l[i] = gkyl_malloc(sizeof(double[3]));
        shift_der_r[i] = gkyl_malloc(sizeof(double[3]));
      }

      double ***spatial_metric_der_l = gkyl_malloc(sizeof(double**[3]));
      double ***spatial_metric_der_r = gkyl_malloc(sizeof(double**[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_der_l[i] = gkyl_malloc(sizeof(double*[3]));
        spatial_metric_der_r[i] = gkyl_malloc(sizeof(double*[3]));

        for (int j = 0; j < 3; j++) {
          spatial_metric_der_l[i][j] = gkyl_malloc(sizeof(double[3]));
          spatial_metric_der_r[i][j] = gkyl_malloc(sizeof(double[3]));
        }
      }

      spacetime->spatial_metric_det_func(spacetime, 0.0, x - 0.1, y, 0.0, &spatial_det_l);
      spacetime->spatial_metric_det_func(spacetime, 0.0, x + 0.1, y, 0.0, &spatial_det_r);
      spacetime->lapse_function_func(spacetime, 0.0, x - 0.1, y, 0.0, &lapse_l);
      spacetime->lapse_function_func(spacetime, 0.0, x + 0.1, y, 0.0, &lapse_r);
      spacetime->shift_vector_func(spacetime, 0.0, x - 0.1, y, 0.0, &shift_l);
      spacetime->shift_vector_func(spacetime, 0.0, x + 0.1, y, 0.0, &shift_r);
      spacetime->excision_region_func(spacetime, 0.0, x - 0.1, y, 0.0, &in_excision_region_l);
      spacetime->excision_region_func(spacetime, 0.0, x + 0.1, y, 0.0, &in_excision_region_r);

      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x - 0.1, y, 0.0, &spatial_metric_l);
      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x + 0.1, y, 0.0, &spatial_metric_r);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x - 0.1, y, 0.0, 0.1, 0.1, 0.1, &extrinsic_curvature_l);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x + 0.1, y, 0.0, 0.1, 0.1, 0.1, &extrinsic_curvature_r);

      spacetime->lapse_function_der_func(spacetime, 0.0, x - 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der_l);
      spacetime->lapse_function_der_func(spacetime, 0.0, x + 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &lapse_der_r);
      spacetime->shift_vector_der_func(spacetime, 0.0, x - 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der_l);
      spacetime->shift_vector_der_func(spacetime, 0.0, x + 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &shift_der_r);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x - 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der_l);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x + 0.1, y, 0.0, pow(10.0, -8.0), pow(10.0, -8.0), pow(10.0, -8.0), &spatial_metric_der_r);

      double *vel_l = gkyl_malloc(sizeof(double[3]));
      double *vel_r = gkyl_malloc(sizeof(double[3]));
      vel_l[0] = u_l; vel_l[1] = v_l; vel_l[2] = w_l;
      vel_r[0] = u_r; vel_r[1] = v_r; vel_r[2] = w_r;

      double v_sq_l = 0.0, v_sq_r = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          v_sq_l += spatial_metric_l[i][j] * vel_l[i] * vel_l[j];
          v_sq_r += spatial_metric_r[i][j] * vel_r[i] * vel_r[j];
        }
      }

      double W_l = 1.0 / sqrt(1.0 - v_sq_l);
      double W_r = 1.0 / sqrt(1.0 - v_sq_r);
      double h_l = 1.0 + ((p_l / rho_l) * (gas_gamma / (gas_gamma - 1.0)));
      double h_r = 1.0 + ((p_r / rho_r) * (gas_gamma / (gas_gamma - 1.0)));

      if (!in_excision_region_l && !in_excision_region_r) {
        double ql[71], qr[71];
        ql[0] = sqrt(spatial_det_l) * rho_l * W_l;
        ql[1] = sqrt(spatial_det_l) * rho_l * h_l * (W_l * W_l) * u_l;
        ql[2] = sqrt(spatial_det_l) * rho_l * h_l * (W_l * W_l) * v_l;
        ql[3] = sqrt(spatial_det_l) * rho_l * h_l * (W_l * W_l) * w_l;
        ql[4] = sqrt(spatial_det_l) * ((rho_l * h_l * (W_l * W_l)) - p_l - (rho_l * W_l));

        ql[5] = lapse_l;
        ql[6] = shift_l[0]; ql[7] = shift_l[1]; ql[8] = shift_l[2];

        ql[9] = spatial_metric_l[0][0]; ql[10] = spatial_metric_l[0][1]; ql[11] = spatial_metric_l[0][2];
        ql[12] = spatial_metric_l[1][0]; ql[13] = spatial_metric_l[1][1]; ql[14] = spatial_metric_l[1][2];
        ql[15] = spatial_metric_l[2][0]; ql[16] = spatial_metric_l[2][1]; ql[17] = spatial_metric_l[2][2];

        ql[18] = extrinsic_curvature_l[0][0]; ql[19] = extrinsic_curvature_l[0][1]; ql[20] = extrinsic_curvature_l[0][2];
        ql[21] = extrinsic_curvature_l[1][0]; ql[22] = extrinsic_curvature_l[1][1]; ql[23] = extrinsic_curvature_l[1][2];
        ql[24] = extrinsic_curvature_l[2][0]; ql[25] = extrinsic_curvature_l[2][1]; ql[26] = extrinsic_curvature_l[2][2];

        ql[27] = 1.0;

        ql[28] = lapse_der_l[0]; ql[29] = lapse_der_l[1]; ql[30] = lapse_der_l[2];
        ql[31] = shift_der_l[0][0]; ql[32] = shift_der_l[0][1]; ql[33] = shift_der_l[0][2];
        ql[34] = shift_der_l[1][0]; ql[35] = shift_der_l[1][1]; ql[36] = shift_der_l[1][2];
        ql[37] = shift_der_l[2][0]; ql[38] = shift_der_l[2][1]; ql[39] = shift_der_l[2][2];

        ql[40] = spatial_metric_der_l[0][0][0]; ql[41] = spatial_metric_der_l[0][0][1]; ql[42] = spatial_metric_der_l[0][0][2];
        ql[43] = spatial_metric_der_l[0][1][0]; ql[44] = spatial_metric_der_l[0][1][1]; ql[45] = spatial_metric_der_l[0][1][2];
        ql[46] = spatial_metric_der_l[0][2][0]; ql[47] = spatial_metric_der_l[0][2][1]; ql[48] = spatial_metric_der_l[0][2][2];

        ql[49] = spatial_metric_der_l[1][0][0]; ql[50] = spatial_metric_der_l[1][0][1]; ql[51] = spatial_metric_der_l[1][0][2];
        ql[52] = spatial_metric_der_l[1][1][0]; ql[53] = spatial_metric_der_l[1][1][1]; ql[54] = spatial_metric_der_l[1][1][2];
        ql[55] = spatial_metric_der_l[1][2][0]; ql[56] = spatial_metric_der_l[1][2][1]; ql[57] = spatial_metric_der_l[1][2][2];

        ql[58] = spatial_metric_der_l[2][0][0]; ql[59] = spatial_metric_der_l[2][0][1]; ql[60] = spatial_metric_der_l[2][0][2];
        ql[61] = spatial_metric_der_l[2][1][0]; ql[62] = spatial_metric_der_l[2][1][1]; ql[63] = spatial_metric_der_l[2][1][2];
        ql[64] = spatial_metric_der_l[2][2][0]; ql[65] = spatial_metric_der_l[2][2][1]; ql[66] = spatial_metric_der_l[2][2][2];

        ql[67] = 0.0;
        ql[68] = x - 0.5; ql[69] = y; ql[70] = 0.0;

        qr[0] = sqrt(spatial_det_r) * rho_r * W_r;
        qr[1] = sqrt(spatial_det_r) * rho_r * h_r * (W_r * W_r) * u_r;
        qr[2] = sqrt(spatial_det_r) * rho_r * h_r * (W_r * W_r) * v_r;
        qr[3] = sqrt(spatial_det_r) * rho_r * h_r * (W_r * W_r) * w_r;
        qr[4] = sqrt(spatial_det_r) * ((rho_r * h_r * (W_r * W_r)) - p_r * (rho_r * W_r));

        qr[5] = lapse_r;
        qr[6] = shift_r[0]; qr[7] = shift_r[1]; qr[8] = shift_r[2];

        qr[9] = spatial_metric_r[0][0]; qr[10] = spatial_metric_r[0][1]; qr[11] = spatial_metric_r[0][2];
        qr[12] = spatial_metric_r[1][0]; qr[13] = spatial_metric_r[1][1]; qr[14] = spatial_metric_r[1][2];
        qr[15] = spatial_metric_r[2][0]; qr[16] = spatial_metric_r[2][1]; qr[17] = spatial_metric_r[2][2];

        qr[18] = extrinsic_curvature_r[0][0]; qr[19] = extrinsic_curvature_r[0][1]; qr[20] = extrinsic_curvature_r[0][2];
        qr[21] = extrinsic_curvature_r[1][0]; qr[22] = extrinsic_curvature_r[1][1]; qr[23] = extrinsic_curvature_r[1][2];
        qr[24] = extrinsic_curvature_r[2][0]; qr[25] = extrinsic_curvature_r[2][1]; qr[26] = extrinsic_curvature_r[2][2];

        qr[27] = 1.0;

        qr[28] = lapse_der_r[0]; qr[29] = lapse_der_r[1]; qr[30] = lapse_der_r[2];
        qr[31] = shift_der_r[0][0]; qr[32] = shift_der_r[0][1]; qr[33] = shift_der_r[0][2];
        qr[34] = shift_der_r[1][0]; qr[35] = shift_der_r[1][1]; qr[36] = shift_der_r[1][2];
        qr[37] = shift_der_r[2][0]; qr[38] = shift_der_r[2][1]; qr[39] = shift_der_r[2][2];

        qr[40] = spatial_metric_der_r[0][0][0]; qr[41] = spatial_metric_der_r[0][0][1]; qr[42] = spatial_metric_der_r[0][0][2];
        qr[43] = spatial_metric_der_r[0][1][0]; qr[44] = spatial_metric_der_r[0][1][1]; qr[45] = spatial_metric_der_r[0][1][2];
        qr[46] = spatial_metric_der_r[0][2][0]; qr[47] = spatial_metric_der_r[0][2][1]; qr[48] = spatial_metric_der_r[0][2][2];

        qr[49] = spatial_metric_der_r[1][0][0]; qr[50] = spatial_metric_der_r[1][0][1]; qr[51] = spatial_metric_der_r[1][0][2];
        qr[52] = spatial_metric_der_r[1][1][0]; qr[53] = spatial_metric_der_r[1][1][1]; qr[54] = spatial_metric_der_r[1][1][2];
        qr[55] = spatial_metric_der_r[1][2][0]; qr[56] = spatial_metric_der_r[1][2][1]; qr[57] = spatial_metric_der_r[1][2][2];

        qr[58] = spatial_metric_der_r[2][0][0]; qr[59] = spatial_metric_der_r[2][0][1]; qr[60] = spatial_metric_der_r[2][0][2];
        qr[61] = spatial_metric_der_r[2][1][0]; qr[62] = spatial_metric_der_r[2][1][1]; qr[63] = spatial_metric_der_r[2][1][2];
        qr[64] = spatial_metric_der_r[2][2][0]; qr[65] = spatial_metric_der_r[2][2][1]; qr[66] = spatial_metric_der_r[2][2][2];

        qr[67] = 0.0;
        qr[68] = x + 0.5; qr[69] = y; qr[70] = 0.0;

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

        for (int d = 0; d < 3; d++) {
          double speeds[3], waves[3 * 71], waves_local[3 * 71];

          double ql_local[71], qr_local[71];
          gkyl_wv_eqn_rotate_to_local(gr_euler, tau1[d], tau2[d], norm[d], ql, ql_local);
          gkyl_wv_eqn_rotate_to_local(gr_euler, tau1[d], tau2[d], norm[d], qr, qr_local);

          double delta[71];
          for (int i = 0; i < 71; i++) {
            delta[i] = qr_local[i] - ql_local[i];
          }

          gkyl_wv_eqn_waves(gr_euler, GKYL_WV_LOW_ORDER_FLUX, delta, ql_local, qr_local, 1.0, 1.0, waves_local, speeds);

          double apdq_local[71], amdq_local[71];
          gkyl_wv_eqn_qfluct(gr_euler, GKYL_WV_LOW_ORDER_FLUX, ql_local, qr_local, 1.0, 1.0, waves_local, speeds, amdq_local, apdq_local);

          for (int i = 0; i < 3; i++) {
            gkyl_wv_eqn_rotate_to_global(gr_euler, tau1[d], tau2[d], norm[d], &waves_local[i * 71], &waves[i * 71]);
          }

          double apdq[71], amdq[71];
          gkyl_wv_eqn_rotate_to_global(gr_euler, tau1[d], tau2[d], norm[d], apdq_local, apdq);
          gkyl_wv_eqn_rotate_to_global(gr_euler, tau1[d], tau2[d], norm[d], amdq_local, amdq);

          double fl_local[71], fr_local[71];
          gkyl_gr_euler_flux(gas_gamma, ql_local, fl_local);
          gkyl_gr_euler_flux(gas_gamma, qr_local, fr_local);

          double fl[71], fr[71];
          gkyl_wv_eqn_rotate_to_global(gr_euler, tau1[d], tau2[d], norm[d], fl_local, fl);
          gkyl_wv_eqn_rotate_to_global(gr_euler, tau1[d], tau2[d], norm[d], fr_local, fr);

          for (int i = 0; i < 71; i++) {
            TEST_CHECK( gkyl_compare(fr[i] - fl[i], amdq[i] + apdq[i], 1e-12) );
          }
        }
      }

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric_l[i]);
        gkyl_free(spatial_metric_r[i]);
        gkyl_free(extrinsic_curvature_l[i]);
        gkyl_free(extrinsic_curvature_r[i]);
        gkyl_free(shift_der_l[i]);
        gkyl_free(shift_der_r[i]);
    
        for (int j = 0; j < 3; j++) {
          gkyl_free(spatial_metric_der_l[i][j]);
          gkyl_free(spatial_metric_der_r[i][j]);
        }
        gkyl_free(spatial_metric_der_l[i]);
        gkyl_free(spatial_metric_der_r[i]);
      }
      gkyl_free(spatial_metric_l);
      gkyl_free(spatial_metric_r);
      gkyl_free(extrinsic_curvature_l);
      gkyl_free(extrinsic_curvature_r);
      gkyl_free(shift_l);
      gkyl_free(shift_r);
      gkyl_free(vel_l);
      gkyl_free(vel_r);
      gkyl_free(lapse_der_l);
      gkyl_free(lapse_der_r);
      gkyl_free(shift_der_l);
      gkyl_free(shift_der_r);
      gkyl_free(spatial_metric_der_l);
      gkyl_free(spatial_metric_der_r);
    }
  }

  gkyl_wv_eqn_release(gr_euler);
  gkyl_gr_spacetime_release(spacetime);
}

// Stringent prim_vars roundtrip stress-test. Sweeps a parameter grid of
// primitive states (varying ρ, p over 5+ orders of magnitude, |v| from low
// to relativistic) at a fixed set of sample positions including explicit
// near-horizon points where γ_ij deviates strongly from δ_ij and β^i is
// large from frame-dragging. For each (state, position) tuple, builds the
// conserved q via Convention-B (q[i+1] = √γ·ρhW²·v^i), calls prim_vars,
// and records the worst-case RELATIVE error per component. Test data lives
// in prim_vars_stringent_data.h so all four GR Euler ctest files exercise
// the same regime.
static void
run_prim_vars_stringent(struct gkyl_gr_spacetime *spacetime, const char *label)
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_wv_eqn *gr_euler = gkyl_wv_gr_euler_new(
    gas_gamma, GKYL_STATIC_GAUGE, 0, spacetime, false);

  double max_rel_rho = 0.0, max_rel_u = 0.0, max_rel_v = 0.0;
  double max_rel_w   = 0.0, max_rel_p = 0.0;
  int    worst_s = -1, worst_p = -1;
  int    n_samples = 0, n_skipped_excise = 0, n_skipped_super = 0;

  const double rel_floor = STRINGENT_REL_FLOOR;

  for (int sx = 0; sx < N_STRINGENT_STATES; sx++) {
    for (int px = 0; px < N_STRINGENT_POSITIONS; px++) {
      double rho_in = stringent_states[sx].rho;
      double p_in   = stringent_states[sx].p;
      double u_in   = stringent_states[sx].vx;
      double v_in   = stringent_states[sx].vy;
      double w_in   = stringent_states[sx].vz;
      double x = stringent_positions[px][0];
      double y = stringent_positions[px][1];
      double z = stringent_positions[px][2];

      double spatial_det, lapse;
      double *shift = gkyl_malloc(sizeof(double[3]));
      bool in_excision_region;
      double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
      double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
      double *lapse_der = gkyl_malloc(sizeof(double[3]));
      double **shift_der = gkyl_malloc(sizeof(double*[3]));
      for (int i = 0; i < 3; i++) shift_der[i] = gkyl_malloc(sizeof(double[3]));
      double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
      for (int i = 0; i < 3; i++) {
        spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));
        for (int j = 0; j < 3; j++) spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
      }

      spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, z, &spatial_det);
      spacetime->lapse_function_func(spacetime, 0.0, x, y, z, &lapse);
      spacetime->shift_vector_func(spacetime, 0.0, x, y, z, &shift);
      spacetime->excision_region_func(spacetime, 0.0, x, y, z, &in_excision_region);
      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, z, &spatial_metric);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, z,
        pow(10.0,-8.0), pow(10.0,-8.0), pow(10.0,-8.0), &extrinsic_curvature);
      spacetime->lapse_function_der_func(spacetime, 0.0, x, y, z,
        pow(10.0,-8.0), pow(10.0,-8.0), pow(10.0,-8.0), &lapse_der);
      spacetime->shift_vector_der_func(spacetime, 0.0, x, y, z,
        pow(10.0,-8.0), pow(10.0,-8.0), pow(10.0,-8.0), &shift_der);
      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, z,
        pow(10.0,-8.0), pow(10.0,-8.0), pow(10.0,-8.0), &spatial_metric_der);

      bool skip = false;
      if (in_excision_region) { n_skipped_excise++; skip = true; }
      double v_sq = 0.0;
      if (!skip) {
        double vel[3] = { u_in, v_in, w_in };
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++)
            v_sq += spatial_metric[i][j] * vel[i] * vel[j];
        if (v_sq >= 1.0 - 1.0e-6) { n_skipped_super++; skip = true; }
      }

      if (!skip) {
        double W = 1.0 / sqrt(1.0 - v_sq);
        double h = 1.0 + ((p_in / rho_in) * (gas_gamma / (gas_gamma - 1.0)));

        double q[71];
        q[0] = sqrt(spatial_det) * rho_in * W;
        q[1] = sqrt(spatial_det) * rho_in * h * (W*W) * u_in;
        q[2] = sqrt(spatial_det) * rho_in * h * (W*W) * v_in;
        q[3] = sqrt(spatial_det) * rho_in * h * (W*W) * w_in;
        q[4] = sqrt(spatial_det) * ((rho_in * h * (W*W)) - p_in - (rho_in * W));
        q[5] = lapse;
        q[6] = shift[0]; q[7] = shift[1]; q[8] = shift[2];
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++) q[9 + 3*i + j] = spatial_metric[i][j];
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++) q[18 + 3*i + j] = extrinsic_curvature[i][j];
        q[27] = 1.0;
        q[28] = lapse_der[0]; q[29] = lapse_der[1]; q[30] = lapse_der[2];
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++) q[31 + 3*i + j] = shift_der[i][j];
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++) q[40 + 9*i + 3*j + k] = spatial_metric_der[i][j][k];
        q[67] = 0.0; q[68] = x; q[69] = y; q[70] = z;

        double prims[71];
        gkyl_gr_euler_prim_vars(gas_gamma, q, prims);

        double rel_rho = fabs(prims[0] - rho_in) / fmax(fabs(rho_in), rel_floor);
        double rel_u   = fabs(prims[1] - u_in)   / fmax(fabs(u_in),   rel_floor);
        double rel_v   = fabs(prims[2] - v_in)   / fmax(fabs(v_in),   rel_floor);
        double rel_w   = fabs(prims[3] - w_in)   / fmax(fabs(w_in),   rel_floor);
        double rel_p   = fabs(prims[4] - p_in)   / fmax(fabs(p_in),   rel_floor);

        if (rel_rho > max_rel_rho) { max_rel_rho = rel_rho; worst_s = sx; worst_p = px; }
        if (rel_u   > max_rel_u)   { max_rel_u   = rel_u;   }
        if (rel_v   > max_rel_v)   { max_rel_v   = rel_v;   }
        if (rel_w   > max_rel_w)   { max_rel_w   = rel_w;   }
        if (rel_p   > max_rel_p)   { max_rel_p   = rel_p;   }
        n_samples++;
      }

      for (int i = 0; i < 3; i++) { gkyl_free(spatial_metric[i]); gkyl_free(extrinsic_curvature[i]); gkyl_free(shift_der[i]); }
      for (int i = 0; i < 3; i++) { for (int j = 0; j < 3; j++) gkyl_free(spatial_metric_der[i][j]); gkyl_free(spatial_metric_der[i]); }
      gkyl_free(spatial_metric); gkyl_free(extrinsic_curvature); gkyl_free(lapse_der);
      gkyl_free(shift_der); gkyl_free(spatial_metric_der); gkyl_free(shift);
    }
  }

  double worst_rel = fmax(fmax(max_rel_rho, fmax(max_rel_u, max_rel_v)),
                          fmax(max_rel_w, max_rel_p));
  TEST_CHECK_( worst_rel < STRINGENT_REL_TOL,
    "[%s] n=%d (skip excise=%d super=%d) "
    "max rel: Δρ=%.3e Δu=%.3e Δv=%.3e Δw=%.3e Δp=%.3e   worst state=%d pos=%d",
    label, n_samples, n_skipped_excise, n_skipped_super,
    max_rel_rho, max_rel_u, max_rel_v, max_rel_w, max_rel_p, worst_s, worst_p);

  gkyl_wv_eqn_release(gr_euler);
}

void
test_prim_vars_stringent_minkowski()
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  run_prim_vars_stringent(spacetime, "Minkowski");
  gkyl_gr_spacetime_release(spacetime);
}

void
test_prim_vars_stringent_schwarzschild()
{
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_prim_vars_stringent(spacetime, "Schwarzschild a=0");
  gkyl_gr_spacetime_release(spacetime);
}

void
test_prim_vars_stringent_kerr_mild()
{
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.5, 0.0, 0.0, 0.0);
  run_prim_vars_stringent(spacetime, "Kerr a=0.5");
  gkyl_gr_spacetime_release(spacetime);
}

void
test_prim_vars_stringent_kerr_extreme()
{
  // a=0.99 puts the outer horizon at r_+ ≈ M·(1 + sqrt(1−0.9801)) ≈ 0.114
  // for M=0.1. Frame-dragging is strong; β^i is large near the horizon.
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.99, 0.0, 0.0, 0.0);
  run_prim_vars_stringent(spacetime, "Kerr a=0.99");
  gkyl_gr_spacetime_release(spacetime);
}

#include <stdio.h>

// Fill q[5..70] with the spacetime data at (x,y,z) from the given spacetime
// callbacks. Returns the excision flag (positive = not excised).
static double
fill_spacetime_block_71(struct gkyl_gr_spacetime *spacetime,
  double x, double y, double z, double q[71])
{
  double spatial_det, lapse;
  bool in_excision_region;
  double *shift = gkyl_malloc(sizeof(double[3]));
  double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
  double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
    extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
  }
  double *lapse_der = gkyl_malloc(sizeof(double[3]));
  double **shift_der = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) shift_der[i] = gkyl_malloc(sizeof(double[3]));
  double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));
    for (int j = 0; j < 3; j++) spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
  }

  spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, z, &spatial_det);
  spacetime->lapse_function_func(spacetime, 0.0, x, y, z, &lapse);
  spacetime->shift_vector_func(spacetime, 0.0, x, y, z, &shift);
  spacetime->excision_region_func(spacetime, 0.0, x, y, z, &in_excision_region);
  spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, z, &spatial_metric);
  spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, z,
    pow(10.0,-8.0), pow(10.0,-8.0), pow(10.0,-8.0), &extrinsic_curvature);
  spacetime->lapse_function_der_func(spacetime, 0.0, x, y, z,
    pow(10.0,-8.0), pow(10.0,-8.0), pow(10.0,-8.0), &lapse_der);
  spacetime->shift_vector_der_func(spacetime, 0.0, x, y, z,
    pow(10.0,-8.0), pow(10.0,-8.0), pow(10.0,-8.0), &shift_der);
  spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, z,
    pow(10.0,-8.0), pow(10.0,-8.0), pow(10.0,-8.0), &spatial_metric_der);

  q[5] = lapse;
  q[6] = shift[0]; q[7] = shift[1]; q[8] = shift[2];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) q[9 + 3*i + j] = spatial_metric[i][j];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) q[18 + 3*i + j] = extrinsic_curvature[i][j];
  q[27] = in_excision_region ? -1.0 : 1.0;
  q[28] = lapse_der[0]; q[29] = lapse_der[1]; q[30] = lapse_der[2];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) q[31 + 3*i + j] = shift_der[i][j];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++) q[40 + 9*i + 3*j + k] = spatial_metric_der[i][j][k];
  q[67] = 0.0; q[68] = x; q[69] = y; q[70] = z;

  double flag = q[27];
  // Also expose spatial_det and metric via the caller using q[9..17].
  (void)spatial_det;  // q[9..17] carries γ_ij; det is derivable.

  for (int i = 0; i < 3; i++) { gkyl_free(spatial_metric[i]); gkyl_free(extrinsic_curvature[i]); gkyl_free(shift_der[i]); }
  for (int i = 0; i < 3; i++) { for (int j = 0; j < 3; j++) gkyl_free(spatial_metric_der[i][j]); gkyl_free(spatial_metric_der[i]); }
  gkyl_free(spatial_metric); gkyl_free(extrinsic_curvature); gkyl_free(lapse_der);
  gkyl_free(shift_der); gkyl_free(spatial_metric_der); gkyl_free(shift);
  return flag;
}

// Combined 3D ultra-stress sweep of prim_vars recovery. Sweeps the full
// cross-product of W × ρ × (p/ρ) from prim_vars_stringent_data.h at a
// specified spacetime cell — Minkowski origin for the clean reference,
// near-horizon Schwarzschild and Kerr to confirm the recovery also works
// in strong-field γ.
//
// Velocity is constructed in the direction (1,1,1) (in coord components)
// and scaled by 1/√(γ_ij·δ^i·δ^j) so the resulting γ-norm hits the target
// |v|² = 1 − 1/W². The scaling stays well-conditioned for high W because
// 1/W² is small even at W=1e4.
//
// W²·ε scaling. The (τ+D)² − |S|² discriminant in the Eulderink-Mellema
// recovery polynomial loses ~2·log₁₀(W) digits of precision because the
// leading W⁴ terms cancel. Empirically: rel_err ≈ W²·machine_ε. At W=1e4
// this is ~1e-8 for moderate p/ρ. The ceiling clamps W to 1e4 in
// flux/max_abs_speed (the 1e-8 v_sq floor) which keeps the production
// code in the reliable regime.
static void
run_ultra_combined(struct gkyl_gr_spacetime *spacetime, const char *label,
  double x, double y, double z)
{
  double gas_gamma = 5.0 / 3.0;
  struct gkyl_wv_eqn *gr_euler = gkyl_wv_gr_euler_new(
    gas_gamma, GKYL_STATIC_GAUGE, 0, spacetime, false);

  // Populate the spacetime block at the test position.
  double q[71] = { 0 };
  double excision_flag = fill_spacetime_block_71(spacetime, x, y, z, q);
  TEST_CHECK_( excision_flag > 0.0,
    "[%s] test position (%.3f, %.3f, %.3f) is inside the excision region; "
    "choose a position outside the horizon", label, x, y, z );
  if (excision_flag <= 0.0) {
    gkyl_wv_eqn_release(gr_euler);
    return;
  }

  // Extract γ_ij at the cell and compute the velocity-direction scaling.
  // We pick the contravariant velocity direction (1,1,1) in coord components,
  // then scale by 1/sqrt(γ_ij·(1,1,1)^i·(1,1,1)^j) so |v|²_γ = 1 − 1/W² for
  // any target W.
  double g[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      g[i][j] = q[9 + 3*i + j];
  double g_sum = 0.0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      g_sum += g[i][j];
  double spatial_det =
    (g[0][0] * (g[1][1]*g[2][2] - g[2][1]*g[1][2]))
  - (g[0][1] * (g[1][0]*g[2][2] - g[1][2]*g[2][0]))
  + (g[0][2] * (g[1][0]*g[2][1] - g[1][1]*g[2][0]));

  const double PRIM_FLOOR = 1.0e-8;

  int n_total = 0, n_skipped_floor = 0;
  double max_rel_overall = 0.0;
  struct { double W, rho, p, rel; } top_failures[5];
  for (int i = 0; i < 5; i++) top_failures[i].rel = 0.0;

  fprintf(stdout,
    "\n[Combined ultra W×ρ×p/ρ sweep — %s @ (%.3f,%.3f,%.3f)]"
    "  %d × %d × %d = %d points\n",
    label, x, y, z,
    N_ULTRA_COMBINED_W, N_ULTRA_COMBINED_RHO, N_ULTRA_COMBINED_PR,
    N_ULTRA_COMBINED_W * N_ULTRA_COMBINED_RHO * N_ULTRA_COMBINED_PR);
  fprintf(stdout,
    "       W         ρ           p          rel ρ      rel v      rel p      rel W      worst\n");

  for (int wi = 0; wi < N_ULTRA_COMBINED_W; wi++) {
    double W_in = ultra_combined_W[wi];
    double v_sq_target = 1.0 - 1.0 / (W_in * W_in);   // well-conditioned
    double scale = sqrt(v_sq_target / g_sum);
    double v_x = scale, v_y = scale, v_z = scale;     // contravariant v^i

    for (int ri = 0; ri < N_ULTRA_COMBINED_RHO; ri++) {
      double rho_in = ultra_combined_rho[ri];
      for (int pi = 0; pi < N_ULTRA_COMBINED_PR; pi++) {
        double p_in = rho_in * ultra_combined_p_over_rho[pi];

        if (rho_in < PRIM_FLOOR || p_in < PRIM_FLOOR) {
          n_skipped_floor++;
          fprintf(stdout, "  %8.1e %10.1e %12.1e    SKIP (below %.0e floor)\n",
            W_in, rho_in, p_in, PRIM_FLOOR);
          continue;
        }

        double h_in = 1.0 + (p_in / rho_in) * (gas_gamma / (gas_gamma - 1.0));
        double rhW2 = rho_in * h_in * (W_in * W_in);
        double sqrt_det = sqrt(spatial_det);

        // Convention B: q[i+1] := √γ · ρhW² · v^i.
        q[0] = sqrt_det * rho_in * W_in;
        q[1] = sqrt_det * rhW2 * v_x;
        q[2] = sqrt_det * rhW2 * v_y;
        q[3] = sqrt_det * rhW2 * v_z;
        q[4] = sqrt_det * (rhW2 - p_in - (rho_in * W_in));

        double prims[71];
        gkyl_gr_euler_prim_vars(gas_gamma, q, prims);

        // Recovered W from recovered velocity via the curved γ-norm.
        double v_sq_rec = 0.0;
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++)
            v_sq_rec += g[i][j] * prims[1 + i] * prims[1 + j];
        if (v_sq_rec >= 1.0) v_sq_rec = 1.0 - 1.0e-300;
        double W_rec = 1.0 / sqrt(1.0 - v_sq_rec);

        double rel_rho = fabs(prims[0] - rho_in) / rho_in;
        double rel_v   = fabs(prims[1] - v_x)    / v_x;
        double rel_p   = fabs(prims[4] - p_in)   / p_in;
        double rel_W   = fabs(W_rec    - W_in)   / W_in;
        double worst   = fmax(fmax(rel_rho, rel_v), fmax(rel_p, rel_W));

        fprintf(stdout,
          "  %8.1e %10.1e %12.1e   %8.2e  %8.2e  %8.2e  %8.2e   %8.2e\n",
          W_in, rho_in, p_in, rel_rho, rel_v, rel_p, rel_W, worst);

        n_total++;
        if (worst > max_rel_overall) max_rel_overall = worst;

        for (int t = 0; t < 5; t++) {
          if (worst > top_failures[t].rel) {
            for (int s = 4; s > t; s--) top_failures[s] = top_failures[s-1];
            top_failures[t].W   = W_in;
            top_failures[t].rho = rho_in;
            top_failures[t].p   = p_in;
            top_failures[t].rel = worst;
            break;
          }
        }
      }
    }
  }

  fprintf(stdout, "\n  Summary (%s):\n", label);
  fprintf(stdout, "    γ_xx=%.3f  γ_yy=%.3f  γ_zz=%.3f  α=%.3f  √γ=%.3f\n",
    g[0][0], g[1][1], g[2][2], q[5], sqrt(spatial_det));
  fprintf(stdout, "    %d total points, %d skipped (below %.0e floor)\n",
    N_ULTRA_COMBINED_W * N_ULTRA_COMBINED_RHO * N_ULTRA_COMBINED_PR,
    n_skipped_floor, PRIM_FLOOR);
  fprintf(stdout, "    Max relative recovery error: %.3e\n", max_rel_overall);
  fprintf(stdout, "    Top failures:\n");
  for (int t = 0; t < 5; t++) {
    if (top_failures[t].rel > 0.0) {
      fprintf(stdout, "      W=%8.1e ρ=%10.1e p=%12.1e  rel=%.3e\n",
        top_failures[t].W, top_failures[t].rho, top_failures[t].p,
        top_failures[t].rel);
    }
  }
  fprintf(stdout, "\n");

  // Empirical reliability at W ≤ 1e4 (with the s² floor at 1e-10):
  //   - "well-conditioned" combos (p/ρ ≥ 1e-3): ≤ 1e-5 worst case
  //   - extreme p/ρ = 1e-6 at W=1e4: ~1.3% (algorithmic-cancellation limit)
  // The previous low-ρ × low-p × low-W corner pathology (rel ~35×) was the
  // s² floor at 1e-8 activating spuriously; tightening to 1e-10 eliminated
  // it across all spacetimes.
  TEST_CHECK( isfinite(max_rel_overall) );
  TEST_CHECK_( max_rel_overall < 2.0e-2,
    "[%s] max relative recovery error %.3e exceeds 2%% threshold", label, max_rel_overall );

  gkyl_wv_eqn_release(gr_euler);
}

void
test_prim_vars_ultra_combined_minkowski()
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);
  run_ultra_combined(spacetime, "Minkowski origin", 0.0, 0.0, 0.0);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_prim_vars_ultra_combined_schwarzschild_near_horizon()
{
  // Schwarzschild M=0.1: r_+ = 2M = 0.2. Sample at r ≈ 0.27 (just outside)
  // off-axis so γ has off-diagonal entries in Cartesian Kerr-Schild.
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);
  run_ultra_combined(spacetime, "Schwarzschild near-horizon", 0.25, 0.10, 0.0);
  gkyl_gr_spacetime_release(spacetime);
}

void
test_prim_vars_ultra_combined_kerr_near_horizon()
{
  // Kerr M=0.1, a=0.9: r_+ ≈ 0.144. Sample at r ≈ 0.27 (close to ergosphere
  // boundary, large β^i from frame-dragging).
  struct gkyl_gr_spacetime *spacetime =
    gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);
  run_ultra_combined(spacetime, "Kerr a=0.9 near-horizon", 0.25, 0.10, 0.0);
  gkyl_gr_spacetime_release(spacetime);
}

TEST_LIST = {
  { "gr_euler_basic_minkowski", test_gr_euler_basic_minkowski },
  { "gr_euler_basic_schwarzschild", test_gr_euler_basic_schwarzschild },
  { "gr_euler_basic_kerr", test_gr_euler_basic_kerr },
  { "gr_euler_waves_minkowski", test_gr_euler_waves_minkowski },
  { "gr_euler_waves_schwarzschild", test_gr_euler_waves_schwarzschild },
  { "gr_euler_waves_kerr", test_gr_euler_waves_kerr },
  { "prim_vars_stringent_minkowski", test_prim_vars_stringent_minkowski },
  { "prim_vars_stringent_schwarzschild", test_prim_vars_stringent_schwarzschild },
  { "prim_vars_stringent_kerr_mild", test_prim_vars_stringent_kerr_mild },
  { "prim_vars_stringent_kerr_extreme", test_prim_vars_stringent_kerr_extreme },
  { "prim_vars_ultra_combined_minkowski",                  test_prim_vars_ultra_combined_minkowski },
  { "prim_vars_ultra_combined_schwarzschild_near_horizon", test_prim_vars_ultra_combined_schwarzschild_near_horizon },
  { "prim_vars_ultra_combined_kerr_near_horizon",          test_prim_vars_ultra_combined_kerr_near_horizon },
  { NULL, NULL },
};
