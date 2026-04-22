#include <acutest.h>

#include <gkyl_util.h>
#include <gkyl_gr_minkowski.h>
#include <gkyl_gr_blackhole.h>
#include <gkyl_gr_neutronstar.h>
#include <gkyl_gr_brill_lindquist.h>
#include "tov_solver.h"
#include <gkyl_gr_tov_spacetime.h>

void
test_gr_minkowski()
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_minkowski_new(false);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
      double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));

      for (int i = 0; i < 3; i++) {
        spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
        inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
      }

      double **spacetime_metric = gkyl_malloc(sizeof(double*[4]));
      double **inv_spacetime_metric = gkyl_malloc(sizeof(double*[4]));

      for (int i = 0; i < 4; i++) {
        spacetime_metric[i] = gkyl_malloc(sizeof(double[4]));
        inv_spacetime_metric[i] = gkyl_malloc(sizeof(double[4]));
      }

      spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);
      spacetime->spacetime_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spacetime_metric);

      spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &inv_spatial_metric);
      spacetime->spacetime_inv_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &inv_spacetime_metric);

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          if (i == j) {
            TEST_CHECK( gkyl_compare(spatial_metric[i][j], 1.0, 1e-10) );
            TEST_CHECK( gkyl_compare(inv_spatial_metric[i][j], 1.0, 1e-10) );
          }
          else {
            TEST_CHECK( gkyl_compare(spatial_metric[i][j], 0.0, 1e-10) );
            TEST_CHECK( gkyl_compare(inv_spatial_metric[i][j], 0.0, 1e-10) );
          }
        }
      }

      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
          if (i == j) {
            if (i == 0) {
              TEST_CHECK( gkyl_compare(spacetime_metric[i][j], -1.0, 1e-10) );
              TEST_CHECK( gkyl_compare(inv_spacetime_metric[i][j], -1.0, 1e-10) );
            }
            else {
              TEST_CHECK( gkyl_compare(spacetime_metric[i][j], 1.0, 1e-10) );
              TEST_CHECK( gkyl_compare(inv_spacetime_metric[i][j], 1.0, 1e-10) );
            }
          }
          else {
            TEST_CHECK( gkyl_compare(spacetime_metric[i][j], 0.0, 1e-10) );
            TEST_CHECK( gkyl_compare(inv_spacetime_metric[i][j], 0.0, 1e-10) );
          }
        }
      }

      double spatial_metric_det;
      double spacetime_metric_det;

      spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &spatial_metric_det);
      spacetime->spacetime_metric_det_func(spacetime, 0.0, x, y, 0.0, &spacetime_metric_det);

      TEST_CHECK( gkyl_compare(spatial_metric_det, 1.0, 1e-10) );
      TEST_CHECK( gkyl_compare(spacetime_metric_det, -1.0, 1e-10) );

      double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
      double ***spatial_christoffel = gkyl_malloc(sizeof(double**[3]));

      for (int i = 0; i < 3; i++) {
        spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));
        spatial_christoffel[i] = gkyl_malloc(sizeof(double*[3]));

        for (int j = 0; j < 3; j++) {
          spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
          spatial_christoffel[i][j] = gkyl_malloc(sizeof(double[3]));
        }
      }

      double ***spacetime_metric_der = gkyl_malloc(sizeof(double**[4]));
      double ***spacetime_christoffel = gkyl_malloc(sizeof(double**[4]));

      for (int i = 0; i < 4; i++) {
        spacetime_metric_der[i] = gkyl_malloc(sizeof(double*[4]));
        spacetime_christoffel[i] = gkyl_malloc(sizeof(double*[4]));

        for (int j = 0; j < 4; j++) {
          spacetime_metric_der[i][j] = gkyl_malloc(sizeof(double[4]));
          spacetime_christoffel[i][j] = gkyl_malloc(sizeof(double[4]));
        }
      }

      spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spatial_metric_der);
      spacetime->spatial_christoffel_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spatial_christoffel);

      spacetime->spacetime_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_metric_der);
      spacetime->spacetime_christoffel_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_christoffel);

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            TEST_CHECK( gkyl_compare(spatial_metric_der[i][j][k], 0.0, 1e-10) );
            TEST_CHECK( gkyl_compare(spatial_christoffel[i][j][k], 0.0, 1e-10) );
          }
        }
      }

      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
          for (int k = 0; k < 4; k++) {
            TEST_CHECK( gkyl_compare(spacetime_metric_der[i][j][k], 0.0, 1e-10) );
            TEST_CHECK( gkyl_compare(spacetime_christoffel[i][j][k], 0.0, 1e-10) );
          }
        }
      }

      double lapse_function;
      double *shift_vector = gkyl_malloc(sizeof(double[3]));

      spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse_function);
      spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift_vector);

      TEST_CHECK( gkyl_compare(lapse_function, 1.0, 1e-10) );
      for (int i = 0; i < 3; i++) {
        TEST_CHECK( gkyl_compare(shift_vector[i], 0.0, 1e-10) );
      }

      double *lapse_function_der = gkyl_malloc(sizeof(double[3]));
      double **shift_vector_der = gkyl_malloc(sizeof(double*[3]));
      double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));

      for (int i = 0; i < 3; i++) {
        shift_vector_der[i] = gkyl_malloc(sizeof(double[3]));
        extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
      }

      spacetime->lapse_function_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &lapse_function_der);
      spacetime->shift_vector_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &shift_vector_der);
      spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &extrinsic_curvature);

      for (int i = 0; i < 3; i++) {
        TEST_CHECK( gkyl_compare(lapse_function_der[i], 0.0, 1e-10) );
      }

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          TEST_CHECK( gkyl_compare(shift_vector_der[i][j], 0.0, 1e-10) );
          TEST_CHECK( gkyl_compare(extrinsic_curvature[i][j], 0.0, 1e-10) );
        }
      }

      double ****spatial_riemann_tensor = gkyl_malloc(sizeof(double***[3]));
      double ****spatial_weyl_tensor = gkyl_malloc(sizeof(double***[3]));
      double **spatial_ricci_tensor = gkyl_malloc(sizeof(double*[3]));
      double spatial_ricci_scalar;

      for (int i = 0; i < 3; i++) {
        spatial_riemann_tensor[i] = gkyl_malloc(sizeof(double**[3]));
        spatial_weyl_tensor[i] = gkyl_malloc(sizeof(double**[3]));
        spatial_ricci_tensor[i] = gkyl_malloc(sizeof(double[3]));

        for (int j = 0; j < 3; j++) {
          spatial_riemann_tensor[i][j] = gkyl_malloc(sizeof(double*[3]));
          spatial_weyl_tensor[i][j] = gkyl_malloc(sizeof(double*[3]));

          for (int k = 0; k < 3; k++) {
            spatial_riemann_tensor[i][j][k] = gkyl_malloc(sizeof(double[3]));
            spatial_weyl_tensor[i][j][k] = gkyl_malloc(sizeof(double[3]));
          }
        }
      }

      spacetime->spatial_riemann_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spatial_riemann_tensor);
      spacetime->spatial_weyl_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spatial_weyl_tensor);
      spacetime->spatial_ricci_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spatial_ricci_tensor);
      spacetime->spatial_ricci_scalar_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spatial_ricci_scalar);

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          TEST_CHECK( gkyl_compare(spatial_ricci_tensor[i][j], 0.0, 1e-10) );

          for (int k = 0; k < 3; k++) {
            for (int l = 0; l < 3; l++) {
              TEST_CHECK( gkyl_compare(spatial_riemann_tensor[i][j][k][l], 0.0, 1e-10) );
              TEST_CHECK( gkyl_compare(spatial_weyl_tensor[i][j][k][l], 0.0, 1e-10) );
            }
          }
        }
      }

      TEST_CHECK( gkyl_compare(spatial_ricci_scalar, 0.0, 1e-10) );

      double ****spacetime_riemann_tensor = gkyl_malloc(sizeof(double***[4]));
      double ****spacetime_weyl_tensor = gkyl_malloc(sizeof(double***[4]));
      double **spacetime_ricci_tensor = gkyl_malloc(sizeof(double*[4]));
      double spacetime_ricci_scalar;

      for (int i = 0; i < 4; i++) {
        spacetime_riemann_tensor[i] = gkyl_malloc(sizeof(double**[4]));
        spacetime_weyl_tensor[i] = gkyl_malloc(sizeof(double**[4]));
        spacetime_ricci_tensor[i] = gkyl_malloc(sizeof(double[4]));

        for (int j = 0; j < 4; j++) {
          spacetime_riemann_tensor[i][j] = gkyl_malloc(sizeof(double*[4]));
          spacetime_weyl_tensor[i][j] = gkyl_malloc(sizeof(double*[4]));

          for (int k = 0; k < 4; k++) {
            spacetime_riemann_tensor[i][j][k] = gkyl_malloc(sizeof(double[4]));
            spacetime_weyl_tensor[i][j][k] = gkyl_malloc(sizeof(double[4]));
          }
        }
      }

      spacetime->spacetime_riemann_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_riemann_tensor);
      spacetime->spacetime_weyl_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_weyl_tensor);
      spacetime->spacetime_ricci_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_ricci_tensor);
      spacetime->spacetime_ricci_scalar_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_ricci_scalar);

      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
          TEST_CHECK( gkyl_compare(spacetime_ricci_tensor[i][j], 0.0, 1e-10) );

          for (int k = 0; k < 4; k++) {
            for (int l = 0; l < 4; l++) {
              TEST_CHECK( gkyl_compare(spacetime_riemann_tensor[i][j][k][l], 0.0, 1e-10) );
              TEST_CHECK( gkyl_compare(spacetime_weyl_tensor[i][j][k][l], 0.0, 1e-10) );
            }
          }
        }
      }

      TEST_CHECK( gkyl_compare(spacetime_ricci_scalar, 0.0, 1e-10) );

      bool in_excision_region;
      spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);

      TEST_CHECK( (in_excision_region == false) );

      for (int i = 0; i < 3; i++) {
        gkyl_free(spatial_metric[i]);
        gkyl_free(inv_spatial_metric[i]);
        gkyl_free(extrinsic_curvature[i]);
        gkyl_free(shift_vector_der[i]);
        gkyl_free(spatial_ricci_tensor[i]);

        for (int j = 0; j < 3; j++) {
          gkyl_free(spatial_metric_der[i][j]);
          gkyl_free(spatial_christoffel[i][j]);

          for (int k = 0; k < 3; k++) {
            gkyl_free(spatial_riemann_tensor[i][j][k]);
            gkyl_free(spatial_weyl_tensor[i][j][k]);
          }
          gkyl_free(spatial_riemann_tensor[i][j]);
          gkyl_free(spatial_weyl_tensor[i][j]);
        }
        gkyl_free(spatial_metric_der[i]);
        gkyl_free(spatial_christoffel[i]);
        gkyl_free(spatial_riemann_tensor[i]);
        gkyl_free(spatial_weyl_tensor[i]);
      }
      gkyl_free(spatial_metric);
      gkyl_free(inv_spatial_metric);
      gkyl_free(spatial_metric_der);
      gkyl_free(spatial_christoffel);
      gkyl_free(extrinsic_curvature);
      gkyl_free(shift_vector_der);
      gkyl_free(lapse_function_der);
      gkyl_free(spatial_riemann_tensor);
      gkyl_free(spatial_weyl_tensor);
      gkyl_free(spatial_ricci_tensor);

      for (int i = 0; i < 4; i++) {
        gkyl_free(spacetime_metric[i]);
        gkyl_free(inv_spacetime_metric[i]);
        gkyl_free(spacetime_ricci_tensor[i]);

        for (int j = 0; j < 4; j++) {
          gkyl_free(spacetime_metric_der[i][j]);
          gkyl_free(spacetime_christoffel[i][j]);
          
          for (int k = 0; k < 4; k++) {
            gkyl_free(spacetime_riemann_tensor[i][j][k]);
            gkyl_free(spacetime_weyl_tensor[i][j][k]);
          }
          gkyl_free(spacetime_riemann_tensor[i][j]);
          gkyl_free(spacetime_weyl_tensor[i][j]);
        }
        gkyl_free(spacetime_metric_der[i]);
        gkyl_free(spacetime_christoffel[i]);
        gkyl_free(spacetime_riemann_tensor[i]);
        gkyl_free(spacetime_weyl_tensor[i]);
      }
      gkyl_free(spacetime_metric);
      gkyl_free(inv_spacetime_metric);
      gkyl_free(spacetime_metric_der);
      gkyl_free(spacetime_christoffel);
      gkyl_free(spacetime_riemann_tensor);
      gkyl_free(spacetime_weyl_tensor);
      gkyl_free(spacetime_ricci_tensor);
    }
  }

  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_schwarzschild()
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.0, 0.0, 0.0, 0.0);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      if (sqrt((x * x) + (y * y)) > 0.2) {
        double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
        double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
        double **spatial_metric_prod = gkyl_malloc(sizeof(double*[3]));

        for (int i = 0; i < 3; i++) {
          spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
          inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
          spatial_metric_prod[i] = gkyl_malloc(sizeof(double[3]));

          for (int j = 0; j < 3; j++) {
            spatial_metric_prod[i][j] = 0.0;
          }
        }

        double **spacetime_metric = gkyl_malloc(sizeof(double*[4]));
        double **inv_spacetime_metric = gkyl_malloc(sizeof(double*[4]));
        double **spacetime_metric_prod = gkyl_malloc(sizeof(double*[4]));

        for (int i = 0; i < 4; i++) {
          spacetime_metric[i] = gkyl_malloc(sizeof(double[4]));
          inv_spacetime_metric[i] = gkyl_malloc(sizeof(double[4]));
          spacetime_metric_prod[i] = gkyl_malloc(sizeof(double[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_metric_prod[i][j] = 0.0;
          }
        }

        spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);
        spacetime->spacetime_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spacetime_metric);

        spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &inv_spatial_metric);
        spacetime->spacetime_inv_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &inv_spacetime_metric);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              spatial_metric_prod[i][j] += spatial_metric[i][k] * inv_spatial_metric[k][j];
            }

            if (i == j) {
              TEST_CHECK( gkyl_compare(spatial_metric_prod[i][j], 1.0, 1e-10) );
            }
            else {
              TEST_CHECK( gkyl_compare(spatial_metric_prod[i][j], 0.0, 1e-10) );
            }
          }
        }

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
              spacetime_metric_prod[i][j] += spacetime_metric[i][k] * inv_spacetime_metric[k][j];
            }

            if (i == j) {
              TEST_CHECK( gkyl_compare(spacetime_metric_prod[i][j], 1.0, 1e-10) );
            }
            else {
              TEST_CHECK( gkyl_compare(spacetime_metric_prod[i][j], 0.0, 1e-10) );
            }
          }
        }

        double spatial_metric_det;
        double spacetime_metric_det;
        double lapse_function;

        spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &spatial_metric_det);
        spacetime->spacetime_metric_det_func(spacetime, 0.0, x, y, 0.0, &spacetime_metric_det);
        spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse_function);

        TEST_CHECK( gkyl_compare(sqrt(-spacetime_metric_det), lapse_function * sqrt(spatial_metric_det), 1e-10) );

        double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
        double ***spatial_christoffel = gkyl_malloc(sizeof(double**[3]));
        double ***spatial_metric_cov_der = gkyl_malloc(sizeof(double**[3]));

        for (int i = 0; i < 3; i++) {
          spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));
          spatial_christoffel[i] = gkyl_malloc(sizeof(double*[3]));
          spatial_metric_cov_der[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
            spatial_christoffel[i][j] = gkyl_malloc(sizeof(double[3]));
            spatial_metric_cov_der[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spatial_metric_der);
        spacetime->spatial_christoffel_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spatial_christoffel);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              spatial_metric_cov_der[i][j][k] = spatial_metric_der[i][j][k];

              for (int l = 0; l < 3; l++) {
                spatial_metric_cov_der[i][j][k] -= spatial_christoffel[l][i][j] * spatial_metric[l][k];
                spatial_metric_cov_der[i][j][k] -= spatial_christoffel[l][i][k] * spatial_metric[j][l];
              }

              TEST_CHECK( gkyl_compare(spatial_metric_cov_der[i][j][k], 0.0, 1e-10) );
            }
          }
        }

        double ***spacetime_metric_der = gkyl_malloc(sizeof(double**[4]));
        double ***spacetime_christoffel = gkyl_malloc(sizeof(double**[4]));
        double ***spacetime_metric_cov_der = gkyl_malloc(sizeof(double**[4]));

        for (int i = 0; i < 4; i++) {
          spacetime_metric_der[i] = gkyl_malloc(sizeof(double*[4]));
          spacetime_christoffel[i] = gkyl_malloc(sizeof(double*[4]));
          spacetime_metric_cov_der[i] = gkyl_malloc(sizeof(double*[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_metric_der[i][j] = gkyl_malloc(sizeof(double[4]));
            spacetime_christoffel[i][j] = gkyl_malloc(sizeof(double[4]));
            spacetime_metric_cov_der[i][j] = gkyl_malloc(sizeof(double[4]));
          }
        }

        spacetime->spacetime_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_metric_der);
        spacetime->spacetime_christoffel_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_christoffel);

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
              spacetime_metric_cov_der[i][j][k] = spacetime_metric_der[i][j][k];

              for (int l = 0; l < 4; l++) {
                spacetime_metric_cov_der[i][j][k] -= spacetime_christoffel[l][i][j] * spacetime_metric[l][k];
                spacetime_metric_cov_der[i][j][k] -= spacetime_christoffel[l][i][k] * spacetime_metric[j][l];
              }

              TEST_CHECK( gkyl_compare(spacetime_metric_cov_der[i][j][k], 0.0, 1e-10) );
            }
          }
        }

        double **shift_vector_der = gkyl_malloc(sizeof(double*[3]));
        double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
        double **shift_vector_cov_der = gkyl_malloc(sizeof(double*[3]));
        double **shift_covector_cov_der = gkyl_malloc(sizeof(double*[3]));

        for (int i = 0; i < 3; i++) {
          shift_vector_der[i] = gkyl_malloc(sizeof(double[3]));
          extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
          shift_vector_cov_der[i] = gkyl_malloc(sizeof(double[3]));
          shift_covector_cov_der[i] = gkyl_malloc(sizeof(double[3]));

          for (int j = 0; j < 3; j++){
            shift_covector_cov_der[i][j] = 0.0;
          }
        }

        spacetime->shift_vector_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &shift_vector_der);
        spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &extrinsic_curvature);

        double *shift_vector = gkyl_malloc(sizeof(double[3]));
        spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift_vector);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            shift_vector_cov_der[i][j] = shift_vector_der[i][j];

            for (int k = 0; k < 3; k++) {
              shift_vector_cov_der[i][j] += spatial_christoffel[j][i][k] * shift_vector[k];
            }
          }
        }

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              shift_covector_cov_der[i][j] += spatial_metric[j][k] * shift_vector_cov_der[i][k];
            }
          }
        }

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            TEST_CHECK( gkyl_compare(2.0 * lapse_function * extrinsic_curvature[i][j], -(shift_covector_cov_der[j][i] + shift_covector_cov_der[i][j]), 1e-6) );
          }
        }

        double **spacetime_ricci_tensor = gkyl_malloc(sizeof(double*[4]));
        double spacetime_ricci_scalar;
        for (int i = 0; i < 4; i++) {
          spacetime_ricci_tensor[i] = gkyl_malloc(sizeof(double[4]));
        }

        spacetime->spacetime_ricci_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_ricci_tensor);
        spacetime->spacetime_ricci_scalar_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_ricci_scalar);

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            TEST_CHECK( gkyl_compare(spacetime_ricci_tensor[i][j], 0.0, 1e-3) );
          }
        }

        TEST_CHECK( gkyl_compare(spacetime_ricci_scalar, 0.0, 1e-2) );

        double ****spacetime_riemann_tensor = gkyl_malloc(sizeof(double***[4]));
        double ****covariant_spacetime_riemann_tensor = gkyl_malloc(sizeof(double***[4]));
        double ****spacetime_weyl_tensor = gkyl_malloc(sizeof(double***[4]));

        for (int i = 0; i < 4; i++) {
          spacetime_riemann_tensor[i] = gkyl_malloc(sizeof(double**[4]));
          covariant_spacetime_riemann_tensor[i] = gkyl_malloc(sizeof(double**[4]));
          spacetime_weyl_tensor[i] = gkyl_malloc(sizeof(double**[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_riemann_tensor[i][j] = gkyl_malloc(sizeof(double*[4]));
            covariant_spacetime_riemann_tensor[i][j] = gkyl_malloc(sizeof(double*[4]));
            spacetime_weyl_tensor[i][j] = gkyl_malloc(sizeof(double*[4]));

            for (int k = 0; k < 4; k++) {
              spacetime_riemann_tensor[i][j][k] = gkyl_malloc(sizeof(double[4]));
              covariant_spacetime_riemann_tensor[i][j][k] = gkyl_malloc(sizeof(double[4]));
              spacetime_weyl_tensor[i][j][k] = gkyl_malloc(sizeof(double[4]));

              for (int l = 0; l < 4; l++) {
                covariant_spacetime_riemann_tensor[i][j][k][l] = 0.0;
              }
            }
          }
        }

        spacetime->spacetime_riemann_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_riemann_tensor);
        spacetime->spacetime_weyl_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_weyl_tensor);

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
              for (int l = 0; l < 4; l++) {
                for (int m = 0; m < 4; m++) {
                  covariant_spacetime_riemann_tensor[i][j][k][l] += spacetime_metric[i][m] * spacetime_riemann_tensor[m][j][k][l];
                }
              }
            }
          }
        }

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
              for (int l = 0; l < 4; l++) {
                TEST_CHECK( gkyl_compare(covariant_spacetime_riemann_tensor[i][j][k][l], spacetime_weyl_tensor[i][j][k][l], 1e-2) );
              }
            }
          }
        }

        bool in_excision_region;
        spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);

        TEST_CHECK( (in_excision_region == false) );

        for (int i = 0; i < 3; i++) {
          gkyl_free(spatial_metric[i]);
          gkyl_free(inv_spatial_metric[i]);
          gkyl_free(spatial_metric_prod[i]);
          gkyl_free(shift_vector_der[i]);
          gkyl_free(extrinsic_curvature[i]);
          gkyl_free(shift_vector_cov_der[i]);
          gkyl_free(shift_covector_cov_der[i]);
          
          for (int j = 0; j < 3; j++) {
            gkyl_free(spatial_metric_der[i][j]);
            gkyl_free(spatial_christoffel[i][j]);
            gkyl_free(spatial_metric_cov_der[i][j]);
          }
          gkyl_free(spatial_metric_der[i]);
          gkyl_free(spatial_christoffel[i]);
          gkyl_free(spatial_metric_cov_der[i]);
        }
        gkyl_free(spatial_metric);
        gkyl_free(inv_spatial_metric);
        gkyl_free(spatial_metric_prod);
        gkyl_free(spatial_metric_der);
        gkyl_free(spatial_christoffel);
        gkyl_free(spatial_metric_cov_der);
        gkyl_free(shift_vector_der);
        gkyl_free(extrinsic_curvature);
        gkyl_free(shift_vector_cov_der);
        gkyl_free(shift_covector_cov_der);
        gkyl_free(shift_vector);

        for (int i = 0; i < 4; i++) {
          gkyl_free(spacetime_metric[i]);
          gkyl_free(inv_spacetime_metric[i]);
          gkyl_free(spacetime_metric_prod[i]);

          for (int j = 0; j < 4; j++) {
            gkyl_free(spacetime_metric_der[i][j]);
            gkyl_free(spacetime_christoffel[i][j]);
            gkyl_free(spacetime_metric_cov_der[i][j]);

            for (int k = 0; k < 4; k++) {
              gkyl_free(spacetime_riemann_tensor[i][j][k]);
              gkyl_free(covariant_spacetime_riemann_tensor[i][j][k]);
              gkyl_free(spacetime_weyl_tensor[i][j][k]);
            }
            gkyl_free(spacetime_riemann_tensor[i][j]);
            gkyl_free(covariant_spacetime_riemann_tensor[i][j]);
            gkyl_free(spacetime_weyl_tensor[i][j]);
          }
          gkyl_free(spacetime_metric_der[i]);
          gkyl_free(spacetime_christoffel[i]);
          gkyl_free(spacetime_metric_cov_der[i]);
          gkyl_free(spacetime_ricci_tensor[i]);
          gkyl_free(spacetime_riemann_tensor[i]);
          gkyl_free(covariant_spacetime_riemann_tensor[i]);
          gkyl_free(spacetime_weyl_tensor[i]);
        }
        gkyl_free(spacetime_metric);
        gkyl_free(inv_spacetime_metric);
        gkyl_free(spacetime_metric_prod);
        gkyl_free(spacetime_metric_der);
        gkyl_free(spacetime_christoffel);
        gkyl_free(spacetime_metric_cov_der);
        gkyl_free(spacetime_ricci_tensor);
        gkyl_free(spacetime_riemann_tensor);
        gkyl_free(covariant_spacetime_riemann_tensor);
        gkyl_free(spacetime_weyl_tensor);
      }
      else {
        bool in_excision_region;
        spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);

        TEST_CHECK( (in_excision_region == true) );
      }
    }
  }

  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_kerr()
{
  struct gkyl_gr_spacetime *spacetime = gkyl_gr_blackhole_new(false, 0.1, 0.9, 0.0, 0.0, 0.0);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      if (sqrt((x * x) + (y * y)) > 0.1 * (1.0 + sqrt(0.19))) {
        double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
        double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
        double **spatial_metric_prod = gkyl_malloc(sizeof(double*[3]));

        for (int i = 0; i < 3; i++) {
          spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
          inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
          spatial_metric_prod[i] = gkyl_malloc(sizeof(double[3]));

          for (int j = 0; j < 3; j++) {
            spatial_metric_prod[i][j] = 0.0;
          }
        }

        double **spacetime_metric = gkyl_malloc(sizeof(double*[4]));
        double **inv_spacetime_metric = gkyl_malloc(sizeof(double*[4]));
        double **spacetime_metric_prod = gkyl_malloc(sizeof(double*[4]));

        for (int i = 0; i < 4; i++) {
          spacetime_metric[i] = gkyl_malloc(sizeof(double[4]));
          inv_spacetime_metric[i] = gkyl_malloc(sizeof(double[4]));
          spacetime_metric_prod[i] = gkyl_malloc(sizeof(double[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_metric_prod[i][j] = 0.0;
          }
        }

        spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);
        spacetime->spacetime_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spacetime_metric);

        spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &inv_spatial_metric);
        spacetime->spacetime_inv_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &inv_spacetime_metric);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              spatial_metric_prod[i][j] += spatial_metric[i][k] * inv_spatial_metric[k][j];
            }

            if (i == j) {
              TEST_CHECK( gkyl_compare(spatial_metric_prod[i][j], 1.0, 1e-10) );
            }
            else {
              TEST_CHECK( gkyl_compare(spatial_metric_prod[i][j], 0.0, 1e-10) );
            }
          }
        }

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
              spacetime_metric_prod[i][j] += spacetime_metric[i][k] * inv_spacetime_metric[k][j];
            }

            if (i == j) {
              TEST_CHECK( gkyl_compare(spacetime_metric_prod[i][j], 1.0, 1e-10) );
            }
            else {
              TEST_CHECK( gkyl_compare(spacetime_metric_prod[i][j], 0.0, 1e-10) );
            }
          }
        }

        double spatial_metric_det;
        double spacetime_metric_det;
        double lapse_function;

        spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &spatial_metric_det);
        spacetime->spacetime_metric_det_func(spacetime, 0.0, x, y, 0.0, &spacetime_metric_det);
        spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse_function);

        TEST_CHECK( gkyl_compare(sqrt(-spacetime_metric_det), lapse_function * sqrt(spatial_metric_det), 1e-10) );

        double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
        double ***spatial_christoffel = gkyl_malloc(sizeof(double**[3]));
        double ***spatial_metric_cov_der = gkyl_malloc(sizeof(double**[3]));

        for (int i = 0; i < 3; i++) {
          spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));
          spatial_christoffel[i] = gkyl_malloc(sizeof(double*[3]));
          spatial_metric_cov_der[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
            spatial_christoffel[i][j] = gkyl_malloc(sizeof(double[3]));
            spatial_metric_cov_der[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spatial_metric_der);
        spacetime->spatial_christoffel_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spatial_christoffel);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              spatial_metric_cov_der[i][j][k] = spatial_metric_der[i][j][k];

              for (int l = 0; l < 3; l++) {
                spatial_metric_cov_der[i][j][k] -= spatial_christoffel[l][i][j] * spatial_metric[l][k];
                spatial_metric_cov_der[i][j][k] -= spatial_christoffel[l][i][k] * spatial_metric[j][l];
              }

              TEST_CHECK( gkyl_compare(spatial_metric_cov_der[i][j][k], 0.0, 1e-9) );
            }
          }
        }

        double ***spacetime_metric_der = gkyl_malloc(sizeof(double**[4]));
        double ***spacetime_christoffel = gkyl_malloc(sizeof(double**[4]));
        double ***spacetime_metric_cov_der = gkyl_malloc(sizeof(double**[4]));

        for (int i = 0; i < 4; i++) {
          spacetime_metric_der[i] = gkyl_malloc(sizeof(double*[4]));
          spacetime_christoffel[i] = gkyl_malloc(sizeof(double*[4]));
          spacetime_metric_cov_der[i] = gkyl_malloc(sizeof(double*[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_metric_der[i][j] = gkyl_malloc(sizeof(double[4]));
            spacetime_christoffel[i][j] = gkyl_malloc(sizeof(double[4]));
            spacetime_metric_cov_der[i][j] = gkyl_malloc(sizeof(double[4]));
          }
        }

        spacetime->spacetime_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_metric_der);
        spacetime->spacetime_christoffel_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_christoffel);

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
              spacetime_metric_cov_der[i][j][k] = spacetime_metric_der[i][j][k];

              for (int l = 0; l < 4; l++) {
                spacetime_metric_cov_der[i][j][k] -= spacetime_christoffel[l][i][j] * spacetime_metric[l][k];
                spacetime_metric_cov_der[i][j][k] -= spacetime_christoffel[l][i][k] * spacetime_metric[j][l];
              }

              TEST_CHECK( gkyl_compare(spacetime_metric_cov_der[i][j][k], 0.0, 1e-9) );
            }
          }
        }

        double **shift_vector_der = gkyl_malloc(sizeof(double*[3]));
        double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
        double **shift_vector_cov_der = gkyl_malloc(sizeof(double*[3]));
        double **shift_covector_cov_der = gkyl_malloc(sizeof(double*[3]));

        for (int i = 0; i < 3; i++) {
          shift_vector_der[i] = gkyl_malloc(sizeof(double[3]));
          extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
          shift_vector_cov_der[i] = gkyl_malloc(sizeof(double[3]));
          shift_covector_cov_der[i] = gkyl_malloc(sizeof(double[3]));

          for (int j = 0; j < 3; j++){
            shift_covector_cov_der[i][j] = 0.0;
          }
        }

        spacetime->shift_vector_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &shift_vector_der);
        spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &extrinsic_curvature);

        double *shift_vector = gkyl_malloc(sizeof(double[3]));
        spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift_vector);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            shift_vector_cov_der[i][j] = shift_vector_der[i][j];

            for (int k = 0; k < 3; k++) {
              shift_vector_cov_der[i][j] += spatial_christoffel[j][i][k] * shift_vector[k];
            }
          }
        }

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              shift_covector_cov_der[i][j] += spatial_metric[j][k] * shift_vector_cov_der[i][k];
            }
          }
        }

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            TEST_CHECK( gkyl_compare(2.0 * lapse_function * extrinsic_curvature[i][j], -(shift_covector_cov_der[j][i] + shift_covector_cov_der[i][j]), 1e-6) );
          }
        }

        double **spacetime_ricci_tensor = gkyl_malloc(sizeof(double*[4]));
        double spacetime_ricci_scalar;
        for (int i = 0; i < 4; i++) {
          spacetime_ricci_tensor[i] = gkyl_malloc(sizeof(double[4]));
        }

        spacetime->spacetime_ricci_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_ricci_tensor);
        spacetime->spacetime_ricci_scalar_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_ricci_scalar);

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            TEST_CHECK( gkyl_compare(spacetime_ricci_tensor[i][j], 0.0, 1e-3) );
          }
        }

        TEST_CHECK( gkyl_compare(spacetime_ricci_scalar, 0.0, 1e-2) );

        double ****spacetime_riemann_tensor = gkyl_malloc(sizeof(double***[4]));
        double ****covariant_spacetime_riemann_tensor = gkyl_malloc(sizeof(double***[4]));
        double ****spacetime_weyl_tensor = gkyl_malloc(sizeof(double***[4]));

        for (int i = 0; i < 4; i++) {
          spacetime_riemann_tensor[i] = gkyl_malloc(sizeof(double**[4]));
          covariant_spacetime_riemann_tensor[i] = gkyl_malloc(sizeof(double**[4]));
          spacetime_weyl_tensor[i] = gkyl_malloc(sizeof(double**[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_riemann_tensor[i][j] = gkyl_malloc(sizeof(double*[4]));
            covariant_spacetime_riemann_tensor[i][j] = gkyl_malloc(sizeof(double*[4]));
            spacetime_weyl_tensor[i][j] = gkyl_malloc(sizeof(double*[4]));

            for (int k = 0; k < 4; k++) {
              spacetime_riemann_tensor[i][j][k] = gkyl_malloc(sizeof(double[4]));
              covariant_spacetime_riemann_tensor[i][j][k] = gkyl_malloc(sizeof(double[4]));
              spacetime_weyl_tensor[i][j][k] = gkyl_malloc(sizeof(double[4]));

              for (int l = 0; l < 4; l++) {
                covariant_spacetime_riemann_tensor[i][j][k][l] = 0.0;
              }
            }
          }
        }

        spacetime->spacetime_riemann_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_riemann_tensor);
        spacetime->spacetime_weyl_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_weyl_tensor);

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
              for (int l = 0; l < 4; l++) {
                for (int m = 0; m < 4; m++) {
                  covariant_spacetime_riemann_tensor[i][j][k][l] += spacetime_metric[i][m] * spacetime_riemann_tensor[m][j][k][l];
                }
              }
            }
          }
        }

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
              for (int l = 0; l < 4; l++) {
                TEST_CHECK( gkyl_compare(covariant_spacetime_riemann_tensor[i][j][k][l], spacetime_weyl_tensor[i][j][k][l], 1e-2) );
              }
            }
          }
        }

        bool in_excision_region;
        spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);

        TEST_CHECK( (in_excision_region == false) );

        for (int i = 0; i < 3; i++) {
          gkyl_free(spatial_metric[i]);
          gkyl_free(inv_spatial_metric[i]);
          gkyl_free(spatial_metric_prod[i]);
          gkyl_free(shift_vector_der[i]);
          gkyl_free(extrinsic_curvature[i]);
          gkyl_free(shift_vector_cov_der[i]);
          gkyl_free(shift_covector_cov_der[i]);
          
          for (int j = 0; j < 3; j++) {
            gkyl_free(spatial_metric_der[i][j]);
            gkyl_free(spatial_christoffel[i][j]);
            gkyl_free(spatial_metric_cov_der[i][j]);
          }
          gkyl_free(spatial_metric_der[i]);
          gkyl_free(spatial_christoffel[i]);
          gkyl_free(spatial_metric_cov_der[i]);
        }
        gkyl_free(spatial_metric);
        gkyl_free(inv_spatial_metric);
        gkyl_free(spatial_metric_prod);
        gkyl_free(spatial_metric_der);
        gkyl_free(spatial_christoffel);
        gkyl_free(spatial_metric_cov_der);
        gkyl_free(shift_vector_der);
        gkyl_free(extrinsic_curvature);
        gkyl_free(shift_vector_cov_der);
        gkyl_free(shift_covector_cov_der);
        gkyl_free(shift_vector);

        for (int i = 0; i < 4; i++) {
          gkyl_free(spacetime_metric[i]);
          gkyl_free(inv_spacetime_metric[i]);
          gkyl_free(spacetime_metric_prod[i]);

          for (int j = 0; j < 4; j++) {
            gkyl_free(spacetime_metric_der[i][j]);
            gkyl_free(spacetime_christoffel[i][j]);
            gkyl_free(spacetime_metric_cov_der[i][j]);

            for (int k = 0; k < 4; k++) {
              gkyl_free(spacetime_riemann_tensor[i][j][k]);
              gkyl_free(covariant_spacetime_riemann_tensor[i][j][k]);
              gkyl_free(spacetime_weyl_tensor[i][j][k]);
            }
            gkyl_free(spacetime_riemann_tensor[i][j]);
            gkyl_free(covariant_spacetime_riemann_tensor[i][j]);
            gkyl_free(spacetime_weyl_tensor[i][j]);
          }
          gkyl_free(spacetime_metric_der[i]);
          gkyl_free(spacetime_christoffel[i]);
          gkyl_free(spacetime_metric_cov_der[i]);
          gkyl_free(spacetime_ricci_tensor[i]);
          gkyl_free(spacetime_riemann_tensor[i]);
          gkyl_free(covariant_spacetime_riemann_tensor[i]);
          gkyl_free(spacetime_weyl_tensor[i]);
        }
        gkyl_free(spacetime_metric);
        gkyl_free(inv_spacetime_metric);
        gkyl_free(spacetime_metric_prod);
        gkyl_free(spacetime_metric_der);
        gkyl_free(spacetime_christoffel);
        gkyl_free(spacetime_metric_cov_der);
        gkyl_free(spacetime_ricci_tensor);
        gkyl_free(spacetime_riemann_tensor);
        gkyl_free(covariant_spacetime_riemann_tensor);
        gkyl_free(spacetime_weyl_tensor);
      }
      else {
        bool in_excision_region;
        spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);

        TEST_CHECK( (in_excision_region == true) );
      }
    }
  }

  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_tov_static()
{
  double gas_gamma = 2.0;
  double K_poly = 100.0; 
  double rho_c = 5e-4;
  double dr_tov = 1e-6;//0.001;
  struct gkyl_tov *tov = gkyl_tov_new(K_poly, gas_gamma, rho_c, dr_tov);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_tov_spacetime_new(false, tov, 0.0, 0.0, 0.0);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind + 0.05;
      double y = 0.1 * y_ind;
      double z = 0.5;

      if (sqrt((x * x) + (y * y) + (z * z)) > 0.2) {
        double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
        double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
        double **spatial_metric_prod = gkyl_malloc(sizeof(double*[3]));

        for (int i = 0; i < 3; i++) {
          spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
          inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
          spatial_metric_prod[i] = gkyl_malloc(sizeof(double[3]));

          for (int j = 0; j < 3; j++) {
            spatial_metric_prod[i][j] = 0.0;
          }
        }

        double **spacetime_metric = gkyl_malloc(sizeof(double*[4]));
        double **inv_spacetime_metric = gkyl_malloc(sizeof(double*[4]));
        double **spacetime_metric_prod = gkyl_malloc(sizeof(double*[4]));

        for (int i = 0; i < 4; i++) {
          spacetime_metric[i] = gkyl_malloc(sizeof(double[4]));
          inv_spacetime_metric[i] = gkyl_malloc(sizeof(double[4]));
          spacetime_metric_prod[i] = gkyl_malloc(sizeof(double[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_metric_prod[i][j] = 0.0;
          }
        }

        spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, z, &spatial_metric);
        spacetime->spacetime_metric_tensor_func(spacetime, 0.0, x, y, z, &spacetime_metric);

        spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x, y, z, &inv_spatial_metric);
        spacetime->spacetime_inv_metric_tensor_func(spacetime, 0.0, x, y, z, &inv_spacetime_metric);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              spatial_metric_prod[i][j] += spatial_metric[i][k] * inv_spatial_metric[k][j];
            }

            if (i == j) {
              TEST_CHECK( gkyl_compare(spatial_metric_prod[i][j], 1.0, 1e-10) );
            }
            else {
              TEST_CHECK( gkyl_compare(spatial_metric_prod[i][j], 0.0, 1e-10) );
            }
          }
        }

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
              spacetime_metric_prod[i][j] += spacetime_metric[i][k] * inv_spacetime_metric[k][j];
            }

            if (i == j) {
              TEST_CHECK( gkyl_compare(spacetime_metric_prod[i][j], 1.0, 1e-10) );
            }
            else {
              TEST_CHECK( gkyl_compare(spacetime_metric_prod[i][j], 0.0, 1e-10) );
            }
          }
        }

        double spatial_metric_det;
        double spacetime_metric_det;
        double lapse_function;

        spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, z, &spatial_metric_det);
        spacetime->spacetime_metric_det_func(spacetime, 0.0, x, y, z, &spacetime_metric_det);
        spacetime->lapse_function_func(spacetime, 0.0, x, y, z, &lapse_function);

        TEST_CHECK( gkyl_compare(sqrt(-spacetime_metric_det), lapse_function * sqrt(spatial_metric_det), 1e-10) );

        double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
        double ***spatial_christoffel = gkyl_malloc(sizeof(double**[3]));
        double ***spatial_metric_cov_der = gkyl_malloc(sizeof(double**[3]));

        for (int i = 0; i < 3; i++) {
          spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));
          spatial_christoffel[i] = gkyl_malloc(sizeof(double*[3]));
          spatial_metric_cov_der[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
            spatial_christoffel[i][j] = gkyl_malloc(sizeof(double[3]));
            spatial_metric_cov_der[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, z, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spatial_metric_der);
        spacetime->spatial_christoffel_func(spacetime, 0.0, x, y, z, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spatial_christoffel);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              spatial_metric_cov_der[i][j][k] = spatial_metric_der[i][j][k];

              for (int l = 0; l < 3; l++) {
                spatial_metric_cov_der[i][j][k] -= spatial_christoffel[l][i][j] * spatial_metric[l][k];
                spatial_metric_cov_der[i][j][k] -= spatial_christoffel[l][i][k] * spatial_metric[j][l];
              }

              TEST_CHECK( gkyl_compare(spatial_metric_cov_der[i][j][k], 0.0, 1e-6) );
            }
          }
        }

        double ***spacetime_metric_der = gkyl_malloc(sizeof(double**[4]));
        double ***spacetime_christoffel = gkyl_malloc(sizeof(double**[4]));
        double ***spacetime_metric_cov_der = gkyl_malloc(sizeof(double**[4]));

        for (int i = 0; i < 4; i++) {
          spacetime_metric_der[i] = gkyl_malloc(sizeof(double*[4]));
          spacetime_christoffel[i] = gkyl_malloc(sizeof(double*[4]));
          spacetime_metric_cov_der[i] = gkyl_malloc(sizeof(double*[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_metric_der[i][j] = gkyl_malloc(sizeof(double[4]));
            spacetime_christoffel[i][j] = gkyl_malloc(sizeof(double[4]));
            spacetime_metric_cov_der[i][j] = gkyl_malloc(sizeof(double[4]));
          }
        }

        spacetime->spacetime_metric_tensor_der_func(spacetime, 0.0, x, y, z, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_metric_der);
        spacetime->spacetime_christoffel_func(spacetime, 0.0, x, y, z, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_christoffel);

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
              spacetime_metric_cov_der[i][j][k] = spacetime_metric_der[i][j][k];

              for (int l = 0; l < 4; l++) {
                spacetime_metric_cov_der[i][j][k] -= spacetime_christoffel[l][i][j] * spacetime_metric[l][k];
                spacetime_metric_cov_der[i][j][k] -= spacetime_christoffel[l][i][k] * spacetime_metric[j][l];
              }

              TEST_CHECK( gkyl_compare(spacetime_metric_cov_der[i][j][k], 0.0, 1e-6) );
            }
          }
        }

        double **shift_vector_der = gkyl_malloc(sizeof(double*[3]));
        double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
        double **shift_vector_cov_der = gkyl_malloc(sizeof(double*[3]));
        double **shift_covector_cov_der = gkyl_malloc(sizeof(double*[3]));

        for (int i = 0; i < 3; i++) {
          shift_vector_der[i] = gkyl_malloc(sizeof(double[3]));
          extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
          shift_vector_cov_der[i] = gkyl_malloc(sizeof(double[3]));
          shift_covector_cov_der[i] = gkyl_malloc(sizeof(double[3]));

          for (int j = 0; j < 3; j++){
            shift_covector_cov_der[i][j] = 0.0;
          }
        }

        spacetime->shift_vector_der_func(spacetime, 0.0, x, y, z, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &shift_vector_der);
        spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, z, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &extrinsic_curvature);

        // // extrinsic_curvature must be symmetric
        // for (int i = 0; i < 3; i++) {
        //     for (int j = 0; j < 3; j++) {
        //         TEST_CHECK( fabs(extrinsic_curvature[i][j] - extrinsic_curvature[j][i]) < 1e-10 );
        //     }
        // }

        // // In vacuum exterior: R + K²- K_ij K^ij = 0 (Hamiltonian constraint, ρ = 0)
        // double extrinsic_curvature_trace = 0.0;
        // double extrinsic_curvature_sq = 0.0;  // K_ij K^ij

        // for (int i = 0; i < 3; i++) {
        //   for (int j = 0; j < 3; j++) {
        //       extrinsic_curvature_trace += inv_spatial_metric[i][j] * extrinsic_curvature[i][j];
        //   }
        // }

        // for (int i = 0; i < 3; i++) {
        //     for (int j = 0; j < 3; j++) {
        //         for (int k = 0; k < 3; k++) {
        //             for (int l = 0; l < 3; l++) {
        //                 extrinsic_curvature_sq += inv_spatial_metric[i][k] * inv_spatial_metric[j][l] * extrinsic_curvature[i][j] * extrinsic_curvature[k][l];
        //             }
        //         }
        //     }
        // }
        // double spatial_ricci_scalar;
        // spacetime->spatial_ricci_scalar_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spatial_ricci_scalar);
        // // Hamiltonian constraint: R + K² - K_ij K^ij = 0
        // TEST_CHECK(fabs(spatial_ricci_scalar + (extrinsic_curvature_trace * extrinsic_curvature_trace) - extrinsic_curvature_sq) < 1e-8);


        double *shift_vector = gkyl_malloc(sizeof(double[3]));
        spacetime->shift_vector_func(spacetime, 0.0, x, y, z, &shift_vector);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            shift_vector_cov_der[i][j] = shift_vector_der[i][j];

            for (int k = 0; k < 3; k++) {
              shift_vector_cov_der[i][j] += spatial_christoffel[j][i][k] * shift_vector[k];
            }
          }
        }

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              shift_covector_cov_der[i][j] += spatial_metric[j][k] * shift_vector_cov_der[i][k];
            }
          }
        }

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            TEST_CHECK( gkyl_compare(2.0 * lapse_function * extrinsic_curvature[i][j], -(shift_covector_cov_der[j][i] + shift_covector_cov_der[i][j]), 1e-6) );
          }
        }

        // Coordinate transformation tests
        double r = sqrt(x * x + y * y + z * z);

        // Lapse must be positive
        TEST_CHECK(lapse_function > 0.0);

        // Spatial metric must be symmetric
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            TEST_CHECK(fabs(spatial_metric[i][j] - spatial_metric[j][i]) < 1e-14);
          }
        }

        // Spatial metric must be positive definite
        for (int i = 0; i < 3; i++) {
          TEST_CHECK(spatial_metric[i][i] > 0.0);
        }
        
        for (int i = 0; i < 3; i++)
        TEST_CHECK(spatial_metric[i][i] > 0.0);
        // Check determinant positive
        double det =
            spatial_metric[0][0] * (spatial_metric[1][1]*spatial_metric[2][2] - spatial_metric[1][2]*spatial_metric[2][1])
          - spatial_metric[0][1] * (spatial_metric[1][0]*spatial_metric[2][2] - spatial_metric[1][2]*spatial_metric[2][0])
          + spatial_metric[0][2] * (spatial_metric[1][0]*spatial_metric[2][1] - spatial_metric[1][1]*spatial_metric[2][0]);
        TEST_CHECK(det > 0.0);

        // positive spacetime det
        TEST_CHECK(gkyl_compare(spacetime_metric_det, -(lapse_function * lapse_function) * det, 1e-10));

        // g_00 = -alpha^2 + gamma_ij beta^i beta^j
        double shift_vector_sq = 0.0;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                shift_vector_sq += spatial_metric[i][j] * shift_vector[i] * shift_vector[j];
            }
        }
        TEST_CHECK(gkyl_compare(spacetime_metric[0][0], -(lapse_function * lapse_function) + shift_vector_sq, 1e-10));

        // Curvature tests

        struct tov_eval_fluid fluid;
        gkyl_tov_eval(tov, x, y, z, &fluid);
        double rho = fluid.rho;
        double p = fluid.P;

        double rho_atm = 1e-2 * rho_c;
        double p_atm = K_poly * pow(rho_atm, gas_gamma);

        double R_star = gkyl_tov_star_radius(tov);

        if (r < R_star) {
          if (rho > rho_atm && p > p_atm) { 

          double h = 1.0 + (gas_gamma / (gas_gamma - 1.0)) * (p / rho);
          
          double **spacetime_ricci_tensor = gkyl_malloc(sizeof(double*[4]));
          for (int i = 0; i < 4; i++) {
            spacetime_ricci_tensor[i] = gkyl_malloc(sizeof(double[4]));
          }

          spacetime->spacetime_ricci_tensor_func(spacetime, 0.0, x, y, z, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_ricci_tensor);

          //For option A:
          double u_con[4] = {1.0 / lapse_function, 0, 0, 0};
          double u_cov[4] = {-lapse_function, 0.0, 0.0, 0.0};

          //double u_cov[4] = {1.0 / sqrt(-spacetime_metric[0][0]), 0, 0, 0};
          // double u_con[4] = {
          //     1.0 / lapse_function,
          //     -shift_vector[0] / lapse_function,
          //     -shift_vector[1] / lapse_function,
          //     -shift_vector[2] / lapse_function
          // };
          // double u_con[4] = {
          //     1.0 / sqrt(-spacetime_metric[0][0]),
          //     0.0,
          //     0.0,
          //     0.0
          // };

          double norm = 0.0;
          for (int mu = 0; mu < 4; mu++) {
            for (int nu = 0; nu < 4; nu++) {
              norm += spacetime_metric[mu][nu] * u_con[mu] * u_con[nu];
            }
          }
          TEST_CHECK(fabs(norm + 1.0) < 1e-12);

          // The BL -> KS/CKS transform only changes the slicing, not the spatial
          // coordinates, so a static TOV fluid element remains coordinate-stationary.
          // double u_con[4] = {
          //     1.0 / sqrt(-spacetime_metric[0][0]),
          //     0.0,
          //     0.0,
          //     0.0
          // };

          
          for (int mu = 0; mu < 4; mu++) {
            for (int nu = 0; nu < 4; nu++) {
                u_cov[mu] += spacetime_metric[mu][nu] * u_con[nu];
            }
          }

          double stress_energy_tensor[4][4];
          for (int mu = 0; mu < 4; mu++) {
              for (int nu = 0; nu < 4; nu++) {
                  stress_energy_tensor[mu][nu] = (rho * h * u_cov[mu] * u_cov[nu]) + (p * spacetime_metric[mu][nu]);
              }
          }

          double spacetime_ricci_scalar = 0.0;
          for (int mu = 0; mu < 4; mu++) {
              for (int nu = 0; nu < 4; nu++) {
                  spacetime_ricci_scalar += inv_spacetime_metric[mu][nu] * spacetime_ricci_tensor[mu][nu];
              }
          }

          // Ricci scalar invariance
          double stress_energy_tensor_trace = 0.0;
          for (int mu = 0; mu < 4; mu++)
              for (int nu = 0; nu < 4; nu++)
                  stress_energy_tensor_trace += inv_spacetime_metric[mu][nu] * stress_energy_tensor[mu][nu];

          double einstein_tensor[4][4];
          for (int mu = 0; mu < 4; mu++) {
              for (int nu = 0; nu < 4; nu++) {
                  einstein_tensor[mu][nu] = spacetime_ricci_tensor[mu][nu] - 0.5 * spacetime_metric[mu][nu] * spacetime_ricci_scalar;
              }
          }
          
          /*
          // // Hamiltonian constraint: R_spatial + K^2 - K_ij K^ij = 16π E,
          // // where E = T_{mu nu} n^mu n^nu is the Eulerian energy density.
          // // K = gamma^ij K_ij
          // double extrinsic_curvature_trace = 0.0;
          // for (int i = 0; i < 3; i++)
          //     for (int j = 0; j < 3; j++)
          //         extrinsic_curvature_trace += inv_spatial_metric[i][j] * extrinsic_curvature[i][j];

          // // K_ij K^ij = gamma^ik gamma^jl K_ij K_kl
          // double extrinsic_curvature_sq = 0.0;
          // for (int i = 0; i < 3; i++)
          //     for (int j = 0; j < 3; j++)
          //         for (int k = 0; k < 3; k++)
          //             for (int l = 0; l < 3; l++)
          //                 extrinsic_curvature_sq += inv_spatial_metric[i][k] * inv_spatial_metric[j][l] * extrinsic_curvature[i][j] * extrinsic_curvature[k][l];

          // double spatial_ricci_scalar; 
          // spacetime->spatial_ricci_scalar_func(spacetime, 0.0, x, y, z,
          //     pow(10.0, -2.0), pow(10.0, -2.0), pow(10.0, -2.0), &spatial_ricci_scalar);
          //     //needs to be e-2 for the Hamiltonian constraint

          // double n_con[4] = {
          //     1.0 / lapse_function,
          //     -shift_vector[0] / lapse_function,
          //     -shift_vector[1] / lapse_function,
          //     -shift_vector[2] / lapse_function
          // };

          // double eulerian_energy_density = 0.0;
          // for (int mu = 0; mu < 4; mu++) {
          //     for (int nu = 0; nu < 4; nu++) {
          //         eulerian_energy_density += stress_energy_tensor[mu][nu] * n_con[mu] * n_con[nu];
          //     }
          // }

          // double lhs = spatial_ricci_scalar
          //     + extrinsic_curvature_trace * extrinsic_curvature_trace 
          //     - extrinsic_curvature_sq;
          // double lhs_alt = -spatial_ricci_scalar
          //     + extrinsic_curvature_trace * extrinsic_curvature_trace
          //     - extrinsic_curvature_sq;
          // double rhs = 16.0 * M_PI * eulerian_energy_density;

          // static bool printed_first_hamiltonian_fail = false;
          // if (fabs(lhs - rhs) > 1e-2 && !printed_first_hamiltonian_fail) {
          //     printf("Hamiltonian constraint FAIL at x=%g y=%g z=%g r=%g\n"
          //           "  R_spatial=%g  K_trace=%g  K_sq=%g\n"
          //           "  lhs=%g  rhs=%g  diff=%g\n",
          //           x, y, z, sqrt(x*x+y*y+z*z),
          //           spatial_ricci_scalar, extrinsic_curvature_trace,
          //           extrinsic_curvature_sq,
          //           lhs, rhs, lhs-rhs);
          //     printf("  alt lhs (-R + K^2 - K_ij K^ij) = %g  alt diff = %g\n",
          //       lhs_alt, lhs_alt - rhs);
          //     printf("At r=%g: rho=%g e=%g p=%g\n", sqrt(x*x+y*y+z*z), rho, fluid.e, p);
          //     printf("rhs = 16pi*E = %g\n", rhs);

          //     printf("shift = %g %g %g\n", shift_vector[0], shift_vector[1], shift_vector[2]);
          //     printf("lapse = %g\n", lapse_function);
          //     printf("K_trace=%g \n", extrinsic_curvature_trace);
          //     printf("Spacetime Ricci scalar = %.12g\n", spacetime_ricci_scalar);
          //     printf("-8pi*T_trace = %.12g\n", -8.0 * M_PI * stress_energy_tensor_trace);
          //     printf("Ricci scalar diff = %.12g\n",
          //       spacetime_ricci_scalar - (-8.0 * M_PI * stress_energy_tensor_trace));
          //     printf("Alt Ricci scalar diff (-R - (-8pi*T_trace)) = %.12g\n",
          //       (-spacetime_ricci_scalar) - (-8.0 * M_PI * stress_energy_tensor_trace));

          //     printf("Einstein equation residuals at failing point:\n");
          //     double max_einstein_diff = 0.0;
          //     double max_alt_einstein_diff = 0.0;
          //     for (int mu = 0; mu < 4; mu++) {
          //       for (int nu = 0; nu < 4; nu++) {
          //         double einstein_rhs = 8.0 * M_PI * stress_energy_tensor[mu][nu];
          //         double einstein_diff = einstein_tensor[mu][nu] - einstein_rhs;
          //         double alt_einstein_diff = -einstein_tensor[mu][nu] - einstein_rhs;
          //         if (fabs(einstein_diff) > max_einstein_diff)
          //           max_einstein_diff = fabs(einstein_diff);
          //         if (fabs(alt_einstein_diff) > max_alt_einstein_diff)
          //           max_alt_einstein_diff = fabs(alt_einstein_diff);
          //         printf("  G[%d][%d] = %.12g  8piT[%d][%d] = %.12g  diff = %.12g\n",
          //           mu, nu, einstein_tensor[mu][nu], mu, nu, einstein_rhs, einstein_diff);
          //       }
          //     }
          //     printf("Max |G - 8piT| = %.12g\n", max_einstein_diff);
          //     printf("Max |-G - 8piT| = %.12g\n", max_alt_einstein_diff);

          //     {
          //       double h_vals[3] = {1e-2, 1e-3, 1e-4};
          //       printf("Hamiltonian convergence :\n");
          //       for (int h_idx = 0; h_idx < 3; h_idx++) {
          //         double h_fd = h_vals[h_idx];
          //         double spatial_ricci_h = 0.0;
          //         double **extrinsic_curvature_h = gkyl_malloc(sizeof(double*[3]));
          //         for (int i = 0; i < 3; i++) {
          //           extrinsic_curvature_h[i] = gkyl_malloc(sizeof(double[3]));
          //         }

          //         spacetime->spatial_ricci_scalar_func(spacetime, 0.0, x, y, z,
          //           h_fd, h_fd, h_fd, &spatial_ricci_h);
          //         spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, z,
          //           h_fd, h_fd, h_fd, &extrinsic_curvature_h);

          //         double K_trace_h = 0.0;
          //         double K_sq_h = 0.0;
          //         for (int i = 0; i < 3; i++) {
          //           for (int j = 0; j < 3; j++) {
          //             K_trace_h += inv_spatial_metric[i][j] * extrinsic_curvature_h[i][j];
          //             for (int k = 0; k < 3; k++) {
          //               for (int l = 0; l < 3; l++) {
          //                 K_sq_h += inv_spatial_metric[i][k] * inv_spatial_metric[j][l]
          //                   * extrinsic_curvature_h[i][j] * extrinsic_curvature_h[k][l];
          //               }
          //             }
          //           }
          //         }

          //         double lhs_h = spatial_ricci_h + K_trace_h * K_trace_h - K_sq_h;
          //         double lhs_alt_h = -spatial_ricci_h + K_trace_h * K_trace_h - K_sq_h;
          //         printf("  h=%g: R=% .12g lhs=% .12g diff=% .12g alt_lhs=% .12g alt_diff=% .12g\n",
          //           h_fd, spatial_ricci_h, lhs_h, lhs_h - rhs, lhs_alt_h, lhs_alt_h - rhs);

          //         for (int i = 0; i < 3; i++) {
          //           gkyl_free(extrinsic_curvature_h[i]);
          //         }
          //         gkyl_free(extrinsic_curvature_h);
          //       }
          //     }
          //     printed_first_hamiltonian_fail = true;
          //     }

          // TEST_CHECK(fabs(lhs - rhs) < 1e-2);

          // //K is good. R_spatial is the problem. 
          */


          for (int mu = 0; mu < 4; mu++) {
              for (int nu = 0; nu < 4; nu++) {
                double curvature_test = fabs(einstein_tensor[mu][nu] - 8.0 * M_PI * stress_energy_tensor[mu][nu]);
                TEST_CHECK(curvature_test < 1e-6);
                if (curvature_test > 1e-6) printf("%e, %e, \n", curvature_test, r);

          // // Time components - well determined - FAIL at 1e-2
          // for (int mu = 0; mu < 4; mu++)
          //     TEST_CHECK(fabs(einstein_tensor[0][mu] - 
          //         8.0*M_PI*stress_energy_tensor[0][mu]) < 1e-2);

          // // Ricci scalar invariant - FAIL at 1e-2
          // double stress_energy_tensor_trace = 0.0;
          // for (int mu = 0; mu < 4; mu++)
          //     for (int nu = 0; nu < 4; nu++)
          //         stress_energy_tensor_trace += inv_spacetime_metric[mu][nu] * stress_energy_tensor[mu][nu];
          // TEST_CHECK(fabs(spacetime_ricci_scalar - (-8.0*M_PI*stress_energy_tensor_trace)) < 1e-2);

          // printf("=== Einstein check at failing point ===\n");
          // printf("Ricci scalar = %g\n", spacetime_ricci_scalar);
          // printf("-8pi*T_trace = %g\n", -8.0*M_PI*stress_energy_tensor_trace);
          // printf("diff = %g\n", spacetime_ricci_scalar - (-8.0*M_PI*stress_energy_tensor_trace));

          // for (int mu = 0; mu < 4; mu++)
          //     for (int nu = 0; nu < 4; nu++) {
          //         double diff = einstein_tensor[mu][nu] - 8.0*M_PI*stress_energy_tensor[mu][nu];
          //         printf("FAIL G[%d][%d]=%g  8piT[%d][%d]=%g  diff=%g\n",
          //                 mu, nu, einstein_tensor[mu][nu],
          //                 mu, nu, 8.0*M_PI*stress_energy_tensor[mu][nu],
          //                 diff);
          //     }
          // printf("x=%g y=%g z=%g r=%g\n", x, y, z, sqrt(x*x+y*y+z*z));

          // //Looks like it might be on the fluid side

              }
          }

          // printf("r at test point = %g\n", sqrt(x*x + y*y));
          // printf("r at x+dx = %g\n", sqrt((x+1e-3)*(x+1e-3) + y*y));

          for (int i = 0; i < 4; i++) {
            gkyl_free(spacetime_ricci_tensor[i]);
          }
        gkyl_free(spacetime_ricci_tensor);
        }
        }

          // double lapse_bl;
          // double shift_bl[3];
          // double gamma_bl[3][3];
          // double gamma_cks_local[3][3];
          // for (int i = 0; i < 3; i++)
          //     for (int j = 0; j < 3; j++)
          //         gamma_cks_local[i][j] = spatial_metric[i][j];

          // cks_to_bl(lapse_function, shift_vector, gamma_cks_local,
          //     x, y, 0.0,
          //     &lapse_bl, shift_bl, gamma_bl);

          // // Verify against gkyl_tov_eval_bl
          // double r = sqrt(x*x + y*y);
          // struct tov_eval_bl bl;
          // gkyl_tov_eval_bl(tov, r, &bl);
          // double f = 1.0 - 2.0*bl.m / r;
          // double A_BL_expected = (f > 1e-10) ? 1.0/f : 1e10;
          // double alpha_BL_expected = exp(bl.Phi);
          // double shift_expected = 

          // printf("lapse_bl=%g  expected=%g  diff=%g\n",
          //     lapse_bl, alpha_BL_expected, lapse_bl - alpha_BL_expected);
          // double gamma_rr_bl_check = 0.0;
          // for (int i = 0; i < 3; i++)
          //     for (int j = 0; j < 3; j++)
          //         gamma_rr_bl_check += gamma_bl[i][j] * n[i] * n[j];

          // printf("gamma_rr_bl=%g  expected=%g  diff=%g\n",
          //     gamma_rr_bl_check, A_BL_expected, gamma_rr_bl_check - A_BL_expected);
          // printf("shift_bl = %g %g %g  expected = 0 0 0\n", shift_bl[0], shift_bl[1], shift_bl[2]);

          // for (int mu = 0; mu < 4; mu++) {
          //     for (int nu = 0; nu < 4; nu++) {
          //         TEST_CHECK(fabs(einstein_tensor[mu][nu] - 8.0 * M_PI * stress_energy_tensor[mu][nu]) < 1e-1);
          //     }
          // }

        // // Curvature tests with BL cooridnates

        // struct tov_eval_fluid fluid;
        // gkyl_tov_eval(tov, x, y, 0.0, &fluid);
        // double rho = fluid.rho;
        // double p = fluid.P;

        // double rho_atm = 1e-2 * rho_c;
        // double p_atm = K_poly * pow(rho_atm, gas_gamma);

        // double R_star = gkyl_tov_star_radius(tov);

        // if (r < R_star) {
        //   if (rho > rho_atm && p > p_atm) { 

        //   double h = 1.0 + (gas_gamma / (gas_gamma - 1.0)) * (p / rho);
          
        //   double **spacetime_ricci_tensor = gkyl_malloc(sizeof(double*[4]));
        //   for (int i = 0; i < 4; i++) {
        //     spacetime_ricci_tensor[i] = gkyl_malloc(sizeof(double[4]));
        //   }

        //   spacetime->spacetime_ricci_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_ricci_tensor);
          
        //   // u^μ = (1/α, -β^i/α) in CKS for fluid at rest in BL
        //   // double u_con[4] = {
        //   //     1.0 / lapse_function,
        //   //     -shift_vector[0] / lapse_function,
        //   //     -shift_vector[1] / lapse_function,
        //   //     -shift_vector[2] / lapse_function
        //   // };
        //   double u_con[4] = {
        //     1.0 / sqrt(-spacetime_metric[0][0]),
        //     0.0,
        //     0.0,
        //     0.0
        //   };

        //   double u_cov[4] = {0.0, 0.0, 0.0, 0.0};
        //   for (int mu = 0; mu < 4; mu++) {
        //     for (int nu = 0; nu < 4; nu++) {
        //         u_cov[mu] += spacetime_metric[mu][nu] * u_con[nu];
        //     }
        //   }

        //   double stress_energy_tensor[4][4];
        //   for (int mu = 0; mu < 4; mu++) {
        //       for (int nu = 0; nu < 4; nu++) {
        //           stress_energy_tensor[mu][nu] = (rho * h * u_cov[mu] * u_cov[nu]) + (p * spacetime_metric[mu][nu]);
        //       }
        //   }

        //   double spacetime_ricci_scalar = 0.0;
        //   for (int mu = 0; mu < 4; mu++) {
        //       for (int nu = 0; nu < 4; nu++) {
        //           spacetime_ricci_scalar += inv_spacetime_metric[mu][nu] * spacetime_ricci_tensor[mu][nu];
        //       }
        //   }

        //   double einstein_tensor[4][4];
        //   for (int mu = 0; mu < 4; mu++) {
        //       for (int nu = 0; nu < 4; nu++) {
        //           einstein_tensor[mu][nu] = spacetime_ricci_tensor[mu][nu] - 0.5 * spacetime_metric[mu][nu] * spacetime_ricci_scalar;
        //       }
        //   }

        //   for (int mu = 0; mu < 4; mu++) {
        //       for (int nu = 0; nu < 4; nu++) {
        //         double curvature_test = fabs(einstein_tensor[mu][nu] - 8.0 * M_PI * stress_energy_tensor[mu][nu]);
        //         TEST_CHECK(curvature_test < 47);
        //         if (curvature_test > 47) printf("%e, %e, \n", curvature_test, r);
        //         // printf("x=%g y=%g\n", x, y);
        //         // printf("rho=%g p=%g h=%g\n", rho, p, h);
        //         // printf("lapse=%g\n", lapse_function);
        //         // printf("g_00=%g\n", spacetime_metric[0][0]);
        //         // printf("u_cov = %g %g %g %g\n", u_cov[0], u_cov[1], u_cov[2], u_cov[3]);
        //         // double norm = 0.0;
        //         // for (int mu = 0; mu < 4; mu++) {
        //         //     norm += u_con[mu] * u_cov[mu];
        //         // }
        //         // printf("u_mu u^mu = %g\n", norm);  // must be -1
        //         // for (int mu = 0; mu < 4; mu++) {
        //         //     for (int nu = 0; nu < 4; nu++) {
        //         //         printf("G[%d][%d]=%g  8piT[%d][%d]=%g  diff=%g\n",
        //         //             mu, nu, einstein_tensor[mu][nu],
        //         //             mu, nu, 8.0*M_PI*stress_energy_tensor[mu][nu],
        //         //             einstein_tensor[mu][nu] - 8.0*M_PI*stress_energy_tensor[mu][nu]);
        //         //         }
        //         // }
        //       }
        //   }

        //   for (int i = 0; i < 4; i++) {
        //     gkyl_free(spacetime_ricci_tensor[i]);
        //   }
        // gkyl_free(spacetime_ricci_tensor);
        // }
        // }

          
        


        //Since the TOV has no excision region, we test several points in the domain
        //For K=100, gamma=2, rho_c=5e-4 the TOV solution gives R_star ~ 9.5
        double test_points[][3] = {
          {0.0, 0.0, 0.0}, // center
          {R_star - 2.0, 0.0, 0.0}, // interior (~0.5 * R_star) R_star
          {R_star,  0.0, 0.0}, // surface  (~R_star)
          {R_star + 2.0, 0.0, 0.0}, // exterior (~2 * R_star)
        };

        for (int p = 0; p < 4; p++) {
          bool in_excision_region;
          spacetime->excision_region_func(spacetime, 0.0, test_points[p][0], test_points[p][1], test_points[p][2], &in_excision_region);
          
          TEST_CHECK(in_excision_region == false);
        }

        for (int i = 0; i < 3; i++) {
          gkyl_free(spatial_metric[i]);
          gkyl_free(inv_spatial_metric[i]);
          gkyl_free(spatial_metric_prod[i]);
          gkyl_free(shift_vector_der[i]);
          gkyl_free(extrinsic_curvature[i]);
          gkyl_free(shift_vector_cov_der[i]);
          gkyl_free(shift_covector_cov_der[i]);
          
          for (int j = 0; j < 3; j++) {
            gkyl_free(spatial_metric_der[i][j]);
            gkyl_free(spatial_christoffel[i][j]);
            gkyl_free(spatial_metric_cov_der[i][j]);
          }
          gkyl_free(spatial_metric_der[i]);
          gkyl_free(spatial_christoffel[i]);
          gkyl_free(spatial_metric_cov_der[i]);
        }
        gkyl_free(spatial_metric);
        gkyl_free(inv_spatial_metric);
        gkyl_free(spatial_metric_prod);
        gkyl_free(spatial_metric_der);
        gkyl_free(spatial_christoffel);
        gkyl_free(spatial_metric_cov_der);
        gkyl_free(shift_vector_der);
        gkyl_free(extrinsic_curvature);
        gkyl_free(shift_vector_cov_der);
        gkyl_free(shift_covector_cov_der);
        gkyl_free(shift_vector);

        for (int i = 0; i < 4; i++) {
          gkyl_free(spacetime_metric[i]);
          gkyl_free(inv_spacetime_metric[i]);
          gkyl_free(spacetime_metric_prod[i]);

          for (int j = 0; j < 4; j++) {
            gkyl_free(spacetime_metric_der[i][j]);
            gkyl_free(spacetime_christoffel[i][j]);
            gkyl_free(spacetime_metric_cov_der[i][j]);
          }
          gkyl_free(spacetime_metric_der[i]);
          gkyl_free(spacetime_christoffel[i]);
          gkyl_free(spacetime_metric_cov_der[i]);
        }
        gkyl_free(spacetime_metric);
        gkyl_free(inv_spacetime_metric);
        gkyl_free(spacetime_metric_prod);
        gkyl_free(spacetime_metric_der);
        gkyl_free(spacetime_christoffel);
        gkyl_free(spacetime_metric_cov_der);

      }

    }
  }  
  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_neutronstar_static()
{
  double mass = 0.1;
  double spin = 0.0;
  double ang_mom = spin * mass;

  double spin_dimensionless = ang_mom / (mass * mass);

  double alpha = 5.0;
  double beta = pow(-0.36 + (1.48 * pow(sqrt(alpha), 0.65)), 3.0);
  double gamma = pow(-4.749 + (0.27613 * pow(sqrt(alpha), 1.5146)) + (5.5168 * pow(sqrt(alpha), 0.22229)), 4.0);

  double mass_quadrupole = -alpha * (spin_dimensionless * spin_dimensionless) * (mass * mass * mass);
  double spin_octupole = -beta * (spin_dimensionless * spin_dimensionless * spin_dimensionless) * (mass * mass * mass * mass);
  double mass_hexadecapole = gamma * (spin_dimensionless * spin_dimensionless * spin_dimensionless * spin_dimensionless) * (mass * mass * mass * mass * mass);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_neutronstar_new(false, mass, spin, mass_quadrupole, spin_octupole, mass_hexadecapole, 0.0, 0.0, 0.0);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      if (sqrt((x * x) + (y * y)) > 0.2) {
        double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
        double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
        double **spatial_metric_prod = gkyl_malloc(sizeof(double*[3]));

        for (int i = 0; i < 3; i++) {
          spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
          inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
          spatial_metric_prod[i] = gkyl_malloc(sizeof(double[3]));

          for (int j = 0; j < 3; j++) {
            spatial_metric_prod[i][j] = 0.0;
          }
        }

        double **spacetime_metric = gkyl_malloc(sizeof(double*[4]));
        double **inv_spacetime_metric = gkyl_malloc(sizeof(double*[4]));
        double **spacetime_metric_prod = gkyl_malloc(sizeof(double*[4]));

        for (int i = 0; i < 4; i++) {
          spacetime_metric[i] = gkyl_malloc(sizeof(double[4]));
          inv_spacetime_metric[i] = gkyl_malloc(sizeof(double[4]));
          spacetime_metric_prod[i] = gkyl_malloc(sizeof(double[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_metric_prod[i][j] = 0.0;
          }
        }

        spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);
        spacetime->spacetime_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spacetime_metric);

        spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &inv_spatial_metric);
        spacetime->spacetime_inv_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &inv_spacetime_metric);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              spatial_metric_prod[i][j] += spatial_metric[i][k] * inv_spatial_metric[k][j];
            }

            if (i == j) {
              TEST_CHECK( gkyl_compare(spatial_metric_prod[i][j], 1.0, 1e-10) );
            }
            else {
              TEST_CHECK( gkyl_compare(spatial_metric_prod[i][j], 0.0, 1e-10) );
            }
          }
        }

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
              spacetime_metric_prod[i][j] += spacetime_metric[i][k] * inv_spacetime_metric[k][j];
            }

            if (i == j) {
              TEST_CHECK( gkyl_compare(spacetime_metric_prod[i][j], 1.0, 1e-10) );
            }
            else {
              TEST_CHECK( gkyl_compare(spacetime_metric_prod[i][j], 0.0, 1e-10) );
            }
          }
        }

        double spatial_metric_det;
        double spacetime_metric_det;
        double lapse_function;

        spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &spatial_metric_det);
        spacetime->spacetime_metric_det_func(spacetime, 0.0, x, y, 0.0, &spacetime_metric_det);
        spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse_function);

        TEST_CHECK( gkyl_compare(sqrt(-spacetime_metric_det), lapse_function * sqrt(spatial_metric_det), 1e-10) );

        double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
        double ***spatial_christoffel = gkyl_malloc(sizeof(double**[3]));
        double ***spatial_metric_cov_der = gkyl_malloc(sizeof(double**[3]));

        for (int i = 0; i < 3; i++) {
          spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));
          spatial_christoffel[i] = gkyl_malloc(sizeof(double*[3]));
          spatial_metric_cov_der[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
            spatial_christoffel[i][j] = gkyl_malloc(sizeof(double[3]));
            spatial_metric_cov_der[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spatial_metric_der);
        spacetime->spatial_christoffel_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spatial_christoffel);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              spatial_metric_cov_der[i][j][k] = spatial_metric_der[i][j][k];

              for (int l = 0; l < 3; l++) {
                spatial_metric_cov_der[i][j][k] -= spatial_christoffel[l][i][j] * spatial_metric[l][k];
                spatial_metric_cov_der[i][j][k] -= spatial_christoffel[l][i][k] * spatial_metric[j][l];
              }

              TEST_CHECK( gkyl_compare(spatial_metric_cov_der[i][j][k], 0.0, 1e-6) );
            }
          }
        }

        double ***spacetime_metric_der = gkyl_malloc(sizeof(double**[4]));
        double ***spacetime_christoffel = gkyl_malloc(sizeof(double**[4]));
        double ***spacetime_metric_cov_der = gkyl_malloc(sizeof(double**[4]));

        for (int i = 0; i < 4; i++) {
          spacetime_metric_der[i] = gkyl_malloc(sizeof(double*[4]));
          spacetime_christoffel[i] = gkyl_malloc(sizeof(double*[4]));
          spacetime_metric_cov_der[i] = gkyl_malloc(sizeof(double*[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_metric_der[i][j] = gkyl_malloc(sizeof(double[4]));
            spacetime_christoffel[i][j] = gkyl_malloc(sizeof(double[4]));
            spacetime_metric_cov_der[i][j] = gkyl_malloc(sizeof(double[4]));
          }
        }

        spacetime->spacetime_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_metric_der);
        spacetime->spacetime_christoffel_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_christoffel);

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
              spacetime_metric_cov_der[i][j][k] = spacetime_metric_der[i][j][k];

              for (int l = 0; l < 4; l++) {
                spacetime_metric_cov_der[i][j][k] -= spacetime_christoffel[l][i][j] * spacetime_metric[l][k];
                spacetime_metric_cov_der[i][j][k] -= spacetime_christoffel[l][i][k] * spacetime_metric[j][l];
              }

              TEST_CHECK( gkyl_compare(spacetime_metric_cov_der[i][j][k], 0.0, 1e-6) );
            }
          }
        }

        double **shift_vector_der = gkyl_malloc(sizeof(double*[3]));
        double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
        double **shift_vector_cov_der = gkyl_malloc(sizeof(double*[3]));
        double **shift_covector_cov_der = gkyl_malloc(sizeof(double*[3]));

        for (int i = 0; i < 3; i++) {
          shift_vector_der[i] = gkyl_malloc(sizeof(double[3]));
          extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
          shift_vector_cov_der[i] = gkyl_malloc(sizeof(double[3]));
          shift_covector_cov_der[i] = gkyl_malloc(sizeof(double[3]));

          for (int j = 0; j < 3; j++){
            shift_covector_cov_der[i][j] = 0.0;
          }
        }

        spacetime->shift_vector_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &shift_vector_der);
        spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &extrinsic_curvature);

        double *shift_vector = gkyl_malloc(sizeof(double[3]));
        spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift_vector);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            shift_vector_cov_der[i][j] = shift_vector_der[i][j];

            for (int k = 0; k < 3; k++) {
              shift_vector_cov_der[i][j] += spatial_christoffel[j][i][k] * shift_vector[k];
            }
          }
        }

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              shift_covector_cov_der[i][j] += spatial_metric[j][k] * shift_vector_cov_der[i][k];
            }
          }
        }

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            TEST_CHECK( gkyl_compare(2.0 * lapse_function * extrinsic_curvature[i][j], -(shift_covector_cov_der[j][i] + shift_covector_cov_der[i][j]), 1e-6) );
          }
        }

        bool in_excision_region;
        spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);

        TEST_CHECK( (in_excision_region == false) );

        for (int i = 0; i < 3; i++) {
          gkyl_free(spatial_metric[i]);
          gkyl_free(inv_spatial_metric[i]);
          gkyl_free(spatial_metric_prod[i]);
          gkyl_free(shift_vector_der[i]);
          gkyl_free(extrinsic_curvature[i]);
          gkyl_free(shift_vector_cov_der[i]);
          gkyl_free(shift_covector_cov_der[i]);
          
          for (int j = 0; j < 3; j++) {
            gkyl_free(spatial_metric_der[i][j]);
            gkyl_free(spatial_christoffel[i][j]);
            gkyl_free(spatial_metric_cov_der[i][j]);
          }
          gkyl_free(spatial_metric_der[i]);
          gkyl_free(spatial_christoffel[i]);
          gkyl_free(spatial_metric_cov_der[i]);
        }
        gkyl_free(spatial_metric);
        gkyl_free(inv_spatial_metric);
        gkyl_free(spatial_metric_prod);
        gkyl_free(spatial_metric_der);
        gkyl_free(spatial_christoffel);
        gkyl_free(spatial_metric_cov_der);
        gkyl_free(shift_vector_der);
        gkyl_free(extrinsic_curvature);
        gkyl_free(shift_vector_cov_der);
        gkyl_free(shift_covector_cov_der);
        gkyl_free(shift_vector);

        for (int i = 0; i < 4; i++) {
          gkyl_free(spacetime_metric[i]);
          gkyl_free(inv_spacetime_metric[i]);
          gkyl_free(spacetime_metric_prod[i]);

          for (int j = 0; j < 4; j++) {
            gkyl_free(spacetime_metric_der[i][j]);
            gkyl_free(spacetime_christoffel[i][j]);
            gkyl_free(spacetime_metric_cov_der[i][j]);
          }
          gkyl_free(spacetime_metric_der[i]);
          gkyl_free(spacetime_christoffel[i]);
          gkyl_free(spacetime_metric_cov_der[i]);
        }
        gkyl_free(spacetime_metric);
        gkyl_free(inv_spacetime_metric);
        gkyl_free(spacetime_metric_prod);
        gkyl_free(spacetime_metric_der);
        gkyl_free(spacetime_christoffel);
        gkyl_free(spacetime_metric_cov_der);
      }
      else {
        bool in_excision_region;
        spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);

        TEST_CHECK( (in_excision_region == true) );
      }
    }
  }

  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_neutronstar_spinning()
{
  double mass = 0.1;
  double spin = -0.12;
  double ang_mom = spin * mass;

  double spin_dimensionless = ang_mom / (mass * mass);

  double alpha = 5.0;
  double beta = pow(-0.36 + (1.48 * pow(sqrt(alpha), 0.65)), 3.0);
  double gamma = pow(-4.749 + (0.27613 * pow(sqrt(alpha), 1.5146)) + (5.5168 * pow(sqrt(alpha), 0.22229)), 4.0);

  double mass_quadrupole = -alpha * (spin_dimensionless * spin_dimensionless) * (mass * mass * mass);
  double spin_octupole = -beta * (spin_dimensionless * spin_dimensionless * spin_dimensionless) * (mass * mass * mass * mass);
  double mass_hexadecapole = gamma * (spin_dimensionless * spin_dimensionless * spin_dimensionless * spin_dimensionless) * (mass * mass * mass * mass * mass);

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_neutronstar_new(false, mass, spin, mass_quadrupole, spin_octupole, mass_hexadecapole, 0.0, 0.0, 0.0);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 0.1 * x_ind;
      double y = 0.1 * y_ind;

      if (sqrt((x * x) + (y * y)) > 0.2) {
        double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
        double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
        double **spatial_metric_prod = gkyl_malloc(sizeof(double*[3]));

        for (int i = 0; i < 3; i++) {
          spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
          inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
          spatial_metric_prod[i] = gkyl_malloc(sizeof(double[3]));

          for (int j = 0; j < 3; j++) {
            spatial_metric_prod[i][j] = 0.0;
          }
        }

        double **spacetime_metric = gkyl_malloc(sizeof(double*[4]));
        double **inv_spacetime_metric = gkyl_malloc(sizeof(double*[4]));
        double **spacetime_metric_prod = gkyl_malloc(sizeof(double*[4]));

        for (int i = 0; i < 4; i++) {
          spacetime_metric[i] = gkyl_malloc(sizeof(double[4]));
          inv_spacetime_metric[i] = gkyl_malloc(sizeof(double[4]));
          spacetime_metric_prod[i] = gkyl_malloc(sizeof(double[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_metric_prod[i][j] = 0.0;
          }
        }

        spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);
        spacetime->spacetime_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spacetime_metric);

        spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &inv_spatial_metric);
        spacetime->spacetime_inv_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &inv_spacetime_metric);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              spatial_metric_prod[i][j] += spatial_metric[i][k] * inv_spatial_metric[k][j];
            }

            if (i == j) {
              TEST_CHECK( gkyl_compare(spatial_metric_prod[i][j], 1.0, 1e-10) );
            }
            else {
              TEST_CHECK( gkyl_compare(spatial_metric_prod[i][j], 0.0, 1e-10) );
            }
          }
        }

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
              spacetime_metric_prod[i][j] += spacetime_metric[i][k] * inv_spacetime_metric[k][j];
            }

            if (i == j) {
              TEST_CHECK( gkyl_compare(spacetime_metric_prod[i][j], 1.0, 1e-10) );
            }
            else {
              TEST_CHECK( gkyl_compare(spacetime_metric_prod[i][j], 0.0, 1e-10) );
            }
          }
        }

        double spatial_metric_det;
        double spacetime_metric_det;
        double lapse_function;

        spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &spatial_metric_det);
        spacetime->spacetime_metric_det_func(spacetime, 0.0, x, y, 0.0, &spacetime_metric_det);
        spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse_function);

        TEST_CHECK( gkyl_compare(sqrt(-spacetime_metric_det), lapse_function * sqrt(spatial_metric_det), 1e-10) );

        double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
        double ***spatial_christoffel = gkyl_malloc(sizeof(double**[3]));
        double ***spatial_metric_cov_der = gkyl_malloc(sizeof(double**[3]));

        for (int i = 0; i < 3; i++) {
          spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));
          spatial_christoffel[i] = gkyl_malloc(sizeof(double*[3]));
          spatial_metric_cov_der[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
            spatial_christoffel[i][j] = gkyl_malloc(sizeof(double[3]));
            spatial_metric_cov_der[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spatial_metric_der);
        spacetime->spatial_christoffel_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spatial_christoffel);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              spatial_metric_cov_der[i][j][k] = spatial_metric_der[i][j][k];

              for (int l = 0; l < 3; l++) {
                spatial_metric_cov_der[i][j][k] -= spatial_christoffel[l][i][j] * spatial_metric[l][k];
                spatial_metric_cov_der[i][j][k] -= spatial_christoffel[l][i][k] * spatial_metric[j][l];
              }

              TEST_CHECK( gkyl_compare(spatial_metric_cov_der[i][j][k], 0.0, 1e-6) );
            }
          }
        }

        double ***spacetime_metric_der = gkyl_malloc(sizeof(double**[4]));
        double ***spacetime_christoffel = gkyl_malloc(sizeof(double**[4]));
        double ***spacetime_metric_cov_der = gkyl_malloc(sizeof(double**[4]));

        for (int i = 0; i < 4; i++) {
          spacetime_metric_der[i] = gkyl_malloc(sizeof(double*[4]));
          spacetime_christoffel[i] = gkyl_malloc(sizeof(double*[4]));
          spacetime_metric_cov_der[i] = gkyl_malloc(sizeof(double*[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_metric_der[i][j] = gkyl_malloc(sizeof(double[4]));
            spacetime_christoffel[i][j] = gkyl_malloc(sizeof(double[4]));
            spacetime_metric_cov_der[i][j] = gkyl_malloc(sizeof(double[4]));
          }
        }

        spacetime->spacetime_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_metric_der);
        spacetime->spacetime_christoffel_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_christoffel);

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
              spacetime_metric_cov_der[i][j][k] = spacetime_metric_der[i][j][k];

              for (int l = 0; l < 4; l++) {
                spacetime_metric_cov_der[i][j][k] -= spacetime_christoffel[l][i][j] * spacetime_metric[l][k];
                spacetime_metric_cov_der[i][j][k] -= spacetime_christoffel[l][i][k] * spacetime_metric[j][l];
              }

              TEST_CHECK( gkyl_compare(spacetime_metric_cov_der[i][j][k], 0.0, 1e-6) );
            }
          }
        }

        double **shift_vector_der = gkyl_malloc(sizeof(double*[3]));
        double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
        double **shift_vector_cov_der = gkyl_malloc(sizeof(double*[3]));
        double **shift_covector_cov_der = gkyl_malloc(sizeof(double*[3]));

        for (int i = 0; i < 3; i++) {
          shift_vector_der[i] = gkyl_malloc(sizeof(double[3]));
          extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
          shift_vector_cov_der[i] = gkyl_malloc(sizeof(double[3]));
          shift_covector_cov_der[i] = gkyl_malloc(sizeof(double[3]));

          for (int j = 0; j < 3; j++){
            shift_covector_cov_der[i][j] = 0.0;
          }
        }

        spacetime->shift_vector_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &shift_vector_der);
        spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &extrinsic_curvature);

        double *shift_vector = gkyl_malloc(sizeof(double[3]));
        spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift_vector);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            shift_vector_cov_der[i][j] = shift_vector_der[i][j];

            for (int k = 0; k < 3; k++) {
              shift_vector_cov_der[i][j] += spatial_christoffel[j][i][k] * shift_vector[k];
            }
          }
        }

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              shift_covector_cov_der[i][j] += spatial_metric[j][k] * shift_vector_cov_der[i][k];
            }
          }
        }

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            TEST_CHECK( gkyl_compare(2.0 * lapse_function * extrinsic_curvature[i][j], -(shift_covector_cov_der[j][i] + shift_covector_cov_der[i][j]), 1e-6) );
          }
        }

        bool in_excision_region;
        spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);

        TEST_CHECK( (in_excision_region == false) );

        for (int i = 0; i < 3; i++) {
          gkyl_free(spatial_metric[i]);
          gkyl_free(inv_spatial_metric[i]);
          gkyl_free(spatial_metric_prod[i]);
          gkyl_free(shift_vector_der[i]);
          gkyl_free(extrinsic_curvature[i]);
          gkyl_free(shift_vector_cov_der[i]);
          gkyl_free(shift_covector_cov_der[i]);
          
          for (int j = 0; j < 3; j++) {
            gkyl_free(spatial_metric_der[i][j]);
            gkyl_free(spatial_christoffel[i][j]);
            gkyl_free(spatial_metric_cov_der[i][j]);
          }
          gkyl_free(spatial_metric_der[i]);
          gkyl_free(spatial_christoffel[i]);
          gkyl_free(spatial_metric_cov_der[i]);
        }
        gkyl_free(spatial_metric);
        gkyl_free(inv_spatial_metric);
        gkyl_free(spatial_metric_prod);
        gkyl_free(spatial_metric_der);
        gkyl_free(spatial_christoffel);
        gkyl_free(spatial_metric_cov_der);
        gkyl_free(shift_vector_der);
        gkyl_free(extrinsic_curvature);
        gkyl_free(shift_vector_cov_der);
        gkyl_free(shift_covector_cov_der);
        gkyl_free(shift_vector);

        for (int i = 0; i < 4; i++) {
          gkyl_free(spacetime_metric[i]);
          gkyl_free(inv_spacetime_metric[i]);
          gkyl_free(spacetime_metric_prod[i]);

          for (int j = 0; j < 4; j++) {
            gkyl_free(spacetime_metric_der[i][j]);
            gkyl_free(spacetime_christoffel[i][j]);
            gkyl_free(spacetime_metric_cov_der[i][j]);
          }
          gkyl_free(spacetime_metric_der[i]);
          gkyl_free(spacetime_christoffel[i]);
          gkyl_free(spacetime_metric_cov_der[i]);
        }
        gkyl_free(spacetime_metric);
        gkyl_free(inv_spacetime_metric);
        gkyl_free(spacetime_metric_prod);
        gkyl_free(spacetime_metric_der);
        gkyl_free(spacetime_christoffel);
        gkyl_free(spacetime_metric_cov_der);
      }
      else {
        bool in_excision_region;
        spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);

        TEST_CHECK( (in_excision_region == true) );
      }
    }
  }

  gkyl_gr_spacetime_release(spacetime);
}

void
test_gr_brill_lindquist()
{
  double mass1 = 0.5;
  double mass2 = 0.5;
  
  double pos_x1 = -3.0;
  double pos_y1 = 0.0;
  double pos_z1 = 0.0;

  double pos_x2 = 3.0;
  double pos_y2 = 0.0;
  double pos_z2 = 0.0;

  struct gkyl_gr_spacetime *spacetime = gkyl_gr_brill_lindquist_new(false, mass1, mass2, pos_x1, pos_y1, pos_z1, pos_x2, pos_y2, pos_z2);

  for (int x_ind = -10; x_ind < 11; x_ind++) {
    for (int y_ind = -10; y_ind < 11; y_ind++) {
      double x = 1.0 * x_ind;
      double y = 1.0 * y_ind;

      if (sqrt(((x - 3.0) * (x - 3.0)) + (y * y)) > 1.0 && sqrt(((x + 3.0) * (x + 3.0)) + (y * y)) > 1.0) {
        double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
        double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
        double **spatial_metric_prod = gkyl_malloc(sizeof(double*[3]));

        for (int i = 0; i < 3; i++) {
          spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
          inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
          spatial_metric_prod[i] = gkyl_malloc(sizeof(double[3]));

          for (int j = 0; j < 3; j++) {
            spatial_metric_prod[i][j] = 0.0;
          }
        }

        double **spacetime_metric = gkyl_malloc(sizeof(double*[4]));
        double **inv_spacetime_metric = gkyl_malloc(sizeof(double*[4]));
        double **spacetime_metric_prod = gkyl_malloc(sizeof(double*[4]));

        for (int i = 0; i < 4; i++) {
          spacetime_metric[i] = gkyl_malloc(sizeof(double[4]));
          inv_spacetime_metric[i] = gkyl_malloc(sizeof(double[4]));
          spacetime_metric_prod[i] = gkyl_malloc(sizeof(double[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_metric_prod[i][j] = 0.0;
          }
        }

        spacetime->spatial_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spatial_metric);
        spacetime->spacetime_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &spacetime_metric);

        spacetime->spatial_inv_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &inv_spatial_metric);
        spacetime->spacetime_inv_metric_tensor_func(spacetime, 0.0, x, y, 0.0, &inv_spacetime_metric);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              spatial_metric_prod[i][j] += spatial_metric[i][k] * inv_spatial_metric[k][j];
            }

            if (i == j) {
              TEST_CHECK( gkyl_compare(spatial_metric_prod[i][j], 1.0, 1e-10) );
            }
            else {
              TEST_CHECK( gkyl_compare(spatial_metric_prod[i][j], 0.0, 1e-10) );
            }
          }
        }

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
              spacetime_metric_prod[i][j] += spacetime_metric[i][k] * inv_spacetime_metric[k][j];
            }

            if (i == j) {
              TEST_CHECK( gkyl_compare(spacetime_metric_prod[i][j], 1.0, 1e-10) );
            }
            else {
              TEST_CHECK( gkyl_compare(spacetime_metric_prod[i][j], 0.0, 1e-10) );
            }
          }
        }

        double spatial_metric_det;
        double spacetime_metric_det;
        double lapse_function;

        spacetime->spatial_metric_det_func(spacetime, 0.0, x, y, 0.0, &spatial_metric_det);
        spacetime->spacetime_metric_det_func(spacetime, 0.0, x, y, 0.0, &spacetime_metric_det);
        spacetime->lapse_function_func(spacetime, 0.0, x, y, 0.0, &lapse_function);

        TEST_CHECK( gkyl_compare(sqrt(-spacetime_metric_det), lapse_function * sqrt(spatial_metric_det), 1e-10) );

        double ***spatial_metric_der = gkyl_malloc(sizeof(double**[3]));
        double ***spatial_christoffel = gkyl_malloc(sizeof(double**[3]));
        double ***spatial_metric_cov_der = gkyl_malloc(sizeof(double**[3]));

        for (int i = 0; i < 3; i++) {
          spatial_metric_der[i] = gkyl_malloc(sizeof(double*[3]));
          spatial_christoffel[i] = gkyl_malloc(sizeof(double*[3]));
          spatial_metric_cov_der[i] = gkyl_malloc(sizeof(double*[3]));

          for (int j = 0; j < 3; j++) {
            spatial_metric_der[i][j] = gkyl_malloc(sizeof(double[3]));
            spatial_christoffel[i][j] = gkyl_malloc(sizeof(double[3]));
            spatial_metric_cov_der[i][j] = gkyl_malloc(sizeof(double[3]));
          }
        }

        spacetime->spatial_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spatial_metric_der);
        spacetime->spatial_christoffel_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spatial_christoffel);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              spatial_metric_cov_der[i][j][k] = spatial_metric_der[i][j][k];

              for (int l = 0; l < 3; l++) {
                spatial_metric_cov_der[i][j][k] -= spatial_christoffel[l][i][j] * spatial_metric[l][k];
                spatial_metric_cov_der[i][j][k] -= spatial_christoffel[l][i][k] * spatial_metric[j][l];
              }

              TEST_CHECK( gkyl_compare(spatial_metric_cov_der[i][j][k], 0.0, 1e-6) );
            }
          }
        }

        double ***spacetime_metric_der = gkyl_malloc(sizeof(double**[4]));
        double ***spacetime_christoffel = gkyl_malloc(sizeof(double**[4]));
        double ***spacetime_metric_cov_der = gkyl_malloc(sizeof(double**[4]));

        for (int i = 0; i < 4; i++) {
          spacetime_metric_der[i] = gkyl_malloc(sizeof(double*[4]));
          spacetime_christoffel[i] = gkyl_malloc(sizeof(double*[4]));
          spacetime_metric_cov_der[i] = gkyl_malloc(sizeof(double*[4]));

          for (int j = 0; j < 4; j++) {
            spacetime_metric_der[i][j] = gkyl_malloc(sizeof(double[4]));
            spacetime_christoffel[i][j] = gkyl_malloc(sizeof(double[4]));
            spacetime_metric_cov_der[i][j] = gkyl_malloc(sizeof(double[4]));
          }
        }

        spacetime->spacetime_metric_tensor_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_metric_der);
        spacetime->spacetime_christoffel_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &spacetime_christoffel);

        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
              spacetime_metric_cov_der[i][j][k] = spacetime_metric_der[i][j][k];

              for (int l = 0; l < 4; l++) {
                spacetime_metric_cov_der[i][j][k] -= spacetime_christoffel[l][i][j] * spacetime_metric[l][k];
                spacetime_metric_cov_der[i][j][k] -= spacetime_christoffel[l][i][k] * spacetime_metric[j][l];
              }

              TEST_CHECK( gkyl_compare(spacetime_metric_cov_der[i][j][k], 0.0, 1e-6) );
            }
          }
        }

        double **shift_vector_der = gkyl_malloc(sizeof(double*[3]));
        double **extrinsic_curvature = gkyl_malloc(sizeof(double*[3]));
        double **shift_vector_cov_der = gkyl_malloc(sizeof(double*[3]));
        double **shift_covector_cov_der = gkyl_malloc(sizeof(double*[3]));

        for (int i = 0; i < 3; i++) {
          shift_vector_der[i] = gkyl_malloc(sizeof(double[3]));
          extrinsic_curvature[i] = gkyl_malloc(sizeof(double[3]));
          shift_vector_cov_der[i] = gkyl_malloc(sizeof(double[3]));
          shift_covector_cov_der[i] = gkyl_malloc(sizeof(double[3]));

          for (int j = 0; j < 3; j++){
            shift_covector_cov_der[i][j] = 0.0;
          }
        }

        spacetime->shift_vector_der_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &shift_vector_der);
        spacetime->extrinsic_curvature_tensor_func(spacetime, 0.0, x, y, 0.0, pow(10.0, -6.0), pow(10.0, -6.0), pow(10.0, -6.0), &extrinsic_curvature);

        double *shift_vector = gkyl_malloc(sizeof(double[3]));
        spacetime->shift_vector_func(spacetime, 0.0, x, y, 0.0, &shift_vector);

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            shift_vector_cov_der[i][j] = shift_vector_der[i][j];

            for (int k = 0; k < 3; k++) {
              shift_vector_cov_der[i][j] += spatial_christoffel[j][i][k] * shift_vector[k];
            }
          }
        }

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
              shift_covector_cov_der[i][j] += spatial_metric[j][k] * shift_vector_cov_der[i][k];
            }
          }
        }

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            TEST_CHECK( gkyl_compare(2.0 * lapse_function * extrinsic_curvature[i][j], -(shift_covector_cov_der[j][i] + shift_covector_cov_der[i][j]), 1e-6) );
          }
        }

        bool in_excision_region;
        spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);

        TEST_CHECK( (in_excision_region == false) );

        for (int i = 0; i < 3; i++) {
          gkyl_free(spatial_metric[i]);
          gkyl_free(inv_spatial_metric[i]);
          gkyl_free(spatial_metric_prod[i]);
          gkyl_free(shift_vector_der[i]);
          gkyl_free(extrinsic_curvature[i]);
          gkyl_free(shift_vector_cov_der[i]);
          gkyl_free(shift_covector_cov_der[i]);
          
          for (int j = 0; j < 3; j++) {
            gkyl_free(spatial_metric_der[i][j]);
            gkyl_free(spatial_christoffel[i][j]);
            gkyl_free(spatial_metric_cov_der[i][j]);
          }
          gkyl_free(spatial_metric_der[i]);
          gkyl_free(spatial_christoffel[i]);
          gkyl_free(spatial_metric_cov_der[i]);
        }
        gkyl_free(spatial_metric);
        gkyl_free(inv_spatial_metric);
        gkyl_free(spatial_metric_prod);
        gkyl_free(spatial_metric_der);
        gkyl_free(spatial_christoffel);
        gkyl_free(spatial_metric_cov_der);
        gkyl_free(shift_vector_der);
        gkyl_free(extrinsic_curvature);
        gkyl_free(shift_vector_cov_der);
        gkyl_free(shift_covector_cov_der);
        gkyl_free(shift_vector);

        for (int i = 0; i < 4; i++) {
          gkyl_free(spacetime_metric[i]);
          gkyl_free(inv_spacetime_metric[i]);
          gkyl_free(spacetime_metric_prod[i]);

          for (int j = 0; j < 4; j++) {
            gkyl_free(spacetime_metric_der[i][j]);
            gkyl_free(spacetime_christoffel[i][j]);
            gkyl_free(spacetime_metric_cov_der[i][j]);
          }
          gkyl_free(spacetime_metric_der[i]);
          gkyl_free(spacetime_christoffel[i]);
          gkyl_free(spacetime_metric_cov_der[i]);
        }
        gkyl_free(spacetime_metric);
        gkyl_free(inv_spacetime_metric);
        gkyl_free(spacetime_metric_prod);
        gkyl_free(spacetime_metric_der);
        gkyl_free(spacetime_christoffel);
        gkyl_free(spacetime_metric_cov_der);
      }
      else {
        bool in_excision_region;
        spacetime->excision_region_func(spacetime, 0.0, x, y, 0.0, &in_excision_region);

        TEST_CHECK( (in_excision_region == true) );
      }
    }
  }

  gkyl_gr_spacetime_release(spacetime);
}

TEST_LIST = {
  { "gr_minkowski", test_gr_minkowski },
  { "gr_schwarzschild", test_gr_schwarzschild },
  //{ "gr_kerr", test_gr_kerr },
  { "gr_tov_static", test_gr_tov_static },
  { "gr_neutronstar_static", test_gr_neutronstar_static },
  { "gr_neutronstar_spinning", test_gr_neutronstar_spinning },
  { "gr_brill_lindquist", test_gr_brill_lindquist },
  { NULL, NULL },
};
