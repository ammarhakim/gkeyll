#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_vacuum_einstein_conformal.h>
#include <gkyl_wv_vacuum_einstein_conformal_priv.h>

void
gkyl_vacuum_einstein_conformal_flux(double excision_threshold, enum gkyl_spacetime_slicing spacetime_slicing, enum gkyl_spacetime_evolution spacetime_evolution,
  const double q[64], double flux[64])
{
  double conformal_fact = 1.0; // Hardcode the conformal factor to 1 everywhere, for now.
  double conformal_fact_der[3];
  for (int i = 0; i < 3; i++) {
    conformal_fact_der[i] = 0.0;
  }

  double conformal_spatial_metric[3][3];
  conformal_spatial_metric[0][0] = q[0]; conformal_spatial_metric[0][1] = q[1]; conformal_spatial_metric[0][2] = q[2];
  conformal_spatial_metric[1][0] = q[3]; conformal_spatial_metric[1][1] = q[4]; conformal_spatial_metric[1][2] = q[5];
  conformal_spatial_metric[2][0] = q[6]; conformal_spatial_metric[2][1] = q[7]; conformal_spatial_metric[2][2] = q[8];

  double conformal_lapse = q[9];

  double conformal_extrinsic_curvature[3][3];
  conformal_extrinsic_curvature[0][0] = q[10]; conformal_extrinsic_curvature[0][1] = q[11]; conformal_extrinsic_curvature[0][2] = q[12];
  conformal_extrinsic_curvature[1][0] = q[13]; conformal_extrinsic_curvature[1][1] = q[14]; conformal_extrinsic_curvature[1][2] = q[15];
  conformal_extrinsic_curvature[2][0] = q[16]; conformal_extrinsic_curvature[2][1] = q[17]; conformal_extrinsic_curvature[2][2] = q[18];

  double conformal_spatial_metric_der[3][3][3];
  conformal_spatial_metric_der[0][0][0] = q[19]; conformal_spatial_metric_der[0][0][1] = q[20]; conformal_spatial_metric_der[0][0][2] = q[21];
  conformal_spatial_metric_der[0][1][0] = q[22]; conformal_spatial_metric_der[0][1][1] = q[23]; conformal_spatial_metric_der[0][1][2] = q[24];
  conformal_spatial_metric_der[0][2][0] = q[25]; conformal_spatial_metric_der[0][2][1] = q[26]; conformal_spatial_metric_der[0][2][2] = q[27];

  conformal_spatial_metric_der[1][0][0] = q[28]; conformal_spatial_metric_der[1][0][1] = q[29]; conformal_spatial_metric_der[1][0][2] = q[30];
  conformal_spatial_metric_der[1][1][0] = q[31]; conformal_spatial_metric_der[1][1][1] = q[32]; conformal_spatial_metric_der[1][1][2] = q[33];
  conformal_spatial_metric_der[1][2][0] = q[34]; conformal_spatial_metric_der[1][2][1] = q[35]; conformal_spatial_metric_der[1][2][2] = q[36];

  conformal_spatial_metric_der[2][0][0] = q[37]; conformal_spatial_metric_der[2][0][1] = q[38]; conformal_spatial_metric_der[2][0][2] = q[39];
  conformal_spatial_metric_der[2][1][0] = q[40]; conformal_spatial_metric_der[2][1][1] = q[41]; conformal_spatial_metric_der[2][1][2] = q[42];
  conformal_spatial_metric_der[2][2][0] = q[43]; conformal_spatial_metric_der[2][2][1] = q[44]; conformal_spatial_metric_der[2][2][2] = q[45];

  double conformal_lapse_der[3];
  conformal_lapse_der[0] = q[46];
  conformal_lapse_der[1] = q[47];
  conformal_lapse_der[2] = q[48];

  double conformal_aux_vect[3];
  conformal_aux_vect[0] = q[49];
  conformal_aux_vect[1] = q[50];
  conformal_aux_vect[2] = q[51];

  double conformal_shift_vect[3];
  conformal_shift_vect[0] = q[52];
  conformal_shift_vect[1] = q[53];
  conformal_shift_vect[2] = q[54];

  double conformal_shift_vect_der[3][3];
  conformal_shift_vect_der[0][0] = q[55]; conformal_shift_vect_der[0][1] = q[56]; conformal_shift_vect_der[0][2] = q[57];
  conformal_shift_vect_der[1][0] = q[58]; conformal_shift_vect_der[1][1] = q[59]; conformal_shift_vect_der[1][2] = q[60];
  conformal_shift_vect_der[2][0] = q[61]; conformal_shift_vect_der[2][1] = q[62]; conformal_shift_vect_der[2][2] = q[63];

  bool in_excision_region = false;
  if (conformal_lapse < excision_threshold) {
    in_excision_region = true;
  }

  if (!in_excision_region) {
    double **inv_conformal_spatial_metric = gkyl_malloc(sizeof(double*[3]));
    for (int i = 0; i < 3; i++) {
      inv_conformal_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
    }

    gkyl_vacuum_einstein_conformal_inv_spatial_metric(q, &inv_conformal_spatial_metric);

    double evolution_func = 0.0;
    if (spacetime_evolution == GKYL_RICCI_EVOLUTION) {
      evolution_func = 0.0;
    }
    else if (spacetime_evolution == GKYL_EINSTEIN_EVOLUTION) {
      evolution_func = 1.0;
    }

    double conformal_extrinsic_curvature_trace = 0.0;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        conformal_extrinsic_curvature_trace += inv_conformal_spatial_metric[i][j] * conformal_extrinsic_curvature[i][j];
      }
    }

    double slicing_func = 0.0;
    if (spacetime_slicing == GKYL_GEODESIC_SLICING) {
      slicing_func = 0.0;
    }
    else if (spacetime_slicing == GKYL_HARMONIC_SLICING) {
      slicing_func = conformal_extrinsic_curvature_trace / (conformal_fact * conformal_fact * conformal_fact * conformal_fact);
    }
    else if (spacetime_slicing == GKYL_1PLUSLOG_SLICING) {
      slicing_func = 2.0 * conformal_extrinsic_curvature_trace / (conformal_lapse * (conformal_fact * conformal_fact * conformal_fact * conformal_fact));
    }

    double conformal_spatial_metric_der_raised1[3][3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          conformal_spatial_metric_der_raised1[k][i][j] = 0.0;
          
          for (int l = 0; l < 3; l++) {
            conformal_spatial_metric_der_raised1[k][i][j] += inv_conformal_spatial_metric[k][l] * conformal_spatial_metric_der[l][i][j];
          }
        }
      }
    }

    double conformal_spatial_metric_der_raised3[3][3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          conformal_spatial_metric_der_raised3[i][j][k] = 0.0;

          for (int l = 0; l < 3; l++) {
            conformal_spatial_metric_der_raised3[i][j][k] += inv_conformal_spatial_metric[l][k] * conformal_spatial_metric_der[i][j][l];
          }
        }
      }
    }

    double conformal_aux_vect_raised[3];
    for (int k = 0; k < 3; k++) {
      conformal_aux_vect_raised[k] = 0.0;
        
      for (int l = 0; l < 3; l++) {
        conformal_aux_vect_raised[k] += inv_conformal_spatial_metric[k][l] * conformal_aux_vect[l];
      }
    }

    double conformal_shift_vect_der_lowered[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        conformal_shift_vect_der_lowered[i][j] = 0.0;

        for (int k = 0; k < 3; k++) {
          conformal_shift_vect_der_lowered[i][j] += conformal_spatial_metric[k][j] * conformal_shift_vect_der[i][k];
        }
      }
    }

    double conformal_shift_vect_der_switched[3][3];
    for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
        conformal_shift_vect_der_switched[i][k] = 0.0;

        for (int l = 0; l < 3; l++) {
          for (int m = 0; m < 3; m++) {
            conformal_shift_vect_der_switched[i][k] += inv_conformal_spatial_metric[i][l] * conformal_spatial_metric[m][k] * conformal_shift_vect_der[l][m];
          }
        }
      }
    }

    double conformal_symmetrized_shift[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        conformal_symmetrized_shift[i][j] = (1.0 / conformal_lapse) * (conformal_shift_vect_der_lowered[i][j] + conformal_shift_vect_der_lowered[j][i]);
      }
    }

    double conformal_extrinsic_curvature_flux[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        conformal_extrinsic_curvature_flux[i][j] = -conformal_shift_vect[0] * conformal_extrinsic_curvature[i][j];
        conformal_extrinsic_curvature_flux[i][j] += conformal_lapse * conformal_spatial_metric_der_raised1[0][i][j];
        conformal_extrinsic_curvature_flux[i][j] -= conformal_lapse * (0.5 * evolution_func) * conformal_aux_vect_raised[0] * conformal_spatial_metric[i][j];

        if (i == 0) {
          conformal_extrinsic_curvature_flux[i][j] += 0.5 * conformal_lapse * conformal_lapse_der[j];
          conformal_extrinsic_curvature_flux[i][j] += conformal_lapse * conformal_aux_vect[j];
          for (int r = 0; r < 3; r++) {
            conformal_extrinsic_curvature_flux[i][j] -= 0.5 * conformal_lapse * conformal_spatial_metric_der_raised3[j][r][r];
          }
        }

        if (j == 0) {
          conformal_extrinsic_curvature_flux[i][j] += 0.5 * conformal_lapse * conformal_lapse_der[i];
          conformal_extrinsic_curvature_flux[i][j] += conformal_lapse * conformal_aux_vect[i];
          for (int r = 0; r < 3; r++) {
            conformal_extrinsic_curvature_flux[i][j] -= 0.5 * conformal_lapse * conformal_spatial_metric_der_raised3[i][r][r];
          }
        }
      }
    }

    double conformal_spatial_metric_der_flux[3][3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          conformal_spatial_metric_der_flux[k][i][j] = 0.0;
        }

        for (int r = 0; r < 3; r++) {
          conformal_spatial_metric_der_flux[0][i][j] -= conformal_shift_vect[r] * conformal_spatial_metric_der[r][i][j];
          conformal_spatial_metric_der_flux[0][i][j] -= 2.0 * conformal_fact_der[r] * conformal_spatial_metric[i][j];
        }

        conformal_spatial_metric_der_flux[0][i][j] += conformal_lapse * (conformal_extrinsic_curvature[i][j] /
          (conformal_fact * conformal_fact * conformal_fact * conformal_fact));
        conformal_spatial_metric_der_flux[0][i][j] -= conformal_lapse * conformal_symmetrized_shift[i][j];
      }
    }

    double conformal_lapse_der_flux[3];
    for (int i = 0; i < 3; i++) {
      conformal_lapse_der_flux[i] = 0.0;
    }
    for (int r = 0; r < 3; r++) {
      conformal_lapse_der_flux[0] -= conformal_shift_vect[r] * conformal_lapse_der[r];
    }
    conformal_lapse_der_flux[0] += conformal_lapse * slicing_func;

    double conformal_aux_vect_flux[3];
    for (int i = 0; i < 3; i++) {
      conformal_aux_vect_flux[i] = -conformal_shift_vect[0] * conformal_aux_vect[i];
      conformal_aux_vect_flux[i] -= 4.0 * conformal_shift_vect[0] * conformal_fact_der[i];
      conformal_aux_vect_flux[i] += conformal_shift_vect_der_switched[0][i];
      conformal_aux_vect_flux[i] -= conformal_shift_vect_der[i][0];
    }
    
    for (int i = 0; i < 10; i++) {
      flux[i] = 0.0;
    }

    flux[10] = conformal_extrinsic_curvature_flux[0][0]; flux[11] = conformal_extrinsic_curvature_flux[0][1]; flux[12] = conformal_extrinsic_curvature_flux[0][2];
    flux[13] = conformal_extrinsic_curvature_flux[1][0]; flux[14] = conformal_extrinsic_curvature_flux[1][1]; flux[15] = conformal_extrinsic_curvature_flux[1][2];
    flux[16] = conformal_extrinsic_curvature_flux[2][0]; flux[17] = conformal_extrinsic_curvature_flux[2][1]; flux[18] = conformal_extrinsic_curvature_flux[2][2];

    flux[19] = conformal_spatial_metric_der_flux[0][0][0]; flux[20] = conformal_spatial_metric_der_flux[0][0][1]; flux[21] = conformal_spatial_metric_der_flux[0][0][2];
    flux[22] = conformal_spatial_metric_der_flux[0][1][0]; flux[23] = conformal_spatial_metric_der_flux[0][1][1]; flux[24] = conformal_spatial_metric_der_flux[0][1][2];
    flux[25] = conformal_spatial_metric_der_flux[0][2][0]; flux[26] = conformal_spatial_metric_der_flux[0][2][1]; flux[27] = conformal_spatial_metric_der_flux[0][2][2];

    flux[28] = conformal_spatial_metric_der_flux[1][0][0]; flux[29] = conformal_spatial_metric_der_flux[1][0][1]; flux[30] = conformal_spatial_metric_der_flux[1][0][2];
    flux[31] = conformal_spatial_metric_der_flux[1][1][0]; flux[32] = conformal_spatial_metric_der_flux[1][1][1]; flux[33] = conformal_spatial_metric_der_flux[1][1][2];
    flux[34] = conformal_spatial_metric_der_flux[1][2][0]; flux[35] = conformal_spatial_metric_der_flux[1][2][1]; flux[36] = conformal_spatial_metric_der_flux[1][2][2];

    flux[37] = conformal_spatial_metric_der_flux[2][0][0]; flux[38] = conformal_spatial_metric_der_flux[2][0][1]; flux[39] = conformal_spatial_metric_der_flux[2][0][2];
    flux[40] = conformal_spatial_metric_der_flux[2][1][0]; flux[41] = conformal_spatial_metric_der_flux[2][1][1]; flux[42] = conformal_spatial_metric_der_flux[2][1][2];
    flux[43] = conformal_spatial_metric_der_flux[2][2][0]; flux[44] = conformal_spatial_metric_der_flux[2][2][1]; flux[45] = conformal_spatial_metric_der_flux[2][2][2];

    flux[46] = conformal_lapse_der_flux[0];
    flux[47] = conformal_lapse_der_flux[1];
    flux[48] = conformal_lapse_der_flux[2];

    flux[49] = conformal_aux_vect_flux[0];
    flux[50] = conformal_aux_vect_flux[1];
    flux[51] = conformal_aux_vect_flux[2];

    for (int i = 52; i < 64; i++) {
      flux[i] = 0.0;
    }
    
    for (int i = 0; i < 3; i++) {
      gkyl_free(inv_conformal_spatial_metric[i]);
    }
    gkyl_free(inv_conformal_spatial_metric);
  }
  else {
    for (int i = 0; i < 64; i++) {
      flux[i] = 0.0;
    }
  }
}

void 
gkyl_vacuum_einstein_conformal_inv_spatial_metric(const double q[64], double ***inv_conformal_spatial_metric)
{
  double conformal_spatial_metric[3][3];
  conformal_spatial_metric[0][0] = q[0]; conformal_spatial_metric[0][1] = q[1]; conformal_spatial_metric[0][2] = q[2];
  conformal_spatial_metric[1][0] = q[3]; conformal_spatial_metric[1][1] = q[4]; conformal_spatial_metric[1][2] = q[5];
  conformal_spatial_metric[2][0] = q[6]; conformal_spatial_metric[2][1] = q[7]; conformal_spatial_metric[2][2] = q[8];

  double conformal_spatial_det =
    (conformal_spatial_metric[0][0] * ((conformal_spatial_metric[1][1] * conformal_spatial_metric[2][2]) - (conformal_spatial_metric[2][1] * conformal_spatial_metric[1][2]))) -
    (conformal_spatial_metric[0][1] * ((conformal_spatial_metric[1][0] * conformal_spatial_metric[2][2]) - (conformal_spatial_metric[1][2] * conformal_spatial_metric[2][0]))) +
    (conformal_spatial_metric[0][2] * ((conformal_spatial_metric[1][0] * conformal_spatial_metric[2][1]) - (conformal_spatial_metric[1][1] * conformal_spatial_metric[2][0])));
  
  double trace = 0.0;
  for (int i = 0; i < 3; i++) {
    trace += conformal_spatial_metric[i][i];
  }

  double conformal_spatial_metric_sq[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      conformal_spatial_metric_sq[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        conformal_spatial_metric_sq[i][j] += conformal_spatial_metric[i][k] * conformal_spatial_metric[k][j];
      }
    }
  }

  double sq_trace = 0.0;
  for (int i = 0; i < 3; i++) {
    sq_trace += conformal_spatial_metric_sq[i][i];
  }

  double euclidean_metric[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        euclidean_metric[i][j] = 1.0;
      }
      else {
        euclidean_metric[i][j] = 0.0;
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*inv_conformal_spatial_metric)[i][j] = (1.0 / conformal_spatial_det) *
        ((0.5 * ((trace * trace) - sq_trace) * euclidean_metric[i][j]) - (trace * conformal_spatial_metric[i][j]) + conformal_spatial_metric_sq[i][j]);
    }
  }
}