#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_vacuum_einstein.h>
#include <gkyl_wv_vacuum_einstein_priv.h>

void
gkyl_vacuum_einstein_flux(enum gkyl_spacetime_slicing spacetime_slicing, enum gkyl_spacetime_evolution spacetime_evolution, const double q[64], double flux[64])
{
  double spatial_metric[3][3];
  spatial_metric[0][0] = q[0]; spatial_metric[0][1] = q[1]; spatial_metric[0][2] = q[2];
  spatial_metric[1][0] = q[3]; spatial_metric[1][1] = q[4]; spatial_metric[1][2] = q[5];
  spatial_metric[2][0] = q[6]; spatial_metric[2][1] = q[7]; spatial_metric[2][2] = q[8];

  double lapse = q[9];

  double extrinsic_curvature[3][3];
  extrinsic_curvature[0][0] = q[10]; extrinsic_curvature[0][1] = q[11]; extrinsic_curvature[0][2] = q[12];
  extrinsic_curvature[1][0] = q[13]; extrinsic_curvature[1][1] = q[14]; extrinsic_curvature[1][2] = q[15];
  extrinsic_curvature[2][0] = q[16]; extrinsic_curvature[2][1] = q[17]; extrinsic_curvature[2][2] = q[18];

  double spatial_metric_der[3][3][3];
  spatial_metric_der[0][0][0] = q[19]; spatial_metric_der[0][0][1] = q[20]; spatial_metric_der[0][0][2] = q[21];
  spatial_metric_der[0][1][0] = q[22]; spatial_metric_der[0][1][1] = q[23]; spatial_metric_der[0][1][2] = q[24];
  spatial_metric_der[0][2][0] = q[25]; spatial_metric_der[0][2][1] = q[26]; spatial_metric_der[0][2][2] = q[27];

  spatial_metric_der[1][0][0] = q[28]; spatial_metric_der[1][0][1] = q[29]; spatial_metric_der[1][0][2] = q[30];
  spatial_metric_der[1][1][0] = q[31]; spatial_metric_der[1][1][1] = q[32]; spatial_metric_der[1][1][2] = q[33];
  spatial_metric_der[1][2][0] = q[34]; spatial_metric_der[1][2][1] = q[35]; spatial_metric_der[1][2][2] = q[36];

  spatial_metric_der[2][0][0] = q[37]; spatial_metric_der[2][0][1] = q[38]; spatial_metric_der[2][0][2] = q[39];
  spatial_metric_der[2][1][0] = q[40]; spatial_metric_der[2][1][1] = q[41]; spatial_metric_der[2][1][2] = q[42];
  spatial_metric_der[2][2][0] = q[43]; spatial_metric_der[2][2][1] = q[44]; spatial_metric_der[2][2][2] = q[45];

  double lapse_der[3];
  lapse_der[0] = q[46];
  lapse_der[1] = q[47];
  lapse_der[2] = q[48];

  double aux_vect[3];
  aux_vect[0] = q[49];
  aux_vect[1] = q[50];
  aux_vect[2] = q[51];

  double shift_vect[3];
  shift_vect[0] = q[52];
  shift_vect[1] = q[53];
  shift_vect[2] = q[54];

  double shift_vect_der[3][3];
  shift_vect_der[0][0] = q[55]; shift_vect_der[0][1] = q[56]; shift_vect_der[0][2] = q[57];
  shift_vect_der[1][0] = q[58]; shift_vect_der[1][1] = q[59]; shift_vect_der[1][2] = q[60];
  shift_vect_der[2][0] = q[61]; shift_vect_der[2][1] = q[62]; shift_vect_der[2][2] = q[63];

  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_vacuum_einstein_inv_spatial_metric(q, &inv_spatial_metric);

  double evolution_func = 0.0;
  if (spacetime_evolution == GKYL_RICCI_EVOLUTION) {
    evolution_func = 0.0;
  }
  else if (spacetime_evolution == GKYL_EINSTEIN_EVOLUTION) {
    evolution_func = 1.0;
  }

  double extrinsic_curvature_trace = 0.0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      extrinsic_curvature_trace += inv_spatial_metric[i][j] * extrinsic_curvature[i][j];
    }
  }

  double slicing_func = 0.0;
  if (spacetime_slicing == GKYL_GEODESIC_SLICING) {
    slicing_func = 0.0;
  }
  else if (spacetime_slicing == GKYL_HARMONIC_SLICING) {
    slicing_func = extrinsic_curvature_trace;
  }
  else if (spacetime_slicing == GKYL_1PLUSLOG_SLICING) {
    slicing_func = extrinsic_curvature_trace / lapse;
  }

  double spatial_metric_der_raised1[3][3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        spatial_metric_der_raised1[k][i][j] = 0.0;

        for (int l = 0; l < 3; l++) {
          spatial_metric_der_raised1[k][i][j] += inv_spatial_metric[k][l] * spatial_metric_der[l][i][j];
        }
      }
    }
  }

  double spatial_metric_der_raised3[3][3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        spatial_metric_der_raised1[i][j][k] = 0.0;

        for (int l = 0; l < 3; l++) {
          spatial_metric_der_raised1[i][j][k] += inv_spatial_metric[l][k] * spatial_metric_der[i][j][l];
        }
      }
    }
  }

  double aux_vect_raised[3];
  for (int k = 0; k < 3; k++) {
    aux_vect_raised[k] = 0.0;
    
    for (int l = 0; l < 3; l++) {
     aux_vect_raised[k] += inv_spatial_metric[k][l] * aux_vect[l];
    }
  }

  double shift_vect_der_lowered[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      shift_vect_der_lowered[i][j] = 0.0;

      for (int k = 0; k < 3; k++) {
        shift_vect_der_lowered[i][j] += spatial_metric[k][j] * shift_vect_der[i][k];
      }
    }
  }

  double shift_vect_der_switched[3][3];
  for (int i = 0; i < 3; i++) {
    for (int k = 0; k < 3; k++) {
      shift_vect_der_switched[i][k] = 0.0;

      for (int l = 0; l < 3; l++) {
        for (int m = 0; m < 3; m++) {
          shift_vect_der_switched[i][k] += inv_spatial_metric[i][l] * spatial_metric[m][k] * shift_vect_der[l][m];
        }
      }
    }
  }

  double symmetrized_shift[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      symmetrized_shift[i][j] = (1.0 / lapse) * (shift_vect_der_lowered[i][j] + shift_vect_der_lowered[j][i]);
    }
  }

  double extrinsic_curvature_flux[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      extrinsic_curvature_flux[i][j] = -shift_vect[0] * extrinsic_curvature[i][j];
      extrinsic_curvature_flux[i][j] += lapse * spatial_metric_der_raised1[0][i][j];
      extrinsic_curvature_flux[i][j] -= lapse * (0.5 * evolution_func) * aux_vect_raised[0] * spatial_metric[i][j];

      if (i == 0) {
        extrinsic_curvature_flux[i][j] += 0.5 * lapse * lapse_der[j];
        extrinsic_curvature_flux[i][j] += lapse * aux_vect[j];
        for (int r = 0; r < 3; r++) {
          extrinsic_curvature_flux[i][j] -= 0.5 * lapse * spatial_metric_der_raised3[j][r][r];
        }
      }

      if (j == 0) {
        extrinsic_curvature_flux[i][j] += 0.5 * lapse * lapse_der[i];
        extrinsic_curvature_flux[i][j] += lapse * aux_vect[i];
        for (int r = 0; r < 3; r++) {
          extrinsic_curvature_flux[i][j] -= 0.5 * lapse * spatial_metric_der_raised3[i][r][r];
        }
      }
    }
  }

  double spatial_metric_der_flux[3][3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        spatial_metric_der_flux[k][i][j] = 0.0;
      }

      for (int r = 0; r < 3; r++) {
        spatial_metric_der_flux[0][i][j] -= shift_vect[r] * spatial_metric_der[r][i][j];
      }

      spatial_metric_der_flux[0][i][j] += lapse * (extrinsic_curvature[i][j] - symmetrized_shift[i][j]);
    }
  }

  double lapse_der_flux[3];
  for (int i = 0; i < 3; i++) {
    lapse_der_flux[i] = 0.0;
  }
  for (int r = 0; r < 3; r++) {
    lapse_der_flux[0] -= shift_vect[r] * lapse_der[r];
  }
  lapse_der_flux[0] += lapse * slicing_func;

  double aux_vect_flux[3];
  for (int i = 0; i < 3; i++) {
    aux_vect_flux[i] = -shift_vect[0] * aux_vect[i];
    aux_vect_flux[i] += shift_vect_der_switched[0][i];
    aux_vect_flux[i] -= shift_vect_der[i][0];
  }
  
  for (int i = 0; i < 10; i++) {
    flux[i] = 0.0;
  }

  flux[10] = extrinsic_curvature_flux[0][0]; flux[11] = extrinsic_curvature_flux[0][1]; flux[12] = extrinsic_curvature_flux[0][2];
  flux[13] = extrinsic_curvature_flux[1][0]; flux[14] = extrinsic_curvature_flux[1][1]; flux[15] = extrinsic_curvature_flux[1][2];
  flux[16] = extrinsic_curvature_flux[2][0]; flux[17] = extrinsic_curvature_flux[2][1]; flux[18] = extrinsic_curvature_flux[2][2];

  flux[19] = spatial_metric_der_flux[0][0][0]; flux[20] = spatial_metric_der_flux[0][0][1]; flux[21] = spatial_metric_der_flux[0][0][2];
  flux[22] = spatial_metric_der_flux[0][1][0]; flux[23] = spatial_metric_der_flux[0][1][1]; flux[24] = spatial_metric_der_flux[0][1][2];
  flux[25] = spatial_metric_der_flux[0][2][0]; flux[26] = spatial_metric_der_flux[0][2][1]; flux[27] = spatial_metric_der_flux[0][2][2];

  flux[28] = spatial_metric_der_flux[1][0][0]; flux[29] = spatial_metric_der_flux[1][0][1]; flux[30] = spatial_metric_der_flux[1][0][2];
  flux[31] = spatial_metric_der_flux[1][1][0]; flux[32] = spatial_metric_der_flux[1][1][1]; flux[33] = spatial_metric_der_flux[1][1][2];
  flux[34] = spatial_metric_der_flux[1][2][0]; flux[35] = spatial_metric_der_flux[1][2][1]; flux[36] = spatial_metric_der_flux[1][2][2];

  flux[37] = spatial_metric_der_flux[2][0][0]; flux[38] = spatial_metric_der_flux[2][0][1]; flux[39] = spatial_metric_der_flux[2][0][2];
  flux[40] = spatial_metric_der_flux[2][1][0]; flux[41] = spatial_metric_der_flux[2][1][1]; flux[42] = spatial_metric_der_flux[2][1][2];
  flux[43] = spatial_metric_der_flux[2][2][0]; flux[44] = spatial_metric_der_flux[2][2][1]; flux[45] = spatial_metric_der_flux[2][2][2];

  flux[46] = lapse_der_flux[0];
  flux[47] = lapse_der_flux[1];
  flux[48] = lapse_der_flux[2];

  flux[49] = aux_vect_flux[0];
  flux[50] = aux_vect_flux[1];
  flux[51] = aux_vect_flux[2];

  for (int i = 52; i < 64; i++) {
    flux[i] = 0.0;
  }
}

void 
gkyl_vacuum_einstein_inv_spatial_metric(const double q[64], double ***inv_spatial_metric)
{
  double spatial_metric[3][3];
  spatial_metric[0][0] = q[0]; spatial_metric[0][1] = q[1]; spatial_metric[0][2] = q[2];
  spatial_metric[1][0] = q[3]; spatial_metric[1][1] = q[4]; spatial_metric[1][2] = q[5];
  spatial_metric[2][0] = q[6]; spatial_metric[2][1] = q[7]; spatial_metric[2][2] = q[8];

  double spatial_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));
  
  double trace = 0.0;
  for (int i = 0; i < 3; i++) {
    trace += spatial_metric[i][i];
  }

  double spatial_metric_sq[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      spatial_metric_sq[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        spatial_metric_sq[i][j] += spatial_metric[i][k] * spatial_metric[k][j];
      }
    }
  }

  double sq_trace = 0.0;
  for (int i = 0; i < 3; i++) {
    sq_trace += spatial_metric_sq[i][i];
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
      (*inv_spatial_metric)[i][j] = (1.0 / spatial_det) *
        ((0.5 * ((trace * trace) - sq_trace) * euclidean_metric[i][j]) - (trace * spatial_metric[i][j]) + spatial_metric_sq[i][j]);
    }
  }
}

static inline double
gkyl_vacuum_einstein_max_abs_speed(enum gkyl_spacetime_slicing spacetime_slicing, const double q[64])
{
  double lapse = q[9];

  double shift_vect[3];
  shift_vect[0] = q[52];
  shift_vect[1] = q[53];
  shift_vect[2] = q[54];

  double **inv_spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    inv_spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  gkyl_vacuum_einstein_inv_spatial_metric(q, &inv_spatial_metric);

  double slicing_func = 0.0;
  if (spacetime_slicing == GKYL_GEODESIC_SLICING) {
    slicing_func = 0.0;
  }
  else if (spacetime_slicing == GKYL_HARMONIC_SLICING) {
    slicing_func = 1.0;
  }
  else if (spacetime_slicing == GKYL_1PLUSLOG_SLICING) {
    slicing_func = 1.0 / lapse;
  }

  double characteristic1 = -shift_vect[0] + (lapse * sqrt(slicing_func * inv_spatial_metric[0][0]));
  double characteristic2 = -shift_vect[0] - (lapse * sqrt(slicing_func * inv_spatial_metric[0][0]));

  return fmax(fabs(characteristic1), fabs(characteristic2));
}

static inline void
cons_to_riem(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* qin, double* wout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 64; i++) {
    wout[i] = qin[i];
  }
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double* qout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 64; i++) {
    qout[i] = win[i];
  }
}

void
gkyl_vacuum_einstein_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct wv_vacuum_einstein *vacuum_einstein = container_of(base->on_dev, struct wv_vacuum_einstein, eqn);
    gkyl_cu_free(vacuum_einstein);
  }

  struct wv_vacuum_einstein *vacuum_einstein = container_of(base, struct wv_vacuum_einstein, eqn);
  gkyl_free(vacuum_einstein);
}

struct gkyl_wv_eqn*
gkyl_wv_vacuum_einstein_new(enum gkyl_spacetime_slicing spacetime_slicing, enum gkyl_spacetime_evolution spacetime_evolution, bool use_gpu)
{
  return gkyl_wv_vacuum_einstein_inew(&(struct gkyl_wv_vacuum_einstein_inp) {
      .spacetime_slicing = spacetime_slicing,
      .spacetime_evolution = spacetime_evolution,
      .rp_type = WV_VACUUM_EINSTEIN_RP_LAX,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_vacuum_einstein_inew(const struct gkyl_wv_vacuum_einstein_inp* inp)
{
  struct wv_vacuum_einstein *vacuum_einstein = gkyl_malloc(sizeof(struct wv_vacuum_einstein));

  vacuum_einstein->eqn.type = GKYL_EQN_VACUUM_EINSTEIN;
  vacuum_einstein->eqn.num_equations = 64;
  vacuum_einstein->eqn.num_diag = 4;

  vacuum_einstein->spacetime_slicing = inp->spacetime_slicing;
  vacuum_einstein->spacetime_evolution = inp->spacetime_evolution;

  if (inp->rp_type == WV_VACUUM_EINSTEIN_RP_LAX) {
    vacuum_einstein->eqn.num_waves = 2;
  }

  vacuum_einstein->eqn.cons_to_riem = cons_to_riem;
  vacuum_einstein->eqn.riem_to_cons = riem_to_cons;

  vacuum_einstein->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(vacuum_einstein->eqn.flags);
  vacuum_einstein->eqn.ref_count = gkyl_ref_count_init(gkyl_vacuum_einstein_free);
  vacuum_einstein->eqn.on_dev = &vacuum_einstein->eqn; // On the CPU, the equation object points to itself.

  return &vacuum_einstein->eqn;
}

enum gkyl_spacetime_slicing
gkyl_wv_vacuum_einstein_spacetime_slicing(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_vacuum_einstein *vacuum_einstein = container_of(eqn, struct wv_vacuum_einstein, eqn);
  enum gkyl_spacetime_slicing spacetime_slicing = vacuum_einstein->spacetime_slicing;

  return spacetime_slicing;
}

enum gkyl_spacetime_evolution
gkyl_wv_vacuum_einstein_spacetime_evolution(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_vacuum_einstein *vacuum_einstein = container_of(eqn, struct wv_vacuum_einstein, eqn);
  enum gkyl_spacetime_evolution spacetime_evolution = vacuum_einstein->spacetime_evolution;

  return spacetime_evolution;
}