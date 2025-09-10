#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_ten_moment_nn_closure.h>
#include <gkyl_moment_non_ideal_priv.h>

// Makes indexing cleaner.
static const unsigned Q111 = 0;
static const unsigned Q112 = 1;
static const unsigned Q113 = 2;
static const unsigned Q122 = 3;
static const unsigned Q123 = 4;
static const unsigned Q133 = 5;
static const unsigned Q222 = 6;
static const unsigned Q223 = 7;
static const unsigned Q233 = 8;
static const unsigned Q333 = 9;

// 1D stencil locations (L: lower, U: upper).
enum loc_1d {
  L_1D, U_1D
};

// 1D, second-order stencil locations (L2: lower, C2: center, U2: upper).
enum loc_1d_second_order {
  L2_1D, C2_1D, U2_1D
};

// 2D stencil locations (L: lower, U: upper).
enum loc_2d {
  LL_2D, LU_2D,
  UL_2D, UU_2D
};

struct gkyl_ten_moment_nn_closure
{
  struct gkyl_rect_grid grid; // Grid on which to solve equations.
  int ndim; // Number of dimensions.
  int poly_order; // Polynomial order of learned DG coefficients.
  double k0; // Damping coefficient.
  kann_t* ann; // Neural network architecture.
};

static void
create_offsets_vertices(const struct gkyl_range* range, long offsets[])
{
  // Box-spanning stencil.
  struct gkyl_range box3;
  gkyl_range_init(&box3, range->ndim, (int[]) { -1, -1, -1 }, (int[]) { 0, 0, 0 });

  struct gkyl_range_iter iter3;
  gkyl_range_iter_init(&iter3, &box3);

  // Construct a list of offsets.
  int count = 0;
  while (gkyl_range_iter_next(&iter3)) {
    offsets[count] = gkyl_range_offset(range, iter3.idx);
    count += 1;
  }
}

static void
create_offsets_centers(const struct gkyl_range *range, long offsets[])
{
  // Box-spanning stencil.
  struct gkyl_range box3;
  gkyl_range_init(&box3, range->ndim, (int[]) { 0, 0, 0 }, (int[]) { 1, 1, 1 });

  struct gkyl_range_iter iter3;
  gkyl_range_iter_init(&iter3, &box3);

  // Construct a list of offsets.
  int count = 0;
  while (gkyl_range_iter_next(&iter3)) {
    offsets[count] = gkyl_range_offset(range, iter3.idx);
    count += 1;
  }
}

static void
var_setup(const gkyl_ten_moment_nn_closure *nnclosure, int start, int end, const double *fluid_d[], double rho[], double p[][6])
{
  for (int j = start; j <= end; j++) {
    rho[j] = fluid_d[j][RHO];
    p[j][0] = fluid_d[j][P11] - ((fluid_d[j][MX] * fluid_d[j][MX]) / fluid_d[j][RHO]);
    p[j][1] = fluid_d[j][P12] - ((fluid_d[j][MX] * fluid_d[j][MY]) / fluid_d[j][RHO]);
    p[j][2] = fluid_d[j][P13] - ((fluid_d[j][MX] * fluid_d[j][MZ]) / fluid_d[j][RHO]);
    p[j][3] = fluid_d[j][P22] - ((fluid_d[j][MY] * fluid_d[j][MY]) / fluid_d[j][RHO]);
    p[j][4] = fluid_d[j][P23] - ((fluid_d[j][MY] * fluid_d[j][MZ]) / fluid_d[j][RHO]);
    p[j][5] = fluid_d[j][P33] - ((fluid_d[j][MZ] * fluid_d[j][MZ]) / fluid_d[j][RHO]);
  }
}

static void
calc_mag_heat_flux(const gkyl_ten_moment_nn_closure *nnclosure, const double *fluid_d[], const double *em_tot_d[], double *heat_flux_d)
{
  const int ndim = nnclosure->ndim;
  const int poly_order = nnclosure->poly_order;
  kann_t* ann = nnclosure->ann;

  double rho_avg = 0.0;
  double drho_dx = 0.0, drho_dy = 0.0, drho_dz = 0.0;
  double drho_dx_dx = 0.0;
  double drho_dx_dy = 0.0;

  double p_avg[6] = { 0.0 };
  double dp_dx[6] = { 0.0 };
  double dp_dy[6] = { 0.0 };
  double dp_dz[6] = { 0.0 };
  double dp_dx_dx[6] = { 0.0 };
  double dp_dx_dy[6] = { 0.0 };

  double B_avg[3] = { 0.0 };
  double dB_dx[3] = { 0.0 };
  double dB_dy[3] = { 0.0 };
  double dB_dz[3] = { 0.0 };
  double dB_dx_dx[3] = { 0.0 };
  double dB_dx_dy[3] = { 0.0 };

  double Q[10] = { 0.0 };

  float *input_data;
  float *output_data;
  if (ndim == 1) {
    if (poly_order == 1) {
      input_data = gkyl_malloc(sizeof(float[6]));
      output_data = gkyl_malloc(sizeof(float[4]));
    }
    else if (poly_order == 2) {
      input_data = gkyl_malloc(sizeof(float[9]));
      output_data = gkyl_malloc(sizeof(float[6]));
    }
  }
  else if (ndim == 2) {
    if (poly_order == 1) {
      input_data = gkyl_malloc(sizeof(float[12]));
      output_data = gkyl_malloc(sizeof(float[8]));
    }
  }

  if (ndim == 1) {
    if (poly_order == 1) {
      const double dx = nnclosure->grid.dx[0];
      double rho[2] = { 0.0 };
      double p[2][6] = { 0.0 };
      var_setup(nnclosure, L_1D, U_1D, fluid_d, rho, p);

      rho_avg = calc_arithm_avg_1D(rho[L_1D], rho[U_1D]);
      p_avg[0] = calc_arithm_avg_1D(p[L_1D][0], p[U_1D][0]);
      p_avg[1] = calc_arithm_avg_1D(p[L_1D][1], p[U_1D][1]);
      p_avg[2] = calc_arithm_avg_1D(p[L_1D][2], p[U_1D][2]);
      p_avg[3] = calc_arithm_avg_1D(p[L_1D][3], p[U_1D][3]);
      p_avg[4] = calc_arithm_avg_1D(p[L_1D][4], p[U_1D][4]);
      p_avg[5] = calc_arithm_avg_1D(p[L_1D][5], p[U_1D][5]);
      B_avg[0] = calc_arithm_avg_1D(em_tot_d[L_1D][BX], em_tot_d[U_1D][BX]);
      B_avg[1] = calc_arithm_avg_1D(em_tot_d[L_1D][BY], em_tot_d[U_1D][BY]);
      B_avg[2] = calc_arithm_avg_1D(em_tot_d[L_1D][BZ], em_tot_d[U_1D][BZ]);
      dB_dx[0] = calc_sym_grad_1D(dx, em_tot_d[L_1D][BX], em_tot_d[U_1D][BX]);
      dB_dx[1] = calc_sym_grad_1D(dx, em_tot_d[L_1D][BY], em_tot_d[U_1D][BY]);
      dB_dx[2] = calc_sym_grad_1D(dx, em_tot_d[L_1D][BZ], em_tot_d[U_1D][BZ]);

      drho_dx = calc_sym_grad_1D(dx, rho[L_1D], rho[U_1D]);
      dp_dx[0] = calc_sym_grad_1D(dx, p[L_1D][0], p[U_1D][0]);
      dp_dx[1] = calc_sym_grad_1D(dx, p[L_1D][1], p[U_1D][1]);
      dp_dx[2] = calc_sym_grad_1D(dx, p[L_1D][2], p[U_1D][2]);
      dp_dx[3] = calc_sym_grad_1D(dx, p[L_1D][3], p[U_1D][3]);
      dp_dx[4] = calc_sym_grad_1D(dx, p[L_1D][4], p[U_1D][4]);
      dp_dx[5] = calc_sym_grad_1D(dx, p[L_1D][5], p[U_1D][5]);

      if (fabs(B_avg[0]) < pow(10.0, -8.0) && fabs(B_avg[1]) < pow(10.0, -8.0) && fabs(B_avg[2]) < pow(10.0, -8.0)) {
        B_avg[0] = 1.0;
        dB_dx[0] = 0.0; dB_dx[1] = 0.0; dB_dx[2] = 0.0;
        dB_dy[0] = 0.0; dB_dy[1] = 0.0; dB_dy[2] = 0.0;
        dB_dz[0] = 0.0; dB_dz[1] = 0.0; dB_dz[2] = 0.0;
      }

      double b_mag = sqrt((B_avg[0] * B_avg[0]) + (B_avg[1] * B_avg[1]) + (B_avg[2] * B_avg[2]));
      double b_mag_dx = ((B_avg[0] * dB_dx[0]) + (B_avg[1] * dB_dx[1]) + (B_avg[2] * dB_dx[2])) / b_mag;

      double local_mag[3];
      for (int i = 0; i < 3; i++) {
        local_mag[i] = B_avg[i] / b_mag;
      }

      double local_mag_dx[3];
      for (int i = 0; i < 3; i++) {
        local_mag_dx[i] = ((b_mag * dB_dx[i]) - (B_avg[i] * b_mag_dx)) / (b_mag * b_mag);
      }

      double p_tensor[3][3];
      double p_tensor_dx[3][3];

      p_tensor[0][0] = p_avg[0]; p_tensor[0][1] = p_avg[1]; p_tensor[0][2] = p_avg[2];
      p_tensor[1][0] = p_avg[1]; p_tensor[1][1] = p_avg[3]; p_tensor[1][2] = p_avg[4];
      p_tensor[2][0] = p_avg[2]; p_tensor[2][1] = p_avg[4]; p_tensor[2][2] = p_avg[5];

      p_tensor_dx[0][0] = dp_dx[0]; p_tensor_dx[0][1] = dp_dx[1]; p_tensor_dx[0][2] = dp_dx[2];
      p_tensor_dx[1][0] = dp_dx[1]; p_tensor_dx[1][1] = dp_dx[3]; p_tensor_dx[1][2] = dp_dx[4];
      p_tensor_dx[2][0] = dp_dx[2]; p_tensor_dx[2][1] = dp_dx[4]; p_tensor_dx[2][2] = dp_dx[5];

      double p_par = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          p_par += p_tensor[i][j] * local_mag[i] * local_mag[j];
        }
      }

      double p_tensor_trace = 0.0;
      for (int i = 0; i < 3; i++) {
        p_tensor_trace += p_tensor[i][i];
      }

      double p_perp = 0.5 * (p_tensor_trace - p_par);

      double p_par_dx = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          p_par_dx += local_mag[i] * local_mag[j] * p_tensor_dx[i][j];
          p_par_dx += p_tensor[i][j] * local_mag[j] * local_mag_dx[i];
          p_par_dx += p_tensor[i][j] * local_mag[i] * local_mag_dx[j];
        }
      }

      double p_tensor_trace_dx = 0.0;
      for (int i = 0; i < 3; i++) {
        p_tensor_trace_dx += p_tensor_dx[i][i];
      }

      double p_perp_dx = 0.5 * (p_tensor_trace_dx - p_par_dx);

      input_data[0] = rho_avg;
      input_data[1] = drho_dx;
      input_data[2] = p_par;
      input_data[3] = p_par_dx;
      input_data[4] = p_perp;
      input_data[5] = p_perp_dx;

      const float *output_data_predicted = kann_apply1(ann, input_data);

      for (int i = 0; i < 4; i++) {
        output_data[i] = output_data_predicted[i];
      }

      double q_par = output_data[0];
      double q_par_dx = output_data[1];
      double q_perp = output_data[2];
      double q_perp_dx = output_data[3];

      double heat_flux_tensor[3][3][3];

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            heat_flux_tensor[i][j][k] = q_par * local_mag[i] * local_mag[j] * local_mag[k];

            if (i == j) {
              heat_flux_tensor[i][j][k] += q_perp * (1.0 - (local_mag[i] * local_mag[j])) * local_mag[k];
            }
            else {
              heat_flux_tensor[i][j][k] -= q_perp * (local_mag[i] * local_mag[j]) * local_mag[k];
            }

            if (i == k) {
              heat_flux_tensor[i][j][k] += q_perp * (1.0 - (local_mag[i] * local_mag[k])) * local_mag[j];
            }
            else {
              heat_flux_tensor[i][j][k] -= q_perp * (local_mag[i] * local_mag[k]) * local_mag[j];
            }

            if (j == k) {
              heat_flux_tensor[i][j][k] += q_perp * (1.0 - (local_mag[j] * local_mag[k])) * local_mag[i];
            }
            else {
              heat_flux_tensor[i][j][k] -= q_perp * (local_mag[j] * local_mag[k]) * local_mag[i];
            }
          }
        }
      }

      Q[0] = B_avg[0] * heat_flux_tensor[0][0][0];
      Q[1] = B_avg[0] * heat_flux_tensor[0][0][1];
      Q[2] = B_avg[0] * heat_flux_tensor[0][0][2];
      Q[3] = B_avg[0] * heat_flux_tensor[0][1][1];
      Q[4] = B_avg[0] * heat_flux_tensor[0][1][2];
      Q[5] = B_avg[0] * heat_flux_tensor[0][2][2];
      Q[6] = B_avg[1] * heat_flux_tensor[1][1][1];
      Q[7] = B_avg[1] * heat_flux_tensor[1][1][2];
      Q[8] = B_avg[1] * heat_flux_tensor[1][2][2];
      Q[9] = B_avg[2] * heat_flux_tensor[2][2][2];
    }
    else if (poly_order == 2) {
      const double dx = nnclosure->grid.dx[0];
      double rho[3] = { 0.0 };
      double p[3][6] = { 0.0 };
      var_setup(nnclosure, L2_1D, U2_1D, fluid_d, rho, p);

      rho_avg = calc_arithm_avg_1D(rho[L2_1D], rho[U2_1D]);
      p_avg[0] = calc_arithm_avg_1D(p[L2_1D][0], p[U2_1D][0]);
      p_avg[1] = calc_arithm_avg_1D(p[L2_1D][1], p[U2_1D][1]);
      p_avg[2] = calc_arithm_avg_1D(p[L2_1D][2], p[U2_1D][2]);
      p_avg[3] = calc_arithm_avg_1D(p[L2_1D][3], p[U2_1D][3]);
      p_avg[4] = calc_arithm_avg_1D(p[L2_1D][4], p[U2_1D][4]);
      p_avg[5] = calc_arithm_avg_1D(p[L2_1D][5], p[U2_1D][5]);
      B_avg[0] = calc_arithm_avg_1D(em_tot_d[L2_1D][BX], em_tot_d[U2_1D][BX]);
      B_avg[1] = calc_arithm_avg_1D(em_tot_d[L2_1D][BY], em_tot_d[U2_1D][BY]);
      B_avg[2] = calc_arithm_avg_1D(em_tot_d[L2_1D][BZ], em_tot_d[U2_1D][BZ]);

      drho_dx = calc_sym_grad_1D(dx, rho[L2_1D], rho[U2_1D]);
      dp_dx[0] = calc_sym_grad_1D(dx, p[L2_1D][0], p[U2_1D][0]);
      dp_dx[1] = calc_sym_grad_1D(dx, p[L2_1D][1], p[U2_1D][1]);
      dp_dx[2] = calc_sym_grad_1D(dx, p[L2_1D][2], p[U2_1D][2]);
      dp_dx[3] = calc_sym_grad_1D(dx, p[L2_1D][3], p[U2_1D][3]);
      dp_dx[4] = calc_sym_grad_1D(dx, p[L2_1D][4], p[U2_1D][4]);
      dp_dx[5] = calc_sym_grad_1D(dx, p[L2_1D][5], p[U2_1D][5]);
      dB_dx[0] = calc_sym_grad_1D(dx, em_tot_d[L2_1D][BX], em_tot_d[U2_1D][BX]);
      dB_dx[1] = calc_sym_grad_1D(dx, em_tot_d[L2_1D][BY], em_tot_d[U2_1D][BY]);
      dB_dx[2] = calc_sym_grad_1D(dx, em_tot_d[L2_1D][BZ], em_tot_d[U2_1D][BZ]);

      drho_dx_dx = calc_sym_grad2_1D(dx, rho[L2_1D], rho[C2_1D], rho[U2_1D]);
      dp_dx_dx[0] = calc_sym_grad2_1D(dx, p[L2_1D][0], p[C2_1D][0], p[U2_1D][0]);
      dp_dx_dx[1] = calc_sym_grad2_1D(dx, p[L2_1D][1], p[C2_1D][1], p[U2_1D][1]);
      dp_dx_dx[2] = calc_sym_grad2_1D(dx, p[L2_1D][2], p[C2_1D][2], p[U2_1D][2]);
      dp_dx_dx[3] = calc_sym_grad2_1D(dx, p[L2_1D][3], p[C2_1D][3], p[U2_1D][3]);
      dp_dx_dx[4] = calc_sym_grad2_1D(dx, p[L2_1D][4], p[C2_1D][4], p[U2_1D][4]);
      dp_dx_dx[5] = calc_sym_grad2_1D(dx, p[L2_1D][5], p[C2_1D][5], p[U2_1D][5]);
      dB_dx_dx[0] = calc_sym_grad2_1D(dx, em_tot_d[L2_1D][BX], em_tot_d[C2_1D][BX], em_tot_d[U2_1D][BX]);
      dB_dx_dx[1] = calc_sym_grad2_1D(dx, em_tot_d[L2_1D][BY], em_tot_d[C2_1D][BY], em_tot_d[U2_1D][BY]);
      dB_dx_dx[2] = calc_sym_grad2_1D(dx, em_tot_d[L2_1D][BZ], em_tot_d[C2_1D][BZ], em_tot_d[U2_1D][BZ]);

      if (fabs(B_avg[0]) < pow(10.0, -8.0) && fabs(B_avg[1]) < pow(10.0, -8.0) && fabs(B_avg[2]) < pow(10.0, -8.0)) {
        B_avg[0] = 1.0;
        dB_dx[0] = 0.0; dB_dx[1] = 0.0; dB_dx[2] = 0.0;
        dB_dy[0] = 0.0; dB_dy[1] = 0.0; dB_dy[2] = 0.0;
        dB_dz[0] = 0.0; dB_dz[1] = 0.0; dB_dz[2] = 0.0;
        dB_dx_dx[0] = 0.0; dB_dx_dx[1] = 0.0; dB_dx_dx[2] = 0.0;
      }

      double b_mag = sqrt((B_avg[0] * B_avg[0]) + (B_avg[1] * B_avg[1]) + (B_avg[2] * B_avg[2]));
      double b_mag_dx = ((B_avg[0] * dB_dx[0]) + (B_avg[1] * dB_dx[1]) + (B_avg[2] * dB_dx[2])) / b_mag;
      double b_mag_dx_dx = b_mag * ((dB_dx[0] * dB_dx[0]) + (B_avg[0] * dB_dx_dx[0]) + (dB_dx[1] * dB_dx[1]) + (B_avg[1] * dB_dx_dx[1]) + (dB_dx[2] * dB_dy[2]) + (B_avg[2] * dB_dx_dy[2]));
      b_mag_dx_dx -= b_mag_dx * b_mag_dx;
      b_mag_dx_dx /= b_mag * b_mag;

      double local_mag[3];
      for (int i = 0; i < 3; i++) {
        local_mag[i] = B_avg[i] / b_mag;
      }

      double local_mag_dx[3];
      double local_mag_dy[3];
      double local_mag_dx_dx[3];
      for (int i = 0; i < 3; i++) {
        local_mag_dx[i] = ((b_mag * dB_dx[i]) - (B_avg[i] * b_mag_dx)) / (b_mag * b_mag);
        local_mag_dx_dx[i] = (b_mag * b_mag) * ((b_mag_dx * dB_dx[i]) + (b_mag * dB_dx_dx[i]) - (dB_dx[i] * b_mag_dx) * (B_avg[i] * b_mag_dx_dx));
        local_mag_dx_dx[i] -= 2.0 * ((b_mag * dB_dx[i]) - (B_avg[i] * b_mag_dx)) * b_mag * b_mag_dx;
        local_mag_dx_dx[i] /= b_mag * b_mag * b_mag * b_mag;
      }

      double p_tensor[3][3];
      double p_tensor_dx[3][3];
      double p_tensor_dx_dx[3][3];

      p_tensor[0][0] = p_avg[0]; p_tensor[0][1] = p_avg[1]; p_tensor[0][2] = p_avg[2];
      p_tensor[1][0] = p_avg[1]; p_tensor[1][1] = p_avg[3]; p_tensor[1][2] = p_avg[4];
      p_tensor[2][0] = p_avg[2]; p_tensor[2][1] = p_avg[4]; p_tensor[2][2] = p_avg[5];

      p_tensor_dx[0][0] = dp_dx[0]; p_tensor_dx[0][1] = dp_dx[1]; p_tensor_dx[0][2] = dp_dx[2];
      p_tensor_dx[1][0] = dp_dx[1]; p_tensor_dx[1][1] = dp_dx[3]; p_tensor_dx[1][2] = dp_dx[4];
      p_tensor_dx[2][0] = dp_dx[2]; p_tensor_dx[2][1] = dp_dx[4]; p_tensor_dx[2][2] = dp_dx[5];

      p_tensor_dx_dx[0][0] = dp_dx_dx[0]; p_tensor_dx_dx[0][1] = dp_dx_dx[1]; p_tensor_dx_dx[0][2] = dp_dx_dx[2];
      p_tensor_dx_dx[1][0] = dp_dx_dx[1]; p_tensor_dx_dx[1][1] = dp_dx_dx[3]; p_tensor_dx_dx[1][2] = dp_dx_dx[4];
      p_tensor_dx_dx[2][0] = dp_dx_dx[2]; p_tensor_dx_dx[2][1] = dp_dx_dx[4]; p_tensor_dx_dx[2][2] = dp_dx_dx[5];

      double p_par = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          p_par += p_tensor[i][j] * local_mag[i] * local_mag[j];
        }
      }

      double p_tensor_trace = 0.0;
      for (int i = 0; i < 3; i++) {
        p_tensor_trace += p_tensor[i][i];
      }

      double p_perp = 0.5 * (p_tensor_trace - p_par);

      double p_par_dx = 0.0;
      double p_par_dx_dx = 0.0;

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          p_par_dx += local_mag[i] * local_mag[j] * p_tensor_dx[i][j];
          p_par_dx += p_tensor[i][j] * local_mag[j] * local_mag_dx[i];
          p_par_dx += p_tensor[i][j] * local_mag[i] * local_mag_dx[j];

          p_par_dx_dx += local_mag[i] * local_mag[i] * p_tensor_dx_dx[i][j];
          p_par_dx_dx += p_tensor_dx[i][j] * local_mag[i] * local_mag_dx[j];
          p_par_dx_dx += p_tensor_dx[i][j] * local_mag[j] * local_mag_dx[i];
          p_par_dx_dx += p_tensor[i][j] * local_mag[i] * local_mag_dx_dx[j];
          p_par_dx_dx += local_mag_dx[j] * local_mag[i] * p_tensor_dx[i][j];
          p_par_dx_dx += local_mag_dx[j] * p_tensor[i][j] * local_mag_dx[i];
          p_par_dx_dx += p_tensor[i][j] * local_mag[j] * local_mag_dx_dx[i];
          p_par_dx_dx += p_tensor[i][j] * local_mag_dx[j] * local_mag_dx[i];
          p_par_dx_dx += local_mag[j] * local_mag_dx[i] * p_tensor_dx[i][j];
        }
      }

      double p_tensor_trace_dx = 0.0;
      double p_tensor_trace_dx_dx = 0.0;
      for (int i = 0; i < 3; i++) {
        p_tensor_trace_dx += p_tensor_dx[i][i];
        p_tensor_trace_dx_dx += p_tensor_dx_dx[i][i];
      }

      double p_perp_dx = 0.5 * (p_tensor_trace_dx - p_par_dx);
      double p_perp_dx_dx = 0.5 * (p_tensor_trace_dx_dx - p_par_dx_dx);

      input_data[0] = rho_avg;
      input_data[1] = drho_dx;
      input_data[2] = drho_dx_dx;
      input_data[3] = p_par;
      input_data[4] = p_par_dx;
      input_data[5] = p_par_dx_dx;
      input_data[6] = p_perp;
      input_data[7] = p_perp_dx;
      input_data[8] = p_perp_dx_dx;

      const float *output_data_predicted = kann_apply1(ann, input_data);

      for (int i = 0; i < 6; i++) {
        output_data[i] = output_data_predicted[i];
      }

      double q_par = output_data[0];
      double q_par_dx = output_data[1];
      double q_perp = output_data[3];
      double q_perp_dx = output_data[4];

      double heat_flux_tensor[3][3][3];

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            heat_flux_tensor[i][j][k] = q_par * local_mag[i] * local_mag[j] * local_mag[k];

            if (i == j) {
              heat_flux_tensor[i][j][k] += q_perp * (1.0 - (local_mag[i] * local_mag[j])) * local_mag[k];
            }
            else {
              heat_flux_tensor[i][j][k] -= q_perp * (local_mag[i] * local_mag[j]) * local_mag[k];
            }

            if (i == k) {
              heat_flux_tensor[i][j][k] += q_perp * (1.0 - (local_mag[i] * local_mag[k])) * local_mag[j];
            }
            else {
              heat_flux_tensor[i][j][k] -= q_perp * (local_mag[i] * local_mag[k]) * local_mag[j];
            }

            if (j == k) {
              heat_flux_tensor[i][j][k] += q_perp * (1.0 - (local_mag[j] * local_mag[k])) * local_mag[i];
            }
            else {
              heat_flux_tensor[i][j][k] -= q_perp * (local_mag[j] * local_mag[k]) * local_mag[i];
            }
          }
        }
      }

      Q[0] = B_avg[0] * heat_flux_tensor[0][0][0];
      Q[1] = B_avg[0] * heat_flux_tensor[0][0][1];
      Q[2] = B_avg[0] * heat_flux_tensor[0][0][2];
      Q[3] = B_avg[0] * heat_flux_tensor[0][1][1];
      Q[4] = B_avg[0] * heat_flux_tensor[0][1][2];
      Q[5] = B_avg[0] * heat_flux_tensor[0][2][2];
      Q[6] = B_avg[1] * heat_flux_tensor[1][1][1];
      Q[7] = B_avg[1] * heat_flux_tensor[1][1][2];
      Q[8] = B_avg[1] * heat_flux_tensor[1][2][2];
      Q[9] = B_avg[2] * heat_flux_tensor[2][2][2];
    }
  }
  else if (ndim == 2) {
    if (poly_order == 1) {
      const double dx = nnclosure->grid.dx[0];
      const double dy = nnclosure->grid.dx[1];
      double rho[4] = { 0.0 };
      double p[4][6] = { 0.0 };
      var_setup(nnclosure, LL_2D, UU_2D, fluid_d, rho, p);
  
      rho_avg = calc_arithm_avg_2D(rho[LL_2D], rho[LU_2D], rho[UL_2D], rho[UU_2D]);
      p_avg[0] = calc_arithm_avg_2D(p[LL_2D][0], p[LU_2D][0], p[UL_2D][0], p[UU_2D][0]);
      p_avg[1] = calc_arithm_avg_2D(p[LL_2D][1], p[LU_2D][1], p[UL_2D][1], p[UU_2D][1]);
      p_avg[2] = calc_arithm_avg_2D(p[LL_2D][2], p[LU_2D][2], p[UL_2D][2], p[UU_2D][2]);
      p_avg[3] = calc_arithm_avg_2D(p[LL_2D][3], p[LU_2D][3], p[UL_2D][3], p[UU_2D][3]);
      p_avg[4] = calc_arithm_avg_2D(p[LL_2D][4], p[LU_2D][4], p[UL_2D][4], p[UU_2D][4]);
      p_avg[5] = calc_arithm_avg_2D(p[LL_2D][5], p[LU_2D][5], p[UL_2D][5], p[UU_2D][5]);
      B_avg[0] = calc_arithm_avg_2D(em_tot_d[LL_2D][BX], em_tot_d[LU_2D][BX], em_tot_d[UL_2D][BX], em_tot_d[UU_2D][BX]);
      B_avg[1] = calc_arithm_avg_2D(em_tot_d[LL_2D][BY], em_tot_d[LU_2D][BY], em_tot_d[UL_2D][BY], em_tot_d[UU_2D][BY]);
      B_avg[2] = calc_arithm_avg_2D(em_tot_d[LL_2D][BZ], em_tot_d[LU_2D][BZ], em_tot_d[UL_2D][BZ], em_tot_d[UU_2D][BZ]);

      drho_dx = calc_sym_gradx_2D(dx, rho[LL_2D], rho[LU_2D], rho[UL_2D], rho[UU_2D]);
      dp_dx[0] = calc_sym_gradx_2D(dx, p[LL_2D][0], p[LU_2D][0], p[UL_2D][0], p[UU_2D][0]);
      dp_dx[1] = calc_sym_gradx_2D(dx, p[LL_2D][1], p[LU_2D][1], p[UL_2D][1], p[UU_2D][1]);
      dp_dx[2] = calc_sym_gradx_2D(dx, p[LL_2D][2], p[LU_2D][2], p[UL_2D][2], p[UU_2D][2]);
      dp_dx[3] = calc_sym_gradx_2D(dx, p[LL_2D][3], p[LU_2D][3], p[UL_2D][3], p[UU_2D][3]);
      dp_dx[4] = calc_sym_gradx_2D(dx, p[LL_2D][4], p[LU_2D][4], p[UL_2D][4], p[UU_2D][4]);
      dp_dx[5] = calc_sym_gradx_2D(dx, p[LL_2D][5], p[LU_2D][5], p[UL_2D][5], p[UU_2D][5]);
      dB_dx[0] = calc_sym_gradx_2D(dx, em_tot_d[LL_2D][BX], em_tot_d[LU_2D][BX], em_tot_d[UL_2D][BX], em_tot_d[UU_2D][BX]);
      dB_dx[1] = calc_sym_gradx_2D(dx, em_tot_d[LL_2D][BY], em_tot_d[LU_2D][BY], em_tot_d[UL_2D][BY], em_tot_d[UU_2D][BY]);
      dB_dx[2] = calc_sym_gradx_2D(dx, em_tot_d[LL_2D][BZ], em_tot_d[LU_2D][BZ], em_tot_d[UL_2D][BZ], em_tot_d[UU_2D][BZ]);

      drho_dy = calc_sym_grady_2D(dy, rho[LL_2D], rho[LU_2D], rho[UL_2D], rho[UU_2D]);
      dp_dy[0] = calc_sym_grady_2D(dy, p[LL_2D][0], p[LU_2D][0], p[UL_2D][0], p[UU_2D][0]);
      dp_dy[1] = calc_sym_grady_2D(dy, p[LL_2D][1], p[LU_2D][1], p[UL_2D][1], p[UU_2D][1]);
      dp_dy[2] = calc_sym_grady_2D(dy, p[LL_2D][2], p[LU_2D][2], p[UL_2D][2], p[UU_2D][2]);
      dp_dy[3] = calc_sym_grady_2D(dy, p[LL_2D][3], p[LU_2D][3], p[UL_2D][3], p[UU_2D][3]);
      dp_dy[4] = calc_sym_grady_2D(dy, p[LL_2D][4], p[LU_2D][4], p[UL_2D][4], p[UU_2D][4]);
      dp_dy[5] = calc_sym_grady_2D(dy, p[LL_2D][5], p[LU_2D][5], p[UL_2D][5], p[UU_2D][5]);
      dB_dy[0] = calc_sym_grady_2D(dy, em_tot_d[LL_2D][BX], em_tot_d[LU_2D][BX], em_tot_d[UL_2D][BX], em_tot_d[UU_2D][BX]);
      dB_dy[1] = calc_sym_grady_2D(dy, em_tot_d[LL_2D][BY], em_tot_d[LU_2D][BY], em_tot_d[UL_2D][BY], em_tot_d[UU_2D][BY]);
      dB_dy[2] = calc_sym_grady_2D(dy, em_tot_d[LL_2D][BZ], em_tot_d[LU_2D][BZ], em_tot_d[UL_2D][BZ], em_tot_d[UU_2D][BZ]);

      drho_dx_dy = calc_sym_gradxy_2D(dx, dy, rho[LL_2D], rho[LU_2D], rho[UL_2D], rho[UU_2D]);
      dp_dx_dy[0] = calc_sym_gradxy_2D(dx, dy, p[LL_2D][0], p[LU_2D][0], p[UL_2D][0], p[UU_2D][0]);
      dp_dx_dy[1] = calc_sym_gradxy_2D(dx, dy, p[LL_2D][1], p[LU_2D][1], p[UL_2D][1], p[UU_2D][1]);
      dp_dx_dy[2] = calc_sym_gradxy_2D(dx, dy, p[LL_2D][2], p[LU_2D][2], p[UL_2D][2], p[UU_2D][2]);
      dp_dx_dy[3] = calc_sym_gradxy_2D(dx, dy, p[LL_2D][3], p[LU_2D][3], p[UL_2D][3], p[UU_2D][3]);
      dp_dx_dy[4] = calc_sym_gradxy_2D(dx, dy, p[LL_2D][4], p[LU_2D][4], p[UL_2D][4], p[UU_2D][4]);
      dp_dx_dy[5] = calc_sym_gradxy_2D(dx, dy, p[LL_2D][5], p[LU_2D][5], p[UL_2D][5], p[UU_2D][5]);
      dB_dx_dy[0] = calc_sym_gradxy_2D(dx, dy, em_tot_d[LL_2D][BX], em_tot_d[LU_2D][BX], em_tot_d[UL_2D][BX], em_tot_d[UU_2D][BX]);
      dB_dx_dy[1] = calc_sym_gradxy_2D(dx, dy, em_tot_d[LL_2D][BY], em_tot_d[LU_2D][BY], em_tot_d[UL_2D][BY], em_tot_d[UU_2D][BY]);
      dB_dx_dy[2] = calc_sym_gradxy_2D(dx, dy, em_tot_d[LL_2D][BZ], em_tot_d[LU_2D][BZ], em_tot_d[UL_2D][BZ], em_tot_d[UU_2D][BZ]);

      if (fabs(B_avg[0]) < pow(10.0, -8.0) && fabs(B_avg[1]) < pow(10.0, -8.0) && fabs(B_avg[2]) < pow(10.0, -8.0)) {
        B_avg[0] = 1.0;
        dB_dx[0] = 0.0; dB_dx[1] = 0.0; dB_dx[2] = 0.0;
        dB_dy[0] = 0.0; dB_dy[1] = 0.0; dB_dy[2] = 0.0;
        dB_dz[0] = 0.0; dB_dz[1] = 0.0; dB_dz[2] = 0.0;
        dB_dx_dy[0] = 0.0; dB_dx_dy[1] = 0.0; dB_dx_dy[2] = 0.0;
      }

      double b_mag = sqrt((B_avg[0] * B_avg[0]) + (B_avg[1] * B_avg[1]) + (B_avg[2] * B_avg[2]));
      double b_mag_dx = ((B_avg[0] * dB_dx[0]) + (B_avg[1] * dB_dx[1]) + (B_avg[2] * dB_dx[2])) / b_mag;
      double b_mag_dy = ((B_avg[0] * dB_dy[0]) + (B_avg[1] * dB_dy[1]) + (B_avg[2] * dB_dy[2])) / b_mag;
      double b_mag_dx_dy = b_mag * ((dB_dx[0] * dB_dy[0]) + (B_avg[0] * dB_dx_dy[0]) + (dB_dx[1] * dB_dy[1]) + (B_avg[1] * dB_dx_dy[1]) + (dB_dx[2] * dB_dy[2]) + (B_avg[2] * dB_dx_dy[2]));
      b_mag_dx_dy -= b_mag_dx * b_mag_dy;
      b_mag_dx_dy /= b_mag * b_mag;

      double local_mag[3];
      for (int i = 0; i < 3; i++) {
        local_mag[i] = B_avg[i] / b_mag;
      }

      double local_mag_dx[3];
      double local_mag_dy[3];
      double local_mag_dx_dy[3];
      for (int i = 0; i < 3; i++) {
        local_mag_dx[i] = ((b_mag * dB_dx[i]) - (B_avg[i] * b_mag_dx)) / (b_mag * b_mag);
        local_mag_dy[i] = ((b_mag * dB_dy[i]) - (B_avg[i] * b_mag_dy)) / (b_mag * b_mag);
        local_mag_dx_dy[i] = (b_mag * b_mag) * ((b_mag_dy * dB_dx[i]) + (b_mag * dB_dx_dy[i]) - (dB_dy[i] * b_mag_dx) * (B_avg[i] * b_mag_dx_dy));
        local_mag_dx_dy[i] -= 2.0 * ((b_mag * dB_dx[i]) - (B_avg[i] * b_mag_dx)) * b_mag * b_mag_dy;
        local_mag_dx_dy[i] /= b_mag * b_mag * b_mag * b_mag;
      }

      double p_tensor[3][3];
      double p_tensor_dx[3][3];
      double p_tensor_dy[3][3];
      double p_tensor_dx_dy[3][3];

      p_tensor[0][0] = p_avg[0]; p_tensor[0][1] = p_avg[1]; p_tensor[0][2] = p_avg[2];
      p_tensor[1][0] = p_avg[1]; p_tensor[1][1] = p_avg[3]; p_tensor[1][2] = p_avg[4];
      p_tensor[2][0] = p_avg[2]; p_tensor[2][1] = p_avg[4]; p_tensor[2][2] = p_avg[5];

      p_tensor_dx[0][0] = dp_dx[0]; p_tensor_dx[0][1] = dp_dx[1]; p_tensor_dx[0][2] = dp_dx[2];
      p_tensor_dx[1][0] = dp_dx[1]; p_tensor_dx[1][1] = dp_dx[3]; p_tensor_dx[1][2] = dp_dx[4];
      p_tensor_dx[2][0] = dp_dx[2]; p_tensor_dx[2][1] = dp_dx[4]; p_tensor_dx[2][2] = dp_dx[5];

      p_tensor_dy[0][0] = dp_dy[0]; p_tensor_dy[0][1] = dp_dy[1]; p_tensor_dy[0][2] = dp_dy[2];
      p_tensor_dy[1][0] = dp_dy[1]; p_tensor_dy[1][1] = dp_dy[3]; p_tensor_dy[1][2] = dp_dy[4];
      p_tensor_dy[2][0] = dp_dy[2]; p_tensor_dy[2][1] = dp_dy[4]; p_tensor_dy[2][2] = dp_dy[5];

      p_tensor_dx_dy[0][0] = dp_dx_dy[0]; p_tensor_dx_dy[0][1] = dp_dx_dy[1]; p_tensor_dx_dy[0][2] = dp_dx_dy[2];
      p_tensor_dx_dy[1][0] = dp_dx_dy[1]; p_tensor_dx_dy[1][1] = dp_dx_dy[3]; p_tensor_dx_dy[1][2] = dp_dx_dy[4];
      p_tensor_dx_dy[2][0] = dp_dx_dy[2]; p_tensor_dx_dy[2][1] = dp_dx_dy[4]; p_tensor_dx_dy[2][2] = dp_dx_dy[5];

      double p_par = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          p_par += p_tensor[i][j] * local_mag[i] * local_mag[j];
        }
      }

      double p_tensor_trace = 0.0;
      for (int i = 0; i < 3; i++) {
        p_tensor_trace += p_tensor[i][i];
      }

      double p_perp = 0.5 * (p_tensor_trace - p_par);

      double p_par_dx = 0.0;
      double p_par_dy = 0.0;
      double p_par_dx_dy = 0.0;

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          p_par_dx += local_mag[i] * local_mag[j] * p_tensor_dx[i][j];
          p_par_dx += p_tensor[i][j] * local_mag[j] * local_mag_dx[i];
          p_par_dx += p_tensor[i][j] * local_mag[i] * local_mag_dx[j];

          p_par_dy += local_mag[i] * local_mag[j] * p_tensor_dy[i][j];
          p_par_dy += p_tensor[i][j] * local_mag[j] * local_mag_dy[i];
          p_par_dy += p_tensor[i][j] * local_mag[i] * local_mag_dy[j];

          p_par_dx_dy += local_mag[i] * local_mag[i] * p_tensor_dx_dy[i][j];
          p_par_dx_dy += p_tensor_dx[i][j] * local_mag[i] * local_mag_dy[j];
          p_par_dx_dy += p_tensor_dx[i][j] * local_mag[j] * local_mag_dy[i];
          p_par_dx_dy += p_tensor[i][j] * local_mag[i] * local_mag_dx_dy[j];
          p_par_dx_dy += local_mag_dx[j] * local_mag[i] * p_tensor_dy[i][j];
          p_par_dx_dy += local_mag_dx[j] * p_tensor[i][j] * local_mag_dy[i];
          p_par_dx_dy += p_tensor[i][j] * local_mag[j] * local_mag_dx_dy[i];
          p_par_dx_dy += p_tensor[i][j] * local_mag_dy[j] * local_mag_dx[i];
          p_par_dx_dy += local_mag[j] * local_mag_dx[i] * p_tensor_dy[i][j];
        }
      }

      double p_tensor_trace_dx = 0.0;
      double p_tensor_trace_dy = 0.0;
      double p_tensor_trace_dx_dy = 0.0;
      for (int i = 0; i < 3; i++) {
        p_tensor_trace_dx += p_tensor_dx[i][i];
        p_tensor_trace_dy += p_tensor_dy[i][i];
        p_tensor_trace_dx_dy += p_tensor_dx_dy[i][i];
      }

      double p_perp_dx = 0.5 * (p_tensor_trace_dx - p_par_dx);
      double p_perp_dy = 0.5 * (p_tensor_trace_dy - p_par_dy);
      double p_perp_dx_dy = 0.5 * (p_tensor_trace_dx_dy - p_par_dx_dy);

      input_data[0] = rho_avg;
      input_data[1] = drho_dx;
      input_data[2] = drho_dy;
      input_data[3] = drho_dx_dy;
      input_data[4] = p_par;
      input_data[5] = p_par_dx;
      input_data[6] = p_par_dy;
      input_data[7] = p_par_dx_dy;
      input_data[8] = p_perp;
      input_data[9] = p_perp_dx;
      input_data[10] = p_perp_dy;
      input_data[11] = p_perp_dx_dy;

      const float *output_data_predicted = kann_apply1(ann, input_data);

      for (int i = 0; i < 8; i++) {
        output_data[i] = output_data_predicted[i];
      }

      double q_par = output_data[0];
      double q_par_dx = output_data[1];
      double q_par_dy = output_data[2];
      double q_perp = output_data[4];
      double q_perp_dx = output_data[5];
      double q_perp_dy = output_data[6];

      double heat_flux_tensor[3][3][3];

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            heat_flux_tensor[i][j][k] = q_par * local_mag[i] * local_mag[j] * local_mag[k];

            if (i == j) {
              heat_flux_tensor[i][j][k] += q_perp * (1.0 - (local_mag[i] * local_mag[j])) * local_mag[k];
            }
            else {
              heat_flux_tensor[i][j][k] -= q_perp * (local_mag[i] * local_mag[j]) * local_mag[k];
            }

            if (i == k) {
              heat_flux_tensor[i][j][k] += q_perp * (1.0 - (local_mag[i] * local_mag[k])) * local_mag[j];
            }
            else {
              heat_flux_tensor[i][j][k] -= q_perp * (local_mag[i] * local_mag[k]) * local_mag[j];
            }

            if (j == k) {
              heat_flux_tensor[i][j][k] += q_perp * (1.0 - (local_mag[j] * local_mag[k])) * local_mag[i];
            }
            else {
              heat_flux_tensor[i][j][k] -= q_perp * (local_mag[j] * local_mag[k]) * local_mag[i];
            }
          }
        }
      }

      Q[0] = B_avg[0] * heat_flux_tensor[0][0][0];
      Q[1] = B_avg[0] * heat_flux_tensor[0][0][1];
      Q[2] = B_avg[0] * heat_flux_tensor[0][0][2];
      Q[3] = B_avg[0] * heat_flux_tensor[0][1][1];
      Q[4] = B_avg[0] * heat_flux_tensor[0][1][2];
      Q[5] = B_avg[0] * heat_flux_tensor[0][2][2];
      Q[6] = B_avg[1] * heat_flux_tensor[1][1][1];
      Q[7] = B_avg[1] * heat_flux_tensor[1][1][2];
      Q[8] = B_avg[1] * heat_flux_tensor[1][2][2];
      Q[9] = B_avg[2] * heat_flux_tensor[2][2][2];
    }
  }

  double k0 = nnclosure->k0;
  double alpha = 1.0 / k0;

  double p_scalar = (p_avg[0] + p_avg[3] + p_avg[5]) / 3.0;
  double vth_avg = sqrt(p_scalar / rho_avg);

  heat_flux_d[Q111] = alpha * rho_avg * vth_avg * Q[0];
  heat_flux_d[Q112] = alpha * rho_avg * vth_avg * Q[1];
  heat_flux_d[Q113] = alpha * rho_avg * vth_avg * Q[2];
  heat_flux_d[Q122] = alpha * rho_avg * vth_avg * Q[3];
  heat_flux_d[Q123] = alpha * rho_avg * vth_avg * Q[4];
  heat_flux_d[Q133] = alpha * rho_avg * vth_avg * Q[5];
  heat_flux_d[Q222] = alpha * rho_avg * vth_avg * Q[6];
  heat_flux_d[Q223] = alpha * rho_avg * vth_avg * Q[7];
  heat_flux_d[Q233] = alpha * rho_avg * vth_avg * Q[8];
  heat_flux_d[Q333] = alpha * rho_avg * vth_avg * Q[9];

  gkyl_free(input_data);
  gkyl_free(output_data);
}

static void
calc_nn_closure_update(const gkyl_ten_moment_nn_closure *nnclosure, const double *fluid_d[], const double *em_tot_d[], double *rhs)
{
  const int ndim = nnclosure->ndim;
  const int poly_order = nnclosure->poly_order;
  kann_t* ann = nnclosure->ann;

  double rho_avg = 0.0;
  double drho_dx = 0.0, drho_dy = 0.0, drho_dz = 0.0;
  double drho_dx_dx = 0.0;
  double drho_dx_dy = 0.0;

  double p_avg[6] = { 0.0 };
  double dp_dx[6] = { 0.0 };
  double dp_dy[6] = { 0.0 };
  double dp_dz[6] = { 0.0 };
  double dp_dx_dx[6] = { 0.0 };
  double dp_dx_dy[6] = { 0.0 };

  double B_avg[3] = { 0.0 };
  double dB_dx[3] = { 0.0 };
  double dB_dy[3] = { 0.0 };
  double dB_dz[3] = { 0.0 };
  double dB_dx_dx[3] = { 0.0 };
  double dB_dx_dy[3] = { 0.0 };

  double divQx[6] = { 0.0 };
  double divQy[6] = { 0.0 };
  double divQz[6] = { 0.0 };

  float *input_data;
  float *output_data;
  if (ndim == 1) {
    if (poly_order == 1) {
      input_data = gkyl_malloc(sizeof(float[6]));
      output_data = gkyl_malloc(sizeof(float[4]));
    }
    else if (poly_order == 2) {
      input_data = gkyl_malloc(sizeof(float[9]));
      output_data = gkyl_malloc(sizeof(float[6]));
    }
  }
  else if (ndim == 2) {
    if (poly_order == 1) {
      input_data = gkyl_malloc(sizeof(float[12]));
      output_data = gkyl_malloc(sizeof(float[8]));
    }
  }

  if (ndim == 1) {
    if (poly_order == 1 ) {
      const double dx = nnclosure->grid.dx[0];
      double rho[2] = { 0.0 };
      double p[2][6] = { 0.0 };
      var_setup(nnclosure, L_1D, U_1D, fluid_d, rho, p);

      rho_avg = calc_arithm_avg_1D(rho[L_1D], rho[U_1D]);
      p_avg[0] = calc_arithm_avg_1D(p[L_1D][0], p[U_1D][0]);
      p_avg[1] = calc_arithm_avg_1D(p[L_1D][1], p[U_1D][1]);
      p_avg[2] = calc_arithm_avg_1D(p[L_1D][2], p[U_1D][2]);
      p_avg[3] = calc_arithm_avg_1D(p[L_1D][3], p[U_1D][3]);
      p_avg[4] = calc_arithm_avg_1D(p[L_1D][4], p[U_1D][4]);
      p_avg[5] = calc_arithm_avg_1D(p[L_1D][5], p[U_1D][5]);
      B_avg[0] = calc_arithm_avg_1D(em_tot_d[L_1D][BX], em_tot_d[U_1D][BX]);
      B_avg[1] = calc_arithm_avg_1D(em_tot_d[L_1D][BY], em_tot_d[U_1D][BY]);
      B_avg[2] = calc_arithm_avg_1D(em_tot_d[L_1D][BZ], em_tot_d[U_1D][BZ]);

      drho_dx = calc_sym_grad_1D(dx, rho[L_1D], rho[U_1D]);
      dp_dx[0] = calc_sym_grad_1D(dx, p[L_1D][0], p[U_1D][0]);
      dp_dx[1] = calc_sym_grad_1D(dx, p[L_1D][1], p[U_1D][1]);
      dp_dx[2] = calc_sym_grad_1D(dx, p[L_1D][2], p[U_1D][2]);
      dp_dx[3] = calc_sym_grad_1D(dx, p[L_1D][3], p[U_1D][3]);
      dp_dx[4] = calc_sym_grad_1D(dx, p[L_1D][4], p[U_1D][4]);
      dp_dx[5] = calc_sym_grad_1D(dx, p[L_1D][5], p[U_1D][5]);
      dB_dx[0] = calc_sym_grad_1D(dx, em_tot_d[L_1D][BX], em_tot_d[U_1D][BX]);
      dB_dx[1] = calc_sym_grad_1D(dx, em_tot_d[L_1D][BY], em_tot_d[U_1D][BY]);
      dB_dx[2] = calc_sym_grad_1D(dx, em_tot_d[L_1D][BZ], em_tot_d[U_1D][BZ]);

      if (fabs(B_avg[0]) < pow(10.0, -8.0) && fabs(B_avg[1]) < pow(10.0, -8.0) && fabs(B_avg[2]) < pow(10.0, -8.0)) {
        B_avg[0] = 1.0;
        dB_dx[0] = 0.0; dB_dx[1] = 0.0; dB_dx[2] = 0.0;
        dB_dy[0] = 0.0; dB_dy[1] = 0.0; dB_dy[2] = 0.0;
        dB_dz[0] = 0.0; dB_dz[1] = 0.0; dB_dz[2] = 0.0;
      }

      double b_mag = sqrt((B_avg[0] * B_avg[0]) + (B_avg[1] * B_avg[1]) + (B_avg[2] * B_avg[2]));
      double b_mag_dx = ((B_avg[0] * dB_dx[0]) + (B_avg[1] * dB_dx[1]) + (B_avg[2] * dB_dx[2])) / b_mag;

      double local_mag[3];
      for (int i = 0; i < 3; i++) {
        local_mag[i] = B_avg[i] / b_mag;
      }

      double local_mag_dx[3];
      for (int i = 0; i < 3; i++) {
        local_mag_dx[i] = ((b_mag * dB_dx[i]) - (B_avg[i] * b_mag_dx)) / (b_mag * b_mag);
      }

      double p_tensor[3][3];
      double p_tensor_dx[3][3];

      p_tensor[0][0] = p_avg[0]; p_tensor[0][1] = p_avg[1]; p_tensor[0][2] = p_avg[2];
      p_tensor[1][0] = p_avg[1]; p_tensor[1][1] = p_avg[3]; p_tensor[1][2] = p_avg[4];
      p_tensor[2][0] = p_avg[2]; p_tensor[2][1] = p_avg[4]; p_tensor[2][2] = p_avg[5];

      p_tensor_dx[0][0] = dp_dx[0]; p_tensor_dx[0][1] = dp_dx[1]; p_tensor_dx[0][2] = dp_dx[2];
      p_tensor_dx[1][0] = dp_dx[1]; p_tensor_dx[1][1] = dp_dx[3]; p_tensor_dx[1][2] = dp_dx[4];
      p_tensor_dx[2][0] = dp_dx[2]; p_tensor_dx[2][1] = dp_dx[4]; p_tensor_dx[2][2] = dp_dx[5];

      double p_par = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          p_par += p_tensor[i][j] * local_mag[i] * local_mag[j];
        }
      }

      double p_tensor_trace = 0.0;
      for (int i = 0; i < 3; i++) {
        p_tensor_trace += p_tensor[i][i];
      }

      double p_perp = 0.5 * (p_tensor_trace - p_par);

      double p_par_dx = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          p_par_dx += local_mag[i] * local_mag[j] * p_tensor_dx[i][j];
          p_par_dx += p_tensor[i][j] * local_mag[j] * local_mag_dx[i];
          p_par_dx += p_tensor[i][j] * local_mag[i] * local_mag_dx[j];
        }
      }

      double p_tensor_trace_dx = 0.0;
      for (int i = 0; i < 3; i++) {
        p_tensor_trace_dx += p_tensor_dx[i][i];
      }

      double p_perp_dx = 0.5 * (p_tensor_trace_dx - p_par_dx);

      input_data[0] = rho_avg;
      input_data[1] = drho_dx;
      input_data[2] = p_par;
      input_data[3] = p_par_dx;
      input_data[4] = p_perp;
      input_data[5] = p_perp_dx;

      const float *output_data_predicted = kann_apply1(ann, input_data);

      for (int i = 0; i < 4; i++) {
        output_data[i] = output_data_predicted[i];
      }

      double q_par = output_data[0];
      double q_par_dx = output_data[1];
      double q_perp = output_data[2];
      double q_perp_dx = output_data[3];

      double heat_flux_tensor[3][3][3];

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            heat_flux_tensor[i][j][k] = q_par * local_mag[i] * local_mag[j] * local_mag[k];

            if (i == j) {
              heat_flux_tensor[i][j][k] += q_perp * (1.0 - (local_mag[i] * local_mag[j])) * local_mag[k];
            }
            else {
              heat_flux_tensor[i][j][k] -= q_perp * (local_mag[i] * local_mag[j]) * local_mag[k];
            }

            if (i == k) {
              heat_flux_tensor[i][j][k] += q_perp * (1.0 - (local_mag[i] * local_mag[k])) * local_mag[j];
            }
            else {
              heat_flux_tensor[i][j][k] -= q_perp * (local_mag[i] * local_mag[k]) * local_mag[j];
            }

            if (j == k) {
              heat_flux_tensor[i][j][k] += q_perp * (1.0 - (local_mag[j] * local_mag[k])) * local_mag[i];
            }
            else {
              heat_flux_tensor[i][j][k] -= q_perp * (local_mag[j] * local_mag[k]) * local_mag[i];
            }
          }
        }
      }

      double heat_flux_tensor_dx[3][3][3];

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            heat_flux_tensor_dx[i][j][k] = local_mag[i] * local_mag[j] * local_mag[k] * q_par_dx;
            heat_flux_tensor_dx[i][j][k] += q_par * local_mag[j] * local_mag[k] * local_mag_dx[i];
            heat_flux_tensor_dx[i][j][k] += q_par * local_mag[i] * local_mag[k] * local_mag_dx[j];
            heat_flux_tensor_dx[i][j][k] += q_par * local_mag[i] * local_mag[j] * local_mag_dx[k];

            if (i == j) {
              heat_flux_tensor_dx[i][j][k] += (1.0 - (local_mag[i] * local_mag[j])) * local_mag[k] * q_perp_dx;
              heat_flux_tensor_dx[i][j][k] += q_perp * (1.0 - (local_mag[i] * local_mag[j])) * local_mag_dx[k];
            }
            else {
              heat_flux_tensor_dx[i][j][k] -= (local_mag[i] * local_mag[j]) * local_mag[k] * q_perp_dx;
              heat_flux_tensor_dx[i][j][k] += q_perp * (local_mag[i] * local_mag[j]) * local_mag_dx[k];
            }

            if (i == k) {
              heat_flux_tensor_dx[i][j][k] += (1.0 - (local_mag[i] * local_mag[k])) * local_mag[j] * q_perp_dx;
              heat_flux_tensor_dx[i][j][k] += q_perp * (1.0 - (local_mag[i] * local_mag[k])) * local_mag_dx[j];
            }
            else {
              heat_flux_tensor_dx[i][j][k] -= (local_mag[i] * local_mag[k]) * local_mag[j] * q_perp_dx;
              heat_flux_tensor_dx[i][j][k] += q_perp * (local_mag[i] * local_mag[k]) * local_mag_dx[j];
            }

            if (j == k) {
              heat_flux_tensor_dx[i][j][k] += (1.0 - (local_mag[j] * local_mag[k])) * local_mag[i] * q_perp_dx;
              heat_flux_tensor_dx[i][j][k] += q_perp * (1.0 - (local_mag[j] * local_mag[k])) * local_mag_dx[i];
            }
            else {
              heat_flux_tensor_dx[i][j][k] -= (local_mag[j] * local_mag[k]) * local_mag[i] * q_perp_dx;
              heat_flux_tensor_dx[i][j][k] += q_perp * (local_mag[j] * local_mag[k]) * local_mag_dx[i];
            }

            heat_flux_tensor_dx[i][j][k] -= ((local_mag[i] * local_mag_dx[j]) + (local_mag[j] * local_mag_dx[i])) * local_mag[k];
            heat_flux_tensor_dx[i][j][k] -= ((local_mag[i] * local_mag_dx[k]) + (local_mag[k] * local_mag_dx[i])) * local_mag[j];
            heat_flux_tensor_dx[i][j][k] -= ((local_mag[j] * local_mag_dx[k]) + (local_mag[k] * local_mag_dx[j])) * local_mag[i];
          }
        }
      }

      divQx[0] = B_avg[0] * heat_flux_tensor_dx[0][0][0];
      divQx[1] = B_avg[0] * heat_flux_tensor_dx[0][0][1];
      divQx[2] = B_avg[0] * heat_flux_tensor_dx[0][0][2];
      divQx[3] = B_avg[0] * heat_flux_tensor_dx[0][1][1];
      divQx[4] = B_avg[0] * heat_flux_tensor_dx[0][1][2];
      divQx[5] = B_avg[0] * heat_flux_tensor_dx[0][2][2];
    }
    else if (poly_order == 2) {
      const double dx = nnclosure->grid.dx[0];
      double rho[3] = { 0.0 };
      double p[3][6] = { 0.0 };
      var_setup(nnclosure, L2_1D, U2_1D, fluid_d, rho, p);

      rho_avg = calc_arithm_avg_1D(rho[L2_1D], rho[U2_1D]);
      p_avg[0] = calc_arithm_avg_1D(p[L2_1D][0], p[U2_1D][0]);
      p_avg[1] = calc_arithm_avg_1D(p[L2_1D][1], p[U2_1D][1]);
      p_avg[2] = calc_arithm_avg_1D(p[L2_1D][2], p[U2_1D][2]);
      p_avg[3] = calc_arithm_avg_1D(p[L2_1D][3], p[U2_1D][3]);
      p_avg[4] = calc_arithm_avg_1D(p[L2_1D][4], p[U2_1D][4]);
      p_avg[5] = calc_arithm_avg_1D(p[L2_1D][5], p[U2_1D][5]);
      B_avg[0] = calc_arithm_avg_1D(em_tot_d[L2_1D][BX], em_tot_d[U2_1D][BX]);
      B_avg[1] = calc_arithm_avg_1D(em_tot_d[L2_1D][BY], em_tot_d[U2_1D][BY]);
      B_avg[2] = calc_arithm_avg_1D(em_tot_d[L2_1D][BZ], em_tot_d[U2_1D][BZ]);

      drho_dx = calc_sym_grad_1D(dx, rho[L2_1D], rho[U2_1D]);
      dp_dx[0] = calc_sym_grad_1D(dx, p[L2_1D][0], p[U2_1D][0]);
      dp_dx[1] = calc_sym_grad_1D(dx, p[L2_1D][1], p[U2_1D][1]);
      dp_dx[2] = calc_sym_grad_1D(dx, p[L2_1D][2], p[U2_1D][2]);
      dp_dx[3] = calc_sym_grad_1D(dx, p[L2_1D][3], p[U2_1D][3]);
      dp_dx[4] = calc_sym_grad_1D(dx, p[L2_1D][4], p[U2_1D][4]);
      dp_dx[5] = calc_sym_grad_1D(dx, p[L2_1D][5], p[U2_1D][5]);
      dB_dx[0] = calc_sym_grad_1D(dx, em_tot_d[L2_1D][BX], em_tot_d[U2_1D][BX]);
      dB_dx[1] = calc_sym_grad_1D(dx, em_tot_d[L2_1D][BY], em_tot_d[U2_1D][BY]);
      dB_dx[2] = calc_sym_grad_1D(dx, em_tot_d[L2_1D][BZ], em_tot_d[U2_1D][BZ]);

      drho_dx_dx = calc_sym_grad2_1D(dx, rho[L2_1D], rho[C2_1D], rho[U2_1D]);
      dp_dx_dx[0] = calc_sym_grad2_1D(dx, p[L2_1D][0], p[C2_1D][0], p[U2_1D][0]);
      dp_dx_dx[1] = calc_sym_grad2_1D(dx, p[L2_1D][1], p[C2_1D][1], p[U2_1D][1]);
      dp_dx_dx[2] = calc_sym_grad2_1D(dx, p[L2_1D][2], p[C2_1D][2], p[U2_1D][2]);
      dp_dx_dx[3] = calc_sym_grad2_1D(dx, p[L2_1D][3], p[C2_1D][3], p[U2_1D][3]);
      dp_dx_dx[4] = calc_sym_grad2_1D(dx, p[L2_1D][4], p[C2_1D][4], p[U2_1D][4]);
      dp_dx_dx[5] = calc_sym_grad2_1D(dx, p[L2_1D][5], p[C2_1D][5], p[U2_1D][5]);
      dB_dx_dx[0] = calc_sym_grad2_1D(dx, em_tot_d[L2_1D][BX], em_tot_d[C2_1D][BX], em_tot_d[U2_1D][BX]);
      dB_dx_dx[1] = calc_sym_grad2_1D(dx, em_tot_d[L2_1D][BY], em_tot_d[C2_1D][BY], em_tot_d[U2_1D][BY]);
      dB_dx_dx[2] = calc_sym_grad2_1D(dx, em_tot_d[L2_1D][BZ], em_tot_d[C2_1D][BZ], em_tot_d[U2_1D][BZ]);

      if (fabs(B_avg[0]) < pow(10.0, -8.0) && fabs(B_avg[1]) < pow(10.0, -8.0) && fabs(B_avg[2]) < pow(10.0, -8.0)) {
        B_avg[0] = 1.0;
        dB_dx[0] = 0.0; dB_dx[1] = 0.0; dB_dx[2] = 0.0;
        dB_dy[0] = 0.0; dB_dy[1] = 0.0; dB_dy[2] = 0.0;
        dB_dz[0] = 0.0; dB_dz[1] = 0.0; dB_dz[2] = 0.0;
        dB_dx_dx[0] = 0.0; dB_dx_dx[1] = 0.0; dB_dx_dx[2] = 0.0;
      }

      double b_mag = sqrt((B_avg[0] * B_avg[0]) + (B_avg[1] * B_avg[1]) + (B_avg[2] * B_avg[2]));
      double b_mag_dx = ((B_avg[0] * dB_dx[0]) + (B_avg[1] * dB_dx[1]) + (B_avg[2] * dB_dx[2])) / b_mag;
      double b_mag_dx_dx = b_mag * ((dB_dx[0] * dB_dx[0]) + (B_avg[0] * dB_dx_dx[0]) + (dB_dx[1] * dB_dx[1]) + (B_avg[1] * dB_dx_dx[1]) + (dB_dx[2] * dB_dy[2]) + (B_avg[2] * dB_dx_dy[2]));
      b_mag_dx_dx -= b_mag_dx * b_mag_dx;
      b_mag_dx_dx /= b_mag * b_mag;

      double local_mag[3];
      for (int i = 0; i < 3; i++) {
        local_mag[i] = B_avg[i] / b_mag;
      }

      double local_mag_dx[3];
      double local_mag_dy[3];
      double local_mag_dx_dx[3];
      for (int i = 0; i < 3; i++) {
        local_mag_dx[i] = ((b_mag * dB_dx[i]) - (B_avg[i] * b_mag_dx)) / (b_mag * b_mag);
        local_mag_dx_dx[i] = (b_mag * b_mag) * ((b_mag_dx * dB_dx[i]) + (b_mag * dB_dx_dx[i]) - (dB_dx[i] * b_mag_dx) * (B_avg[i] * b_mag_dx_dx));
        local_mag_dx_dx[i] -= 2.0 * ((b_mag * dB_dx[i]) - (B_avg[i] * b_mag_dx)) * b_mag * b_mag_dx;
        local_mag_dx_dx[i] /= b_mag * b_mag * b_mag * b_mag;
      }

      double p_tensor[3][3];
      double p_tensor_dx[3][3];
      double p_tensor_dx_dx[3][3];

      p_tensor[0][0] = p_avg[0]; p_tensor[0][1] = p_avg[1]; p_tensor[0][2] = p_avg[2];
      p_tensor[1][0] = p_avg[1]; p_tensor[1][1] = p_avg[3]; p_tensor[1][2] = p_avg[4];
      p_tensor[2][0] = p_avg[2]; p_tensor[2][1] = p_avg[4]; p_tensor[2][2] = p_avg[5];

      p_tensor_dx[0][0] = dp_dx[0]; p_tensor_dx[0][1] = dp_dx[1]; p_tensor_dx[0][2] = dp_dx[2];
      p_tensor_dx[1][0] = dp_dx[1]; p_tensor_dx[1][1] = dp_dx[3]; p_tensor_dx[1][2] = dp_dx[4];
      p_tensor_dx[2][0] = dp_dx[2]; p_tensor_dx[2][1] = dp_dx[4]; p_tensor_dx[2][2] = dp_dx[5];

      p_tensor_dx_dx[0][0] = dp_dx_dx[0]; p_tensor_dx_dx[0][1] = dp_dx_dx[1]; p_tensor_dx_dx[0][2] = dp_dx_dx[2];
      p_tensor_dx_dx[1][0] = dp_dx_dx[1]; p_tensor_dx_dx[1][1] = dp_dx_dx[3]; p_tensor_dx_dx[1][2] = dp_dx_dx[4];
      p_tensor_dx_dx[2][0] = dp_dx_dx[2]; p_tensor_dx_dx[2][1] = dp_dx_dx[4]; p_tensor_dx_dx[2][2] = dp_dx_dx[5];

      double p_par = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          p_par += p_tensor[i][j] * local_mag[i] * local_mag[j];
        }
      }

      double p_tensor_trace = 0.0;
      for (int i = 0; i < 3; i++) {
        p_tensor_trace += p_tensor[i][i];
      }

      double p_perp = 0.5 * (p_tensor_trace - p_par);

      double p_par_dx = 0.0;
      double p_par_dx_dx = 0.0;

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          p_par_dx += local_mag[i] * local_mag[j] * p_tensor_dx[i][j];
          p_par_dx += p_tensor[i][j] * local_mag[j] * local_mag_dx[i];
          p_par_dx += p_tensor[i][j] * local_mag[i] * local_mag_dx[j];

          p_par_dx_dx += local_mag[i] * local_mag[i] * p_tensor_dx_dx[i][j];
          p_par_dx_dx += p_tensor_dx[i][j] * local_mag[i] * local_mag_dx[j];
          p_par_dx_dx += p_tensor_dx[i][j] * local_mag[j] * local_mag_dx[i];
          p_par_dx_dx += p_tensor[i][j] * local_mag[i] * local_mag_dx_dx[j];
          p_par_dx_dx += local_mag_dx[j] * local_mag[i] * p_tensor_dx[i][j];
          p_par_dx_dx += local_mag_dx[j] * p_tensor[i][j] * local_mag_dx[i];
          p_par_dx_dx += p_tensor[i][j] * local_mag[j] * local_mag_dx_dx[i];
          p_par_dx_dx += p_tensor[i][j] * local_mag_dx[j] * local_mag_dx[i];
          p_par_dx_dx += local_mag[j] * local_mag_dx[i] * p_tensor_dx[i][j];
        }
      }

      double p_tensor_trace_dx = 0.0;
      double p_tensor_trace_dx_dx = 0.0;
      for (int i = 0; i < 3; i++) {
        p_tensor_trace_dx += p_tensor_dx[i][i];
        p_tensor_trace_dx_dx += p_tensor_dx_dx[i][i];
      }

      double p_perp_dx = 0.5 * (p_tensor_trace_dx - p_par_dx);
      double p_perp_dx_dx = 0.5 * (p_tensor_trace_dx_dx - p_par_dx_dx);

      input_data[0] = rho_avg;
      input_data[1] = drho_dx;
      input_data[2] = drho_dx_dx;
      input_data[3] = p_par;
      input_data[4] = p_par_dx;
      input_data[5] = p_par_dx_dx;
      input_data[6] = p_perp;
      input_data[7] = p_perp_dx;
      input_data[8] = p_perp_dx_dx;

      const float *output_data_predicted = kann_apply1(ann, input_data);

      for (int i = 0; i < 6; i++) {
        output_data[i] = output_data_predicted[i];
      }

      double q_par = output_data[0];
      double q_par_dx = output_data[1];
      double q_perp = output_data[3];
      double q_perp_dx = output_data[4];

      double heat_flux_tensor[3][3][3];

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            heat_flux_tensor[i][j][k] = q_par * local_mag[i] * local_mag[j] * local_mag[k];

            if (i == j) {
              heat_flux_tensor[i][j][k] += q_perp * (1.0 - (local_mag[i] * local_mag[j])) * local_mag[k];
            }
            else {
              heat_flux_tensor[i][j][k] -= q_perp * (local_mag[i] * local_mag[j]) * local_mag[k];
            }

            if (i == k) {
              heat_flux_tensor[i][j][k] += q_perp * (1.0 - (local_mag[i] * local_mag[k])) * local_mag[j];
            }
            else {
              heat_flux_tensor[i][j][k] -= q_perp * (local_mag[i] * local_mag[k]) * local_mag[j];
            }

            if (j == k) {
              heat_flux_tensor[i][j][k] += q_perp * (1.0 - (local_mag[j] * local_mag[k])) * local_mag[i];
            }
            else {
              heat_flux_tensor[i][j][k] -= q_perp * (local_mag[j] * local_mag[k]) * local_mag[i];
            }
          }
        }
      }

      double heat_flux_tensor_dx[3][3][3];

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            heat_flux_tensor_dx[i][j][k] = local_mag[i] * local_mag[j] * local_mag[k] * q_par_dx;
            heat_flux_tensor_dx[i][j][k] += q_par * local_mag[j] * local_mag[k] * local_mag_dx[i];
            heat_flux_tensor_dx[i][j][k] += q_par * local_mag[i] * local_mag[k] * local_mag_dx[j];
            heat_flux_tensor_dx[i][j][k] += q_par * local_mag[i] * local_mag[j] * local_mag_dx[k];

            if (i == j) {
              heat_flux_tensor_dx[i][j][k] += (1.0 - (local_mag[i] * local_mag[j])) * local_mag[k] * q_perp_dx;
              heat_flux_tensor_dx[i][j][k] += q_perp * (1.0 - (local_mag[i] * local_mag[j])) * local_mag_dx[k];
            }
            else {
              heat_flux_tensor_dx[i][j][k] -= (local_mag[i] * local_mag[j]) * local_mag[k] * q_perp_dx;
              heat_flux_tensor_dx[i][j][k] += q_perp * (local_mag[i] * local_mag[j]) * local_mag_dx[k];
            }

            if (i == k) {
              heat_flux_tensor_dx[i][j][k] += (1.0 - (local_mag[i] * local_mag[k])) * local_mag[j] * q_perp_dx;
              heat_flux_tensor_dx[i][j][k] += q_perp * (1.0 - (local_mag[i] * local_mag[k])) * local_mag_dx[j];
            }
            else {
              heat_flux_tensor_dx[i][j][k] -= (local_mag[i] * local_mag[k]) * local_mag[j] * q_perp_dx;
              heat_flux_tensor_dx[i][j][k] += q_perp * (local_mag[i] * local_mag[k]) * local_mag_dx[j];
            }

            if (j == k) {
              heat_flux_tensor_dx[i][j][k] += (1.0 - (local_mag[j] * local_mag[k])) * local_mag[i] * q_perp_dx;
              heat_flux_tensor_dx[i][j][k] += q_perp * (1.0 - (local_mag[j] * local_mag[k])) * local_mag_dx[i];
            }
            else {
              heat_flux_tensor_dx[i][j][k] -= (local_mag[j] * local_mag[k]) * local_mag[i] * q_perp_dx;
              heat_flux_tensor_dx[i][j][k] += q_perp * (local_mag[j] * local_mag[k]) * local_mag_dx[i];
            }

            heat_flux_tensor_dx[i][j][k] -= ((local_mag[i] * local_mag_dx[j]) + (local_mag[j] * local_mag_dx[i])) * local_mag[k];
            heat_flux_tensor_dx[i][j][k] -= ((local_mag[i] * local_mag_dx[k]) + (local_mag[k] * local_mag_dx[i])) * local_mag[j];
            heat_flux_tensor_dx[i][j][k] -= ((local_mag[j] * local_mag_dx[k]) + (local_mag[k] * local_mag_dx[j])) * local_mag[i];
          }
        }
      }

      divQx[0] = B_avg[0] * heat_flux_tensor_dx[0][0][0];
      divQx[1] = B_avg[0] * heat_flux_tensor_dx[0][0][1];
      divQx[2] = B_avg[0] * heat_flux_tensor_dx[0][0][2];
      divQx[3] = B_avg[0] * heat_flux_tensor_dx[0][1][1];
      divQx[4] = B_avg[0] * heat_flux_tensor_dx[0][1][2];
      divQx[5] = B_avg[0] * heat_flux_tensor_dx[0][2][2];
    }
  }
  else if (ndim == 2) {
    if (poly_order == 1) {
      const double dx = nnclosure->grid.dx[0];
      const double dy = nnclosure->grid.dx[1];
      double rho[4] = { 0.0 };
      double p[4][6] = { 0.0 };
      var_setup(nnclosure, LL_2D, UU_2D, fluid_d, rho, p);

      rho_avg = calc_arithm_avg_2D(rho[LL_2D], rho[LU_2D], rho[UL_2D], rho[UU_2D]);
      p_avg[0] = calc_arithm_avg_2D(p[LL_2D][0], p[LU_2D][0], p[UL_2D][0], p[UU_2D][0]);
      p_avg[1] = calc_arithm_avg_2D(p[LL_2D][1], p[LU_2D][1], p[UL_2D][1], p[UU_2D][1]);
      p_avg[2] = calc_arithm_avg_2D(p[LL_2D][2], p[LU_2D][2], p[UL_2D][2], p[UU_2D][2]);
      p_avg[3] = calc_arithm_avg_2D(p[LL_2D][3], p[LU_2D][3], p[UL_2D][3], p[UU_2D][3]);
      p_avg[4] = calc_arithm_avg_2D(p[LL_2D][4], p[LU_2D][4], p[UL_2D][4], p[UU_2D][4]);
      p_avg[5] = calc_arithm_avg_2D(p[LL_2D][5], p[LU_2D][5], p[UL_2D][5], p[UU_2D][5]);
      B_avg[0] = calc_arithm_avg_2D(em_tot_d[LL_2D][BX], em_tot_d[LU_2D][BX], em_tot_d[UL_2D][BX], em_tot_d[UU_2D][BX]);
      B_avg[1] = calc_arithm_avg_2D(em_tot_d[LL_2D][BY], em_tot_d[LU_2D][BY], em_tot_d[UL_2D][BY], em_tot_d[UU_2D][BY]);
      B_avg[2] = calc_arithm_avg_2D(em_tot_d[LL_2D][BZ], em_tot_d[LU_2D][BZ], em_tot_d[UL_2D][BZ], em_tot_d[UU_2D][BZ]);

      drho_dx = calc_sym_gradx_2D(dx, rho[LL_2D], rho[LU_2D], rho[UL_2D], rho[UU_2D]);
      dp_dx[0] = calc_sym_gradx_2D(dx, p[LL_2D][0], p[LU_2D][0], p[UL_2D][0], p[UU_2D][0]);
      dp_dx[1] = calc_sym_gradx_2D(dx, p[LL_2D][1], p[LU_2D][1], p[UL_2D][1], p[UU_2D][1]);
      dp_dx[2] = calc_sym_gradx_2D(dx, p[LL_2D][2], p[LU_2D][2], p[UL_2D][2], p[UU_2D][2]);
      dp_dx[3] = calc_sym_gradx_2D(dx, p[LL_2D][3], p[LU_2D][3], p[UL_2D][3], p[UU_2D][3]);
      dp_dx[4] = calc_sym_gradx_2D(dx, p[LL_2D][4], p[LU_2D][4], p[UL_2D][4], p[UU_2D][4]);
      dp_dx[5] = calc_sym_gradx_2D(dx, p[LL_2D][5], p[LU_2D][5], p[UL_2D][5], p[UU_2D][5]);
      dB_dx[0] = calc_sym_gradx_2D(dx, em_tot_d[LL_2D][BX], em_tot_d[LU_2D][BX], em_tot_d[UL_2D][BX], em_tot_d[UU_2D][BX]);
      dB_dx[1] = calc_sym_gradx_2D(dx, em_tot_d[LL_2D][BY], em_tot_d[LU_2D][BY], em_tot_d[UL_2D][BY], em_tot_d[UU_2D][BY]);
      dB_dx[2] = calc_sym_gradx_2D(dx, em_tot_d[LL_2D][BZ], em_tot_d[LU_2D][BZ], em_tot_d[UL_2D][BZ], em_tot_d[UU_2D][BZ]);

      drho_dy = calc_sym_grady_2D(dy, rho[LL_2D], rho[LU_2D], rho[UL_2D], rho[UU_2D]);
      dp_dy[0] = calc_sym_grady_2D(dy, p[LL_2D][0], p[LU_2D][0], p[UL_2D][0], p[UU_2D][0]);
      dp_dy[1] = calc_sym_grady_2D(dy, p[LL_2D][1], p[LU_2D][1], p[UL_2D][1], p[UU_2D][1]);
      dp_dy[2] = calc_sym_grady_2D(dy, p[LL_2D][2], p[LU_2D][2], p[UL_2D][2], p[UU_2D][2]);
      dp_dy[3] = calc_sym_grady_2D(dy, p[LL_2D][3], p[LU_2D][3], p[UL_2D][3], p[UU_2D][3]);
      dp_dy[4] = calc_sym_grady_2D(dy, p[LL_2D][4], p[LU_2D][4], p[UL_2D][4], p[UU_2D][4]);
      dp_dy[5] = calc_sym_grady_2D(dy, p[LL_2D][5], p[LU_2D][5], p[UL_2D][5], p[UU_2D][5]);
      dB_dy[0] = calc_sym_grady_2D(dy, em_tot_d[LL_2D][BX], em_tot_d[LU_2D][BX], em_tot_d[UL_2D][BX], em_tot_d[UU_2D][BX]);
      dB_dy[1] = calc_sym_grady_2D(dy, em_tot_d[LL_2D][BY], em_tot_d[LU_2D][BY], em_tot_d[UL_2D][BY], em_tot_d[UU_2D][BY]);
      dB_dy[2] = calc_sym_grady_2D(dy, em_tot_d[LL_2D][BZ], em_tot_d[LU_2D][BZ], em_tot_d[UL_2D][BZ], em_tot_d[UU_2D][BZ]);

      drho_dx_dy = calc_sym_gradxy_2D(dx, dy, rho[LL_2D], rho[LU_2D], rho[UL_2D], rho[UU_2D]);
      dp_dx_dy[0] = calc_sym_gradxy_2D(dx, dy, p[LL_2D][0], p[LU_2D][0], p[UL_2D][0], p[UU_2D][0]);
      dp_dx_dy[1] = calc_sym_gradxy_2D(dx, dy, p[LL_2D][1], p[LU_2D][1], p[UL_2D][1], p[UU_2D][1]);
      dp_dx_dy[2] = calc_sym_gradxy_2D(dx, dy, p[LL_2D][2], p[LU_2D][2], p[UL_2D][2], p[UU_2D][2]);
      dp_dx_dy[3] = calc_sym_gradxy_2D(dx, dy, p[LL_2D][3], p[LU_2D][3], p[UL_2D][3], p[UU_2D][3]);
      dp_dx_dy[4] = calc_sym_gradxy_2D(dx, dy, p[LL_2D][4], p[LU_2D][4], p[UL_2D][4], p[UU_2D][4]);
      dp_dx_dy[5] = calc_sym_gradxy_2D(dx, dy, p[LL_2D][5], p[LU_2D][5], p[UL_2D][5], p[UU_2D][5]);
      dB_dx_dy[0] = calc_sym_gradxy_2D(dx, dy, em_tot_d[LL_2D][BX], em_tot_d[LU_2D][BX], em_tot_d[UL_2D][BX], em_tot_d[UU_2D][BX]);
      dB_dx_dy[1] = calc_sym_gradxy_2D(dx, dy, em_tot_d[LL_2D][BY], em_tot_d[LU_2D][BY], em_tot_d[UL_2D][BY], em_tot_d[UU_2D][BY]);
      dB_dx_dy[2] = calc_sym_gradxy_2D(dx, dy, em_tot_d[LL_2D][BZ], em_tot_d[LU_2D][BZ], em_tot_d[UL_2D][BZ], em_tot_d[UU_2D][BZ]);

      if (fabs(B_avg[0]) < pow(10.0, -8.0) && fabs(B_avg[1]) < pow(10.0, -8.0) && fabs(B_avg[2]) < pow(10.0, -8.0)) {
        B_avg[0] = 1.0;
        dB_dx[0] = 0.0; dB_dx[1] = 0.0; dB_dx[2] = 0.0;
        dB_dy[0] = 0.0; dB_dy[1] = 0.0; dB_dy[2] = 0.0;
        dB_dz[0] = 0.0; dB_dz[1] = 0.0; dB_dz[2] = 0.0;
        dB_dx_dy[0] = 0.0; dB_dx_dy[1] = 0.0; dB_dx_dy[2] = 0.0;
      }

      double b_mag = sqrt((B_avg[0] * B_avg[0]) + (B_avg[1] * B_avg[1]) + (B_avg[2] * B_avg[2]));
      double b_mag_dx = ((B_avg[0] * dB_dx[0]) + (B_avg[1] * dB_dx[1]) + (B_avg[2] * dB_dx[2])) / b_mag;
      double b_mag_dy = ((B_avg[0] * dB_dy[0]) + (B_avg[1] * dB_dy[1]) + (B_avg[2] * dB_dy[2])) / b_mag;
      double b_mag_dx_dy = b_mag * ((dB_dx[0] * dB_dy[0]) + (B_avg[0] * dB_dx_dy[0]) + (dB_dx[1] * dB_dy[1]) + (B_avg[1] * dB_dx_dy[1]) + (dB_dx[2] * dB_dy[2]) + (B_avg[2] * dB_dx_dy[2]));
      b_mag_dx_dy -= b_mag_dx * b_mag_dy;
      b_mag_dx_dy /= b_mag * b_mag;

      double local_mag[3];
      for (int i = 0; i < 3; i++) {
        local_mag[i] = B_avg[i] / b_mag;
      }

      double local_mag_dx[3];
      double local_mag_dy[3];
      double local_mag_dx_dy[3];
      for (int i = 0; i < 3; i++) {
        local_mag_dx[i] = ((b_mag * dB_dx[i]) - (B_avg[i] * b_mag_dx)) / (b_mag * b_mag);
        local_mag_dy[i] = ((b_mag * dB_dy[i]) - (B_avg[i] * b_mag_dy)) / (b_mag * b_mag);
        local_mag_dx_dy[i] = (b_mag * b_mag) * ((b_mag_dy * dB_dx[i]) + (b_mag * dB_dx_dy[i]) - (dB_dy[i] * b_mag_dx) * (B_avg[i] * b_mag_dx_dy));
        local_mag_dx_dy[i] -= 2.0 * ((b_mag * dB_dx[i]) - (B_avg[i] * b_mag_dx)) * b_mag * b_mag_dy;
        local_mag_dx_dy[i] /= b_mag * b_mag * b_mag * b_mag;
      }

      double p_tensor[3][3];
      double p_tensor_dx[3][3];
      double p_tensor_dy[3][3];
      double p_tensor_dx_dy[3][3];

      p_tensor[0][0] = p_avg[0]; p_tensor[0][1] = p_avg[1]; p_tensor[0][2] = p_avg[2];
      p_tensor[1][0] = p_avg[1]; p_tensor[1][1] = p_avg[3]; p_tensor[1][2] = p_avg[4];
      p_tensor[2][0] = p_avg[2]; p_tensor[2][1] = p_avg[4]; p_tensor[2][2] = p_avg[5];

      p_tensor_dx[0][0] = dp_dx[0]; p_tensor_dx[0][1] = dp_dx[1]; p_tensor_dx[0][2] = dp_dx[2];
      p_tensor_dx[1][0] = dp_dx[1]; p_tensor_dx[1][1] = dp_dx[3]; p_tensor_dx[1][2] = dp_dx[4];
      p_tensor_dx[2][0] = dp_dx[2]; p_tensor_dx[2][1] = dp_dx[4]; p_tensor_dx[2][2] = dp_dx[5];

      p_tensor_dy[0][0] = dp_dy[0]; p_tensor_dy[0][1] = dp_dy[1]; p_tensor_dy[0][2] = dp_dy[2];
      p_tensor_dy[1][0] = dp_dy[1]; p_tensor_dy[1][1] = dp_dy[3]; p_tensor_dy[1][2] = dp_dy[4];
      p_tensor_dy[2][0] = dp_dy[2]; p_tensor_dy[2][1] = dp_dy[4]; p_tensor_dy[2][2] = dp_dy[5];

      p_tensor_dx_dy[0][0] = dp_dx_dy[0]; p_tensor_dx_dy[0][1] = dp_dx_dy[1]; p_tensor_dx_dy[0][2] = dp_dx_dy[2];
      p_tensor_dx_dy[1][0] = dp_dx_dy[1]; p_tensor_dx_dy[1][1] = dp_dx_dy[3]; p_tensor_dx_dy[1][2] = dp_dx_dy[4];
      p_tensor_dx_dy[2][0] = dp_dx_dy[2]; p_tensor_dx_dy[2][1] = dp_dx_dy[4]; p_tensor_dx_dy[2][2] = dp_dx_dy[5];

      double p_par = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          p_par += p_tensor[i][j] * local_mag[i] * local_mag[j];
        }
      }

      double p_tensor_trace = 0.0;
      for (int i = 0; i < 3; i++) {
        p_tensor_trace += p_tensor[i][i];
      }

      double p_perp = 0.5 * (p_tensor_trace - p_par);

      double p_par_dx = 0.0;
      double p_par_dy = 0.0;
      double p_par_dx_dy = 0.0;

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          p_par_dx += local_mag[i] * local_mag[j] * p_tensor_dx[i][j];
          p_par_dx += p_tensor[i][j] * local_mag[j] * local_mag_dx[i];
          p_par_dx += p_tensor[i][j] * local_mag[i] * local_mag_dx[j];

          p_par_dy += local_mag[i] * local_mag[j] * p_tensor_dy[i][j];
          p_par_dy += p_tensor[i][j] * local_mag[j] * local_mag_dy[i];
          p_par_dy += p_tensor[i][j] * local_mag[i] * local_mag_dy[j];

          p_par_dx_dy += local_mag[i] * local_mag[i] * p_tensor_dx_dy[i][j];
          p_par_dx_dy += p_tensor_dx[i][j] * local_mag[i] * local_mag_dy[j];
          p_par_dx_dy += p_tensor_dx[i][j] * local_mag[j] * local_mag_dy[i];
          p_par_dx_dy += p_tensor[i][j] * local_mag[i] * local_mag_dx_dy[j];
          p_par_dx_dy += local_mag_dx[j] * local_mag[i] * p_tensor_dy[i][j];
          p_par_dx_dy += local_mag_dx[j] * p_tensor[i][j] * local_mag_dy[i];
          p_par_dx_dy += p_tensor[i][j] * local_mag[j] * local_mag_dx_dy[i];
          p_par_dx_dy += p_tensor[i][j] * local_mag_dy[j] * local_mag_dx[i];
          p_par_dx_dy += local_mag[j] * local_mag_dx[i] * p_tensor_dy[i][j];
        }
      }

      double p_tensor_trace_dx = 0.0;
      double p_tensor_trace_dy = 0.0;
      double p_tensor_trace_dx_dy = 0.0;
      for (int i = 0; i < 3; i++) {
        p_tensor_trace_dx += p_tensor_dx[i][i];
        p_tensor_trace_dy += p_tensor_dy[i][i];
        p_tensor_trace_dx_dy += p_tensor_dx_dy[i][i];
      }

      double p_perp_dx = 0.5 * (p_tensor_trace_dx - p_par_dx);
      double p_perp_dy = 0.5 * (p_tensor_trace_dy - p_par_dy);
      double p_perp_dx_dy = 0.5 * (p_tensor_trace_dx_dy - p_par_dx_dy);

      input_data[0] = rho_avg;
      input_data[1] = drho_dx;
      input_data[2] = drho_dy;
      input_data[3] = drho_dx_dy;
      input_data[4] = p_par;
      input_data[5] = p_par_dx;
      input_data[6] = p_par_dy;
      input_data[7] = p_par_dx_dy;
      input_data[8] = p_perp;
      input_data[9] = p_perp_dx;
      input_data[10] = p_perp_dy;
      input_data[11] = p_perp_dx_dy;

      const float *output_data_predicted = kann_apply1(ann, input_data);

      for (int i = 0; i < 8; i++) {
        output_data[i] = output_data_predicted[i];
      }

      double q_par = output_data[0];
      double q_par_dx = output_data[1];
      double q_par_dy = output_data[2];
      double q_perp = output_data[4];
      double q_perp_dx = output_data[5];
      double q_perp_dy = output_data[6];

      double heat_flux_tensor[3][3][3];

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            heat_flux_tensor[i][j][k] = q_par * local_mag[i] * local_mag[j] * local_mag[k];

            if (i == j) {
              heat_flux_tensor[i][j][k] += q_perp * (1.0 - (local_mag[i] * local_mag[j])) * local_mag[k];
            }
            else {
              heat_flux_tensor[i][j][k] -= q_perp * (local_mag[i] * local_mag[j]) * local_mag[k];
            }

            if (i == k) {
              heat_flux_tensor[i][j][k] += q_perp * (1.0 - (local_mag[i] * local_mag[k])) * local_mag[j];
            }
            else {
              heat_flux_tensor[i][j][k] -= q_perp * (local_mag[i] * local_mag[k]) * local_mag[j];
            }

            if (j == k) {
              heat_flux_tensor[i][j][k] += q_perp * (1.0 - (local_mag[j] * local_mag[k])) * local_mag[i];
            }
            else {
              heat_flux_tensor[i][j][k] -= q_perp * (local_mag[j] * local_mag[k]) * local_mag[i];
            }
          }
        }
      }

      double heat_flux_tensor_dx[3][3][3];
      double heat_flux_tensor_dy[3][3][3];

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            heat_flux_tensor_dx[i][j][k] = local_mag[i] * local_mag[j] * local_mag[k] * q_par_dx;
            heat_flux_tensor_dx[i][j][k] += q_par * local_mag[j] * local_mag[k] * local_mag_dx[i];
            heat_flux_tensor_dx[i][j][k] += q_par * local_mag[i] * local_mag[k] * local_mag_dx[j];
            heat_flux_tensor_dx[i][j][k] += q_par * local_mag[i] * local_mag[j] * local_mag_dx[k];

            heat_flux_tensor_dy[i][j][k] = local_mag[i] * local_mag[j] * local_mag[k] * q_par_dy;
            heat_flux_tensor_dy[i][j][k] += q_par * local_mag[j] * local_mag[k] * local_mag_dy[i];
            heat_flux_tensor_dy[i][j][k] += q_par * local_mag[i] * local_mag[k] * local_mag_dy[j];
            heat_flux_tensor_dy[i][j][k] += q_par * local_mag[i] * local_mag[j] * local_mag_dy[k];

            if (i == j) {
              heat_flux_tensor_dx[i][j][k] += (1.0 - (local_mag[i] * local_mag[j])) * local_mag[k] * q_perp_dx;
              heat_flux_tensor_dx[i][j][k] += q_perp * (1.0 - (local_mag[i] * local_mag[j])) * local_mag_dx[k];

              heat_flux_tensor_dy[i][j][k] += (1.0 - (local_mag[i] * local_mag[j])) * local_mag[k] * q_perp_dy;
              heat_flux_tensor_dy[i][j][k] += q_perp * (1.0 - (local_mag[i] * local_mag[j])) * local_mag_dy[k];
            }
            else {
              heat_flux_tensor_dx[i][j][k] -= (local_mag[i] * local_mag[j]) * local_mag[k] * q_perp_dx;
              heat_flux_tensor_dx[i][j][k] += q_perp * (local_mag[i] * local_mag[j]) * local_mag_dx[k];

              heat_flux_tensor_dy[i][j][k] -= (local_mag[i] * local_mag[j]) * local_mag[k] * q_perp_dy;
              heat_flux_tensor_dy[i][j][k] += q_perp * (local_mag[i] * local_mag[j]) * local_mag_dy[k];
            }

            if (i == k) {
              heat_flux_tensor_dx[i][j][k] += (1.0 - (local_mag[i] * local_mag[k])) * local_mag[j] * q_perp_dx;
              heat_flux_tensor_dx[i][j][k] += q_perp * (1.0 - (local_mag[i] * local_mag[k])) * local_mag_dx[j];

              heat_flux_tensor_dy[i][j][k] += (1.0 - (local_mag[i] * local_mag[k])) * local_mag[j] * q_perp_dy;
              heat_flux_tensor_dy[i][j][k] += q_perp * (1.0 - (local_mag[i] * local_mag[k])) * local_mag_dy[j];
            }
            else {
              heat_flux_tensor_dx[i][j][k] -= (local_mag[i] * local_mag[k]) * local_mag[j] * q_perp_dx;
              heat_flux_tensor_dx[i][j][k] += q_perp * (local_mag[i] * local_mag[k]) * local_mag_dx[j];

              heat_flux_tensor_dy[i][j][k] -= (local_mag[i] * local_mag[k]) * local_mag[j] * q_perp_dy;
              heat_flux_tensor_dy[i][j][k] += q_perp * (local_mag[i] * local_mag[k]) * local_mag_dy[j];
            }

            if (j == k) {
              heat_flux_tensor_dx[i][j][k] += (1.0 - (local_mag[j] * local_mag[k])) * local_mag[i] * q_perp_dx;
              heat_flux_tensor_dx[i][j][k] += q_perp * (1.0 - (local_mag[j] * local_mag[k])) * local_mag_dx[i];

              heat_flux_tensor_dy[i][j][k] += (1.0 - (local_mag[j] * local_mag[k])) * local_mag[i] * q_perp_dy;
              heat_flux_tensor_dy[i][j][k] += q_perp * (1.0 - (local_mag[j] * local_mag[k])) * local_mag_dy[i];
            }
            else {
              heat_flux_tensor_dx[i][j][k] -= (local_mag[j] * local_mag[k]) * local_mag[i] * q_perp_dx;
              heat_flux_tensor_dx[i][j][k] += q_perp * (local_mag[j] * local_mag[k]) * local_mag_dx[i];

              heat_flux_tensor_dy[i][j][k] -= (local_mag[j] * local_mag[k]) * local_mag[i] * q_perp_dy;
              heat_flux_tensor_dy[i][j][k] += q_perp * (local_mag[j] * local_mag[k]) * local_mag_dy[i];
            }

            heat_flux_tensor_dx[i][j][k] -= ((local_mag[i] * local_mag_dx[j]) + (local_mag[j] * local_mag_dx[i])) * local_mag[k];
            heat_flux_tensor_dx[i][j][k] -= ((local_mag[i] * local_mag_dx[k]) + (local_mag[k] * local_mag_dx[i])) * local_mag[j];
            heat_flux_tensor_dx[i][j][k] -= ((local_mag[j] * local_mag_dx[k]) + (local_mag[k] * local_mag_dx[j])) * local_mag[i];

            heat_flux_tensor_dy[i][j][k] -= ((local_mag[i] * local_mag_dy[j]) + (local_mag[j] * local_mag_dy[i])) * local_mag[k];
            heat_flux_tensor_dy[i][j][k] -= ((local_mag[i] * local_mag_dy[k]) + (local_mag[k] * local_mag_dy[i])) * local_mag[j];
            heat_flux_tensor_dy[i][j][k] -= ((local_mag[j] * local_mag_dy[k]) + (local_mag[k] * local_mag_dy[j])) * local_mag[i];
          }
        }
      }

      divQx[0] = B_avg[0] * heat_flux_tensor_dx[0][0][0];
      divQx[1] = B_avg[0] * heat_flux_tensor_dx[0][0][1];
      divQx[2] = B_avg[0] * heat_flux_tensor_dx[0][0][2];
      divQx[3] = B_avg[0] * heat_flux_tensor_dx[0][1][1];
      divQx[4] = B_avg[0] * heat_flux_tensor_dx[0][1][2];
      divQx[5] = B_avg[0] * heat_flux_tensor_dx[0][2][2];

      divQy[0] = B_avg[1] * heat_flux_tensor_dy[0][0][1];
      divQy[1] = B_avg[1] * heat_flux_tensor_dy[0][1][1];
      divQy[2] = B_avg[1] * heat_flux_tensor_dy[0][1][2];
      divQy[3] = B_avg[1] * heat_flux_tensor_dy[1][1][1];
      divQy[4] = B_avg[1] * heat_flux_tensor_dy[1][1][2];
      divQy[5] = B_avg[1] * heat_flux_tensor_dy[1][2][2];
    }
  }

  double k0 = nnclosure->k0;
  double alpha = 1.0 / k0;

  double p_scalar = (p_avg[0] + p_avg[3] + p_avg[5]) / 3.0;
  double vth_avg = sqrt(p_scalar / rho_avg);

  rhs[RHO] = 0.0;
  rhs[MX] = 0.0;
  rhs[MY] = 0.0;
  rhs[MZ] = 0.0;
  rhs[P11] = alpha * rho_avg * vth_avg * (divQx[0] + divQy[0] + divQz[0]);
  rhs[P12] = alpha * rho_avg * vth_avg * (divQx[1] + divQy[1] + divQz[1]);
  rhs[P13] = alpha * rho_avg * vth_avg * (divQx[2] + divQy[2] + divQz[2]);
  rhs[P22] = alpha * rho_avg * vth_avg * (divQx[3] + divQy[3] + divQz[3]);
  rhs[P23] = alpha * rho_avg * vth_avg * (divQx[4] + divQy[4] + divQz[4]);
  rhs[P33] = alpha * rho_avg * vth_avg * (divQx[5] + divQy[5] + divQz[5]);

  gkyl_free(input_data);
  gkyl_free(output_data);
}

void
gkyl_ten_moment_nn_closure_advance(const gkyl_ten_moment_nn_closure *nnclosure, const struct gkyl_range *heat_flux_rng, const struct gkyl_range *update_rng,
  const struct gkyl_array *fluid, const struct gkyl_array *em_tot, struct gkyl_array *heat_flux, struct gkyl_array *rhs)
{
  int poly_order = nnclosure->poly_order;
  int ndim = update_rng->ndim;
  long sz[] = { 2, 4, 8 };
  long sz_p2[] = { 3, 9, 27 };

  long offsets_vertices[sz[ndim - 1]];
  long offsets_vertices_p2[sz_p2[ndim - 1]];
  if (poly_order == 1) {
    create_offsets_vertices(update_rng, offsets_vertices);
  }
  else {
    create_offsets_vertices(update_rng, offsets_vertices_p2);
  }

  long offsets_centers[sz[ndim - 1]];
  long offsets_centers_p2[sz_p2[ndim - 1]];
  if (poly_order == 1) {
    create_offsets_centers(heat_flux_rng, offsets_centers);
  }
  else {
    create_offsets_centers(heat_flux_rng, offsets_centers_p2);
  }

  const double *fluid_d[sz[ndim - 1]];
  const double *fluid_d_p2[sz_p2[ndim - 1]];
  const double *em_tot_d[sz[ndim - 1]];
  const double *em_tot_d_p2[sz_p2[ndim - 1]];
  double *heat_flux_d;
  const double *heat_flux_up[sz[ndim - 1]];
  const double *heat_flux_up_p2[sz_p2[ndim - 1]];
  double *rhs_d;

  struct gkyl_range_iter iter_vertex;
  gkyl_range_iter_init(&iter_vertex, heat_flux_rng);
  while (gkyl_range_iter_next(&iter_vertex)) {
    long linc_vertex = gkyl_range_idx(heat_flux_rng, iter_vertex.idx);
    long linc_center = gkyl_range_idx(update_rng, iter_vertex.idx);

    if (poly_order == 1) {
      for (int i = 0; i < sz[ndim - 1]; i++) {
        em_tot_d[i] = gkyl_array_cfetch(em_tot, linc_center + offsets_vertices[i]);
        fluid_d[i] = gkyl_array_cfetch(fluid, linc_center + offsets_vertices[i]);
      }
    }
    else {
      for (int i = 0; i < sz_p2[ndim - 1]; i++) {
        em_tot_d_p2[i] = gkyl_array_cfetch(em_tot, linc_center + offsets_vertices_p2[i]);
        fluid_d_p2[i] = gkyl_array_cfetch(fluid, linc_center + offsets_vertices_p2[i]);
      }
    }

    heat_flux_d = gkyl_array_fetch(heat_flux, linc_vertex);

    // Unnecessary to calculate heat flux separately for now, so commenting this out for efficiency.
    // calc_mag_heat_flux(nnclosure, fluid_d, em_tot_d, heat_flux_d);
  }

  struct gkyl_range_iter iter_center;
  gkyl_range_iter_init(&iter_center, update_rng);
  while (gkyl_range_iter_next(&iter_center)) {
    long linc_vertex = gkyl_range_idx(heat_flux_rng, iter_center.idx);
    long linc_center = gkyl_range_idx(update_rng, iter_center.idx);

    if (poly_order == 1) {
      for (int i = 0; i < sz[ndim - 1]; i++) {
        em_tot_d[i] = gkyl_array_cfetch(em_tot, linc_vertex + offsets_centers[i]);
        fluid_d[i] = gkyl_array_cfetch(fluid, linc_vertex + offsets_centers[i]);
      }
    }
    else {
      for (int i = 0; i < sz_p2[ndim - 1]; i++) {
        em_tot_d_p2[i] = gkyl_array_cfetch(em_tot, linc_vertex + offsets_centers_p2[i]);
        fluid_d_p2[i] = gkyl_array_cfetch(fluid, linc_vertex + offsets_centers_p2[i]);
      }
    }

    rhs_d = gkyl_array_fetch(rhs, linc_center);

    if (poly_order == 1) {
      calc_nn_closure_update(nnclosure, fluid_d, em_tot_d, rhs_d);
    }
    else {
      calc_nn_closure_update(nnclosure, fluid_d_p2, em_tot_d_p2, rhs_d);
    }
  }
}

gkyl_ten_moment_nn_closure*
gkyl_ten_moment_nn_closure_new(struct gkyl_ten_moment_nn_closure_inp inp)
{
  gkyl_ten_moment_nn_closure *up = gkyl_malloc(sizeof(gkyl_ten_moment_nn_closure));

  up->grid = *(inp.grid);
  up->ndim = up->grid.ndim;
  up->k0 = inp.k0;
  up->poly_order = inp.poly_order;
  up->ann = inp.ann;

  return up;
}

void
gkyl_ten_moment_nn_closure_release(gkyl_ten_moment_nn_closure *nnclosure)
{
  free(nnclosure);
}