#include <math.h>
#include <stdlib.h>

#include <gkyl_alloc.h>
#include <gkyl_ref_count.h>
#include <gkyl_gr_blackhole.h>
#include <gkyl_gr_spacetime.h>
#include <gkyl_gr_spacetime_diff.h>

#include "gkyl_gr_tov_spacetime.h"
#include "tov_solver.h"


struct gr_tov_spacetime {
  struct gkyl_gr_spacetime spacetime;

  struct gkyl_tov *tov; 
  struct gkyl_gr_spacetime *bh;

  double pos_c[3];  // star centre
  double R_star;
  double M_star;
};

static inline struct gr_tov_spacetime*
tov_st(const struct gkyl_gr_spacetime *spacetime)
{
  return (struct gr_tov_spacetime*)spacetime;
}


static inline double
star_radius(const struct gr_tov_spacetime *ts, double x, double y, double z)
{
  double dx = x - ts->pos_c[0];
  double dy = y - ts->pos_c[1];
  double dz = z - ts->pos_c[2];
  return sqrt(dx * dx + dy * dy + dz * dz);
}


//  BL -> Spherical KS
//  Time transformation:  dt_KS = dt_BL + h'(r) dr
//  with  h'(r) = 2m(r) / [ r (1 - 2m(r)/r) ]
//  Spherical KS 3+1:
//  alpha_KS = sqrt( alpha_BL^2 + beta_r^2 / spatial_metric_rr )
//  beta_r = alpha_BL^2  h'
//  spatial_metric_rr = A_BL - alpha_BL^2  h'^2


struct sks_data {
  double alpha;
  double beta_r;
  double spatial_metric_rr;
};

static bool
eval_cks_interior(const struct gkyl_tov *tov,
    double dx, double dy, double dz,
    double *lapse, double shift_vector[3], double spatial_metric[3][3])
{
  double r = sqrt(dx*dx + dy*dy + dz*dz);
  double r_safe = fmax(r, 1e-6);

  // Step 1: BL metric from TOV arrays
  struct tov_eval_bl bl;
  if (!gkyl_tov_eval_bl(tov, r_safe, &bl)) return false;

  double alpha_BL = exp(bl.Phi);
  double f        = 1.0 - 2.0*bl.m / r_safe;
  double A_BL     = (f > 1e-10) ? 1.0/f : 1e10;

  // Step 2: BL -> SKS
  double h_prime       = 2.0*bl.m * A_BL / r_safe;
  double beta_r        = alpha_BL*alpha_BL * h_prime;
  double gamma_rr      = A_BL - alpha_BL*alpha_BL * h_prime*h_prime;
  double lapse_SKS     = (gamma_rr > 1e-10) ?
    sqrt(alpha_BL*alpha_BL + beta_r*beta_r / gamma_rr) : 0.0;

  // Step 3: SKS -> CKS
  double n[3] = {dx/r_safe, dy/r_safe, dz/r_safe};

  *lapse = lapse_SKS;

  for (int i = 0; i < 3; i++) {
    shift_vector[i] = beta_r * n[i];
    for (int j = 0; j < 3; j++)
      spatial_metric[i][j] = (i==j ? 1.0 : 0.0) + (gamma_rr - 1.0)*n[i]*n[j];
  }

  return true;
}

static void
tov_lapse_function(const struct gkyl_gr_spacetime *spacetime, const double t,
    const double x, const double y, const double z, double *lapse)
{
    struct gr_tov_spacetime *ts = tov_st(spacetime);
    double r = star_radius(ts, x, y, z);

    if (r < ts->R_star) {
        double shift_vector[3], spatial_metric[3][3];
        eval_cks_interior(ts->tov,
            x - ts->pos_c[0], y - ts->pos_c[1], z - ts->pos_c[2],
            lapse, shift_vector, spatial_metric);
    }
    else {
        ts->bh->lapse_function_func(ts->bh, t, x, y, z, lapse);
    }
}

static void
tov_shift_vector(const struct gkyl_gr_spacetime *spacetime, const double t,
    const double x, const double y, const double z, double **shift)
{
  struct gr_tov_spacetime *ts = tov_st(spacetime);
  double r = star_radius(ts, x, y, z);

  if (r < ts->R_star) {
    double lapse;
    double shift_vector[3];
    double spatial_metric[3][3];
    eval_cks_interior(ts->tov,
        x - ts->pos_c[0], y - ts->pos_c[1], z - ts->pos_c[2],
        &lapse, shift_vector, spatial_metric);
    (*shift)[0] = shift_vector[0];
    (*shift)[1] = shift_vector[1];
    (*shift)[2] = shift_vector[2];
  }
  else {
      ts->bh->shift_vector_func(ts->bh, t, x, y, z, shift);
  }
}

static void
tov_spatial_metric_tensor(const struct gkyl_gr_spacetime *spacetime, const double t,
    const double x, const double y, const double z, double ***spatial_metric)
{
  struct gr_tov_spacetime *ts = tov_st(spacetime);
  double r = star_radius(ts, x, y, z);

  if (r < ts->R_star) {
    double lapse;
    double shift_vector[3];
    double spatial_metric_local[3][3];
    eval_cks_interior(ts->tov,
        x - ts->pos_c[0], y - ts->pos_c[1], z - ts->pos_c[2],
        &lapse, shift_vector, spatial_metric_local);
    for (int i = 0; i < 3; i++) 
        for (int j = 0; j < 3; j++)
            (*spatial_metric)[i][j] = spatial_metric_local[i][j];
  }
  else {
      ts->bh->spatial_metric_tensor_func(ts->bh, t, x, y, z, spatial_metric);
  }
}

static void
tov_spatial_metric_det(const struct gkyl_gr_spacetime* spacetime, const double t,
const double x, const double y, const double z, double* spatial_metric_det)
{
  struct gr_tov_spacetime *ts = tov_st(spacetime);
  double r = star_radius(ts, x, y, z);
  if (r < ts->R_star) {
    double lapse, shift_vector[3], spatial_metric_local[3][3];
    eval_cks_interior(ts->tov, x - ts->pos_c[0], y - ts->pos_c[1], z - ts->pos_c[2], &lapse, shift_vector, spatial_metric_local);
    double dx = x - ts->pos_c[0], dy = y - ts->pos_c[1], dz = z - ts->pos_c[2];
    double r_safe = fmax(r, 1e-6);
    double n[3] = {dx/r_safe, dy/r_safe, dz/r_safe};
    double spatial_metric_local_rr = 0.0;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        spatial_metric_local_rr += spatial_metric_local[i][j] * n[i] * n[j];
      }
    }
    *spatial_metric_det = spatial_metric_local_rr;
  }
  else {
    ts->bh->spatial_metric_det_func(ts->bh, t, x, y, z, spatial_metric_det);
  }
}

static void
tov_spatial_inv_metric_tensor(const struct gkyl_gr_spacetime *spacetime, const double t,
    const double x, const double y, const double z, double ***spatial_inv_metric_tensor)
{
    struct gr_tov_spacetime *ts = tov_st(spacetime);
    double r = star_radius(ts, x, y, z);

    if (r < ts->R_star) {
        double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++)
            spatial_metric[i] = gkyl_malloc(sizeof(double[3]));

        tov_spatial_metric_tensor(spacetime, t, x, y, z, &spatial_metric);

        double spatial_metric_det =
            (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
            (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
            (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));

        double trace = 0.0;
        for (int i = 0; i < 3; i++)
            trace += spatial_metric[i][i];

        double **spatial_metric_sq = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
            spatial_metric_sq[i] = gkyl_malloc(sizeof(double[3]));
            for (int j = 0; j < 3; j++)
                spatial_metric_sq[i][j] = 0.0;
        }

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                    spatial_metric_sq[i][j] += spatial_metric[i][k] * spatial_metric[k][j];

        double sq_trace = 0.0;
        for (int i = 0; i < 3; i++)
            sq_trace += spatial_metric_sq[i][i];

        double **euclidean_metric = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) {
            euclidean_metric[i] = gkyl_malloc(sizeof(double[3]));
            for (int j = 0; j < 3; j++)
                euclidean_metric[i][j] = (i == j) ? 1.0 : 0.0;
        }

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                (*spatial_inv_metric_tensor)[i][j] = (1.0 / spatial_metric_det) *
                    ((0.5 * ((trace * trace) - sq_trace) * euclidean_metric[i][j]) -
                     (trace * spatial_metric[i][j]) + spatial_metric_sq[i][j]);

        for (int i = 0; i < 3; i++) {
            gkyl_free(spatial_metric[i]);
            gkyl_free(spatial_metric_sq[i]);
            gkyl_free(euclidean_metric[i]);
        }
        gkyl_free(spatial_metric);
        gkyl_free(spatial_metric_sq);
        gkyl_free(euclidean_metric);
    }
    else {
        ts->bh->spatial_inv_metric_tensor_func(ts->bh, t, x, y, z, spatial_inv_metric_tensor);
    }
}

static void
tov_lapse_function_der(const struct gkyl_gr_spacetime* spacetime, const double t,
  const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double** lapse_function_der)
{
  
  gkyl_gr_lapse_function_diff(spacetime, t, x, y, z, dx, dy, dz, lapse_function_der);
}

static void
tov_shift_vector_der(const struct gkyl_gr_spacetime* spacetime, const double t,
  const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** shift_vector_der)
{
  gkyl_gr_shift_vector_diff(spacetime, t, x, y, z, dx, dy, dz, shift_vector_der);
}

static void
tov_spatial_metric_tensor_der(const struct gkyl_gr_spacetime* spacetime, const double t,
  const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double**** spatial_metric_tensor_der)
{
  gkyl_gr_spatial_metric_tensor_diff(spacetime, t, x, y, z, dx, dy, dz, spatial_metric_tensor_der);
}

static void
tov_spatial_christoffel(const struct gkyl_gr_spacetime* spacetime, const double t,
  const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double**** spatial_christoffel)
{
  gkyl_gr_spatial_christoffel_fd(spacetime, t, x, y, z, dx, dy, dz, spatial_christoffel);
}

static void
tov_extrinsic_curvature_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** extrinsic_curvature_tensor)
{
  double lapse_function;
  double *shift_vector = gkyl_malloc(sizeof(double[3]));

  double **spatial_metric = gkyl_malloc(sizeof(double*[3]));
  double **shift_vector_der = gkyl_malloc(sizeof(double*[3]));
  double **shift_vector_cov_der = gkyl_malloc(sizeof(double*[3]));
  double **shift_covector_cov_der = gkyl_malloc(sizeof(double*[3]));

  double ***spatial_christoffel = gkyl_malloc(sizeof(double**[3]));

  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
    shift_vector_der[i] = gkyl_malloc(sizeof(double[3]));
    shift_vector_cov_der[i] = gkyl_malloc(sizeof(double[3]));
    shift_covector_cov_der[i] = gkyl_malloc(sizeof(double[3]));

    spatial_christoffel[i] = gkyl_malloc(sizeof(double*[3]));
    for (int j = 0; j < 3; j++) {
      spatial_christoffel[i][j] = gkyl_malloc(sizeof(double[3]));

      shift_covector_cov_der[i][j] = 0.0;
    }
  }

  tov_lapse_function(spacetime, t, x, y, z, &lapse_function);
  tov_shift_vector(spacetime, t, x, y, z, &shift_vector);
  tov_spatial_metric_tensor(spacetime, t, x, y, z, &spatial_metric);
  tov_shift_vector_der(spacetime, t, x, y, z, dx, dy, dz, &shift_vector_der);
  tov_spatial_christoffel(spacetime, t, x, y, z, dx, dy, dz, &spatial_christoffel);

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
      (*extrinsic_curvature_tensor)[i][j] = -(1.0 / (2.0 * lapse_function)) * (shift_covector_cov_der[j][i] + shift_covector_cov_der[i][j]);
    }
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric[i]);
    gkyl_free(shift_vector_der[i]);
    gkyl_free(shift_vector_cov_der[i]);
    gkyl_free(shift_covector_cov_der[i]);

    for (int j = 0; j < 3; j++) {
      gkyl_free(spatial_christoffel[i][j]);
    }
    gkyl_free(spatial_christoffel[i]);
  }
  gkyl_free(shift_vector);
  gkyl_free(spatial_metric);
  gkyl_free(shift_vector_der);
  gkyl_free(shift_vector_cov_der);
  gkyl_free(shift_covector_cov_der);
  gkyl_free(spatial_christoffel);
}

static void
tov_excision_region(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  bool* in_excision_region)
{
  *in_excision_region = false;
}

static void
tov_verify_metric(const struct gkyl_gr_spacetime *spacetime)
{
    struct gr_tov_spacetime *ts = tov_st(spacetime);
    double M = gkyl_tov_star_mass(ts->tov);
    double R = gkyl_tov_star_radius(ts->tov);

    printf("\n=== METRIC VERIFICATION (on x-axis) ===\n");
    printf("%-8s %-10s %-10s %-10s %-10s %-10s %-10s\n",
        "r", "lapse", "shift_x", "gam_xx", "gam_yy", "gam_xy", "det");

    double r_tests[5] = {1e-4, 0.25*R, 0.5*R, R-0.05, R+0.05};
    for (int p = 0; p < 5; p++) {
        double r = r_tests[p];
        double xp = ts->pos_c[0] + r;

        double lapse;
        double *shift = gkyl_malloc(sizeof(double[3]));
        double **gam = gkyl_malloc(sizeof(double*[3]));
        for (int i = 0; i < 3; i++) gam[i] = gkyl_malloc(sizeof(double[3]));
        double det;

        tov_lapse_function(spacetime, 0.0, xp, ts->pos_c[1], ts->pos_c[2], &lapse);
        tov_shift_vector(spacetime, 0.0, xp, ts->pos_c[1], ts->pos_c[2], &shift);
        tov_spatial_metric_tensor(spacetime, 0.0, xp, ts->pos_c[1], ts->pos_c[2], &gam);
        tov_spatial_metric_det(spacetime, 0.0, xp, ts->pos_c[1], ts->pos_c[2], &det);

        // on x-axis: n=(1,0,0), so gamma_xx=gamma_rr, gamma_yy=gamma_zz=1, gamma_xy=0
        // det = gamma_rr = gamma_xx
        // shift_x = beta_r, shift_y=shift_z=0
        printf("%-8.4f %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f\n",
            r, lapse, shift[0], gam[0][0], gam[1][1], gam[0][1], det);

        // check det = gamma_xx on x-axis
        if (fabs(det - gam[0][0]) > 1e-10)
            printf("  !! det != gamma_xx: det=%.6f gamma_xx=%.6f\n", det, gam[0][0]);
        // check gamma_yy = 1 on x-axis
        if (fabs(gam[1][1] - 1.0) > 1e-10)
            printf("  !! gamma_yy != 1: %.6f\n", gam[1][1]);
        // check gamma_xy = 0
        if (fabs(gam[0][1]) > 1e-10)
            printf("  !! gamma_xy != 0: %.6f\n", gam[0][1]);
        // check shift_y = shift_z = 0
        if (fabs(shift[1]) > 1e-10 || fabs(shift[2]) > 1e-10)
            printf("  !! shift_y=%.2e shift_z=%.2e (should be 0)\n", shift[1], shift[2]);

        gkyl_free(shift);
        for (int i = 0; i < 3; i++) gkyl_free(gam[i]);
        gkyl_free(gam);
    }

    // surface continuity check
    printf("\n--- Surface continuity ---\n");
    double r_in  = R - 0.05;
    double r_out = R + 0.05;
    double lapse_in, lapse_out;
    tov_lapse_function(spacetime, 0.0, ts->pos_c[0]+r_in,  ts->pos_c[1], ts->pos_c[2], &lapse_in);
    tov_lapse_function(spacetime, 0.0, ts->pos_c[0]+r_out, ts->pos_c[1], ts->pos_c[2], &lapse_out);
    printf("lapse just inside =%.6f, just outside=%.6f, jump=%.2e\n",
        lapse_in, lapse_out, lapse_out - lapse_in);
}

static void
tov_verify_lapse_derivative(const struct gkyl_gr_spacetime *spacetime)
{
    struct gr_tov_spacetime *ts = tov_st(spacetime);
    double R = gkyl_tov_star_radius(ts->tov);

    printf("\n=== LAPSE DERIVATIVE VERIFICATION ===\n");

    // check at r=R/2, on x-axis
    double xp = ts->pos_c[0] + 0.5*R;
    double h = 1e-4;

    // finite difference by hand
    double lapse_p, lapse_m;
    tov_lapse_function(spacetime, 0.0, xp+h, ts->pos_c[1], ts->pos_c[2], &lapse_p);
    tov_lapse_function(spacetime, 0.0, xp-h, ts->pos_c[1], ts->pos_c[2], &lapse_m);
    double dalpha_dx_fd = (lapse_p - lapse_m) / (2.0*h);

    // what the derivative function returns
    double **lapse_der = gkyl_malloc(sizeof(double*[1]));
    lapse_der[0] = gkyl_malloc(sizeof(double[3]));
    tov_lapse_function_der(spacetime, 0.0, xp, ts->pos_c[1], ts->pos_c[2],
        h, h, h, lapse_der);

    printf("  d(lapse)/dx at r=R/2:\n");
    printf("    manual FD = %.6e\n", dalpha_dx_fd);
    printf("    from func = %.6e\n", lapse_der[0][0]);
    printf("    diff      = %.2e\n", lapse_der[0][0] - dalpha_dx_fd);

    gkyl_free(lapse_der[0]);
    gkyl_free(lapse_der);
}

struct gkyl_gr_spacetime*
gkyl_gr_tov_spacetime_new(bool use_gpu, struct gkyl_tov *tov,
  double pos_c_x, double pos_c_y, double pos_c_z)
{
  struct gr_tov_spacetime *ts = gkyl_malloc(sizeof(struct gr_tov_spacetime));
  ts->tov    = tov;
  ts->M_star = gkyl_tov_star_mass(tov);
  ts->R_star = gkyl_tov_star_radius(tov);
  ts->pos_c[0] = pos_c_x;
  ts->pos_c[1] = pos_c_y;
  ts->pos_c[2] = pos_c_z;

  ts->bh = gkyl_gr_blackhole_new(use_gpu, ts->M_star, 0.0, pos_c_x, pos_c_y, pos_c_z);

  ts->spacetime.lapse_function_func = tov_lapse_function;
  ts->spacetime.spatial_metric_tensor_func = tov_spatial_metric_tensor;
  ts->spacetime.shift_vector_func = tov_shift_vector;

  ts->spacetime.excision_region_func = tov_excision_region;

  ts->spacetime.spatial_metric_det_func = tov_spatial_metric_det;
  ts->spacetime.spatial_inv_metric_tensor_func = tov_spatial_inv_metric_tensor;
  ts->spacetime.spatial_metric_tensor_der_func  = tov_spatial_metric_tensor_der;
  ts->spacetime.lapse_function_der_func = tov_lapse_function_der;
  ts->spacetime.shift_vector_der_func = tov_shift_vector_der;
  ts->spacetime.spatial_christoffel_func = tov_spatial_christoffel;
  ts->spacetime.extrinsic_curvature_tensor_func = tov_extrinsic_curvature_tensor;

  //ts->spacetime.vacuum_einstein = NULL;
  tov_verify_metric(&ts->spacetime);
  tov_verify_lapse_derivative(&ts->spacetime);
  return &ts->spacetime;
}
