#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_gr_blackhole_isotropic.h>
#include <gkyl_gr_spacetime_diff.h>

double**
blackhole_isotropic_spatial_transformation_tensor(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z)
{
  const struct gr_blackhole_isotropic *blackhole_isotropic = container_of(spacetime, struct gr_blackhole_isotropic, spacetime);

  double pos_x = blackhole_isotropic->pos_x;
  double pos_y = blackhole_isotropic->pos_y;
  double pos_z = blackhole_isotropic->pos_z;

  double eta = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)) + ((z - pos_z) * (z - pos_z)));
  double theta = atan2(sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y))), (z - pos_z));
  double phi = atan2((y - pos_y), (x - pos_x));

  double **spatial_transformation_tensor = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_transformation_tensor[i] = gkyl_malloc(sizeof(double[3]));
  }

  spatial_transformation_tensor[0][0] = sin(theta) * cos(phi);
  spatial_transformation_tensor[0][1] = sin(theta) * sin(phi);
  spatial_transformation_tensor[0][2] = cos(theta);

  spatial_transformation_tensor[1][0] = (cos(theta) * cos(phi)) / eta;
  spatial_transformation_tensor[1][1] = (cos(theta) * sin(phi)) / eta;
  spatial_transformation_tensor[1][2] = -sin(theta) / eta;

  spatial_transformation_tensor[2][0] = -sin(phi) / (eta * sin(theta));
  spatial_transformation_tensor[2][1] = cos(phi) / (eta * sin(theta));
  spatial_transformation_tensor[2][2] = 0.0;

  return spatial_transformation_tensor;
}

double**
blackhole_isotropic_spacetime_transformation_tensor(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z)
{
  const struct gr_blackhole_isotropic *blackhole_isotropic = container_of(spacetime, struct gr_blackhole_isotropic, spacetime);

  double pos_x = blackhole_isotropic->pos_x;
  double pos_y = blackhole_isotropic->pos_y;
  double pos_z = blackhole_isotropic->pos_z;

  double eta = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)) + ((z - pos_z) * (z - pos_z)));
  double theta = atan2(sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y))), (z - pos_z));
  double phi = atan2((y - pos_y), (x - pos_x));

  double **spacetime_transformation_tensor = gkyl_malloc(sizeof(double*[4]));
  for (int i = 0; i < 4; i++) {
    spacetime_transformation_tensor[i] = gkyl_malloc(sizeof(double[4]));
  }

  spacetime_transformation_tensor[0][0] = 1.0;
  spacetime_transformation_tensor[0][1] = 0.0;
  spacetime_transformation_tensor[0][2] = 0.0;
  spacetime_transformation_tensor[0][3] = 0.0;
  
  spacetime_transformation_tensor[1][0] = 0.0;
  spacetime_transformation_tensor[1][1] = sin(theta) * cos(phi);
  spacetime_transformation_tensor[1][2] = sin(theta) * sin(phi);
  spacetime_transformation_tensor[1][3] = cos(theta);

  spacetime_transformation_tensor[2][0] = 0.0;
  spacetime_transformation_tensor[2][1] = (cos(theta) * cos(phi)) / eta;
  spacetime_transformation_tensor[2][2] = (cos(theta) * sin(phi)) / eta;
  spacetime_transformation_tensor[2][3] = -sin(theta) / eta;

  spacetime_transformation_tensor[3][0] = 0.0;
  spacetime_transformation_tensor[3][1] = -sin(phi) / (eta * sin(theta));
  spacetime_transformation_tensor[3][2] = cos(phi) / (eta * sin(theta));
  spacetime_transformation_tensor[3][3] = 0.0;

  return spacetime_transformation_tensor;
}

double**
blackhole_isotropic_spatial_inv_transformation_tensor(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z)
{
  const struct gr_blackhole_isotropic *blackhole_isotropic = container_of(spacetime, struct gr_blackhole_isotropic, spacetime);

  double pos_x = blackhole_isotropic->pos_x;
  double pos_y = blackhole_isotropic->pos_y;
  double pos_z = blackhole_isotropic->pos_z;

  double eta = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)) + ((z - pos_z) * (z - pos_z)));
  double theta = atan2(sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y))), (z - pos_z));
  double phi = atan2((y - pos_y), (x - pos_x));

  double **spatial_inv_transformation_tensor = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_inv_transformation_tensor[i] = gkyl_malloc(sizeof(double[3]));
  }

  spatial_inv_transformation_tensor[0][0] = sin(theta) * cos(phi);
  spatial_inv_transformation_tensor[0][1] = eta * cos(theta) * cos(phi);
  spatial_inv_transformation_tensor[0][2] = -eta * sin(theta) * sin(phi);

  spatial_inv_transformation_tensor[1][0] = sin(theta) * sin(phi);
  spatial_inv_transformation_tensor[1][1] = eta * cos(theta) * sin(phi);
  spatial_inv_transformation_tensor[1][2] = eta * sin(theta) * cos(phi);

  spatial_inv_transformation_tensor[2][0] = cos(theta);
  spatial_inv_transformation_tensor[2][1] = -eta * sin(theta);
  spatial_inv_transformation_tensor[2][2] = 0.0;

  return spatial_inv_transformation_tensor;
}

double**
blackhole_isotropic_spacetime_inv_transformation_tensor(const struct gkyl_gr_spacetime* spacetime, const double x, const double y, const double z)
{
  const struct gr_blackhole_isotropic *blackhole_isotropic = container_of(spacetime, struct gr_blackhole_isotropic, spacetime);

  double pos_x = blackhole_isotropic->pos_x;
  double pos_y = blackhole_isotropic->pos_y;
  double pos_z = blackhole_isotropic->pos_z;

  double eta = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)) + ((z - pos_z) * (z - pos_z)));
  double theta = atan2(sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y))), (z - pos_z));
  double phi = atan2((y - pos_y), (x - pos_x));

  double **spacetime_inv_transformation_tensor = gkyl_malloc(sizeof(double*[4]));
  for (int i = 0; i < 4; i++) {
    spacetime_inv_transformation_tensor[i] = gkyl_malloc(sizeof(double[4]));
  }

  spacetime_inv_transformation_tensor[0][0] = 1.0;
  spacetime_inv_transformation_tensor[0][1] = 0.0;
  spacetime_inv_transformation_tensor[0][2] = 0.0;
  spacetime_inv_transformation_tensor[0][3] = 0.0;
  
  spacetime_inv_transformation_tensor[1][0] = 0.0;
  spacetime_inv_transformation_tensor[1][1] = sin(theta) * cos(phi);
  spacetime_inv_transformation_tensor[1][2] = eta * cos(theta) * cos(phi);
  spacetime_inv_transformation_tensor[1][3] = -eta * sin(theta) * sin(phi);

  spacetime_inv_transformation_tensor[2][0] = 0.0;
  spacetime_inv_transformation_tensor[2][1] = sin(theta) * sin(phi);
  spacetime_inv_transformation_tensor[2][2] = eta * cos(theta) * sin(phi);
  spacetime_inv_transformation_tensor[2][3] = eta * sin(theta) * cos(phi);

  spacetime_inv_transformation_tensor[3][0] = 0.0;
  spacetime_inv_transformation_tensor[3][1] = cos(theta);
  spacetime_inv_transformation_tensor[3][2] = -eta * sin(theta);
  spacetime_inv_transformation_tensor[3][3] = 0.0;

  return spacetime_inv_transformation_tensor;
}

static void
blackhole_isotropic_spatial_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spatial_metric_tensor)
{
  const struct gr_blackhole_isotropic *blackhole_isotropic = container_of(spacetime, struct gr_blackhole_isotropic, spacetime);

  double mass = blackhole_isotropic->mass;
  double spin = blackhole_isotropic->spin;

  double pos_x = blackhole_isotropic->pos_x;
  double pos_y = blackhole_isotropic->pos_y;
  double pos_z = blackhole_isotropic->pos_z;

  double Mspin = mass * spin;

  double eta = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)) + ((z - pos_z) * (z - pos_z)));
  double r_plus = mass + sqrt((mass * mass) - (Mspin * Mspin));
  double r_minus = mass - sqrt((mass * mass) - (Mspin * Mspin));
  double rho = (r_plus * (cosh(0.5 * eta) * cosh(0.5 * eta))) - (r_minus * (sinh(0.5 * eta) * sinh(0.5 * eta)));

  double theta = atan2(sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y))), (z - pos_z));
  
  double rho_kerr_sq = (rho * rho) + ((Mspin * Mspin) * (cos(theta) * cos(theta)));
  double delta = (rho * rho) - (2.0 * mass * rho) + (Mspin * Mspin);

  double **spatial_metric_tensor_spherical = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric_tensor_spherical[i] = gkyl_malloc(sizeof(double[3]));
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      spatial_metric_tensor_spherical[i][j] = 0.0;
      (*spatial_metric_tensor)[i][j] = 0.0;
    }
  }

  spatial_metric_tensor_spherical[0][0] = rho_kerr_sq;
  spatial_metric_tensor_spherical[1][1] = rho_kerr_sq;
  spatial_metric_tensor_spherical[2][2] = ((((rho * rho) + (Mspin * Mspin)) * ((rho * rho) + (Mspin * Mspin))) -
    (delta * (Mspin * Mspin) * (sin(theta) * sin(theta)))) * (sin(theta) * sin(theta)) / rho_kerr_sq;
  
  double **spatial_transformation_tensor = blackhole_isotropic_spatial_transformation_tensor(spacetime, x, y, z);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          (*spatial_metric_tensor)[i][j] += spatial_transformation_tensor[k][i] * spatial_transformation_tensor[l][j] * spatial_metric_tensor_spherical[k][l];
        }
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric_tensor_spherical[i]);
    gkyl_free(spatial_transformation_tensor[i]);
  }
  gkyl_free(spatial_metric_tensor_spherical);
  gkyl_free(spatial_transformation_tensor);
}

static void
blackhole_isotropic_spacetime_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spacetime_metric_tensor)
{
  const struct gr_blackhole_isotropic *blackhole_isotropic = container_of(spacetime, struct gr_blackhole_isotropic, spacetime);

  double mass = blackhole_isotropic->mass;
  double spin = blackhole_isotropic->spin;

  double pos_x = blackhole_isotropic->pos_x;
  double pos_y = blackhole_isotropic->pos_y;
  double pos_z = blackhole_isotropic->pos_z;

  double Mspin = mass * spin;

  double eta = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)) + ((z - pos_z) * (z - pos_z)));
  double r_plus = mass + sqrt((mass * mass) - (Mspin * Mspin));
  double r_minus = mass - sqrt((mass * mass) - (Mspin * Mspin));
  double rho = (r_plus * (cosh(0.5 * eta) * cosh(0.5 * eta))) - (r_minus * (sinh(0.5 * eta) * sinh(0.5 * eta)));

  double theta = atan2(sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y))), (z - pos_z));
  
  double rho_kerr_sq = (rho * rho) + ((Mspin * Mspin) * (cos(theta) * cos(theta)));
  double delta = (rho * rho) - (2.0 * mass * rho) + (Mspin * Mspin);

  double **spacetime_metric_tensor_spherical = gkyl_malloc(sizeof(double*[4]));
  for (int i = 0; i < 4; i++) {
    spacetime_metric_tensor_spherical[i] = gkyl_malloc(sizeof(double[4]));
  }

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      spacetime_metric_tensor_spherical[i][j] = 0.0;
      (*spacetime_metric_tensor)[i][j] = 0.0;
    }
  }

  spacetime_metric_tensor_spherical[0][0] = -(1.0 / (1.0 + (((2.0 * mass * rho) * ((Mspin * Mspin) + (rho * rho))) / (delta * rho_kerr_sq)))) +
    ((4.0 * (Mspin * Mspin) * (mass * mass) * (rho * rho) * (sin(theta) * sin(theta))) / ((rho_kerr_sq * ((Mspin * Mspin) + (rho * rho)) * ((Mspin * Mspin) + (rho * rho))) -
      ((Mspin * Mspin) * delta * (sin(theta) * sin(theta)))));
  spacetime_metric_tensor_spherical[1][1] = rho_kerr_sq;
  spacetime_metric_tensor_spherical[2][2] = rho_kerr_sq;
  spacetime_metric_tensor_spherical[3][3] = ((((rho * rho) + (Mspin * Mspin)) * ((rho * rho) + (Mspin * Mspin))) - (delta * (Mspin * Mspin) * (sin(theta) * sin(theta)))) *
    (sin(theta) * sin(theta)) / rho_kerr_sq;

  spacetime_metric_tensor_spherical[0][3] = -(2.0 * Mspin * mass * rho * (sin(theta) * sin(theta))) / rho_kerr_sq;
  spacetime_metric_tensor_spherical[3][0] = -(2.0 * Mspin * mass * rho * (sin(theta) * sin(theta))) / rho_kerr_sq;
  
  double **spacetime_transformation_tensor = blackhole_isotropic_spacetime_transformation_tensor(spacetime, x, y, z);

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        for (int l = 0; l < 4; l++) {
          (*spacetime_metric_tensor)[i][j] += spacetime_transformation_tensor[k][i] * spacetime_transformation_tensor[l][j] * spacetime_metric_tensor_spherical[k][l];
        }
      }
    }
  }

  for (int i = 0; i < 4; i++) {
    gkyl_free(spacetime_metric_tensor_spherical[i]);
    gkyl_free(spacetime_transformation_tensor[i]);
  }
  gkyl_free(spacetime_metric_tensor_spherical);
  gkyl_free(spacetime_transformation_tensor);
}

static void
blackhole_isotropic_spatial_inv_metric_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double*** spatial_inv_metric_tensor)
{
  const struct gr_blackhole_isotropic *blackhole_isotropic = container_of(spacetime, struct gr_blackhole_isotropic, spacetime);

  double mass = blackhole_isotropic->mass;
  double spin = blackhole_isotropic->spin;

  double pos_x = blackhole_isotropic->pos_x;
  double pos_y = blackhole_isotropic->pos_y;
  double pos_z = blackhole_isotropic->pos_z;

  double Mspin = mass * spin;

  double eta = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)) + ((z - pos_z) * (z - pos_z)));
  double r_plus = mass + sqrt((mass * mass) - (Mspin * Mspin));
  double r_minus = mass - sqrt((mass * mass) - (Mspin * Mspin));
  double rho = (r_plus * (cosh(0.5 * eta) * cosh(0.5 * eta))) - (r_minus * (sinh(0.5 * eta) * sinh(0.5 * eta)));

  double theta = atan2(sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y))), (z - pos_z));
  
  double rho_kerr_sq = (rho * rho) + ((Mspin * Mspin) * (cos(theta) * cos(theta)));
  double delta = (rho * rho) - (2.0 * mass * rho) + (Mspin * Mspin);

  double **spatial_inv_metric_tensor_spherical = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_inv_metric_tensor_spherical[i] = gkyl_malloc(sizeof(double[3]));
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      spatial_inv_metric_tensor_spherical[i][j] = 0.0;
      (*spatial_inv_metric_tensor)[i][j] = 0.0;
    }
  }

  spatial_inv_metric_tensor_spherical[0][0] = 1.0 / rho_kerr_sq;
  spatial_inv_metric_tensor_spherical[1][1] = 1.0 / rho_kerr_sq;
  spatial_inv_metric_tensor_spherical[2][2] = rho_kerr_sq / (((((rho * rho) + (Mspin * Mspin)) * ((rho * rho) + (Mspin * Mspin))) - (delta * (Mspin * Mspin) *
    (sin(theta) * sin(theta)))) * (sin(theta) * sin(theta)));
  
  double **spatial_inv_transformation_tensor = blackhole_isotropic_spatial_inv_transformation_tensor(spacetime, x, y, z);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          (*spatial_inv_metric_tensor)[i][j] += spatial_inv_transformation_tensor[k][i] * spatial_inv_transformation_tensor[l][j] * spatial_inv_metric_tensor_spherical[k][l];
        }
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_inv_metric_tensor_spherical[i]);
    gkyl_free(spatial_inv_transformation_tensor[i]);
  }
  gkyl_free(spatial_inv_metric_tensor_spherical);
  gkyl_free(spatial_inv_transformation_tensor);
}

static void
blackhole_isotropic_spatial_metric_det(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* spatial_metric_det)
{
  double** spatial_metric = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_metric[i] = gkyl_malloc(sizeof(double[3]));
  }

  blackhole_isotropic_spatial_metric_tensor(spacetime, t, x, y, z, &spatial_metric);

  *spatial_metric_det = (spatial_metric[0][0] * ((spatial_metric[1][1] * spatial_metric[2][2]) - (spatial_metric[2][1] * spatial_metric[1][2]))) -
    (spatial_metric[0][1] * ((spatial_metric[1][0] * spatial_metric[2][2]) - (spatial_metric[1][2] * spatial_metric[2][0]))) +
    (spatial_metric[0][2] * ((spatial_metric[1][0] * spatial_metric[2][1]) - (spatial_metric[1][1] * spatial_metric[2][0])));
  
  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_metric[i]);
  }
  gkyl_free(spatial_metric);
}

static void
blackhole_isotropic_spacetime_metric_det(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* spacetime_metric_det)
{
  double spatial_metric_det;
  double lapse_function;
  blackhole_isotropic_spatial_metric_det(spacetime, t, x, y, z, &spatial_metric_det);
  blackhole_isotropic_lapse_function(spacetime, t, x, y, z, &lapse_function);

  *spacetime_metric_det = - (lapse_function * lapse_function) * spatial_metric_det;
}

static void
blackhole_isotropic_spatial_metric_tensor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
   const double dx, const double dy, const double dz, double**** spatial_metric_tensor_der)
{
  gkyl_gr_spatial_metric_tensor_diff(spacetime, t, x, y, z, dx, dy, dz, spatial_metric_tensor_der);
}

static void
blackhole_isotropic_spacetime_metric_tensor_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double**** spacetime_metric_tensor_der)
{
  gkyl_gr_spacetime_metric_tensor_diff(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_metric_tensor_der);
}

static void
blackhole_isotropic_lapse_function(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double* lapse_function)
{
  const struct gr_blackhole_isotropic *blackhole_isotropic = container_of(spacetime, struct gr_blackhole_isotropic, spacetime);

  double mass = blackhole_isotropic->mass;
  double spin = blackhole_isotropic->spin;

  double pos_x = blackhole_isotropic->pos_x;
  double pos_y = blackhole_isotropic->pos_y;
  double pos_z = blackhole_isotropic->pos_z;

  double Mspin = mass * spin;

  double eta = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)) + ((z - pos_z) * (z - pos_z)));
  double r_plus = mass + sqrt((mass * mass) - (Mspin * Mspin));
  double r_minus = mass - sqrt((mass * mass) - (Mspin * Mspin));
  double rho = (r_plus * (cosh(0.5 * eta) * cosh(0.5 * eta))) - (r_minus * (sinh(0.5 * eta) * sinh(0.5 * eta)));

  double theta = atan2(sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y))), (z - pos_z));
  
  double rho_kerr_sq = (rho * rho) + ((Mspin * Mspin) * (cos(theta) * cos(theta)));
  double delta = (rho * rho) - (2.0 * mass * rho) + (Mspin * Mspin);

  *lapse_function = sqrt(1.0 / (1.0 + ((2.0 * mass * rho) * ((rho * rho) + (Mspin * Mspin)) / (delta * rho_kerr_sq))));
}

static void
blackhole_isotropic_shift_vector(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  double** shift_vector)
{
  const struct gr_blackhole_isotropic *blackhole_isotropic = container_of(spacetime, struct gr_blackhole_isotropic, spacetime);

  double mass = blackhole_isotropic->mass;
  double spin = blackhole_isotropic->spin;

  double pos_x = blackhole_isotropic->pos_x;
  double pos_y = blackhole_isotropic->pos_y;
  double pos_z = blackhole_isotropic->pos_z;

  double Mspin = mass * spin;

  double eta = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)) + ((z - pos_z) * (z - pos_z)));
  double r_plus = mass + sqrt((mass * mass) - (Mspin * Mspin));
  double r_minus = mass - sqrt((mass * mass) - (Mspin * Mspin));
  double rho = (r_plus * (cosh(0.5 * eta) * cosh(0.5 * eta))) - (r_minus * (sinh(0.5 * eta) * sinh(0.5 * eta)));

  double theta = atan2(sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y))), (z - pos_z));
  
  double rho_kerr_sq = (rho * rho) + ((Mspin * Mspin) * (cos(theta) * cos(theta)));
  double delta = (rho * rho) - (2.0 * mass * rho) + (Mspin * Mspin);

  double *shift_vector_spherical = gkyl_malloc(sizeof(double[3]));
  for (int i = 0; i < 3; i++) {
    (*shift_vector)[i] = 0.0;
  }

  shift_vector_spherical[0] = 0.0;
  shift_vector_spherical[1] = 0.0;
  shift_vector_spherical[2] = -(2.0 * Mspin * mass * rho) / ((((rho * rho) + (Mspin * Mspin)) * ((rho * rho) + (Mspin * Mspin))) - ((Mspin * Mspin) * delta * (sin(theta) * sin(theta))));

  double **spatial_inv_transformation_tensor = blackhole_isotropic_spatial_inv_transformation_tensor(spacetime, x, y, z);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      (*shift_vector)[i] += spatial_inv_transformation_tensor[i][j] * shift_vector_spherical[j];
    }
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(spatial_inv_transformation_tensor[i]);
  }
  gkyl_free(shift_vector_spherical);
  gkyl_free(spatial_inv_transformation_tensor);
}

static void
blackhole_isotropic_lapse_function_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double** lapse_function_der)
{
  gkyl_gr_lapse_function_diff(spacetime, t, x, y, z, dx, dy, dz, lapse_function_der);
}

static void
blackhole_isotropic_shift_vector_der(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** shift_vector_der)
{
  gkyl_gr_shift_vector_diff(spacetime, t, x, y, z, dx, dy, dz, shift_vector_der);
}

static void
blackhole_isotropic_spatial_christoffel(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double**** spatial_christoffel)
{
  gkyl_gr_spatial_christoffel_fd(spacetime, t, x, y, z, dx, dy, dz, spatial_christoffel);
}

static void
blackhole_isotropic_spacetime_christoffel(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double**** spacetime_christoffel)
{
  gkyl_gr_spacetime_christoffel_fd(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_christoffel);
}

static void
blackhole_isotropic_spatial_riemann_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double***** spatial_riemann_tensor)
{
  gkyl_gr_spatial_riemann_tensor_fd(spacetime, t, x, y, z, dx, dy, dz, spatial_riemann_tensor);
}

static void
blackhole_isotropic_spacetime_riemann_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double***** spacetime_riemann_tensor)
{
  gkyl_gr_spacetime_riemann_tensor_fd(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_riemann_tensor);
}

static void
blackhole_isotropic_spatial_ricci_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** spatial_ricci_tensor)
{
  gkyl_gr_spatial_ricci_tensor_fd(spacetime, t, x, y, z, dx, dy, dz, spatial_ricci_tensor);
}

static void
blackhole_isotropic_spacetime_ricci_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double*** spacetime_ricci_tensor)
{
  gkyl_gr_spacetime_ricci_tensor_fd(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_ricci_tensor);
}

static void
blackhole_isotropic_spatial_ricci_scalar(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double* spatial_ricci_scalar)
{
  gkyl_gr_spatial_ricci_scalar_fd(spacetime, t, x, y, z, dx, dy, dz, spatial_ricci_scalar);
}

static void
blackhole_isotropic_spacetime_ricci_scalar(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double* spacetime_ricci_scalar)
{
  gkyl_gr_spacetime_ricci_scalar_fd(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_ricci_scalar);
}

static void
blackhole_isotropic_spatial_weyl_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double***** spatial_weyl_tensor)
{
  gkyl_gr_spatial_weyl_tensor_fd(spacetime, t, x, y, z, dx, dy, dx, spatial_weyl_tensor);
}

static void
blackhole_isotropic_spacetime_weyl_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dt, const double dx, const double dy, const double dz, double***** spacetime_weyl_tensor)
{
  gkyl_gr_spacetime_weyl_tensor_fd(spacetime, t, x, y, z, dt, dx, dy, dz, spacetime_weyl_tensor);
}

static void
blackhole_isotropic_extrinsic_curvature_tensor(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  const double dx, const double dy, const double dz, double*** extrinsic_curvature_tensor)
{
  const struct gr_blackhole_isotropic *blackhole_isotropic = container_of(spacetime, struct gr_blackhole_isotropic, spacetime);

  double mass = blackhole_isotropic->mass;
  double spin = blackhole_isotropic->spin;

  double pos_x = blackhole_isotropic->pos_x;
  double pos_y = blackhole_isotropic->pos_y;
  double pos_z = blackhole_isotropic->pos_z;

  double Mspin = mass * spin;

  double eta = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)) + ((z - pos_z) * (z - pos_z)));
  double r_plus = mass + sqrt((mass * mass) - (Mspin * Mspin));
  double r_minus = mass - sqrt((mass * mass) - (Mspin * Mspin));
  double rho = (r_plus * (cosh(0.5 * eta) * cosh(0.5 * eta))) - (r_minus * (sinh(0.5 * eta) * sinh(0.5 * eta)));
  
  double theta = atan2(sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y))), (z - pos_z));
  
  double rho_kerr_sq = (rho * rho) + ((Mspin * Mspin) * (cos(theta) * cos(theta)));
  double delta = (rho * rho) - (2.0 * mass * rho) + (Mspin * Mspin);

  double **extrinsic_curvature_tensor_spherical = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    extrinsic_curvature_tensor_spherical[i] = gkyl_malloc(sizeof(double[3]));
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      extrinsic_curvature_tensor_spherical[i][j] = 0.0;
      (*extrinsic_curvature_tensor)[i][j] = 0.0;
    }
  }

  double conformal_factor_sq = sin(theta) / sqrt(((((rho * rho) + (Mspin * Mspin)) * ((rho * rho) + (Mspin * Mspin))) -
    (delta * (Mspin * Mspin) * (sin(theta) * sin(theta)))) * (sin(theta) * sin(theta)) / rho_kerr_sq);

  extrinsic_curvature_tensor_spherical[0][2] = conformal_factor_sq * ((Mspin * mass * ((2.0 * (rho * rho) * ((rho * rho) + (Mspin * Mspin))) +
    (rho_kerr_sq * ((rho * rho) - (Mspin * Mspin)))) * (sin(theta) * sin(theta))) / (rho_kerr_sq * rho_kerr_sq));
  extrinsic_curvature_tensor_spherical[2][0] = conformal_factor_sq * ((Mspin * mass * ((2.0 * (rho * rho) * ((rho * rho) + (Mspin * Mspin))) +
    (rho_kerr_sq * ((rho * rho) - (Mspin * Mspin)))) * (sin(theta) * sin(theta))) / (rho_kerr_sq * rho_kerr_sq));

  extrinsic_curvature_tensor_spherical[1][2] = -conformal_factor_sq * ((2.0 * (Mspin * Mspin * Mspin) * mass * rho * sqrt(delta) * cos(theta) *
    (sin(theta) * sin(theta) * sin(theta))) / (rho_kerr_sq * rho_kerr_sq));
  extrinsic_curvature_tensor_spherical[2][1] = -conformal_factor_sq * ((2.0 * (Mspin * Mspin * Mspin) * mass * rho * sqrt(delta) * cos(theta) *
    (sin(theta) * sin(theta) * sin(theta))) / (rho_kerr_sq * rho_kerr_sq));

  double **spatial_transformation_tensor = blackhole_isotropic_spatial_transformation_tensor(spacetime, x, y, z);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          (*extrinsic_curvature_tensor)[i][j] += spatial_transformation_tensor[k][i] * spatial_transformation_tensor[l][j] * extrinsic_curvature_tensor_spherical[k][l];
        }
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    gkyl_free(extrinsic_curvature_tensor_spherical[i]);
    gkyl_free(spatial_transformation_tensor[i]);
  }
  gkyl_free(extrinsic_curvature_tensor_spherical);
  gkyl_free(spatial_transformation_tensor);
}

static void
blackhole_isotropic_excision_region(const struct gkyl_gr_spacetime* spacetime, const double t, const double x, const double y, const double z,
  bool* in_excision_region)
{
  const struct gr_blackhole_isotropic *blackhole_isotropic = container_of(spacetime, struct gr_blackhole_isotropic, spacetime);

  double mass = blackhole_isotropic->mass;

  double pos_x = blackhole_isotropic->pos_x;
  double pos_y = blackhole_isotropic->pos_y;
  double pos_z = blackhole_isotropic->pos_z;

  double eta = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)) + ((z - pos_z) * (z - pos_z)));

  if (eta <= mass * 2.0) {
    *in_excision_region = true;
  }
  else {
    *in_excision_region = false;
  }
}

void
gkyl_gr_blackhole_isotropic_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_gr_spacetime* base = container_of(ref, struct gkyl_gr_spacetime, ref_count);

  if (gkyl_gr_spacetime_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct gr_blackhole_isotropic *gr_blackhole_isotropic = container_of(base->on_dev, struct gr_blackhole_isotropic, spacetime);
    gkyl_cu_free(gr_blackhole_isotropic);
  }

  struct gr_blackhole_isotropic *gr_blackhole_isotropic = container_of(base, struct gr_blackhole_isotropic, spacetime);
  gkyl_free(gr_blackhole_isotropic);
}

struct gkyl_gr_spacetime*
gkyl_gr_blackhole_isotropic_new(bool use_gpu, double mass, double spin, double pos_x, double pos_y, double pos_z)
{
  return gkyl_gr_blackhole_isotropic_inew(&(struct gkyl_gr_blackhole_isotropic_inp) {
      .use_gpu = use_gpu,
      .mass = mass,
      .spin = spin,
      .pos_x = pos_x,
      .pos_y = pos_y,
      .pos_z = pos_z,
    }
  );
}

struct gkyl_gr_spacetime*
gkyl_gr_blackhole_isotropic_inew(const struct gkyl_gr_blackhole_isotropic_inp* inp)
{
  struct gr_blackhole_isotropic *gr_blackhole_isotropic = gkyl_malloc(sizeof(struct gr_blackhole_isotropic));

  gr_blackhole_isotropic->mass = inp->mass;
  gr_blackhole_isotropic->spin = inp->spin;

  gr_blackhole_isotropic->pos_x = inp->pos_x;
  gr_blackhole_isotropic->pos_y = inp->pos_y;
  gr_blackhole_isotropic->pos_z = inp->pos_z;

  gr_blackhole_isotropic->spacetime.spatial_metric_tensor_func = blackhole_isotropic_spatial_metric_tensor;
  gr_blackhole_isotropic->spacetime.spacetime_metric_tensor_func = blackhole_isotropic_spacetime_metric_tensor;

  gr_blackhole_isotropic->spacetime.spatial_inv_metric_tensor_func = blackhole_isotropic_spatial_inv_metric_tensor;
  //gr_blackhole_isotropic->spacetime.spacetime_inv_metric_tensor_func = blackhole_isotropic_spacetime_inv_metric_tensor;

  gr_blackhole_isotropic->spacetime.spatial_metric_det_func = blackhole_isotropic_spatial_metric_det;
  gr_blackhole_isotropic->spacetime.spacetime_metric_det_func = blackhole_isotropic_spacetime_metric_det;

  gr_blackhole_isotropic->spacetime.spatial_metric_tensor_der_func = blackhole_isotropic_spatial_metric_tensor_der;
  gr_blackhole_isotropic->spacetime.spacetime_metric_tensor_der_func = blackhole_isotropic_spacetime_metric_tensor_der;

  gr_blackhole_isotropic->spacetime.lapse_function_func = blackhole_isotropic_lapse_function;
  gr_blackhole_isotropic->spacetime.shift_vector_func = blackhole_isotropic_shift_vector;

  gr_blackhole_isotropic->spacetime.lapse_function_der_func = blackhole_isotropic_lapse_function_der;
  gr_blackhole_isotropic->spacetime.shift_vector_der_func = blackhole_isotropic_shift_vector_der;

  gr_blackhole_isotropic->spacetime.spatial_christoffel_func = blackhole_isotropic_spatial_christoffel;
  gr_blackhole_isotropic->spacetime.spacetime_christoffel_func = blackhole_isotropic_spacetime_christoffel;

  gr_blackhole_isotropic->spacetime.spatial_riemann_tensor_func = blackhole_isotropic_spatial_riemann_tensor;
  gr_blackhole_isotropic->spacetime.spacetime_riemann_tensor_func = blackhole_isotropic_spacetime_riemann_tensor;

  gr_blackhole_isotropic->spacetime.spatial_ricci_tensor_func = blackhole_isotropic_spatial_ricci_tensor;
  gr_blackhole_isotropic->spacetime.spacetime_ricci_tensor_func = blackhole_isotropic_spacetime_ricci_tensor;

  gr_blackhole_isotropic->spacetime.spatial_ricci_scalar_func = blackhole_isotropic_spatial_ricci_scalar;
  gr_blackhole_isotropic->spacetime.spacetime_ricci_scalar_func = blackhole_isotropic_spacetime_ricci_scalar;

  gr_blackhole_isotropic->spacetime.spatial_weyl_tensor_func = blackhole_isotropic_spatial_weyl_tensor;
  gr_blackhole_isotropic->spacetime.spacetime_weyl_tensor_func = blackhole_isotropic_spacetime_weyl_tensor;

  gr_blackhole_isotropic->spacetime.extrinsic_curvature_tensor_func = blackhole_isotropic_extrinsic_curvature_tensor;

  gr_blackhole_isotropic->spacetime.excision_region_func = blackhole_isotropic_excision_region;

  gr_blackhole_isotropic->spacetime.flags = 0;
  GKYL_CLEAR_CU_ALLOC(gr_blackhole_isotropic->spacetime.flags);
  gr_blackhole_isotropic->spacetime.ref_count = gkyl_ref_count_init(gkyl_gr_blackhole_isotropic_free);
  gr_blackhole_isotropic->spacetime.on_dev = &gr_blackhole_isotropic->spacetime; // On the CPU, the spacetime object points to itself.

  return &gr_blackhole_isotropic->spacetime;
}