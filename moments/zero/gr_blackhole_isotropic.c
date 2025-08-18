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

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)) + ((z - pos_z) * (z - pos_z)));
  double theta = atan2(sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y))), (z - pos_z));
  double phi = atan2((y - pos_y), (x - pos_x));

  double **spatial_transformation_tensor = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_transformation_tensor[i] = gkyl_malloc(sizeof(double[3]));
  }

  spatial_transformation_tensor[0][0] = sin(theta) * cos(phi);
  spatial_transformation_tensor[0][1] = sin(theta) * sin(phi);
  spatial_transformation_tensor[0][2] = cos(theta);

  spatial_transformation_tensor[1][0] = (cos(theta) * cos(phi)) / rho;
  spatial_transformation_tensor[1][1] = (cos(theta) * sin(phi)) / rho;
  spatial_transformation_tensor[1][2] = -sin(theta) / rho;

  spatial_transformation_tensor[2][0] = -sin(phi) / (rho * sin(theta));
  spatial_transformation_tensor[2][1] = cos(phi) / (rho * sin(theta));
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

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)) + ((z - pos_z) * (z - pos_z)));
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
  spacetime_transformation_tensor[2][1] = (cos(theta) * cos(phi)) / rho;
  spacetime_transformation_tensor[2][2] = (cos(theta) * sin(phi)) / rho;
  spacetime_transformation_tensor[2][3] = -sin(theta) / rho;

  spacetime_transformation_tensor[3][0] = 0.0;
  spacetime_transformation_tensor[3][1] = -sin(phi) / (rho * sin(theta));
  spacetime_transformation_tensor[3][2] = cos(phi) / (rho * sin(theta));
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

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)) + ((z - pos_z) * (z - pos_z)));
  double theta = atan2(sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y))), (z - pos_z));
  double phi = atan2((y - pos_y), (x - pos_x));

  double **spatial_inv_transformation_tensor = gkyl_malloc(sizeof(double*[3]));
  for (int i = 0; i < 3; i++) {
    spatial_inv_transformation_tensor[i] = gkyl_malloc(sizeof(double[3]));
  }

  spatial_inv_transformation_tensor[0][0] = sin(theta) * cos(phi);
  spatial_inv_transformation_tensor[0][1] = rho * cos(theta) * cos(phi);
  spatial_inv_transformation_tensor[0][2] = -rho * sin(theta) * sin(phi);

  spatial_inv_transformation_tensor[1][0] = sin(theta) * sin(phi);
  spatial_inv_transformation_tensor[1][1] = rho * cos(theta) * sin(phi);
  spatial_inv_transformation_tensor[1][2] = rho * sin(theta) * cos(phi);

  spatial_inv_transformation_tensor[2][0] = cos(theta);
  spatial_inv_transformation_tensor[2][1] = -rho * sin(theta);
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

  double rho = sqrt(((x - pos_x) * (x - pos_x)) + ((y - pos_y) * (y - pos_y)) + ((z - pos_z) * (z - pos_z)));
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
  spacetime_inv_transformation_tensor[1][2] = rho * cos(theta) * cos(phi);
  spacetime_inv_transformation_tensor[1][3] = -rho * sin(theta) * sin(phi);

  spacetime_inv_transformation_tensor[2][0] = 0.0;
  spacetime_inv_transformation_tensor[2][1] = sin(theta) * sin(phi);
  spacetime_inv_transformation_tensor[2][2] = rho * cos(theta) * sin(phi);
  spacetime_inv_transformation_tensor[2][3] = rho * sin(theta) * cos(phi);

  spacetime_inv_transformation_tensor[3][0] = 0.0;
  spacetime_inv_transformation_tensor[3][1] = cos(theta);
  spacetime_inv_transformation_tensor[3][2] = -rho * sin(theta);
  spacetime_inv_transformation_tensor[3][3] = 0.0;

  return spacetime_inv_transformation_tensor;
}