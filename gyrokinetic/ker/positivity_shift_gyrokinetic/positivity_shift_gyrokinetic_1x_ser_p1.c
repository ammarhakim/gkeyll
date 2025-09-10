#include <gkyl_positivity_shift_gyrokinetic_kernels.h> 

GKYL_CU_DH bool positivity_shift_gyrokinetic_conf_pos_check_1x_ser_p1(const double *fld) 
{ 
  // fld: conf-space field.

  bool is_positive = true;

  if (0.7071067811865475*fld[0]-0.7071067811865475*fld[1] < 0.) is_positive = false;
  if (0.7071067811865475*fld[1]+0.7071067811865475*fld[0] < 0.) is_positive = false;

  return is_positive;

}
