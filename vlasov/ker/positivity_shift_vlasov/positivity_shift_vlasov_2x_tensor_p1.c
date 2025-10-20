#include <gkyl_positivity_shift_vlasov_kernels.h> 

GKYL_CU_DH bool positivity_shift_vlasov_conf_pos_check_2x_tensor_p1(const double *fld) 
{ 
  // fld: conf-space field.

  bool is_positive = true;

  if (0.5*fld[3]-0.5*fld[2]-0.5*fld[1]+0.5*fld[0] < 0.) is_positive = false;
  if (-(0.5*fld[3])+0.5*fld[2]-0.5*fld[1]+0.5*fld[0] < 0.) is_positive = false;
  if (-(0.5*fld[3])-0.5*fld[2]+0.5*fld[1]+0.5*fld[0] < 0.) is_positive = false;
  if (0.5*fld[3]+0.5*fld[2]+0.5*fld[1]+0.5*fld[0] < 0.) is_positive = false;

  return is_positive;

}
