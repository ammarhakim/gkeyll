#include <gkyl_array_integrate_kernels.h>
GKYL_CU_DH void gkyl_array_integrate_op_sq_weighted_1x2v_gkhyb_p1(double *dxSq, double vol, int num_comp, int num_basis, const double *weight, const double *fIn, double *out) 
{ 

  const double volFac = 0.04714045207910316*vol;

  out[0] += (15.0*weight[0]*(fIn[11]*fIn[11])+30.000000000000004*weight[1]*fIn[10]*fIn[11]+15.0*weight[0]*(fIn[10]*fIn[10])+15.0*weight[0]*(fIn[9]*fIn[9])+30.000000000000004*weight[1]*fIn[8]*fIn[9]+15.0*weight[0]*(fIn[8]*fIn[8])+15.0*weight[0]*(fIn[7]*fIn[7])+30.0*weight[1]*fIn[6]*fIn[7]+15.0*weight[0]*(fIn[6]*fIn[6])+15.0*weight[0]*(fIn[5]*fIn[5])+30.0*weight[1]*fIn[3]*fIn[5]+15.0*weight[0]*(fIn[4]*fIn[4])+30.0*weight[1]*fIn[2]*fIn[4]+15.0*weight[0]*(fIn[3]*fIn[3])+15.0*weight[0]*(fIn[2]*fIn[2])+30.0*fIn[0]*fIn[1]*weight[1]+15.0*weight[0]*(fIn[1]*fIn[1])+15.0*(fIn[0]*fIn[0])*weight[0])*volFac;
}

