#include <gkyl_positivity_shift_vlasov_kernels.h> 
#include <math.h> 
#include <float.h> 

GKYL_CU_DH bool positivity_shift_vlasov_shift_only_1x3v_tensor_p1(double ffloor, double *distf) 
{ 
  // ffloor: Distribution function floor to shift to when f<0.
  // distf: distribution function.

  bool shifted = false;

  double fnod[16];
  fnod[0] = 0.25*distf[15]-0.25*distf[14]-0.25*distf[13]-0.25*distf[12]-0.25*distf[11]+0.25*distf[10]+0.25*distf[9]+0.25*distf[8]+0.25*distf[7]+0.25*distf[6]+0.25*distf[5]-0.25*distf[4]-0.25*distf[3]-0.25*distf[2]-0.25*distf[1]+0.25*distf[0]; 
  fnod[1] = -(0.25*distf[15])+0.25*distf[14]+0.25*distf[13]+0.25*distf[12]-0.25*distf[11]-0.25*distf[10]-0.25*distf[9]-0.25*distf[8]+0.25*distf[7]+0.25*distf[6]+0.25*distf[5]+0.25*distf[4]-0.25*distf[3]-0.25*distf[2]-0.25*distf[1]+0.25*distf[0]; 
  fnod[2] = -(0.25*distf[15])+0.25*distf[14]+0.25*distf[13]-0.25*distf[12]+0.25*distf[11]-0.25*distf[10]+0.25*distf[9]+0.25*distf[8]-0.25*distf[7]-0.25*distf[6]+0.25*distf[5]-0.25*distf[4]+0.25*distf[3]-0.25*distf[2]-0.25*distf[1]+0.25*distf[0]; 
  fnod[3] = 0.25*distf[15]-0.25*distf[14]-0.25*distf[13]+0.25*distf[12]+0.25*distf[11]+0.25*distf[10]-0.25*distf[9]-0.25*distf[8]-0.25*distf[7]-0.25*distf[6]+0.25*distf[5]+0.25*distf[4]+0.25*distf[3]-0.25*distf[2]-0.25*distf[1]+0.25*distf[0]; 
  fnod[4] = -(0.25*distf[15])+0.25*distf[14]-0.25*distf[13]+0.25*distf[12]+0.25*distf[11]+0.25*distf[10]-0.25*distf[9]+0.25*distf[8]-0.25*distf[7]+0.25*distf[6]-0.25*distf[5]-0.25*distf[4]-0.25*distf[3]+0.25*distf[2]-0.25*distf[1]+0.25*distf[0]; 
  fnod[5] = 0.25*distf[15]-0.25*distf[14]+0.25*distf[13]-0.25*distf[12]+0.25*distf[11]-0.25*distf[10]+0.25*distf[9]-0.25*distf[8]-0.25*distf[7]+0.25*distf[6]-0.25*distf[5]+0.25*distf[4]-0.25*distf[3]+0.25*distf[2]-0.25*distf[1]+0.25*distf[0]; 
  fnod[6] = 0.25*distf[15]-0.25*distf[14]+0.25*distf[13]+0.25*distf[12]-0.25*distf[11]-0.25*distf[10]-0.25*distf[9]+0.25*distf[8]+0.25*distf[7]-0.25*distf[6]-0.25*distf[5]-0.25*distf[4]+0.25*distf[3]+0.25*distf[2]-0.25*distf[1]+0.25*distf[0]; 
  fnod[7] = -(0.25*distf[15])+0.25*distf[14]-0.25*distf[13]-0.25*distf[12]-0.25*distf[11]+0.25*distf[10]+0.25*distf[9]-0.25*distf[8]+0.25*distf[7]-0.25*distf[6]-0.25*distf[5]+0.25*distf[4]+0.25*distf[3]+0.25*distf[2]-0.25*distf[1]+0.25*distf[0]; 
  fnod[8] = -(0.25*distf[15])-0.25*distf[14]+0.25*distf[13]+0.25*distf[12]+0.25*distf[11]+0.25*distf[10]+0.25*distf[9]-0.25*distf[8]+0.25*distf[7]-0.25*distf[6]-0.25*distf[5]-0.25*distf[4]-0.25*distf[3]-0.25*distf[2]+0.25*distf[1]+0.25*distf[0]; 
  fnod[9] = 0.25*distf[15]+0.25*distf[14]-0.25*distf[13]-0.25*distf[12]+0.25*distf[11]-0.25*distf[10]-0.25*distf[9]+0.25*distf[8]+0.25*distf[7]-0.25*distf[6]-0.25*distf[5]+0.25*distf[4]-0.25*distf[3]-0.25*distf[2]+0.25*distf[1]+0.25*distf[0]; 
  fnod[10] = 0.25*distf[15]+0.25*distf[14]-0.25*distf[13]+0.25*distf[12]-0.25*distf[11]-0.25*distf[10]+0.25*distf[9]-0.25*distf[8]-0.25*distf[7]+0.25*distf[6]-0.25*distf[5]-0.25*distf[4]+0.25*distf[3]-0.25*distf[2]+0.25*distf[1]+0.25*distf[0]; 
  fnod[11] = -(0.25*distf[15])-0.25*distf[14]+0.25*distf[13]-0.25*distf[12]-0.25*distf[11]+0.25*distf[10]-0.25*distf[9]+0.25*distf[8]-0.25*distf[7]+0.25*distf[6]-0.25*distf[5]+0.25*distf[4]+0.25*distf[3]-0.25*distf[2]+0.25*distf[1]+0.25*distf[0]; 
  fnod[12] = 0.25*distf[15]+0.25*distf[14]+0.25*distf[13]-0.25*distf[12]-0.25*distf[11]+0.25*distf[10]-0.25*distf[9]-0.25*distf[8]-0.25*distf[7]-0.25*distf[6]+0.25*distf[5]-0.25*distf[4]-0.25*distf[3]+0.25*distf[2]+0.25*distf[1]+0.25*distf[0]; 
  fnod[13] = -(0.25*distf[15])-0.25*distf[14]-0.25*distf[13]+0.25*distf[12]-0.25*distf[11]-0.25*distf[10]+0.25*distf[9]+0.25*distf[8]-0.25*distf[7]-0.25*distf[6]+0.25*distf[5]+0.25*distf[4]-0.25*distf[3]+0.25*distf[2]+0.25*distf[1]+0.25*distf[0]; 
  fnod[14] = -(0.25*distf[15])-0.25*distf[14]-0.25*distf[13]-0.25*distf[12]+0.25*distf[11]-0.25*distf[10]-0.25*distf[9]-0.25*distf[8]+0.25*distf[7]+0.25*distf[6]+0.25*distf[5]-0.25*distf[4]+0.25*distf[3]+0.25*distf[2]+0.25*distf[1]+0.25*distf[0]; 
  fnod[15] = 0.25*distf[15]+0.25*distf[14]+0.25*distf[13]+0.25*distf[12]+0.25*distf[11]+0.25*distf[10]+0.25*distf[9]+0.25*distf[8]+0.25*distf[7]+0.25*distf[6]+0.25*distf[5]+0.25*distf[4]+0.25*distf[3]+0.25*distf[2]+0.25*distf[1]+0.25*distf[0]; 

  // If f < 0. at check nodes, set it to ffloor.
  if (fnod[0] < 0.) {
    fnod[0] = ffloor;
    shifted = true;
  }
  if (fnod[1] < 0.) {
    fnod[1] = ffloor;
    shifted = true;
  }
  if (fnod[2] < 0.) {
    fnod[2] = ffloor;
    shifted = true;
  }
  if (fnod[3] < 0.) {
    fnod[3] = ffloor;
    shifted = true;
  }
  if (fnod[4] < 0.) {
    fnod[4] = ffloor;
    shifted = true;
  }
  if (fnod[5] < 0.) {
    fnod[5] = ffloor;
    shifted = true;
  }
  if (fnod[6] < 0.) {
    fnod[6] = ffloor;
    shifted = true;
  }
  if (fnod[7] < 0.) {
    fnod[7] = ffloor;
    shifted = true;
  }
  if (fnod[8] < 0.) {
    fnod[8] = ffloor;
    shifted = true;
  }
  if (fnod[9] < 0.) {
    fnod[9] = ffloor;
    shifted = true;
  }
  if (fnod[10] < 0.) {
    fnod[10] = ffloor;
    shifted = true;
  }
  if (fnod[11] < 0.) {
    fnod[11] = ffloor;
    shifted = true;
  }
  if (fnod[12] < 0.) {
    fnod[12] = ffloor;
    shifted = true;
  }
  if (fnod[13] < 0.) {
    fnod[13] = ffloor;
    shifted = true;
  }
  if (fnod[14] < 0.) {
    fnod[14] = ffloor;
    shifted = true;
  }
  if (fnod[15] < 0.) {
    fnod[15] = ffloor;
    shifted = true;
  }

  if (shifted) {
  distf[0] = 0.25*fnod[15]+0.25*fnod[14]+0.25*fnod[13]+0.25*fnod[12]+0.25*fnod[11]+0.25*fnod[10]+0.25*fnod[9]+0.25*fnod[8]+0.25*fnod[7]+0.25*fnod[6]+0.25*fnod[5]+0.25*fnod[4]+0.25*fnod[3]+0.25*fnod[2]+0.25*fnod[1]+0.25*fnod[0]; 
  distf[1] = 0.25*fnod[15]+0.25*fnod[14]+0.25*fnod[13]+0.25*fnod[12]+0.25*fnod[11]+0.25*fnod[10]+0.25*fnod[9]+0.25*fnod[8]-0.25*fnod[7]-0.25*fnod[6]-0.25*fnod[5]-0.25*fnod[4]-0.25*fnod[3]-0.25*fnod[2]-0.25*fnod[1]-0.25*fnod[0]; 
  distf[2] = 0.25*fnod[15]+0.25*fnod[14]+0.25*fnod[13]+0.25*fnod[12]-0.25*fnod[11]-0.25*fnod[10]-0.25*fnod[9]-0.25*fnod[8]+0.25*fnod[7]+0.25*fnod[6]+0.25*fnod[5]+0.25*fnod[4]-0.25*fnod[3]-0.25*fnod[2]-0.25*fnod[1]-0.25*fnod[0]; 
  distf[3] = 0.25*fnod[15]+0.25*fnod[14]-0.25*fnod[13]-0.25*fnod[12]+0.25*fnod[11]+0.25*fnod[10]-0.25*fnod[9]-0.25*fnod[8]+0.25*fnod[7]+0.25*fnod[6]-0.25*fnod[5]-0.25*fnod[4]+0.25*fnod[3]+0.25*fnod[2]-0.25*fnod[1]-0.25*fnod[0]; 
  distf[4] = 0.25*fnod[15]-0.25*fnod[14]+0.25*fnod[13]-0.25*fnod[12]+0.25*fnod[11]-0.25*fnod[10]+0.25*fnod[9]-0.25*fnod[8]+0.25*fnod[7]-0.25*fnod[6]+0.25*fnod[5]-0.25*fnod[4]+0.25*fnod[3]-0.25*fnod[2]+0.25*fnod[1]-0.25*fnod[0]; 
  distf[5] = 0.25*fnod[15]+0.25*fnod[14]+0.25*fnod[13]+0.25*fnod[12]-0.25*fnod[11]-0.25*fnod[10]-0.25*fnod[9]-0.25*fnod[8]-0.25*fnod[7]-0.25*fnod[6]-0.25*fnod[5]-0.25*fnod[4]+0.25*fnod[3]+0.25*fnod[2]+0.25*fnod[1]+0.25*fnod[0]; 
  distf[6] = 0.25*fnod[15]+0.25*fnod[14]-0.25*fnod[13]-0.25*fnod[12]+0.25*fnod[11]+0.25*fnod[10]-0.25*fnod[9]-0.25*fnod[8]-0.25*fnod[7]-0.25*fnod[6]+0.25*fnod[5]+0.25*fnod[4]-0.25*fnod[3]-0.25*fnod[2]+0.25*fnod[1]+0.25*fnod[0]; 
  distf[7] = 0.25*fnod[15]+0.25*fnod[14]-0.25*fnod[13]-0.25*fnod[12]-0.25*fnod[11]-0.25*fnod[10]+0.25*fnod[9]+0.25*fnod[8]+0.25*fnod[7]+0.25*fnod[6]-0.25*fnod[5]-0.25*fnod[4]-0.25*fnod[3]-0.25*fnod[2]+0.25*fnod[1]+0.25*fnod[0]; 
  distf[8] = 0.25*fnod[15]-0.25*fnod[14]+0.25*fnod[13]-0.25*fnod[12]+0.25*fnod[11]-0.25*fnod[10]+0.25*fnod[9]-0.25*fnod[8]-0.25*fnod[7]+0.25*fnod[6]-0.25*fnod[5]+0.25*fnod[4]-0.25*fnod[3]+0.25*fnod[2]-0.25*fnod[1]+0.25*fnod[0]; 
  distf[9] = 0.25*fnod[15]-0.25*fnod[14]+0.25*fnod[13]-0.25*fnod[12]-0.25*fnod[11]+0.25*fnod[10]-0.25*fnod[9]+0.25*fnod[8]+0.25*fnod[7]-0.25*fnod[6]+0.25*fnod[5]-0.25*fnod[4]-0.25*fnod[3]+0.25*fnod[2]-0.25*fnod[1]+0.25*fnod[0]; 
  distf[10] = 0.25*fnod[15]-0.25*fnod[14]-0.25*fnod[13]+0.25*fnod[12]+0.25*fnod[11]-0.25*fnod[10]-0.25*fnod[9]+0.25*fnod[8]+0.25*fnod[7]-0.25*fnod[6]-0.25*fnod[5]+0.25*fnod[4]+0.25*fnod[3]-0.25*fnod[2]-0.25*fnod[1]+0.25*fnod[0]; 
  distf[11] = 0.25*fnod[15]+0.25*fnod[14]-0.25*fnod[13]-0.25*fnod[12]-0.25*fnod[11]-0.25*fnod[10]+0.25*fnod[9]+0.25*fnod[8]-0.25*fnod[7]-0.25*fnod[6]+0.25*fnod[5]+0.25*fnod[4]+0.25*fnod[3]+0.25*fnod[2]-0.25*fnod[1]-0.25*fnod[0]; 
  distf[12] = 0.25*fnod[15]-0.25*fnod[14]+0.25*fnod[13]-0.25*fnod[12]-0.25*fnod[11]+0.25*fnod[10]-0.25*fnod[9]+0.25*fnod[8]-0.25*fnod[7]+0.25*fnod[6]-0.25*fnod[5]+0.25*fnod[4]+0.25*fnod[3]-0.25*fnod[2]+0.25*fnod[1]-0.25*fnod[0]; 
  distf[13] = 0.25*fnod[15]-0.25*fnod[14]-0.25*fnod[13]+0.25*fnod[12]+0.25*fnod[11]-0.25*fnod[10]-0.25*fnod[9]+0.25*fnod[8]-0.25*fnod[7]+0.25*fnod[6]+0.25*fnod[5]-0.25*fnod[4]-0.25*fnod[3]+0.25*fnod[2]+0.25*fnod[1]-0.25*fnod[0]; 
  distf[14] = 0.25*fnod[15]-0.25*fnod[14]-0.25*fnod[13]+0.25*fnod[12]-0.25*fnod[11]+0.25*fnod[10]+0.25*fnod[9]-0.25*fnod[8]+0.25*fnod[7]-0.25*fnod[6]-0.25*fnod[5]+0.25*fnod[4]-0.25*fnod[3]+0.25*fnod[2]+0.25*fnod[1]-0.25*fnod[0]; 
  distf[15] = 0.25*fnod[15]-0.25*fnod[14]-0.25*fnod[13]+0.25*fnod[12]-0.25*fnod[11]+0.25*fnod[10]+0.25*fnod[9]-0.25*fnod[8]-0.25*fnod[7]+0.25*fnod[6]+0.25*fnod[5]-0.25*fnod[4]+0.25*fnod[3]-0.25*fnod[2]-0.25*fnod[1]+0.25*fnod[0]; 
  }

  return shifted;

}

GKYL_CU_DH bool positivity_shift_vlasov_MRS_limiter_1x3v_tensor_p1(double ffloor, double *distf) 
{ 
  // ffloor: Distribution function floor to shift to when f<0.
  // distf: distribution function.

  double fnod[16];
  fnod[0] = 0.25*distf[15]-0.25*distf[14]-0.25*distf[13]-0.25*distf[12]-0.25*distf[11]+0.25*distf[10]+0.25*distf[9]+0.25*distf[8]+0.25*distf[7]+0.25*distf[6]+0.25*distf[5]-0.25*distf[4]-0.25*distf[3]-0.25*distf[2]-0.25*distf[1]+0.25*distf[0]; 
  fnod[1] = -(0.25*distf[15])+0.25*distf[14]+0.25*distf[13]+0.25*distf[12]-0.25*distf[11]-0.25*distf[10]-0.25*distf[9]-0.25*distf[8]+0.25*distf[7]+0.25*distf[6]+0.25*distf[5]+0.25*distf[4]-0.25*distf[3]-0.25*distf[2]-0.25*distf[1]+0.25*distf[0]; 
  fnod[2] = -(0.25*distf[15])+0.25*distf[14]+0.25*distf[13]-0.25*distf[12]+0.25*distf[11]-0.25*distf[10]+0.25*distf[9]+0.25*distf[8]-0.25*distf[7]-0.25*distf[6]+0.25*distf[5]-0.25*distf[4]+0.25*distf[3]-0.25*distf[2]-0.25*distf[1]+0.25*distf[0]; 
  fnod[3] = 0.25*distf[15]-0.25*distf[14]-0.25*distf[13]+0.25*distf[12]+0.25*distf[11]+0.25*distf[10]-0.25*distf[9]-0.25*distf[8]-0.25*distf[7]-0.25*distf[6]+0.25*distf[5]+0.25*distf[4]+0.25*distf[3]-0.25*distf[2]-0.25*distf[1]+0.25*distf[0]; 
  fnod[4] = -(0.25*distf[15])+0.25*distf[14]-0.25*distf[13]+0.25*distf[12]+0.25*distf[11]+0.25*distf[10]-0.25*distf[9]+0.25*distf[8]-0.25*distf[7]+0.25*distf[6]-0.25*distf[5]-0.25*distf[4]-0.25*distf[3]+0.25*distf[2]-0.25*distf[1]+0.25*distf[0]; 
  fnod[5] = 0.25*distf[15]-0.25*distf[14]+0.25*distf[13]-0.25*distf[12]+0.25*distf[11]-0.25*distf[10]+0.25*distf[9]-0.25*distf[8]-0.25*distf[7]+0.25*distf[6]-0.25*distf[5]+0.25*distf[4]-0.25*distf[3]+0.25*distf[2]-0.25*distf[1]+0.25*distf[0]; 
  fnod[6] = 0.25*distf[15]-0.25*distf[14]+0.25*distf[13]+0.25*distf[12]-0.25*distf[11]-0.25*distf[10]-0.25*distf[9]+0.25*distf[8]+0.25*distf[7]-0.25*distf[6]-0.25*distf[5]-0.25*distf[4]+0.25*distf[3]+0.25*distf[2]-0.25*distf[1]+0.25*distf[0]; 
  fnod[7] = -(0.25*distf[15])+0.25*distf[14]-0.25*distf[13]-0.25*distf[12]-0.25*distf[11]+0.25*distf[10]+0.25*distf[9]-0.25*distf[8]+0.25*distf[7]-0.25*distf[6]-0.25*distf[5]+0.25*distf[4]+0.25*distf[3]+0.25*distf[2]-0.25*distf[1]+0.25*distf[0]; 
  fnod[8] = -(0.25*distf[15])-0.25*distf[14]+0.25*distf[13]+0.25*distf[12]+0.25*distf[11]+0.25*distf[10]+0.25*distf[9]-0.25*distf[8]+0.25*distf[7]-0.25*distf[6]-0.25*distf[5]-0.25*distf[4]-0.25*distf[3]-0.25*distf[2]+0.25*distf[1]+0.25*distf[0]; 
  fnod[9] = 0.25*distf[15]+0.25*distf[14]-0.25*distf[13]-0.25*distf[12]+0.25*distf[11]-0.25*distf[10]-0.25*distf[9]+0.25*distf[8]+0.25*distf[7]-0.25*distf[6]-0.25*distf[5]+0.25*distf[4]-0.25*distf[3]-0.25*distf[2]+0.25*distf[1]+0.25*distf[0]; 
  fnod[10] = 0.25*distf[15]+0.25*distf[14]-0.25*distf[13]+0.25*distf[12]-0.25*distf[11]-0.25*distf[10]+0.25*distf[9]-0.25*distf[8]-0.25*distf[7]+0.25*distf[6]-0.25*distf[5]-0.25*distf[4]+0.25*distf[3]-0.25*distf[2]+0.25*distf[1]+0.25*distf[0]; 
  fnod[11] = -(0.25*distf[15])-0.25*distf[14]+0.25*distf[13]-0.25*distf[12]-0.25*distf[11]+0.25*distf[10]-0.25*distf[9]+0.25*distf[8]-0.25*distf[7]+0.25*distf[6]-0.25*distf[5]+0.25*distf[4]+0.25*distf[3]-0.25*distf[2]+0.25*distf[1]+0.25*distf[0]; 
  fnod[12] = 0.25*distf[15]+0.25*distf[14]+0.25*distf[13]-0.25*distf[12]-0.25*distf[11]+0.25*distf[10]-0.25*distf[9]-0.25*distf[8]-0.25*distf[7]-0.25*distf[6]+0.25*distf[5]-0.25*distf[4]-0.25*distf[3]+0.25*distf[2]+0.25*distf[1]+0.25*distf[0]; 
  fnod[13] = -(0.25*distf[15])-0.25*distf[14]-0.25*distf[13]+0.25*distf[12]-0.25*distf[11]-0.25*distf[10]+0.25*distf[9]+0.25*distf[8]-0.25*distf[7]-0.25*distf[6]+0.25*distf[5]+0.25*distf[4]-0.25*distf[3]+0.25*distf[2]+0.25*distf[1]+0.25*distf[0]; 
  fnod[14] = -(0.25*distf[15])-0.25*distf[14]-0.25*distf[13]-0.25*distf[12]+0.25*distf[11]-0.25*distf[10]-0.25*distf[9]-0.25*distf[8]+0.25*distf[7]+0.25*distf[6]+0.25*distf[5]-0.25*distf[4]+0.25*distf[3]+0.25*distf[2]+0.25*distf[1]+0.25*distf[0]; 
  fnod[15] = 0.25*distf[15]+0.25*distf[14]+0.25*distf[13]+0.25*distf[12]+0.25*distf[11]+0.25*distf[10]+0.25*distf[9]+0.25*distf[8]+0.25*distf[7]+0.25*distf[6]+0.25*distf[5]+0.25*distf[4]+0.25*distf[3]+0.25*distf[2]+0.25*distf[1]+0.25*distf[0]; 

  bool shifted_node = false;

  if (distf[0] > 0.) {
    // Apply Moe–Rossmanith–Seal limiter.
    double fnod_min = DBL_MAX;
    fnod_min = fmin(fnod_min, fnod[0]);
    fnod_min = fmin(fnod_min, fnod[1]);
    fnod_min = fmin(fnod_min, fnod[2]);
    fnod_min = fmin(fnod_min, fnod[3]);
    fnod_min = fmin(fnod_min, fnod[4]);
    fnod_min = fmin(fnod_min, fnod[5]);
    fnod_min = fmin(fnod_min, fnod[6]);
    fnod_min = fmin(fnod_min, fnod[7]);
    fnod_min = fmin(fnod_min, fnod[8]);
    fnod_min = fmin(fnod_min, fnod[9]);
    fnod_min = fmin(fnod_min, fnod[10]);
    fnod_min = fmin(fnod_min, fnod[11]);
    fnod_min = fmin(fnod_min, fnod[12]);
    fnod_min = fmin(fnod_min, fnod[13]);
    fnod_min = fmin(fnod_min, fnod[14]);
    fnod_min = fmin(fnod_min, fnod[15]);

    if (fnod_min < 0.0) {
      double f_cellav = distf[0]/4.000000000000001;
      double denom = f_cellav - fnod_min;
      double theta = denom > 1.0e-12*f_cellav? fmin(1.0, f_cellav/denom) : 1.0;

  distf[0] = -(4.0*f_cellav*theta)+distf[0]*theta+4.0*f_cellav; 
  distf[1] = distf[1]*theta; 
  distf[2] = distf[2]*theta; 
  distf[3] = distf[3]*theta; 
  distf[4] = distf[4]*theta; 
  distf[5] = distf[5]*theta; 
  distf[6] = distf[6]*theta; 
  distf[7] = distf[7]*theta; 
  distf[8] = distf[8]*theta; 
  distf[9] = distf[9]*theta; 
  distf[10] = distf[10]*theta; 
  distf[11] = distf[11]*theta; 
  distf[12] = distf[12]*theta; 
  distf[13] = distf[13]*theta; 
  distf[14] = distf[14]*theta; 
  distf[15] = distf[15]*theta; 

      shifted_node = true;
    }
  }

  else {

    // If f < 0. at check nodes, set it to ffloor.
    if (fnod[0] < 0.) fnod[0] = ffloor;
    if (fnod[1] < 0.) fnod[1] = ffloor;
    if (fnod[2] < 0.) fnod[2] = ffloor;
    if (fnod[3] < 0.) fnod[3] = ffloor;
    if (fnod[4] < 0.) fnod[4] = ffloor;
    if (fnod[5] < 0.) fnod[5] = ffloor;
    if (fnod[6] < 0.) fnod[6] = ffloor;
    if (fnod[7] < 0.) fnod[7] = ffloor;
    if (fnod[8] < 0.) fnod[8] = ffloor;
    if (fnod[9] < 0.) fnod[9] = ffloor;
    if (fnod[10] < 0.) fnod[10] = ffloor;
    if (fnod[11] < 0.) fnod[11] = ffloor;
    if (fnod[12] < 0.) fnod[12] = ffloor;
    if (fnod[13] < 0.) fnod[13] = ffloor;
    if (fnod[14] < 0.) fnod[14] = ffloor;
    if (fnod[15] < 0.) fnod[15] = ffloor;

  distf[0] = 0.25*fnod[15]+0.25*fnod[14]+0.25*fnod[13]+0.25*fnod[12]+0.25*fnod[11]+0.25*fnod[10]+0.25*fnod[9]+0.25*fnod[8]+0.25*fnod[7]+0.25*fnod[6]+0.25*fnod[5]+0.25*fnod[4]+0.25*fnod[3]+0.25*fnod[2]+0.25*fnod[1]+0.25*fnod[0]; 
  distf[1] = 0.25*fnod[15]+0.25*fnod[14]+0.25*fnod[13]+0.25*fnod[12]+0.25*fnod[11]+0.25*fnod[10]+0.25*fnod[9]+0.25*fnod[8]-0.25*fnod[7]-0.25*fnod[6]-0.25*fnod[5]-0.25*fnod[4]-0.25*fnod[3]-0.25*fnod[2]-0.25*fnod[1]-0.25*fnod[0]; 
  distf[2] = 0.25*fnod[15]+0.25*fnod[14]+0.25*fnod[13]+0.25*fnod[12]-0.25*fnod[11]-0.25*fnod[10]-0.25*fnod[9]-0.25*fnod[8]+0.25*fnod[7]+0.25*fnod[6]+0.25*fnod[5]+0.25*fnod[4]-0.25*fnod[3]-0.25*fnod[2]-0.25*fnod[1]-0.25*fnod[0]; 
  distf[3] = 0.25*fnod[15]+0.25*fnod[14]-0.25*fnod[13]-0.25*fnod[12]+0.25*fnod[11]+0.25*fnod[10]-0.25*fnod[9]-0.25*fnod[8]+0.25*fnod[7]+0.25*fnod[6]-0.25*fnod[5]-0.25*fnod[4]+0.25*fnod[3]+0.25*fnod[2]-0.25*fnod[1]-0.25*fnod[0]; 
  distf[4] = 0.25*fnod[15]-0.25*fnod[14]+0.25*fnod[13]-0.25*fnod[12]+0.25*fnod[11]-0.25*fnod[10]+0.25*fnod[9]-0.25*fnod[8]+0.25*fnod[7]-0.25*fnod[6]+0.25*fnod[5]-0.25*fnod[4]+0.25*fnod[3]-0.25*fnod[2]+0.25*fnod[1]-0.25*fnod[0]; 
  distf[5] = 0.25*fnod[15]+0.25*fnod[14]+0.25*fnod[13]+0.25*fnod[12]-0.25*fnod[11]-0.25*fnod[10]-0.25*fnod[9]-0.25*fnod[8]-0.25*fnod[7]-0.25*fnod[6]-0.25*fnod[5]-0.25*fnod[4]+0.25*fnod[3]+0.25*fnod[2]+0.25*fnod[1]+0.25*fnod[0]; 
  distf[6] = 0.25*fnod[15]+0.25*fnod[14]-0.25*fnod[13]-0.25*fnod[12]+0.25*fnod[11]+0.25*fnod[10]-0.25*fnod[9]-0.25*fnod[8]-0.25*fnod[7]-0.25*fnod[6]+0.25*fnod[5]+0.25*fnod[4]-0.25*fnod[3]-0.25*fnod[2]+0.25*fnod[1]+0.25*fnod[0]; 
  distf[7] = 0.25*fnod[15]+0.25*fnod[14]-0.25*fnod[13]-0.25*fnod[12]-0.25*fnod[11]-0.25*fnod[10]+0.25*fnod[9]+0.25*fnod[8]+0.25*fnod[7]+0.25*fnod[6]-0.25*fnod[5]-0.25*fnod[4]-0.25*fnod[3]-0.25*fnod[2]+0.25*fnod[1]+0.25*fnod[0]; 
  distf[8] = 0.25*fnod[15]-0.25*fnod[14]+0.25*fnod[13]-0.25*fnod[12]+0.25*fnod[11]-0.25*fnod[10]+0.25*fnod[9]-0.25*fnod[8]-0.25*fnod[7]+0.25*fnod[6]-0.25*fnod[5]+0.25*fnod[4]-0.25*fnod[3]+0.25*fnod[2]-0.25*fnod[1]+0.25*fnod[0]; 
  distf[9] = 0.25*fnod[15]-0.25*fnod[14]+0.25*fnod[13]-0.25*fnod[12]-0.25*fnod[11]+0.25*fnod[10]-0.25*fnod[9]+0.25*fnod[8]+0.25*fnod[7]-0.25*fnod[6]+0.25*fnod[5]-0.25*fnod[4]-0.25*fnod[3]+0.25*fnod[2]-0.25*fnod[1]+0.25*fnod[0]; 
  distf[10] = 0.25*fnod[15]-0.25*fnod[14]-0.25*fnod[13]+0.25*fnod[12]+0.25*fnod[11]-0.25*fnod[10]-0.25*fnod[9]+0.25*fnod[8]+0.25*fnod[7]-0.25*fnod[6]-0.25*fnod[5]+0.25*fnod[4]+0.25*fnod[3]-0.25*fnod[2]-0.25*fnod[1]+0.25*fnod[0]; 
  distf[11] = 0.25*fnod[15]+0.25*fnod[14]-0.25*fnod[13]-0.25*fnod[12]-0.25*fnod[11]-0.25*fnod[10]+0.25*fnod[9]+0.25*fnod[8]-0.25*fnod[7]-0.25*fnod[6]+0.25*fnod[5]+0.25*fnod[4]+0.25*fnod[3]+0.25*fnod[2]-0.25*fnod[1]-0.25*fnod[0]; 
  distf[12] = 0.25*fnod[15]-0.25*fnod[14]+0.25*fnod[13]-0.25*fnod[12]-0.25*fnod[11]+0.25*fnod[10]-0.25*fnod[9]+0.25*fnod[8]-0.25*fnod[7]+0.25*fnod[6]-0.25*fnod[5]+0.25*fnod[4]+0.25*fnod[3]-0.25*fnod[2]+0.25*fnod[1]-0.25*fnod[0]; 
  distf[13] = 0.25*fnod[15]-0.25*fnod[14]-0.25*fnod[13]+0.25*fnod[12]+0.25*fnod[11]-0.25*fnod[10]-0.25*fnod[9]+0.25*fnod[8]-0.25*fnod[7]+0.25*fnod[6]+0.25*fnod[5]-0.25*fnod[4]-0.25*fnod[3]+0.25*fnod[2]+0.25*fnod[1]-0.25*fnod[0]; 
  distf[14] = 0.25*fnod[15]-0.25*fnod[14]-0.25*fnod[13]+0.25*fnod[12]-0.25*fnod[11]+0.25*fnod[10]+0.25*fnod[9]-0.25*fnod[8]+0.25*fnod[7]-0.25*fnod[6]-0.25*fnod[5]+0.25*fnod[4]-0.25*fnod[3]+0.25*fnod[2]+0.25*fnod[1]-0.25*fnod[0]; 
  distf[15] = 0.25*fnod[15]-0.25*fnod[14]-0.25*fnod[13]+0.25*fnod[12]-0.25*fnod[11]+0.25*fnod[10]+0.25*fnod[9]-0.25*fnod[8]-0.25*fnod[7]+0.25*fnod[6]+0.25*fnod[5]-0.25*fnod[4]+0.25*fnod[3]-0.25*fnod[2]-0.25*fnod[1]+0.25*fnod[0]; 

    shifted_node = true;
  }

  return shifted_node;

}
