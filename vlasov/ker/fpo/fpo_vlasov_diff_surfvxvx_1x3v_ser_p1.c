#include <gkyl_fpo_vlasov_kernels.h>

GKYL_CU_DH double
fpo_vlasov_diff_surfvxvx_1x3v_ser_p1(const double* w, const double* dx,
  const double* g[], const double* f[], double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // g: 
  // f: 
  // out: Incremented output

  const double Jvxvx = 4/dx[1]/dx[1];

  const double* gllc = g[1];
  const double* glcl = g[3];
  const double* glcc = g[4];
  const double* glcu = g[5];
  const double* gluc = g[7];
  const double* gcll = g[9];
  const double* gclc = g[10];
  const double* gclu = g[11];
  const double* gccl = g[12];
  const double* gccc = g[13];
  const double* gccu = g[14];
  const double* gcul = g[15];
  const double* gcuc = g[16];
  const double* gcuu = g[17];
  const double* gulc = g[19];
  const double* gucl = g[21];
  const double* gucc = g[22];
  const double* gucu = g[23];
  const double* guuc = g[25];

  const double* fllc = f[1];
  const double* flcl = f[3];
  const double* flcc = f[4];
  const double* flcu = f[5];
  const double* fluc = f[7];
  const double* fcll = f[9];
  const double* fclc = f[10];
  const double* fclu = f[11];
  const double* fccl = f[12];
  const double* fccc = f[13];
  const double* fccu = f[14];
  const double* fcul = f[15];
  const double* fcuc = f[16];
  const double* fcuu = f[17];
  const double* fulc = f[19];
  const double* fucl = f[21];
  const double* fucc = f[22];
  const double* fucu = f[23];
  const double* fuuc = f[25];

  double D_proj1_l[8];
  D_proj1_l[0] = 0.6123724356957944*gccc[2]-0.6123724356957944*glcc[2];
  D_proj1_l[1] = 0.6123724356957944*gccc[5]-0.6123724356957944*glcc[5];
  D_proj1_l[2] = 0.6123724356957944*gccc[7]-0.6123724356957944*glcc[7];
  D_proj1_l[3] = 0.6123724356957944*gccc[9]-0.6123724356957944*glcc[9];
  D_proj1_l[4] = 0.6123724356957944*gccc[11]-0.6123724356957944*glcc[11];
  D_proj1_l[5] = 0.6123724356957944*gccc[12]-0.6123724356957944*glcc[12];
  D_proj1_l[6] = 0.6123724356957944*gccc[14]-0.6123724356957944*glcc[14];
  D_proj1_l[7] = 0.6123724356957944*gccc[15]-0.6123724356957944*glcc[15];

  double D_proj1_u[8];
  D_proj1_u[0] = 0.6123724356957944*gucc[2]-0.6123724356957944*gccc[2];
  D_proj1_u[1] = 0.6123724356957944*gucc[5]-0.6123724356957944*gccc[5];
  D_proj1_u[2] = 0.6123724356957944*gucc[7]-0.6123724356957944*gccc[7];
  D_proj1_u[3] = 0.6123724356957944*gucc[9]-0.6123724356957944*gccc[9];
  D_proj1_u[4] = 0.6123724356957944*gucc[11]-0.6123724356957944*gccc[11];
  D_proj1_u[5] = 0.6123724356957944*gucc[12]-0.6123724356957944*gccc[12];
  D_proj1_u[6] = 0.6123724356957944*gucc[14]-0.6123724356957944*gccc[14];
  D_proj1_u[7] = 0.6123724356957944*gucc[15]-0.6123724356957944*gccc[15];

  double df_proj1_l[8];
  df_proj1_l[0] = (-0.7654655446197428*flcc[2])-0.7654655446197428*fccc[2]-0.7954951288348656*flcc[0]+0.7954951288348656*fccc[0];
  df_proj1_l[1] = (-0.7654655446197428*flcc[5])-0.7654655446197428*fccc[5]-0.7954951288348656*flcc[1]+0.7954951288348656*fccc[1];
  df_proj1_l[2] = (-0.7654655446197428*flcc[7])-0.7654655446197428*fccc[7]-0.7954951288348656*flcc[3]+0.7954951288348656*fccc[3];
  df_proj1_l[3] = (-0.7654655446197428*flcc[9])-0.7654655446197428*fccc[9]-0.7954951288348656*flcc[4]+0.7954951288348656*fccc[4];
  df_proj1_l[4] = (-0.7654655446197428*flcc[11])-0.7654655446197428*fccc[11]-0.7954951288348656*flcc[6]+0.7954951288348656*fccc[6];
  df_proj1_l[5] = (-0.7654655446197428*flcc[12])-0.7654655446197428*fccc[12]-0.7954951288348656*flcc[8]+0.7954951288348656*fccc[8];
  df_proj1_l[6] = (-0.7654655446197428*flcc[14])-0.7654655446197428*fccc[14]-0.7954951288348656*flcc[10]+0.7954951288348656*fccc[10];
  df_proj1_l[7] = (-0.7654655446197428*flcc[15])-0.7654655446197428*fccc[15]-0.7954951288348656*flcc[13]+0.7954951288348656*fccc[13];

  double df_proj1_u[8];
  df_proj1_u[0] = (-0.7654655446197428*fucc[2])-0.7654655446197428*fccc[2]+0.7954951288348656*fucc[0]-0.7954951288348656*fccc[0];
  df_proj1_u[1] = (-0.7654655446197428*fucc[5])-0.7654655446197428*fccc[5]+0.7954951288348656*fucc[1]-0.7954951288348656*fccc[1];
  df_proj1_u[2] = (-0.7654655446197428*fucc[7])-0.7654655446197428*fccc[7]+0.7954951288348656*fucc[3]-0.7954951288348656*fccc[3];
  df_proj1_u[3] = (-0.7654655446197428*fucc[9])-0.7654655446197428*fccc[9]+0.7954951288348656*fucc[4]-0.7954951288348656*fccc[4];
  df_proj1_u[4] = (-0.7654655446197428*fucc[11])-0.7654655446197428*fccc[11]+0.7954951288348656*fucc[6]-0.7954951288348656*fccc[6];
  df_proj1_u[5] = (-0.7654655446197428*fucc[12])-0.7654655446197428*fccc[12]+0.7954951288348656*fucc[8]-0.7954951288348656*fccc[8];
  df_proj1_u[6] = (-0.7654655446197428*fucc[14])-0.7654655446197428*fccc[14]+0.7954951288348656*fucc[10]-0.7954951288348656*fccc[10];
  df_proj1_u[7] = (-0.7654655446197428*fucc[15])-0.7654655446197428*fccc[15]+0.7954951288348656*fucc[13]-0.7954951288348656*fccc[13];

  double f_proj2_l[8];
  f_proj2_l[0] = 0.408248290463863*flcc[2]-0.408248290463863*fccc[2]+0.3535533905932737*flcc[0]+0.3535533905932737*fccc[0];
  f_proj2_l[1] = 0.408248290463863*flcc[5]-0.408248290463863*fccc[5]+0.3535533905932737*flcc[1]+0.3535533905932737*fccc[1];
  f_proj2_l[2] = 0.408248290463863*flcc[7]-0.408248290463863*fccc[7]+0.3535533905932737*flcc[3]+0.3535533905932737*fccc[3];
  f_proj2_l[3] = 0.408248290463863*flcc[9]-0.408248290463863*fccc[9]+0.3535533905932737*flcc[4]+0.3535533905932737*fccc[4];
  f_proj2_l[4] = 0.408248290463863*flcc[11]-0.408248290463863*fccc[11]+0.3535533905932737*flcc[6]+0.3535533905932737*fccc[6];
  f_proj2_l[5] = 0.408248290463863*flcc[12]-0.408248290463863*fccc[12]+0.3535533905932737*flcc[8]+0.3535533905932737*fccc[8];
  f_proj2_l[6] = 0.408248290463863*flcc[14]-0.408248290463863*fccc[14]+0.3535533905932737*flcc[10]+0.3535533905932737*fccc[10];
  f_proj2_l[7] = 0.408248290463863*flcc[15]-0.408248290463863*fccc[15]+0.3535533905932737*flcc[13]+0.3535533905932737*fccc[13];

  double f_proj2_u[8];
  f_proj2_u[0] = (-0.408248290463863*fucc[2])+0.408248290463863*fccc[2]+0.3535533905932737*fucc[0]+0.3535533905932737*fccc[0];
  f_proj2_u[1] = (-0.408248290463863*fucc[5])+0.408248290463863*fccc[5]+0.3535533905932737*fucc[1]+0.3535533905932737*fccc[1];
  f_proj2_u[2] = (-0.408248290463863*fucc[7])+0.408248290463863*fccc[7]+0.3535533905932737*fucc[3]+0.3535533905932737*fccc[3];
  f_proj2_u[3] = (-0.408248290463863*fucc[9])+0.408248290463863*fccc[9]+0.3535533905932737*fucc[4]+0.3535533905932737*fccc[4];
  f_proj2_u[4] = (-0.408248290463863*fucc[11])+0.408248290463863*fccc[11]+0.3535533905932737*fucc[6]+0.3535533905932737*fccc[6];
  f_proj2_u[5] = (-0.408248290463863*fucc[12])+0.408248290463863*fccc[12]+0.3535533905932737*fucc[8]+0.3535533905932737*fccc[8];
  f_proj2_u[6] = (-0.408248290463863*fucc[14])+0.408248290463863*fccc[14]+0.3535533905932737*fucc[10]+0.3535533905932737*fccc[10];
  f_proj2_u[7] = (-0.408248290463863*fucc[15])+0.408248290463863*fccc[15]+0.3535533905932737*fucc[13]+0.3535533905932737*fccc[13];

  out[0] +=  Jvxvx*(0.125*D_proj1_u[7]*df_proj1_u[7]-0.125*D_proj1_l[7]*df_proj1_l[7]+0.125*D_proj1_u[6]*df_proj1_u[6]-0.125*D_proj1_l[6]*df_proj1_l[6]+0.125*D_proj1_u[5]*df_proj1_u[5]-0.125*D_proj1_l[5]*df_proj1_l[5]+0.125*D_proj1_u[4]*df_proj1_u[4]-0.125*D_proj1_l[4]*df_proj1_l[4]+0.125*D_proj1_u[3]*df_proj1_u[3]-0.125*D_proj1_l[3]*df_proj1_l[3]+0.125*D_proj1_u[2]*df_proj1_u[2]-0.125*D_proj1_l[2]*df_proj1_l[2]+0.125*D_proj1_u[1]*df_proj1_u[1]-0.125*D_proj1_l[1]*df_proj1_l[1]+0.125*D_proj1_u[0]*df_proj1_u[0]-0.125*D_proj1_l[0]*df_proj1_l[0]);
  out[1] +=  Jvxvx*(0.125*D_proj1_u[6]*df_proj1_u[7]-0.125*D_proj1_l[6]*df_proj1_l[7]+0.125*df_proj1_u[6]*D_proj1_u[7]-0.125*df_proj1_l[6]*D_proj1_l[7]+0.125*D_proj1_u[3]*df_proj1_u[5]-0.125*D_proj1_l[3]*df_proj1_l[5]+0.125*df_proj1_u[3]*D_proj1_u[5]-0.125*df_proj1_l[3]*D_proj1_l[5]+0.125*D_proj1_u[2]*df_proj1_u[4]-0.125*D_proj1_l[2]*df_proj1_l[4]+0.125*df_proj1_u[2]*D_proj1_u[4]-0.125*df_proj1_l[2]*D_proj1_l[4]+0.125*D_proj1_u[0]*df_proj1_u[1]-0.125*D_proj1_l[0]*df_proj1_l[1]+0.125*df_proj1_u[0]*D_proj1_u[1]-0.125*df_proj1_l[0]*D_proj1_l[1]);
  out[2] +=  Jvxvx*((-0.2165063509461096*D_proj1_u[7]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[7]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[7]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[7]*df_proj1_l[7]-0.2165063509461096*D_proj1_u[6]*f_proj2_u[6]+0.2165063509461096*D_proj1_l[6]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[6]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[6]*df_proj1_l[6]-0.2165063509461096*D_proj1_u[5]*f_proj2_u[5]+0.2165063509461096*D_proj1_l[5]*f_proj2_l[5]+0.2165063509461096*D_proj1_u[5]*df_proj1_u[5]+0.2165063509461096*D_proj1_l[5]*df_proj1_l[5]-0.2165063509461096*D_proj1_u[4]*f_proj2_u[4]+0.2165063509461096*D_proj1_l[4]*f_proj2_l[4]+0.2165063509461096*D_proj1_u[4]*df_proj1_u[4]+0.2165063509461096*D_proj1_l[4]*df_proj1_l[4]-0.2165063509461096*D_proj1_u[3]*f_proj2_u[3]+0.2165063509461096*D_proj1_l[3]*f_proj2_l[3]+0.2165063509461096*D_proj1_u[3]*df_proj1_u[3]+0.2165063509461096*D_proj1_l[3]*df_proj1_l[3]-0.2165063509461096*D_proj1_u[2]*f_proj2_u[2]+0.2165063509461096*D_proj1_l[2]*f_proj2_l[2]+0.2165063509461096*D_proj1_u[2]*df_proj1_u[2]+0.2165063509461096*D_proj1_l[2]*df_proj1_l[2]-0.2165063509461096*D_proj1_u[1]*f_proj2_u[1]+0.2165063509461096*D_proj1_l[1]*f_proj2_l[1]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[1]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[1]-0.2165063509461096*D_proj1_u[0]*f_proj2_u[0]+0.2165063509461096*D_proj1_l[0]*f_proj2_l[0]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[0]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[0]);
  out[3] +=  Jvxvx*(0.125*D_proj1_u[5]*df_proj1_u[7]-0.125*D_proj1_l[5]*df_proj1_l[7]+0.125*df_proj1_u[5]*D_proj1_u[7]-0.125*df_proj1_l[5]*D_proj1_l[7]+0.125*D_proj1_u[3]*df_proj1_u[6]-0.125*D_proj1_l[3]*df_proj1_l[6]+0.125*df_proj1_u[3]*D_proj1_u[6]-0.125*df_proj1_l[3]*D_proj1_l[6]+0.125*D_proj1_u[1]*df_proj1_u[4]-0.125*D_proj1_l[1]*df_proj1_l[4]+0.125*df_proj1_u[1]*D_proj1_u[4]-0.125*df_proj1_l[1]*D_proj1_l[4]+0.125*D_proj1_u[0]*df_proj1_u[2]-0.125*D_proj1_l[0]*df_proj1_l[2]+0.125*df_proj1_u[0]*D_proj1_u[2]-0.125*df_proj1_l[0]*D_proj1_l[2]);
  out[4] +=  Jvxvx*(0.125*D_proj1_u[4]*df_proj1_u[7]-0.125*D_proj1_l[4]*df_proj1_l[7]+0.125*df_proj1_u[4]*D_proj1_u[7]-0.125*df_proj1_l[4]*D_proj1_l[7]+0.125*D_proj1_u[2]*df_proj1_u[6]-0.125*D_proj1_l[2]*df_proj1_l[6]+0.125*df_proj1_u[2]*D_proj1_u[6]-0.125*df_proj1_l[2]*D_proj1_l[6]+0.125*D_proj1_u[1]*df_proj1_u[5]-0.125*D_proj1_l[1]*df_proj1_l[5]+0.125*df_proj1_u[1]*D_proj1_u[5]-0.125*df_proj1_l[1]*D_proj1_l[5]+0.125*D_proj1_u[0]*df_proj1_u[3]-0.125*D_proj1_l[0]*df_proj1_l[3]+0.125*df_proj1_u[0]*D_proj1_u[3]-0.125*df_proj1_l[0]*D_proj1_l[3]);
  out[5] +=  Jvxvx*((-0.2165063509461096*D_proj1_u[6]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[6]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[6]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[6]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[6]*D_proj1_u[7]+0.2165063509461096*df_proj1_u[6]*D_proj1_u[7]+0.2165063509461096*f_proj2_l[6]*D_proj1_l[7]+0.2165063509461096*df_proj1_l[6]*D_proj1_l[7]-0.2165063509461096*D_proj1_u[3]*f_proj2_u[5]+0.2165063509461096*D_proj1_l[3]*f_proj2_l[5]+0.2165063509461096*D_proj1_u[3]*df_proj1_u[5]+0.2165063509461096*D_proj1_l[3]*df_proj1_l[5]-0.2165063509461096*f_proj2_u[3]*D_proj1_u[5]+0.2165063509461096*df_proj1_u[3]*D_proj1_u[5]+0.2165063509461096*f_proj2_l[3]*D_proj1_l[5]+0.2165063509461096*df_proj1_l[3]*D_proj1_l[5]-0.2165063509461096*D_proj1_u[2]*f_proj2_u[4]+0.2165063509461096*D_proj1_l[2]*f_proj2_l[4]+0.2165063509461096*D_proj1_u[2]*df_proj1_u[4]+0.2165063509461096*D_proj1_l[2]*df_proj1_l[4]-0.2165063509461096*f_proj2_u[2]*D_proj1_u[4]+0.2165063509461096*df_proj1_u[2]*D_proj1_u[4]+0.2165063509461096*f_proj2_l[2]*D_proj1_l[4]+0.2165063509461096*df_proj1_l[2]*D_proj1_l[4]-0.2165063509461096*D_proj1_u[0]*f_proj2_u[1]+0.2165063509461096*D_proj1_l[0]*f_proj2_l[1]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[1]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[1]-0.2165063509461096*f_proj2_u[0]*D_proj1_u[1]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[1]+0.2165063509461096*f_proj2_l[0]*D_proj1_l[1]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[1]);
  out[6] +=  Jvxvx*(0.125*D_proj1_u[3]*df_proj1_u[7]-0.125*D_proj1_l[3]*df_proj1_l[7]+0.125*df_proj1_u[3]*D_proj1_u[7]-0.125*df_proj1_l[3]*D_proj1_l[7]+0.125*D_proj1_u[5]*df_proj1_u[6]-0.125*D_proj1_l[5]*df_proj1_l[6]+0.125*df_proj1_u[5]*D_proj1_u[6]-0.125*df_proj1_l[5]*D_proj1_l[6]+0.125*D_proj1_u[0]*df_proj1_u[4]-0.125*D_proj1_l[0]*df_proj1_l[4]+0.125*df_proj1_u[0]*D_proj1_u[4]-0.125*df_proj1_l[0]*D_proj1_l[4]+0.125*D_proj1_u[1]*df_proj1_u[2]-0.125*D_proj1_l[1]*df_proj1_l[2]+0.125*df_proj1_u[1]*D_proj1_u[2]-0.125*df_proj1_l[1]*D_proj1_l[2]);
  out[7] +=  Jvxvx*((-0.2165063509461096*D_proj1_u[5]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[5]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[5]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[5]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[5]*D_proj1_u[7]+0.2165063509461096*df_proj1_u[5]*D_proj1_u[7]+0.2165063509461096*f_proj2_l[5]*D_proj1_l[7]+0.2165063509461096*df_proj1_l[5]*D_proj1_l[7]-0.2165063509461096*D_proj1_u[3]*f_proj2_u[6]+0.2165063509461096*D_proj1_l[3]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[3]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[3]*df_proj1_l[6]-0.2165063509461096*f_proj2_u[3]*D_proj1_u[6]+0.2165063509461096*df_proj1_u[3]*D_proj1_u[6]+0.2165063509461096*f_proj2_l[3]*D_proj1_l[6]+0.2165063509461096*df_proj1_l[3]*D_proj1_l[6]-0.2165063509461096*D_proj1_u[1]*f_proj2_u[4]+0.2165063509461096*D_proj1_l[1]*f_proj2_l[4]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[4]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[4]-0.2165063509461096*f_proj2_u[1]*D_proj1_u[4]+0.2165063509461096*df_proj1_u[1]*D_proj1_u[4]+0.2165063509461096*f_proj2_l[1]*D_proj1_l[4]+0.2165063509461096*df_proj1_l[1]*D_proj1_l[4]-0.2165063509461096*D_proj1_u[0]*f_proj2_u[2]+0.2165063509461096*D_proj1_l[0]*f_proj2_l[2]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[2]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[2]-0.2165063509461096*f_proj2_u[0]*D_proj1_u[2]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[2]+0.2165063509461096*f_proj2_l[0]*D_proj1_l[2]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[2]);
  out[8] +=  Jvxvx*(0.125*D_proj1_u[2]*df_proj1_u[7]-0.125*D_proj1_l[2]*df_proj1_l[7]+0.125*df_proj1_u[2]*D_proj1_u[7]-0.125*df_proj1_l[2]*D_proj1_l[7]+0.125*D_proj1_u[4]*df_proj1_u[6]-0.125*D_proj1_l[4]*df_proj1_l[6]+0.125*df_proj1_u[4]*D_proj1_u[6]-0.125*df_proj1_l[4]*D_proj1_l[6]+0.125*D_proj1_u[0]*df_proj1_u[5]-0.125*D_proj1_l[0]*df_proj1_l[5]+0.125*df_proj1_u[0]*D_proj1_u[5]-0.125*df_proj1_l[0]*D_proj1_l[5]+0.125*D_proj1_u[1]*df_proj1_u[3]-0.125*D_proj1_l[1]*df_proj1_l[3]+0.125*df_proj1_u[1]*D_proj1_u[3]-0.125*df_proj1_l[1]*D_proj1_l[3]);
  out[9] +=  Jvxvx*((-0.2165063509461096*D_proj1_u[4]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[4]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[4]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[4]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[4]*D_proj1_u[7]+0.2165063509461096*df_proj1_u[4]*D_proj1_u[7]+0.2165063509461096*f_proj2_l[4]*D_proj1_l[7]+0.2165063509461096*df_proj1_l[4]*D_proj1_l[7]-0.2165063509461096*D_proj1_u[2]*f_proj2_u[6]+0.2165063509461096*D_proj1_l[2]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[2]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[2]*df_proj1_l[6]-0.2165063509461096*f_proj2_u[2]*D_proj1_u[6]+0.2165063509461096*df_proj1_u[2]*D_proj1_u[6]+0.2165063509461096*f_proj2_l[2]*D_proj1_l[6]+0.2165063509461096*df_proj1_l[2]*D_proj1_l[6]-0.2165063509461096*D_proj1_u[1]*f_proj2_u[5]+0.2165063509461096*D_proj1_l[1]*f_proj2_l[5]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[5]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[5]-0.2165063509461096*f_proj2_u[1]*D_proj1_u[5]+0.2165063509461096*df_proj1_u[1]*D_proj1_u[5]+0.2165063509461096*f_proj2_l[1]*D_proj1_l[5]+0.2165063509461096*df_proj1_l[1]*D_proj1_l[5]-0.2165063509461096*D_proj1_u[0]*f_proj2_u[3]+0.2165063509461096*D_proj1_l[0]*f_proj2_l[3]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[3]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[3]-0.2165063509461096*f_proj2_u[0]*D_proj1_u[3]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[3]+0.2165063509461096*f_proj2_l[0]*D_proj1_l[3]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[3]);
  out[10] +=  Jvxvx*(0.125*D_proj1_u[1]*df_proj1_u[7]-0.125*D_proj1_l[1]*df_proj1_l[7]+0.125*df_proj1_u[1]*D_proj1_u[7]-0.125*df_proj1_l[1]*D_proj1_l[7]+0.125*D_proj1_u[0]*df_proj1_u[6]-0.125*D_proj1_l[0]*df_proj1_l[6]+0.125*df_proj1_u[0]*D_proj1_u[6]-0.125*df_proj1_l[0]*D_proj1_l[6]+0.125*D_proj1_u[4]*df_proj1_u[5]-0.125*D_proj1_l[4]*df_proj1_l[5]+0.125*df_proj1_u[4]*D_proj1_u[5]-0.125*df_proj1_l[4]*D_proj1_l[5]+0.125*D_proj1_u[2]*df_proj1_u[3]-0.125*D_proj1_l[2]*df_proj1_l[3]+0.125*df_proj1_u[2]*D_proj1_u[3]-0.125*df_proj1_l[2]*D_proj1_l[3]);
  out[11] +=  Jvxvx*((-0.2165063509461096*D_proj1_u[3]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[3]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[3]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[3]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[3]*D_proj1_u[7]+0.2165063509461096*df_proj1_u[3]*D_proj1_u[7]+0.2165063509461096*f_proj2_l[3]*D_proj1_l[7]+0.2165063509461096*df_proj1_l[3]*D_proj1_l[7]-0.2165063509461096*D_proj1_u[5]*f_proj2_u[6]+0.2165063509461096*D_proj1_l[5]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[5]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[5]*df_proj1_l[6]-0.2165063509461096*f_proj2_u[5]*D_proj1_u[6]+0.2165063509461096*df_proj1_u[5]*D_proj1_u[6]+0.2165063509461096*f_proj2_l[5]*D_proj1_l[6]+0.2165063509461096*df_proj1_l[5]*D_proj1_l[6]-0.2165063509461096*D_proj1_u[0]*f_proj2_u[4]+0.2165063509461096*D_proj1_l[0]*f_proj2_l[4]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[4]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[4]-0.2165063509461096*f_proj2_u[0]*D_proj1_u[4]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[4]+0.2165063509461096*f_proj2_l[0]*D_proj1_l[4]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[4]-0.2165063509461096*D_proj1_u[1]*f_proj2_u[2]+0.2165063509461096*D_proj1_l[1]*f_proj2_l[2]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[2]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[2]-0.2165063509461096*f_proj2_u[1]*D_proj1_u[2]+0.2165063509461096*df_proj1_u[1]*D_proj1_u[2]+0.2165063509461096*f_proj2_l[1]*D_proj1_l[2]+0.2165063509461096*df_proj1_l[1]*D_proj1_l[2]);
  out[12] +=  Jvxvx*((-0.2165063509461096*D_proj1_u[2]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[2]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[2]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[2]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[2]*D_proj1_u[7]+0.2165063509461096*df_proj1_u[2]*D_proj1_u[7]+0.2165063509461096*f_proj2_l[2]*D_proj1_l[7]+0.2165063509461096*df_proj1_l[2]*D_proj1_l[7]-0.2165063509461096*D_proj1_u[4]*f_proj2_u[6]+0.2165063509461096*D_proj1_l[4]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[4]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[4]*df_proj1_l[6]-0.2165063509461096*f_proj2_u[4]*D_proj1_u[6]+0.2165063509461096*df_proj1_u[4]*D_proj1_u[6]+0.2165063509461096*f_proj2_l[4]*D_proj1_l[6]+0.2165063509461096*df_proj1_l[4]*D_proj1_l[6]-0.2165063509461096*D_proj1_u[0]*f_proj2_u[5]+0.2165063509461096*D_proj1_l[0]*f_proj2_l[5]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[5]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[5]-0.2165063509461096*f_proj2_u[0]*D_proj1_u[5]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[5]+0.2165063509461096*f_proj2_l[0]*D_proj1_l[5]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[5]-0.2165063509461096*D_proj1_u[1]*f_proj2_u[3]+0.2165063509461096*D_proj1_l[1]*f_proj2_l[3]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[3]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[3]-0.2165063509461096*f_proj2_u[1]*D_proj1_u[3]+0.2165063509461096*df_proj1_u[1]*D_proj1_u[3]+0.2165063509461096*f_proj2_l[1]*D_proj1_l[3]+0.2165063509461096*df_proj1_l[1]*D_proj1_l[3]);
  out[13] +=  Jvxvx*(0.125*D_proj1_u[0]*df_proj1_u[7]-0.125*D_proj1_l[0]*df_proj1_l[7]+0.125*df_proj1_u[0]*D_proj1_u[7]-0.125*df_proj1_l[0]*D_proj1_l[7]+0.125*D_proj1_u[1]*df_proj1_u[6]-0.125*D_proj1_l[1]*df_proj1_l[6]+0.125*df_proj1_u[1]*D_proj1_u[6]-0.125*df_proj1_l[1]*D_proj1_l[6]+0.125*D_proj1_u[2]*df_proj1_u[5]-0.125*D_proj1_l[2]*df_proj1_l[5]+0.125*df_proj1_u[2]*D_proj1_u[5]-0.125*df_proj1_l[2]*D_proj1_l[5]+0.125*D_proj1_u[3]*df_proj1_u[4]-0.125*D_proj1_l[3]*df_proj1_l[4]+0.125*df_proj1_u[3]*D_proj1_u[4]-0.125*df_proj1_l[3]*D_proj1_l[4]);
  out[14] +=  Jvxvx*((-0.2165063509461096*D_proj1_u[1]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[1]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[1]*D_proj1_u[7]+0.2165063509461096*df_proj1_u[1]*D_proj1_u[7]+0.2165063509461096*f_proj2_l[1]*D_proj1_l[7]+0.2165063509461096*df_proj1_l[1]*D_proj1_l[7]-0.2165063509461096*D_proj1_u[0]*f_proj2_u[6]+0.2165063509461096*D_proj1_l[0]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[6]-0.2165063509461096*f_proj2_u[0]*D_proj1_u[6]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[6]+0.2165063509461096*f_proj2_l[0]*D_proj1_l[6]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[6]-0.2165063509461096*D_proj1_u[4]*f_proj2_u[5]+0.2165063509461096*D_proj1_l[4]*f_proj2_l[5]+0.2165063509461096*D_proj1_u[4]*df_proj1_u[5]+0.2165063509461096*D_proj1_l[4]*df_proj1_l[5]-0.2165063509461096*f_proj2_u[4]*D_proj1_u[5]+0.2165063509461096*df_proj1_u[4]*D_proj1_u[5]+0.2165063509461096*f_proj2_l[4]*D_proj1_l[5]+0.2165063509461096*df_proj1_l[4]*D_proj1_l[5]-0.2165063509461096*D_proj1_u[2]*f_proj2_u[3]+0.2165063509461096*D_proj1_l[2]*f_proj2_l[3]+0.2165063509461096*D_proj1_u[2]*df_proj1_u[3]+0.2165063509461096*D_proj1_l[2]*df_proj1_l[3]-0.2165063509461096*f_proj2_u[2]*D_proj1_u[3]+0.2165063509461096*df_proj1_u[2]*D_proj1_u[3]+0.2165063509461096*f_proj2_l[2]*D_proj1_l[3]+0.2165063509461096*df_proj1_l[2]*D_proj1_l[3]);
  out[15] +=  Jvxvx*((-0.2165063509461096*D_proj1_u[0]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[0]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[0]*D_proj1_u[7]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[7]+0.2165063509461096*f_proj2_l[0]*D_proj1_l[7]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[7]-0.2165063509461096*D_proj1_u[1]*f_proj2_u[6]+0.2165063509461096*D_proj1_l[1]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[6]-0.2165063509461096*f_proj2_u[1]*D_proj1_u[6]+0.2165063509461096*df_proj1_u[1]*D_proj1_u[6]+0.2165063509461096*f_proj2_l[1]*D_proj1_l[6]+0.2165063509461096*df_proj1_l[1]*D_proj1_l[6]-0.2165063509461096*D_proj1_u[2]*f_proj2_u[5]+0.2165063509461096*D_proj1_l[2]*f_proj2_l[5]+0.2165063509461096*D_proj1_u[2]*df_proj1_u[5]+0.2165063509461096*D_proj1_l[2]*df_proj1_l[5]-0.2165063509461096*f_proj2_u[2]*D_proj1_u[5]+0.2165063509461096*df_proj1_u[2]*D_proj1_u[5]+0.2165063509461096*f_proj2_l[2]*D_proj1_l[5]+0.2165063509461096*df_proj1_l[2]*D_proj1_l[5]-0.2165063509461096*D_proj1_u[3]*f_proj2_u[4]+0.2165063509461096*D_proj1_l[3]*f_proj2_l[4]+0.2165063509461096*D_proj1_u[3]*df_proj1_u[4]+0.2165063509461096*D_proj1_l[3]*df_proj1_l[4]-0.2165063509461096*f_proj2_u[3]*D_proj1_u[4]+0.2165063509461096*df_proj1_u[3]*D_proj1_u[4]+0.2165063509461096*f_proj2_l[3]*D_proj1_l[4]+0.2165063509461096*df_proj1_l[3]*D_proj1_l[4]);
  return 0.;

}

