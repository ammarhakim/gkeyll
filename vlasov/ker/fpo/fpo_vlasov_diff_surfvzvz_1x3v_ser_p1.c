#include <gkyl_fpo_vlasov_kernels.h>

GKYL_CU_DH double
fpo_vlasov_diff_surfvzvz_1x3v_ser_p1(const double* w, const double* dx,
  const double* g[], const double* f[], double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinates
  // dxv[NDIM]: Cell spacing
  // g: 
  // f: 
  // out: Incremented output

  const double Jvzvz = 4/dx[3]/dx[3];

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
  D_proj1_l[0] = 0.6123724356957944*gccc[4]-0.6123724356957944*gccl[4];
  D_proj1_l[1] = 0.6123724356957944*gccc[8]-0.6123724356957944*gccl[8];
  D_proj1_l[2] = 0.6123724356957944*gccc[9]-0.6123724356957944*gccl[9];
  D_proj1_l[3] = 0.6123724356957944*gccc[10]-0.6123724356957944*gccl[10];
  D_proj1_l[4] = 0.6123724356957944*gccc[12]-0.6123724356957944*gccl[12];
  D_proj1_l[5] = 0.6123724356957944*gccc[13]-0.6123724356957944*gccl[13];
  D_proj1_l[6] = 0.6123724356957944*gccc[14]-0.6123724356957944*gccl[14];
  D_proj1_l[7] = 0.6123724356957944*gccc[15]-0.6123724356957944*gccl[15];

  double D_proj1_u[8];
  D_proj1_u[0] = 0.6123724356957944*gccu[4]-0.6123724356957944*gccc[4];
  D_proj1_u[1] = 0.6123724356957944*gccu[8]-0.6123724356957944*gccc[8];
  D_proj1_u[2] = 0.6123724356957944*gccu[9]-0.6123724356957944*gccc[9];
  D_proj1_u[3] = 0.6123724356957944*gccu[10]-0.6123724356957944*gccc[10];
  D_proj1_u[4] = 0.6123724356957944*gccu[12]-0.6123724356957944*gccc[12];
  D_proj1_u[5] = 0.6123724356957944*gccu[13]-0.6123724356957944*gccc[13];
  D_proj1_u[6] = 0.6123724356957944*gccu[14]-0.6123724356957944*gccc[14];
  D_proj1_u[7] = 0.6123724356957944*gccu[15]-0.6123724356957944*gccc[15];

  double df_proj1_l[8];
  df_proj1_l[0] = (-0.7654655446197428*fccl[4])-0.7654655446197428*fccc[4]-0.7954951288348656*fccl[0]+0.7954951288348656*fccc[0];
  df_proj1_l[1] = (-0.7654655446197428*fccl[8])-0.7654655446197428*fccc[8]-0.7954951288348656*fccl[1]+0.7954951288348656*fccc[1];
  df_proj1_l[2] = (-0.7654655446197428*fccl[9])-0.7654655446197428*fccc[9]-0.7954951288348656*fccl[2]+0.7954951288348656*fccc[2];
  df_proj1_l[3] = (-0.7654655446197428*fccl[10])-0.7654655446197428*fccc[10]-0.7954951288348656*fccl[3]+0.7954951288348656*fccc[3];
  df_proj1_l[4] = (-0.7654655446197428*fccl[12])-0.7654655446197428*fccc[12]-0.7954951288348656*fccl[5]+0.7954951288348656*fccc[5];
  df_proj1_l[5] = (-0.7654655446197428*fccl[13])-0.7654655446197428*fccc[13]-0.7954951288348656*fccl[6]+0.7954951288348656*fccc[6];
  df_proj1_l[6] = (-0.7654655446197428*fccl[14])-0.7654655446197428*fccc[14]-0.7954951288348656*fccl[7]+0.7954951288348656*fccc[7];
  df_proj1_l[7] = (-0.7654655446197428*fccl[15])-0.7654655446197428*fccc[15]-0.7954951288348656*fccl[11]+0.7954951288348656*fccc[11];

  double df_proj1_u[8];
  df_proj1_u[0] = (-0.7654655446197428*fccu[4])-0.7654655446197428*fccc[4]+0.7954951288348656*fccu[0]-0.7954951288348656*fccc[0];
  df_proj1_u[1] = (-0.7654655446197428*fccu[8])-0.7654655446197428*fccc[8]+0.7954951288348656*fccu[1]-0.7954951288348656*fccc[1];
  df_proj1_u[2] = (-0.7654655446197428*fccu[9])-0.7654655446197428*fccc[9]+0.7954951288348656*fccu[2]-0.7954951288348656*fccc[2];
  df_proj1_u[3] = (-0.7654655446197428*fccu[10])-0.7654655446197428*fccc[10]+0.7954951288348656*fccu[3]-0.7954951288348656*fccc[3];
  df_proj1_u[4] = (-0.7654655446197428*fccu[12])-0.7654655446197428*fccc[12]+0.7954951288348656*fccu[5]-0.7954951288348656*fccc[5];
  df_proj1_u[5] = (-0.7654655446197428*fccu[13])-0.7654655446197428*fccc[13]+0.7954951288348656*fccu[6]-0.7954951288348656*fccc[6];
  df_proj1_u[6] = (-0.7654655446197428*fccu[14])-0.7654655446197428*fccc[14]+0.7954951288348656*fccu[7]-0.7954951288348656*fccc[7];
  df_proj1_u[7] = (-0.7654655446197428*fccu[15])-0.7654655446197428*fccc[15]+0.7954951288348656*fccu[11]-0.7954951288348656*fccc[11];

  double f_proj2_l[8];
  f_proj2_l[0] = 0.408248290463863*fccl[4]-0.408248290463863*fccc[4]+0.3535533905932737*fccl[0]+0.3535533905932737*fccc[0];
  f_proj2_l[1] = 0.408248290463863*fccl[8]-0.408248290463863*fccc[8]+0.3535533905932737*fccl[1]+0.3535533905932737*fccc[1];
  f_proj2_l[2] = 0.408248290463863*fccl[9]-0.408248290463863*fccc[9]+0.3535533905932737*fccl[2]+0.3535533905932737*fccc[2];
  f_proj2_l[3] = 0.408248290463863*fccl[10]-0.408248290463863*fccc[10]+0.3535533905932737*fccl[3]+0.3535533905932737*fccc[3];
  f_proj2_l[4] = 0.408248290463863*fccl[12]-0.408248290463863*fccc[12]+0.3535533905932737*fccl[5]+0.3535533905932737*fccc[5];
  f_proj2_l[5] = 0.408248290463863*fccl[13]-0.408248290463863*fccc[13]+0.3535533905932737*fccl[6]+0.3535533905932737*fccc[6];
  f_proj2_l[6] = 0.408248290463863*fccl[14]-0.408248290463863*fccc[14]+0.3535533905932737*fccl[7]+0.3535533905932737*fccc[7];
  f_proj2_l[7] = 0.408248290463863*fccl[15]-0.408248290463863*fccc[15]+0.3535533905932737*fccl[11]+0.3535533905932737*fccc[11];

  double f_proj2_u[8];
  f_proj2_u[0] = (-0.408248290463863*fccu[4])+0.408248290463863*fccc[4]+0.3535533905932737*fccu[0]+0.3535533905932737*fccc[0];
  f_proj2_u[1] = (-0.408248290463863*fccu[8])+0.408248290463863*fccc[8]+0.3535533905932737*fccu[1]+0.3535533905932737*fccc[1];
  f_proj2_u[2] = (-0.408248290463863*fccu[9])+0.408248290463863*fccc[9]+0.3535533905932737*fccu[2]+0.3535533905932737*fccc[2];
  f_proj2_u[3] = (-0.408248290463863*fccu[10])+0.408248290463863*fccc[10]+0.3535533905932737*fccu[3]+0.3535533905932737*fccc[3];
  f_proj2_u[4] = (-0.408248290463863*fccu[12])+0.408248290463863*fccc[12]+0.3535533905932737*fccu[5]+0.3535533905932737*fccc[5];
  f_proj2_u[5] = (-0.408248290463863*fccu[13])+0.408248290463863*fccc[13]+0.3535533905932737*fccu[6]+0.3535533905932737*fccc[6];
  f_proj2_u[6] = (-0.408248290463863*fccu[14])+0.408248290463863*fccc[14]+0.3535533905932737*fccu[7]+0.3535533905932737*fccc[7];
  f_proj2_u[7] = (-0.408248290463863*fccu[15])+0.408248290463863*fccc[15]+0.3535533905932737*fccu[11]+0.3535533905932737*fccc[11];

  out[0] +=  Jvzvz*(0.125*D_proj1_u[7]*df_proj1_u[7]-0.125*D_proj1_l[7]*df_proj1_l[7]+0.125*D_proj1_u[6]*df_proj1_u[6]-0.125*D_proj1_l[6]*df_proj1_l[6]+0.125*D_proj1_u[5]*df_proj1_u[5]-0.125*D_proj1_l[5]*df_proj1_l[5]+0.125*D_proj1_u[4]*df_proj1_u[4]-0.125*D_proj1_l[4]*df_proj1_l[4]+0.125*D_proj1_u[3]*df_proj1_u[3]-0.125*D_proj1_l[3]*df_proj1_l[3]+0.125*D_proj1_u[2]*df_proj1_u[2]-0.125*D_proj1_l[2]*df_proj1_l[2]+0.125*D_proj1_u[1]*df_proj1_u[1]-0.125*D_proj1_l[1]*df_proj1_l[1]+0.125*D_proj1_u[0]*df_proj1_u[0]-0.125*D_proj1_l[0]*df_proj1_l[0]);
  out[1] +=  Jvzvz*(0.125*D_proj1_u[6]*df_proj1_u[7]-0.125*D_proj1_l[6]*df_proj1_l[7]+0.125*df_proj1_u[6]*D_proj1_u[7]-0.125*df_proj1_l[6]*D_proj1_l[7]+0.125*D_proj1_u[3]*df_proj1_u[5]-0.125*D_proj1_l[3]*df_proj1_l[5]+0.125*df_proj1_u[3]*D_proj1_u[5]-0.125*df_proj1_l[3]*D_proj1_l[5]+0.125*D_proj1_u[2]*df_proj1_u[4]-0.125*D_proj1_l[2]*df_proj1_l[4]+0.125*df_proj1_u[2]*D_proj1_u[4]-0.125*df_proj1_l[2]*D_proj1_l[4]+0.125*D_proj1_u[0]*df_proj1_u[1]-0.125*D_proj1_l[0]*df_proj1_l[1]+0.125*df_proj1_u[0]*D_proj1_u[1]-0.125*df_proj1_l[0]*D_proj1_l[1]);
  out[2] +=  Jvzvz*(0.125*D_proj1_u[5]*df_proj1_u[7]-0.125*D_proj1_l[5]*df_proj1_l[7]+0.125*df_proj1_u[5]*D_proj1_u[7]-0.125*df_proj1_l[5]*D_proj1_l[7]+0.125*D_proj1_u[3]*df_proj1_u[6]-0.125*D_proj1_l[3]*df_proj1_l[6]+0.125*df_proj1_u[3]*D_proj1_u[6]-0.125*df_proj1_l[3]*D_proj1_l[6]+0.125*D_proj1_u[1]*df_proj1_u[4]-0.125*D_proj1_l[1]*df_proj1_l[4]+0.125*df_proj1_u[1]*D_proj1_u[4]-0.125*df_proj1_l[1]*D_proj1_l[4]+0.125*D_proj1_u[0]*df_proj1_u[2]-0.125*D_proj1_l[0]*df_proj1_l[2]+0.125*df_proj1_u[0]*D_proj1_u[2]-0.125*df_proj1_l[0]*D_proj1_l[2]);
  out[3] +=  Jvzvz*(0.125*D_proj1_u[4]*df_proj1_u[7]-0.125*D_proj1_l[4]*df_proj1_l[7]+0.125*df_proj1_u[4]*D_proj1_u[7]-0.125*df_proj1_l[4]*D_proj1_l[7]+0.125*D_proj1_u[2]*df_proj1_u[6]-0.125*D_proj1_l[2]*df_proj1_l[6]+0.125*df_proj1_u[2]*D_proj1_u[6]-0.125*df_proj1_l[2]*D_proj1_l[6]+0.125*D_proj1_u[1]*df_proj1_u[5]-0.125*D_proj1_l[1]*df_proj1_l[5]+0.125*df_proj1_u[1]*D_proj1_u[5]-0.125*df_proj1_l[1]*D_proj1_l[5]+0.125*D_proj1_u[0]*df_proj1_u[3]-0.125*D_proj1_l[0]*df_proj1_l[3]+0.125*df_proj1_u[0]*D_proj1_u[3]-0.125*df_proj1_l[0]*D_proj1_l[3]);
  out[4] +=  Jvzvz*((-0.2165063509461096*D_proj1_u[7]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[7]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[7]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[7]*df_proj1_l[7]-0.2165063509461096*D_proj1_u[6]*f_proj2_u[6]+0.2165063509461096*D_proj1_l[6]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[6]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[6]*df_proj1_l[6]-0.2165063509461096*D_proj1_u[5]*f_proj2_u[5]+0.2165063509461096*D_proj1_l[5]*f_proj2_l[5]+0.2165063509461096*D_proj1_u[5]*df_proj1_u[5]+0.2165063509461096*D_proj1_l[5]*df_proj1_l[5]-0.2165063509461096*D_proj1_u[4]*f_proj2_u[4]+0.2165063509461096*D_proj1_l[4]*f_proj2_l[4]+0.2165063509461096*D_proj1_u[4]*df_proj1_u[4]+0.2165063509461096*D_proj1_l[4]*df_proj1_l[4]-0.2165063509461096*D_proj1_u[3]*f_proj2_u[3]+0.2165063509461096*D_proj1_l[3]*f_proj2_l[3]+0.2165063509461096*D_proj1_u[3]*df_proj1_u[3]+0.2165063509461096*D_proj1_l[3]*df_proj1_l[3]-0.2165063509461096*D_proj1_u[2]*f_proj2_u[2]+0.2165063509461096*D_proj1_l[2]*f_proj2_l[2]+0.2165063509461096*D_proj1_u[2]*df_proj1_u[2]+0.2165063509461096*D_proj1_l[2]*df_proj1_l[2]-0.2165063509461096*D_proj1_u[1]*f_proj2_u[1]+0.2165063509461096*D_proj1_l[1]*f_proj2_l[1]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[1]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[1]-0.2165063509461096*D_proj1_u[0]*f_proj2_u[0]+0.2165063509461096*D_proj1_l[0]*f_proj2_l[0]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[0]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[0]);
  out[5] +=  Jvzvz*(0.125*D_proj1_u[3]*df_proj1_u[7]-0.125*D_proj1_l[3]*df_proj1_l[7]+0.125*df_proj1_u[3]*D_proj1_u[7]-0.125*df_proj1_l[3]*D_proj1_l[7]+0.125*D_proj1_u[5]*df_proj1_u[6]-0.125*D_proj1_l[5]*df_proj1_l[6]+0.125*df_proj1_u[5]*D_proj1_u[6]-0.125*df_proj1_l[5]*D_proj1_l[6]+0.125*D_proj1_u[0]*df_proj1_u[4]-0.125*D_proj1_l[0]*df_proj1_l[4]+0.125*df_proj1_u[0]*D_proj1_u[4]-0.125*df_proj1_l[0]*D_proj1_l[4]+0.125*D_proj1_u[1]*df_proj1_u[2]-0.125*D_proj1_l[1]*df_proj1_l[2]+0.125*df_proj1_u[1]*D_proj1_u[2]-0.125*df_proj1_l[1]*D_proj1_l[2]);
  out[6] +=  Jvzvz*(0.125*D_proj1_u[2]*df_proj1_u[7]-0.125*D_proj1_l[2]*df_proj1_l[7]+0.125*df_proj1_u[2]*D_proj1_u[7]-0.125*df_proj1_l[2]*D_proj1_l[7]+0.125*D_proj1_u[4]*df_proj1_u[6]-0.125*D_proj1_l[4]*df_proj1_l[6]+0.125*df_proj1_u[4]*D_proj1_u[6]-0.125*df_proj1_l[4]*D_proj1_l[6]+0.125*D_proj1_u[0]*df_proj1_u[5]-0.125*D_proj1_l[0]*df_proj1_l[5]+0.125*df_proj1_u[0]*D_proj1_u[5]-0.125*df_proj1_l[0]*D_proj1_l[5]+0.125*D_proj1_u[1]*df_proj1_u[3]-0.125*D_proj1_l[1]*df_proj1_l[3]+0.125*df_proj1_u[1]*D_proj1_u[3]-0.125*df_proj1_l[1]*D_proj1_l[3]);
  out[7] +=  Jvzvz*(0.125*D_proj1_u[1]*df_proj1_u[7]-0.125*D_proj1_l[1]*df_proj1_l[7]+0.125*df_proj1_u[1]*D_proj1_u[7]-0.125*df_proj1_l[1]*D_proj1_l[7]+0.125*D_proj1_u[0]*df_proj1_u[6]-0.125*D_proj1_l[0]*df_proj1_l[6]+0.125*df_proj1_u[0]*D_proj1_u[6]-0.125*df_proj1_l[0]*D_proj1_l[6]+0.125*D_proj1_u[4]*df_proj1_u[5]-0.125*D_proj1_l[4]*df_proj1_l[5]+0.125*df_proj1_u[4]*D_proj1_u[5]-0.125*df_proj1_l[4]*D_proj1_l[5]+0.125*D_proj1_u[2]*df_proj1_u[3]-0.125*D_proj1_l[2]*df_proj1_l[3]+0.125*df_proj1_u[2]*D_proj1_u[3]-0.125*df_proj1_l[2]*D_proj1_l[3]);
  out[8] +=  Jvzvz*((-0.2165063509461096*D_proj1_u[6]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[6]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[6]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[6]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[6]*D_proj1_u[7]+0.2165063509461096*df_proj1_u[6]*D_proj1_u[7]+0.2165063509461096*f_proj2_l[6]*D_proj1_l[7]+0.2165063509461096*df_proj1_l[6]*D_proj1_l[7]-0.2165063509461096*D_proj1_u[3]*f_proj2_u[5]+0.2165063509461096*D_proj1_l[3]*f_proj2_l[5]+0.2165063509461096*D_proj1_u[3]*df_proj1_u[5]+0.2165063509461096*D_proj1_l[3]*df_proj1_l[5]-0.2165063509461096*f_proj2_u[3]*D_proj1_u[5]+0.2165063509461096*df_proj1_u[3]*D_proj1_u[5]+0.2165063509461096*f_proj2_l[3]*D_proj1_l[5]+0.2165063509461096*df_proj1_l[3]*D_proj1_l[5]-0.2165063509461096*D_proj1_u[2]*f_proj2_u[4]+0.2165063509461096*D_proj1_l[2]*f_proj2_l[4]+0.2165063509461096*D_proj1_u[2]*df_proj1_u[4]+0.2165063509461096*D_proj1_l[2]*df_proj1_l[4]-0.2165063509461096*f_proj2_u[2]*D_proj1_u[4]+0.2165063509461096*df_proj1_u[2]*D_proj1_u[4]+0.2165063509461096*f_proj2_l[2]*D_proj1_l[4]+0.2165063509461096*df_proj1_l[2]*D_proj1_l[4]-0.2165063509461096*D_proj1_u[0]*f_proj2_u[1]+0.2165063509461096*D_proj1_l[0]*f_proj2_l[1]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[1]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[1]-0.2165063509461096*f_proj2_u[0]*D_proj1_u[1]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[1]+0.2165063509461096*f_proj2_l[0]*D_proj1_l[1]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[1]);
  out[9] +=  Jvzvz*((-0.2165063509461096*D_proj1_u[5]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[5]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[5]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[5]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[5]*D_proj1_u[7]+0.2165063509461096*df_proj1_u[5]*D_proj1_u[7]+0.2165063509461096*f_proj2_l[5]*D_proj1_l[7]+0.2165063509461096*df_proj1_l[5]*D_proj1_l[7]-0.2165063509461096*D_proj1_u[3]*f_proj2_u[6]+0.2165063509461096*D_proj1_l[3]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[3]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[3]*df_proj1_l[6]-0.2165063509461096*f_proj2_u[3]*D_proj1_u[6]+0.2165063509461096*df_proj1_u[3]*D_proj1_u[6]+0.2165063509461096*f_proj2_l[3]*D_proj1_l[6]+0.2165063509461096*df_proj1_l[3]*D_proj1_l[6]-0.2165063509461096*D_proj1_u[1]*f_proj2_u[4]+0.2165063509461096*D_proj1_l[1]*f_proj2_l[4]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[4]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[4]-0.2165063509461096*f_proj2_u[1]*D_proj1_u[4]+0.2165063509461096*df_proj1_u[1]*D_proj1_u[4]+0.2165063509461096*f_proj2_l[1]*D_proj1_l[4]+0.2165063509461096*df_proj1_l[1]*D_proj1_l[4]-0.2165063509461096*D_proj1_u[0]*f_proj2_u[2]+0.2165063509461096*D_proj1_l[0]*f_proj2_l[2]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[2]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[2]-0.2165063509461096*f_proj2_u[0]*D_proj1_u[2]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[2]+0.2165063509461096*f_proj2_l[0]*D_proj1_l[2]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[2]);
  out[10] +=  Jvzvz*((-0.2165063509461096*D_proj1_u[4]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[4]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[4]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[4]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[4]*D_proj1_u[7]+0.2165063509461096*df_proj1_u[4]*D_proj1_u[7]+0.2165063509461096*f_proj2_l[4]*D_proj1_l[7]+0.2165063509461096*df_proj1_l[4]*D_proj1_l[7]-0.2165063509461096*D_proj1_u[2]*f_proj2_u[6]+0.2165063509461096*D_proj1_l[2]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[2]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[2]*df_proj1_l[6]-0.2165063509461096*f_proj2_u[2]*D_proj1_u[6]+0.2165063509461096*df_proj1_u[2]*D_proj1_u[6]+0.2165063509461096*f_proj2_l[2]*D_proj1_l[6]+0.2165063509461096*df_proj1_l[2]*D_proj1_l[6]-0.2165063509461096*D_proj1_u[1]*f_proj2_u[5]+0.2165063509461096*D_proj1_l[1]*f_proj2_l[5]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[5]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[5]-0.2165063509461096*f_proj2_u[1]*D_proj1_u[5]+0.2165063509461096*df_proj1_u[1]*D_proj1_u[5]+0.2165063509461096*f_proj2_l[1]*D_proj1_l[5]+0.2165063509461096*df_proj1_l[1]*D_proj1_l[5]-0.2165063509461096*D_proj1_u[0]*f_proj2_u[3]+0.2165063509461096*D_proj1_l[0]*f_proj2_l[3]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[3]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[3]-0.2165063509461096*f_proj2_u[0]*D_proj1_u[3]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[3]+0.2165063509461096*f_proj2_l[0]*D_proj1_l[3]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[3]);
  out[11] +=  Jvzvz*(0.125*D_proj1_u[0]*df_proj1_u[7]-0.125*D_proj1_l[0]*df_proj1_l[7]+0.125*df_proj1_u[0]*D_proj1_u[7]-0.125*df_proj1_l[0]*D_proj1_l[7]+0.125*D_proj1_u[1]*df_proj1_u[6]-0.125*D_proj1_l[1]*df_proj1_l[6]+0.125*df_proj1_u[1]*D_proj1_u[6]-0.125*df_proj1_l[1]*D_proj1_l[6]+0.125*D_proj1_u[2]*df_proj1_u[5]-0.125*D_proj1_l[2]*df_proj1_l[5]+0.125*df_proj1_u[2]*D_proj1_u[5]-0.125*df_proj1_l[2]*D_proj1_l[5]+0.125*D_proj1_u[3]*df_proj1_u[4]-0.125*D_proj1_l[3]*df_proj1_l[4]+0.125*df_proj1_u[3]*D_proj1_u[4]-0.125*df_proj1_l[3]*D_proj1_l[4]);
  out[12] +=  Jvzvz*((-0.2165063509461096*D_proj1_u[3]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[3]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[3]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[3]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[3]*D_proj1_u[7]+0.2165063509461096*df_proj1_u[3]*D_proj1_u[7]+0.2165063509461096*f_proj2_l[3]*D_proj1_l[7]+0.2165063509461096*df_proj1_l[3]*D_proj1_l[7]-0.2165063509461096*D_proj1_u[5]*f_proj2_u[6]+0.2165063509461096*D_proj1_l[5]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[5]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[5]*df_proj1_l[6]-0.2165063509461096*f_proj2_u[5]*D_proj1_u[6]+0.2165063509461096*df_proj1_u[5]*D_proj1_u[6]+0.2165063509461096*f_proj2_l[5]*D_proj1_l[6]+0.2165063509461096*df_proj1_l[5]*D_proj1_l[6]-0.2165063509461096*D_proj1_u[0]*f_proj2_u[4]+0.2165063509461096*D_proj1_l[0]*f_proj2_l[4]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[4]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[4]-0.2165063509461096*f_proj2_u[0]*D_proj1_u[4]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[4]+0.2165063509461096*f_proj2_l[0]*D_proj1_l[4]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[4]-0.2165063509461096*D_proj1_u[1]*f_proj2_u[2]+0.2165063509461096*D_proj1_l[1]*f_proj2_l[2]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[2]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[2]-0.2165063509461096*f_proj2_u[1]*D_proj1_u[2]+0.2165063509461096*df_proj1_u[1]*D_proj1_u[2]+0.2165063509461096*f_proj2_l[1]*D_proj1_l[2]+0.2165063509461096*df_proj1_l[1]*D_proj1_l[2]);
  out[13] +=  Jvzvz*((-0.2165063509461096*D_proj1_u[2]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[2]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[2]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[2]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[2]*D_proj1_u[7]+0.2165063509461096*df_proj1_u[2]*D_proj1_u[7]+0.2165063509461096*f_proj2_l[2]*D_proj1_l[7]+0.2165063509461096*df_proj1_l[2]*D_proj1_l[7]-0.2165063509461096*D_proj1_u[4]*f_proj2_u[6]+0.2165063509461096*D_proj1_l[4]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[4]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[4]*df_proj1_l[6]-0.2165063509461096*f_proj2_u[4]*D_proj1_u[6]+0.2165063509461096*df_proj1_u[4]*D_proj1_u[6]+0.2165063509461096*f_proj2_l[4]*D_proj1_l[6]+0.2165063509461096*df_proj1_l[4]*D_proj1_l[6]-0.2165063509461096*D_proj1_u[0]*f_proj2_u[5]+0.2165063509461096*D_proj1_l[0]*f_proj2_l[5]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[5]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[5]-0.2165063509461096*f_proj2_u[0]*D_proj1_u[5]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[5]+0.2165063509461096*f_proj2_l[0]*D_proj1_l[5]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[5]-0.2165063509461096*D_proj1_u[1]*f_proj2_u[3]+0.2165063509461096*D_proj1_l[1]*f_proj2_l[3]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[3]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[3]-0.2165063509461096*f_proj2_u[1]*D_proj1_u[3]+0.2165063509461096*df_proj1_u[1]*D_proj1_u[3]+0.2165063509461096*f_proj2_l[1]*D_proj1_l[3]+0.2165063509461096*df_proj1_l[1]*D_proj1_l[3]);
  out[14] +=  Jvzvz*((-0.2165063509461096*D_proj1_u[1]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[1]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[1]*D_proj1_u[7]+0.2165063509461096*df_proj1_u[1]*D_proj1_u[7]+0.2165063509461096*f_proj2_l[1]*D_proj1_l[7]+0.2165063509461096*df_proj1_l[1]*D_proj1_l[7]-0.2165063509461096*D_proj1_u[0]*f_proj2_u[6]+0.2165063509461096*D_proj1_l[0]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[6]-0.2165063509461096*f_proj2_u[0]*D_proj1_u[6]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[6]+0.2165063509461096*f_proj2_l[0]*D_proj1_l[6]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[6]-0.2165063509461096*D_proj1_u[4]*f_proj2_u[5]+0.2165063509461096*D_proj1_l[4]*f_proj2_l[5]+0.2165063509461096*D_proj1_u[4]*df_proj1_u[5]+0.2165063509461096*D_proj1_l[4]*df_proj1_l[5]-0.2165063509461096*f_proj2_u[4]*D_proj1_u[5]+0.2165063509461096*df_proj1_u[4]*D_proj1_u[5]+0.2165063509461096*f_proj2_l[4]*D_proj1_l[5]+0.2165063509461096*df_proj1_l[4]*D_proj1_l[5]-0.2165063509461096*D_proj1_u[2]*f_proj2_u[3]+0.2165063509461096*D_proj1_l[2]*f_proj2_l[3]+0.2165063509461096*D_proj1_u[2]*df_proj1_u[3]+0.2165063509461096*D_proj1_l[2]*df_proj1_l[3]-0.2165063509461096*f_proj2_u[2]*D_proj1_u[3]+0.2165063509461096*df_proj1_u[2]*D_proj1_u[3]+0.2165063509461096*f_proj2_l[2]*D_proj1_l[3]+0.2165063509461096*df_proj1_l[2]*D_proj1_l[3]);
  out[15] +=  Jvzvz*((-0.2165063509461096*D_proj1_u[0]*f_proj2_u[7])+0.2165063509461096*D_proj1_l[0]*f_proj2_l[7]+0.2165063509461096*D_proj1_u[0]*df_proj1_u[7]+0.2165063509461096*D_proj1_l[0]*df_proj1_l[7]-0.2165063509461096*f_proj2_u[0]*D_proj1_u[7]+0.2165063509461096*df_proj1_u[0]*D_proj1_u[7]+0.2165063509461096*f_proj2_l[0]*D_proj1_l[7]+0.2165063509461096*df_proj1_l[0]*D_proj1_l[7]-0.2165063509461096*D_proj1_u[1]*f_proj2_u[6]+0.2165063509461096*D_proj1_l[1]*f_proj2_l[6]+0.2165063509461096*D_proj1_u[1]*df_proj1_u[6]+0.2165063509461096*D_proj1_l[1]*df_proj1_l[6]-0.2165063509461096*f_proj2_u[1]*D_proj1_u[6]+0.2165063509461096*df_proj1_u[1]*D_proj1_u[6]+0.2165063509461096*f_proj2_l[1]*D_proj1_l[6]+0.2165063509461096*df_proj1_l[1]*D_proj1_l[6]-0.2165063509461096*D_proj1_u[2]*f_proj2_u[5]+0.2165063509461096*D_proj1_l[2]*f_proj2_l[5]+0.2165063509461096*D_proj1_u[2]*df_proj1_u[5]+0.2165063509461096*D_proj1_l[2]*df_proj1_l[5]-0.2165063509461096*f_proj2_u[2]*D_proj1_u[5]+0.2165063509461096*df_proj1_u[2]*D_proj1_u[5]+0.2165063509461096*f_proj2_l[2]*D_proj1_l[5]+0.2165063509461096*df_proj1_l[2]*D_proj1_l[5]-0.2165063509461096*D_proj1_u[3]*f_proj2_u[4]+0.2165063509461096*D_proj1_l[3]*f_proj2_l[4]+0.2165063509461096*D_proj1_u[3]*df_proj1_u[4]+0.2165063509461096*D_proj1_l[3]*df_proj1_l[4]-0.2165063509461096*f_proj2_u[3]*D_proj1_u[4]+0.2165063509461096*df_proj1_u[3]*D_proj1_u[4]+0.2165063509461096*f_proj2_l[3]*D_proj1_l[4]+0.2165063509461096*df_proj1_l[3]*D_proj1_l[4]);
  return 0.;

}

