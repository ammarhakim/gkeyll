GKYL_CU_DH static inline double 
ser_3x_p1_surfx3_eval_quad_node_0_r(const double* GKYL_RESTRICT f) { 
  return 0.6123724356957944*f[7]-0.6123724356957944*(f[6]+f[5])+0.3535533905932737*f[4]+0.6123724356957944*f[3]-0.3535533905932737*(f[2]+f[1])+0.3535533905932737*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_3x_p1_surfx3_eval_quad_node_0_l(const double* GKYL_RESTRICT f) { 
    return (-0.6123724356957944*f[7])+0.6123724356957944*(f[6]+f[5])+0.3535533905932737*f[4]-0.6123724356957944*f[3]-0.3535533905932737*(f[2]+f[1])+0.3535533905932737*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_3x_p1_surfx3_eval_quad_node_1_r(const double* GKYL_RESTRICT f) { 
  return (-0.6123724356957944*f[7])+0.6123724356957944*f[6]-0.6123724356957944*f[5]-0.3535533905932737*f[4]+0.6123724356957944*f[3]+0.3535533905932737*f[2]-0.3535533905932737*f[1]+0.3535533905932737*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_3x_p1_surfx3_eval_quad_node_1_l(const double* GKYL_RESTRICT f) { 
    return 0.6123724356957944*f[7]-0.6123724356957944*f[6]+0.6123724356957944*f[5]-0.3535533905932737*f[4]-0.6123724356957944*f[3]+0.3535533905932737*f[2]-0.3535533905932737*f[1]+0.3535533905932737*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_3x_p1_surfx3_eval_quad_node_2_r(const double* GKYL_RESTRICT f) { 
  return (-0.6123724356957944*(f[7]+f[6]))+0.6123724356957944*f[5]-0.3535533905932737*f[4]+0.6123724356957944*f[3]-0.3535533905932737*f[2]+0.3535533905932737*(f[1]+f[0]); 
} 
GKYL_CU_DH static inline double 
ser_3x_p1_surfx3_eval_quad_node_2_l(const double* GKYL_RESTRICT f) { 
    return 0.6123724356957944*(f[7]+f[6])-0.6123724356957944*f[5]-0.3535533905932737*f[4]-0.6123724356957944*f[3]-0.3535533905932737*f[2]+0.3535533905932737*(f[1]+f[0]); 
} 
GKYL_CU_DH static inline double 
ser_3x_p1_surfx3_eval_quad_node_3_r(const double* GKYL_RESTRICT f) { 
  return 0.6123724356957944*(f[7]+f[6]+f[5])+0.3535533905932737*f[4]+0.6123724356957944*f[3]+0.3535533905932737*(f[2]+f[1]+f[0]); 
} 
GKYL_CU_DH static inline double 
ser_3x_p1_surfx3_eval_quad_node_3_l(const double* GKYL_RESTRICT f) { 
    return (-0.6123724356957944*(f[7]+f[6]+f[5]))+0.3535533905932737*f[4]-0.6123724356957944*f[3]+0.3535533905932737*(f[2]+f[1]+f[0]); 
} 
