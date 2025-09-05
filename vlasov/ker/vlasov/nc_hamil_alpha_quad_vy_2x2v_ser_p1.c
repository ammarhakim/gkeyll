#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void nc_hamil_alpha_quad_vy_2x2v_ser_p1(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double dx11 = 2.0/dxv[1]; 

  const double wx1 = w[2]; 
  const double wx2 = w[3]; 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[16]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[20]; 
  alpha_quad[0] += -dx11*((0.6339745962155614*hamil[3]+0.8660254037844386*hamil[2]-0.8660254037844386*hamil[1])*((0.5*poisson_tensor_conf_1[3]-0.5*poisson_tensor_conf_1[2]-0.5*poisson_tensor_conf_1[1]+0.5*poisson_tensor_conf_1[0])*(wx2-0.5*dv11)+(0.5*poisson_tensor_conf_0[3]-0.5*poisson_tensor_conf_0[2]-0.5*poisson_tensor_conf_0[1]+0.5*poisson_tensor_conf_0[0])*(wx1-0.2886751345948129*dv10))); 
  alpha_quad[1] += -dx11*((2.366025403784438*hamil[3]+0.8660254037844386*hamil[2]-0.8660254037844386*hamil[1])*((0.5*poisson_tensor_conf_1[3]-0.5*poisson_tensor_conf_1[2]-0.5*poisson_tensor_conf_1[1]+0.5*poisson_tensor_conf_1[0])*(wx2-0.5*dv11)+(0.5*poisson_tensor_conf_0[3]-0.5*poisson_tensor_conf_0[2]-0.5*poisson_tensor_conf_0[1]+0.5*poisson_tensor_conf_0[0])*(wx1+0.2886751345948129*dv10))); 
  alpha_quad[2] += -dx11*((0.6339745962155614*hamil[3]+0.8660254037844386*hamil[2]-0.8660254037844386*hamil[1])*(((-0.5*poisson_tensor_conf_1[3])+0.5*poisson_tensor_conf_1[2]-0.5*poisson_tensor_conf_1[1]+0.5*poisson_tensor_conf_1[0])*(wx2-0.5*dv11)+((-0.5*poisson_tensor_conf_0[3])+0.5*poisson_tensor_conf_0[2]-0.5*poisson_tensor_conf_0[1]+0.5*poisson_tensor_conf_0[0])*(wx1-0.2886751345948129*dv10))); 
  alpha_quad[3] += -dx11*((2.366025403784438*hamil[3]+0.8660254037844386*hamil[2]-0.8660254037844386*hamil[1])*(((-0.5*poisson_tensor_conf_1[3])+0.5*poisson_tensor_conf_1[2]-0.5*poisson_tensor_conf_1[1]+0.5*poisson_tensor_conf_1[0])*(wx2-0.5*dv11)+((-0.5*poisson_tensor_conf_0[3])+0.5*poisson_tensor_conf_0[2]-0.5*poisson_tensor_conf_0[1]+0.5*poisson_tensor_conf_0[0])*(wx1+0.2886751345948129*dv10))); 
  alpha_quad[4] += -dx11*((0.6339745962155614*hamil[3]+0.8660254037844386*hamil[2]-0.8660254037844386*hamil[1])*(((-0.5*poisson_tensor_conf_1[3])-0.5*poisson_tensor_conf_1[2]+0.5*poisson_tensor_conf_1[1]+0.5*poisson_tensor_conf_1[0])*(wx2-0.5*dv11)+((-0.5*poisson_tensor_conf_0[3])-0.5*poisson_tensor_conf_0[2]+0.5*poisson_tensor_conf_0[1]+0.5*poisson_tensor_conf_0[0])*(wx1-0.2886751345948129*dv10))); 
  alpha_quad[5] += -dx11*((2.366025403784438*hamil[3]+0.8660254037844386*hamil[2]-0.8660254037844386*hamil[1])*(((-0.5*poisson_tensor_conf_1[3])-0.5*poisson_tensor_conf_1[2]+0.5*poisson_tensor_conf_1[1]+0.5*poisson_tensor_conf_1[0])*(wx2-0.5*dv11)+((-0.5*poisson_tensor_conf_0[3])-0.5*poisson_tensor_conf_0[2]+0.5*poisson_tensor_conf_0[1]+0.5*poisson_tensor_conf_0[0])*(wx1+0.2886751345948129*dv10))); 
  alpha_quad[6] += -dx11*((0.6339745962155614*hamil[3]+0.8660254037844386*hamil[2]-0.8660254037844386*hamil[1])*((0.5*poisson_tensor_conf_1[3]+0.5*poisson_tensor_conf_1[2]+0.5*poisson_tensor_conf_1[1]+0.5*poisson_tensor_conf_1[0])*(wx2-0.5*dv11)+(0.5*poisson_tensor_conf_0[3]+0.5*poisson_tensor_conf_0[2]+0.5*poisson_tensor_conf_0[1]+0.5*poisson_tensor_conf_0[0])*(wx1-0.2886751345948129*dv10))); 
  alpha_quad[7] += -dx11*((2.366025403784438*hamil[3]+0.8660254037844386*hamil[2]-0.8660254037844386*hamil[1])*((0.5*poisson_tensor_conf_1[3]+0.5*poisson_tensor_conf_1[2]+0.5*poisson_tensor_conf_1[1]+0.5*poisson_tensor_conf_1[0])*(wx2-0.5*dv11)+(0.5*poisson_tensor_conf_0[3]+0.5*poisson_tensor_conf_0[2]+0.5*poisson_tensor_conf_0[1]+0.5*poisson_tensor_conf_0[0])*(wx1+0.2886751345948129*dv10))); 

} 
