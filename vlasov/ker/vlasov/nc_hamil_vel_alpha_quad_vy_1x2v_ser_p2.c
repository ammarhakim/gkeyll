#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void nc_hamil_vel_alpha_quad_vy_1x2v_ser_p2(const double *w, const double *dxv, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 

  const double wx1 = w[1]; 
  const double wx2 = w[2]; 
  const double dv0 = dxv[1]; 
  const double dv1 = dxv[2]; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[12]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[15]; 
  alpha_quad[0] += (1.936491673103709*hamil[7]+4.499999999999999*hamil[6]-2.598076211353316*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1])*dv10*((-1.0*(0.6324555320336759*poisson_tensor_conf_1[2]-0.9486832980505137*poisson_tensor_conf_1[1]+0.7071067811865475*poisson_tensor_conf_1[0])*(wx2-0.5*dv1))-1.0*(0.6324555320336759*poisson_tensor_conf_0[2]-0.9486832980505137*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(wx1-0.3872983346207416*dv0)); 
  alpha_quad[1] += (1.936491673103709*hamil[7]-1.5*hamil[3]+0.8660254037844386*hamil[1])*dv10*((-1.0*(0.6324555320336759*poisson_tensor_conf_1[2]-0.9486832980505137*poisson_tensor_conf_1[1]+0.7071067811865475*poisson_tensor_conf_1[0])*(wx2-0.5*dv1))-1.0*(0.6324555320336759*poisson_tensor_conf_0[2]-0.9486832980505137*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*wx1); 
  alpha_quad[2] += (1.936491673103709*hamil[7]-4.499999999999999*hamil[6]+2.598076211353316*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1])*dv10*((-1.0*(0.6324555320336759*poisson_tensor_conf_1[2]-0.9486832980505137*poisson_tensor_conf_1[1]+0.7071067811865475*poisson_tensor_conf_1[0])*(wx2-0.5*dv1))-1.0*(0.6324555320336759*poisson_tensor_conf_0[2]-0.9486832980505137*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(wx1+0.3872983346207416*dv0)); 
  alpha_quad[3] += (1.936491673103709*hamil[7]+4.499999999999999*hamil[6]-2.598076211353316*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1])*dv10*((-1.0*(0.7071067811865475*poisson_tensor_conf_1[0]-0.7905694150420947*poisson_tensor_conf_1[2])*(wx2-0.5*dv1))-1.0*(0.7071067811865475*poisson_tensor_conf_0[0]-0.7905694150420947*poisson_tensor_conf_0[2])*(wx1-0.3872983346207416*dv0)); 
  alpha_quad[4] += (1.936491673103709*hamil[7]-1.5*hamil[3]+0.8660254037844386*hamil[1])*dv10*((-1.0*(0.7071067811865475*poisson_tensor_conf_1[0]-0.7905694150420947*poisson_tensor_conf_1[2])*(wx2-0.5*dv1))-1.0*(0.7071067811865475*poisson_tensor_conf_0[0]-0.7905694150420947*poisson_tensor_conf_0[2])*wx1); 
  alpha_quad[5] += (1.936491673103709*hamil[7]-4.499999999999999*hamil[6]+2.598076211353316*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1])*dv10*((-1.0*(0.7071067811865475*poisson_tensor_conf_1[0]-0.7905694150420947*poisson_tensor_conf_1[2])*(wx2-0.5*dv1))-1.0*(0.7071067811865475*poisson_tensor_conf_0[0]-0.7905694150420947*poisson_tensor_conf_0[2])*(wx1+0.3872983346207416*dv0)); 
  alpha_quad[6] += (1.936491673103709*hamil[7]+4.499999999999999*hamil[6]-2.598076211353316*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1])*dv10*((-1.0*(0.6324555320336759*poisson_tensor_conf_1[2]+0.9486832980505137*poisson_tensor_conf_1[1]+0.7071067811865475*poisson_tensor_conf_1[0])*(wx2-0.5*dv1))-1.0*(0.6324555320336759*poisson_tensor_conf_0[2]+0.9486832980505137*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(wx1-0.3872983346207416*dv0)); 
  alpha_quad[7] += (1.936491673103709*hamil[7]-1.5*hamil[3]+0.8660254037844386*hamil[1])*dv10*((-1.0*(0.6324555320336759*poisson_tensor_conf_1[2]+0.9486832980505137*poisson_tensor_conf_1[1]+0.7071067811865475*poisson_tensor_conf_1[0])*(wx2-0.5*dv1))-1.0*(0.6324555320336759*poisson_tensor_conf_0[2]+0.9486832980505137*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*wx1); 
  alpha_quad[8] += (1.936491673103709*hamil[7]-4.499999999999999*hamil[6]+2.598076211353316*hamil[4]-1.5*hamil[3]+0.8660254037844386*hamil[1])*dv10*((-1.0*(0.6324555320336759*poisson_tensor_conf_1[2]+0.9486832980505137*poisson_tensor_conf_1[1]+0.7071067811865475*poisson_tensor_conf_1[0])*(wx2-0.5*dv1))-1.0*(0.6324555320336759*poisson_tensor_conf_0[2]+0.9486832980505137*poisson_tensor_conf_0[1]+0.7071067811865475*poisson_tensor_conf_0[0])*(wx1+0.3872983346207416*dv0)); 

} 
