#include <gkyl_nc_hamil_vol_comps_kernels.h> 
#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH static void dot_pb_1x2v_tensor_p1(const double *list_1, const double *list_2, double* GKYL_RESTRICT out) 
{ 
  for (int i = 0; i<18; ++i) out[0] += list_1[i]*list_2[i]; 
} 
  
GKYL_CU_DH static inline void zero_pb_list(double *pb_list) 
{ 
  for (int i = 0; i<18; ++i) pb_list[i] = 0.0; 
} 
  
GKYL_CU_DH void vlasov_nc_hamil_gen_vol_1x2v_tensor_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  double alphai[18] = {0.0}; 
  double bracket_w_Hi[18] = {0.0}; 
  for (int dir = 0; dir<3; ++dir) { 
    // Compute the per-direction charateristic 
    zero_pb_list(alphai); 
    vlasov_nc_hamil_gen_vol_alpha_1x2v_tensor_p1(w, dxv, dir, poisson_tensor_conf, hamil, alphai); 
    
    // Compute the modal basis backet coeffcients, accumulate the outputs 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp0_1x2v_tensor_p1(w, dxv, dir, alphai, bracket_w_Hi); 
    dot_pb_1x2v_tensor_p1(f, bracket_w_Hi, &out[0]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp1_1x2v_tensor_p1(w, dxv, dir, alphai, bracket_w_Hi); 
    dot_pb_1x2v_tensor_p1(f, bracket_w_Hi, &out[1]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp2_1x2v_tensor_p1(w, dxv, dir, alphai, bracket_w_Hi); 
    dot_pb_1x2v_tensor_p1(f, bracket_w_Hi, &out[2]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp3_1x2v_tensor_p1(w, dxv, dir, alphai, bracket_w_Hi); 
    dot_pb_1x2v_tensor_p1(f, bracket_w_Hi, &out[3]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp4_1x2v_tensor_p1(w, dxv, dir, alphai, bracket_w_Hi); 
    dot_pb_1x2v_tensor_p1(f, bracket_w_Hi, &out[4]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp5_1x2v_tensor_p1(w, dxv, dir, alphai, bracket_w_Hi); 
    dot_pb_1x2v_tensor_p1(f, bracket_w_Hi, &out[5]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp6_1x2v_tensor_p1(w, dxv, dir, alphai, bracket_w_Hi); 
    dot_pb_1x2v_tensor_p1(f, bracket_w_Hi, &out[6]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp7_1x2v_tensor_p1(w, dxv, dir, alphai, bracket_w_Hi); 
    dot_pb_1x2v_tensor_p1(f, bracket_w_Hi, &out[7]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp8_1x2v_tensor_p1(w, dxv, dir, alphai, bracket_w_Hi); 
    dot_pb_1x2v_tensor_p1(f, bracket_w_Hi, &out[8]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp9_1x2v_tensor_p1(w, dxv, dir, alphai, bracket_w_Hi); 
    dot_pb_1x2v_tensor_p1(f, bracket_w_Hi, &out[9]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp10_1x2v_tensor_p1(w, dxv, dir, alphai, bracket_w_Hi); 
    dot_pb_1x2v_tensor_p1(f, bracket_w_Hi, &out[10]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp11_1x2v_tensor_p1(w, dxv, dir, alphai, bracket_w_Hi); 
    dot_pb_1x2v_tensor_p1(f, bracket_w_Hi, &out[11]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp12_1x2v_tensor_p1(w, dxv, dir, alphai, bracket_w_Hi); 
    dot_pb_1x2v_tensor_p1(f, bracket_w_Hi, &out[12]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp13_1x2v_tensor_p1(w, dxv, dir, alphai, bracket_w_Hi); 
    dot_pb_1x2v_tensor_p1(f, bracket_w_Hi, &out[13]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp14_1x2v_tensor_p1(w, dxv, dir, alphai, bracket_w_Hi); 
    dot_pb_1x2v_tensor_p1(f, bracket_w_Hi, &out[14]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp15_1x2v_tensor_p1(w, dxv, dir, alphai, bracket_w_Hi); 
    dot_pb_1x2v_tensor_p1(f, bracket_w_Hi, &out[15]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp16_1x2v_tensor_p1(w, dxv, dir, alphai, bracket_w_Hi); 
    dot_pb_1x2v_tensor_p1(f, bracket_w_Hi, &out[16]); 
    zero_pb_list(bracket_w_Hi); 
    vlasov_nc_hamil_gen_vol_comp17_1x2v_tensor_p1(w, dxv, dir, alphai, bracket_w_Hi); 
    dot_pb_1x2v_tensor_p1(f, bracket_w_Hi, &out[17]); 
    
  } 
} 
