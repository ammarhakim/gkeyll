#pragma once 
// Generated nodal surface-evaluation tables (see surf-eval-tables.mac). 
// Header-static so both host and device instantiations of GKYL_CU_DH 
// kernels can read them; unused arrays are elided per TU. 
#include <gkyl_util.h> 
#include <gkyl_vlasov_surf_meta.h> 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p2_conf_ev[12] = { 
  0.7071067811865475, 
  -1.0546722811938851, 
  0.9681844646844028, 
  0.7071067811865475, 
  -0.41639003950091297, 
  -0.5164305132317774, 
  0.7071067811865475, 
  0.41639003950091297, 
  -0.5164305132317774, 
  0.7071067811865475, 
  1.0546722811938851, 
  0.9681844646844028, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p2_conf_dx0[12] = { 
  0.0, 
  1.224744871391589, 
  -4.084728180770504, 
  0.0, 
  1.224744871391589, 
  -1.6126716885136845, 
  0.0, 
  1.224744871391589, 
  1.6126716885136845, 
  0.0, 
  1.224744871391589, 
  4.084728180770504, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p2_vel_ev_v0[3] = { 
  0.7071067811865475, 
  -1.224744871391589, 
  1.5811388300841895, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_1x1v_tensor_p2_vel_sparse_idx[3] = { 
  0, 
  1, 
  2, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p2_ph_v0_Cm[12] = { 
  1.0, 
  -2.1528407789851314, 
  -2.355129310456592, 
  1.0, 
  -0.8499526089621406, 
  1.2562282115554928, 
  1.0, 
  0.8499526089621406, 
  1.2562282115554928, 
  1.0, 
  2.1528407789851314, 
  -2.355129310456592, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p2_ph_v0_CmD[12] = { 
  0.0, 
  2.5, 
  9.936188210700607, 
  0.0, 
  2.5, 
  3.92285819520988, 
  0.0, 
  2.5, 
  -3.92285819520988, 
  0.0, 
  2.5, 
  -9.936188210700607, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p2_ph_v0_CmDx0[12] = { 
  0.0, 
  2.5, 
  9.936188210700607, 
  0.0, 
  2.5, 
  3.92285819520988, 
  0.0, 
  2.5, 
  -3.92285819520988, 
  0.0, 
  2.5, 
  -9.936188210700607, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p2_ph_v0_V[3] = { 
  1.0, 
  2.5, 
  -3.8461538461538463, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p2_ph_v0_coefl[9] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  -0.2906888370749726, 
  -0.20139513400278566, 
  -0.20139513400278566, 
  0.169, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p2_ph_v0_coefr[9] = { 
  0.5, 
  0.34641016151377546, 
  -0.34641016151377546, 
  -0.24, 
  -0.2906888370749726, 
  -0.2906888370749726, 
  0.20139513400278566, 
  -0.20139513400278566, 
  0.169, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_1x1v_tensor_p2_ph_v0_cmap[9] = { 
  0, 
  1, 
  0, 
  1, 
  2, 
  0, 
  2, 
  1, 
  2, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_1x1v_tensor_p2_ph_v0_vlmap[9] = { 
  0, 
  0, 
  1, 
  1, 
  0, 
  2, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_1x1v_tensor_p2_ph_v0_vrmap[9] = { 
  0, 
  0, 
  1, 
  1, 
  0, 
  2, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p2_ph_v0_Vd0[3] = { 
  0.0, 
  2.5, 
  -11.538461538461538, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_1x1v_tensor_p2_ph_v0_vld0map[9] = { 
  0, 
  0, 
  1, 
  1, 
  0, 
  2, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_1x1v_tensor_p2_ph_v0_vrd0map[9] = { 
  0, 
  0, 
  1, 
  1, 
  0, 
  2, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p2_ph_v0_dcoefl0[9] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  -0.2906888370749726, 
  -0.20139513400278566, 
  -0.20139513400278566, 
  0.169, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p2_ph_v0_dcoefr0[9] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  0.2906888370749726, 
  -0.20139513400278566, 
  0.20139513400278566, 
  -0.169, 
}; 
GKYL_CU_D __attribute__((unused)) static const struct gkyl_vlasov_surf_meta vst_1x1v_tensor_p2_ph_v0_meta = { .nb = 9, .na = 3, .nmv = 3, .no = 4, .ni = 1, .ns = 3, .fhat_off = 0 }; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p2_prj_v0_Cw[12] = { 
  0.3478548451374539, 
  -0.7488760957794683, 
  -0.8192431415675561, 
  0.6521451548625461, 
  -0.5542924757974403, 
  0.8192431415675563, 
  0.6521451548625461, 
  0.5542924757974403, 
  0.8192431415675563, 
  0.3478548451374539, 
  0.7488760957794683, 
  -0.8192431415675561, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p2_prj_v0_Vw[1] = { 
  1.0, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_1x1v_tensor_p2_prj_v0_kamap[3] = { 
  0, 
  1, 
  2, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_1x1v_tensor_p2_prj_v0_kbmap[3] = { 
  0, 
  0, 
  0, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_1x1v_tensor_p2_prj_v0_out_off[4] = { 
  0, 
  3, 
  6, 
  9, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_1x1v_tensor_p2_prj_v0_out_mode[9] = { 
  0, 
  2, 
  5, 
  1, 
  3, 
  7, 
  4, 
  6, 
  8, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p2_prj_v0_out_cl[9] = { 
  0.5, 
  -0.8660254037844386, 
  1.118033988749895, 
  0.34641016151377546, 
  -0.6, 
  0.7745966692414833, 
  -0.2906888370749726, 
  0.503487835006964, 
  -0.65, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p2_prj_v0_out_cr[9] = { 
  -0.5, 
  -0.8660254037844386, 
  -1.118033988749895, 
  -0.34641016151377546, 
  -0.6, 
  -0.7745966692414833, 
  0.2906888370749726, 
  0.503487835006964, 
  0.65, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_1x1v_tensor_p2_prj_v0_pm_a[9] = { 
  0, 
  1, 
  0, 
  1, 
  2, 
  0, 
  2, 
  1, 
  2, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_1x1v_tensor_p2_prj_v0_pm_b[9] = { 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p2_prj_v0_pm_cl[9] = { 
  0.5, 
  0.34641016151377546, 
  -0.8660254037844386, 
  -0.6, 
  -0.2906888370749726, 
  1.118033988749895, 
  0.503487835006964, 
  0.7745966692414833, 
  -0.65, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p2_prj_v0_pm_cr[9] = { 
  -0.5, 
  -0.34641016151377546, 
  -0.8660254037844386, 
  -0.6, 
  0.2906888370749726, 
  -1.118033988749895, 
  0.503487835006964, 
  -0.7745966692414833, 
  0.65, 
}; 
GKYL_CU_D __attribute__((unused)) static const struct gkyl_vlasov_prj_meta vst_1x1v_tensor_p2_prj_v0_meta = { .nk = 3, .np = 9, .na = 3, .nmv = 1, .no = 4, .ni = 1, .fhat_off = 0 }; 
