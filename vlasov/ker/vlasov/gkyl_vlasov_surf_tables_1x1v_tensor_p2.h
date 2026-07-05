#pragma once 
// Generated nodal surface-evaluation tables (see surf-eval-tables.mac). 
// Header-static so both host and device instantiations of GKYL_CU_DH 
// kernels can read them; unused arrays are elided per TU. 
__attribute__((unused)) static const double vst_1x1v_tensor_p2_conf_ev[12] = { 
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
__attribute__((unused)) static const double vst_1x1v_tensor_p2_conf_dx0[12] = { 
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
__attribute__((unused)) static const double vst_1x1v_tensor_p2_vel_ev_v0[3] = { 
  0.7071067811865475, 
  -1.224744871391589, 
  1.5811388300841895, 
}; 
__attribute__((unused)) static const int vst_1x1v_tensor_p2_vel_sparse_idx[3] = { 
  0, 
  1, 
  2, 
}; 
__attribute__((unused)) static const double vst_1x1v_tensor_p2_ph_v0_Cm[12] = { 
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
__attribute__((unused)) static const double vst_1x1v_tensor_p2_ph_v0_CmD[12] = { 
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
__attribute__((unused)) static const double vst_1x1v_tensor_p2_ph_v0_V[3] = { 
  1.0, 
  2.5, 
  -3.8461538461538463, 
}; 
__attribute__((unused)) static const double vst_1x1v_tensor_p2_ph_v0_coefl[9] = { 
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
__attribute__((unused)) static const double vst_1x1v_tensor_p2_ph_v0_coefr[9] = { 
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
__attribute__((unused)) static const int vst_1x1v_tensor_p2_ph_v0_cmap[9] = { 
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
__attribute__((unused)) static const int vst_1x1v_tensor_p2_ph_v0_vlmap[9] = { 
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
__attribute__((unused)) static const int vst_1x1v_tensor_p2_ph_v0_vrmap[9] = { 
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
