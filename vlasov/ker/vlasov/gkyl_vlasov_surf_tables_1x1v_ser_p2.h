#pragma once 
// Generated nodal surface-evaluation tables (see surf-eval-tables.mac). 
// Header-static so both host and device instantiations of GKYL_CU_DH 
// kernels can read them; unused arrays are elided per TU. 
__attribute__((unused)) static const double vst_1x1v_ser_p2_conf_ev[9] = { 
  0.7071067811865475, 
  -0.9486832980505137, 
  0.6324555320336759, 
  0.7071067811865475, 
  0.0, 
  -0.7905694150420947, 
  0.7071067811865475, 
  0.9486832980505137, 
  0.6324555320336759, 
}; 
__attribute__((unused)) static const double vst_1x1v_ser_p2_conf_dx0[9] = { 
  0.0, 
  1.224744871391589, 
  -3.6742346141747664, 
  0.0, 
  1.224744871391589, 
  0.0, 
  0.0, 
  1.224744871391589, 
  3.6742346141747664, 
}; 
__attribute__((unused)) static const double vst_1x1v_ser_p2_vel_ev_v0[3] = { 
  0.7071067811865475, 
  -1.224744871391589, 
  1.5811388300841895, 
}; 
__attribute__((unused)) static const int vst_1x1v_ser_p2_vel_sparse_idx[3] = { 
  0, 
  1, 
  2, 
}; 
__attribute__((unused)) static const double vst_1x1v_ser_p2_ho_conf_ev[12] = { 
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
__attribute__((unused)) static const double vst_1x1v_ser_p2_ho_conf_dx0[12] = { 
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
__attribute__((unused)) static const double vst_1x1v_ser_p2_ho_vel_ev_v0[3] = { 
  0.7071067811865475, 
  -1.224744871391589, 
  1.5811388300841895, 
}; 
__attribute__((unused)) static const int vst_1x1v_ser_p2_ho_vel_sparse_idx[3] = { 
  0, 
  1, 
  2, 
}; 
