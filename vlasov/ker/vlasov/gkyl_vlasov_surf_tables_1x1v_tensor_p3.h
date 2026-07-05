#pragma once 
// Generated nodal surface-evaluation tables (see surf-eval-tables.mac). 
// Header-static so both host and device instantiations of GKYL_CU_DH 
// kernels can read them; unused arrays are elided per TU. 
__attribute__((unused)) static const double vst_1x1v_tensor_p3_conf_ev[20] = { 
  0.7071067811865475, 
  -1.1098391188717989, 
  1.1569870650434417, 
  -0.9373434910716291, 
  0.7071067811865475, 
  -0.6594875259537024, 
  -0.10289451165398211, 
  0.7808504157580675, 
  0.7071067811865475, 
  0.0, 
  -0.7905694150420947, 
  0.0, 
  0.7071067811865475, 
  0.6594875259537024, 
  -0.10289451165398211, 
  -0.7808504157580675, 
  0.7071067811865475, 
  1.1098391188717989, 
  1.1569870650434417, 
  0.9373434910716291, 
}; 
__attribute__((unused)) static const double vst_1x1v_tensor_p3_conf_dx0[20] = { 
  0.0, 
  1.224744871391589, 
  -4.298388424359992, 
  8.715656478198236, 
  0.0, 
  1.224744871391589, 
  -2.554184205050222, 
  1.2620965531989388, 
  0.0, 
  1.224744871391589, 
  0.0, 
  -2.806243040080455, 
  0.0, 
  1.224744871391589, 
  2.554184205050222, 
  1.2620965531989388, 
  0.0, 
  1.224744871391589, 
  4.298388424359992, 
  8.715656478198236, 
}; 
__attribute__((unused)) static const double vst_1x1v_tensor_p3_vel_ev_v0[4] = { 
  0.7071067811865475, 
  -1.224744871391589, 
  1.5811388300841895, 
  -1.8708286933869707, 
}; 
__attribute__((unused)) static const int vst_1x1v_tensor_p3_vel_sparse_idx[4] = { 
  0, 
  1, 
  2, 
  3, 
}; 
