#pragma once 
// Generated nodal surface-evaluation tables (see surf-eval-tables.mac). 
// Header-static so both host and device instantiations of GKYL_CU_DH 
// kernels can read them; unused arrays are elided per TU. 
__attribute__((unused)) static const double vst_2x2v_ser_p1_conf_ev[16] = { 
  0.5, 
  -0.5, 
  -0.5, 
  0.5, 
  0.5, 
  -0.5, 
  0.5, 
  -0.5, 
  0.5, 
  0.5, 
  -0.5, 
  -0.5, 
  0.5, 
  0.5, 
  0.5, 
  0.5, 
}; 
__attribute__((unused)) static const double vst_2x2v_ser_p1_conf_dx0[16] = { 
  0.0, 
  0.8660254037844386, 
  0.0, 
  -0.8660254037844386, 
  0.0, 
  0.8660254037844386, 
  0.0, 
  0.8660254037844386, 
  0.0, 
  0.8660254037844386, 
  0.0, 
  -0.8660254037844386, 
  0.0, 
  0.8660254037844386, 
  0.0, 
  0.8660254037844386, 
}; 
__attribute__((unused)) static const double vst_2x2v_ser_p1_conf_dx1[16] = { 
  0.0, 
  0.0, 
  0.8660254037844386, 
  -0.8660254037844386, 
  0.0, 
  0.0, 
  0.8660254037844386, 
  -0.8660254037844386, 
  0.0, 
  0.0, 
  0.8660254037844386, 
  0.8660254037844386, 
  0.0, 
  0.0, 
  0.8660254037844386, 
  0.8660254037844386, 
}; 
__attribute__((unused)) static const double vst_2x2v_ser_p1_vel_ev_v0[8] = { 
  0.5, 
  -0.8660254037844386, 
  -0.5, 
  0.8660254037844386, 
  0.5, 
  -0.8660254037844386, 
  0.5, 
  -0.8660254037844386, 
}; 
__attribute__((unused)) static const double vst_2x2v_ser_p1_vel_dv1_v0[8] = { 
  0.0, 
  0.0, 
  0.8660254037844386, 
  -1.5, 
  0.0, 
  0.0, 
  0.8660254037844386, 
  -1.5, 
}; 
__attribute__((unused)) static const double vst_2x2v_ser_p1_vel_ev_v1[8] = { 
  0.5, 
  -0.5, 
  -0.8660254037844386, 
  0.8660254037844386, 
  0.5, 
  0.5, 
  -0.8660254037844386, 
  -0.8660254037844386, 
}; 
__attribute__((unused)) static const double vst_2x2v_ser_p1_vel_dv0_v1[8] = { 
  0.0, 
  0.8660254037844386, 
  0.0, 
  -1.5, 
  0.0, 
  0.8660254037844386, 
  0.0, 
  -1.5, 
}; 
__attribute__((unused)) static const int vst_2x2v_ser_p1_vel_sparse_idx[3] = { 
  0, 
  1, 
  2, 
}; 
