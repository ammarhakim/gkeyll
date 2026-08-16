#pragma once 
// Generated nodal surface-evaluation (MODAL->NODAL) tables for the 
// alpha_quad producers and lax_flux_nodal kernels (see 
// surf-eval-tables.mac). Header-static so both host and device 
// instantiations of GKYL_CU_DH kernels can read them; unused arrays 
// are elided per TU. 
#include <gkyl_util.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_meta.h> 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p1_conf_ev[4] = { 
  0.7071067811865475, 
  -0.7071067811865475, 
  0.7071067811865475, 
  0.7071067811865475, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p1_conf_dx0[4] = { 
  0.0, 
  1.224744871391589, 
  0.0, 
  1.224744871391589, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p1_vel_ev_v0[3] = { 
  0.7071067811865475, 
  -1.224744871391589, 
  1.5811388300841895, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_1x1v_tensor_p1_vel_sparse_idx[3] = { 
  0, 
  1, 
  2, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p1_ph_v0_Cm[4] = { 
  1.0, 
  -1.4433756729740645, 
  1.0, 
  1.4433756729740645, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p1_ph_v0_CmD[4] = { 
  0.0, 
  2.5, 
  0.0, 
  2.5, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p1_ph_v0_CmDx0[4] = { 
  0.0, 
  2.5, 
  0.0, 
  2.5, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p1_ph_v0_V[3] = { 
  1.0, 
  2.5, 
  -3.8461538461538463, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p1_ph_v0_coefl[6] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  -0.20139513400278566, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p1_ph_v0_coefr[6] = { 
  0.5, 
  0.34641016151377546, 
  -0.34641016151377546, 
  -0.24, 
  -0.2906888370749726, 
  -0.20139513400278566, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_1x1v_tensor_p1_ph_v0_cmap[6] = { 
  0, 
  1, 
  0, 
  1, 
  0, 
  1, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_1x1v_tensor_p1_ph_v0_vlmap[6] = { 
  0, 
  0, 
  1, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_1x1v_tensor_p1_ph_v0_vrmap[6] = { 
  0, 
  0, 
  1, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p1_ph_v0_Vd0[3] = { 
  0.0, 
  2.5, 
  -11.538461538461538, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_1x1v_tensor_p1_ph_v0_vld0map[6] = { 
  0, 
  0, 
  1, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_1x1v_tensor_p1_ph_v0_vrd0map[6] = { 
  0, 
  0, 
  1, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p1_ph_v0_dcoefl0[6] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  -0.20139513400278566, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p1_ph_v0_dcoefr0[6] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  0.2906888370749726, 
  0.20139513400278566, 
}; 
GKYL_CU_D __attribute__((unused)) static const struct gkyl_vlasov_surf_meta vst_1x1v_tensor_p1_ph_v0_meta = { .nb = 6, .na = 2, .nmv = 3, .no = 2, .ni = 1, .ns = 3, .fhat_off = 0 }; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_conf_ev[4] = { 
  0.7071067811865475, 
  -0.7071067811865475, 
  0.7071067811865475, 
  0.7071067811865475, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_conf_dx0[4] = { 
  0.0, 
  1.224744871391589, 
  0.0, 
  1.224744871391589, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_vel_ev_v0[3] = { 
  0.7071067811865475, 
  -1.224744871391589, 
  1.5811388300841895, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_1x1v_tensor_p1_ho_vel_sparse_idx[3] = { 
  0, 
  1, 
  2, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_ph_v0_Cm[4] = { 
  1.0, 
  -1.4433756729740645, 
  1.0, 
  1.4433756729740645, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_ph_v0_CmD[4] = { 
  0.0, 
  2.5, 
  0.0, 
  2.5, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_ph_v0_CmDx0[4] = { 
  0.0, 
  2.5, 
  0.0, 
  2.5, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_ph_v0_V[3] = { 
  1.0, 
  2.5, 
  -3.8461538461538463, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_ph_v0_coefl[6] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  -0.20139513400278566, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_ph_v0_coefr[6] = { 
  0.5, 
  0.34641016151377546, 
  -0.34641016151377546, 
  -0.24, 
  -0.2906888370749726, 
  -0.20139513400278566, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_1x1v_tensor_p1_ho_ph_v0_cmap[6] = { 
  0, 
  1, 
  0, 
  1, 
  0, 
  1, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_1x1v_tensor_p1_ho_ph_v0_vlmap[6] = { 
  0, 
  0, 
  1, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_1x1v_tensor_p1_ho_ph_v0_vrmap[6] = { 
  0, 
  0, 
  1, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_ph_v0_Vd0[3] = { 
  0.0, 
  2.5, 
  -11.538461538461538, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_1x1v_tensor_p1_ho_ph_v0_vld0map[6] = { 
  0, 
  0, 
  1, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_1x1v_tensor_p1_ho_ph_v0_vrd0map[6] = { 
  0, 
  0, 
  1, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_ph_v0_dcoefl0[6] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  -0.2906888370749726, 
  -0.20139513400278566, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_1x1v_tensor_p1_ho_ph_v0_dcoefr0[6] = { 
  0.5, 
  0.34641016151377546, 
  0.34641016151377546, 
  0.24, 
  0.2906888370749726, 
  0.20139513400278566, 
}; 
GKYL_CU_D __attribute__((unused)) static const struct gkyl_vlasov_surf_meta vst_1x1v_tensor_p1_ho_ph_v0_meta = { .nb = 6, .na = 2, .nmv = 3, .no = 2, .ni = 1, .ns = 3, .fhat_off = 0 }; 
