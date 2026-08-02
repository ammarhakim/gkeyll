#pragma once 
// Generated nodal surface-evaluation tables (see surf-eval-tables.mac). 
// Header-static so both host and device instantiations of GKYL_CU_DH 
// kernels can read them; unused arrays are elided per TU. 
#include <gkyl_util.h> 
#include <gkyl_vlasov_surf_meta.h> 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_conf_ev[16] = { 
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
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_conf_dx0[16] = { 
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
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_vel_ev_v0[2] = { 
  0.7071067811865475, 
  -1.224744871391589, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_ser_p1_vel_sparse_idx[2] = { 
  0, 
  1, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_v0_Cm[16] = { 
  1.0, 
  -1.4433756729740645, 
  -1.4433756729740645, 
  2.0833333333333335, 
  1.0, 
  -1.4433756729740645, 
  1.4433756729740645, 
  -2.0833333333333335, 
  1.0, 
  1.4433756729740645, 
  -1.4433756729740645, 
  -2.0833333333333335, 
  1.0, 
  1.4433756729740645, 
  1.4433756729740645, 
  2.0833333333333335, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_v0_CmD[16] = { 
  0.0, 
  2.5, 
  0.0, 
  -3.6084391824351614, 
  0.0, 
  2.5, 
  0.0, 
  3.6084391824351614, 
  0.0, 
  2.5, 
  0.0, 
  -3.6084391824351614, 
  0.0, 
  2.5, 
  0.0, 
  3.6084391824351614, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_v0_CmDx0[16] = { 
  0.0, 
  2.5, 
  0.0, 
  -3.6084391824351614, 
  0.0, 
  2.5, 
  0.0, 
  3.6084391824351614, 
  0.0, 
  2.5, 
  0.0, 
  -3.6084391824351614, 
  0.0, 
  2.5, 
  0.0, 
  3.6084391824351614, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_v0_CmDx1[16] = { 
  0.0, 
  0.0, 
  2.5, 
  -3.6084391824351614, 
  0.0, 
  0.0, 
  2.5, 
  -3.6084391824351614, 
  0.0, 
  0.0, 
  2.5, 
  3.6084391824351614, 
  0.0, 
  0.0, 
  2.5, 
  3.6084391824351614, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_v0_V[2] = { 
  1.0, 
  2.5, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_v0_coefl[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.11757550765359256, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_v0_coefr[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  -0.2449489742783178, 
  0.16970562748477142, 
  -0.16970562748477142, 
  -0.16970562748477142, 
  -0.11757550765359256, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_ser_p1_ph_v0_cmap[8] = { 
  0, 
  1, 
  2, 
  0, 
  3, 
  1, 
  2, 
  3, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_ser_p1_ph_v0_vlmap[8] = { 
  0, 
  0, 
  0, 
  1, 
  0, 
  1, 
  1, 
  1, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_ser_p1_ph_v0_vrmap[8] = { 
  0, 
  0, 
  0, 
  1, 
  0, 
  1, 
  1, 
  1, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_v0_Vd0[2] = { 
  0.0, 
  2.5, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_ser_p1_ph_v0_vld0map[8] = { 
  0, 
  0, 
  0, 
  1, 
  0, 
  1, 
  1, 
  1, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_ser_p1_ph_v0_vrd0map[8] = { 
  0, 
  0, 
  0, 
  1, 
  0, 
  1, 
  1, 
  1, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_v0_dcoefl0[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.11757550765359256, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_v0_dcoefr0[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.11757550765359256, 
}; 
GKYL_CU_D __attribute__((unused)) static const struct gkyl_vlasov_surf_meta vst_2x1v_ser_p1_ph_v0_meta = { .nb = 8, .na = 4, .nmv = 2, .no = 4, .ni = 1, .ns = 2, .fhat_off = 0 }; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_prj_v0_Cw[16] = { 
  1.0, 
  -1.4433756729740645, 
  -1.4433756729740645, 
  2.0833333333333335, 
  1.0, 
  -1.4433756729740645, 
  1.4433756729740645, 
  -2.0833333333333335, 
  1.0, 
  1.4433756729740645, 
  -1.4433756729740645, 
  -2.0833333333333335, 
  1.0, 
  1.4433756729740645, 
  1.4433756729740645, 
  2.0833333333333335, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_prj_v0_Vw[1] = { 
  1.0, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_ser_p1_prj_v0_kamap[4] = { 
  0, 
  1, 
  2, 
  3, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_ser_p1_prj_v0_kbmap[4] = { 
  0, 
  0, 
  0, 
  0, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_ser_p1_prj_v0_out_off[5] = { 
  0, 
  2, 
  4, 
  6, 
  8, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_ser_p1_prj_v0_out_mode[8] = { 
  0, 
  3, 
  1, 
  5, 
  2, 
  6, 
  4, 
  7, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_prj_v0_out_cl[8] = { 
  0.3535533905932737, 
  -0.6123724356957944, 
  0.2449489742783178, 
  -0.42426406871192845, 
  0.2449489742783178, 
  -0.42426406871192845, 
  0.16970562748477142, 
  -0.29393876913398137, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_prj_v0_out_cr[8] = { 
  -0.3535533905932737, 
  -0.6123724356957944, 
  -0.2449489742783178, 
  -0.42426406871192845, 
  -0.2449489742783178, 
  -0.42426406871192845, 
  -0.16970562748477142, 
  -0.29393876913398137, 
}; 
GKYL_CU_D __attribute__((unused)) static const struct gkyl_vlasov_prj_meta vst_2x1v_ser_p1_prj_v0_meta = { .nk = 4, .np = 8, .na = 4, .nmv = 1, .no = 4, .ni = 1, .fhat_off = 0 }; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_vel_vol_dv0[4] = { 
  0.0, 
  1.224744871391589, 
  0.0, 
  1.224744871391589, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_confsurf_x0_ev_l[8] = { 
  0.5, 
  0.8660254037844386, 
  -0.5, 
  -0.8660254037844386, 
  0.5, 
  0.8660254037844386, 
  0.5, 
  0.8660254037844386, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_confsurf_x0_ev_r[8] = { 
  0.5, 
  -0.8660254037844386, 
  -0.5, 
  0.8660254037844386, 
  0.5, 
  -0.8660254037844386, 
  0.5, 
  -0.8660254037844386, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_x0_Cm[4] = { 
  1.0, 
  -1.4433756729740645, 
  1.0, 
  1.4433756729740645, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_x0_V[8] = { 
  1.0, 
  2.5, 
  -1.4433756729740645, 
  -3.6084391824351614, 
  1.0, 
  2.5, 
  1.4433756729740645, 
  3.6084391824351614, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_x0_coefl[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.11757550765359256, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_x0_coefr[8] = { 
  0.3535533905932737, 
  -0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  -0.16970562748477142, 
  -0.16970562748477142, 
  0.16970562748477142, 
  -0.11757550765359256, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_ser_p1_ph_x0_cmap[8] = { 
  0, 
  0, 
  1, 
  0, 
  1, 
  0, 
  1, 
  1, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_ser_p1_ph_x0_vlmap[8] = { 
  0, 
  1, 
  0, 
  2, 
  1, 
  3, 
  2, 
  3, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_ser_p1_ph_x0_vrmap[8] = { 
  0, 
  1, 
  0, 
  2, 
  1, 
  3, 
  2, 
  3, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_x0_Vd0[6] = { 
  0.0, 
  2.5, 
  6.25, 
  0.0, 
  2.5, 
  6.25, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_ser_p1_ph_x0_vld0map[8] = { 
  0, 
  0, 
  0, 
  1, 
  0, 
  2, 
  1, 
  2, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_ser_p1_ph_x0_vrd0map[8] = { 
  0, 
  0, 
  0, 
  1, 
  0, 
  2, 
  1, 
  2, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_x0_dcoefl0[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.11757550765359256, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_x0_dcoefr0[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  -0.16970562748477142, 
  0.16970562748477142, 
  -0.11757550765359256, 
}; 
GKYL_CU_D __attribute__((unused)) static const struct gkyl_vlasov_surf_meta vst_2x1v_ser_p1_ph_x0_meta = { .nb = 8, .na = 2, .nmv = 4, .no = 2, .ni = 2, .ns = 2, .fhat_off = 0 }; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_prj_x0_Cw[4] = { 
  1.0, 
  -1.4433756729740645, 
  1.0, 
  1.4433756729740645, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_prj_x0_Vw[4] = { 
  1.0, 
  -1.4433756729740645, 
  1.0, 
  1.4433756729740645, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_ser_p1_prj_x0_kamap[4] = { 
  0, 
  1, 
  0, 
  1, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_ser_p1_prj_x0_kbmap[4] = { 
  0, 
  0, 
  1, 
  1, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_ser_p1_prj_x0_out_off[5] = { 
  0, 
  2, 
  4, 
  6, 
  8, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_ser_p1_prj_x0_out_mode[8] = { 
  0, 
  1, 
  2, 
  4, 
  3, 
  5, 
  6, 
  7, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_prj_x0_out_cl[8] = { 
  0.3535533905932737, 
  -0.6123724356957944, 
  0.2449489742783178, 
  -0.42426406871192845, 
  0.2449489742783178, 
  -0.42426406871192845, 
  0.16970562748477142, 
  -0.29393876913398137, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_prj_x0_out_cr[8] = { 
  -0.3535533905932737, 
  -0.6123724356957944, 
  -0.2449489742783178, 
  -0.42426406871192845, 
  -0.2449489742783178, 
  -0.42426406871192845, 
  -0.16970562748477142, 
  -0.29393876913398137, 
}; 
GKYL_CU_D __attribute__((unused)) static const struct gkyl_vlasov_prj_meta vst_2x1v_ser_p1_prj_x0_meta = { .nk = 4, .np = 8, .na = 2, .nmv = 2, .no = 2, .ni = 2, .fhat_off = 0 }; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_confsurf_x1_ev_l[8] = { 
  0.5, 
  -0.5, 
  0.8660254037844386, 
  -0.8660254037844386, 
  0.5, 
  0.5, 
  0.8660254037844386, 
  0.8660254037844386, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_confsurf_x1_ev_r[8] = { 
  0.5, 
  -0.5, 
  -0.8660254037844386, 
  0.8660254037844386, 
  0.5, 
  0.5, 
  -0.8660254037844386, 
  -0.8660254037844386, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_x1_Cm[4] = { 
  1.0, 
  -1.4433756729740645, 
  1.0, 
  1.4433756729740645, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_x1_V[8] = { 
  1.0, 
  2.5, 
  -1.4433756729740645, 
  -3.6084391824351614, 
  1.0, 
  2.5, 
  1.4433756729740645, 
  3.6084391824351614, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_x1_coefl[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.11757550765359256, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_x1_coefr[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  -0.2449489742783178, 
  0.2449489742783178, 
  -0.16970562748477142, 
  0.16970562748477142, 
  -0.16970562748477142, 
  -0.11757550765359256, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_ser_p1_ph_x1_cmap[8] = { 
  0, 
  1, 
  0, 
  0, 
  1, 
  1, 
  0, 
  1, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_ser_p1_ph_x1_vlmap[8] = { 
  0, 
  0, 
  1, 
  2, 
  1, 
  2, 
  3, 
  3, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_ser_p1_ph_x1_vrmap[8] = { 
  0, 
  0, 
  1, 
  2, 
  1, 
  2, 
  3, 
  3, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_x1_Vd0[6] = { 
  0.0, 
  2.5, 
  6.25, 
  0.0, 
  2.5, 
  6.25, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_ser_p1_ph_x1_vld0map[8] = { 
  0, 
  0, 
  0, 
  1, 
  0, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_ser_p1_ph_x1_vrd0map[8] = { 
  0, 
  0, 
  0, 
  1, 
  0, 
  1, 
  2, 
  2, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_x1_dcoefl0[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.16970562748477142, 
  0.11757550765359256, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_ph_x1_dcoefr0[8] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.2449489742783178, 
  0.16970562748477142, 
  0.16970562748477142, 
  -0.16970562748477142, 
  -0.11757550765359256, 
}; 
GKYL_CU_D __attribute__((unused)) static const struct gkyl_vlasov_surf_meta vst_2x1v_ser_p1_ph_x1_meta = { .nb = 8, .na = 2, .nmv = 4, .no = 2, .ni = 2, .ns = 2, .fhat_off = 4 }; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_prj_x1_Cw[4] = { 
  1.0, 
  -1.4433756729740645, 
  1.0, 
  1.4433756729740645, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_prj_x1_Vw[4] = { 
  1.0, 
  -1.4433756729740645, 
  1.0, 
  1.4433756729740645, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_ser_p1_prj_x1_kamap[4] = { 
  0, 
  1, 
  0, 
  1, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_ser_p1_prj_x1_kbmap[4] = { 
  0, 
  0, 
  1, 
  1, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_ser_p1_prj_x1_out_off[5] = { 
  0, 
  2, 
  4, 
  6, 
  8, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_ser_p1_prj_x1_out_mode[8] = { 
  0, 
  2, 
  1, 
  4, 
  3, 
  6, 
  5, 
  7, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_prj_x1_out_cl[8] = { 
  0.3535533905932737, 
  -0.6123724356957944, 
  0.2449489742783178, 
  -0.42426406871192845, 
  0.2449489742783178, 
  -0.42426406871192845, 
  0.16970562748477142, 
  -0.29393876913398137, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_ser_p1_prj_x1_out_cr[8] = { 
  -0.3535533905932737, 
  -0.6123724356957944, 
  -0.2449489742783178, 
  -0.42426406871192845, 
  -0.2449489742783178, 
  -0.42426406871192845, 
  -0.16970562748477142, 
  -0.29393876913398137, 
}; 
GKYL_CU_D __attribute__((unused)) static const struct gkyl_vlasov_prj_meta vst_2x1v_ser_p1_prj_x1_meta = { .nk = 4, .np = 8, .na = 2, .nmv = 2, .no = 2, .ni = 2, .fhat_off = 4 }; 
