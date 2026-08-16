#pragma once 
// Generated projection (NODAL->MODAL) tables for the from-flux 
// vlasov_(boundary_)(ho_)surf* consumer kernels (see 
// surf-eval-tables.mac). Header-static so both host and device 
// instantiations of GKYL_CU_DH kernels can read them; unused arrays 
// are elided per TU. 
#include <gkyl_util.h> 
#include <gkyl_vlasov_surf_nod2mod_meta.h> 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_tensor_p1_prj_v0_Cw[16] = { 
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
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_tensor_p1_prj_v0_Vw[1] = { 
  1.0, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_tensor_p1_prj_v0_kamap[4] = { 
  0, 
  1, 
  2, 
  3, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_tensor_p1_prj_v0_kbmap[4] = { 
  0, 
  0, 
  0, 
  0, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_tensor_p1_prj_v0_out_off[5] = { 
  0, 
  3, 
  6, 
  9, 
  12, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_tensor_p1_prj_v0_out_mode[12] = { 
  0, 
  3, 
  7, 
  1, 
  5, 
  9, 
  2, 
  6, 
  10, 
  4, 
  8, 
  11, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_tensor_p1_prj_v0_out_cl[12] = { 
  0.3535533905932737, 
  -0.6123724356957944, 
  0.7905694150420947, 
  0.2449489742783178, 
  -0.42426406871192845, 
  0.5477225575051661, 
  0.2449489742783178, 
  -0.42426406871192845, 
  0.5477225575051661, 
  0.16970562748477142, 
  -0.29393876913398137, 
  0.3794733192202055, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_tensor_p1_prj_v0_out_cr[12] = { 
  -0.3535533905932737, 
  -0.6123724356957944, 
  -0.7905694150420947, 
  -0.2449489742783178, 
  -0.42426406871192845, 
  -0.5477225575051661, 
  -0.2449489742783178, 
  -0.42426406871192845, 
  -0.5477225575051661, 
  -0.16970562748477142, 
  -0.29393876913398137, 
  -0.3794733192202055, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_tensor_p1_prj_v0_pm_a[12] = { 
  0, 
  1, 
  2, 
  0, 
  3, 
  1, 
  2, 
  0, 
  3, 
  1, 
  2, 
  3, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_tensor_p1_prj_v0_pm_b[12] = { 
  0, 
  0, 
  0, 
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
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_tensor_p1_prj_v0_pm_cl[12] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  -0.6123724356957944, 
  0.16970562748477142, 
  -0.42426406871192845, 
  -0.42426406871192845, 
  0.7905694150420947, 
  -0.29393876913398137, 
  0.5477225575051661, 
  0.5477225575051661, 
  0.3794733192202055, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_tensor_p1_prj_v0_pm_cr[12] = { 
  -0.3535533905932737, 
  -0.2449489742783178, 
  -0.2449489742783178, 
  -0.6123724356957944, 
  -0.16970562748477142, 
  -0.42426406871192845, 
  -0.42426406871192845, 
  -0.7905694150420947, 
  -0.29393876913398137, 
  -0.5477225575051661, 
  -0.5477225575051661, 
  -0.3794733192202055, 
}; 
GKYL_CU_D __attribute__((unused)) static const struct gkyl_vlasov_prj_meta vst_2x1v_tensor_p1_prj_v0_meta = { .nk = 4, .np = 12, .na = 4, .nmv = 1, .no = 4, .ni = 1, .fhat_off = 0 }; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_tensor_p1_ho_prj_v0_Cw[16] = { 
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
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_tensor_p1_ho_prj_v0_Vw[1] = { 
  1.0, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_tensor_p1_ho_prj_v0_kamap[4] = { 
  0, 
  1, 
  2, 
  3, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_tensor_p1_ho_prj_v0_kbmap[4] = { 
  0, 
  0, 
  0, 
  0, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_tensor_p1_ho_prj_v0_out_off[5] = { 
  0, 
  3, 
  6, 
  9, 
  12, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_tensor_p1_ho_prj_v0_out_mode[12] = { 
  0, 
  3, 
  7, 
  1, 
  5, 
  9, 
  2, 
  6, 
  10, 
  4, 
  8, 
  11, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_tensor_p1_ho_prj_v0_out_cl[12] = { 
  0.3535533905932737, 
  -0.6123724356957944, 
  0.7905694150420947, 
  0.2449489742783178, 
  -0.42426406871192845, 
  0.5477225575051661, 
  0.2449489742783178, 
  -0.42426406871192845, 
  0.5477225575051661, 
  0.16970562748477142, 
  -0.29393876913398137, 
  0.3794733192202055, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_tensor_p1_ho_prj_v0_out_cr[12] = { 
  -0.3535533905932737, 
  -0.6123724356957944, 
  -0.7905694150420947, 
  -0.2449489742783178, 
  -0.42426406871192845, 
  -0.5477225575051661, 
  -0.2449489742783178, 
  -0.42426406871192845, 
  -0.5477225575051661, 
  -0.16970562748477142, 
  -0.29393876913398137, 
  -0.3794733192202055, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_tensor_p1_ho_prj_v0_pm_a[12] = { 
  0, 
  1, 
  2, 
  0, 
  3, 
  1, 
  2, 
  0, 
  3, 
  1, 
  2, 
  3, 
}; 
GKYL_CU_D __attribute__((unused)) static const int vst_2x1v_tensor_p1_ho_prj_v0_pm_b[12] = { 
  0, 
  0, 
  0, 
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
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_tensor_p1_ho_prj_v0_pm_cl[12] = { 
  0.3535533905932737, 
  0.2449489742783178, 
  0.2449489742783178, 
  -0.6123724356957944, 
  0.16970562748477142, 
  -0.42426406871192845, 
  -0.42426406871192845, 
  0.7905694150420947, 
  -0.29393876913398137, 
  0.5477225575051661, 
  0.5477225575051661, 
  0.3794733192202055, 
}; 
GKYL_CU_D __attribute__((unused)) static const double vst_2x1v_tensor_p1_ho_prj_v0_pm_cr[12] = { 
  -0.3535533905932737, 
  -0.2449489742783178, 
  -0.2449489742783178, 
  -0.6123724356957944, 
  -0.16970562748477142, 
  -0.42426406871192845, 
  -0.42426406871192845, 
  -0.7905694150420947, 
  -0.29393876913398137, 
  -0.5477225575051661, 
  -0.5477225575051661, 
  -0.3794733192202055, 
}; 
GKYL_CU_D __attribute__((unused)) static const struct gkyl_vlasov_prj_meta vst_2x1v_tensor_p1_ho_prj_v0_meta = { .nk = 4, .np = 12, .na = 4, .nmv = 1, .no = 4, .ni = 1, .fhat_off = 0 }; 
