#pragma once
#include <math.h>
#include <gkyl_util.h>
EXTERN_C_BEG

GKYL_CU_DH double fpo_vlasov_drag_vol_1x3v_ser_p1(const double* w, const double* dx, const double* h, const double* f, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_vol_1x3v_ser_p1(const double* w, const double* dx, const double* g, const double* f, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_drag_surfvxvx_1x3v_ser_p1(const double* w, const double* dx, const double* hl, const double* hc, const double* hu, const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvxvx_1x3v_ser_p1(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_drag_surfvxvy_1x3v_ser_p1(const double* w, const double* dx, const double* hl, const double* hc, const double* hu, const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvxvy_1x3v_ser_p1(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_drag_surfvxvz_1x3v_ser_p1(const double* w, const double* dx, const double* hl, const double* hc, const double* hu, const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvxvz_1x3v_ser_p1(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_drag_surfvyvx_1x3v_ser_p1(const double* w, const double* dx, const double* hl, const double* hc, const double* hu, const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvyvx_1x3v_ser_p1(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_drag_surfvyvy_1x3v_ser_p1(const double* w, const double* dx, const double* hl, const double* hc, const double* hu, const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvyvy_1x3v_ser_p1(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_drag_surfvyvz_1x3v_ser_p1(const double* w, const double* dx, const double* hl, const double* hc, const double* hu, const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvyvz_1x3v_ser_p1(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_drag_surfvzvx_1x3v_ser_p1(const double* w, const double* dx, const double* hl, const double* hc, const double* hu, const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvzvx_1x3v_ser_p1(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_drag_surfvzvy_1x3v_ser_p1(const double* w, const double* dx, const double* hl, const double* hc, const double* hu, const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvzvy_1x3v_ser_p1(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_drag_surfvzvz_1x3v_ser_p1(const double* w, const double* dx, const double* hl, const double* hc, const double* hu, const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvzvz_1x3v_ser_p1(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_drag_vol_1x3v_ser_p2(const double* w, const double* dx, const double* h, const double* f, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_vol_1x3v_ser_p2(const double* w, const double* dx, const double* g, const double* f, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_drag_surfvxvx_1x3v_ser_p2(const double* w, const double* dx, const double* hl, const double* hc, const double* hu, const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvxvx_1x3v_ser_p2(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_drag_surfvxvy_1x3v_ser_p2(const double* w, const double* dx, const double* hl, const double* hc, const double* hu, const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvxvy_1x3v_ser_p2(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_drag_surfvxvz_1x3v_ser_p2(const double* w, const double* dx, const double* hl, const double* hc, const double* hu, const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvxvz_1x3v_ser_p2(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_drag_surfvyvx_1x3v_ser_p2(const double* w, const double* dx, const double* hl, const double* hc, const double* hu, const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvyvx_1x3v_ser_p2(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_drag_surfvyvy_1x3v_ser_p2(const double* w, const double* dx, const double* hl, const double* hc, const double* hu, const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvyvy_1x3v_ser_p2(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_drag_surfvyvz_1x3v_ser_p2(const double* w, const double* dx, const double* hl, const double* hc, const double* hu, const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvyvz_1x3v_ser_p2(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_drag_surfvzvx_1x3v_ser_p2(const double* w, const double* dx, const double* hl, const double* hc, const double* hu, const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvzvx_1x3v_ser_p2(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_drag_surfvzvy_1x3v_ser_p2(const double* w, const double* dx, const double* hl, const double* hc, const double* hu, const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvzvy_1x3v_ser_p2(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_drag_surfvzvz_1x3v_ser_p2(const double* w, const double* dx, const double* hl, const double* hc, const double* hu, const double* fl, const double* fc, const double* fu, double* GKYL_RESTRICT out);
GKYL_CU_DH double fpo_vlasov_diff_surfvzvz_1x3v_ser_p2(const double* w, const double* dx, const double* g[], const double* f[], double* GKYL_RESTRICT out);

EXTERN_C_END