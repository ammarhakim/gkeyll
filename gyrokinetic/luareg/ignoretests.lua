-- Auto-updated by runregression: tests that exceeded the per-test timeout.
-- Remove an entry manually to re-enable the test.
return {
   lua = {
      "rt_gk_sheath_1x2v_p1",
      "rt_gk_sheath_1x2v_p1_cons",
      "rt_gk_sheath_2x2v_p1",
      "rt_gk_sheath_3x2v_p1",
      "rt_gk_sheath_nonuniformx_1x2v_p1",
      "rt_gk_sheath_nonuniformx_2x2v_p1",
      "rt_gk_sheath_nonuniformx_3x2v_p1",
      "rt_gk_sheath_nonuniformv_1x2v_p1",
      "rt_gk_sheath_nonuniformv_2x2v_p1",
      "rt_gk_sheath_nonuniformv_3x2v_p1",
      "rt_gk_bgk_relax_1x2v_p1",
      "rt_gk_bgk_cross_relax_1x2v_p1",
      "rt_gk_lbo_cross_relax_1x2v_p1",
      "rt_gk_ar_react_nonuniformv_1x2v_p1",
      "rt_gk_rad_1x2v_p1",
      "rt_gk_rad_nonuniformv_1x2v_p1",
      "rt_gk_li_react_nonuniformv_3x2v_p1",
      "rt_gk_li_react_3x2v_p1",
   },
   c = {
      "rt_gk_mirror_boltz_elc_poa_1x2v_p1",
      "rt_gk_sheath_nonuniformx_3x2v_p1",
      "rt_gk_multib_step_sol_2x2v_p1",
      "rt_gk_wham_3x2v_p1",
      "rt_gk_multib_step_nonuniform_2x2v_p1",
      "rt_gk_wham_boltz_elc_poa_1x2v_p1",
      "rt_gk_multib_step_2x2v_p1",
      "rt_gk_wham_nonuniformx_1x2v_p1",
      "rt_gk_wham_1x2v_p1_static_field",
      "rt_gk_wham_nonuniformx_3x2v_p1",
      -- runtime segfault (3x2v test, 24k output files generated during create):
      "rt_gk_wham_2xIC_3x2v_p1",
   },
}
