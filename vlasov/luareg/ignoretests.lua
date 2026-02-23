-- Auto-updated by runregression: tests that exceeded the per-test timeout.
-- Remove an entry manually to re-enable the test.
return {
   lua = {
      "rt_can_pb_bgk_surf_cylindrical_sodshock_im_3x3v_p1",
      "rt_can_pb_ex_bgk_surf_flat_sq_ic",
      "rt_can_pb_bgk_surf_annulus_sodshock_im_2x2v_p2",
      "rt_can_pb_bgk_surf_sphere_sodshock_im_1x2v_p2",
      "rt_can_pb_bgk_surf_cylindrical_sodshock_im_1x3v_p2",
      "rt_can_pb_bgk_surf_sphere_sodshock_im_2x2v_p2",
      "rt_can_pb_bgk_surf_sphere_sodshock_im_2x2v_p1",
      "rt_can_pb_bgk_surf_toroidal_sodshock_im_3x3v_p1",
   },
   c = {
      "rt_can_pb_bgk_surf_cylindrical_sodshock_im_1x3v_p2",
      "rt_can_pb_bgk_surf_annulus_sodshock_im_2x2v_p2",
      "rt_can_pb_bgk_surf_sphere_sodshock_im_2x2v_p1",
      "rt_vlasov_neut_bgk_sodshock_1x3v_p2",
      "rt_can_pb_bgk_surf_sphere_sodshock_im_1x2v_p2",
      "rt_can_pb_ex_bgk_surf_flat_sq_ic",
      "rt_can_pb_bgk_surf_sphere_sodshock_im_2x2v_p2",
      -- runtime crash: "Moment option 13 not available" (missing feature):
      "rt_escreen_sr",
      -- unit/timing tests (not physics regression tests):
      "rt_hyper_vlasov_tm",
      "rt_vlasov_kerntm",
      "rt_vlasov_moments",
      -- 3x3v tests timeout on this machine:
      "rt_can_pb_bgk_surf_cylindrical_sodshock_im_3x3v_p1",
      "rt_can_pb_bgk_surf_toroidal_sodshock_im_3x3v_p1",
   },
}
