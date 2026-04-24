-- MOAT (Mother Of All Tests) for the vlasov C suite.
-- These tests run when --moat is passed and form the core acceptance set.
-- Format: bare test basenames (without .c extension).
return {
   "rt_vlasov_landau_damping_1x1v_p2",
   "rt_vlasov_twostream_p2",
   "rt_vlasov_es_shock",
   "rt_vlasov_bgk_relax_1x1v_p2",
}
