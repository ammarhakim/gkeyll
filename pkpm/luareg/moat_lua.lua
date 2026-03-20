-- MOAT (Mother Of All Tests) for the pkpm Lua suite.
-- These tests run when --moat is passed and form the core acceptance set.
-- Format: bare test basenames (without .lua extension).
return {
   "rt_pkpm_landau_damping_p1",
   "rt_pkpm_em_advect_p1",
   "rt_pkpm_neut_sodshock_p1",
   "rt_pkpm_travel_pulse_p1",
}
