-- MOAT (Mother Of All Tests) for the moments Lua suite.
-- These tests run when --moat is passed and form the core acceptance set.
-- Format: bare test basenames (without .lua extension).
return {
   "rt_euler_sodshock",
   "rt_mhd_brio_wu",
   "rt_10m_riem",
   "rt_5m_gem",
   "rt_maxwell_plane_wave_1d",
   "rt_gr_einstein_plane_shock",
}
