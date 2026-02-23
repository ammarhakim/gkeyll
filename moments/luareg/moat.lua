-- moat.lua for the moments layer
--
-- Lists the Main-On-Acceptance-Tests (MOAT) for the moments layer.
-- These are the key regression tests that should always pass before
-- accepting a change. Running with --moat runs only this subset.
--
-- Format: return a list of absolute paths to .lua test files.
--
-- Example:
--   return {
--     "/path/to/gkeyll/moments/luareg/rt_euler_sodshock.lua",
--     "/path/to/gkeyll/moments/luareg/rt_5m_gem.lua",
--   }

return {}
