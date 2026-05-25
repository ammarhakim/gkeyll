-- 2D Bondi-Hoyle-Lyttleton accretion onto a static (Schwarzschild) black
-- hole, modular tetrad-basis GR Euler with the RYU-CHATTOPADHYAY equation
-- of state. Structural clone of rt_gr_bhl_static_tetrad_mod_tm.lua with
-- the EOS switched from Mathews-Taub to Ryu-Chattopadhyay (Ryu+ 2006).
--
-- EOS selection: explicit `eos = "rc"` dispatches the equation object to
-- the RC closure
--   h(θ) = 2(6θ² + 4θ + 1)/(3θ + 2),  θ = p/ρ
-- which has the same non-rel (Γ→5/3) and ultra-rel (Γ→4/3) asymptotic
-- limits as TM but fits the Synge (single-component perfect relativistic
-- gas) enthalpy tighter: ~0.8% maximum h error vs TM's ~2%. Useful for
-- BHL/shock-heated flows where the trans-rel post-shock thermodynamics
-- depends on the enthalpy curve away from the asymptotic limits.
--
-- Recovery cost: TM has a closed-form analytic cubic. RC requires Newton
-- iteration on a degree-8 polynomial (Ryu+ eq 29), warm-started from the
-- TM cubic for fast convergence (~3-5 iterations in practice). The cold-
-- flow fallback to EM Newton at γ=5/3 (θ < 1e-4) is shared with TM.
--
-- Roe is incompatible with non-IDEAL EOSs at the equation-object level
-- (SR-Roe eigenstructure uses the ideal-gas Jacobian). rpType="hll".
--
-- See rt_gr_bhl_static_tetrad_mod.lua for the geometric / spacetime setup.

local Moments = G0.Moments
local GREulerTetradMod = G0.Moments.Eq.GREulerTetradMod
local BlackHole = G0.Moments.Spacetime.BlackHole

rhol = 3.0
ul   = 0.3
pl   = 0.05

rhor = 0.01
ur   = 0.0
pr   = 0.01

mass = 0.3
spin = 0.0

pos_x = 2.5
pos_y = 2.5
pos_z = 0.0

Nx = 256
Ny = 256
Lx = 5.0
Ly = 5.0
cfl_frac = 0.5  -- Strict 2D directional-splitting CFL (cfla_x + cfla_y ≤ 1.0).

reinit_freq = 100

t_end = 15.0
num_frames = 30
field_energy_calcs = GKYL_MAX_INT
integrated_mom_calcs = GKYL_MAX_INT
dt_failure_tol = 1.0e-4
num_failures_max = 20

x_loc = 1.0

-- Ryu-Chattopadhyay specific enthalpy: h(θ) = 2(6θ² + 4θ + 1)/(3θ + 2),
-- θ = p/ρ. Same IC convention as TM — the C-side closure type is what
-- changes; this lua-side helper only needs to seed the conservative
-- state from the initial primitive state, which uses the same enthalpy
-- formula the equation object will compute internally.
local function rc_enthalpy(rho, p)
  local theta = p / rho
  return 2.0 * (6.0*theta*theta + 4.0*theta + 1.0) / (3.0*theta + 2.0)
end

momentApp = Moments.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { 0.0, 0.0 },
  upper = { Lx, Ly },
  cells = { Nx, Ny },
  cflFrac = cfl_frac,

  decompCuts = { 1, 1 },
  periodicDirs = { },

  spacetime = Moments.Spacetime.new {
    evolve = false,
    equation = BlackHole.new(mass, spin, pos_x, pos_y, pos_z),
    spacetimeGauge = G0.SpacetimeGauge.Static,
    reinitFreq = reinit_freq,
  },

  fluid = Moments.Species.new {
    -- Explicit EOS selection — "rc" dispatches to Ryu-Chattopadhyay.
    equation = GREulerTetradMod.new {
      eos = "rc",
      rpType = "hll",
    },

    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local rho = 0.0
      local u = 0.0
      local p = 0.0

      if x < x_loc then
        rho = rhol
        u = ul
        p = pl
      else
        rho = rhor
        u = ur
        p = pr
      end

      local spatial_metric = BlackHole.spatialMetricTensor(mass, spin, pos_x, pos_y, pos_z, 0.0, x, y, 0.0)
      local spatial_det    = BlackHole.spatialMetricDeterminant(mass, spin, pos_x, pos_y, pos_z, 0.0, x, y, 0.0)
      local in_excision_region = BlackHole.excisionRegion(mass, spin, pos_x, pos_y, pos_z, 0.0, x, y, 0.0)

      local vel = { u, 0.0, 0.0 }
      local v_sq = 0.0
      for i = 1, 3 do
        for j = 1, 3 do
          v_sq = v_sq + spatial_metric[i][j] * vel[i] * vel[j]
        end
      end

      local W = 1.0 / math.sqrt(1.0 - v_sq)
      if v_sq > 1.0 - math.pow(10.0, -8.0) then
        W = 1.0 / math.sqrt(1.0 - math.pow(10.0, -8.0))
      end

      -- Ryu-Chattopadhyay specific enthalpy at the IC primitive state.
      local h = rc_enthalpy(rho, p)

      -- Convention A: S_i = γ_ij·ρhW²·v^j (genuine covariant momentum).
      local v_lower = { 0.0, 0.0, 0.0 }
      for i = 1, 3 do
        for j = 1, 3 do
          v_lower[i] = v_lower[i] + spatial_metric[i][j] * vel[j]
        end
      end

      local rho_rel = math.sqrt(spatial_det) * rho * W
      local mom_x   = math.sqrt(spatial_det) * rho * h * (W * W) * v_lower[1]
      local mom_y   = math.sqrt(spatial_det) * rho * h * (W * W) * v_lower[2]
      local mom_z   = math.sqrt(spatial_det) * rho * h * (W * W) * v_lower[3]
      local Etot    = math.sqrt(spatial_det) * ((rho * h * (W * W)) - p - (rho * W))

      if in_excision_region then
        rho_rel, mom_x, mom_y, mom_z, Etot = 0.0, 0.0, 0.0, 0.0, 0.0
      end

      return rho_rel, mom_x, mom_y, mom_z, Etot
    end,

    evolve = true,
    forceLowOrderFlux = false,
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy },
    bcy = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy },
  }
}

momentApp:run()
