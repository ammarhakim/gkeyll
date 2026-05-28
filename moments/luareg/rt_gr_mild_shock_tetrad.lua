-- Mildly relativistic blast wave test for the modular tetrad-basis GR
-- Euler equations. Structural mirror of rt_gr_mild_shock_tetrad.lua
-- (packed); the spacetime is supplied as a separate Moments.Spacetime
-- component instead of being bundled into the species conserved state.
-- Initial conditions from Del Zanna & Bucciantini 2002, A&A 390, 1177,
-- Section 4.1 (blast wave 1). https://arxiv.org/abs/astro-ph/0205290

local Moments = G0.Moments
local GREulerTetrad = G0.Moments.Eq.GREulerTetrad
local Minkowski = G0.Moments.Spacetime.Minkowski

gas_gamma = 5.0 / 3.0

rhol = 10.0
ul   = 0.0
pl   = 13.3

rhor = 1.0
ur   = 0.0
pr   = math.pow(10.0, -6.0)

Nx = 4096
Lx = 1.0
cfl_frac = 0.85

reinit_freq = 100

t_end = 0.4
num_frames = 1
field_energy_calcs = GKYL_MAX_INT
integrated_mom_calcs = GKYL_MAX_INT
dt_failure_tol = 1.0e-4
num_failures_max = 20

momentApp = Moments.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { 0.0 },
  upper = { Lx },
  cells = { Nx },
  cflFrac = cfl_frac,

  decompCuts = { 1 },
  periodicDirs = { },

  spacetime = Moments.Spacetime.new {
    evolve = false,
    equation = Minkowski.new(),
    spacetimeGauge = G0.SpacetimeGauge.Static,
    reinitFreq = reinit_freq,
  },

  fluid = Moments.Species.new {
    equation = GREulerTetrad.new {
      gasGamma = gas_gamma,
      rpType = "hll",
    },

    init = function (t, xn)
      local x = xn[1]

      local rho = 0.0
      local u = 0.0
      local p = 0.0

      if x < 0.5 then
        rho = rhol
        u = ul
        p = pl
      else
        rho = rhor
        u = ur
        p = pr
      end

      local spatial_metric = Minkowski.spatialMetricTensor(0.0, x, 0.0, 0.0)
      local spatial_det    = Minkowski.spatialMetricDeterminant(0.0, x, 0.0, 0.0)

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

      local h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)))

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

      return rho_rel, mom_x, mom_y, mom_z, Etot
    end,

    evolve = true,
    limiter = G0.WaveLimiter.MinMod,
    forceLowOrderFlux = false,
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy },
  }
}

momentApp:run()
