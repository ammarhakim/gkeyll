-- 2D quadrants test for the general relativistic Euler equations.
-- Input parameters taken from the initial conditions in Section 4.2 (Riemann 2-D), from the article:
-- L. Del Zanna and N. Bucciantini (2002), "An efficient shock-capturing central-type scheme for multdimensional flows. I. Hydrodynamics",
-- Astronomy and Astrophysics, Volume 390 (3): 1177-1186.
-- https://arxiv.org/abs/astro-ph/0205290

local Moments = G0.Moments
local GREuler = G0.Moments.Eq.GREuler
local Minkowski = G0.Moments.Spacetime.Minkowski

-- Physical constants (using normalized code units).
gas_gamma = 5.0 / 3.0 -- Adiabatic index.

rho_ul = 0.1 -- Upper-left fluid mass density.
u_ul = 0.99 -- Upper-left fluid x-velocity.
v_ul = 0.0 -- Upper-left fluid y-velocity.
p_ul = 1.0 -- Upper-left fluid pressure.

rho_ur = 0.1 -- Upper-right fluid mass density.
u_ur = 0.0 -- Upper-right fluid x-velocity.
v_ur = 0.0 -- Upper-right fluid y-velocity.
p_ur = 0.01 -- Upper-right fluid pressure.

rho_ll = 0.5 -- Lower-left fluid mass density.
u_ll = 0.0 -- Lower-left fluid x-velocity.
v_ll = 0.0 -- Lower-left fluid y-velocity.
p_ll = 1.0 -- Lower-left fluid pressure.

rho_lr = 0.1 -- Lower-right fluid mass density.
u_lr = 0.0 -- Lower-right fluid x-velocity.
v_lr = 0.99 -- Lower-right fluid y-velocity.
p_lr = 1.0 -- Lower-right fluid pressure.

-- Simulation parameters.
Nx = 256 -- Cell count (x-direction).
Ny = 256 -- Cell count (y-direction).
Lx = 1.0 -- Domain size (x-direction).
Ly = 1.0 -- Domain size (y-direction).
cfl_frac = 0.8 -- CFL coefficient.

reinit_freq = GKYL_MAX_INT -- Spacetime reinitialization frequency.

t_end = 0.4 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

loc = 0.5 -- Fluid boundaries (both x and y coordinates).

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

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coodinate direction (x- and y-directions).
  
  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).

  -- Fluid.
  fluid = Moments.Species.new {
    equation = GREuler.new {
      gasGamma = gas_gamma,
      reinitFreq = reinit_freq
    },

    hasGREuler = true,
    GREulerGasGamma = gas_gamma,
  
    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]
      
      local rho = 0.0
      local u = 0.0
      local v = 0.0
      local p = 0.0
    
      if y > loc then
        if x < loc then
          rho = rho_ul -- Fluid mass density (upper-left).
          u = u_ul -- Fluid x-velocity (upper-left).
          v = v_ul -- Fluid y-velocity (upper-left).
          p = p_ul -- Fluid pressure (upper-left).
        else
          rho = rho_ur -- Fluid mass density (upper-right).
          u = u_ur -- Fluid x-velocity (upper-right).
          v = v_ur -- Fluid y-velocity (upper-right).
          p = p_ur -- Fluid pressure (upper-right).
        end
      else
        if x < loc then
          rho = rho_ll -- Fluid mass density (lower-left).
          u = u_ll -- Fluid x-velocity (lower-left).
          v = v_ll -- Fluid y-velocity (lower-left).
          p = p_ll -- Fluid pressure (lower-left).
        else
          rho = rho_lr -- Fluid mass density (lower-right).
          u = u_lr -- Fluid x-velocity (lower-right).
          v = v_lr -- Fluid y-velocity (lower-right).
          p = p_lr -- Fluid pressure (lower-right).
        end
      end

      local lapse = Minkowski.lapseFunction(0.0, x, y, 0.0)
      local shift = Minkowski.shiftVector(0.0, x, y, 0.0)
      local spatial_metric = Minkowski.spatialMetricTensor(0.0, x, y, 0.0)
      local spatial_det = Minkowski.spatialMetricDeterminant(0.0, x, y, 0.0)
      local extrinsic_curvature = Minkowski.extrinsicCurvatureTensor(0.0, x, y, 0.0, 1.0, 1.0, 1.0)
      local in_excision_region = Minkowski.excisionRegion(0.0, x, y, 0.0)

      local lapse_der = Minkowski.lapseFunctionDer(0.0, x, y, 0.0, 1.0, 1.0, 1.0)
      local shift_der = Minkowski.shiftVectorDer(0.0, x, y, 0.0, 1.0, 1.0, 1.0)
      local spatial_metric_der = Minkowski.spatialMetricTensorDer(0.0, x, y, 0.0, 1.0, 1.0, 1.0)

      local vel = { u, v, 0.0 }
      local v_sq = 0.0

      for i = 1, 3 do
        for j = 1, 3 do
          v_sq = v_sq + spatial_metric[i][j] * vel[i] * vel[j]
        end
      end

      local W = 1.0 / math.sqrt(1.0 - v_sq)
      if v_sq > 1.0 - math.pow(10.0, -8.0) then
        W = 1.0 / math.sqrt(1.0 - pow(10.0, -8.0))
      end

      local h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)))

      local rho_rel = math.sqrt(spatial_det) * rho * W -- Fluid relativistic mass density.
      local mom_x = math.sqrt(spatial_det) * rho * h * (W * W) * u -- Fluid momentum density (x-direction).
      local mom_y = math.sqrt(spatial_det) * rho * h * (W * W) * v -- Fluid momentum density (y-direction).
      local mom_z = 0.0 -- Fluid momentum density (z-direction).
      local Etot = math.sqrt(spatial_det) * ((rho * h * (W * W)) - p - (rho * W)) -- Fluid total energy density.

      local excision = 0.0
      if in_excision_region then
        rho_rel, mom_x, mom_y, mom_z, Etot, lapse = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        for i = 1, 3 do
          shift[i] = 0.0
          lapse_der[i] = 0.0
          for j = 1, 3 do
            spatial_metric[i][j] = 0.0
            extrinsic_curvature[i][j] = 0.0
            shift_der[i][j] = 0.0
            for k = 1, 3 do
              spatial_metric_der[i][j][k] = 0.0
            end
          end  
        end
        
        excision = -1.0
      else
        excision = 1.0
      end
    
      return rho_rel, mom_x, mom_y, mom_z, Etot,
        lapse,
        shift[1], shift[2], shift[3],
        spatial_metric[1][1], spatial_metric[1][2], spatial_metric[1][3],
        spatial_metric[2][1], spatial_metric[2][2], spatial_metric[2][3],
        spatial_metric[3][1], spatial_metric[3][2], spatial_metric[3][3],
        extrinsic_curvature[1][1], extrinsic_curvature[1][2], extrinsic_curvature[1][3],
        extrinsic_curvature[2][1], extrinsic_curvature[2][2], extrinsic_curvature[2][3],
        extrinsic_curvature[3][1], extrinsic_curvature[3][2], extrinsic_curvature[3][3],
        excision,
        lapse_der[1], lapse_der[2], lapse_der[3],
        shift_der[1][1], shift_der[1][2], shift_der[1][3],
        shift_der[2][1], shift_der[2][2], shift_der[2][3],
        shift_der[3][1], shift_der[3][2], shift_der[3][3],
        spatial_metric_der[1][1][1], spatial_metric_der[1][1][2], spatial_metric_der[1][1][3],
        spatial_metric_der[1][2][1], spatial_metric_der[1][2][2], spatial_metric_der[1][2][3],
        spatial_metric_der[1][3][1], spatial_metric_der[1][3][2], spatial_metric_der[1][3][3],
        spatial_metric_der[2][1][1], spatial_metric_der[2][1][2], spatial_metric_der[2][1][3],
        spatial_metric_der[2][2][1], spatial_metric_der[2][2][2], spatial_metric_der[2][2][3],
        spatial_metric_der[2][3][1], spatial_metric_der[2][3][2], spatial_metric_der[2][3][3],
        spatial_metric_der[3][1][1], spatial_metric_der[3][1][2], spatial_metric_der[3][1][3],
        spatial_metric_der[3][2][1], spatial_metric_der[3][2][2], spatial_metric_der[3][2][3],
        spatial_metric_der[3][3][1], spatial_metric_der[3][3][2], spatial_metric_der[3][3][3],
        0.0,
        x, y, 0.0
    end,

    evolve = true, -- Evolve species?
    limiter = G0.WaveLimiter.MinMod,
    forceLowOrderFlux = false, -- Use HLL fluxes.
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy }, -- Copy boundary conditions (x-direction).
    bcy = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy } -- Copy boundary conditions (y-direction).
  }
}

-- Run application.
momentApp:run()
