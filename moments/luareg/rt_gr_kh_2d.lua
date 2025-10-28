local Moments = G0.Moments
local GREuler = G0.Moments.Eq.GREuler
local Minkowski = G0.Moments.Spacetime.Minkowski

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
gas_gamma = 4.0 / 3.0 -- Adiabatic index.

rho_u = 1.5 -- Upper fluid mass density.
u_u = 0.25 -- Upper fluid x-velocity.
v_u = 1.0 / 400.0 -- Upper fluid y-velocity.
p_u = 20.0 -- Upper fluid pressure.

rho_l = 0.5 -- Lower fluid mass density.
u_l = 0.25 -- Lower fluid x-velocity.
v_l = 1.0 / 400.0 -- Lower fluid y-velocity.
p_l = 20.0 -- Lower fluid pressure.

-- Simulation parameters.
Nx = 256 -- Cell count (x-direction).
Ny = 128 -- Cell count (y-direction).
Lx = 1.0 -- Domain size (x-direction).
Ly = 0.5 -- Domain size (y-direction).
cfl_frac = 0.95 -- CFL coefficient.

reinit_freq = GKYL_MAX_INT -- Spacetime reinitialization frequency.

t_end = 5.0 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

y_loc = 0.0 -- Fluid boundary (y-direction).

momentApp = Moments.App.new {
  
  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { 0.0, -0.5 * Ly },
  upper = { Lx, 0.5 * Ly },
  cells = { Nx, Ny },
  cflFrac = cfl_frac,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coodinate direction (x- and y-directions).
  
  -- Boundary conditions for configuration space.
  periodicDirs = { 1 }, -- Periodic directions (x-direction).

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
    
      if y > y_loc then
        rho = rho_u -- Fluid mass density (upper).
        u = u_u * math.tanh(100.0 * y) -- Fluid x-velocity (upper).
        v = v_u * math.sin(2.0 * pi * x) * math.exp(-100.0 * y * y) -- Fluid y-velocity (upper).
        p = p_u -- Fluid pressure (upper).
      else
        rho = rho_l -- Fluid mass density (lower).
        u = u_l * math.tanh(100.0 * y) -- Fluid x-velocity (lower).
        v = v_l * math.sin(2.0 * pi * x) * math.exp(-100.0 * y * y) -- Fluid y-velocity (lower).
        p = p_l -- Fluid pressure (lower).
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
    bcy = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy } -- Copy boundary conditions (y-direction).
  }
}

-- Run application.
momentApp:run()
