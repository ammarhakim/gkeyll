local Moments = G0.Moments
local GRMHDTetrad = G0.Moments.Eq.GRMHDTetrad
local Minkowski = G0.Moments.Spacetime.Minkowski

-- Physical constants (using normalized code units).
gas_gamma = 5.0 / 3.0 -- Adiabatic index.

rhol = 1.0 -- Left fluid mass density.
ul = 0.0 -- Left fluid velocity.
pl = 30.0 -- Left fluid pressure.

Bxl = 5.0 -- Left fluid magnetic field (x-direction).
Byl = 6.0 -- Left fluid magnetic field (y-direction).
Bzl = 6.0 -- Left fluid magnetic field (z-direction).

rhor = 1.0 -- Right fluid mass density.
ur = 0.0 -- Right fluid velocity.
pr = 1.0 -- Right fluid pressure.

Bxr = 5.0 -- Right fluid magnetic field (x-direction).
Byr = 0.7 -- Right fluid magnetic field (y-direction).
Bzr = 0.7 -- Right fluid magnetic field (z-direction).

light_speed = 1.0 -- Speed of light.
b_fact = 0.8 -- Factor of speed of light for magnetic field correction.

-- Simulation parameters.
Nx = 4096 -- Cell count (x-direction).
Lx = 1.0 -- Domain size (x-direction).
cfl_frac = 0.85 -- CFL coefficient.

spacetime_gauge = G0.SpacetimeGauge.Static -- Spacetime gauge choice.
reinit_freq = 10 -- Spacetime reinitialization frequency.

t_end = 0.4 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

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

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).
  
  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).

  -- Fluid.
  fluid = Moments.Species.new {
    equation = GRMHDTetrad.new {
      gasGamma = gas_gamma,
      lightSpeed = light_speed,
      mgnErrorSpeedFactor = b_fact,
      spacetimeGauge = spacetime_gauge,
      reinitFreq = reinit_freq
    },

    hasGRMHD = true,
    GRMHDGasGamma = gas_gamma,
  
    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]
      
      local rho = 0.0
      local u = 0.0
      local p = 0.0

      local Bx = 0.0
      local By = 0.0
      local Bz = 0.0

      if x < 0.5 then
        rho = rhol -- Fluid mass density (left).
        u = ul -- Fluid velocity (left).
        p = pl -- Fluid pressure (left).

        Bx = Bxl -- Fluid magnetic field (x-direction, left).
        By = Byl -- Fluid magnetic field (y-direction, left).
        Bz = Bzl -- Fluid magnetic field (z-direction, left).
      else
        rho = rhor -- Fluid mass density (right).
        u = ur -- Fluid velocity (right).
        p = pr -- Fluid pressure (right).

        Bx = Bxr -- Fluid magnetic field (x-direction, right).
        By = Byr -- Fluid magnetic field (y-direction, right).
        Bz = Bzr -- Fluid magnetic field (z-direction, right).
      end

      local lapse = Minkowski.lapseFunction(0.0, x, 0.0, 0.0)
      local shift = Minkowski.shiftVector(0.0, x, 0.0, 0.0)
      local spatial_metric = Minkowski.spatialMetricTensor(0.0, x, 0.0, 0.0)
      local spatial_det = Minkowski.spatialMetricDeterminant(0.0, x, 0.0, 0.0)
      local extrinsic_curvature = Minkowski.extrinsicCurvatureTensor(0.0, x, 0.0, 0.0, 1.0, 1.0, 1.0)
      local in_excision_region = Minkowski.excisionRegion(0.0, x, 0.0, 0.0)

      local lapse_der = Minkowski.lapseFunctionDer(0.0, x, 0.0, 0.0, 1.0, 1.0, 1.0)
      local shift_der = Minkowski.shiftVectorDer(0.0, x, 0.0, 0.0, 1.0, 1.0, 1.0)
      local spatial_metric_der = Minkowski.spatialMetricTensorDer(0.0, x, 0.0, 0.0, 1.0, 1.0, 1.0)

      local vel = { u, 0.0, 0.0 }
      local v_sq = 0.0

      for i = 1, 3 do
        for j = 1, 3 do
          v_sq = v_sq + (spatial_metric[i][j] * vel[i] * vel[j])
        end
      end

      local W = 1.0 / math.sqrt(1.0 - v_sq)
      if v_sq > 1.0 - math.pow(10.0, -8.0) then
        W = 1.0 / math.sqrt(1.0 - pow(10.0, -8.0))
      end

      local mag = { Bx, By, Bz }

      local cov_mag = { 0.0, 0.0, 0.0 }
      for i = 1, 3 do
        for j = 1, 3 do
          cov_mag[i] = cov_mag[i] + (spatial_metric[i][j] * mag[j])
        end
      end

      local cov_vel = { 0.0, 0.0, 0.0 }
      for i = 1, 3 do
        for j = 1, 3 do
          cov_vel[i] = cov_vel[i] + (spatial_metric[i][j] * vel[j])
        end
      end

      local b0 = 0.0
      for i = 1, 3 do
        b0 = b0 + (W * mag[i] * (cov_vel[i] / lapse))
      end

      local spacetime_vel = { W / lapse, 0.0, 0.0, 0.0 }
      for i = 1, 3 do
        spacetime_vel[i + 1] = (W * vel[i]) - (shift[i] * (W / lapse))
      end

      local b = { 0.0, 0.0, 0.0 }
      for i = 1, 3 do
        b[i] = (mag[i] + (lapse * b0 * spacetime_vel[i + 1])) / W
      end

      local b_sq = 0.0
      for i = 1, 3 do
        b_sq = b_sq + ((mag[i] * cov_mag[i]) / (W * W))
      end
      b_sq = b_sq + (((lapse * lapse) * (b0 * b0)) / (W * W))

      local cov_b = { 0.0, 0.0, 0.0 }
      for i = 1, 3 do
        for j = 1, 3 do
          cov_b[i] = cov_b[i] + (spatial_metric[i][j] * b[j])
        end
      end

      local h_star = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0))) + (b_sq / rho)
      local p_star = p + (0.5 * b_sq)

      local rho_rel = math.sqrt(spatial_det) * rho * W -- Fluid relativistic mass density.
      local mom_x = math.sqrt(spatial_det) * ((rho * h_star * (W * W) * cov_vel[1]) - (lapse * b0 * cov_b[1])) -- Fluid momentum density (x-direction).
      local mom_y = math.sqrt(spatial_det) * ((rho * h_star * (W * W) * cov_vel[2]) - (lapse * b0 * cov_b[2])) -- Fluid momentum density (y-direction).
      local mom_z = math.sqrt(spatial_det) * ((rho * h_star * (W * W) * cov_vel[3]) - (lapse * b0 * cov_b[3])) -- Fluid momentum density (z-direction).
      local Etot = math.sqrt(spatial_det) * ((rho * h_star * (W * W)) - p_star - ((lapse * lapse) * (b0 * b0)) - (rho * W)) -- Fluid total energy density.

      local Bx_rel = math.sqrt(spatial_det) * Bx -- Fluid relativistic magnetic field (x-direction).
      local By_rel = math.sqrt(spatial_det) * By -- Fluid relativistic magnetic field (y-direction).
      local Bz_rel = math.sqrt(spatial_det) * Bz -- Fluid relativistic magnetic field (z-direction).

      local excision = 0.0
      if in_excision_region then
        rho_rel, mom_x, mom_y, mom_z, Etot, Bx_rel, By_rel, Bz_rel, lapse = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
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
    
      return rho_rel, mom_x, mom_y, mom_z, Etot, Bx_rel, By_rel, Bz_rel, 0.0,
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
        x, 0.0, 0.0
    end,

    evolve = true, -- Evolve species?
    limiter = G0.WaveLimiter.MinMod,
    forceLowOrderFlux = false, -- Use HLL fluxes.
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy }, -- Copy boundary conditions (x-direction).
  }
}

-- Run application.
momentApp:run()
