-- Non-zero tangential velocity Riemann problem for the general relativistic Euler equations in the tetrad basis, assuming an ultra-relativistic equation of state.
-- Input parameters taken from the initial conditions in Figure 3 (Riemann problem 1), from the article:
-- P. Mach and M. Piętka (2010), "Exact solution of the hydrodynamical Riemann problem with nonzero tangential velocities and the ultrarelativistic equation of state",
-- Physical Review E, Volume 81 (4): 046313.
-- https://arxiv.org/abs/0905.0349

local Moments = G0.Moments
local GRUltraRelEulerTetrad = G0.Moments.Eq.GRUltraRelEulerTetrad
local Minkowski = G0.Moments.Spacetime.Minkowski

-- Physical constants (using normalized code units).
gas_gamma = 4.0 / 3.0 -- Adiabatic index.

rhol = 1.0 -- Left fluid mass density.
ul = 0.5 -- Left fluid normal velocity.
vl = 1.0 / 3.0 -- Left fluid tangential velocity.

rhor = 20.0 -- Right fluid mass density.
ur = 0.5 -- Right fluid normal velocity.
vr = 0.5 -- Right fluid tangential velocity.

-- Simulation parameters.
Nx = 4096 -- Cell count (x-direction).
Lx = 2.0 -- Domain size (x-direction).
cfl_frac = 0.95 -- CFL coefficient.

reinit_freq = 100 -- Spacetime reinitialization frequency.

t_end = 1.0 -- Final simulation time.
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
  lower = { -0.5 * Lx },
  upper = { 0.5 * Lx },
  cells = { Nx },
  cflFrac = cfl_frac,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).
  
  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).

  -- Fluid.
  fluid = Moments.Species.new {
    equation = GRUltraRelEulerTetrad.new {
      gasGamma = gas_gamma,
      reinitFreq = reinit_freq
    },

    hasGRUltraRel = true,
    GRUltraRelGasGamma = gas_gamma,
  
    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]
      
      local rho = 0.0
      local u = 0.0
      local v = 0.0

      if x < 0.0 then
        rho = rhol -- Fluid mass density (left).
        u = ul -- Fluid normal velocity (left).
        v = vl -- Fluid tangential velocity (left).
      else
        rho = rhor -- Fluid mass density (right).
        u = ur -- Fluid normal velocity (right).
        v = vr -- Fluid tangential velocity (right).
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

      local p = (gas_gamma - 1.0) * rho

      local Etot = math.sqrt(spatial_det) * (((rho + p) * (W * W)) - p) -- Fluid total energy density.
      local mom_x = math.sqrt(spatial_det) * (rho + p) * (W * W) * u -- Fluid momentum density (x-direction).
      local mom_y = math.sqrt(spatial_det) * (rho + p) * (W * W) * v -- Fluid momentum density (y-direction).
      local mom_z = 0.0 -- Fluid momentum density (z-direction).

      local excision = 0.0
      if in_excision_region then
        Etot, mom_x, mom_y, mom_z, lapse = 0.0, 0.0, 0.0, 0.0, 0.0
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
    
      return Etot, mom_x, mom_y, mom_z,
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
    forceLowOrderFlux = true, -- Use Lax fluxes.
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy } -- Copy boundary conditions (x-direction).
  }
}

-- Run application.
momentApp:run()
