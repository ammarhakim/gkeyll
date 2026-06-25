-- Auto-generated modular twin: the SAME vacuum-Einstein case driven through
-- Moments.Spacetime.new{ evolve = true, ... } (moment_spacetime flux +
-- moment_spacetime_coupling source) with no fluid species. Cross-checked
-- against the species-driven original (rt_vacuum_einstein_gowdywave.lua).
local Moments = G0.Moments
local VacuumEinstein = G0.Moments.Eq.VacuumEinstein
local Minkowski = G0.Moments.Spacetime.Minkowski
local FFI = require("ffi")

-- Define Bessel functions.
FFI.cdef[[
  double jn(int n, double x);
]]

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
tau0 = 9.8753205829098 -- Reparameterized initial time.

-- Evolution parameters.
excision_threshold = 0.3 -- Excision threshold (lapse).
spacetime_slicing = G0.SpacetimeSlicing.Harmonic -- Spacetime slicing condition.
spacetime_evolution = G0.SpacetimeEvolution.Einstein -- Spacetime evolution system.

-- Simulation parameters.
Nx = 50 -- Cell count (x-direction).
Lx = 1.0 -- Domain size (x-direction).
cfl_frac = 0.25 -- CFL coefficient.

t_end = 2.0 -- Final simulation time.
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
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x- and y-directions).
  
  -- Boundary conditions for configuration space.
  periodicDirs = { 1 }, -- Periodic directions (x-direction only).

  -- Fluid.
  spacetime = Moments.Spacetime.new {
    evolve = true,
    equation = VacuumEinstein.new {
      excisionThreshold = excision_threshold,
      spacetimeSlicing = spacetime_slicing,
      spacetimeEvolution = spacetime_evolution
    },

  
    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local lapse = Minkowski.lapseFunction(0.0, x, 0.0, 0.0)
      local shift = Minkowski.shiftVector(0.0, x, 0.0, 0.0)
      local spatial_metric = Minkowski.spatialMetricTensor(0.0, x, 0.0, 0.0)
      local inv_spatial_metric = Minkowski.invSpatialMetricTensor(0.0, x, 0.0, 0.0)
      local spatial_det = Minkowski.spatialMetricDeterminant(0.0, x, 0.0, 0.0)
      local extrinsic_curvature = Minkowski.extrinsicCurvatureTensor(0.0, x, 0.0, 0.0, 1.0, 1.0, 1.0)
      local in_excision_region = Minkowski.excisionRegion(0.0, x, 0.0, 0.0)

      local lapse_der = Minkowski.lapseFunctionDer(0.0, x, 0.0, 0.0, 1.0, 1.0, 1.0)
      local shift_der = Minkowski.shiftVectorDer(0.0, x, 0.0, 0.0, 1.0, 1.0, 1.0)
      local spatial_metric_der = Minkowski.spatialMetricTensorDer(0.0, x, 0.0, 0.0, 1.0, 1.0, 1.0)

      lambda = (-2.0 * pi * tau0 * FFI.C.jn(0, 2.0 * pi * tau0) * FFI.C.jn(1, 2.0 * pi * tau0) * (math.cos(2.0 * pi * x) * math.cos(2.0 * pi * x)))
        + (2.0 * (pi * pi) * (tau0 * tau0) * ((FFI.C.jn(0, 2.0 * pi * tau0) * FFI.C.jn(0, 2.0 * pi * tau0)) + (FFI.C.jn(1, 2.0 * pi * tau0) * FFI.C.jn(1, 2.0 * pi * tau0))))
        - (0.5 * (((2.0 * pi) * (2.0 * pi)) * ((FFI.C.jn(0, 2.0 * pi) * FFI.C.jn(0, 2.0 * pi)) + (FFI.C.jn(1, 2.0 * pi) * FFI.C.jn(1, 2.0 * pi)))))
        + (pi * FFI.C.jn(0, 2.0 * pi) * FFI.C.jn(1, 2.0 * pi))
      lambda_dt = 2.0 * (pi * pi) * tau0 * ((FFI.C.jn(1, 2.0 * pi * tau0) * FFI.C.jn(1, 2.0 * pi * tau0)) * (1.0 + math.cos(4.0 * pi * x))
        + (2.0 * (FFI.C.jn(0, 2.0 * pi * tau0) * FFI.C.jn(0, 2.0 * pi * tau0)) * (math.sin(2.0 * pi * x) * math.sin(2.0 * pi * x))))
      lambda_dx = 4.0 * (pi * pi * tau0) * FFI.C.jn(0, 2.0 * pi * tau0) * FFI.C.jn(1, 2.0 * pi * tau0) * math.sin(4.0 * pi * x)

      P = FFI.C.jn(0, 2.0 * pi * tau0) * math.cos(2.0 * pi * x)
      P_dt = -2.0 * pi * FFI.C.jn(1, 2.0 * pi * tau0) * math.cos(2.0 * pi * x)
      P_dx = -2.0 * pi * FFI.C.jn(0, 2.0 * pi * tau0) * math.sin(2.0 * pi * x)

      spatial_metric[1][1] = (1.0 / math.sqrt(tau0)) * math.exp(0.5 * lambda)
      spatial_metric[2][2] = tau0 * math.exp(P)
      spatial_metric[3][3] = tau0 * math.exp(-P)

      lapse = (1.0 / math.pow(tau0, 0.25)) * math.exp(0.25 * lambda)
      lapse_der[0] = (math.exp(0.25 * lambda) * lambda_dx) / (4.0 * math.pow(tau0, 0.25))

      extrinsic_curvature[1][1] = 0.25 * (1.0 / math.pow(tau0, 0.25)) * math.exp(0.25 * lambda) * ((1.0 / tau0) - lambda_dt)
      extrinsic_curvature[2][2] = -0.5 * math.pow(tau0, 0.25) * math.exp(-0.25 * lambda) * math.exp(P) * (1.0 + (tau0 * P_dt))
      extrinsic_curvature[3][3] = -0.5 * math.pow(tau0, 0.25) * math.exp(-0.25 * lambda) * math.exp(-P) * (1.0 - (tau0 * P_dt))

      spatial_metric_der[1][1][1] = (math.exp(0.25 * lambda) * lambda_dx) / (2.0 * math.sqrt(tau0))
      spatial_metric_der[1][2][2] = math.exp(P) * tau0 * P_dx
      spatial_metric_der[1][3][3] = -math.exp(-P) * tau0 * P_dx

      spatial_det = math.exp(0.5 * lambda) * math.pow(tau0, 1.5)
      inv_spatial_metric[1][1] = math.exp(-0.5 * lambda) * math.sqrt(tau0)
      inv_spatial_metric[2][2] = math.exp(-P) / tau0
      inv_spatial_metric[3][3] = math.exp(P) / tau0

      for i = 1, 3 do
        for j = 1, 3 do
          for k = 1, 3 do
            spatial_metric_der[i][j][k] = 0.5 * spatial_metric_der[i][j][k]
          end

        shift_der[i][j] = 0.5 * shift_der[i][j]
        end
      end

      for i = 1, 3 do
        lapse_der[i] = lapse_der[i] / lapse
      end

      local spatial_metric_der_raised1 = { }
      for i = 1, 3 do
        spatial_metric_der_raised1[i] = { }
        for j = 1, 3 do
          spatial_metric_der_raised1[i][j] = { }
        end
      end
      for i = 1, 3 do
        for j = 1, 3 do
          for k = 1, 3 do
            spatial_metric_der_raised1[k][i][j] = 0.0
        
            for l = 1, 3 do
              spatial_metric_der_raised1[k][i][j] = spatial_metric_der_raised1[k][i][j] + inv_spatial_metric[k][l] * spatial_metric_der[l][i][j]
            end
          end
        end
      end
  
      local spatial_metric_der_raised3 = { }
      for i = 1, 3 do
        spatial_metric_der_raised3[i] = { }
        for j = 1, 3 do
          spatial_metric_der_raised3[i][j] = { }
        end
      end
      for i = 1, 3 do
        for j = 1, 3 do
          for k = 1, 3 do
            spatial_metric_der_raised3[i][j][k] = 0.0

            for l = 1, 3 do
              spatial_metric_der_raised3[i][j][k] = spatial_metric_der_raised3[i][j][k] + inv_spatial_metric[l][k] * spatial_metric_der[i][j][l]
            end
          end
        end
      end

      local aux_vect = { }
      for i = 1, 3 do
        aux_vect[i] = 0.0

        for s = 1, 3 do
          aux_vect[i] = aux_vect[i] + spatial_metric_der_raised3[i][s][s]
          aux_vect[i] = aux_vect[i] - spatial_metric_der_raised1[s][s][i]
        end
      end

      local excision = 0.0
      if in_excision_region then
        lapse = 0.0
        for i = 1, 3 do
          shift[i] = 0.0
          lapse_der[i] = 0.0
          aux_vect[i] = 0.0
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
    
      return spatial_metric[1][1], spatial_metric[1][2], spatial_metric[1][3],
        spatial_metric[2][1], spatial_metric[2][2], spatial_metric[2][3],
        spatial_metric[3][1], spatial_metric[3][2], spatial_metric[3][3],
        lapse,
        extrinsic_curvature[1][1], extrinsic_curvature[1][2], extrinsic_curvature[1][3],
        extrinsic_curvature[2][1], extrinsic_curvature[2][2], extrinsic_curvature[2][3],
        extrinsic_curvature[3][1], extrinsic_curvature[3][2], extrinsic_curvature[3][3],
        spatial_metric_der[1][1][1], spatial_metric_der[1][1][2], spatial_metric_der[1][1][3],
        spatial_metric_der[1][2][1], spatial_metric_der[1][2][2], spatial_metric_der[1][2][3],
        spatial_metric_der[1][3][1], spatial_metric_der[1][3][2], spatial_metric_der[1][3][3],
        spatial_metric_der[2][1][1], spatial_metric_der[2][1][2], spatial_metric_der[2][1][3],
        spatial_metric_der[2][2][1], spatial_metric_der[2][2][2], spatial_metric_der[2][2][3],
        spatial_metric_der[2][3][1], spatial_metric_der[2][3][2], spatial_metric_der[2][3][3],
        spatial_metric_der[3][1][1], spatial_metric_der[3][1][2], spatial_metric_der[3][1][3],
        spatial_metric_der[3][2][1], spatial_metric_der[3][2][2], spatial_metric_der[3][2][3],
        spatial_metric_der[3][3][1], spatial_metric_der[3][3][2], spatial_metric_der[3][3][3],
        lapse_der[1], lapse_der[2], lapse_der[3],
        aux_vect[1], aux_vect[2], aux_vect[3],
        shift[1], shift[2], shift[3],
        shift_der[1][1], shift_der[1][2], shift_der[1][3],
        shift_der[2][1], shift_der[2][2], shift_der[2][3],
        shift_der[3][1], shift_der[3][2], shift_der[3][3]
    end,

    forceLowOrderFlux = false -- Use HLL fluxes.
  }
}

-- Run application.
momentApp:run()
