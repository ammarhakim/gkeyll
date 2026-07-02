-- Auto-generated modular twin: the SAME vacuum-Einstein case driven through
-- Moments.Spacetime.new{ evolve = true, ... } (moment_spacetime flux +
-- moment_spacetime_coupling source) with no fluid species. Cross-checked
-- against the species-driven original (rt_vacuum_einstein_conformal_brill_lindquist.lua).
local Moments = G0.Moments
local VacuumEinsteinConformal = G0.Moments.Eq.VacuumEinsteinConformal
local BrillLindquist = G0.Moments.Spacetime.BrillLindquist

-- Spacetime parameters (using geometric units).
mass1 = 0.5 -- Mass of the first black hole.
mass2 = 0.5 -- Mass of the second black hole.

pos_x1 = 2.5 -- Position of the first black hole (x-direction).
pos_y1 = 5.0 -- Position of the first black hole (y-direction).
pos_z1 = 0.0 -- Position of the first black hole (z-direction).

pos_x2 = 7.5 -- Position of the second black hole (x-direction).
pos_y2 = 5.0 -- Position of the second black hole (y-direction).
pos_z2 = 0.0 -- Position of the second black hole (z-direction).

-- Evolution parameters.
excision_threshold = 0.3 -- Excision threshold (lapse).
spacetime_slicing = G0.SpacetimeSlicing.OnePlusLog -- Spacetime slicing condition.
spacetime_evolution = G0.SpacetimeEvolution.Einstein -- Spacetime evolution system.

-- Simulation parameters.
Nx = 256 -- Cell count (x-direction).
Ny = 256 -- Cell count (y-direction).
Lx = 10.0 -- Domain size (x-direction).
Ly = 10.0 -- Domain size (y-direction).
cfl_frac = 0.8 -- CFL coefficient.

t_end = 9.0 -- Final simulation time.
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
  lower = { 0.0, 0.0 },
  upper = { Lx, Ly },
  cells = { Nx, Ny },
  cflFrac = cfl_frac,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coodinate direction (x- and y-directions).
  
  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).

  -- Fluid.
  spacetime = Moments.Spacetime.new {
    evolve = true,
    equation = VacuumEinsteinConformal.new {
      excisionThreshold = excision_threshold,
      spacetimeSlicing = spacetime_slicing,
      spacetimeEvolution = spacetime_evolution
    },

  
    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local conformal_lapse = BrillLindquist.lapseFunction(mass1, mass2, pos_x1, pos_y1, pos_z1, pos_x2, pos_y2, pos_z2, 0.0, x, y, 0.0)
      local conformal_shift = BrillLindquist.shiftVector(mass1, mass2, pos_x1, pos_y1, pos_z1, pos_x2, pos_y2, pos_z2, 0.0, x, y, 0.0)
      local conformal_spatial_metric = BrillLindquist.spatialMetricTensor(mass1, mass2, pos_x1, pos_y1, pos_z1, pos_x2, pos_y2, pos_z2, 0.0, x, y, 0.0)
      local inv_conformal_spatial_metric = BrillLindquist.invSpatialMetricTensor(mass1, mass2, pos_x1, pos_y1, pos_z1, pos_x2, pos_y2, pos_z2, 0.0, x, y, 0.0)
      local conformal_spatial_det = BrillLindquist.spatialMetricDeterminant(mass1, mass2, pos_x1, pos_y1, pos_z1, pos_x2, pos_y2, pos_z2, 0.0, x, y, 0.0)
      local conformal_extrinsic_curvature = BrillLindquist.extrinsicCurvatureTensor(mass1, mass2, pos_x1, pos_y1, pos_z1, pos_x2, pos_y2, pos_z2, 0.0, x, y, 0.0,
        math.pow(10.0, -8.0), math.pow(10.0, -8.0), math.pow(10.0, -8.0))
      local in_excision_region = BrillLindquist.excisionRegion(mass1, mass2, pos_x1, pos_y1, pos_z1, pos_x2, pos_y2, pos_z2, 0.0, x, y, 0.0)
      
      local conformal_lapse_der = BrillLindquist.lapseFunctionDer(mass1, mass2, pos_x1, pos_y1, pos_z1, pos_x2, pos_y2, pos_z2, 0.0, x, y, 0.0,
        math.pow(10.0, -8.0), math.pow(10.0, -8.0), math.pow(10.0, -8.0))
      local conformal_shift_der = BrillLindquist.shiftVectorDer(mass1, mass2, pos_x1, pos_y1, pos_z1, pos_x2, pos_y2, pos_z2, 0.0, x, y, 0.0,
        math.pow(10.0, -8.0), math.pow(10.0, -8.0), math.pow(10.0, -8.0))
      local conformal_spatial_metric_der = BrillLindquist.spatialMetricTensorDer(mass1, mass2, pos_x1, pos_y1, pos_z1, pos_x2, pos_y2, pos_z2, 0.0, x, y, 0.0,
        math.pow(10.0, -8.0), math.pow(10.0, -8.0), math.pow(10.0, -8.0))
      
      local conformal_fact = BrillLindquist.conformalFactor(mass1, mass2, pos_x1, pos_y1, pos_z1, pos_x2, pos_y2, pos_z2, 0.0, x, y, 0.0)
      local bssn_conformal_fact = BrillLindquist.bssnConformalFactor(mass1, mass2, pos_x1, pos_y1, pos_z1, pos_x2, pos_y2, pos_z2, 0.0, x, y, 0.0)

      local conformal_fact_der = BrillLindquist.conformalFactorDer(mass1, mass2, pos_x1, pos_y1, pos_z1, pos_x2, pos_y2, pos_z2, 0.0, x, y, 0.0,
        math.pow(10.0, -8.0), math.pow(10.0, -8.0), math.pow(10.0, -8.0))
      local bssn_conformal_fact_der = BrillLindquist.bssnConformalFactorDer(mass1, mass2, pos_x1, pos_y1, pos_z1, pos_x2, pos_y2, pos_z2, 0.0, x, y, 0.0,
        math.pow(10.0, -8.0), math.pow(10.0, -8.0), math.pow(10.0, -8.0))
      local bssn_conformal_fact_der2 = BrillLindquist.bssnConformalFactorDer2(mass1, mass2, pos_x1, pos_y1, pos_z1, pos_x2, pos_y2, pos_z2, 0.0, x, y, 0.0,
        math.pow(10.0, -8.0), math.pow(10.0, -8.0), math.pow(10.0, -8.0))
    
      for i = 1, 3 do
        conformal_fact_der[i] = conformal_fact_der[i] / conformal_fact
      end

      -- Set first and second conformal derivatives to zero, to improve stability.
      for i = 1, 3 do
        conformal_fact_der[i] = 0.0
        bssn_conformal_fact_der[i] = 0.0

        for j = 1, 3 do
          bssn_conformal_fact_der2[i][j] = 0.0
        end
      end

      for i = 1, 3 do
        for j = 1, 3 do
          conformal_spatial_metric[i][j] = conformal_spatial_metric[i][j] / (conformal_fact * conformal_fact * conformal_fact * conformal_fact)
          inv_conformal_spatial_metric[i][j] = inv_conformal_spatial_metric[i][j] * (conformal_fact * conformal_fact * conformal_fact * conformal_fact)
        end
      end

      for i = 1, 3 do
        for j = 1, 3 do
          for k = 1, 3 do
            conformal_spatial_metric_der[i][j][k] = 0.5 * conformal_spatial_metric_der[i][j][k]
            conformal_spatial_metric_der[i][j][k] = conformal_spatial_metric_der[i][j][k] / (conformal_fact * conformal_fact * conformal_fact * conformal_fact)
            conformal_spatial_metric_der[i][j][k] = conformal_spatial_metric_der[i][j][k] - (2.0 * conformal_fact_der[i] * conformal_spatial_metric[j][k]);

          end

        conformal_shift_der[i][j] = 0.5 * conformal_shift_der[i][j]
        end
      end

      for i = 1, 3 do
        conformal_lapse_der[i] = conformal_lapse_der[i] / conformal_lapse
      end

      local conformal_spatial_metric_der_raised1 = { }
      for i = 1, 3 do
        conformal_spatial_metric_der_raised1[i] = { }
        for j = 1, 3 do
          conformal_spatial_metric_der_raised1[i][j] = { }
        end
      end
      for i = 1, 3 do
        for j = 1, 3 do
          for k = 1, 3 do
            conformal_spatial_metric_der_raised1[k][i][j] = 0.0
        
            for l = 1, 3 do
              conformal_spatial_metric_der_raised1[k][i][j] = conformal_spatial_metric_der_raised1[k][i][j] + inv_conformal_spatial_metric[k][l] * conformal_spatial_metric_der[l][i][j]
            end
          end
        end
      end
  
      local conformal_spatial_metric_der_raised3 = { }
      for i = 1, 3 do
        conformal_spatial_metric_der_raised3[i] = { }
        for j = 1, 3 do
          conformal_spatial_metric_der_raised3[i][j] = { }
        end
      end
      for i = 1, 3 do
        for j = 1, 3 do
          for k = 1, 3 do
            conformal_spatial_metric_der_raised3[i][j][k] = 0.0

            for l = 1, 3 do
              conformal_spatial_metric_der_raised3[i][j][k] = conformal_spatial_metric_der_raised3[i][j][k] + inv_conformal_spatial_metric[l][k] * conformal_spatial_metric_der[i][j][l]
            end
          end
        end
      end

      local conformal_aux_vect = { }
      for i = 1, 3 do
        conformal_aux_vect[i] = 0.0

        for s = 1, 3 do
          conformal_aux_vect[i] = conformal_aux_vect[i] + conformal_spatial_metric_der_raised3[i][s][s]
          conformal_aux_vect[i] = conformal_aux_vect[i] - conformal_spatial_metric_der_raised1[s][s][i]
        end

        conformal_aux_vect[i] = conformal_aux_vect[i] - 4.0 * conformal_fact_der[i]
      end

      local excision = 0.0
      if in_excision_region then
        conformal_lapse = 0.0
        conformal_fact = 0.0
        bssn_conformal_fact = 0.0
        for i = 1, 3 do
          conformal_shift[i] = 0.0
          conformal_lapse_der[i] = 0.0
          conformal_aux_vect[i] = 0.0
          conformal_fact_der[i] = 0.0
          bssn_conformal_fact_der[i] = 0.0
          for j = 1, 3 do
            conformal_spatial_metric[i][j] = 0.0
            conformal_extrinsic_curvature[i][j] = 0.0
            conformal_shift_der[i][j] = 0.0
            bssn_conformal_fact_der2[i][j] = 0.0
            for k = 1, 3 do
              conformal_spatial_metric_der[i][j][k] = 0.0
            end
          end  
        end
        
        excision = -1.0
      else
        excision = 1.0
      end
    
      return conformal_spatial_metric[1][1], conformal_spatial_metric[1][2], conformal_spatial_metric[1][3],
        conformal_spatial_metric[2][1], conformal_spatial_metric[2][2], conformal_spatial_metric[2][3],
        conformal_spatial_metric[3][1], conformal_spatial_metric[3][2], conformal_spatial_metric[3][3],
        conformal_lapse,
        conformal_extrinsic_curvature[1][1], conformal_extrinsic_curvature[1][2], conformal_extrinsic_curvature[1][3],
        conformal_extrinsic_curvature[2][1], conformal_extrinsic_curvature[2][2], conformal_extrinsic_curvature[2][3],
        conformal_extrinsic_curvature[3][1], conformal_extrinsic_curvature[3][2], conformal_extrinsic_curvature[3][3],
        conformal_spatial_metric_der[1][1][1], conformal_spatial_metric_der[1][1][2], conformal_spatial_metric_der[1][1][3],
        conformal_spatial_metric_der[1][2][1], conformal_spatial_metric_der[1][2][2], conformal_spatial_metric_der[1][2][3],
        conformal_spatial_metric_der[1][3][1], conformal_spatial_metric_der[1][3][2], conformal_spatial_metric_der[1][3][3],
        conformal_spatial_metric_der[2][1][1], conformal_spatial_metric_der[2][1][2], conformal_spatial_metric_der[2][1][3],
        conformal_spatial_metric_der[2][2][1], conformal_spatial_metric_der[2][2][2], conformal_spatial_metric_der[2][2][3],
        conformal_spatial_metric_der[2][3][1], conformal_spatial_metric_der[2][3][2], conformal_spatial_metric_der[2][3][3],
        conformal_spatial_metric_der[3][1][1], conformal_spatial_metric_der[3][1][2], conformal_spatial_metric_der[3][1][3],
        conformal_spatial_metric_der[3][2][1], conformal_spatial_metric_der[3][2][2], conformal_spatial_metric_der[3][2][3],
        conformal_spatial_metric_der[3][3][1], conformal_spatial_metric_der[3][3][2], conformal_spatial_metric_der[3][3][3],
        conformal_lapse_der[1], conformal_lapse_der[2], conformal_lapse_der[3],
        conformal_aux_vect[1], conformal_aux_vect[2], conformal_aux_vect[3],
        conformal_shift[1], conformal_shift[2], conformal_shift[3],
        conformal_shift_der[1][1], conformal_shift_der[1][2], conformal_shift_der[1][3],
        conformal_shift_der[2][1], conformal_shift_der[2][2], conformal_shift_der[2][3],
        conformal_shift_der[3][1], conformal_shift_der[3][2], conformal_shift_der[3][3],
        bssn_conformal_fact,
        bssn_conformal_fact_der[1], bssn_conformal_fact_der[2], bssn_conformal_fact_der[3],
        bssn_conformal_fact_der2[1][1], bssn_conformal_fact_der2[1][2], bssn_conformal_fact_der2[1][3],
        bssn_conformal_fact_der2[2][1], bssn_conformal_fact_der2[2][2], bssn_conformal_fact_der2[2][3],
        bssn_conformal_fact_der2[3][1], bssn_conformal_fact_der2[3][2], bssn_conformal_fact_der2[3][3]
    end,

    forceLowOrderFlux = false, -- Use HLL fluxes.
  }
}

-- Run application.
momentApp:run()
