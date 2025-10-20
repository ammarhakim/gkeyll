local Moments = G0.Moments
local VacuumEinstein = G0.Moments.Eq.VacuumEinstein
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

t_end = 10.0 -- Final simulation time.
num_frames = 100 -- Number of output frames.
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
  fluid = Moments.Species.new {
    equation = VacuumEinstein.new {
      excisionThreshold = excision_threshold,
      spacetimeSlicing = spacetime_slicing,
      spacetimeEvolution = spacetime_evolution
    },

    hasVacuumEinstein = true,
    vacuumEinsteinExcisionThreshold = excision_threshold,
    vacuumEinsteinSpacetimeSlicing = spacetime_slicing,
    vacuumEinsteinSpacetimeEvolution = spacetime_evolution,
  
    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]

      local lapse = BrillLindquist.lapseFunction(mass1, mass2, pos_x1, pos_y1, pos_z1, pos_x2, pos_y2, pos_z2, 0.0, x, y, 0.0)
      local shift = BrillLindquist.shiftVector(mass1, mass2, pos_x1, pos_y1, pos_z1, pos_x2, pos_y2, pos_z2, 0.0, x, y, 0.0)
      local spatial_metric = BrillLindquist.spatialMetricTensor(mass1, mass2, pos_x1, pos_y1, pos_z1, pos_x2, pos_y2, pos_z2, 0.0, x, y, 0.0)
      local inv_spatial_metric = BrillLindquist.invSpatialMetricTensor(mass1, mass2, pos_x1, pos_y1, pos_z1, pos_x2, pos_y2, pos_z2, 0.0, x, y, 0.0)
      local spatial_det = BrillLindquist.spatialMetricDeterminant(mass1, mass2, pos_x1, pos_y1, pos_z1, pos_x2, pos_y2, pos_z2, 0.0, x, y, 0.0)
      local extrinsic_curvature = BrillLindquist.extrinsicCurvatureTensor(mass1, mass2, pos_x1, pos_y1, pos_z1, pos_x2, pos_y2, pos_z2, 0.0, x, y, 0.0,
        math.pow(10.0, -8.0), math.pow(10.0, -8.0), math.pow(10.0, -8.0))
      local in_excision_region = BrillLindquist.excisionRegion(mass1, mass2, pos_x1, pos_y1, pos_z1, pos_x2, pos_y2, pos_z2, 0.0, x, y, 0.0)
      
      local lapse_der = BrillLindquist.lapseFunctionDer(mass1, mass2, pos_x1, pos_y1, pos_z1, pos_x2, pos_y2, pos_z2, 0.0, x, y, 0.0,
        math.pow(10.0, -8.0), math.pow(10.0, -8.0), math.pow(10.0, -8.0))
      local shift_der = BrillLindquist.shiftVectorDer(mass1, mass2, pos_x1, pos_y1, pos_z1, pos_x2, pos_y2, pos_z2, 0.0, x, y, 0.0,
        math.pow(10.0, -8.0), math.pow(10.0, -8.0), math.pow(10.0, -8.0))
      local spatial_metric_der = BrillLindquist.spatialMetricTensorDer(mass1, mass2, pos_x1, pos_y1, pos_z1, pos_x2, pos_y2, pos_z2, 0.0, x, y, 0.0,
        math.pow(10.0, -8.0), math.pow(10.0, -8.0), math.pow(10.0, -8.0))

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

    evolve = true, -- Evolve species?
    forceLowOrderFlux = false, -- Use HLL fluxes.
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy }, -- Copy boundary conditions (x-direction).
    bcy = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy } -- Copy boundary conditions (y-direction).
  }
}

-- Run application.
momentApp:run()
