-- 1D current sheet Riemann problem for the general relativistic Maxwell equations.
-- Input parameters taken from the initial conditions in Section C3.2 (current sheet), from the article:
-- S. S. Komissarov (2004), "Electrodynamics of black hole magnetospheres",
-- Monthly Notices of the Royal Astronomical Society, Volume 350 (2): 427-448.
-- https://arxiv.org/abs/astro-ph/0402403

local Moments = G0.Moments
local GRMaxwell = G0.Moments.Eq.GRMaxwell
local Minkowski = G0.Moments.Spacetime.Minkowski

-- Physical constants (using normalized code units).
light_speed = 1.0 -- Speed of light.
e_fact = 0.0 -- Factor of speed of light for electric field correction.
b_fact = 0.0 -- Factor of speed of light for magnetic field correction.

Bx = 1.0 -- Total magnetic field (x-direction).
Bz = 0.0 -- Total magnetic field (z-direction).
Dx = 0.0 -- Total electric field (x-direction).
Dy = 0.0 -- Total electric field (y-direction).
Dz = 0.0 -- Total electric field (z-direction).

B0 = 0.5 -- Reference magnetic field strength.

-- Derived physical quantities (using normalized code units).
By_l = B0 -- Left magnetic field strength (y-direction).
By_r = -B0 -- Right magnetic field strength (z-direction).

-- Simulation parameters.
Nx = 100 -- Cell count (x-direction).
Lx = 3.0 -- Domain size (x-direction).
cfl_frac = 1.0 -- CFL coefficient.

reinit_freq = 10 -- Spacetime reinitialization frequency.

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

  -- Field.
  field = Moments.Species.new {
    equation = GRMaxwell.new {
      lightSpeed = light_speed,
      elcErrorSpeedFactor = e_fact,
      mgnErrorSpeedFactor = b_fact,
      reinitFreq = reinit_freq
    },
  
    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      local By = 0.0

      if x < 0.0 then
        By = By_l -- Left total magnetic field (y-direction).
      else
        By = By_r -- Right total magnetic field (y-direction).
      end

      local lapse = Minkowski.lapseFunction(0.0, x, 0.0, 0.0)
      local shift = Minkowski.shiftVector(0.0, x, 0.0, 0.0)
      local spatial_metric = Minkowski.spatialMetricTensor(0.0, x, 0.0, 0.0)
      local in_excision_region = Minkowski.excisionRegion(0.0, x, 0.0, 0.0)

      local excision = 0.0
      if in_excision_region then
        Dx, Dy, Dz, Bx, By, Bz, lapse = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        for i = 1, 3 do
          shift[i] = 0.0
          for j = 1, 3 do
            spatial_metric[i][j] = 0.0
          end  
        end
        
        excision = -1.0
      else
        excision = 1.0
      end
    
      return Dx, Dy, Dz, Bx, By, Bz, 0.0, 0.0,
        lapse,
        shift[1], shift[2], shift[3],
        spatial_metric[1][1], spatial_metric[1][2], spatial_metric[1][3],
        spatial_metric[2][1], spatial_metric[2][2], spatial_metric[2][3],
        spatial_metric[3][1], spatial_metric[3][2], spatial_metric[3][3],
        excision,
        0.0,
        x, 0.0, 0.0
    end,

    evolve = true, -- Evolve species?
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy } -- Copy boundary conditions (x-direction).
  }
}

-- Run application.
momentApp:run()
