-- Square wave initial data test for the inviscid Burgers' equation.
-- Input parameters correspond to a square wave, equivalent to the initial conditions in Section 6, but with the wave amplitude set to 3 between x = 2 and x = 4, and -1 everywhere else, from the article:
-- J. W. L. Wan and A. Jameson (2008), "Monotonicity-preserving multigrid time stepping schemes for conservation laws",
-- Computing and Visualization in Science, Volume 11: 41-58.
-- https://link.springer.com/article/10.1007/s00791-006-0056-3

local Moments = G0.Moments
local Burgers = G0.Moments.Eq.Burgers

-- Simulation parameters.
Nx = 128 -- Cell count (x-direction).
Lx = 6.0 -- Domain size (x-direction).
cfl_frac = 0.9 -- CFL coefficient.

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
  
  -- Advected quantity.
  fluid = Moments.Species.new {
    equation = Burgers.new { },
    
    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      local f = 0.0

      if x > 2.0 and x < 4.0 then
        f = 3.0 -- Advected quantity (between 2 and 4).
      else
        f = -1.0 -- Advected quantity (elsewhere).
      end
      
      return f
    end,
  
    evolve = true, -- Evolve species?
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy } -- Copy boundary conditions (x-direction).
  }
}

-- Run application.
momentApp:run()