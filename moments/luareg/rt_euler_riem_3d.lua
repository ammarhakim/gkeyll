-- 3D (spherical) Riemann problem for the 5-moment (Euler) equations.
-- Input parameters match the initial conditions found in Entry JE23 of Ammar's Simulation Journal (https://www.ammar-hakim.org/sj/je/je23/je23-euler-3d.html#spherical-riemann-problem),
-- adapted from Section 3.2, from the article:
-- J. O. Langseth and R. J. LeVeque (2000), "A Wave Propagation Method for Three-Dimensional Hyperbolic Conservation Laws",
-- Journal of Computational Physics, Volume 165 (1): 126-166.
-- https://www.sciencedirect.com/science/article/pii/S0021999100966063

local Moments = G0.Moments
local Euler = G0.Moments.Eq.Euler

-- Physical constants (using normalized code units).
gas_gamma = 1.4 -- Adiabatic index.

rhol = 1.0 -- Left/inner fluid mass density.
pl = 5.0 -- Left/inner fluid pressure.

rhor = 1.0 -- Right/outer fluid mass density.
pr = 1.0 -- Right/outer fluid pressure.

-- Simulation parameters.
Nx = 37 -- Cell count (x-direction).
Ny = 37 -- Cell count (y-direction).
Nz = 25 -- Cell count (z-direction).
Lx = 1.5 -- Domain size (x-direction).
Ly = 1.5 -- Domain size (y-direction).
Lz = 1.0 -- Domain size (z-direction).
cfl_frac = 0.9 -- CFL coefficient.

t_end = 0.7 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

rloc = 0.2 -- Fluid boundary (radial coordinate).

momentApp = Moments.App.new {
  
  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { 0.0, 0.0, 0.0 },
  upper = { Lx, Ly, Lz },
  cells = { Nx, Ny, Nz },
  cflFrac = cfl_frac,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1, 1 }, -- Cuts in each coodinate direction (x-, y- and z-directions).
    
  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).
  
  -- Fluid.
  fluid = Moments.Species.new {
    equation = Euler.new { gasGamma = gas_gamma },
    
    -- Initial conditions function.
    init = function (t, xn)
      local x, y, z = xn[1], xn[2], xn[3]

      local rho = 0.0
      local p = 0.0
    
      local r = math.sqrt(x * x + y * y + (z - 0.4) * (z - 0.4))
    
      if r < rloc then
        rho = rhol -- Fluid mass density (left/inner).
        p = pl -- Fluid pressure (left/inner).
      else
        rho = rhor -- Fluid mass density (right/outer).
        p = pr -- Fluid pressure (right/outer).
      end
      
      local mom_x = 0.0 -- Fluid momentum density (x-direction).
      local mom_y = 0.0 -- Fluid momentum density (y-direction).
      local mom_z = 0.0 -- Fluid momentum density (z-direction).
      local Etot = p / (gas_gamma - 1.0) -- Fluid total energy density.
      
      return rho, mom_x, mom_y, mom_z, Etot
    end,
  
    evolve = true, -- Evolve species?
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy }, -- Copy boundary conditions (x-direction).
    bcy = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy }, -- Copy boundary conditions (y-direction).
    bcz = { G0.SpeciesBc.bcWall, G0.SpeciesBc.bcWall } -- Wall boundary conditions (y-direction).
  }
}

-- Run application.
momentApp:run()