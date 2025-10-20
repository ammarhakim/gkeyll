local Moments = G0.Moments

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.

E0 = 1.0 -- Reference electric field strength.

-- Derived physical quantities (using normalized code units).
light_speed = 1.0 / math.sqrt(mu0 * epsilon0) -- Speed of light.

-- Simulation parameters.
Nr = 80 -- Cell count (radial direction).
Ntheta = 360 -- Cell count (angular direction).
Lr = 1.0 -- Domain size (radial direction).
Ltheta = 2.0 * pi -- Domain size (angular direction).
cfl_frac = 1.0 -- CFL coefficient.

t_end = 0.25 -- Final simulation time.
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
  lower = { 0.25, 0.0 },
  upper = { 0.25 + Lr, Ltheta },
  cells = { Nr, Ntheta },
  cflFrac = cfl_frac,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coodinate direction (radial and angular directions).
      
  -- Boundary conditions for configuration space.
  periodicDirs = { 2 }, -- Periodic directions (angular-direction only).

  -- Computational coordinates (r, theta) to physical coordinates (x, y).
  mapc2p = function (t, zc)
    local r, theta = zc[1], zc[2]

    local xp = { }
    xp[1] = r * math.cos(theta)
    xp[2] = r * math.sin(theta)

    return xp[1], xp[2]
  end,
      
  -- Field.
  field = Moments.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local r, theta = xn[1], xn[2]

      local E_over_r = E0 / r
      local B_over_r = E0 / light_speed / r
    
      local Ex = E_over_r * math.cos(theta) -- Total electric field (x-direction).
      local Ey = E_over_r * math.sin(theta) -- Total electric field (y-direction).
      local Ez = 0.0 -- Total electric field (z-direction).
      
      local Bx = -B_over_r * math.sin(theta) -- Total magnetic field (x-direction).
      local By = B_over_r * math.cos(theta) -- Total magnetic field (y-direction).
      local Bz = 0.0 -- Total magnetic field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

    evolve = true, -- Evolve field?
    limiter = G0.WaveLimiter.NoLimiter,
    bcx = { G0.FieldBc.bcPECWall, G0.FieldBc.bcPECWall } -- PEC wall boundary conditions (radial direction).
  }
}
  
-- Run application.
momentApp:run()
  