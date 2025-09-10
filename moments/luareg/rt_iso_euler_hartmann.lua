local Moments = G0.Moments
local IsoEuler = G0.Moments.Eq.IsoEuler

-- Physical constants (using normalized code units).
epsilon0 = 8.854187817620389850536563031710750260608e-12 -- Permittivity of free space.
mu0 = 12.56637061435917295385057353311801153679e-7 -- Permeability of free space.
mass_ion = 1.672621637e-27 -- Ion mass.
charge_ion = 1.602176487e-19 -- Ion charge.
mass_elc = 9.10938215e-31 -- Electron mass.
charge_elc = -1.602176487e-19 -- Electron charge.

vte = 3.0e6 -- Electron thermal velocity.
vti = 6.0e4 -- Ion thermal velocity.
B0 = 10.0 -- Reference magnetic field strength.
Bx = 1.0 -- Magnetic field (x-direction).

beta = 0.01 -- Plasma beta.
grav = 9.8 -- Gravitational acceleration.
coll_fac = 1.0 -- Collision factor.

-- Derived physical quantities (using normalized code units).
T_elc = mass_elc * vte * vte / 2.0 -- Electron temperature.
T_ion = mass_ion * vti * vti / 2.0 -- Ion temperature.

n_elc = beta * B0 * B0 / (2.0 * mu0 * T_elc) -- Electron number density. 
n_ion = beta * B0 * B0 / (2.0 * mu0 * T_elc) -- Ion number density.

rho_elc = n_elc * mass_elc -- Electron mass density.
rho_ion = n_ion * mass_ion -- Ion mass density.

-- Simulation parameters.
Nx = 64 -- Cell count (x-direction).
Lx = 1.0 -- Domain size (x-direction).
cfl_frac = 1.0 -- CFL coefficient.

t_end = 1.0e-6 -- Final simulation time.
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

  hasBraginskii = true,
  collisionFactor = coll_fac,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).
    
  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).
  
  -- Electrons.
  elc = Moments.Species.new {
    charge = charge_elc, mass = mass_elc,
    equation = IsoEuler.new { vThermal = vte },
    braginskiiType = G0.Braginskii.MagFull,
    
    -- Initial conditions function.
    init = function (t, xn)
      local rhoe = rho_elc -- Electron mass density.
      local mome_x = 0.0 -- Electron momentum density (x-direction).
      local mome_y = 0.0 -- Electron momentum density (y-direction).
      local mome_z = 0.0 -- Electron momentum density (z-direction).
      
      return rhoe, mome_x, mome_y, mome_z
    end,

    -- Applied acceleration function.
    appliedAcceleration = function (t, xn)
      local accele_x = 0.0 -- Electron applied acceleration (x-direction).
      local accele_y = grav -- Electron applied acceleration (y-direction).
      local accele_z = 0.0 -- Electron applied acceleration (z-direction).

      return accele_x, accele_y, accele_z
    end,
  
    evolve = true, -- Evolve species?
    bcx = { G0.SpeciesBc.bcNoSlip, G0.SpeciesBc.bcNoSlip } -- No-slip boundary conditions (x-direction).
  },

  -- Ions.
  ion = Moments.Species.new {
    charge = charge_ion, mass = mass_ion,
    equation = IsoEuler.new { vThermal = vti },
    braginskiiType = G0.Braginskii.MagFull,
    
    -- Initial conditions function.
    init = function (t, xn)
      local rhoi = rho_ion -- Ion mass density.
      local momi_x = 0.0 -- Ion momentum density (x-direction).
      local momi_y = 0.0 -- Ion momentum density (y-direction).
      local momi_z = 0.0 -- Ion momentum density (z-direction).
      
      return rhoi, momi_x, momi_y, momi_z
    end,

    -- Applied acceleration function.
    appliedAcceleration = function (t, xn)
      local acceli_x = 0.0 -- Ion applied acceleration (x-direction).
      local acceli_y = grav -- Ion applied acceleration (y-direction).
      local acceli_z = 0.0 -- Ion applied acceleration (z-direction).
  
      return acceli_x, acceli_y, acceli_z
    end,
  
    evolve = true, -- Evolve species?
    bcx = { G0.SpeciesBc.bcNoSlip, G0.SpeciesBc.bcNoSlip } -- No-slip boundary conditions (x-direction).
  },
    
  -- Field.
  field = Moments.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local Ex = 0.0 -- Total electric field (x-direction).
      local Ey = 0.0 -- Total electric field (y-direction).
      local Ez = 0.0 -- Total electric field (z-direction).

      local Bx = Bx -- Total magnetic field (x-direction).
      local By = 0.0 -- Total magnetic field (y-direction).
      local Bz = 0.0 -- Total magnetic field (z-direction).
      
      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

    evolve = false, -- Evolve field?
    bcx = { G0.FieldBc.bcPECWall, G0.FieldBc.bcPECWall } -- PEC wall boundary conditions (x-direction).
  }
}

-- Run application.
momentApp:run()
