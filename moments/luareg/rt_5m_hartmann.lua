local Moments = G0.Moments
local Euler = G0.Moments.Eq.Euler

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
gas_gamma = 5.0 / 3.0 -- Adiabatic index.
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_ion = 1.0 -- Ion mass.
charge_ion = 1.0 -- Ion charge.
mass_elc = 1.0 / math.sqrt(1836.153) -- Electron mass.
charge_elc = -1.0 -- Electron charge.

n0 = 1.0 -- Reference number density.
coll_fac = 1e5 -- Collision factor.

beta = 0.1 -- Plasma beta.
grav = 0.01 -- Gravitational acceleration.

-- Derived physical quantities (using normalized code units).
light_speed = 1.0 / math.sqrt(epsilon0 * mu0) -- Speed of light.
vAe = light_speed -- Electron Alfven velocity.
B0 = vAe * math.sqrt(mu0 * n0 * mass_elc) -- Reference magnetic field strength.
omega_pi = math.sqrt(charge_ion * charge_ion * n0 / (epsilon0 * mass_ion)) -- Ion plasma frequency.
di = light_speed / omega_pi -- Ion skin depth.

T_elc = beta * (B0 * B0) / (2.0 * n0 * mu0) -- Electron temperature.
T_ion = beta * (B0 * B0) / (2.0 * n0 * mu0) -- Ion temperature.

rho_elc = n0 * mass_elc -- Electron mass density.
rho_ion = n0 * mass_ion -- Ion mass density.

E_elc = n0 * T_elc / (gas_gamma - 1.0) -- Electron total energy density.
E_ion = n0 * T_ion / (gas_gamma - 1.0) -- Ion total energy density.

tau = coll_fac * 6.0 * math.sqrt(2.0 * pi * mass_elc * T_elc * pi * T_elc * pi * T_elc) * epsilon0 * epsilon0 /
  (charge_ion * charge_ion * charge_ion * charge_ion * n0) -- Collision time.
lambda = 1.0 / mu0 * (mass_elc / (charge_ion * charge_ion * n0 * tau)) -- Collision wavelength.

-- Simulation parameters.
Nx = 256 -- Cell count (x-direction).
Lx = 256.0 * di -- Domain size (x-direction).
cfl_frac = 0.001 -- CFL coefficient.

t_end = 1.0 * tau -- Final simulation time.
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
    equation = Euler.new { gasGamma = gas_gamma },
    braginskiiType = G0.Braginskii.MagFull,
    
    -- Initial conditions function.
    init = function (t, xn)
      local rhoe = rho_elc -- Electron mass density.
      local mome_x = 0.0 -- Electron momentum density (x-direction).
      local mome_y = 0.0 -- Electron momentum density (y-direction).
      local mome_z = 0.0 -- Electron momentum density (z-direction).
      local Ee_tot = E_elc -- Electron total energy density.
      
      return rhoe, mome_x, mome_y, mome_z, Ee_tot
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
    equation = Euler.new { gasGamma = gas_gamma },
    braginskiiType = G0.Braginskii.MagFull,
    
    -- Initial conditions function.
    init = function (t, xn)
      local rhoi = rho_ion -- Ion mass density.
      local momi_x = 0.0 -- Ion momentum density (x-direction).
      local momi_y = 0.0 -- Ion momentum density (y-direction).
      local momi_z = 0.0 -- Ion momentum density (z-direction).
      local Ei_tot = E_ion -- Ion total energy density.
      
      return rhoi, momi_x, momi_y, momi_z, Ei_tot
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

      local Bx = B0 -- Total magnetic field (x-direction).
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
