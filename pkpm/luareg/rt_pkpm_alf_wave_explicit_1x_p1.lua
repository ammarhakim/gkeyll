local PKPM = G0.PKPM

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.
mass_ion = 100.0 -- Ion mass.
charge_ion = 1.0 -- Ion charge.

n0 = 1.0 -- Reference number density.
vAe = 0.1 -- Electron Alfven velocity.
beta = 0.1 -- Plasma beta.

-- Derived physical quantities (using normalized code units).
B0 = vAe * math.sqrt(mu0 * n0 * mass_elc) -- Reference magnetic field strength.
vAi = vAe / math.sqrt(mass_ion) -- Ion Alfven velocity.

vte = vAe * math.sqrt(beta / 2.0) -- Electron thermal velocity.
vti = vte / math.sqrt(mass_ion) -- Ion thermal velocity.

omega_ci = charge_ion * B0 / mass_ion -- Ion cyclotron frequency.
rho_s = math.sqrt(2.0) * vti / omega_ci -- Ion sound gyroradius.

k_perp = 1.005412 / rho_s -- Perpendicular wavenumber.
k_par = 0.087627 / rho_s -- Parallel wavenumber.
k = math.sqrt((k_par * k_par) + (k_perp * k_perp)) -- Total wavenumber.
theta = -math.atan(k_par / k_perp) -- Angle between wave and x-direction.

B0_perp = (1.0 / 3.0) * B0 -- Reference perpendicular magnetic field strength.
u0_perp = (1.0 / 3.0) * vAi -- Reference perpendicular fluid velocity.

nu_elc = 0.01 * omega_ci -- Electron collision frequency.
nu_ion = 0.01 * omega_ci / math.sqrt(mass_ion) -- Ion collision frequency.

-- Simulation parameters.
Nx = 16 -- Cell count (configuration space: x-direction).
Nvx = 32 -- Cell count (velocity space: vx-direction).
Lx = 2.0 * pi / k -- Domain size (configuration space: x-direction).
vx_max_elc = 6.0 * vte -- Domain boundary (electron velocity space: vx-direction).
vx_max_ion = 6.0 * vti -- Domain boundary (ion velocity space: vx-direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 100.0 / omega_ci -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

-- Training parameters.
train_nn = false -- Train neural network on simulation data?
train_ab_initio = true -- Train neural network ab initio?
nn_width = 256 -- Number of neurons to use per layer.
nn_depth = 5 -- Number of layers to use.
train_nn_file = "rt_pkpm_alf_wave_explicit_1x_p1_moms_nn_1" -- File path of neural network to train.
num_trains = GKYL_MAX_INT -- Number of times to train neural network.
num_nn_writes = 1 -- Number of times to write out neural network.
input_moms = { 1, 3, 4 } -- Array of "input" moments to train on.
output_moms = { 5, 6 } -- Array of "output" moments to train on.
test_nn = false -- Test neural network on simulation data?
test_nn_file = "rt_pkpm_alf_wave_explicit_1x_p1_moms_nn_1" -- File path of neural network to test.
num_tests = 1 -- Number of times to test neural network.

pkpmApp = PKPM.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedL2fCalcs = integrated_L2_f_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { 0.0 },
  upper = { Lx },
  cells = { Nx },
  cflFrac = cfl_frac,

  useExplicitSource = true,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1 }, -- Periodic directions (x-direction only).

  -- Electrons.
  elc = PKPM.Species.new {
    modelID = G0.Model.Default,
    charge = charge_elc, mass = mass_elc,
    
    -- Velocity space grid.
    lower = { -vx_max_elc },
    upper = { vx_max_elc },
    cells = { Nvx },

    -- Initial conditions (distribution function).
    initDist = function (t, xn)
      local vx = xn[2]

      local F0 = (n0 / math.sqrt(2.0 * pi * vte * vte)) * (math.exp(-(vx * vx) / (2.0 * vte * vte))) -- Electron distribution function (F0).
      local G = (vte * vte) * F0 -- Electron distribution function (G).
      
      return F0, G
    end,

    -- Initial conditions (fluid).
    initFluid = function (t, xn)
      local x = xn[1]

      local Jy = -B0_perp * k * math.sin(k * x) / mu0

      local vx_drift = 0.0 -- Electron drift velocity (x-direction).
      local vy_drift = -Jy / charge_ion -- Electron drift velocity (y-direction).
      local vz_drift = -u0_perp * math.cos(k * x) -- Electron drift velocity (z-direction).
    
      local mom_x = mass_elc * vx_drift -- Electron total momentum density (x-direction).
      local mom_y = mass_elc * vy_drift -- Electron total momentum density (y-direction).
      local mom_z = mass_elc * vz_drift -- Electron total momentum density (z-direction).

      return mom_x, mom_y, mom_z
    end,

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu_elc -- Electron collision frequency.
      end
    },

    evolve = true -- Evolve species?
  },

  -- Ions.
  ion = PKPM.Species.new {
    modelID = G0.Model.Default,
    charge = charge_ion, mass = mass_ion,
    
    -- Velocity space grid.
    lower = { -vx_max_ion },
    upper = { vx_max_ion },
    cells = { Nvx },

    -- Initial conditions (distribution function).
    initDist = function (t, xn)
      local vx = xn[2]

      local F0 = (n0 / math.sqrt(2.0 * pi * vti * vti)) * (math.exp(-(vx * vx) / (2.0 * vti * vti))) -- Ion distribution function (F0).
      local G = (vti * vti) * F0 -- Ion distribution function (G).
      
      return F0, G
    end,

    -- Initial conditions (fluid).
    initFluid = function (t, xn)
      local x = xn[1]

      local vx_drift = 0.0 -- Ion drift velocity (x-direction).
      local vy_drift = 0.0 -- Ion drift velocity (y-direction).
      local vz_drift = -u0_perp * math.cos(k * x) -- Ion drift velocity (z-direction).
    
      local mom_x = mass_ion * vx_drift -- Ion total momentum density (x-direction).
      local mom_y = mass_ion * vy_drift -- Ion total momentum density (y-direction).
      local mom_z = mass_ion * vz_drift -- Ion total momentum density (z-direction).

      return mom_x, mom_y, mom_z
    end,

    collisions = {
      collisionID = G0.Collisions.LBO,

      selfNu = function (t, xn)
        return nu_ion -- Ion collision frequency.
      end
    },

    evolve = true -- Evolve species?
  },

  -- Field.
  field = PKPM.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]
      
      local Jy = -B0_perp * k * math.sin(k * x) / mu0

      local Bx = B0 * math.sin(theta) -- Total magnetic field (x-direction).
      local By = -B0 * math.cos(theta) -- Total magnetic field (y-direction).
      local Bz = -B0_perp * math.cos(k * x) -- Total magnetic field (z-direction).
    
      local vx_drift = 0.0 -- Electron drift velocity (x-direction).
      local vy_drift = -Jy / charge_ion -- Electron drift velocity (y-direction).
      local vz_drift = -u0_perp * math.cos(k * x) -- Electron drift velocity (z-direction).
    
      local Ex = -((vy_drift * Bz) - (vz_drift * By)) -- Total electric field (x-direction).
      local Ey = -((vz_drift * Bx) - (vx_drift * Bz)) -- Total electric field (y-direction).
      local Ez = -((vx_drift * By) - (vy_drift * Bx)) -- Total electric field (z-direction).

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

    evolve = true, -- Evolve field?
    
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  },

  -- Training parameters.
  trainNN = train_nn,
  trainAbInitio = train_ab_initio,
  NNWidth = nn_width,
  NNDepth = nn_depth,
  trainNNFile = train_nn_file,
  numTrains = num_trains,
  numNNWrites = num_nn_writes,
  inputMoms = input_moms,
  outputMoms = output_moms,
  testNN = test_nn,
  testNNFile = test_nn_file,
  numTests = num_tests
}

pkpmApp:run()