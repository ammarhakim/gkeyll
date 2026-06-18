-- Relativistic pair-plasma force-free current sheet (GEM-style reconnection),
-- using the MODULAR GR-Euler-tetrad fluids coupled to a standard perfectly
-- hyperbolic Maxwell field through the explicit relativistic EM-coupling source
-- solver (SI units, no 4*pi, no D = eps0 E constitutive relation).
--
-- This is the modular replacement for the packed GRTwoFluid setup: instead of a
-- single species bundling electrons + positrons + Maxwell + geometry into one
-- 84-component state, we use
--   * two GREulerTetrad fluids (electrons, positrons), each a 5-component
--     hydro state [D, S_x, S_y, S_z, tau] in Convention A (covariant momentum),
--   * a separate 8-component Maxwell field [E, B, phi, psi] (E and B as-is),
--   * a Minkowski spacetime supplied as a separate Moments.Spacetime component.
-- The fluid <-> field coupling is done by the explicit SSP-RK3 GR-EM source
-- solver, enabled with field.useGRExplicitEmCoupling = true.

local Moments = G0.Moments
local GREulerTetrad = G0.Moments.Eq.GREulerTetrad
local Minkowski = G0.Moments.Spacetime.Minkowski

-- Mathematical constants (dimensionless).
pi = math.pi

local function noise_init(noise_amp, noise_index, k_init, k_final, Lx, Ly, x, y)
  math.randomseed(0)
  local kindex = (noise_index + 1.0) / 2.0
  local B_amp = 0.0
  local B_phase = 0.0
  local Bx_noise = 0.0
  local By_noise = 0.0
  local Jz_noise = 0.0
  for i = k_init, k_final do
    B_amp = math.random()
    B_phase = math.random()

    Bx_noise = Bx_noise - 2.0*(2.0*pi/Ly)*(Lx/(i*2.0*pi))*B_amp*math.sin(2.0*pi*y/Ly)*(math.cos(2.0*pi*y/Ly)+1)*math.cos(i*2.0*pi*x/Lx +  2.0*pi*B_phase)*i^kindex
    By_noise = By_noise + B_amp*(math.cos(2.0*pi*y/Ly) + 1.0)*(math.cos(2.0*pi*y/Ly) + 1.0)*math.sin(i*2.0*pi*x/Lx + 2.0*pi*B_phase)*i^kindex
    Jz_noise = Jz_noise + (2.0*pi*i/Lx)*B_amp*(math.cos(2.0*pi*y/Ly) + 1.0)*(math.cos(2.0*pi*y/Ly) + 1.0)*math.cos(i*2.0*pi*x/Lx + 2*pi*B_phase)*i^kindex +
                 2.0*(2.0*pi/Ly)*(2.0*pi/Ly)*(Lx/(i*2.0*pi))*B_amp*(math.sin(2.0*pi*y/Ly)*math.sin(2.0*pi*y/Ly) - math.cos(2.0*pi*y/Ly)*(math.cos(2.0*pi*y/Ly)+1.0))*math.cos(i*2.0*pi*x/Lx +  2.0*pi*B_phase)*i^kindex
  end

  local kdiff = k_final - k_init + 1.0;
  Bx_noise = noise_amp*Bx_noise/math.sqrt(2.0*kdiff*kdiff/3.0)
  By_noise = noise_amp*By_noise/math.sqrt(2.0*kdiff*kdiff/3.0)
  Jz_noise = noise_amp*Jz_noise/math.sqrt(2.0*kdiff*kdiff/3.0)

  return Bx_noise, By_noise, Jz_noise
end

-- Physical constants (using normalized code units).
gas_gamma_elc = 5.0/3.0
gas_gamma_ion = 5.0/3.0
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
light_speed = 1.0/math.sqrt(epsilon0*mu0) -- Speed of light.
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.
mass_ion = 1.0 -- Positron mass.
charge_ion = 1.0 -- Positron charge.

sigma = 10.0 -- B^2/(mu0*n0*m*c^2)
n0 = 1.0 -- reference density
T0 = 0.01 -- Reference temperature in units of mc^2

-- Derived parameters
wpe = math.sqrt(charge_ion^2*n0/(epsilon0*mass_elc))
B0 = math.sqrt(sigma*mu0*n0*mass_elc*light_speed^2) -- in-plane magnetic field strength
vA = math.sqrt(sigma/(sigma + 1.0))*light_speed -- in-plane Alfven velocity
omegaCi = charge_ion*B0/mass_ion
de = light_speed/wpe

-- Reconnection parameters
guide = 0.01 -- Guide-field strength.
w0 = de -- Layer width.
psi0 = 0.1*B0*de -- Layer perturbation.
noise_amp = 0.001*B0 -- Noise amplitude.
k_init = 1 -- First wave mode to perturb with noise, 1.0 correspond to box size.
k_final = 20 -- Last wave mode to perturb with noise.
noise_index = -1.0 -- Spectral index of the noise.

-- Simulation parameters.
Nx = 512 -- Cell count (configuration space: x-direction).
Ny = 512 -- Cell count (configuration space: y-direction).
Lx = 8.0*pi*de  -- Domain size (configuration space: x-direction).
Ly = 8.0*pi*de  -- Domain size (configuration space: y-direction).
cfl_frac = 0.95 -- CFL coefficient.

reinit_freq = 100 -- Spacetime reinitialization frequency.

t_end = 100.0 -- Final simulation time.
num_frames = 10 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

-- Build the modular GR-Euler-tetrad conserved hydro state for one species with
-- rest density rho, pressure p, spatial four-velocity u^i (= W v^i), adiabatic
-- index gas_gamma, in Convention A: D = sqrt(gamma) rho W,
-- S_i = sqrt(gamma) rho h W^2 v_i, tau = sqrt(gamma)(rho h W^2 - p - rho W).
local function gr_euler_state(rho, p, gas_gamma, ux, uy, uz, spatial_metric, spatial_det)
  -- W from the spatial four-velocity: W = sqrt(1 + g_ij u^i u^j).
  local u = { ux, uy, uz }
  local u_sq = 0.0
  for i = 1, 3 do
    for j = 1, 3 do
      u_sq = u_sq + spatial_metric[i][j] * u[i] * u[j]
    end
  end
  local W = math.sqrt(1.0 + u_sq)

  local vel = { ux / W, uy / W, uz / W }

  local h = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0)))

  local v_lower = { 0.0, 0.0, 0.0 }
  for i = 1, 3 do
    for j = 1, 3 do
      v_lower[i] = v_lower[i] + spatial_metric[i][j] * vel[j]
    end
  end

  local D = math.sqrt(spatial_det) * rho * W
  local mom_x = math.sqrt(spatial_det) * rho * h * (W * W) * v_lower[1]
  local mom_y = math.sqrt(spatial_det) * rho * h * (W * W) * v_lower[2]
  local mom_z = math.sqrt(spatial_det) * rho * h * (W * W) * v_lower[3]
  local Etot = math.sqrt(spatial_det) * ((rho * h * (W * W)) - p - (rho * W))

  return D, mom_x, mom_y, mom_z, Etot
end

-- Per-species spatial four-velocity that carries (1/2) of the equilibrium
-- current J in the physically correct direction: contribution_s = q_s n u_s,
-- so u_s = 0.5 J / (q_s n). With a pair plasma this gives electrons and
-- positrons counter-streaming, each carrying half the total current.
local function species_u(charge, Jx, Jy, Jz)
  return (0.5 * Jx) / (charge * n0), (0.5 * Jy) / (charge * n0), (0.5 * Jz) / (charge * n0)
end

momentApp = Moments.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { -0.5 * Lx, -0.5 * Ly },
  upper = { 0.5 * Lx, 0.5 * Ly },
  cells = { Nx, Ny },
  cflFrac = cfl_frac,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coordinate direction (x- and y-directions).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1 }, -- Periodic directions (x-direction).

  -- Minkowski spacetime supplied as a separate (static) component; the
  -- GR-Euler-tetrad fluids acquire the geometry through it.
  spacetime = Moments.Spacetime.new {
    evolve = false,
    equation = Minkowski.new(),
    spacetimeGauge = G0.SpacetimeGauge.Static,
    reinitFreq = reinit_freq,
  },

  -- Electron fluid (modular GR Euler, tetrad basis).
  elc = Moments.Species.new {
    charge = charge_elc, mass = mass_elc,

    equation = GREulerTetrad.new {
      rpType = "hllc",
      eos = "approx_synge",
      useRcc = true
    },

    init = function (t, xn)
      local x, y = xn[1], xn[2]
      local sech_sq = (1.0 / math.cosh(y / w0))^2

      local _, _, Jz_noise = noise_init(noise_amp, noise_index, k_init, k_final, Lx, Ly, x, y)

      local Jx = -B0/w0*math.tanh(y/w0)*sech_sq/(math.sqrt(sech_sq + guide*guide));
      local Jy = 0.0;
      local Jz = B0/w0*sech_sq - psi0*math.cos(2.0*pi*x/Lx)*math.cos(pi*y/Ly)*((2.0*pi/Lx)^2 + (pi/Ly)^2) + Jz_noise;

      local ux, uy, uz = species_u(charge_elc, Jx, Jy, Jz)

      local rho = n0 * mass_elc
      local p = n0 * T0

      local spatial_metric = Minkowski.spatialMetricTensor(0.0, x, y, 0.0)
      local spatial_det = Minkowski.spatialMetricDeterminant(0.0, x, y, 0.0)

      return gr_euler_state(rho, p, gas_gamma_elc, ux, uy, uz, spatial_metric, spatial_det)
    end,

    evolve = true,
    bcy = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy }, -- Copy boundary conditions (y-direction).
  },

  -- Positron fluid (modular GR Euler, tetrad basis).
  ion = Moments.Species.new {
    charge = charge_ion, mass = mass_ion,

    equation = GREulerTetrad.new {
      rpType = "hllc",
      eos = "approx_synge",
      useRcc = true
    },

    init = function (t, xn)
      local x, y = xn[1], xn[2]
      local sech_sq = (1.0 / math.cosh(y / w0))^2

      local _, _, Jz_noise = noise_init(noise_amp, noise_index, k_init, k_final, Lx, Ly, x, y)

      local Jx = -B0/w0*math.tanh(y/w0)*sech_sq/(math.sqrt(sech_sq + guide*guide));
      local Jy = 0.0;
      local Jz = B0/w0*sech_sq - psi0*math.cos(2.0*pi*x/Lx)*math.cos(pi*y/Ly)*((2.0*pi/Lx)^2 + (pi/Ly)^2) + Jz_noise;

      local ux, uy, uz = species_u(charge_ion, Jx, Jy, Jz)

      local rho = n0 * mass_ion
      local p = n0 * T0

      local spatial_metric = Minkowski.spatialMetricTensor(0.0, x, y, 0.0)
      local spatial_det = Minkowski.spatialMetricDeterminant(0.0, x, y, 0.0)

      return gr_euler_state(rho, p, gas_gamma_ion, ux, uy, uz, spatial_metric, spatial_det)
    end,

    evolve = true,
    bcy = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy }, -- Copy boundary conditions (y-direction).
  },

  -- Standard perfectly hyperbolic Maxwell field; E and B are used directly by
  -- the relativistic EM-coupling source solver (no D = eps0 E relation).
  field = Moments.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 1.0,

    -- Enable the explicit relativistic multi-fluid + Maxwell coupling.
    useGRExplicitEmCoupling = true,

    init = function (t, xn)
      local x, y = xn[1], xn[2]
      local sech_sq = (1.0 / math.cosh(y / w0))^2

      local Bx_noise, By_noise, _ = noise_init(noise_amp, noise_index, k_init, k_final, Lx, Ly, x, y)

      local b1x = -B0*math.tanh(y/w0)
      local b1y = 0.0;
      local b1z = B0*math.sqrt(guide*guide + sech_sq);

      local Ex = 0.0
      local Ey = 0.0
      local Ez = 0.0
      local Bx = b1x + psi0 * (pi / Ly) * math.cos(2 * pi * x / Lx) * math.sin(pi * y / Ly) + Bx_noise;
      local By = b1y - psi0 * (2 * pi / Lx) * math.sin(2 * pi * x / Lx) * math.cos(pi * y / Ly) + By_noise;
      local Bz = b1z;

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,

    evolve = true,
    bcy = { G0.FieldBc.bcCopy, G0.FieldBc.bcCopy }, -- Copy boundary conditions (y-direction).
  }
}

-- Run application.
momentApp:run()
