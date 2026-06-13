-- Lua version of vlasov/creg/rt_vlasov_lbo_cross_relax_1x3v_p2.c (the physical analog of
-- the gyrokinetic 1x2v test rt_gk_lbo_cross_relax_1x2v_p1): a Vlasov 1x3v cross-species
-- LBO relaxation with anisotropic bi-Maxwellians (parallel = v_z, perpendicular = v_x,v_y).

local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using non-normalized physical units; values match gkyl_const.h).
epsilon0 = 8.854187817620389850536563031710750260608e-12 -- Permittivity of free space.
planck_h = 6.62606896e-34 -- Planck's constant.
hbar = planck_h / (2.0 * pi) -- Planck's constant / 2 pi.
mass_elc = 9.10938215e-31 -- Electron mass.
mass_ion = 2.014 * 1.672621637e-27 -- Ion (deuterium) mass.
elem_charge = 1.602176487e-19 -- Elementary charge.
charge_elc = -elem_charge -- Electron charge.
charge_ion = elem_charge -- Ion charge.

T_par_elc = 300.0 * elem_charge -- Parallel electron temperature.
T_par_ion = 200.0 * elem_charge -- Parallel ion temperature.
alpha = 1.3 -- Ratio of perpendicular to parallel temperatures.

n0 = 7.0e19 -- Reference number density (1 / m^3).
B0 = 1.0 -- Reference magnetic field strength (Tesla).

-- Derived physical quantities (using non-normalized physical units).
T_perp_elc = alpha * T_par_elc -- Perpendicular electron temperature.
T_perp_ion = alpha * T_par_ion -- Perpendicular ion temperature.

Te = (T_par_elc + (2.0 * T_perp_elc)) / 3.0 -- Electron temperature.
Ti = (T_par_ion + (2.0 * T_perp_ion)) / 3.0 -- Ion temperature.

log_lambda_elc = 6.6 - 0.5 * math.log(n0 / 1.0e20) + 1.5 * math.log(T_par_elc / charge_ion) -- Electron Coulomb logarithm.
log_lambda_ion = 6.6 - 0.5 * math.log(n0 / 1.0e20) + 1.5 * math.log(T_par_ion / charge_ion) -- Ion Coulomb logarithm.
nu_elc = log_lambda_elc * (charge_ion ^ 4.0) * n0 /
  (6.0 * math.sqrt(2.0) * (pi ^ (3.0 / 2.0)) * (epsilon0 ^ 2.0) * math.sqrt(mass_elc) * (T_par_elc ^ (3.0 / 2.0))) -- Electron collision frequency.
nu_ion = log_lambda_ion * (charge_ion ^ 4.0) * n0 /
  (12.0 * (pi ^ (3.0 / 2.0)) * (epsilon0 ^ 2.0) * math.sqrt(mass_ion) * (T_par_ion ^ (3.0 / 2.0))) -- Ion collision frequency.

vte_par = math.sqrt(T_par_elc / mass_elc) -- Parallel electron thermal velocity.
vti_par = math.sqrt(T_par_ion / mass_ion) -- Parallel ion thermal velocity.
vte = math.sqrt(Te / mass_elc) -- Electron thermal velocity.
vti = math.sqrt(Ti / mass_ion) -- Ion thermal velocity.

upar_elc = 0.5 * math.sqrt(mass_elc / mass_ion) * vte -- Parallel electron velocity.
upar_ion = 50.0 * (mass_elc / mass_ion) * vti -- Parallel ion velocity.

-- Simulation parameters.
Nz = 1 -- Cell count (configuration space: z-direction).
Nvx = 8 -- Cell count (velocity space: vx-direction).
Nvy = 8 -- Cell count (velocity space: vy-direction).
Nvz = 8 -- Cell count (velocity space: vz-direction).
Lz = 4.0 -- Domain size (configuration space: z-direction).
vx_max_elc = 5.0 * vte_par -- Electron domain boundary (velocity space: perpendicular vx, vy).
vz_max_elc = 5.0 * vte_par -- Electron domain boundary (velocity space: parallel vz).
vx_max_ion = 5.0 * vti_par -- Ion domain boundary (velocity space: perpendicular vx, vy).
vz_max_ion = 5.0 * vti_par -- Ion domain boundary (velocity space: parallel vz).
poly_order = 2 -- Polynomial order (Vlasov LBO requires p=2; GK uses its p=1 hybrid basis).
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 0.1 / nu_ion -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

vlasovApp = Vlasov.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedL2fCalcs = integrated_L2_f_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { -0.5 * Lz },
  upper = { 0.5 * Lz },
  cells = { Nz },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (z-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1 }, -- Periodic directions (z-direction only).

  -- Electrons.
  elc = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_elc, mass = mass_elc,

    -- Velocity space grid (perpendicular = vx, vy; parallel = vz).
    lower = { -vx_max_elc, -vx_max_elc, -vz_max_elc },
    upper = { vx_max_elc, vx_max_elc, vz_max_elc },
    cells = { Nvx, Nvy, Nvz },

    -- 0D-in-configuration collision relaxation: streaming and velocity-space EM
    -- advection are identically zero, so skip the collisionless terms entirely.
    noCollisionlessTerms = true,

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,

        init = function (t, xn)
          local vx, vy, vz = xn[2], xn[3], xn[4]

          local vtpar_sq = T_par_elc / mass_elc
          local vtperp_sq = T_perp_elc / mass_elc

          -- Anisotropic bi-Maxwellian: parallel = vz (drift u_par, T_par), perp = (vx, vy) (T_perp).
          local norm = n0 / ((2.0 * pi) ^ (3.0 / 2.0) * vtperp_sq * math.sqrt(vtpar_sq))
          local n = norm * math.exp(-(vx * vx + vy * vy) / (2.0 * vtperp_sq)
            - (vz - upar_elc) * (vz - upar_elc) / (2.0 * vtpar_sq)) -- Distribution function.

          return n
        end
      }
    },

    collisions = {
      collisionID = G0.Collisions.LBO,

      -- No selfNu provided: use the Spitzer (normNu) collisionality computed from the
      -- reference density and temperature, matching the gyrokinetic reference.
      n0 = n0, -- Reference density.
      T0 = Te, -- Reference temperature.
      epsilon0 = epsilon0,
      hbar = hbar,
      eV = elem_charge,

      numCrossCollisions = 1,
      collideWith = { "ion" }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2, G0.Moment.M2ij }
  },

  -- Ions.
  ion = Vlasov.Species.new {
    modelID = G0.Model.Default,
    charge = charge_ion, mass = mass_ion,

    -- Velocity space grid (perpendicular = vx, vy; parallel = vz).
    lower = { -vx_max_ion, -vx_max_ion, -vz_max_ion },
    upper = { vx_max_ion, vx_max_ion, vz_max_ion },
    cells = { Nvx, Nvy, Nvz },

    noCollisionlessTerms = true,

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,

        init = function (t, xn)
          local vx, vy, vz = xn[2], xn[3], xn[4]

          local vtpar_sq = T_par_ion / mass_ion
          local vtperp_sq = T_perp_ion / mass_ion

          -- Anisotropic bi-Maxwellian: parallel = vz (drift u_par, T_par), perp = (vx, vy) (T_perp).
          local norm = n0 / ((2.0 * pi) ^ (3.0 / 2.0) * vtperp_sq * math.sqrt(vtpar_sq))
          local n = norm * math.exp(-(vx * vx + vy * vy) / (2.0 * vtperp_sq)
            - (vz - upar_ion) * (vz - upar_ion) / (2.0 * vtpar_sq)) -- Distribution function.

          return n
        end
      }
    },

    collisions = {
      collisionID = G0.Collisions.LBO,

      n0 = n0, -- Reference density.
      T0 = Ti, -- Reference temperature.
      epsilon0 = epsilon0,
      hbar = hbar,
      eV = elem_charge,

      numCrossCollisions = 1,
      collideWith = { "elc" }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.M2, G0.Moment.M2ij }
  },

  skipField = true,

  -- Field.
  field = Vlasov.Field.new {
    epsilon0 = 1.0, mu0 = 1.0,

    -- Initial conditions function.
    init = function (t, xn)
      return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    end,

    evolve = false, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0,

    isStatic = true
  }
}

vlasovApp:run()
