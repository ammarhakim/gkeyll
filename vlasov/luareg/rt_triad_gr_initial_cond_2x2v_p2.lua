local Vlasov = G0.Vlasov

-- Physical constants (using normalized code units).
mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.

nl = 1.0 -- Left/inner number density.
Tl = 1.0 -- Left/inner temperature.
V_r_drift_l = 0.25 -- Left/inner drift velocity (radial direction).
V_theta_drift_l = 0.25 -- Left/inner drift velocity (angular direction).

vt = 1.0 -- Thermal velocity.
nu = 15000.0 -- Collision frequency.

-- Simulation parameters.
Nr = 24 -- Cell count (configuration space: radial direction).
Ntheta = 24 -- Cell count (configuration space: radial direction).
Nvr = 12 -- Cell count (velocity space: radial direction).
Nvtheta = 12 -- Cell count (velocity space: angular direction).
Lr = 8.7 -- Domain size (configuration space: radial direction).
Ltheta = 2 * math.pi -- Domain size (configuration space: radial direction).
vr_max = 8 * vt -- Domain boundary (velocity space: radial direction).
vtheta_max = 8 * vt -- Domain boundary (velocity space: angular direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 0.1 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

massBH = 1.0 -- Mass of the Kerr-Schild black hole
spinBH = 0.8 -- Spin parameter, a = J/M, of the Kerr-Schild black hole

vlasovApp = Vlasov.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedL2fCalcs = integrated_L2_f_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { 0.3 , 0.0 },
  upper = { 0.3 + Lr, Ltheta },
  cells = { Nr, Ntheta },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 2 }, -- Periodic directions (phi-dir).

  -- Vlasov Geometry
  geom = Vlasov.Geom.new {
    
    -- Black hole parameters
    massBH = massBH,
    spinBH = spinBH
    
  },

  -- Neutral species.
  neut = Vlasov.Species.new {
    modelID = G0.Model.TriadGR,
    charge = charge, mass = mass,

    -- Use the vierbein inputs to contruct the Poisson Tensor
    useLo = false,
    usePresetGeom = true,
    triadPresetGeomType = G0.TriadGeom.GR_KS_rphi,

    -- Velocity space grid.
    lower = { -vr_max, -vtheta_max },
    upper = { vr_max, vtheta_max },
    cells = { Nvr, Nvtheta },

  -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)

          local n = 0.0
          n = nl -- Total number density (left/inner).

          local metric_det = 1.0

          return metric_det * n
        end,
        temperatureInit = function (t, xn)

          local T = 0.0
          T = Tl -- Isotropic temperature (left/inner)

          return T
        end,
        driftVelocityInit = function (t, xn)

          local V_r_drift = 0.0
          local V_r_drift = 0.0

          V_r_drift = V_r_drift_l -- x drift velocity (left/inner).
          V_r_drift = V_r_drift_l -- y drift velocity (left/inner).

          return V_r_drift, V_r_drift
        end,

        correctAllMoments = true,
        iterationEpsilon = 0.0,
        maxIterations = 0,
        useLastConverged = false
      }
    },

    bcx = {
      lower = {
        type = G0.SpeciesBc.bcAbsorb
      },
      upper = {
        type = G0.SpeciesBc.bcReflect
      }
    },

    evolve = true, -- Evolve species?
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.LTEMoments, G0.Moment.EnergyMoment }
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