local Vlasov = G0.Vlasov
-- Check the case with no blackhole as a diagnostic for the 1x3v neutral relativistic Sod shock test.

-- Physical constants (using normalized code units).
mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.

nl = 1.0 -- Left number density.
Tl = 0.1 -- Left temperature.
V_r_drift_l = 0.0 -- Left drift velocity (x direction).
V_y_drift_l = 0.0 -- Left drift velocity (y direction).
V_z_drift_l = 0.0 -- Left drift velocity (z direction).

nr = 0.1 -- Right number density.
Tr = 0.08 -- Right temperature.
V_r_drift_r = 0.0 -- Right drift velocity (x direction).
V_y_drift_r = 0.0 -- Right drift velocity (y direction).
V_z_drift_r = 0.0 -- Right drift velocity (z direction).

vt = 1.0 -- Thermal velocity.
nu = 15000.0 -- Collision frequency.

-- Simulation parameters.
Nx = 64 --64 -- Cell count (configuration space: x direction).
Nvx = 16 -- 64 -- Cell count (velocity space: x direction).
Nvy =  16 -- 64 --Cell count (velocity space: y direction).
Nvz = 16 -- 64 --Cell count (velocity space: z direction).
Lx = 1.0 -- Domain size (configuration space: x direction).
vx_max = 2.0 * vt -- Domain boundary (velocity space: x direction).
vy_max = 2.0 * vt -- Domain boundary (velocity space: y direction).
vz_max = 2.0 * vt -- Domain boundary (velocity space: z direction).
poly_order = 1 -- Polynomial order.
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

-- Check the case with no blackhole as a diagnostic
massBH = 0.0 -- Mass of the Kerr-Schild black hole
spinBH = 0.0 -- Spin parameter, a = J/M, of the Kerr-Schild black hole

vlasovApp = Vlasov.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedL2fCalcs = integrated_L2_f_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { 1.0 },
  upper = { 1.0 + Lx },
  cells = { Nx },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (z-dir).

  -- Vlasov Geometry
  geom = Vlasov.Geom.new {
    
    -- Black hole parameters
    massBH = massBH,
    spinBH = spinBH,

    -- Use the vierbein inputs to contruct the Poisson Tensor
    usePresetGeom = true,
    triadPresetGeomType = G0.TriadGeom.GR_KS_Cart_3V

  },

  -- Neutral species.
  neut = Vlasov.Species.new {
    modelID = G0.Model.TriadGR,
    charge = charge, mass = mass,

    -- Velocity space grid.
    lower = { -vx_max, -vy_max , -vz_max },
    upper = { vx_max, vy_max, vz_max },
    cells = { Nvx, Nvy, Nvz },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local x = xn[1]

          local n = 0.0
          if x < 1.5 then
            n = nl -- Total number density (left).
          else
            n = nr -- Total number density (right).
          end

          return n
        end,
        temperatureInit = function (t, xn)
          local x = xn[1]

          local T = 0.0
          if x < 1.5 then
            T = Tl -- Isotropic temperature (left).
          else
            T = Tr -- Isotropic temperature (right).
          end

          return T
        end,
        driftVelocityInit = function (t, xn)
          local x = xn[1]

          local Vx_drift = 0.0

          if x < 1.5 then
            Vx_drift = Vx_drift_l -- Total drift velocity (left).
            Vy_drift = Vy_drift_l -- Total drift velocity (left).
            Vz_drift = Vz_drift_l -- Total drift velocity (left).
          else
            Vx_drift = Vx_drift_r -- Total drift velocity (right).
            Vy_drift = Vy_drift_r -- Total drift velocity (right).
            Vz_drift = Vz_drift_r -- Total drift velocity (right).
          end

          return Vx_drift, Vy_drift, Vz_drift
        end,


        correctAllMoments = true,
        iterationEpsilon = 0.0,
        maxIterations = 0,
        useLastConverged = false
      }
    },

    bcx = {
      lower = {
        type = G0.SpeciesBc.bcCopy
      },
      upper = {
        type = G0.SpeciesBc.bcCopy
      }
    },

    collisions = {
      collisionID = G0.Collisions.BGK,

      selfNu = function (t, xn)
        return nu -- Collision frequency.
      end,
      
      useImplicitCollisionScheme = true
    },

    correct = {
      correctAllMoments = true,
      iterationEpsilon = 1e-12,
      maxIterations = 100,
      useLastConverged = false
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
