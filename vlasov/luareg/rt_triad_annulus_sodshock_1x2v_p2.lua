local Vlasov = G0.Vlasov

-- Physical constants (using normalized code units).
mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.

nl = 1.0 -- Left/inner number density.
Tl = 1.0 -- Left/inner temperature.
V_r_drift_l = 0.0 -- Left/inner drift velocity (radial direction).
V_theta_drift_l = 0.0 -- Left/inner drift velocity (angular direction).

nr = 0.125 -- Right/outer number density.
Tr = math.sqrt(0.1 / 0.125) -- Right/outer temperature.
V_r_drift_r = 0.0 -- Right/outer drift velocity (radial direction).
V_theta_drift_r = 0.0 -- Right/outer drift velocity (angular direction).

vt = 1.0 -- Thermal velocity.
nu = 15000.0 -- Collision frequency.

-- Simulation parameters.
Nr = 128 -- Cell count (configuration space: radial direction).
Nvr = 12 -- Cell count (velocity space: radial direction).
Nvtheta = 12 -- Cell count (velocity space: angular direction).
Lr = 1.0 -- Domain size (configuration space: radial direction).
vr_max = 8.0 * vt -- Domain boundary (velocity space: radial direction).
vtheta_max = 8.0 * vt -- Domain boundary (velocity space: angular direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 0.1 -- Final simulation time.
num_frames = 20 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

midplane = 1.0 -- Radial midplane location designating jump in quantities.

vlasovApp = Vlasov.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedL2fCalcs = integrated_L2_f_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { 0.5 },
  upper = { 0.5 + Lr },
  cells = { Nr },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).

  -- Neutral species.
  neut = Vlasov.Species.new {
    modelID = G0.Model.Triad,
    charge = charge, mass = mass,

    -- vals_ij = e_{r,\theta} . \sigma_{x,y}
    covTangentBasis = function (t, xn)
      local q_r = xn[1]
      local q_theta = 0.0

      local e_rx = math.cos(q_theta) -- Covariant Tangent Basis Coefficients (r-x coefficient).
      local e_ry = math.sin(q_theta) -- Covariant Tangent Basis Coefficients (r-y coefficient).
      local e_tx = - q_r * math.sin(q_theta) -- Covariant Tangent Basis Coefficients (theta-x coefficient).
      local e_ty = q_r * math.cos(q_theta)-- Covariant Tangent Basis Coefficients (theta-y coefficient).

      return e_rx, e_ry, e_tx, e_ty
    end,

    -- vals_ij = \sigma_{r,\theta} . \sigma_{x,y}
    triadBasis = function (t, xn)
      local q_r = xn[1]
      local q_theta = 0.0

      local sigma_rx = math.cos(q_theta) -- Triad Basis Coefficients (r-x coefficient).
      local sigma_ry = math.sin(q_theta) -- Triad Basis Coefficients (r-y coefficient).
      local sigma_tx = - math.sin(q_theta) -- Triad Basis Coefficients (theta-x coefficient).
      local sigma_ty = math.cos(q_theta)-- Triad Basis Coefficients (theta-y coefficient).

      return sigma_rx, sigma_ry, sigma_tx, sigma_ty
    end,

    -- d(vals_ij)/dx^k = d(\sigma_{r, \theta} . \sigma_{x,y})/dx^k
    -- x^k = { r, \theta }
    triadBasisGradient = function (t, xn)
      local q_r = xn[1]
      local q_theta = 0.0

      -- d/dr components
      local d_sigma_rx_dr = 0.0 -- Triad Basis Gradient Coefficients (r-x coefficient).
      local d_sigma_ry_dr = 0.0 -- Triad Basis Gradient Coefficients (r-y coefficient).
      local d_sigma_tx_dr = 0.0 -- Triad Basis Gradient Coefficients (theta-x coefficient).
      local d_sigma_ty_dr = 0.0-- Triad Basis Gradient Coefficients (theta-y coefficient).

      -- d/dtheta components
      local d_sigma_rx_dtheta = - math.sin(q_theta) -- Triad Basis Gradient Coefficients (r-x coefficient).
      local d_sigma_ry_dtheta = math.cos(q_theta) -- Triad Basis Gradient Coefficients (r-y coefficient).
      local d_sigma_tx_dtheta = - math.cos(q_theta) -- Triad Basis Gradient Coefficients (theta-x coefficient).
      local d_sigma_ty_dtheta = - math.sin(q_theta)-- Triad Basis Gradient Coefficients (theta-y coefficient).

      return d_sigma_rx_dr, d_sigma_ry_dr, d_sigma_tx_dr, d_sigma_ty_dr, d_sigma_rx_dtheta, d_sigma_ry_dtheta, d_sigma_tx_dtheta, d_sigma_ty_dtheta
    end,

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
          local r = xn[1]

          local n = 0.0
          if r < midplane then
            n = nl -- Total number density (left/inner).
          else
            n = nr -- Total number density (right/outer).
          end

          return n
        end,
        temperatureInit = function (t, xn)
          local r = xn[1]

          local T = 0.0
          if r < midplane then
            T = Tl -- Isotropic temperature (left/inner).
          else
            T = Tr -- Isotropic temperature (right/outer).
          end

          return T
        end,
        driftVelocityInit = function (t, xn)
          local r = xn[1]

          local V_r_drift = 0.0
          local V_theta_drift = 0.0

          if r < midplane then
            V_r_drift = V_r_drift_l -- Radial drift velocity (left/inner).
            V_theta_drift = V_theta_drift_l -- Angular drift velocity (left/inner).
          else
            V_r_drift = V_r_drift_r -- Radial drift velocity (right/outer).
            V_theta_drift = V_theta_drift_r --Angular drift velocity (right/outer).
          end

          return V_r_drift, V_theta_drift
        end,

        correctAllMoments = true,
        iterationEpsilon = 0.0,
        maxIterations = 0,
        useLastConverged = false
      }
    },

    collisions = {
      collisionID = G0.Collisions.BGK,

      selfNu = function (t, xn)
        return nu -- Collision frequency.
      end,
      
      useImplicitCollisionScheme = true,
      correctAllMoments = true,
      iterationEpsilon = 0.0,
      maxIterations = 0,
      useLastConverged = false
    },

    bcx = {
      lower = {
        type = G0.SpeciesBc.bcReflect
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