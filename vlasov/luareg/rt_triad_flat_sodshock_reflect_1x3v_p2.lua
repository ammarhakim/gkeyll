local Vlasov = G0.Vlasov

-- Physical constants (using normalized code units).
mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.

nl = 1.0 -- Left/inner number density.
Tl = 1.0 -- Left/inner temperature.
V_x_drift_l = 0.0 -- Left/inner drift velocity (x direction).
V_y_drift_l = 0.0 -- Left/inner drift velocity (y direction).
V_z_drift_l = 0.0 -- Left/inner drift velocity (y direction).

nx = 0.125 -- Right/outer number density.
Tr = math.sqrt(0.1 / 0.125) -- Right/outer temperature.
V_x_drift_r = 0.0 -- Right/outer drift velocity (x direction).
V_y_drift_r = 0.0 -- Right/outer drift velocity (y direction).
V_z_drift_r = 0.0 -- Right/outer drift velocity (y direction).

vt = 1.0 -- Thermal velocity.
nu = 15000.0 -- Collision frequency.

-- Simulation parameters.
Nx = 48 -- Cell count (configuration space: x direction).
Nvx = 8 -- Cell count (velocity space: x direction).
Nvy = 8 -- Cell count (velocity space: y direction).
Nvz = 8 -- Cell count (velocity space: y direction).
Lx = 1.0 -- Domain size (configuration space: x direction).
vx_max = 8.0 * vt -- Domain boundary (velocity space: x direction).
vy_max = 8.0 * vt -- Domain boundary (velocity space: y direction).
vz_max = 8.0 * vt -- Domain boundary (velocity space: y direction).
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

midplane = 0.5 -- x midplane location designating jump in quantities.

vlasovApp = Vlasov.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedL2fCalcs = integrated_L2_f_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { 0 },
  upper = { Lx },
  cells = { Nx },
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

    -- vals_ij = e_{x,y} . \sigma_{x,y}
    covTangentBasis = function (t, xn)

      local e_xx = 1.0 -- Covariant Tangent Basis Coefficients (x-x coefficient).
      local e_xy = 0.0 -- Covariant Tangent Basis Coefficients (x-y coefficient).
      local e_xz = 0.0 -- Covariant Tangent Basis Coefficients (x-z coefficient).
      local e_yx = 0.0 -- Covariant Tangent Basis Coefficients (y-x coefficient).
      local e_yy = 1.0 -- Covariant Tangent Basis Coefficients (y-y coefficient).
      local e_yz = 0.0 -- Covariant Tangent Basis Coefficients (y-z coefficient).
      local e_zx = 0.0 -- Covariant Tangent Basis Coefficients (z-x coefficient).
      local e_zy = 0.0 -- Covariant Tangent Basis Coefficients (z-y coefficient).
      local e_zz = 1.0 -- Covariant Tangent Basis Coefficients (z-z coefficient).

      return e_xx, e_xy, e_xz, e_yx, e_yy, e_yz, e_zx, e_zy, e_zz
    end,

    -- vals_ij = \sigma_{x, y} . \sigma_{x,y}
    triadBasis = function (t, xn)

      local sigma_xx = 1.0 -- Triad Basis Coefficients (x-x coefficient).
      local sigma_xy = 0.0 -- Triad Basis Coefficients (x-y coefficient).
      local sigma_xz = 0.0 -- Triad Basis Coefficients (x-z coefficient).
      local sigma_yx = 0.0 -- Triad Basis Coefficients (y-x coefficient).
      local sigma_yy = 1.0 -- Triad Basis Coefficients (y-y coefficient).
      local sigma_yz = 0.0 -- Triad Basis Coefficients (y-z coefficient).
      local sigma_zx = 0.0 -- Triad Basis Coefficients (z-x coefficient).
      local sigma_zy = 0.0 -- Triad Basis Coefficients (z-y coefficient).
      local sigma_zz = 1.0 -- Triad Basis Coefficients (z-z coefficient).

      return sigma_xx, sigma_xy, sigma_xz, sigma_yx, sigma_yy, sigma_yz, sigma_zx, sigma_zy, sigma_zz
    end,

    -- d(vals_ij)/dx^k = d(\sigma_{x, y} . \sigma_{x,y})/dx^k
    -- x^k = { x, y }
    triadBasisGradient = function (t, xn)

      -- d/dx components
      local d_sigma_xx_dx = 0.0 -- Triad Basis Gradient Coefficients (x-x coefficient).
      local d_sigma_xy_dx = 0.0 -- Triad Basis Gradient Coefficients (x-y coefficient).
      local d_sigma_xz_dx = 0.0 -- Triad Basis Gradient Coefficients (x-z coefficient).
      local d_sigma_yx_dx = 0.0 -- Triad Basis Gradient Coefficients (y-x coefficient).
      local d_sigma_yy_dx = 0.0 -- Triad Basis Gradient Coefficients (y-y coefficient).
      local d_sigma_yz_dx = 0.0 -- Triad Basis Gradient Coefficients (y-z coefficient).
      local d_sigma_zx_dx = 0.0 -- Triad Basis Gradient Coefficients (z-x coefficient).
      local d_sigma_zy_dx = 0.0 -- Triad Basis Gradient Coefficients (z-y coefficient).
      local d_sigma_zz_dx = 0.0 -- Triad Basis Gradient Coefficients (z-z coefficient).

      -- d/dy components
      local d_sigma_xx_dy = 0.0 -- Triad Basis Gradient Coefficients (x-x coefficient).
      local d_sigma_xy_dy = 0.0 -- Triad Basis Gradient Coefficients (x-y coefficient).
      local d_sigma_xz_dy = 0.0 -- Triad Basis Gradient Coefficients (x-z coefficient).
      local d_sigma_yx_dy = 0.0 -- Triad Basis Gradient Coefficients (y-x coefficient).
      local d_sigma_yy_dy = 0.0 -- Triad Basis Gradient Coefficients (y-y coefficient).
      local d_sigma_yz_dy = 0.0 -- Triad Basis Gradient Coefficients (y-z coefficient).
      local d_sigma_zx_dy = 0.0 -- Triad Basis Gradient Coefficients (z-x coefficient).
      local d_sigma_zy_dy = 0.0 -- Triad Basis Gradient Coefficients (z-y coefficient).
      local d_sigma_zz_dy = 0.0 -- Triad Basis Gradient Coefficients (z-z coefficient).

      -- d/dz components
      local d_sigma_xx_dz = 0.0 -- Triad Basis Gradient Coefficients (x-x coefficient).
      local d_sigma_xy_dz = 0.0 -- Triad Basis Gradient Coefficients (x-y coefficient).
      local d_sigma_xz_dz = 0.0 -- Triad Basis Gradient Coefficients (x-z coefficient).
      local d_sigma_yx_dz = 0.0 -- Triad Basis Gradient Coefficients (y-x coefficient).
      local d_sigma_yy_dz = 0.0 -- Triad Basis Gradient Coefficients (y-y coefficient).
      local d_sigma_yz_dz = 0.0 -- Triad Basis Gradient Coefficients (y-z coefficient).
      local d_sigma_zx_dz = 0.0 -- Triad Basis Gradient Coefficients (z-x coefficient).
      local d_sigma_zy_dz = 0.0 -- Triad Basis Gradient Coefficients (z-y coefficient).
      local d_sigma_zz_dz = 0.0 -- Triad Basis Gradient Coefficients (z-z coefficient).

      return d_sigma_xx_dx, d_sigma_xy_dx, d_sigma_xz_dx, d_sigma_yx_dx, d_sigma_yy_dx, d_sigma_yz_dx, d_sigma_zx_dx, d_sigma_zy_dx, d_sigma_zz_dx, 
             d_sigma_xx_dy, d_sigma_xy_dy, d_sigma_xz_dy, d_sigma_yx_dy, d_sigma_yy_dy, d_sigma_yz_dy, d_sigma_zx_dy, d_sigma_zy_dy, d_sigma_zz_dy, 
             d_sigma_xx_dz, d_sigma_xy_dz, d_sigma_xz_dz, d_sigma_yx_dz, d_sigma_yy_dz, d_sigma_yz_dz, d_sigma_zx_dz, d_sigma_zy_dz, d_sigma_zz_dz 
    end,

    -- Velocity space grid.
    lower = { -vx_max, -vy_max, -vz_max },
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
          if math.abs(x) < midplane then
            n = nl -- Total number density (left/inner).
          else
            n = nx -- Total number density (right/outer).
          end

          local metric_det = 1.0

          return metric_det * n
        end,
        temperatureInit = function (t, xn)
          local x = xn[1]

          local T = 0.0
          if math.abs(x) < midplane then
            T = Tl -- Isotropic temperature (left/inner).
          else
            T = Tr -- Isotropic temperature (right/outer).
          end

          return T
        end,
        driftVelocityInit = function (t, xn)
          local x = xn[1]

          local V_x_drift = 0.0
          local V_y_drift = 0.0
          local V_z_drift = 0.0

          if math.abs(x) < midplane then
            V_x_drift = V_x_drift_l -- x drift velocity (left/inner).
            V_y_drift = V_y_drift_l -- y drift velocity (left/inner).
            V_z_drift = V_z_drift_l -- y drift velocity (left/inner).
          else
            V_x_drift = V_x_drift_r -- x drift velocity (right/outer).
            V_y_drift = V_y_drift_r --y drift velocity (right/outer).
            V_z_drift = V_z_drift_r --y drift velocity (right/outer).
          end

          return V_x_drift, V_y_drift, V_z_drift
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