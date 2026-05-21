local Vlasov = G0.Vlasov

-- 2D kinetic Bondi-Hoyle-Lyttleton-like accretion setup on a Schwarzschild
-- Kerr-Schild background. This uses the same physical box size as
-- rt_gr_bhl_static_tetrad, shifted so the black hole is at the origin.

-- Physical constants (using normalized code units).
mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.

-- Fluid states from moments/luareg/rt_gr_bhl_static_tetrad.lua.
rhol = 3.0 -- Left mass density.
ul = 0.3 -- Left x drift velocity.
pl = 0.05 -- Left pressure.

rhor = 0.01 -- Right mass density.
ur = 0.0 -- Right x drift velocity.
pr = 0.01 -- Right pressure.

-- Spacetime parameters (using geometric units).
massBH = 0.3 -- Mass of the Kerr-Schild black hole.
spinBH = 0.0 -- Spin parameter, a = J/M, of the Kerr-Schild black hole.
excision_radius = 1.9 * massBH -- Kinetic excision region, inside the horizon.
excision_density = 1.0e-4 -- Avoid a singular LTE projection in excised cells.

-- Simulation parameters.
Nx = 12-- 32 -- Cell count (configuration space: x direction).
Ny = 12 -- 32 -- Cell count (configuration space: y direction).
Nvx = 6 -- 12 -- Cell count (velocity space: x direction).
Nvy = 6 -- 12 -- Cell count (velocity space: y direction).
Nvz = 4 -- 12 -- Cell count (velocity space: z direction).
Lx = 5.0 -- Domain size (configuration space: x direction).
Ly = 5.0 -- Domain size (configuration space: y direction).
vx_max = 2.0 -- Domain boundary (velocity space: x direction).
vy_max = 2.0 -- Domain boundary (velocity space: y direction).
vz_max = 2.0 -- Domain boundary (velocity space: z direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 15.0 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
integrated_L2_f_calcs = GKYL_MAX_INT -- Number of times to calculate L2 norm of distribution function.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

x_loc = -1.5 -- Shock location shifted from x=1.0 on [0,5] to [-2.5,2.5].

function in_excision_region(x, y)
  return math.sqrt(x * x + y * y) < excision_radius
end

function sqrt_spatial_det(x, y)
  local r = math.sqrt(x * x + y * y)
  if r < excision_radius then
    r = excision_radius
  end

  return math.sqrt(1.0 + (2.0 * massBH / r))
end

function primitive_state(x, y)
  local rho, ux, p = 0.0, 0.0, 0.0

  if x < x_loc then
    rho = rhol
    ux = ul
    p = pl
  else
    rho = rhor
    ux = ur
    p = pr
  end

  return rho, ux, p
end

vlasovApp = Vlasov.App.new {

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedL2fCalcs = integrated_L2_f_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { -0.5 * Lx, -0.5 * Ly },
  upper = { 0.5 * Lx, 0.5 * Ly },
  cells = { Nx, Ny },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coordinate direction.

  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).

  -- Vlasov Geometry.
  geom = Vlasov.Geom.new {
    massBH = massBH,
    spinBH = spinBH,

    usePresetGeom = true,
    triadPresetGeomType = G0.TriadGeom.GR_KS_Cart_3V
  },

  -- Neutral species.
  neut = Vlasov.Species.new {
    modelID = G0.Model.TriadGR,
    charge = charge, mass = mass,

    -- Velocity space grid.
    lower = { -vx_max, -vy_max, -vz_max },
    upper = { vx_max, vy_max, vz_max },
    cells = { Nvx, Nvy, Nvz },

    -- Skip effectively excised cells.
    skipCellThresh = 1.0e-5,

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          local x, y = xn[1], xn[2]
          if in_excision_region(x, y) then
            return excision_density
          end

          local rho, ux, p = primitive_state(x, y)
          return sqrt_spatial_det(x, y) * rho / mass
        end,
        temperatureInit = function (t, xn)
          local x, y = xn[1], xn[2]
          if in_excision_region(x, y) then
            return pr / rhor
          end

          local rho, ux, p = primitive_state(x, y)
          return p / rho
        end,
        driftVelocityInit = function (t, xn)
          local x, y = xn[1], xn[2]
          if in_excision_region(x, y) then
            return 0.0, 0.0, 0.0
          end

          local rho, ux, p = primitive_state(x, y)
          return ux, 0.0, 0.0
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

    bcy = {
      lower = {
        type = G0.SpeciesBc.bcCopy
      },
      upper = {
        type = G0.SpeciesBc.bcCopy
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
