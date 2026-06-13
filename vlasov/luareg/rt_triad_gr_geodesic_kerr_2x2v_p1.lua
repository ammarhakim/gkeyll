local Vlasov = G0.Vlasov

-- Physical constants (using normalized code units).
mass = 1.0 -- Neutral mass.
charge = 0.0 -- Neutral charge.

nl = 1.0 -- Left/inner number density.
Tl = 1.0 -- Left/inner temperature.
V_r_drift_l = 0.0 -- Left/inner drift velocity (radial direction).
V_theta_drift_l = 0.0 -- Left/inner drift velocity (angular direction).

vt = 1.0 -- Thermal velocity.
nu = 15000.0 -- Collision frequency.

-- Simulation parameters.
Nr = 12 -- 96 -- Cell count (configuration space: radial direction).
Ntheta = 12 -- 64 -- Cell count (configuration space: radial direction).
Nvr = 16 -- 200 -- Cell count (velocity space: radial direction).
Nvtheta = 16 -- 86 -- Cell count (velocity space: angular direction).
Lr = 7.5 -- Domain size (configuration space: radial direction).
Ltheta = 2 * math.pi -- Domain size (configuration space: radial direction).
vr_max = 1.0 * vt -- Domain boundary (velocity space: radial direction).
vtheta_max = 0.6 * vt -- Domain boundary (velocity space: angular direction).
poly_order = 1 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 20.0 -- Final simulation time.
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
  lower = { 1.5 , 0.0 },
  upper = { 1.5 + Lr, 3.0 * math.pi/4.0 },
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
    spinBH = spinBH,

    -- Use the vierbein inputs to contruct the Poisson Tensor
    usePresetGeom = true,
    triadPresetGeomType = G0.TriadGeom.GR_KS_rphi
    
  },

  -- Neutral species.
  neut = Vlasov.Species.new {
    modelID = G0.Model.TriadGR,
    charge = charge, mass = mass,

    -- Velocity space grid.
    lower = { -vr_max, -vtheta_max },
    upper = { vr_max, vtheta_max/4.0 },
    cells = { Nvr, Nvtheta },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.Func,

        init = function (t, xn)
          local r = xn[1]
          local phi = xn[2]
          local pr = xn[3]
          local pphi = xn[4]

          -- Widths in each direction
          local sigma_r = 0.005 
          local sigma_phi = 0.0025 
          local sigma_pr = 0.0002
          local sigma_pphi = 0.0002

          -- Peak location (at the geodesic we aim to follow)
          local r0 = 8.0
          local phi0 = math.pi/2
          local pr0 = -0.342762596284762
          local pphi0 = -0.163322023271852

          -- Initial Gaussian bump in phase-space with a small thermal spread
          local f = math.exp( - (r - r0) * (r - r0) / (2 * sigma_r) ) *
                    math.exp( - (phi - phi0) * (phi - phi0) / (2 * sigma_phi) ) *
                    math.exp( - (pr - pr0) * (pr - pr0) / (2 * sigma_pr) ) *
                    math.exp( - (pphi - pphi0) * (pphi - pphi0) / (2 * sigma_pphi) ) 
	  
	  f = f + 1e-10 * math.exp( - (pr) * (pr) / (2 * sigma_pr) ) *
                    math.exp( - (pphi) * (pphi) / (2 * sigma_pphi) )		
          return f
        end
      }
    },

    bcx = {
      lower = {
        type = G0.SpeciesBc.bcAbsorb
      },
      upper = {
        type = G0.SpeciesBc.bcAbsorb
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
