local Vlasov = G0.Vlasov
-- Wald solution referenced from Section 3.2 of the Tetrad-First paper: https://arxiv.org/pdf/2410.02549
-- Axisymmetry in phi, startign from a uniform backgrond B-field
-- Also see: K. Parfery et al. Introducing PHAEDRA: a new spectral code for simulations of relativistic magnetospheres

-- Mathematical constants (dimensionless).
pi = math.pi

-- Simulation parameters.
Nr = 48 -- Cell count (r-direction).
Ntheta = 48 -- Cell count (theta-direction).
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

massBH = 0.3 -- Mass of the black hole
spinBH = 0.0 -- Spin parameter, a = J/M, of the black hole ( Kerr-Schild, coordinates , |a| <= M).

vlasovApp = Vlasov.App.new {
  -- Model ID sets Maxwell's Eq type
  modelID = G0.Model.TriadGR, -- General Relativity Maxwell's used here

  -- Vlasov Geometry
  -- NOTE: Kerr-Schild spherical assumed for DG-GR-Kinetic and DG-GR-Maxwell
  geom = Vlasov.Geom.new {
    
    -- Black hole parameters
    massBH = massBH,
    spinBH = spinBH,

    -- Specify the geometry rtheta choice
    usePresetGeom = true,
    triadPresetGeomType = G0.TriadGeom.GR_KS_rtheta
    
  },

  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedL2fCalcs = integrated_L2_f_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { 1.8 * massBH, math.pi / 8.0 },
  upper = { 5.0 * massBH, 7.0 * math.pi / 8.0 },
  cells = { Nr, Ntheta },
  cflFrac = cfl_frac,
    
  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (r-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = {  }, -- Periodic directions (x- and theta-directions only).

  -- Field.
  field = Vlasov.Field.new {

    -- Use GR field ID
    fieldID = G0.FieldModel.GR,

    -- Initial conditions function.
    init = function (t, xn)
      local r, theta = xn[1], xn[2]

      -- D^i, and B^i (contravaraint components)
      -- (TODO: Compute the spherical KS initial conditions)
      local Dr = 0.0 -- Total electric field (r-direction).
      local Dtheta = 0.0 -- Total electric field (theta-direction).
      local Dphi = 0.0 -- Total electric field (phi-direction).
    
      local Br = 0.0 -- Total magnetic field (r-direction).
      local Btheta = 0.0 -- Total magnetic field (theta-direction).
      local Bphi = 0.0 -- Total magnetic field (phi-direction).

      -- Must return conserved variables
      local rho = math.sqrt(r * r + spinBH * spinBH * math.cos(theta) * math.cos(theta) )
      local metric_det = rho * math.sqrt(2*massBH*r + rho * rho) * math.sin(theta) 
      
      -- Compute Jc * D^i and Jc * B^i
      local JDr = metric_det * Dr
      local JDtheta = metric_det * Dtheta
      local JDphi = metric_det * Dphi
      local JBr = metric_det * Br
      local JBtheta = metric_det * Btheta
      local JBphi = metric_det * Bphi

      -- Hand off the conserved varaibles (J * Q^\xi)
      return JDr, JDtheta, JDphi, JBr, JBtheta, JBphi, 0.0, 0.0
    end,

    -- (TODO: Will need outflow / excision BCs at r = 0, but Copy may be fine since no charateristics escape.)
    -- (TODO: Theta boundaries need an added (reflecting?) BC that doesn't currently exist )
    bcx = { G0.FieldBc.bcCopy, G0.FieldBc.bcCopy }, -- boundary conditions (r-direction).
    bcy = { G0.FieldBc.bcCopy, G0.FieldBc.bcCopy } -- boundary conditions (theta-direction).
  }
}

-- Run application.
vlasovApp:run()