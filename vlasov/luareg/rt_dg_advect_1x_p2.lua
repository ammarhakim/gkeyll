-- Constant advection of a sine wave using a DG discretization of the advection equation. 

local Vlasov = G0.Vlasov
local Advect = G0.Vlasov.Eq.Advect

-- Simulation parameters.
Nx = 16 -- Cell count (x-direction).
Lx = 2.0*math.pi -- Domain size (x-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 20.0*math.pi -- Final simulation time.
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
  lower = { 0.0 },
  upper = { Lx },
  cells = { Nx },
  cflFrac = cfl_frac,
    
  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1 }, -- Periodic directions.
  
  -- Fluid.
  fluid = Vlasov.FluidSpecies.new {
    equation = Advect.new { },

    -- Constant advection function.
    appAdvect = function (t, xn)
      local x = xn[1]
      local ux = 1.0
      local uy = 0.0
      local uz = 0.0

      return 1.0, 0.0, 0.0
    end,
    
    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]
      
      return math.sin(x)
    end,
  },

  skipField = true
}

-- Run application.
vlasovApp:run()