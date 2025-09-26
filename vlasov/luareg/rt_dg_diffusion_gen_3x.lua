-- General diffusion (with constant diffusion tensor) of a 3D square wave using a p2 DG discretization of the advection-diffusion equation.

local Vlasov = G0.Vlasov
local LinearAdvection = G0.Vlasov.Eq.LinearAdvection

-- Mathematical constants (dimensionless).
pi = math.pi

-- Physical constants (using normalized code units).
v_advect = 1.0 -- Advection velocity.
diffusion_coeff = 1.0 -- Diffusion coefficient.

-- Simulation parameters.
Nx = 16 -- Cell count (configuration space: x-direction).
Ny = 16 -- Cell count (configuration space: y-direction).
Nz = 16 -- Cell count (configuration space: z-direction).
Lx = 4.0 -- Domain size (configuration space: x-direction).
Ly = 4.0 -- Domain size (configuration space: y-direction).
Lz = 4.0 -- Domain size (configuration space: z-direction).
poly_order = 2 -- Polynomial order.
basis_type = "serendipity" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 0.9 -- CFL coefficient.

t_end = 0.01 -- Final simulation time.
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
  lower = { -0.5 * Lx, -0.5 * Ly, -0.5 * Lz },
  upper = { 0.5 * Lx, 0.5 * Ly, 0.5 * Lz },
  cells = { Nx, Ny, Nz },
  cflFrac = cfl_frac,
    
  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1, 1 }, -- Cuts in each coordinate direction (x-, y- and z-directions only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1, 2, 3 }, -- Periodic directions (x-, y- and z-directions only).
  
  -- Fluid.
  fluid = Vlasov.FluidSpecies.new {
    equation = LinearAdvection.new { },

    -- Constant advection function.
    appAdvect = function (t, xn)
      local ux = 0.0 -- Advection velocity (x-direction).
      local uy = 0.0 -- Advection velocity (y-direction).
      local uz = 0.0 -- Advection velocity (z-direction).

      return ux, uy, uz
    end,
    
    -- Initial conditions function.
    init = function (t, xn)
      local x, y, z = xn[1], xn[2], xn[3]

      local f = 0.0
      if math.abs(x) < 1.0 and math.abs(y) < 1.0 and math.abs(z) < 1.0 then
        f = 1.0 -- Advected quantity (interior).
      else
        f = 0.0 -- Advected quantity (exterior).
      end

      return f
    end,

    -- Diffusion.
    diffusion = {
      diffusionTensor = function (t, xn)
        local diffusion_xx = diffusion_coeff -- Diffusion tensor (xx-component).
        local diffusion_xy = diffusion_coeff -- Diffusion tensor (xy-component).
        local diffusion_xz = diffusion_coeff -- Diffusion tensor (xz-component).
        local diffusion_yy = diffusion_coeff -- Diffusion tensor (yy-component).
        local diffusion_yz = diffusion_coeff -- Diffusion tensor (yz-component).
        local diffusion_zz = diffusion_coeff -- Diffusion tensor (zz-component).
  
        return diffusion_xx, diffusion_xy, diffusion_xz, diffusion_yy, diffusion_yz, diffusion_zz
      end
    },

    evolve = true -- Evolve species?
  },

  skipField = true
}

-- Run application.
vlasovApp:run()