local Moments = G0.Moments
local MHDRGFM = G0.Moments.Eq.MHDRGFM

-- Physical constants (using normalized code units).
gas_gamma1 = 1.67 -- First species adiabatic index.
gas_gamma2 = 2.0 -- Second species adiabatic index.

rhol = 1.0 -- Left fluid mass density.
ul = 0.0 -- Left fluid velocity.
pl = 1.0 -- Left fluid pressure.
Bx_l = 0.75 -- Left magnetic field (x-direction).
By_l = 1.0 -- Left magnetic field (y-direction).
phi1_l = 0.99999 -- Left fluid level set value (first species).

rhor = 0.125 -- Right fluid mass density.
ur = 0.0 -- Right fluid velocity.
pr = 0.1 -- Right fluid pressure.
Bx_r = 0.75 -- Right magnetic field (x-direction).
By_r = -1.0 -- Right magnetic field (y-direction).
phi1_r = 0.00001 -- Right fluid level set value (first species).

light_speed = 1.0 -- Speed of light.
b_fact = 0.0 -- Factor of speed of light for magnetic field correction.

-- Simulation parameters.
Nx = 2048 -- Cell count (x-direction).
Lx = 1.0 -- Domain size (x-direction).
cfl_frac = 0.95 -- CFL coefficient.
reinit_freq = 3 -- Reinitialization frequency (for level set).
surface_tension = 0.0 -- Surface tension (for level set).

t_end = 0.1 -- Final simulation time.
num_frames = 1 -- Number of output frames.
field_energy_calcs = GKYL_MAX_INT -- Number of times to calculate field energy.
integrated_mom_calcs = GKYL_MAX_INT -- Number of times to calculate integrated moments.
dt_failure_tol = 1.0e-4 -- Minimum allowable fraction of initial time-step.
num_failures_max = 20 -- Maximum allowable number of consecutive small time-steps.

momentApp = Moments.App.new {
  
  tEnd = t_end,
  nFrame = num_frames,
  fieldEnergyCalcs = field_energy_calcs,
  integratedMomentCalcs = integrated_mom_calcs,
  dtFailureTol = dt_failure_tol,
  numFailuresMax = num_failures_max,
  lower = { 0.0 },
  upper = { Lx },
  cells = { Nx },
  cflFrac = cfl_frac,

  -- Decomposition for configuration space.
  decompCuts = { 1 }, -- Cuts in each coodinate direction (x-direction only).
    
  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).
  
  -- Fluid.
  fluid = Moments.Species.new {
    equation = MHDRGFM.new {
      numComponents = 2,
      gasGamma = {
        gas_gamma1,
        gas_gamma2
      },
      lightSpeed = light_speed,
      mgnErrorSpeedFactor = b_fact,
      reinitFreq = reinit_freq,
      surfaceTension = surface_tension
    },
    
    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]

      local rho1 = 0.0
      local rho2 = 0.0
      local phi1 = 0.0
    
      local vx_total = 0.0
      local vy_total = 0.0
      local vz_total = 0.0
      local p_total = 0.0
      local Bx_total = 0.0
      local By_total = 0.0
      local Bz_total = 0.0
    
      if x < 0.5 then
        rho1 = rhol -- First species fluid mass density (left).
        rho2 = rhor -- Second species fluid mass density (right).
        phi1 = phi1_l -- First species level set value (left).
    
        vx_total = ul -- Total fluid velocity (left).
        p_total = pl -- Total fluid pressure (left).

        Bx_total = Bx_l -- Magnetic field (x-direction, left).
        By_total = By_l -- Magnetic field (y-direction, left).
      else
        rho1 = rhol -- First species fluid mass density (left).
        rho2 = rhor -- Second species fluid mass density (right).
        phi1 = phi1_r -- First species level set value (right).
    
        vx_total = ur -- Total fluid velocity (right).
        p_total = pr -- Total fluid pressure (right).

        Bx_total = Bx_r -- Magnetic field (x-direction, right).
        By_total = By_r -- Magnetic field (y-direction, right).
      end

      local rho_total = (phi1 * rho1) + ((1.0 - phi1) * rho2) -- Total fluid density.

      local momx_total = rho_total * vx_total -- Total fluid momentum density (x-direction).
      local momy_total = rho_total * vy_total -- Total fluid momentum density (y-direction).
      local momz_total = rho_total * vz_total -- Total fluid momentum density (z-direction).
    
      local E1 = (p_total / (gas_gamma1 - 1.0)) + (0.5 * rho1 * (vx_total * vx_total)) + (0.5 * ((Bx_total * Bx_total) + (By_total * By_total) + (Bz_total * Bz_total))) -- First species total energy.
      local E2 = (p_total / (gas_gamma2 - 1.0)) + (0.5 * rho2 * (vx_total * vx_total)) + (0.5 * ((Bx_total * Bx_total) + (By_total * By_total) + (Bz_total * Bz_total))) -- Second species total energy.
      local E_total = (phi1 * E1) + ((1.0 - phi1) * E2) -- Total fluid energy.

      local level_set1 = rho_total * phi1 -- Conserved level set value (first species).
      local mass_frac1 = phi1 * rho1 -- Conserved mass density (first species).
      local mass_frac2 = (1.0 - phi1) * rho2 -- Conserved mass density (second species).
      
      return rho_total, momx_total, momy_total, momz_total, E_total, Bx_total, By_total, Bz_total, 0.0, level_set1, mass_frac1, mass_frac2, 0.0
    end,
  
    evolve = true, -- Evolve species?
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy } -- Copy boundary conditions (x-direction).
  }
}

-- Run application.
momentApp:run()
