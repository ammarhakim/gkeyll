local Moments = G0.Moments
local GRTwoFluid = G0.Moments.Eq.GRTwoFluid
local Minkowski = G0.Moments.Spacetime.Minkowski

-- Physical constants (using normalized code units).
gas_gamma_elc = 2.0 -- Adiabatic index (electrons).
gas_gamma_ion = 2.0 -- Adiabatic index (ions).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
mass_ion = 1.0 -- Proton mass.
charge_ion = 1.0 -- Proton charge.
mass_elc = 1.0 / 1836.2 -- Electron mass.
charge_elc = -1.0 -- Electron charge.

rhol_ion = 1.0 -- Left ion mass density.
rhor_ion = 0.125 -- Right ion mass density
pl = 1.0 -- Left electron/ion pressure.
pr = 0.1 -- Right electron/ion pressure.

Bx = 0.5e-2 -- Total magnetic field (x-direction).
Byl = 1.0e-2 -- Left total magneic field (y-direction).
Byr = -1.0e-2 -- Right total magnetic field (y-direction).

has_collision = false -- Whether to include collisions.
nu_base_ei = 0.5 -- Base electron-ion collision frequency.

light_speed = 1.0 -- Speed of light.
e_fact = 0.0 -- Factor of speed of light for electric field correction.
b_fact = 0.0 -- Factor of speed of light for magnetic field correction.

-- Derived physical quantities (using normalized code units).
rhol_elc = rhol_ion * mass_elc / mass_ion -- Left electron mass density.
rhor_elc = rhor_ion * mass_elc / mass_ion -- Right electron mass density.

-- Simulation parameters.
Nx = 4096 -- Cell count (x-direction).
Lx = 10.0 -- Domain size (x-direction).
cfl_frac = 0.95 -- CFL coefficient.

reinit_freq = 100 -- Spacetime reinitialization frequency.

t_end = 2.0 -- Final simulation time.
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
    equation = GRTwoFluid.new {
      massElc = mass_elc,
      massIon = mass_ion,
      chargeElc = charge_elc,
      chargeIon = charge_ion,
      gasGammaElc = gas_gamma_elc,
      gasGammaIon = gas_gamma_ion,
      reinitFreq = reinit_freq
    },

    hasGRTwoFluid = true,
    GRTwoFluidMassElc = mass_elc,
    GRTwoFluidMassIon = mass_ion,
    GRTwoFluidChargeElc = charge_elc,
    GRTwoFluidChargeIon = charge_ion,
    GRTwoFluidGasGammaElc = gas_gamma_elc,
    GRTwoFluidGasGammaIon = gas_gamma_ion,
  
    -- Initial conditions function.
    init = function (t, xn)
      local x = xn[1]
      
      local Dx = 0.0 -- Total electric field (x-direction).
      local Dy = 0.0 -- Total electric field (y-direction).
      local Dz = 0.0 -- Total electric field (z-direction).

      local By = 0.0
      local Bz = 0.0 -- Total magnetic field (z-direction).

      local rhoe = 0.0
      local rhoi = 0.0
      local p = 0.0

      if x < 0.5 * Lx then
        rhoe = rhol_elc -- Electron mass density (left).
        rhoi = rhol_ion -- Ion mass density (left).
        p = pl -- Electron/ion pressure (left).
      else
        rhoe = rhor_elc -- Electron mass density (right).
        rhoi = rhor_ion -- Ion mass density (right).
        p = pr -- Electron/ion pressure (right).
      end

      if x < 0.5 * Lx then
        By = Byl -- Total magnetic field (y-direction, left).
      else
        By = Byr -- Total magnetic field (y-direction, right).
      end

      local lapse = Minkowski.lapseFunction(0.0, x, 0.0, 0.0)
      local shift = Minkowski.shiftVector(0.0, x, 0.0, 0.0)
      local spatial_metric = Minkowski.spatialMetricTensor(0.0, x, 0.0, 0.0)
      local spatial_det = Minkowski.spatialMetricDeterminant(0.0, x, 0.0, 0.0)
      local extrinsic_curvature = Minkowski.extrinsicCurvatureTensor(0.0, x, 0.0, 0.0, 1.0, 1.0, 1.0)
      local in_excision_region = Minkowski.excisionRegion(0.0, x, 0.0, 0.0)

      local lapse_der = Minkowski.lapseFunctionDer(0.0, x, 0.0, 0.0, 1.0, 1.0, 1.0)
      local shift_der = Minkowski.shiftVectorDer(0.0, x, 0.0, 0.0, 1.0, 1.0, 1.0)
      local spatial_metric_der = Minkowski.spatialMetricTensorDer(0.0, x, 0.0, 0.0, 1.0, 1.0, 1.0)

      local We = 1.0
      local Wi = 1.0

      local he = 1.0 + ((p / rhoe) * (gas_gamma_elc / (gas_gamma_elc - 1.0)))
      local hi = 1.0 + ((p / rhoi) * (gas_gamma_ion / (gas_gamma_ion - 1.0)))

      local rhoe_rel = math.sqrt(spatial_det) * rhoe * We -- Electron relativistic mass density.
      local mome_x = 0.0 -- Electron momentum density (x-direction).
      local mome_y = 0.0 -- Electron momentum density (y-direction).
      local mome_z = 0.0 -- Electron momentum density (z-direction).
      local Ee_tot = math.sqrt(spatial_det) * ((rhoe * he * (We * We)) - p - (rhoe * We)) -- Electron total energy density.

      local rhoi_rel = math.sqrt(spatial_det) * rhoi * Wi -- Ion relativistic mass density.
      local momi_x = 0.0 -- Ion momentum density (x-direction).
      local momi_y = 0.0 -- Ion momentum density (y-direction).
      local momi_z = 0.0 -- Ion momentum density (z-direction).
      local Ei_tot = math.sqrt(spatial_det) * ((rhoi * hi * (Wi * Wi)) - p - (rhoi * Wi)) -- Ion total energy density.

      local excision = 0.0
      if in_excision_region then
        rhoe_rel, mome_x, mome_y, mome_z, Ee_tot, rhoi_rel, momi_x, momi_y, momi_z, Ei_tot = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        Dx, Dy, Dz, Bx, By, Bz, lapse = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        for i = 1, 3 do
          shift[i] = 0.0
          lapse_der[i] = 0.0
          for j = 1, 3 do
            spatial_metric[i][j] = 0.0
            extrinsic_curvature[i][j] = 0.0
            shift_der[i][j] = 0.0
            for k = 1, 3 do
              spatial_metric_der[i][j][k] = 0.0
            end
          end  
        end
        
        excision = -1.0
      else
        excision = 1.0
      end
    
      return rhoe_rel, mome_x, mome_y, mome_z, Ee_tot,
        rhoi_rel, momi_x, momi_y, momi_z, Ei_tot,
        Dx, Dy, Dz, Bx, By, Bz, 0.0, 0.0,
        lapse,
        shift[1], shift[2], shift[3],
        spatial_metric[1][1], spatial_metric[1][2], spatial_metric[1][3],
        spatial_metric[2][1], spatial_metric[2][2], spatial_metric[2][3],
        spatial_metric[3][1], spatial_metric[3][2], spatial_metric[3][3],
        extrinsic_curvature[1][1], extrinsic_curvature[1][2], extrinsic_curvature[1][3],
        extrinsic_curvature[2][1], extrinsic_curvature[2][2], extrinsic_curvature[2][3],
        extrinsic_curvature[3][1], extrinsic_curvature[3][2], extrinsic_curvature[3][3],
        excision,
        lapse_der[1], lapse_der[2], lapse_der[3],
        shift_der[1][1], shift_der[1][2], shift_der[1][3],
        shift_der[2][1], shift_der[2][2], shift_der[2][3],
        shift_der[3][1], shift_der[3][2], shift_der[3][3],
        spatial_metric_der[1][1][1], spatial_metric_der[1][1][2], spatial_metric_der[1][1][3],
        spatial_metric_der[1][2][1], spatial_metric_der[1][2][2], spatial_metric_der[1][2][3],
        spatial_metric_der[1][3][1], spatial_metric_der[1][3][2], spatial_metric_der[1][3][3],
        spatial_metric_der[2][1][1], spatial_metric_der[2][1][2], spatial_metric_der[2][1][3],
        spatial_metric_der[2][2][1], spatial_metric_der[2][2][2], spatial_metric_der[2][2][3],
        spatial_metric_der[2][3][1], spatial_metric_der[2][3][2], spatial_metric_der[2][3][3],
        spatial_metric_der[3][1][1], spatial_metric_der[3][1][2], spatial_metric_der[3][1][3],
        spatial_metric_der[3][2][1], spatial_metric_der[3][2][2], spatial_metric_der[3][2][3],
        spatial_metric_der[3][3][1], spatial_metric_der[3][3][2], spatial_metric_der[3][3][3],
        0.0,
        x, 0.0, 0.0
    end,

    evolve = true, -- Evolve species?
    limiter = G0.WaveLimiter.MinMod,
    forceLowOrderFlux = false, -- Use HLL fluxes.
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy }, -- Copy boundary conditions (x-direction).
  }
}

-- Run application.
momentApp:run()
