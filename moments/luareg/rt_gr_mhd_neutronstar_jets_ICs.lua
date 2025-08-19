local Moments = G0.Moments
local GRMHD = G0.Moments.Eq.GRMHD
local NeutronStar = G0.Moments.Spacetime.NeutronStar

-- Physical constants (using normalized code units).
gas_gamma = 4.0 / 3.0 -- Adiabatic index.

light_speed = 1.0 -- Speed of light.
b_fact = 0.8 -- Factor of speed of light for magnetic field correction.

-- Spacetime parameters (using geometric units).
mass = 1.0 -- Mass of the neutron star.
spin = -0.069 -- Spin of the neutron star.

pos_x = 50.0 -- Position of the neutron star (x-direction).
pos_y = 5.0 -- Position of the neutron star (y-direction).
pos_z = 50.0 -- Position of the neutron star (z-direction).

-- Derived spacetime parameters.
ang_mom = spin * mass -- Angular momentum of the neutron star.
spin_dimensionless = ang_mom / (mass * mass) -- Dimensionless spin of the neutron star.

alpha = 5.0 -- Neutron star alpha parameter (mass quadrupole moment coefficient).
beta = math.pow(-0.36 + (1.48 * math.pow(math.sqrt(alpha), 0.65)), 3.0) -- Neutron star beta parameter (spin octupole moment coefficient).
gamma = math.pow(-4.749 + (0.27613 * math.pow(math.sqrt(alpha), 1.5146)) + (5.5168 * math.pow(math.sqrt(alpha), 0.22229)), 4.0) -- Neutron star gamma parameter (mass hexadecapole moment coefficient).

mass_quadrupole = -alpha * (spin_dimensionless * spin_dimensionless) * (mass * mass * mass) -- Neutron star mass quadrupole.
spin_octupole = -beta * (spin_dimensionless * spin_dimensionless * spin_dimensionless) * (mass * mass * mass * mass) -- Neutron star spin octupole.
mass_hexadecapole = gamma * (spin_dimensionless * spin_dimensionless * spin_dimensionless * spin_dimensionless) * (mass * mass * mass * mass * mass) -- Neutron star mass hexadecapole.

-- Simulation parameters.
Nx = 32 -- Cell count (x-direction).
Ny = 32 -- Cell count (y-direction).
Nz = 32 -- Cell count (z-direction).
Lx = 100.0 -- Domain size (x-direction).
Ly = 10.0 -- Domain size (y-direction).
Lz = 100.0 -- Domain size (z-direction).
cfl_frac = 0.9 -- CFL coefficient.

spacetime_gauge = G0.SpacetimeGauge.Static -- Spacetime gauge choice.
reinit_freq = 10 -- Spacetime reinitialization frequency.

t_end = 50.0 -- Final simulation time.
num_frames = 100 -- Number of output frames.
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
  lower = { 0.0, 0.0, 0.0 },
  upper = { Lx, Ly, Lz },
  cells = { Nx, Ny, Nz },
  cflFrac = cfl_frac,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1, 1 }, -- Cuts in each coodinate direction (x- and y-directions).
  
  -- Boundary conditions for configuration space.
  periodicDirs = { }, -- Periodic directions (none).

  -- Fluid.
  fluid = Moments.Species.new {
    equation = GRMHD.new {
      gasGamma = gas_gamma,
      lightSpeed = light_speed,
      mgnErrorSpeedFactor = b_fact,
      neutronStarParameters = {
        mass = mass,
        spin = spin,
        massQuadrupole = mass_quadrupole,
        spinOctupole = spin_octupole,
        massHexadecapole = mass_hexadecapole,
        posX = pos_x,
        posY = pos_y,
        posZ = pos_z
      },
      spacetimeGauge = spacetime_gauge,
      reinitFreq = reinit_freq
    },

    hasGRMHD = true,
    GRMHDGasGamma = gas_gamma,
  
    -- Initial conditions function.
    init = function (t, xn)
      local x, y, z = xn[1], xn[2], xn[3]
      local r = math.sqrt(((x - 50.0) * (x - 50.0)) + ((y - 5.0) * (y - 5.0)) + ((z - 50.0) * (z - 50.0)))
      local theta = math.atan((y - 5.0) / (x - 5.0))
      local phi = math.atan((((x - 50.0) * (x - 50.0)) + ((y - 5.0) * (y - 5.0))) / (z - 50.0))
      --print(theta)
      local rho = 0.0
      local p = 0.0
      local rho0 = 10
      local r_star = 150740.0
      local omega0 = -6.08113e-6
      local rc = 1507.4
      local Rhe = 21103.5
      local ux = 0.0
      local uy = 0.0

      local Bx = 0.0
      local By = 0.0
      local Bz = 0.0

      if r < r_star then
        rho = rho0 * math.pow(1.21052e38 / (math.pow(r, 2) + 3.63155e37) , 1.25) * math.pow((1 - (r / r_star)) , 5)
        --rho = rho0 * math.pow(1e16 / (math.pow(r, 2) + 1e16) , 1.25) * math.pow((1 - (r / r_star)) , 5)
        p = 3.64945e10 * rho / r
        --p = 1e-6 * rho / r

        if r < rc then
          ux = - omega0 * (r * r * math.pow(math.sin(theta), 2)) * math.sin(theta) * r * math.sin(phi)
          uy = omega0 * (r * r * math.pow(math.sin(theta), 2)) * r * math.cos(theta) * math.sin(phi)
          --ux = 0.0
          --uy = 0.0
        elseif r < Rhe then
          ux = - omega0 * (rc * rc * math.pow(math.sin(theta), 2)) * r * math.sin(theta) * math.sin(phi)
          uy = omega0 * (rc * rc * math.pow(math.sin(theta), 2)) * r * math.cos(theta) * math.sin(phi)
        end

        if r < 2.4 then
          Bz = 8.54651e-8
          --Bz = 1e6
        end
      else
        --rho = 1e-20
        rho = 121.052
        ux = 0.0
        uy = 0.0
        -- p = rho / 1e6
        p = rho / 3.64945e22
      end

      local lapse = NeutronStar.lapseFunction(mass, spin, mass_quadrupole, spin_octupole, mass_hexadecapole, pos_x, pos_y, pos_z, 0.0, x, y, z)
      local shift = NeutronStar.shiftVector(mass, spin, mass_quadrupole, spin_octupole, mass_hexadecapole, pos_x, pos_y, pos_z, 0.0, x, y, z)
      local spatial_metric = NeutronStar.spatialMetricTensor(mass, spin, mass_quadrupole, spin_octupole, mass_hexadecapole, pos_x, pos_y, pos_z, 0.0, x, y, z)
      local spatial_det = NeutronStar.spatialMetricDeterminant(mass, spin, mass_quadrupole, spin_octupole, mass_hexadecapole, pos_x, pos_y, pos_z, 0.0, x, y, z)
      local extrinsic_curvature = NeutronStar.extrinsicCurvatureTensor(mass, spin, mass_quadrupole, spin_octupole, mass_hexadecapole, pos_x, pos_y, pos_z, 0.0, x, y, z,
        math.pow(10.0, -8.0), math.pow(10.0, -8.0), math.pow(10.0, -8.0))
      local in_excision_region = NeutronStar.excisionRegion(mass, spin, mass_quadrupole, spin_octupole, mass_hexadecapole, pos_x, pos_y, pos_z, 0.0, x, y, z)

      local lapse_der = NeutronStar.lapseFunctionDer(mass, spin, mass_quadrupole, spin_octupole, mass_hexadecapole, pos_x, pos_y, pos_z, 0.0, x, y, z,
        math.pow(10.0, -8.0), math.pow(10.0, -8.0), math.pow(10.0, -8.0))
      local shift_der = NeutronStar.shiftVectorDer(mass, spin, mass_quadrupole, spin_octupole, mass_hexadecapole, pos_x, pos_y, pos_z, 0.0, x, y, z,
        math.pow(10.0, -8.0), math.pow(10.0, -8.0), math.pow(10.0, -8.0))
      local spatial_metric_der = NeutronStar.spatialMetricTensorDer(mass, spin, mass_quadrupole, spin_octupole, mass_hexadecapole, pos_x, pos_y, pos_z, 0.0, x, y, z,
        math.pow(10.0, -8.0), math.pow(10.0, -8.0), math.pow(10.0, -8.0))

      local vel = { ux, uy, 0 }
      local v_sq = 0.0

      for i = 1, 3 do
        for j = 1, 3 do
          v_sq = v_sq + (spatial_metric[i][j] * vel[i] * vel[j])
        end
      end

      local W = 1.0 / math.sqrt(1.0 - v_sq)
      if v_sq > 1.0 - math.pow(10.0, -8.0) then
        W = 1.0 / math.sqrt(1.0 - math.pow(10.0, -8.0))
      end

      local mag = { Bx, By, Bz }

      local cov_mag = { 0.0, 0.0, 0.0 }
      for i = 1, 3 do
        for j = 1, 3 do
          cov_mag[i] = cov_mag[i] + (spatial_metric[i][j] * mag[j])
        end
      end

      local cov_vel = { 0.0, 0.0, 0.0 }
      for i = 1, 3 do
        for j = 1, 3 do
          cov_vel[i] = cov_vel[i] + (spatial_metric[i][j] * vel[j])
        end
      end

      local b0 = 0.0
      for i = 1, 3 do
        b0 = b0 + (W * mag[i] * (cov_vel[i] / lapse))
      end

      local spacetime_vel = { W / lapse, 0.0, 0.0, 0.0 }
      for i = 1, 3 do
        spacetime_vel[i + 1] = (W * vel[i]) - (shift[i] * (W / lapse))
      end

      local b = { 0.0, 0.0, 0.0 }
      for i = 1, 3 do
        b[i] = (mag[i] + (lapse * b0 * spacetime_vel[i + 1])) / W
      end

      local b_sq = 0.0
      for i = 1, 3 do
        b_sq = b_sq + ((mag[i] * cov_mag[i]) / (W * W))
      end
      b_sq = b_sq + (((lapse * lapse) * (b0 * b0)) / (W * W))

      local cov_b = { 0.0, 0.0, 0.0 }
      for i = 1, 3 do
        for j = 1, 3 do
          cov_b[i] = cov_b[i] + (spatial_metric[i][j] * b[j])
        end
      end

      local h_star = 1.0 + ((p / rho) * (gas_gamma / (gas_gamma - 1.0))) + (b_sq / rho)
      local p_star = p + (0.5 * b_sq)

      local rho_rel = math.sqrt(spatial_det) * rho * W -- Fluid relativistic mass density.
      local mom_x = math.sqrt(spatial_det) * ((rho * h_star * (W * W) * cov_vel[1]) - (lapse * b0 * cov_b[1])) -- Fluid momentum density (x-direction).
      local mom_y = math.sqrt(spatial_det) * ((rho * h_star * (W * W) * cov_vel[2]) - (lapse * b0 * cov_b[2])) -- Fluid momentum density (y-direction).
      local mom_z = math.sqrt(spatial_det) * ((rho * h_star * (W * W) * cov_vel[3]) - (lapse * b0 * cov_b[3])) -- Fluid momentum density (z-direction).
      local Etot = math.sqrt(spatial_det) * ((rho * h_star * (W * W)) - p_star - ((lapse * lapse) * (b0 * b0)) - (rho * W)) -- Fluid total energy density.

      local Bx_rel = math.sqrt(spatial_det) * Bx -- Fluid relativistic magnetic field (x-direction).
      local By_rel = math.sqrt(spatial_det) * By -- Fluid relativistic magnetic field (y-direction).
      local Bz_rel = math.sqrt(spatial_det) * Bz -- Fluid relativistic magnetic field (z-direction).

      local excision = 0.0
      if in_excision_region then
        rho_rel, mom_x, mom_y, mom_z, Etot, Bx_rel, By_rel, Bz_rel, lapse = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
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
    
      return rho_rel, mom_x, mom_y, mom_z, Etot, Bx_rel, By_rel, Bz_rel, 0.0,
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
        x, y, z
    end,

    evolve = true, -- Evolve species?
    forceLowOrderFlux = true, -- Use Lax fluxes.
    bcx = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy }, -- Copy boundary conditions (x-direction).
    bcy = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy }, -- Copy boundary conditions (y-direction).
    bcz = { G0.SpeciesBc.bcCopy, G0.SpeciesBc.bcCopy } -- Copy boundary conditions (z-direction).
  }
}

-- Run application.
momentApp:run()
