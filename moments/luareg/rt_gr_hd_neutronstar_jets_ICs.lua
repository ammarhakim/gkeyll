-- ==========================================
-- Load MESA profile .data file
-- ==========================================

-- Change this to your MESA profile filename
local mesa_file = "/Users/nikolabukowiecka/Desktop/Nikola/Gkyell/MESAfiles/CHE_test.data"

-- Debug print: make sure this path is exactly what you expect
print("Attempting to open MESA file:", mesa_file)

-- Try to open file
local f, err = io.open(mesa_file, "r")
if not f then
    io.stderr:write("ERROR: Could not open MESA profile file!\n")
    io.stderr:write("Path checked: " .. mesa_file .. "\n")
    io.stderr:write("System error message: " .. tostring(err) .. "\n")
    os.exit(1)
end

-- Skip header lines until we find the column names
local header_found = false
local col_names = {}
while true do
   local line = f:read("*l")
   if not line then break end
   if string.find(line, "logR") then
      for name in string.gmatch(line, "[^%s]+") do
         table.insert(col_names, name)
      end
      header_found = true
      break
   end
end

if not header_found then
    io.stderr:write("ERROR: Could not find column names in MESA profile file!\n")
    os.exit(1)
end

-- Build a map from column names to their indices
local col_index = {}
for i, name in ipairs(col_names) do
   col_index[name] = i
end

-- Tables for variables we care about
local radius = {}
local rhoStar = {}
local pStar = {}
local vphiStar = {}

if not f then
   io.stderr:write("INTERNAL ERROR: f was nil before reading lines.\n")
   os.exit(1)
end

-- Read rest of file line by line
for line in f:lines() do
   local values = {}
   for value in string.gmatch(line, "[^%s]+") do
      table.insert(values, tonumber(value))
   end

   -- Convert radius: MESA gives logR in log10(cm)
   local r = 10^(values[col_index["logR"]]) -- cm
   table.insert(radius, r)

   -- Density: logRho is log10(g/cm^3)
   local d = 10^(values[col_index["logRho"]]) -- g/cm^3
   table.insert(rhoStar, d)

   -- Pressure: logP is log10(dyn/cm^2)
   local p = 10^(values[col_index["logP"]]) -- dyn/cm^2
   table.insert(pStar, p)

   -- Azimuthal velocity: assuming MESA has v_rot column (cm/s)
   local v = values[col_index["v_rot"]] or 0.0
   table.insert(vphiStar, v)
end

f:close()



-- ==========================================
-- Your existing interpolation functions below
-- ==========================================
local findNearestIndex = function(x, x0)
   local nx = #x
   local idxLo = 1
   while (idxLo+1 < nx and x[idxLo+1] <= x0) do
      idxLo = idxLo + 1
   end
   local idxUp = math.min(idxLo + 1, nx)
   return idxLo, idxUp
end

local linearInterp = function(fun, x, x0)
   local idxLo, idxUp = findNearestIndex(x, x0)
   local x0l = math.min(x0, x[#x])
   local y = fun[idxLo] + (x0l - x[idxLo]) * (fun[idxUp] - fun[idxLo]) / (x[idxUp] - x[idxLo])
   return y
end

-- Now radius, rho, press, vphi are 1D tables ready for interpolation

-- Example: nearest neighbor/linear interpolation (reuse your function below)
-- local density_at_r = linearInterp(rho, radius, 1.0e10)

local Moments = G0.Moments
local GRMHD = G0.Moments.Eq.GRMHD
local NeutronStar = G0.Moments.Spacetime.NeutronStar

-- Physical constants (using normalized code units).
gas_gamma = 4.0 / 3.0 -- Adiabatic index.

light_speed = 1.0 -- Speed of light.
b_fact = 0.85 -- Factor of speed of light for magnetic field correction.

-- Spacetime parameters (using geometric units).
mass = 1.0 -- Mass of the neutron star.
spin = -0.11 -- Spin of the neutron star.

pos_x = 20.0 -- Position of the neutron star (x-direction).
pos_y = 2.0 -- Position of the neutron star (y-direction).
pos_z = 20.0 -- Position of the neutron star (z-direction).

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
Nx = 100 -- Cell count (x-direction).
Ny = 32 -- Cell count (y-direction).
Nz = 100 -- Cell count (z-direction).
Lx = 40.0 -- Domain size (x-direction).
Ly = 4.0 -- Domain size (y-direction).
Lz = 40.0 -- Domain size (z-direction).
cfl_frac = 0.9 -- CFL coefficient.

spacetime_gauge = G0.SpacetimeGauge.Static -- Spacetime gauge choice.
reinit_freq = 10 -- Spacetime reinitialization frequency.

t_end = 600.0 -- Final simulation time.
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
      local r = math.sqrt(((x - (0.5 * Lx)) * (x - (0.5 * Lx))) + ((y - (0.5 * Ly)) * (y - (0.5 * Ly))) + ((z - (0.5 * Lz)) * (z - (0.5 * Lz))))
      local theta = math.atan2((y - (0.5 * Ly)), (x - (0.5 * Lx)))
      local phi = math.atan2(math.sqrt(((x - (0.5 * Lx)) * (x - (0.5 * Lx))) + ((y - (0.5 * Ly)) * (y - (0.5 * Ly)))), (z - (0.5 * Lz)))
      
      local function velocity_cartesian(x, y, z)
        local r = math.sqrt(((x - (0.5 * Lx)) * (x - (0.5 * Lx))) + ((y - (0.5 * Ly)) * (y - (0.5 * Ly))) + ((z - (0.5 * Lz)) * (z - (0.5 * Lz))))
        if r == 0 then 
          return 0.0, 0.0, 0.0 
        end

        local phi = math.atan2((y - (0.5 * Ly)), (x - (0.5 * Lx)))
        local vphi_val = linearInterp(vphiStar, radius, r) -- interpolate v_phi(r)
          -- radius: the table of radii you read from the MESA profile file. It is an array of sampled radii where you know v_phi
          -- r: the specific radius at the grid point you’re evaluating the velocity for

          -- Convert to Cartesian
        local vx = -vphi_val * math.sin(phi)
        local vy =  vphi_val * math.cos(phi)
        local vz = 0.0

        return vx, vy, vz
      end
      
      --print(theta)

      local r_star = radius[#radius]
      local ux = 0.0
      local uy = 0.0

      local Bx = 0.0
      local By = 0.0
      local Bz = 0.0

      local vx, vy, vz = velocity_cartesian(x, y, z)

      if r < r_star then
        rho = linearInterp(rhoStar, radius, r)
        p = linearInterp(pStar, radius, r)
        ux = vx
        uy = vy
        uz = vz
    
      else
        rho = 121.052
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
