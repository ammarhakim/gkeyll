local Vlasov = G0.Vlasov

-- Mathematical constants (dimensionless).
pi = math.pi

local function noise_init(noise_amp, noise_index, k_init, k_final, Lx, Ly, x, y)
  math.randomseed(0)
  local kindex = (noise_index + 1.0) / 2.0
  local B_amp = 0.0
  local B_phase = 0.0
  local Bx_noise = 0.0
  local By_noise = 0.0
  local Jz_noise = 0.0
  for i = k_init, k_final do 
    B_amp = math.random()
    B_phase = math.random()

    Bx_noise = Bx_noise - 2.0*(2.0*pi/Ly)*(Lx/(i*2.0*pi))*B_amp*math.sin(2.0*pi*y/Ly)*(math.cos(2.0*pi*y/Ly)+1)*math.cos(i*2.0*pi*x/Lx +  2.0*pi*B_phase)*i^kindex
    By_noise = By_noise + B_amp*(math.cos(2.0*pi*y/Ly) + 1.0)*(math.cos(2.0*pi*y/Ly) + 1.0)*math.sin(i*2.0*pi*x/Lx + 2.0*pi*B_phase)*i^kindex
    Jz_noise = Jz_noise + (2.0*pi*i/Lx)*B_amp*(math.cos(2.0*pi*y/Ly) + 1.0)*(math.cos(2.0*pi*y/Ly) + 1.0)*math.cos(i*2.0*pi*x/Lx + 2*pi*B_phase)*i^kindex + 
                 2.0*(2.0*pi/Ly)*(2.0*pi/Ly)*(Lx/(i*2.0*pi))*B_amp*(math.sin(2.0*pi*y/Ly)*math.sin(2.0*pi*y/Ly) - math.cos(2.0*pi*y/Ly)*(math.cos(2.0*pi*y/Ly)+1.0))*math.cos(i*2.0*pi*x/Lx +  2.0*pi*B_phase)*i^kindex
  end

  local kdiff = k_final - k_init + 1.0;
  Bx_noise = noise_amp*Bx_noise/math.sqrt(2.0*kdiff*kdiff/3.0)
  By_noise = noise_amp*By_noise/math.sqrt(2.0*kdiff*kdiff/3.0)
  Jz_noise = noise_amp*Jz_noise/math.sqrt(2.0*kdiff*kdiff/3.0)

  return Bx_noise, By_noise, Jz_noise
end

-- Physical constants (using normalized code units).
epsilon0 = 1.0 -- Permittivity of free space.
mu0 = 1.0 -- Permeability of free space.
light_speed = 1.0/math.sqrt(epsilon0*mu0) -- Speed of light. 
mass_elc = 1.0 -- Electron mass.
charge_elc = -1.0 -- Electron charge.
mass_ion = 1.0 -- Positron mass.
charge_ion = 1.0 -- Positron charge.

sigma = 1.0 -- B^2/(mu0*n0*m*c^2)
n0 = 1.0 -- reference density
T0 = 1.0 -- Reference temperature in units of mc^2

-- Derived parameters
wpe = math.sqrt(charge_ion^2*n0/(epsilon0*mass_elc))
B0 = math.sqrt(sigma*mu0*n0*mass_elc*light_speed^2) -- in-plane magnetic field strength
vA = math.sqrt(sigma/(sigma + 1.0))*light_speed -- in-plane Alfven velocity
omegaCi = charge_ion*B0/mass_ion
de = light_speed/wpe

-- Reconnection parameters
guide = 0.01 -- Guide-field strength. 
w0 = de -- Layer width.
psi0 = 0.1*B0*de -- Layer perturbation.
noise_amp = 0.001*B0 -- Noise amplitude.
k_init = 1 -- First wave mode to perturb with noise, 1.0 correspond to box size.
k_final = 20 -- Last wave mode to perturb with noise.
noise_index = -1.0 -- Spectral index of the noise.

-- Simulation parameters.
Nx = 8 -- Cell count (configuration space: x-direction).
Ny = 4 -- Cell count (configuration space: y-direction).
Nvx = 8 -- Cell count (velocity space: vx-direction).
Nvy = 8 -- Cell count (velocity space: vy-direction).
Nvz = 8 -- Cell count (velocity space: vz-direction).
Lx = 8.0*pi*de  -- Domain size (configuration space: x-direction).
Ly = 4.0*pi*de  -- Domain size (configuration space: x-direction).
vx_max = 8.0 -- Domain boundary (velocity space: vx-direction).
vy_max = 8.0 -- Domain boundary (velocity space: vy-direction).
vz_max = 8.0 -- Domain boundary (velocity space: vz-direction).
nonuniform_v_pow = 2.0 -- Quadratic velocity map. 
vx_linear_res = 1.0/32.0 -- Transition from linear to quadratic velocity map for vx. 
vy_linear_res = 1.0/32.0 -- Transition from linear to quadratic velocity map for vy. 
vz_linear_res = 1.0/32.0 -- Transition from linear to quadratic velocity map for vz. 
poly_order = 2 -- Polynomial order.
basis_type = "tensor" -- Basis function set.
time_stepper = "rk3" -- Time integrator.
cfl_frac = 1.0 -- CFL coefficient.

t_end = 0.1 -- Final simulation time.
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
  lower = { -0.5 * Lx, -0.5 * Ly },
  upper = { 0.5 * Lx, 0.5 * Ly },
  cells = { Nx, Ny },
  cflFrac = cfl_frac,

  basis = basis_type,
  polyOrder = poly_order,
  timeStepper = time_stepper,

  -- Decomposition for configuration space.
  decompCuts = { 1, 1 }, -- Cuts in each coodinate direction (x-direction only).

  -- Boundary conditions for configuration space.
  periodicDirs = { 1 }, -- Periodic directions (x-direction only).

  -- Electrons.
  elc = Vlasov.Species.new {
    modelID = G0.Model.SR,
    charge = charge_elc, mass = mass_elc,
    
    -- Velocity space grid.
    lower = { -1.0, -1.0, -1.0 },
    upper = { 1.0, 1.0, 1.0 },
    cells = { Nvx, Nvy, Nvz },

    mapc2pVel = { 
      -- vx mapping 
      { 
        vmap = function (t, xn)
          local vc = xn[1]
          local vp = 0.0
          local ncells_linear = Nvx/2

          if (vc < 0.0) then 
            vp = vx_linear_res*ncells_linear*vc - vx_max*vc^nonuniform_v_pow
          else
            vp = vx_linear_res*ncells_linear*vc + vx_max*vc^nonuniform_v_pow
          end
          return vp
        end
      },
      -- vy mapping 
      { 
        vmap = function (t, xn)
          local vc = xn[1]
          local vp = 0.0
          local ncells_linear = Nvy/2

          if (vc < 0.0) then 
            vp = vy_linear_res*ncells_linear*vc - vy_max*vc^nonuniform_v_pow
          else
            vp = vy_linear_res*ncells_linear*vc + vy_max*vc^nonuniform_v_pow
          end
          return vp
        end
      },  
      -- vz mapping 
      { 
        vmap = function (t, xn)
          local vc = xn[1]
          local vp = 0.0
          local ncells_linear = Nvz/2

          if (vc < 0.0) then 
            vp = vz_linear_res*ncells_linear*vc - vz_max*vc^nonuniform_v_pow
          else
            vp = vz_linear_res*ncells_linear*vc + vz_max*vc^nonuniform_v_pow
          end
          return vp
        end
      },      
    },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          return n0
        end,
        temperatureInit = function (t, xn)
          return T0
        end,
        driftVelocityInit = function (t, xn)
          local x, y = xn[1], xn[2]
          local sech_sq = (1.0 / math.cosh(y / w0))^2

          local Bx_noise, By_noise, Jz_noise = noise_init(noise_amp, noise_index, k_init, k_final, Lx, Ly, x, y)

          -- Half the current goes to electrons and half to positrons
          local Jx = -B0/w0*math.tanh(y/w0)*sech_sq/(math.sqrt(sech_sq + guide*guide));
          local Jy = 0.0;
          local Jz = B0/w0*sech_sq - psi0*math.cos(2.0*pi*x/Lx)*math.cos(pi*y/Ly)*((2.0*pi/Lx)^2 + (pi/Ly)^2) + Jz_noise;

          local vdrift_x = 0.5*Jx/charge_elc
          local vdrift_y = 0.5*Jy/charge_elc
          local vdrift_z = 0.5*Jz/charge_elc
          return vdrift_x, vdrift_y, vdrift_z 
        end,

        correctAllMoments = true,
        useLastConverged = true
      },
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
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.LTEMoments }
  },

  -- Positrons.
  ion = Vlasov.Species.new {
    modelID = G0.Model.SR,
    charge = charge_ion, mass = mass_ion,
    
    -- Velocity space grid.
    lower = { -1.0, -1.0, -1.0 },
    upper = { 1.0, 1.0, 1.0 },
    cells = { Nvx, Nvy, Nvz },

    mapc2pVel = { 
      -- vx mapping 
      { 
        vmap = function (t, xn)
          local vc = xn[1]
          local vp = 0.0
          local ncells_linear = Nvx/2

          if (vc < 0.0) then 
            vp = vx_linear_res*ncells_linear*vc - vx_max*vc^nonuniform_v_pow
          else
            vp = vx_linear_res*ncells_linear*vc + vx_max*vc^nonuniform_v_pow
          end
          return vp
        end
      },
      -- vy mapping 
      { 
        vmap = function (t, xn)
          local vc = xn[1]
          local vp = 0.0
          local ncells_linear = Nvy/2

          if (vc < 0.0) then 
            vp = vy_linear_res*ncells_linear*vc - vy_max*vc^nonuniform_v_pow
          else
            vp = vy_linear_res*ncells_linear*vc + vy_max*vc^nonuniform_v_pow
          end
          return vp
        end
      },  
      -- vz mapping 
      { 
        vmap = function (t, xn)
          local vc = xn[1]
          local vp = 0.0
          local ncells_linear = Nvz/2

          if (vc < 0.0) then 
            vp = vz_linear_res*ncells_linear*vc - vz_max*vc^nonuniform_v_pow
          else
            vp = vz_linear_res*ncells_linear*vc + vz_max*vc^nonuniform_v_pow
          end
          return vp
        end
      },      
    },

    -- Initial conditions.
    numInit = 1,
    projections = {
      {
        projectionID = G0.Projection.LTE,

        densityInit = function (t, xn)
          return n0
        end,
        temperatureInit = function (t, xn)
          return T0
        end,
        driftVelocityInit = function (t, xn)
          local x, y = xn[1], xn[2]
          local sech_sq = (1.0 / math.cosh(y / w0))^2

          local Bx_noise, By_noise, Jz_noise = noise_init(noise_amp, noise_index, k_init, k_final, Lx, Ly, x, y)

          -- Half the current goes to electrons and half to positrons
          local Jx = -B0/w0*math.tanh(y/w0)*sech_sq/(math.sqrt(sech_sq + guide*guide));
          local Jy = 0.0;
          local Jz = B0/w0*sech_sq - psi0*math.cos(2.0*pi*x/Lx)*math.cos(pi*y/Ly)*((2.0*pi/Lx)^2 + (pi/Ly)^2) + Jz_noise;

          local vdrift_x = 0.5*Jx/charge_ion
          local vdrift_y = 0.5*Jy/charge_ion
          local vdrift_z = 0.5*Jz/charge_ion
          return vdrift_x, vdrift_y, vdrift_z 
        end,

        correctAllMoments = true,
        useLastConverged = true
      },
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
    diagnostics = { G0.Moment.M0, G0.Moment.M1, G0.Moment.LTEMoments }
  },

  field = Vlasov.Field.new {
    epsilon0 = epsilon0, mu0 = mu0,

    -- Initial conditions function.
    init = function (t, xn)
      local x, y = xn[1], xn[2]
      local sech_sq = (1.0 / math.cosh(y / w0))^2
      local Bx_noise, By_noise, Jz_noise = noise_init(noise_amp, noise_index, k_init, k_final, Lx, Ly, x, y)

      local b1x = -B0*math.tanh(y/w0)
      local b1y = 0.0;
      local b1z = B0*math.sqrt(guide*guide + sech_sq);

      local Ex = 0.0;
      local Ey = 0.0;
      local Ez = 0.0;
      local Bx = b1x + psi0 * (pi / Ly) * math.cos(2 * pi * x / Lx) * math.sin(pi * y / Ly) + Bx_noise;
      local By = b1y - psi0 * (2 * pi / Lx) * math.sin(2 * pi * x / Lx) * math.cos(pi * y / Ly) + By_noise;
      local Bz = b1z;

      return Ex, Ey, Ez, Bx, By, Bz, 0.0, 0.0
    end,
    bcy = {
      lower = {
        type = G0.FieldBc.bcCopy
      },
      upper = {
        type = G0.FieldBc.bcCopy
      }
    },
    
    evolve = true, -- Evolve field?
    elcErrorSpeedFactor = 0.0,
    mgnErrorSpeedFactor = 0.0
  }
}

vlasovApp:run()