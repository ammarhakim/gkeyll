-- Gkyl ------------------------------------------------------------------------
--
-- Generate grid for use in simulations of tokamaks
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local argparse = require "Lib.argparse"
local ffi = require "ffi"
local lfs = require "Lib.lfs"
local xsys = require "xsys"
local Time = require "Lib.Time"
local Grid = require "Grid"
local DataStruct = require "DataStruct"

-- Create CLI parser to handle commands and options
local parser = argparse()
   :name("tokagridgen")
   :description [[

Generates a grid for use in tokamak simulations. This tool takes as
input the solution to the Grad-Shafranov equation and constructs nodes
the RZ plane along with a minimum set of geometric quantities needed
in various solvers.

]]

-- add flags and options
parser:option("-e --examples", "Fetch example input file (nstxu_DN, asdex, step)")
parser:option("-i --input", "Input file")
parser:option("-n --name", "Application name (used for app_inp.name)")
parser:option("-t --template", "Print blank template (blankSN or blankDN)")

-- parse command line parameters
local args = parser:parse(GKYL_COMMANDS_L)

if args.examples then
   -- write out example to terminal based on argument
   if args.examples == "nstxu_DN" then
      io.write([[
-- Name of tokamak
name = "nstxu_DN"
-- Input file for tokagridgen Tool. File path must be specified relative to where the tool is called
-- "gyrokinetic/data/eqdsk/nstxu_DN.geqdsk" for double null
--"gyrokinetic/data/eqdsk/asdex.geqdsk" for single null
geqdsk_path = "gyrokinetic/data/eqdsk/nstxu_DN.geqdsk"

--Number of null points. one of:
-- single_null
-- double_null
toka_type = double_null

-- Whether to generate only half the poloidal domain (DN only)
use_half_domain = true

-- Lower Private Flux region parameters
lower_PF = {
    rleft = 0.2,      -- rleft_PF
    rright = 1.6,     -- rright_PF  
    rmin = 0.2,       -- rmin_PF
    rmax = 1.6,       -- rmax_PF
    zmin_left = -1.22, -- zmin_left_PF
    zmin_right = -1.5  -- zmin_right_PF
}

-- Outer Scrape-Off Layer region parameters
outer_SOL = {
    rclose = 1.6,     -- rclose_outer_SOL
    rleft = 0.2,      -- rleft_outer_SOL
    rright = 1.6,     -- rright_outer_SOL
    rmin = 0.2,       -- rmin_outer_SOL
    rmax = 1.6,       -- rmax_outer_SOL
    zmin = -1.5,      -- zmin_outer_SOL
    zmax = 1.5        -- zmax_outer_SOL
}

-- Inner Scrape-Off Layer region parameters
inner_SOL = {
    rleft = 0.2,      -- rleft_inner_SOL
    rright = 1.6,     -- rright_inner_SOL
    rmin = 0.2,       -- rmin_inner_SOL
    rmax = 1.6,       -- rmax_inner_SOL
    zmin = -1.22,     -- zmin_inner_SOL
    zmax = 1.22       -- zmax_inner_SOL
}

-- Core region parameters
core = {
    rclose_inner = 0.2,  -- rclose_inner_core
    rclose_outer = 1.6,  -- rclose_outer_core
    rleft = 0.2,         -- rleft_core
    rright = 1.6,        -- rright_core
    rmin = 0.2,          -- rmin_core
    rmax = 1.6           -- rmax_core
}

-- Relative widths of psi for different regions (fraction of the total psi difference between center and separatrix, from 0 to 1)
-- For DN: (PF Width, Outer SOL Width, Inner SOL Width, Core Width)
-- For SN: (PF Width, SOL Width, Core Width) 
wpsi = {0.1, 0.25, 0.1, 0.9}


-- Number of cells to include along the directions of different blocks
-- (npsi_pf, npsi_sol, npsi_core, ntheta_lower_outer, ntheta_middle_SOL, ntheta_lower_inner)
cell_count = {4, 4, 4, 6, 6, 6}


-- Plate function definitions
-- These functions define the divertor plate geometry as functions of parameter s in [0,1]
function plate_func_lower_outer(s, RZ)
    RZ[0] = 0.3 + 0.7*s  -- R coordinate
    RZ[1] = -1.5         -- Z coordinate
  end
  
  function plate_func_upper_outer(s, RZ)
    RZ[0] = 0.3 + 0.7*s  -- R coordinate
    RZ[1] = 1.5          -- Z coordinate
  end
  
  function plate_func_upper_inner(s, RZ)
    RZ[0] = 0.3          -- R coordinate
    RZ[1] = 1.0 + (1.5 - 1.0)*s  -- Z coordinate
  end
  
  function plate_func_lower_inner(s, RZ)
    RZ[0] = 0.3          -- R coordinate
    RZ[1] = -(1.0 + (1.5 - 1.0)*s)  -- Z coordinate
  end

]]
      )
   elseif args.examples == "asdex" then
      io.write([[
-- Name of tokamak
name = "asdex"
-- Input file for tokagridgen Tool.
-- "gyrokinetic/data/eqdsk/nstxu_DN.geqdsk" for double null
--"gyrokinetic/data/eqdsk/asdex.geqdsk" for single null
geqdsk_path = "gyrokinetic/data/eqdsk/asdex.geqdsk"

--Number of null points. one of:
-- single_null
-- double_null
toka_type = single_null

-- Lower Private Flux region parameters
lower_PF = {
    rleft = 1.1,      -- rleft_PF
    rright = 1.7,     -- rright_PF  
    rmin = 1.1,       -- rmin_PF
    rmax = 1.7,       -- rmax_PF
    zmin = -1.3,      -- zmin_PF
    zmax = -0.9,      -- zmax_PF
    zmin_left = -1.18, -- zmin_left_PF
    zmin_right = -1.18 -- zmin_right_PF
}

-- Scrape-Off Layer region parameters
SOL = {
    rclose = 2.5,     -- rclose_SOL
    rleft = 0.8,      -- rleft_SOL
    rright = 2.5,     -- rright_SOL
    rmin = 0.8,       -- rmin_SOL
    rmax = 2.5,       -- rmax_SOL
    zmin = -1.3,      -- zmin_SOL
    zmax = 1.0        -- zmax_SOL
}

-- Core region parameters
core = {
    rclose = 2.0,     -- rclose_core
    rleft = 0.8,      -- rleft_core
    rright = 2.5,     -- rright_core
    rmin = 0.8,       -- rmin_core
    rmax = 2.5,       -- rmax_core
    zmin = -1.3,      -- zmin_core
    zmax = 1.0        -- zmax_core
}

-- Relative widths of psi for different regions (fraction of the total psi difference between center and separatrix, from 0 to 1)
-- For DN: (PF Width, Outer SOL Width, Inner SOL Width, Core Width)
-- For SN: (PF Width, SOL Width, Core Width) 
wpsi = {0.01, 0.04, 0.1}


-- Number of cells to include along the directions of different blocks
-- (npsi_pf, npsi_sol, npsi_core, ntheta_lower_outer, ntheta_middle_SOL, ntheta_lower_inner)
cell_count = {2, 4, 4, 4, 18, 4}

-- Plate function definitions
-- These functions define the divertor plate geometry as functions of parameter s in [0,1]
function plate_func_lower_outer(s, RZ)
    local npts = 100

    local t = {
    0.        , 0.01010101, 0.02020202, 0.03030303, 0.04040404, 0.05050505,
    0.06060606, 0.07070707, 0.08080808, 0.09090909, 0.1010101 , 0.11111111,
    0.12121212, 0.13131313, 0.14141414, 0.15151515, 0.16161616, 0.17171717,
    0.18181818, 0.19191919, 0.2020202 , 0.21212121, 0.22222222, 0.23232323,
    0.24242424, 0.25252525, 0.26262626, 0.27272727, 0.28282828, 0.29292929,
    0.3030303 , 0.31313131, 0.32323232, 0.33333333, 0.34343434, 0.35353535,
    0.36363636, 0.37373737, 0.38383838, 0.39393939, 0.4040404 , 0.41414141,
    0.42424242, 0.43434343, 0.44444444, 0.45454545, 0.46464646, 0.47474747,
    0.48484848, 0.49494949, 0.50505051, 0.51515152, 0.52525253, 0.53535354,
    0.54545455, 0.55555556, 0.56565657, 0.57575758, 0.58585859, 0.5959596 ,
    0.60606061, 0.61616162, 0.62626263, 0.63636364, 0.64646465, 0.65656566,
    0.66666667, 0.67676768, 0.68686869, 0.6969697 , 0.70707071, 0.71717172,
    0.72727273, 0.73737374, 0.74747475, 0.75757576, 0.76767677, 0.77777778,
    0.78787879, 0.7979798 , 0.80808081, 0.81818182, 0.82828283, 0.83838384,
    0.84848485, 0.85858586, 0.86868687, 0.87878788, 0.88888889, 0.8989899 ,
    0.90909091, 0.91919192, 0.92929293, 0.93939394, 0.94949495, 0.95959596,
    0.96969697, 0.97979798, 0.98989899, 1.        ,
    }
    local R = {
    1.57598284, 1.57749872, 1.57719396, 1.57796683, 1.57991018, 1.58038241,
    1.5809759 , 1.58399201, 1.58425793, 1.58426253, 1.58499896, 1.58614281,
    1.58753427, 1.58902876, 1.5904817 , 1.59179306, 1.59302565, 1.59427949,
    1.59565459, 1.59723624, 1.59893569, 1.60055062, 1.60187675, 1.60274568,
    1.60335474, 1.60411485, 1.60542574, 1.60728998, 1.60926429, 1.61088179,
    1.61186578, 1.61252968, 1.61330099, 1.61451591, 1.61611502, 1.61793178,
    1.61980027, 1.62157904, 1.62315876, 1.62443226, 1.6252949 , 1.62583334,
    1.62645555, 1.62760027, 1.62957002, 1.63208393, 1.63470503, 1.63700839,
    1.63888443, 1.64056222, 1.6422873 , 1.64422197, 1.6463157 , 1.64848453,
    1.65066913, 1.65286514, 1.65507565, 1.65730495, 1.65958062, 1.66195139,
    1.66446661, 1.66714803, 1.6699525 , 1.6728276 , 1.67572691, 1.67866399,
    1.68168688, 1.68484367, 1.68816064, 1.69163081, 1.69524436, 1.69899094,
    1.70285788, 1.70683192, 1.71089965, 1.7150371 , 1.71920156, 1.72334843,
    1.72743844, 1.73147165, 1.73546566, 1.7394382 , 1.7434051 , 1.74737814,
    1.7513686 , 1.75538652, 1.75942234, 1.76345072, 1.76744589, 1.77138913,
    1.77529078, 1.7791686 , 1.78304031, 1.78692   , 1.79081551, 1.79473402,
    1.79868275, 1.80266888, 1.80669962, 1.81078217,
    }
    local Z = {
    -1.23058861, -1.22520737, -1.21958678, -1.21420037, -1.2089721 , -1.20333331,
    -1.19784428, -1.19319607, -1.18774139, -1.18215984, -1.17672523, -1.171336  ,
    -1.1659685 , -1.16060623, -1.15523274, -1.14983804, -1.14443584, -1.13904528,
    -1.1336855 , -1.12837252, -1.12308551, -1.11777961, -1.11240953, -1.10693817,
    -1.10141171, -1.09592495, -1.09057007, -1.08534602, -1.08014713, -1.07486217,
    -1.069425  , -1.06390939, -1.05841614, -1.05302462, -1.0477214 , -1.04246792,
    -1.03722574, -1.03196133, -1.02664761, -1.02125793, -1.01576636, -1.01019997,
    -1.00467492, -0.99931589, -0.99420628, -0.98925267, -0.98431435, -0.97925461,
    -0.97404232, -0.96875963, -0.96349423, -0.95830444, -0.95317355, -0.94807301,
    -0.94298211, -0.93789756, -0.9328184 , -0.92774445, -0.92269006, -0.91768278,
    -0.91275056, -0.90790634, -0.90312784, -0.8983877 , -0.89366215, -0.88896287,
    -0.88432198, -0.87977162, -0.8753329 , -0.87101007, -0.86680593, -0.86272171,
    -0.85875197, -0.85488954, -0.85112709, -0.84744285, -0.84378957, -0.84011741,
    -0.8363825 , -0.83258511, -0.82874522, -0.82488291, -0.8210157 , -0.81715554,
    -0.81331368, -0.80950011, -0.80570491, -0.80190214, -0.79806541, -0.79417556,
    -0.79024315, -0.78628637, -0.78232333, -0.77836847, -0.77442978, -0.77051462,
    -0.76663036, -0.76278437, -0.75898399, -0.7552366 ,
    }
    local idx_tlo, idx_tup
    if s < 1e-8 then
        idx_tlo = 1  -- Lua arrays are 1-indexed
        idx_tup = 1
    elseif math.abs(s-1.0) < 1e-8 then
        idx_tlo = npts
        idx_tup = npts
    else
        for i = 1, npts-1 do
            if t[i] <= s and s < t[i+1] then
                idx_tlo = i
                idx_tup = i+1
                break
            end
        end
    end
    local Dt = t[idx_tup] - t[idx_tlo]
    if idx_tlo == idx_tup then
        RZ[0] = R[idx_tlo]  -- R coordinate (0-indexed for C compatibility)
        RZ[1] = Z[idx_tlo]  -- Z coordinate (0-indexed for C compatibility)
    else
        RZ[0] = ((s-t[idx_tlo])/Dt)*R[idx_tup] + ((t[idx_tup]-s)/Dt)*R[idx_tlo]
        RZ[1] = ((s-t[idx_tlo])/Dt)*Z[idx_tup] + ((t[idx_tup]-s)/Dt)*Z[idx_tlo]
    end
end

function plate_func_lower_inner(s, RZ)
    local npts = 100
  
    local t = {
      0.        , 0.01010101, 0.02020202, 0.03030303, 0.04040404, 0.05050505,
      0.06060606, 0.07070707, 0.08080808, 0.09090909, 0.1010101 , 0.11111111,
      0.12121212, 0.13131313, 0.14141414, 0.15151515, 0.16161616, 0.17171717,
      0.18181818, 0.19191919, 0.2020202 , 0.21212121, 0.22222222, 0.23232323,
      0.24242424, 0.25252525, 0.26262626, 0.27272727, 0.28282828, 0.29292929,
      0.3030303 , 0.31313131, 0.32323232, 0.33333333, 0.34343434, 0.35353535,
      0.36363636, 0.37373737, 0.38383838, 0.39393939, 0.4040404 , 0.41414141,
      0.42424242, 0.43434343, 0.44444444, 0.45454545, 0.46464646, 0.47474747,
      0.48484848, 0.49494949, 0.50505051, 0.51515152, 0.52525253, 0.53535354,
      0.54545455, 0.55555556, 0.56565657, 0.57575758, 0.58585859, 0.5959596 ,
      0.60606061, 0.61616162, 0.62626263, 0.63636364, 0.64646465, 0.65656566,
      0.66666667, 0.67676768, 0.68686869, 0.6969697 , 0.70707071, 0.71717172,
      0.72727273, 0.73737374, 0.74747475, 0.75757576, 0.76767677, 0.77777778,
      0.78787879, 0.7979798 , 0.80808081, 0.81818182, 0.82828283, 0.83838384,
      0.84848485, 0.85858586, 0.86868687, 0.87878788, 0.88888889, 0.8989899 ,
      0.90909091, 0.91919192, 0.92929293, 0.93939394, 0.94949495, 0.95959596,
      0.96969697, 0.97979798, 0.98989899, 1.        ,
    }
    local R = {
      1.23400147, 1.23596982, 1.23692438, 1.23760647, 1.23873964, 1.24046625,
      1.24223112, 1.24357829, 1.24473493, 1.24609806, 1.24768076, 1.24927123,
      1.25072248, 1.25216181, 1.25378298, 1.25571142, 1.25800448, 1.26071586,
      1.26375963, 1.26674464, 1.26923947, 1.270992  , 1.27246039, 1.27427511,
      1.27662659, 1.27883727, 1.28012797, 1.27985427, 1.27885712, 1.2788943 ,
      1.28076395, 1.28304269, 1.28410951, 1.28383   , 1.28313288, 1.28288408,
      1.28294436, 1.28269267, 1.28188045, 1.28045658, 1.27875059, 1.27715353,
      1.27600109, 1.27507859, 1.27385662, 1.27199829, 1.26963499, 1.26696543,
      1.26414096, 1.26124713, 1.2583607 , 1.25547296, 1.25249254, 1.24935156,
      1.24611617, 1.24289388, 1.23972915, 1.23652666, 1.23317353, 1.22961774,
      1.22589572, 1.22205141, 1.21817485, 1.21443592, 1.21101198, 1.20792728,
      1.20483227, 1.20134014, 1.19741473, 1.19333087, 1.18934956, 1.18552628,
      1.18181437, 1.17817041, 1.17456403, 1.17096752, 1.16735816, 1.16372798,
      1.1600751 , 1.15645747, 1.15298346, 1.14971821, 1.14660428, 1.14356418,
      1.1405544 , 1.13757905, 1.1346467 , 1.13179194, 1.12907914, 1.12655402,
      1.12413273, 1.12168382, 1.11916745, 1.11664274, 1.11417392, 1.11183193,
      1.10969647, 1.10784783, 1.10636628, 1.10533213,
    }
    local Z = {
      -1.15621288, -1.15003049, -1.14365673, -1.13723457, -1.13090344, -1.12468749,
      -1.11847249, -1.11216563, -1.10581805, -1.09951633, -1.09326264, -1.08700972,
      -1.080725  , -1.07443785, -1.06819323, -1.06203269, -1.05599439, -1.05011605,
      -1.04437945, -1.03864409, -1.03275332, -1.02660506, -1.02031329, -1.01404485,
      -1.00788248, -1.00174305, -0.99552402, -0.98913856, -0.98267291, -0.97632012,
      -0.97016517, -0.96404269, -0.95776537, -0.95131811, -0.94480896, -0.93833879,
      -0.93189358, -0.92543579, -0.91904105, -0.91275869, -0.90652847, -0.90027248,
      -0.89392434, -0.88753274, -0.88121219, -0.87504236, -0.86901791, -0.86312083,
      -0.85731909, -0.85156118, -0.84579663, -0.84002768, -0.83430755, -0.8286761 ,
      -0.82309595, -0.81750594, -0.81188219, -0.80628355, -0.80077929, -0.7954017 ,
      -0.79012948, -0.78493671, -0.77976693, -0.77451069, -0.76905359, -0.76338046,
      -0.75771843, -0.75231866, -0.74720265, -0.74218826, -0.73710256, -0.73190848,
      -0.72663647, -0.72131539, -0.71596723, -0.71061258, -0.70526774, -0.69993643,
      -0.69461799, -0.68927613, -0.68384451, -0.67828401, -0.67263256, -0.66694068,
      -0.66123602, -0.65551422, -0.64976824, -0.6439824 , -0.63813111, -0.63219691,
      -0.6262177 , -0.62025164, -0.61431515, -0.60837962, -0.60241435, -0.59639522,
      -0.59030675, -0.584134  , -0.57786206, -0.57147602,
    }
    local idx_tlo, idx_tup
    if s < 1e-8 then
        idx_tlo = 1  -- Lua arrays are 1-indexed
        idx_tup = 1
    elseif math.abs(s-1.0) < 1e-8 then
        idx_tlo = npts
        idx_tup = npts
    else
        for i = 1, npts-1 do
            if t[i] <= s and s < t[i+1] then
                idx_tlo = i
                idx_tup = i+1
                break
            end
        end
    end
    -- Interpolate the value of R and Z.
    local Dt = t[idx_tup] - t[idx_tlo]
    if idx_tlo == idx_tup then
        RZ[0] = R[idx_tlo]  -- R coordinate (0-indexed for C compatibility)
        RZ[1] = Z[idx_tlo]  -- Z coordinate (0-indexed for C compatibility)
    else
        RZ[0] = ((s-t[idx_tlo])/Dt)*R[idx_tup] + ((t[idx_tup]-s)/Dt)*R[idx_tlo]
        RZ[1] = ((s-t[idx_tlo])/Dt)*Z[idx_tup] + ((t[idx_tup]-s)/Dt)*Z[idx_tlo]
    end
end

]]
      )
   elseif args.examples == "step" then
      io.write([[
-- Name of tokamak
name = "step"
-- Input file for tokagridgen Tool.
-- "gyrokinetic/data/eqdsk/nstxu_DN.geqdsk" for double null
--"gyrokinetic/data/eqdsk/asdex.geqdsk" for single null
geqdsk_path = "gyrokinetic/data/eqdsk/step.geqdsk"

--Number of null points. one of:
-- single_null
-- double_null
toka_type = double_null

-- Whether to generate only half the poloidal domain (DN only)
use_half_domain = false

-- Lower Private Flux region parameters
lower_PF = {
    rleft = 0.0,      -- rleft_PF
    rright = 5.2,     -- rright_PF  
    rmin = 1.7,       -- rmin_PF
    rmax = 6.2,       -- rmax_PF
    zmin_left = -8.29, -- zmin_left_PF
    zmin_right = -6.34 -- zmin_right_PF
}

-- Outer Scrape-Off Layer region parameters
outer_SOL = {
    rclose = 6.2,     -- rclose_outer_SOL
    rleft = 0.0,      -- rleft_outer_SOL
    rright = 5.2,     -- rright_outer_SOL
    rmin = 0.7,       -- rmin_outer_SOL
    rmax = 6.2,       -- rmax_outer_SOL
    zmin = -8.29,     -- zmin_outer_SOL
    zmax = 8.29       -- zmax_outer_SOL
}

-- Inner Scrape-Off Layer region parameters
inner_SOL = {
    rleft = 2.0,      -- rleft_inner_SOL
    rright = 5.2,     -- rright_inner_SOL
    rmin = 0.0,       -- rmin_inner_SOL
    rmax = 6.2,       -- rmax_inner_SOL
    zmin = -6.34,     -- zmin_inner_SOL
    zmax = 6.34       -- zmax_inner_SOL
}

-- Core region parameters
core = {
    rclose_inner = 0.0,  -- rclose_inner_core
    rclose_outer = 6.2,  -- rclose_outer_core
    rleft = 2.0,         -- rleft_core
    rright = 5.2,        -- rright_core
    rmin = 1.58,         -- rmin_core
    rmax = 6.2           -- rmax_core
}

-- Relative widths of psi for different regions (fraction of the total psi difference between center and separatrix, from 0 to 1)
-- For DN: (PF Width, Outer SOL Width, Inner SOL Width, Core Width)
-- For SN: (PF Width, SOL Width, Core Width) 
wpsi = {0.01, 0.025, 0.01, 0.04}


-- Number of cells to include along the directions of different blocks
-- (npsi_pf, npsi_sol, npsi_core, ntheta_lower_outer, ntheta_middle_SOL, ntheta_lower_inner)
cell_count = {2, 5, 8, 8, 24, 4}

-- Plate function definitions
-- These functions define the divertor plate geometry as functions of parameter s in [0,1]
function plate_func_lower_outer(s, RZ)
    local p0 = {5.488-0.6, -8.600}
    local p1 = {5.855-0.6, -8.52318}
    p1[1] = (p1[1] - p0[1])*2 + p1[1]
    p1[2] = (p1[2] - p0[2])*2 + p1[2]
    RZ[0] = (1-s)*p0[1] + s*p1[1]  -- R coordinate
    RZ[1] = (1-s)*p0[2] + s*p1[2]  -- Z coordinate
  end
  
  function plate_func_upper_outer(s, RZ)
    local p0 = {5.488-0.6, 8.600}
    local p1 = {5.855-0.6, 8.52318}
    p1[1] = (p1[1] - p0[1])*2 + p1[1]
    p1[2] = (p1[2] - p0[2])*2 + p1[2]
    RZ[0] = (1-s)*p0[1] + s*p1[1]  -- R coordinate
    RZ[1] = (1-s)*p0[2] + s*p1[2]  -- Z coordinate
  end
  
  function plate_func_upper_inner(s, RZ)
    local p0 = {1.65, 6.1}
    local p1 = {1.95, 7.2}
    p1[1] = (p1[1] - p0[1])*2 + p1[1]
    p1[2] = (p1[2] - p0[2])*2 + p1[2]
    RZ[0] = (1-s)*p0[1] + s*p1[1]  -- R coordinate
    RZ[1] = (1-s)*p0[2] + s*p1[2]  -- Z coordinate
  end
  
  function plate_func_lower_inner(s, RZ)
    local p0 = {1.65, -6.1}
    local p1 = {1.95, -7.2}
    p1[1] = (p1[1] - p0[1])*2 + p1[1]
    p1[2] = (p1[2] - p0[2])*2 + p1[2]
    RZ[0] = (1-s)*p0[1] + s*p1[1]  -- R coordinate
    RZ[1] = (1-s)*p0[2] + s*p1[2]  -- Z coordinate
  end

]]
      )
   else
      print("Error: Invalid example name. Use 'nstxu_DN', 'asdex', or 'step'")
      print("Usage: --examples nstxu_DN  (for double null example)")
      print("       --examples step      (for double null example)")
      print("       --examples asdex     (for single null example)")
      return
   end
   
   return
end

if args.template then
   if args.template == "SN" then
      io.write([[
-- Name of tokamak
name = ""
-- Input file for tokagridgen Tool.
-- "gyrokinetic/data/eqdsk/nstxu_DN.geqdsk" for double null
--"gyrokinetic/data/eqdsk/asdex.geqdsk" for single null
geqdsk_path = "gyrokinetic/data/eqdsk/"

--Number of null points. one of:
-- single_null
-- double_null
toka_type = single_null

-- Lower Private Flux region parameters
lower_PF = {
    rleft = ,      -- rleft_PF
    rright = ,     -- rright_PF  
    rmin = ,       -- rmin_PF
    rmax = ,       -- rmax_PF
    zmin = ,      -- zmin_PF
    zmax = ,      -- zmax_PF
    zmin_left = , -- zmin_left_PF
    zmin_right =  -- zmin_right_PF
}

-- Scrape-Off Layer region parameters
SOL = {
    rclose = ,     -- rclose_SOL
    rleft = ,      -- rleft_SOL
    rright = ,     -- rright_SOL
    rmin = ,       -- rmin_SOL
    rmax = ,       -- rmax_SOL
    zmin = ,      -- zmin_SOL
    zmax =         -- zmax_SOL
}

-- Core region parameters
core = {
    rclose = ,     -- rclose_core
    rleft = ,      -- rleft_core
    rright = ,     -- rright_core
    rmin = ,       -- rmin_core
    rmax = ,       -- rmax_core
    zmin = ,      -- zmin_core
    zmax =         -- zmax_core
}

-- Relative widths of psi for different regions (fraction of the total psi difference between center and separatrix, from 0 to 1)
-- For DN: (PF Width, Outer SOL Width, Inner SOL Width, Core Width)
-- For SN: (PF Width, SOL Width, Core Width) 
wpsi = { , , }

-- Number of cells to include along the directions of different blocks
-- (npsi_pf, npsi_sol, npsi_core, ntheta_lower_outer, ntheta_middle_SOL, ntheta_lower_inner)
cell_count = { , , , , , }

-- Plate function definitions
-- These functions define the divertor plate geometry as functions of parameter s in [0,1]
function plate_func_lower_outer(s, RZ)
    RZ[0] =   -- R coordinate
    RZ[1] =          -- Z coordinate
end

function plate_func_lower_inner(s, RZ)
    RZ[0] =           -- R coordinate
    RZ[1] =   -- Z coordinate
end
]]
      )
   elseif args.template == "DN" then
      io.write([[
-- Name of tokamak
name = ""
-- Input file for tokagridgen Tool.
-- "gyrokinetic/data/eqdsk/nstxu_DN.geqdsk" for double null
--"gyrokinetic/data/eqdsk/asdex.geqdsk" for single null
geqdsk_path = "gyrokinetic/data/eqdsk/"

--Number of null points. one of:
-- single_null
-- double_null
toka_type = double_null

-- Whether to generate only half the poloidal domain (DN only) (true or false)
use_half_domain = 

-- Lower Private Flux region parameters
lower_PF = {
    rleft = ,      -- rleft_PF
    rright = ,     -- rright_PF  
    rmin = ,       -- rmin_PF
    rmax = ,       -- rmax_PF
    zmin_left = , -- zmin_left_PF
    zmin_right =   -- zmin_right_PF
}

-- Outer Scrape-Off Layer region parameters
outer_SOL = {
    rclose = ,     -- rclose_outer_SOL
    rleft = ,      -- rleft_outer_SOL
    rright = ,     -- rright_outer_SOL
    rmin = ,       -- rmin_outer_SOL
    rmax = ,       -- rmax_outer_SOL
    zmin = ,      -- zmin_outer_SOL
    zmax =         -- zmax_outer_SOL
}

-- Inner Scrape-Off Layer region parameters
inner_SOL = {
    rleft = ,      -- rleft_inner_SOL
    rright = ,     -- rright_inner_SOL
    rmin = ,       -- rmin_inner_SOL
    rmax = ,       -- rmax_inner_SOL
    zmin = ,     -- zmin_inner_SOL
    zmax =        -- zmax_inner_SOL
}

-- Core region parameters
core = {
    rclose_inner = ,  -- rclose_inner_core
    rclose_outer = ,  -- rclose_outer_core
    rleft = ,         -- rleft_core
    rright = ,        -- rright_core
    rmin = ,          -- rmin_core
    rmax =            -- rmax_core
}

-- Relative widths of psi for different regions (fraction of the total psi difference between center and separatrix, from 0 to 1)
-- For DN: (PF Width, Outer SOL Width, Inner SOL Width, Core Width)
-- For SN: (PF Width, SOL Width, Core Width) 
wpsi = { , , , }

-- Number of cells to include along the directions of different blocks
-- (npsi_pf, npsi_sol, npsi_core, ntheta_lower_outer, ntheta_middle_SOL, ntheta_lower_inner)
cell_count = { , , , , , }

-- Plate function definitions
-- These functions define the divertor plate geometry as functions of parameter s in [0,1]
function plate_func_lower_outer(s, RZ)
    RZ[0] =   -- R coordinate
    RZ[1] =          -- Z coordinate
end

function plate_func_upper_outer(s, RZ)
    RZ[0] =  -- R coordinate
    RZ[1] =           -- Z coordinate
end

function plate_func_upper_inner(s, RZ)
    RZ[0] =           -- R coordinate
    RZ[1] =   -- Z coordinate
end

function plate_func_lower_inner(s, RZ)
    RZ[0] =           -- R coordinate
    RZ[1] =   -- Z coordinate
end
]]
      )
   else
      print("Error: Invalid template name. Use 'blankSN' or 'blankDN'")
      print("Usage: --template SN   (for single null blank template)")
      print("       --template DN   (for double null blank template)")
      return
   end
   return
end

-- FFI call prototypes 
ffi.cdef [[

// flag to indicate what tokamak configuration
enum gkyl_toka_grid_gen_null_points {
    GKYL_TOKA_GRID_GEN_SINGLE_NULL, // Single Null  Configuration
    GKYL_TOKA_GRID_GEN_DOUBLE_NULL, // Double Null Configuration
  }; 

struct gkylt_tokagridgen_inp {
  enum gkyl_toka_grid_gen_null_points toka_type; //Tokamak configuration

  double wpsi[4]; // Width of psi for different regions
  int cell_count[6]; // number of cells in computational space
  
  // Plate function pointers - these will be set from Lua input
  void (*plate_func_lower_outer)(double s, double* RZ);
  void (*plate_func_upper_outer)(double s, double* RZ);
  void (*plate_func_lower_inner)(double s, double* RZ);
  void (*plate_func_upper_inner)(double s, double* RZ);
  
  // Region-specific parameters (conditional based on toka_type)
  // For double_null: lower_PF[6], outer_SOL[7], inner_SOL[6], core[5]
  // For single_null: lower_PF[8], SOL[7], core[8]
  double lower_PF[8]; // lower private flux region parameters (max size for both types)
  double outer_SOL[7]; // outer scrape-off layer parameters (DN only)
  double inner_SOL[6]; // inner scrape-off layer parameters (DN only)
  double SOL[7]; // scrape-off layer parameters (SN only)
  double core[7]; // core region parameters (max size for both types)
  
  const char *geqdsk_path; // name for file with psi(R,Z)
  const char *out_prefix; // output prefix
  const char *name; // application name
  bool use_half_domain; // flag to indicate if half domain should be used (DN only)

};

// Settings struct to return metadata
struct gkylt_tokagridgen_settings {
  int toka_type;
  bool use_half_domain;
  double psisep;
};

// Return struct containing both bgeom and settings
struct gkylt_tokagridgen_output {
  struct gkyl_gk_block_geom* bgeom;
  struct gkylt_tokagridgen_settings settings;
};

struct gkylt_tokagridgen_output gkylt_tokagridgen(const struct gkylt_tokagridgen_inp *tginp);


// Lite struct and accessor to safely read nested geometry info from Lua
struct gkyl_gk_block_geom_info_lite {
  double lower[3], upper[3];
  int cells[3], cuts[3];

  int geometry_id;

  struct {
    char filepath[1024];
    int rz_poly_order;
    int flux_poly_order;
    int reflect;
    int use_gpu;
  } efit_info;

  struct {
    int ftype;
    int half_domain;
    double rclose, rleft, rright, rmin, rmax;
    double zmin, zmax, zmin_left, zmin_right, zmax_left, zmax_right;
    int plate_spec;
    int inexact_roots;
    int use_cubics;
    int use_hyperbolic_numbers;
  } tok_grid_info;

  struct { int bid; int dir; int edge; } connections[3][2];
};

void gkyl_gk_block_geom_get_block_lite(const struct gkyl_gk_block_geom *bgeom,
  int bidx, struct gkyl_gk_block_geom_info_lite *out);

// Block geometry access functions
int gkyl_gk_block_geom_ndim(const struct gkyl_gk_block_geom *bgeom);
int gkyl_gk_block_geom_num_blocks(const struct gkyl_gk_block_geom *bgeom);
]]

single_null = ffi.C.GKYL_TOKA_GRID_GEN_SINGLE_NULL
double_null = ffi.C.GKYL_TOKA_GRID_GEN_DOUBLE_NULL
-- open input file and read contents: this loads 'speciesList' and
-- 'field' into this module's namespace
if args.input then
   local inpFile = assert(loadfile(args.input))
   inpFile() -- load contents of input file
else
   print("Must specify an input file to run!")
   return
end

-- set input struct ...
local inp = ffi.new("struct gkylt_tokagridgen_inp")
inp.toka_type = toka_type
inp.geqdsk_path = geqdsk_path
inp.out_prefix = outPrefix
inp.name = name
-- Default behavior: full domain unless specified; DN may want half
if use_half_domain ~= nil then
   inp.use_half_domain = use_half_domain
else
   inp.use_half_domain = false
end

-- Set plate function pointers from Lua input
if plate_func_lower_outer then
   inp.plate_func_lower_outer = plate_func_lower_outer
end
if plate_func_upper_outer then
   inp.plate_func_upper_outer = plate_func_upper_outer
end
if plate_func_lower_inner then
   inp.plate_func_lower_inner = plate_func_lower_inner
end
if plate_func_upper_inner then
   inp.plate_func_upper_inner = plate_func_upper_inner
end

for d = 1, #cell_count do
   inp.cell_count[d-1] = cell_count[d]
end
for d = 1, #wpsi do
    inp.wpsi[d-1] = wpsi[d]
 end

 -- Initialize all arrays to zero
for d = 0, 7 do
   inp.lower_PF[d] = 0.0
end
for d = 0, 6 do
   inp.outer_SOL[d] = 0.0
   inp.SOL[d] = 0.0
   inp.core[d] = 0.0
end
for d = 0, 5 do
   inp.inner_SOL[d] = 0.0
end

-- Handle region-specific parameters based on toka_type
if toka_type == double_null then
   -- Double null configuration: lower_PF[6], outer_SOL[7], inner_SOL[6], core[6]
   if lower_PF then
      -- Convert dictionary to array format for C interface
      if type(lower_PF) == "table" and lower_PF.rleft then
         -- Dictionary format: convert to array
         inp.lower_PF[0] = lower_PF.rleft or 0.0
         inp.lower_PF[1] = lower_PF.rright or 0.0
         inp.lower_PF[2] = lower_PF.rmin or 0.0
         inp.lower_PF[3] = lower_PF.rmax or 0.0
         inp.lower_PF[4] = lower_PF.zmin_left or 0.0
         inp.lower_PF[5] = lower_PF.zmin_right or 0.0
      else
         -- Array format: use existing logic
      for d = 1, math.min(#lower_PF, 6) do
         inp.lower_PF[d-1] = lower_PF[d]
         end
      end
   end
   
   if outer_SOL then
      -- Convert dictionary to array format for C interface
      if type(outer_SOL) == "table" and outer_SOL.rclose then
         -- Dictionary format: convert to array
         inp.outer_SOL[0] = outer_SOL.rclose or 0.0
         inp.outer_SOL[1] = outer_SOL.rleft or 0.0
         inp.outer_SOL[2] = outer_SOL.rright or 0.0
         inp.outer_SOL[3] = outer_SOL.rmin or 0.0
         inp.outer_SOL[4] = outer_SOL.rmax or 0.0
         inp.outer_SOL[5] = outer_SOL.zmin or 0.0
         inp.outer_SOL[6] = outer_SOL.zmax or 0.0
      else
         -- Array format: use existing logic
      for d = 1, math.min(#outer_SOL, 7) do
         inp.outer_SOL[d-1] = outer_SOL[d]
         end
      end
   end
   
   if inner_SOL then
      -- Convert dictionary to array format for C interface
      if type(inner_SOL) == "table" and inner_SOL.rleft then
         -- Dictionary format: convert to array
         inp.inner_SOL[0] = inner_SOL.rleft or 0.0
         inp.inner_SOL[1] = inner_SOL.rright or 0.0
         inp.inner_SOL[2] = inner_SOL.rmin or 0.0
         inp.inner_SOL[3] = inner_SOL.rmax or 0.0
         inp.inner_SOL[4] = inner_SOL.zmin or 0.0
         inp.inner_SOL[5] = inner_SOL.zmax or 0.0
      else
         -- Array format: use existing logic
      for d = 1, math.min(#inner_SOL, 6) do
         inp.inner_SOL[d-1] = inner_SOL[d]
         end
      end
   end
   
   if core then
      -- Convert dictionary to array format for C interface
      if type(core) == "table" and core.rclose_inner then
         -- Dictionary format: convert to array
         inp.core[0] = core.rclose_inner or 0.0
         inp.core[1] = core.rclose_outer or 0.0
         inp.core[2] = core.rleft or 0.0
         inp.core[3] = core.rright or 0.0
         inp.core[4] = core.rmin or 0.0
         inp.core[5] = core.rmax or 0.0
      else
         -- Array format: use existing logic
      for d = 1, math.min(#core, 6) do
         inp.core[d-1] = core[d]
         end
      end
   end
   
elseif toka_type == single_null then
   -- Single null configuration: lower_PF[8], SOL[7], core[8]
   if lower_PF then
      -- Convert dictionary to array format for C interface
      if type(lower_PF) == "table" and lower_PF.rleft then
         -- Dictionary format: convert to array
         inp.lower_PF[0] = lower_PF.rleft or 0.0
         inp.lower_PF[1] = lower_PF.rright or 0.0
         inp.lower_PF[2] = lower_PF.rmin or 0.0
         inp.lower_PF[3] = lower_PF.rmax or 0.0
         inp.lower_PF[4] = lower_PF.zmin or 0.0
         inp.lower_PF[5] = lower_PF.zmax or 0.0
         inp.lower_PF[6] = lower_PF.zmin_left or 0.0
         inp.lower_PF[7] = lower_PF.zmin_right or 0.0
      else
         -- Array format: use existing logic
      for d = 1, math.min(#lower_PF, 8) do
         inp.lower_PF[d-1] = lower_PF[d]
         end
      end
   end
   
   if SOL then
      -- Convert dictionary to array format for C interface
      if type(SOL) == "table" and SOL.rclose then
         -- Dictionary format: convert to array
         inp.SOL[0] = SOL.rclose or 0.0
         inp.SOL[1] = SOL.rleft or 0.0
         inp.SOL[2] = SOL.rright or 0.0
         inp.SOL[3] = SOL.rmin or 0.0
         inp.SOL[4] = SOL.rmax or 0.0
         inp.SOL[5] = SOL.zmin or 0.0
         inp.SOL[6] = SOL.zmax or 0.0
      else
         -- Array format: use existing logic
      for d = 1, math.min(#SOL, 7) do
         inp.SOL[d-1] = SOL[d]
         end
      end
   end
   
   if core then
      -- Convert dictionary to array format for C interface
      if type(core) == "table" and core.rclose then
         -- Dictionary format: convert to array
         inp.core[0] = core.rclose or 0.0
         inp.core[1] = core.rleft or 0.0
         inp.core[2] = core.rright or 0.0
         inp.core[3] = core.rmin or 0.0
         inp.core[4] = core.rmax or 0.0
         inp.core[5] = core.zmin or 0.0
         inp.core[6] = core.zmax or 0.0
      else
         -- Array format: use existing logic
      for d = 1, math.min(#core, 7) do
         inp.core[d-1] = core[d]
         end
      end
   end
end

-- Create output directory named <name>_geometry
local output_dir = name .. "_geometry"
lfs.mkdir(output_dir)

-- Change to the output directory
local original_dir = lfs.currentdir()
lfs.chdir(output_dir)

-- Adjust geqdsk_path to point to parent directory since we're now in a subdirectory
inp.geqdsk_path = "../" .. geqdsk_path
 
local tmStart = Time.clock()
-- generate grid with settings
local result = ffi.C.gkylt_tokagridgen(inp)
local bgeom = result.bgeom
local settings = result.settings

-- Extract block geometry data for metadata
local function extractBlockGeometryData(bgeom)
   local ndim = ffi.C.gkyl_gk_block_geom_ndim(bgeom)
   local num_blocks = ffi.C.gkyl_gk_block_geom_num_blocks(bgeom)
  
  -- Enum normalizers
  local ftype_names = {
    [0] = "GKYL_DN_SOL_OUT",
    [1] = "GKYL_DN_SOL_IN",
    [2] = "GKYL_LSN_SOL",
    [3] = "GKYL_USN_UP",
    [4] = "GKYL_CORE",
    [5] = "GKYL_DN_SOL_OUT_LO",
    [6] = "GKYL_DN_SOL_OUT_MID",
    [7] = "GKYL_DN_SOL_OUT_UP",
    [8] = "GKYL_DN_SOL_IN_LO",
    [9] = "GKYL_DN_SOL_IN_MID",
    [10] = "GKYL_DN_SOL_IN_UP",
    [11] = "GKYL_LSN_SOL_LO",
    [12] = "GKYL_LSN_SOL_MID",
    [13] = "GKYL_LSN_SOL_UP",
    [14] = "GKYL_PF_UP_L",
    [15] = "GKYL_PF_UP_R",
    [16] = "GKYL_PF_LO_L",
    [17] = "GKYL_PF_LO_R",
    [18] = "GKYL_CORE_L",
    [19] = "GKYL_CORE_R",
    [20] = "GKYL_IWL",
  }
  local edge_names = {
    [1] = "GKYL_LOWER_POSITIVE",
    [2] = "GKYL_LOWER_NEGATIVE",
    [3] = "GKYL_UPPER_POSITIVE",
    [4] = "GKYL_UPPER_NEGATIVE",
    [5] = "GKYL_PHYSICAL",
  }
  local geometry_id_names = {
    [0] = "GKYL_TOKAMAK",
  }

  -- Sets of ftypes for applicability rules
  local PF_UP = { [14]=true, [15]=true }
  local PF_LO = { [16]=true, [17]=true }
  local DN_SOL_OUT = { [5]=true, [6]=true, [7]=true }
  local DN_SOL_IN  = { [8]=true, [9]=true, [10]=true }
  local LSN_SOL    = { [11]=true, [12]=true, [13]=true }
  local CORE_SPLIT = { [18]=true, [19]=true }

  local blocks_data = {}
  local topo_connections = {}
   
   -- Extract data for each block
   for bidx = 0, num_blocks-1 do
     local lite = ffi.new("struct gkyl_gk_block_geom_info_lite[1]")
     ffi.C.gkyl_gk_block_geom_get_block_lite(bgeom, bidx, lite)
     local info = lite[0]

     local block_data = {
         lower = {},
         upper = {},
         cells = {},
         cuts = {},
         connections = {}
      }
      
      -- Extract extents, cells, and cuts
      for d = 0, ndim-1 do
        block_data.lower[d+1] = info.lower[d]
        block_data.upper[d+1] = info.upper[d]
        block_data.cells[d+1] = info.cells[d]
        block_data.cuts[d+1] = info.cuts[d]
         
         -- Extract connections for each direction and edge
         block_data.connections[d+1] = {}
         for e = 0, 1 do -- 0: lower edge, 1: upper edge
           local conn = info.connections[d][e]
            block_data.connections[d+1][e+1] = {
              bid = conn.bid,
              dir = conn.dir,
              edge = conn.edge,
              edge_name = edge_names[conn.edge]
            }
         end
      end
      
      -- Extract geometry information from lite struct
      local ftype = info.tok_grid_info.ftype -- numeric enum
      local geom = {
        geometry_id = info.geometry_id,
        geometry_id_name = geometry_id_names[info.geometry_id],
        efit_info = {
           filepath = ffi.string(info.efit_info.filepath),
           rz_poly_order = info.efit_info.rz_poly_order,
           flux_poly_order = info.efit_info.flux_poly_order,
           reflect = info.efit_info.reflect ~= 0,
           use_gpu = info.efit_info.use_gpu ~= 0
        },
        tok_grid_info = {
           ftype = ftype_names[ftype],
           half_domain = info.tok_grid_info.half_domain ~= 0,
           rclose = info.tok_grid_info.rclose,
           rleft = info.tok_grid_info.rleft,
           rright = info.tok_grid_info.rright,
           rmin = info.tok_grid_info.rmin,
           rmax = info.tok_grid_info.rmax,
           zmin = info.tok_grid_info.zmin,
           zmax = info.tok_grid_info.zmax,
           zmin_left = info.tok_grid_info.zmin_left,
           zmin_right = info.tok_grid_info.zmin_right,
           zmax_left = info.tok_grid_info.zmax_left,
           zmax_right = info.tok_grid_info.zmax_right,
           plate_spec = info.tok_grid_info.plate_spec ~= 0,
           inexact_roots = info.tok_grid_info.inexact_roots ~= 0,
           use_cubics = info.tok_grid_info.use_cubics ~= 0,
           use_hyperbolic_numbers = info.tok_grid_info.use_hyperbolic_numbers ~= 0
        }
      }

      -- ftype-based omission of inapplicable fields (remove rclose, zmin_left, etc from blocks where that isn't a parameter)
      local tgi = geom.tok_grid_info
      if not (DN_SOL_OUT[ftype] or LSN_SOL[ftype] or CORE_SPLIT[ftype]) then
        tgi.rclose = nil
      end
      if not (PF_LO[ftype] or PF_UP[ftype]) then
        tgi.zmin_left = nil; tgi.zmin_right = nil; tgi.zmax_left = nil; tgi.zmax_right = nil
      else
        if PF_LO[ftype] then
          -- lower PF uses zmin_left/right; zmax_* not applicable
          tgi.zmax_left = nil; tgi.zmax_right = nil
        end
        if PF_UP[ftype] then
          -- upper PF uses zmax_left/right; zmin_* not applicable
          tgi.zmin_left = nil; tgi.zmin_right = nil
        end
      end

      block_data.geometry = geom
      -- duplicate block info connections into btopo 
      topo_connections[bidx+1] = block_data.connections
      
      
      blocks_data[bidx+1] = block_data
   end
   
   return {
      ndim = ndim,
      num_blocks = num_blocks,
      blocks = blocks_data,
      btopo = {
        ndim = ndim,
        num_blocks = num_blocks,
        connections = topo_connections,
      },
   }
end

-- Extract the geometry data
local geometry_metadata = extractBlockGeometryData(bgeom)

-- Function to extract settings metadata from the returned settings struct
local function extractSettingsMetadata(settings)
   local toka_type_names = {
      [0] = "GKYL_TOKA_GRID_GEN_SINGLE_NULL",
      [1] = "GKYL_TOKA_GRID_GEN_DOUBLE_NULL"
   }
   
   return {
      toka_type = toka_type_names[settings.toka_type],
      use_half_domain = settings.use_half_domain,
      psisep = settings.psisep
   }
end

-- Extract the settings metadata
local settings_metadata = extractSettingsMetadata(settings)

-- Set GKYL_OUT_PREFIX to empty string to avoid path concatenation
GKYL_OUT_PREFIX = name

-- Create a simple grid for writing the metadata
local grid1d = Grid.RectCart {
   lower = { 0.0 },
   upper = { 1.0 },
   cells = { 1 }
}

-- Create a field with the geometry metadata
local geometry_field = DataStruct.Field {
   onGrid = grid1d,
   numComponents = 1,
   metaData = geometry_metadata,
}

-- Create a field with the settings metadata
local settings_field = DataStruct.Field {
   onGrid = grid1d,
   numComponents = 1,
   metaData = settings_metadata,
}

-- Write the geometry metadata to file with a simple filename
-- DataStruct automatically prepends the name, so we just need "_bgeom.gkyl"
local output_filename = "bgeom.gkyl"
geometry_field:write(output_filename)

-- Write the settings metadata to file
local settings_filename = "settings.gkyl"
settings_field:write(settings_filename)

-- Change back to original directory
lfs.chdir(original_dir)

io.write(string.format("Completed in %g seconds\n", Time.clock()-tmStart))

