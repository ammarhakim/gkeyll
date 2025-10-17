-- Gkyl ------------------------------------------------------------------------
--
-- Generate grid for use in simulations of magnetic mirrors
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local argparse = require "Lib.argparse"
local ffi = require "ffi"
local lfs = require "Lib.lfs"
local xsys = require "xsys"
local Time = require "Lib.Time"

-- Create CLI parser to handle commands and options
local parser = argparse()
   :name("mirrorgridgen")
   :description [[

Generates a grid for use in mirror simulations. This tool takes as
input the solution to the Grad-Shafranov equation and constructs nodes
the RZ plane along with a minimum set of geometric quantities needed
in various solvers.

]]

-- add flags and options
parser:flag("-e --example", "Fetch example input file", false)
parser:option("-i --input", "Input file")

-- parse command line parameters
local args = parser:parse(GKYL_COMMANDS_L)

if args.example then
   -- write out example to terminal
   io.write([[
-- Input file for mirrorgridgen Tool.

psiRZ_fname = "wham_hires.geqdsk_psi.gkyl"
include_axis = true -- should we include the r=0 axis?
write_psi_cubic = true -- write the bicubic interpolation to psi

-- field-line coordinate to use. one of:
-- sqrt_psi_cart_z
-- psi_cart_z
field_line_coordinate = sqrt_psi_cart_z

-- lower and upper extents of computational space grid
-- (psi, phi, z)
lower = { 2.0e-6, 0.0, -2.0 }
upper = { 1.0e-3, 2*math.pi, 2.0 }
cells = { 10, 16, 64 }

]]
   )
   
   return
end

-- FFI call prototypes 
ffi.cdef [[

// flag to indicate what field-line coordinate to use
enum gkyl_mirror_grid_gen_field_line_coord {
  GKYL_MIRROR_GRID_GEN_PSI_CART_Z, // use psi and Cartesian Z coordinate
  GKYL_MIRROR_GRID_GEN_SQRT_PSI_CART_Z, // use sqrt(psi) and Cartesian Z coordinate
  GKYL_MIRROR_GRID_GEN_FL_LENGTH, // use field-line length NYI
};  

struct gkylt_mirrorgridgen_inp {
  enum gkyl_mirror_grid_gen_field_line_coord fl_coord; // field-line coordinate to use
  bool include_axis; // add nodes on r=0 axis (the axis is assumed be psi=0)

  double lower[3], upper[3]; // lower and upper bounds of computational space
  int cells[3]; // number of cells in computational space
  
  bool write_psi_cubic; // set to true to write the cubic fit to file
  const char *psiRZ_fname; // name for file with psi(R,Z)
  const char *out_prefix; // output prefix
};

bool gkylt_mirrorgridgen(const struct gkylt_mirrorgridgen_inp *mginp);
]]

psi_cart_z = ffi.C.GKYL_MIRROR_GRID_GEN_PSI_CART_Z
sqrt_psi_cart_z = ffi.C.GKYL_MIRROR_GRID_GEN_SQRT_PSI_CART_Z

-- open input file and read contents: this loads 'speciesList' and
-- 'field' into this module's namespace
if args.input then
   local inpFile = assert(loadfile(args.input))
   inpFile() -- load contents of input file
else
   print("Must specify an input file to run!")
   return
end

-- compute output prefix
local outPrefix = lfs.currentdir() .. "/" .. args.input:sub(1, args.input:len()-4)

-- set input struct ...
local inp = ffi.new("struct gkylt_mirrorgridgen_inp")
inp.fl_coord = field_line_coordinate
inp.include_axis = include_axis 
inp.write_psi_cubic = write_psi_cubic
inp.psiRZ_fname = psiRZ_fname
inp.out_prefix = outPrefix

for d = 1, #lower do
   inp.lower[d-1] = lower[d]
end
for d = 1, #upper do
   inp.upper[d-1] = upper[d]
end
for d = 1, #cells do
   inp.cells[d-1] = cells[d]
end

local tmStart = Time.clock()
-- generate grid
ffi.C.gkylt_mirrorgridgen(inp)
io.write(string.format("Completed in %g seconds\n", Time.clock()-tmStart))
