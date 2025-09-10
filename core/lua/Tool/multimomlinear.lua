-- Gkyl ------------------------------------------------------------------------
--
-- Linear dispersion solver for multi-moment multifluid equations.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local DataStruct = require "DataStruct"
local ZeroDynVector = require "DataStruct.ZeroDynVector"
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"
local Species = require "Tool.LinearSpecies"
local argparse = require "Lib.argparse"
local complex = require "sci.complex"
local ffi = require "ffi"
local lfs = require "Lib.lfs"
local xsys = require "xsys"
local Time = require "Lib.Time"
local Mpi = require "Comm.Mpi"

-- Create CLI parser to handle commands and options
local parser = argparse()
   :name("multimomlinear")
   :description [[
Solves the linear dispersion relation for multi-moment multifluid
equations. Arbitrary number of species are supported. Each plasma
species can be either an isothermal fluid (including cold fluid),
five-moment or a ten-moment fluid. Some support for Hammett-Perkin
type collisionless closures is provided. The field equations can
either be full Maxwell equations or, for electrostatic problems,
Poisson equations.

For full details on how the dispersion solver works see Developer
Documentation on Gkeyll RTD website.
]]

-- add flags and options
parser:flag("-e --example", "Fetch example input file", false)
parser:option("-i --input", "Input file")

-- parse command line parameters
local args = parser:parse(GKYL_COMMANDS_L)

if args.example then
   -- write out example to terminal
   io.write([[
-- Input file for multimomlinear Tool.

local Species = require "Tool.LinearSpecies"

-- Electrons
elc = Species.Euler {
   mass = 1.0, -- mass
   charge = -1.0, -- charge
   density = 1.0, -- number density
   velocity = {1.0, 0.0, 0.0}, -- velocity vector
   pressure = 0.1, -- pressure
}

-- Ions
ion = Species.Euler {
   mass = 25.0, -- mass
   charge = 1.0, -- charge
   density = 1.0, -- number density
   velocity = {0.0, 0.0, 0.0}, -- velocity vector
   pressure = 0.1, -- pressure
}

-- EM field
field = Species.Maxwell {
   epsilon0 = 1.0, mu0 = 1.0,

   electricField = {0.0, 0.0, 0.0}, -- background electric field
   magneticField = {0.0, 0.0, 1.0}, -- background magnetic field
}

-- list of species to include in dispersion relation
speciesList = { elc, ion }

-- Wave vectors for which to compute dispersion 
kvectors = { {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0} }

-- flag to indicate if eigenvectors should be computed
calcEigenVectors = false
]]
   )
   
   return
end

-- FFI call prototypes 
ffi.cdef [[
  bool gkyl_multi_mom_eigensolve(double _Complex *A, double _Complex *x, double _Complex *vl, double _Complex *vr, int N, int eig_vec);
]]

-- clear out contents of matrix
local function matrixClear(m, val)
   for i = 1, m:numRows() do
      for j = 1, m:numCols() do
	 m[i][j] = val
      end
   end
end

-- Increment src matrix in sub-region by dest matrix starting from
-- row,col location
local function matrixSubIncr(dest, src, row, col)
   for r = 1, src:numRows() do
      for c = 1, src:numCols() do
	 dest[r+row-1][c+col-1] = dest[r+row-1][c+col-1] + src[r][c]
      end
   end
end

-- Computes eigensystem
local function computeEigenSystem(A, N, calcVec)
   local ere, eim = Lin.Vec(N), Lin.Vec(N)
   local ecomp = Lin.ComplexVec(N)
   local vre, vim = Lin.ComplexMat(N, N), Lin.ComplexMat(N, N)

   local info = ffi.C.gkyl_multi_mom_eigensolve(A:data(), ecomp:data(), vre:data(), vim:data(), N, calcVec == true and 1 or 0)

   for i = 1, N do
      ere[i], eim[i] = ecomp[i].re, ecomp[i].im
   end

   return ere, eim
end

local count = 1 -- just a counter for number of wave-vectors

-- Solves the dispersion relation in the EM case. Inputs are the wave
-- vector, the list of species objects and the field
-- object. Eigenvalues corresponding to kvec are appended to the
-- eigValues DynVector.
local function solveDispEM(kvec, speciesList, field, eigValues, calcEigVec)
   -- total number of equations 
   local numEqn = field:numEquations()
   for _, s in ipairs(speciesList) do
      numEqn = numEqn + s:numEquations()
   end
   -- Allocate space for dispersion matrix
   local dispMat = Lin.ComplexMat(numEqn, numEqn)
   matrixClear(dispMat, 0.0)

   -- compute starting index into each species sub-matrices in dispMat
   local speciesIdx = {}
   local currStart = 1
   for sidx, s in ipairs(speciesList) do
      speciesIdx[sidx] = currStart
      currStart = currStart+s:numEquations()
   end
   local fieldIdx = currStart -- starting index into EM field

   -- loop over each species, inserting contribution to dispersion matrix
   for sidx, s in ipairs(speciesList) do
      -- species' contribution to moment equations
      local Dmom = s:calcMomDispMat(kvec, field)
      matrixSubIncr(dispMat, Dmom, speciesIdx[sidx], speciesIdx[sidx])

      -- species' contribution to field equations
      local Dfld = s:calcFieldDispMat(kvec, field)
      matrixSubIncr(dispMat, Dfld, speciesIdx[sidx], fieldIdx)

      if field:isa(Species.Poisson) then
	 -- for ES problems we need an additional contribution to this
	 -- species from every other species, including itself
	 for spidx, sp in ipairs(speciesList) do
	    local Des = s:calcElectrostaticDispMat(kvec, sp, field)
	    matrixSubIncr(dispMat, Des, speciesIdx[sidx], speciesIdx[spidx])
	 end
      end
   end

   -- get field contribution to field equation
   local Dfld = field:calcFieldDispMat(kvec)
   matrixSubIncr(dispMat, Dfld, fieldIdx, fieldIdx)

   -- for each species compute contribution to moment equations (this is
   -- from current sources.)
   for sidx, s in ipairs(speciesList) do
      local Dmom = field:calcMomDispMat(kvec, s.charge, s.density, s.velocity)
      matrixSubIncr(dispMat, Dmom, fieldIdx, speciesIdx[sidx])
   end

   -- create eigenvalue solver and solve system
   local evr, evi = computeEigenSystem(dispMat, dispMat:numRows(), calcEigVec)

   -- append to output field
   local kw = { kvec[1], kvec[2], kvec[3] }
   for i = 1, #evr do
      kw[3+2*i-1], kw[3+2*i] = evr[i], evi[i]
   end
   eigValues:appendData(count, kw )
   count = count+1

   --deleteEigenSolver(eigSolver)
end

-- open input file and read contents: this loads 'speciesList' and
-- 'field' into this module's namespace
if args.input then
   local inpFile = assert(loadfile(args.input))
   inpFile() -- load contents of input file
else
   print("Must specify an input file to run!")
   return
end

assert(speciesList, "Input file must have a 'speciesList' table!")
assert(field, "Input file must have a 'field'!")

-- check if we should compute eigenvectors
local calcEigVec = xsys.pickBool(calcEigenVectors, false)

-- print some info about run
io.write(string.format("Number of species %d:", #speciesList))
for _, s in ipairs(speciesList) do
   io.write(string.format(" %s", s:type()))
end
io.write("\n")
io.write(string.format("Field type is %s\n", field:type()))
io.write(string.format("Computing dispersion relation for %d wave-vectors\n", #kvectors))
if calcEigVec then
   io.write("Eigenvectors will be computed\n")
else
   io.write("Eigenvectors will NOT be computed\n")
end

-- modify the output prefix to ensure files are written to the proper
-- place
GKYL_OUT_PREFIX = lfs.currentdir() .. "/" .. args.input:sub(1, args.input:len()-4)

-- total number of equations 
local numEqn = field:numEquations()
for _, s in ipairs(speciesList) do
   numEqn = numEqn + s:numEquations()
end
-- allocate dynvector for output. Arranged as:
-- kx, ky, kz, w1r, w1i, w2r, w2i, ...
local eigValues = ZeroDynVector.DynVector { numComponents = 3 + 2*numEqn }

local tmStart = Time.clock()

-- solve dispersion relation for each kvector
for _, kv in ipairs(kvectors) do
   solveDispEM(kv, speciesList, field, eigValues, calcEigVec)
end

-- write out eigenvaules
eigValues:write("frequencies.gkyl", 0.0, 0)

io.write(string.format("Completed in %g seconds\n", Time.clock()-tmStart))


