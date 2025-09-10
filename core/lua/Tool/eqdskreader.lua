-- Gkyl ------------------------------------------------------------------------
--
-- Reads eqdsk files and converts to Gkeyll files
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Alloc = require "Lib.Alloc"
local DataStruct = require "DataStruct"
local ZeroDynVector = require "DataStruct.ZeroDynVector"
local Grid = require "Grid"
local argparse = require "Lib.argparse"
local Proto = require "Lib.Proto"
local Mpi = require "Comm.Mpi"

-- Create CLI parser to handle commands and options
local parser = argparse()
   :name("eqdskreader")
   :description [[
Reads eqdsk files and converts the data into Gkyell file so they can be
read in gkeyll and postgkyl. The eqdsk file format is somewhat strange
and it is not clear if this tool can read every eqdsk file out
there. Some caution is needed.
]]

-- add flags and options
parser:option("-i --input", "Input eqdsk file to read")

-- parse command line parameters
local args = parser:parse(GKYL_COMMANDS_L)

-- modify prefix to write files based on eqdsk file name
GKYL_OUT_PREFIX = lfs.currentdir() .. "/" .. args.input

-- Reads and returns header information
local function readHeader(fh)
   local ln1 = fh:read()
   local n = ln1:len()
   
   local s,e = ln1:find("%d+")
   local nx = tonumber(ln1:sub(s,e))

   ln1 = ln1:sub(e+1, n)
   local s,e = ln1:find("%d+")
   local ny = tonumber(ln1:sub(s,e))

   return nx, ny
end

-- reads eqdsk file all into memory at once. Further splitting into
-- pieces is not done in this function.
--
-- fh : Handle to open file
-- returns nx, ny, data
local function readEqdsk(fh)
   local nx, ny = readHeader(fh)
   local data = Alloc.Double()

   numMatch = "[+-]?%d+%.%d+[eE][+-]%d+"
   intNumMatch = "[+-]?%d+"
   for ln in fh:lines() do
      local n = ln:len()
      local ln1 = ln:sub(1, n) -- copy so we can mod
   
      while true do
	 local s,e = ln1:find(numMatch)
	 local ints,inte = ln1:find(intNumMatch)
	 if s ~= nil and e ~= nil then
	    data:push( tonumber(ln1:sub(s,e)) )
	    ln1 = ln1:sub(e+1, n)
	 elseif ints ~= nil and inte ~= nil then
	    data:push( tonumber(ln1:sub(ints,inte)) )
	    ln1 = ln1:sub(inte+1, n)
	 else
	    break
	 end
      end
   end

   return nx, ny, data
end

-- open file and read data
fh = io.open(args.input, "r")
local nx, ny, data = readEqdsk(fh)
fh:close()

-- function to allocate fields
local function makeField(meta, grid)
   return DataStruct.Field {
      onGrid = grid,
      numComponents = 1,
      metaData = meta,
   }
end

-- 1D iterator
local Iter1D = Proto()
function Iter1D:init(arr)
   self.arr = arr
   self.currIdx = 1
end
function Iter1D:seek(idx)
   self.currIdx = idx
end
function Iter1D:next()
   local val = 0
   if self.currIdx <= self.arr:size() then
      val = self.arr[self.currIdx]
   end
   self.currIdx = self.currIdx+1
   return val
end


dataItr = Iter1D(data) -- so we can read number by number
-- extract various pieces of from read data
local rdim, zdim = dataItr:next(), dataItr:next()
local rcenter, rleft = dataItr:next(), dataItr:next()
local zmid = dataItr:next()

local rmaxis, zmaxis = dataItr:next(), dataItr:next()
local simag, sibry = dataItr:next(), dataItr:next()
local bcenter = dataItr:next()
local current = dataItr:next()

-- store scalars in a table to put into meta-data file
local metaData = {
   rdim = rdim, zdim = zdim,
   rcenter = rcenter, rleft = rleft,
   zmid = zmid,

   rmaxis = rmaxis, zmaxis = zmaxis,
   simag = simag, sibry = sibry,
   bcenter = bcenter,

   current = current,
}

-- create grids
local grid1d = Grid.RectCart {
   lower = { simag },
   upper = { sibry },
   cells = { nx }
}
local grid2d = Grid.RectCart {
   lower = { rleft, zmid-zdim/2},
   upper = { rleft+rdim, zmid+zdim/2 },
   cells = { nx, ny }
}

-- allocate various fields
fpol = makeField(metaData, grid1d)
pres = makeField(metaData, grid1d)
ffprime = makeField(metaData, grid1d)
pprime = makeField(metaData, grid1d)
psi = makeField(metaData, grid2d)
qpsi = makeField(metaData, grid1d)

dataItr:seek(21) -- need to skip forward
-- construct 1D arrays
for i = 1, nx do
   fpol:get(i)[1] = dataItr:next()
end
for i = 1, nx do
   pres:get(i)[1] = dataItr:next()
end
for i = 1, nx do
   ffprime:get(i)[1] = dataItr:next()
end
for i = 1, nx do
   pprime:get(i)[1] = dataItr:next()
end

local indexer = psi:indexer()
for j = 1, ny do
   for i = 1, nx do
      psi:get(indexer(i,j))[1] = dataItr:next()
   end
end

for i = 1, nx do
   qpsi:get(i)[1] = dataItr:next()
end

-- read (R,Z) coordinates of LCFS
local sepRZ = ZeroDynVector.DynVector { numComponents = 2 }

local nbbbs = dataItr:next()
local limtr = dataItr:next()

for i = 1, nbbbs do
   sepRZ:appendData(i, { dataItr:next(), dataItr:next() } )
end

-- write data to Gkeyll files
fpol:write("fpol.gkyl")
pres:write("pres.gkyl")
ffprime:write("ffprime.gkyl")
pprime:write("pprime.gkyl")
psi:write("psi.gkyl")
qpsi:write("qpsi.gkyl")
sepRZ:write("sepRZ.gkyl")
