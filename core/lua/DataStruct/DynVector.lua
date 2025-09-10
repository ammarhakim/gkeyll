-- Gkyl ------------------------------------------------------------------------
--
-- Dynamically growing 1D vector.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- System libraries.
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- Gkyl libraries.
local Proto = require "Lib.Proto"
local Alloc = require "Lib.Alloc"
local Lin = require "Lib.Linalg"
local Mpi = require "Comm.Mpi"

-- Template to copy from table/vector.
local copyTempl = xsys.template [[
return function (src, dest)
| for i = 1, NCOMP do
  dest[${i}] = src[${i}]
| end
end 	    
]]

-- DynVector -------------------------------------------------------------------
--
-- Dynamically growing 1D vector. Used to store diagnostics
--------------------------------------------------------------------------------

local DynVector = Proto()

-- Constructor for DynVector.
function DynVector:init(tbl)
   self._numComponents = tbl.numComponents and tbl.numComponents or 1

   -- We store 1 extra element than requested to allow for 1-based
   -- indexing of returned values.
   local allocator = Alloc.createAllocator(
      string.format("double[%d]", self._numComponents+1))
   self._timeMesh = Alloc.Double()
   self._data = allocator()
   -- Temp storage for single entry.
   self._tmpData = new(string.format("double[%d]", self._numComponents+1))
   self._tmpTable = {}
   for i = 1, self._numComponents do self._tmpTable[i] = nil end

   -- Construct various functions from template representations.
   self._copyToTempData = loadstring( copyTempl {NCOMP=self._numComponents} )()

   -- Used to save the last data, in case we need it later.
   self.tLast, self.dLast = -1., Lin.Vec(self._numComponents)
   self.flushed = false

end

function DynVector:appendData(t, v)
   self._timeMesh:push(t)
   self._copyToTempData(v, self._tmpData)
   self._data:push(self._tmpData)
end

function DynVector:appendLast(dynV)
   local tm, lv = dynV:lastData()
   self:appendData(tm, lv)
end

function DynVector:assignLastTime(tm)
   self._timeMesh:assignLast(tm)
end

function DynVector:assignLastVal(fact, val)
   for i = 1, self._numComponents do self._tmpData[i] = val[i]*fact end
   self._data:assignLast(self._tmpData)
end

function DynVector:assignLast(fact, dynV)
   local _, lv = dynV:lastData()
   for i = 1, self._numComponents do self._tmpData[i] = lv[i]*fact end
   self._data:assignLast(self._tmpData)
end

function DynVector:accumulateLastOne(fact, dynV)
   local _, lv = dynV:lastData()
   local selflv
   if self:size() == 0 and self.flushed then
      selflv = self.dLast
   else
      selflv = self._data:last()
   end
   for i = 1, self._numComponents do self._tmpData[i] = selflv[i]+lv[i]*fact end
   self._data:assignLast(self._tmpData)
end

function DynVector:accumulateLast(c1, dynV1, ...)
   local args = {...}                  -- Package up rest of args as table.
   self:accumulateLastOne(c1, dynV1)   -- Accumulate first factor*vec:last() pair.
   for i = 1, #args/2 do               -- Accumulate rest of the factor*vec:last() pairs.
      self:accumulateLastOne(args[2*i-1], args[2*i])
   end
end

function DynVector:combineLast(c1, dynV1, ...)
   local args = {...}           -- Package up rest of args as table.
   self:assignLast(c1, dynV1)   -- Assign first factor*vec:last() pair.
   for i = 1, #args/2 do        -- Accumulate rest of the factor*vec:last() pairs.
      self:accumulateLastOne(args[2*i-1], args[2*i])
   end
end

function DynVector:_copy_from_dynvector(buff)
   local c = 1
   for i = 1, self._data:size() do
      local v = self._data[i]
      for n = 1, self._numComponents do
	 buff[c] = v[n]
	 c = c+1
      end
   end
end

function DynVector:_copy_to_dynvector(buff)
   local c = 1
   for i = 1, self._data:size() do
      local v = self._data[i]
      for n = 1, self._numComponents do
	 v[n] = buff[c]
	 c = c+1
      end
   end
end

function DynVector:clear()
   self._data:clear()
   self._timeMesh:clear()
end

function DynVector:copy(dynV)
   -- Copy over time-mesh and data. Assumes number of components is the same.
   local dataIn = dynV:data()
   local timeMeshIn = dynV:timeMesh()
   local nValIn = dataIn:size()

   self._timeMesh:expand(nValIn-(self:size()-1))
   for i = 1, self._timeMesh:size() do
      self._timeMesh[i] = timeMeshIn[i]
   end

   self._data:expand(nValIn-(self:size()-1))
   for i = 1, self._data:size() do
      for n = 1, self._numComponents do
	 self._data[i][n] = dataIn[i][n]
      end
   end
end

function DynVector:copyLast(dynV)
   local tm, lv = dynV:lastData()
   self._timeMesh:assignLast(tm)
   self._data:assignLast(lv)
end

function DynVector:data() return self._data end
function DynVector:lastTime() return self._timeMesh:last() end

function DynVector:lastTime()
   if self:size() == 0 and self.flushed then 
      return self.tLast
   else
      return self._timeMesh:last()
   end
end

function DynVector:lastData()
   if self:size() == 0 and self.flushed then 
      return self.tLast, self.dLast
   else
      return self._timeMesh:last(), self._data:last() 
   end
end

function DynVector:numComponents() return self._numComponents end

-- Returns time-stamp and frame number.
function DynVector:read(fName, flushData)
   assert(self.adiosIo, "DataStruct.DynVector: 'GKYL_ADIOS2_MPI' must be defined in order to do I/O.")
   self.adiosIo:read(self, fName)

   local tLast, dLast = self:lastData()

   if flushData then 
      self.tLast = tLast
      for i = 1, self._numComponents do self.dLast[i] = dLast[i] end
      self:clear()
      self.flushed = true
   end
   return tLast, dLast
end

function DynVector:removeLast()
   local tm = self._timeMesh:popLast()
   local v = self._data:popLast()
   return tm, v
end

function DynVector:size() return self._data:size() end

local function split (inputstr, sep)
   if sep == nil then
      sep = "%s"
   end
   local t={}
   for str in string.gmatch(inputstr, "([^"..sep.."]+)") do
      table.insert(t, str)
   end
   return t
end

function DynVector:timeMesh() return self._timeMesh end

function DynVector:write(outNm, tmStamp, frNum, flushData, appendData)
   local flushData = xsys.pickBool(flushData, true)  -- Default flush data on write.

   if appendData and (frNum and frNum>=0) then 
      frNum = frNum 
   else 
      frNum = frNum or 0
   end

   -- IO should go here
   io.write("Write method NYI!\n")

   -- Clear data for next round of IO.
   if flushData then 
      local tLast, dLast = self:lastData() -- Save the last data, in case we need it later.
      self.tLast = tLast
      for i = 1, self._numComponents do self.dLast[i] = dLast[i] end
      self:clear()
      self.flushed = true
   end
end

return { DynVector = DynVector }
