-- Gkyl ------------------------------------------------------------------------
--
-- Multi-component fields on cartesian grids
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- System libraries.
local ffi = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"
local new, sizeof, typeof, metatype = xsys.from(ffi,
     "new, sizeof, typeof, metatype")

-- Gkyl libraries.
local Alloc = require "Lib.Alloc"
local CartDecompNeigh = require "Lib.CartDecompNeigh"
local Grid = require "Grid.RectCart"
local Lin = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local Mpi = require "Comm.Mpi"
local Range = require "Lib.Range"
local ZeroArray = require "DataStruct.ZeroArray"
local ZeroArrayRio = require "DataStruct.ZeroArrayRio"

local lume = require "Lib.lume"
local mpack = require "Lib.MessagePack"

local RangeVec = Lin.new_vec_ct(ffi.typeof("struct gkyl_range"))


-- C interfaces
ffi.cdef [[
/**
 * Create ghost and skin sub-ranges given parent range. The skin and
 * ghost ranges are sub-ranges of the parent range and DO NOT include
 * corners.
 *
 * @param skin On output, skin range
 * @param ghost On outout, ghost range
 * @param dir Direction in which skin/ghost are computed
 * @param edge Edge on which skin/ghost are computed
 * @param parent Range for which skin/ghost are computed
 * @param nghost Number of ghost cells in 'dir' are nghost[dir]
 */
void gkyl_skin_ghost_ranges(struct gkyl_range *skin, struct gkyl_range *ghost,
  int dir, int edge, const struct gkyl_range *parent, const int *nghost);
]]

-- Local definitions.
local rowMajLayout, colMajLayout = Range.rowMajor, Range.colMajor -- data layout
local indexerMakerFuncs = {} -- list of functions that make indexers
indexerMakerFuncs[rowMajLayout] = Range.makeRowMajorIndexer
indexerMakerFuncs[colMajLayout] = Range.makeColMajorIndexer
-- Default layout.
local defaultLayout = rowMajLayout

local genIndexerMakerFuncs = {} -- List of functions that make generic indexers.
genIndexerMakerFuncs[rowMajLayout] = Range.makeRowMajorGenIndexer
genIndexerMakerFuncs[colMajLayout] = Range.makeColMajorGenIndexer

-- Helper function to check if two fields are compatible.
local function field_compatible(y, x)
   return y:localRange() == x:localRange() and y:numComponents() == x:numComponents()
end
local function field_compatible_ext(y, x)
   return y:localExtRange() == x:localExtRange() and y:numComponents() == x:numComponents()
end
-- Helper function to check if two fields have the same range.
-- Useful when manipulating two fields with different number of components, but same range.
local function field_check_range(y, x)
   return y:localRange() == x:localRange()
end

-- Turn numbers in a table into strings and concatenate them.
local tblToStr = function(tblIn)
   local strOut = ""
   for _, v in ipairs(tblIn) do 
      strOut = v<0 and (strOut .. math.abs(v) .. 1) or (strOut .. math.abs(v) .. 2) 
   end
   return strOut
end

-- Field accessor object: allows access to field values in cell.
local function new_field_comp_ct(elct)
   local field_comp_mf = {
      data = function(self)
	 return self._cdata
      end
   }
   local field_comp_mt = {
      __index = function(self, k)
	 if type(k) == "number" then
	    return self._cdata[k-1]
	 else
	    return field_comp_mf[k]
	 end
      end,
      __newindex = function(self, k, v)
	 self._cdata[k-1] = v
      end,
   }
   return metatype(typeof("struct { int numComponents; $* _cdata; }", elct), field_comp_mt)
end


-- A function to create constructors for Field objects
local function Field_meta_ctor(elct)
   -- ctor for component data
   local fcompct = new_field_comp_ct(elct)

   local elctSize = sizeof(elct)
   
   -- Meta-data for type
   local isNumberType = false   
   local elctCommType = nil
   if ffi.istype(new(elct), new("double")) then
      elctCommType = Mpi.DOUBLE
      isNumberType = true
   elseif ffi.istype(new(elct), new("float")) then
      elctCommType = Mpi.FLOAT
      isNumberType = true
   elseif ffi.istype(new(elct), new("int")) then
      elctCommType = Mpi.INT
      isNumberType = true
   elseif ffi.istype(new(elct), new("long")) then
      elctCommType = Mpi.LONG
      isNumberType = true
   else
      elctCommType = Mpi.BYTE -- by default, send stuff as byte array
   end

   -- Binary operation functions and reduce MPI types (used in reduce method).
   local reduceOpsMPI = {max = Mpi.MAX, min = Mpi.MIN, sum = Mpi.SUM}
   
   -- Make constructor for Field.
   local Field = {}
   function Field:new(tbl)
      local self = setmetatable({}, Field)

      -- Read data from input table.
      local grid = tbl.onGrid
      local nc = tbl.numComponents and tbl.numComponents or 1 -- Default numComponents=1.
      local ghost = tbl.ghost and tbl.ghost or {0, 0} -- No ghost cells by default.

      self._syncCorners = xsys.pickBool(tbl.syncCorners, false) -- Don't sync corners by default.
      self._syncPeriodicDirs = xsys.pickBool(tbl.syncPeriodicDirs, true) -- Sync periodic BCs by default.

      -- Local and global ranges.
      local globalRange = grid:globalRange()
      local localRange = grid:localRange()

      local sz = localRange:extend(ghost[1], ghost[2]):volume()*nc -- Amount of data in field.
      -- Setup object.
      self._grid = grid
      self._ndim = grid:ndim()
      self._lowerGhost, self._upperGhost = ghost[1], ghost[2]
      self._numComponents = nc
      self._size = sz

      -- Underlying ZeroArray data structure, which handles memory allocation
      self._zero = ZeroArray.Array(ZeroArray.double, self._numComponents, self._size/self._numComponents, 0)
      -- get pointer to data in ZeroArray
      self._data = self._zero:data()

      -- Create a device copy if needed.
      if xsys.pickBool(tbl.useDevice, GKYL_USE_GPU) then
         self.useDevice = 1
         self._zeroDevice = ZeroArray.Array(ZeroArray.double, self._numComponents, self._size/self._numComponents, 1)
         self._devAllocData = self._zeroDevice:data()

         self._zeroForOps = self._zeroDevice
      else
         self._zeroDevice = nil
         self._devAllocData = nil

         self._zeroForOps = self._zero
      end

      -- Extended ranges include ghost cells.
      self._globalExtRange = globalRange:extend(self._lowerGhost, self._upperGhost)
      self._localExtRange = localRange:extend(self._lowerGhost, self._upperGhost)

      -- re-initialize localRange and globalRange as sub-ranges of localExtRange and globalExtRange
      self._localRange = self._localExtRange:subRange(localRange:lowerAsVec(), localRange:upperAsVec())
      self._globalRange = self._globalExtRange:subRange(globalRange:lowerAsVec(), globalRange:upperAsVec())

      self._lowerGhostVec, self._upperGhostVec = Lin.IntVec(self._ndim), Lin.IntVec(self._ndim)
      for d=1, self._ndim do
         self._lowerGhostVec[d], self._upperGhostVec[d] = self._lowerGhost, self._upperGhost
      end
      local lowerEdge, upperEdge = 0, 1 -- Match gkyl_edge_loc in gkylzero/zero/gkyl_range.h.
      -- Create local skin/ghost ranges in each direction.
      self._localGhostRngLo, self._localSkinRngLo = RangeVec(self._ndim), RangeVec(self._ndim)
      self._localGhostRngUp, self._localSkinRngUp = RangeVec(self._ndim), RangeVec(self._ndim)
      for d=1, self._ndim do
         ffiC.gkyl_skin_ghost_ranges(self._localSkinRngLo[d], self._localGhostRngLo[d], d-1, lowerEdge, self._localExtRange, self._lowerGhostVec:data())
         ffiC.gkyl_skin_ghost_ranges(self._localSkinRngUp[d], self._localGhostRngUp[d], d-1, upperEdge, self._localExtRange, self._upperGhostVec:data())
      end

      -- Create global skin/ghost ranges in each direction.
      self._globalGhostRngLo, self._globalSkinRngLo = RangeVec(self._ndim), RangeVec(self._ndim)
      self._globalGhostRngUp, self._globalSkinRngUp = RangeVec(self._ndim), RangeVec(self._ndim)
      for d=1, self._ndim do
         ffiC.gkyl_skin_ghost_ranges(self._globalSkinRngLo[d], self._globalGhostRngLo[d], d-1, lowerEdge, self._globalExtRange, self._lowerGhostVec:data())
         ffiC.gkyl_skin_ghost_ranges(self._globalSkinRngUp[d], self._globalGhostRngUp[d], d-1, upperEdge, self._globalExtRange, self._upperGhostVec:data())
      end

      -- Create ranges covering the global ghost range owned by this MPI rank.
      -- Ranks without part of the global ghost range return a range with 0 volume.
      -- Use: when we want every rank to pass a range to an operation,
      -- but only want the boundary ranks to actually perform the operation.
      self._locGloGhostRngInterLo, self._locGloSkinRngInterLo = RangeVec(self._ndim), RangeVec(self._ndim)
      self._locGloGhostRngInterUp, self._locGloSkinRngInterUp = RangeVec(self._ndim), RangeVec(self._ndim)
      for d=1, self._ndim do
         local ghostLo = self._globalGhostRngLo[d]:intersect(self._localGhostRngLo[d])
         local ghostUp = self._globalGhostRngUp[d]:intersect(self._localGhostRngUp[d])
         local skinLo = self._globalSkinRngLo[d]:intersect(self._localSkinRngLo[d])
         local skinUp = self._globalSkinRngUp[d]:intersect(self._localSkinRngUp[d])
         self._locGloGhostRngInterLo[d] = self._localExtRange:subRange(ghostLo:lowerAsVec(),ghostLo:upperAsVec())
         self._locGloGhostRngInterUp[d] = self._localExtRange:subRange(ghostUp:lowerAsVec(),ghostUp:upperAsVec())
         self._locGloSkinRngInterLo[d] = self._localExtRange:subRange(skinLo:lowerAsVec(),skinLo:upperAsVec())
         self._locGloSkinRngInterUp[d] = self._localExtRange:subRange(skinUp:lowerAsVec(),skinUp:upperAsVec())
      end

      -- All real-cell edges.
      self._localEdgeRange = self._localRange:extend(1, 0) -- Or (1, 0)?

      -- All cell-cell edges, including those of a ghost cell.
      self._localExtEdgeRange = self._localRange:extend(self._lowerGhost-1, self._upperGhost)

      -- Local and (MPI) global values of a reduction (reduce method).
      if self.useDevice then
	 local on_gpu = 1
         self.localReductionVal = ZeroArray.Array(ZeroArray.double, self._numComponents, 1, on_gpu)
         self.globalReductionVal = ZeroArray.Array(ZeroArray.double, self._numComponents, 1)
         self.localReductionVal_h = ZeroArray.Array(ZeroArray.double, self._numComponents, 1)
      else
         self.localReductionVal = ZeroArray.Array(ZeroArray.double, self._numComponents, 1)
         self.globalReductionVal = ZeroArray.Array(ZeroArray.double, self._numComponents, 1)
	 self.localReductionVal_h = self.localReductionVal
      end
      
      self._layout = defaultLayout -- Default layout is column-major.
      if tbl.layout then
         self._layout = tbl.layout=="row-major" and rowMajLayout or colMajLayout
      end

      -- Store start index and size handled by local rank for local and extended range.
      self._localStartIdx, self._localNumBump = self._localRange:lowerAsVec(), self._localRange:volume()
      self._localExtStartIdx, self._localExtNumBump = self._localExtRange:lowerAsVec(), self._localExtRange:volume()

      -- Compute communication neighbors.
      self._decompNeigh = CartDecompNeigh(grid:decomposedRange())
      if self._syncCorners then
         self._decompNeigh:calcAllCommNeigh(ghost[1], ghost[2])
      else
         self._decompNeigh:calcFaceCommNeigh(ghost[1], ghost[2])
      end

      local localExtRange = self._localExtRange
      local indexer = self:genIndexer()
      local decomposedRange = self._grid:decomposedRange()
      local myId = self._grid:subGridId() -- Grid ID on this processor.
      local neigIds = self._decompNeigh:neighborData(myId) -- List of neighbors.

      -- Create buffers for periodic copy if Mpi.Comm_size(nodeComm) = 1.
      for dir = 1, self._ndim do
         if self._grid:cuts(dir) == 1 then
            self._lowerPeriodicBuff, self._upperPeriodicBuff = {}, {}
            break
         end
      end

      -- Following loop creates buffers for periodic directions.
      -- This is complicated as one needs to treat lower -> upper
      -- transfers differently than upper -> lower as the number of
      -- ghost cells may be different on each lower/upper side. (AHH)
      for dir = 1, self._ndim do
         -- Set up periodic-sync buffers for all dirs.
         if self._lowerGhost > 0 and self._upperGhost > 0 then
            local skelIds = decomposedRange:boundarySubDomainIds(dir)
            for i = 1, #skelIds do
               local loId, upId = skelIds[i].lower, skelIds[i].upper

               if myId == loId then  -- Redundant for cuts=1.
                  local rgnSend = decomposedRange:subDomain(loId):lowerSkin(dir, self._upperGhost)
                  if self._grid:cuts(dir) == 1 then
                     self._lowerPeriodicBuff[dir] = ZeroArray.Array(ZeroArray.double, self._numComponents, rgnSend:volume(), self.useDevice)
                  end
               end
               if myId == upId then  -- Redundant for cuts=1.
                  local rgnSend = decomposedRange:subDomain(upId):upperSkin(dir, self._lowerGhost)
                  if self._grid:cuts(dir) == 1 then
                     self._upperPeriodicBuff[dir] = ZeroArray.Array(ZeroArray.double, self._numComponents, rgnSend:volume(), self.useDevice)
                  end
               end	       
            end
         end
      end
      

      -- Get info for syncs without MPI datatypes. We'll loop over blocks and post a send/recv pair ourselves.
      self._syncSendRng, self._syncRecvRng = {}, {}
      self._syncSendBufVol, self._syncRecvBufVol = 0, 0
      for _, sendId in ipairs(neigIds) do
         local neighRgn = decomposedRange:subDomain(sendId)
         local sendRgn = self._localRange:intersect(
            neighRgn:extend(self._lowerGhost, self._upperGhost))
         self._syncSendRng[sendId] = self._localExtRange:subRange(sendRgn:lowerAsVec(), sendRgn:upperAsVec())
         self._syncSendBufVol = math.max(self._syncSendBufVol, sendRgn:volume())
      end
      for _, recvId in ipairs(neigIds) do
         local neighRgn = decomposedRange:subDomain(recvId)
         local recvRgn = self._localExtRange:intersect(neighRgn)
         self._syncRecvRng[recvId] = self._localExtRange:subRange(recvRgn:lowerAsVec(), recvRgn:upperAsVec())
         self._syncRecvBufVol = math.max(self._syncRecvBufVol, recvRgn:volume())
      end

      -- Get info for periodic syncs without MPI datatypes.
      self._syncPerNeigh = {}
      self._syncPerSendRng, self._syncPerRecvRng = {}, {}
      self._onBoundary = {}
      for dir = 1, self._ndim do
         -- set up periodic-sync Datatypes for all dirs, in case we want to change periodicDirs later
         if self._lowerGhost > 0 and self._upperGhost > 0 and decomposedRange:numSubDomains() > 1 then
            self._onBoundary[dir] = false
            local skelIds = decomposedRange:boundarySubDomainIds(dir)
            for i = 1, #skelIds do
               local loId, upId = skelIds[i].lower, skelIds[i].upper
	       self._onBoundary[dir] = self._onBoundary[dir] or ((myId==skelIds[i].lower) or (myId==skelIds[i].upper))
               -- Only create if we are on proper ranks.
               -- Note that if the node communicator has rank size of 1, then we can access all the
               -- memory needed for periodic boundary conditions and no communication is needed.
               local rgnSend, rgnRecv, oppId
               if myId == loId and self._grid:cuts(dir) > 1 then
                  rgnSend = decomposedRange:subDomain(loId):lowerSkin(dir, self._upperGhost)
                  rgnRecv = decomposedRange:subDomain(loId):lowerGhost(dir, self._lowerGhost)
                  oppId = upId
               end
               if myId == upId and self._grid:cuts(dir) > 1 then
                  rgnSend = decomposedRange:subDomain(upId):upperSkin(dir, self._lowerGhost)
                  rgnRecv = decomposedRange:subDomain(upId):upperGhost(dir, self._upperGhost)
                  oppId = loId
               end
               if oppId ~= nil then
                  self._syncPerNeigh[dir] = oppId
                  self._syncPerSendRng[dir] = self._localExtRange:subRange(rgnSend:lowerAsVec(), rgnSend:upperAsVec())

                  self._syncPerRecvRng[dir] = self._localExtRange:subRange(rgnRecv:lowerAsVec(), rgnRecv:upperAsVec())
               end
            end
         end
      end

      self._metaData = tbl.metaData
      mpack.set_number("double")
      -- tag to identify basis used to set this field
      self._basisId = "none"

      return self
   end
   setmetatable(Field, { __call = function (self, o) return self.new(self, o) end })

   -- Set callable methods.
   Field.__index = {
      data = function(self)
         return self._zeroForOps:data()
      end,
      elemType = function(self)
         return elct
      end,
      elemSize = function(self)
         return elctSize
      end,
      elemCommType = function(self)
         return elctCommType
      end,
      ndim = function(self)
         return self._ndim
      end,
      grid = function(self)
         return self._grid
      end,
      numComponents = function(self)
         return self._numComponents
      end,
      hasCuDev = function(self)
         if self._zeroDevice then return self._zeroDevice:is_cu_dev() else return false end
      end,
      setHostOps = function(self)
         self._zeroForOps = self._zero
      end,
      setDeviceOps = function(self)
         self._zeroForOps = self._zeroDevice
      end,
      copy = function(self, fIn)
         self._zeroForOps:copy(fIn._zeroForOps)
      end,
      deviceCopy = function(self, fIn)
         self._zeroDevice:copy(fIn._zeroDevice)
      end,
      copyHostToDevice = function(self)
         if self._zeroDevice then self._zeroDevice:copy(self._zero) end
      end,
      copyDeviceToHost = function(self)
         if self._zeroDevice then self._zero:copy(self._zeroDevice) end
      end,
      copyDeviceToHostAsync = function(self)
         self._zero:copyAsync(self._zeroDevice)
      end,
      copyHostToDeviceAsync = function(self)
         self._zero:copyAsync(self._zeroDevice)
      end,
      copyRangeToRange = function(self, fIn, outRange, inRange)
         self._zeroForOps:copyRangeToRange(fIn._zeroForOps, outRange, inRange)
      end,
      copyRange = function(self, fIn, inRange)
         self._zeroForOps:copyRange(fIn._zeroForOps, inRange)
      end,
      deviceDataPointer = function(self)
         return self._devAllocData
      end,
      dataPointer = function(self)
         return self._data
      end,
      clear = function(self, val)
         self._zeroForOps:clear(val)
      end,
      clearRange = function (self, val, rng)
	 self._zeroForOps:clearRange(val, rng)
      end,
      fill = function(self, k, fc)
         local loc = (k - 1) * self._numComponents -- (k-1) as k is 1-based index	
         fc._cdata = self._data + loc
      end,
      _localLower = function(self)
         return 0
      end,
      _localShape = function(self)
         return self._localExtRange:volume() * self:numComponents()
      end,
      _assign = function(self, fact, fld)
         self._zeroForOps:set(fact, fld._zeroForOps)
      end,
      _assignRange = function(self, fact, fld, rng)
         self._zeroForOps:setRange(fact, fld._zeroForOps, rng)
      end,
      -- assignOffsetOneFld assumes that one of the input or output fields have fewer components than the other.
      --   a) nCompOut > nCompIn: assigns all of the input field w/ part of the output field, the (0-based)
      --                          offset indicating which is the 1st component in the output field to assign.
      --   b) nCompOut < nCompIn: assigns nCompOut components of the input field onto the output field, the
      --                          (0-based) offset indicates which is the first component in the input field to read.
      -- It assumes that the components being summed are continuous.
      _assignOffsetOneFld = function(self, fact, fld, compStart)
	 assert(field_check_range(self, fld),
		"CartField:assignOffsetOneFld: Can only assign fields with the same range")
	 assert(type(fact) == "number",
		"CartField:assignOffsetOneFld: Factor not a number")
         assert(self:layout() == fld:layout(),
		"CartField:assignOffsetOneFld: Fields should have same layout for sums to make sense")
         self._zeroForOps:setOffset(fact, fld._zeroForOps, compStart)
      end,
      _assignOffsetOneFldRange = function(self, fact, fld, compStart, rng)
	 assert(field_check_range(self, fld),
		"CartField:assignOffsetOneFld: Can only assign fields with the same range")
	 assert(type(fact) == "number",
		"CartField:assignOffsetOneFld: Factor not a number")
         assert(self:layout() == fld:layout(),
		"CartField:assignOffsetOneFld: Fields should have same layout for sums to make sense")
         self._zeroForOps:setOffsetRange(fact, fld._zeroForOps, compStart, rng)
      end,
      _accumulateOneFld = function(self, fact, fld)
         self._zeroForOps:accumulate(fact, fld._zeroForOps)
      end,
      _accumulateOneFldRange = function(self, fact, fld, rng)
         self._zeroForOps:accumulateRange(fact, fld._zeroForOps, rng)
      end,
      -- accumulateOffsetOneFld assumes that one of the input or output fields have fewer components than the other.
      --   a) nCompOut > nCompIn: accumulates all of the input field w/ part of the output field, the (0-based)
      --                          offset indicating which is the 1st component in the output field to accumulate to.
      --   b) nCompOut < nCompIn: accumulates nCompIn components of the input field onto the output field, the
      --                          (0-based) offset indicates which is the first component in the input field to accumulate. 
      -- It assumes that the components being summed are continuous.
      _accumulateOffsetOneFld = function(self, fact, fld, compStart)
         assert(field_check_range(self, fld),
            "CartField:accumulateOffsetOneFld: Can only accumulate fields with the same range")
         assert(type(fact) == "number",
            "CartField:accumulateOffsetOneFld: Factor not a number")
         assert(self:layout() == fld:layout(),
            "CartField:accumulateOffsetOneFld: Fields should have same layout for sums to make sense")
         self._zeroForOps:accumulateOffset(fact, fld._zeroForOps, compStart)
      end,
      _accumulateOffsetOneFldRange = function(self, fact, fld, compStart, rng)
         assert(field_check_range(self, fld),
            "CartField:accumulateOffsetOneFld: Can only accumulate fields with the same range")
         assert(type(fact) == "number",
            "CartField:accumulateOffsetOneFld: Factor not a number")
         assert(self:layout() == fld:layout(),
            "CartField:accumulateOffsetOneFld: Fields should have same layout for sums to make sense")
         self._zeroForOps:accumulateOffsetRange(fact, fld._zeroForOps, compStart, rng)
      end,
      accumulate = isNumberType and
	 function (self, c1, fld1, ...)
	    local args = {...} -- Package up rest of args as table.
	    local nFlds = #args/2
	    self:_accumulateOneFld(c1, fld1) -- Accumulate first field.
	    for i = 1, nFlds do -- Accumulate rest of the fields.
	       self:_accumulateOneFld(args[2*i-1], args[2*i])
	    end
          end or
          function(self, c1, fld1, ...)
             assert(false, "CartField:accumulate: Accumulate only works on numeric fields")
          end,
      accumulateOffset = isNumberType and
          function(self, c1, fld1, compStart1, ...)
             local args = { ... }                       -- Package up rest of args as table.
             local nFlds = #args / 3
             self:_accumulateOffsetOneFld(c1, fld1, compStart1) -- Accumulate first field.
             for i = 1, nFlds do                        -- Accumulate rest of the fields
                self:_accumulateOffsetOneFld(args[3 * i - 2], args[3 * i - 1], args[3 * i])
             end
          end or
          function(self, c1, fld1, compStart1, ...)
             assert(false, "CartField:accumulateOffset: Accumulate only works on numeric fields")
          end,
      combine = isNumberType and
          function(self, c1, fld1, ...)
             local args = { ... }  -- Package up rest of args as table.
             local nFlds = #args / 2
             self:_assign(c1, fld1) -- Assign first field.
             for i = 1, nFlds do   -- Accumulate rest of the fields.
                self:_accumulateOneFld(args[2 * i - 1], args[2 * i])
             end
          end or
          function(self, c1, fld1, ...)
             assert(false, "CartField:combine: Combine only works on numeric fields")
          end,
      combineOffset = isNumberType and
          function(self, c1, fld1, compStart1, ...)
             local args = { ... } -- Package up rest of args as table.
             local nFlds = #args / 3
             local notAssigned = {}
             for i = 1, self:numComponents() do table.insert(notAssigned, true) end -- Boolean indicates if already assigned.
             self:_assignOffsetOneFld(c1, fld1, compStart1)                       -- Assign first field.
             notAssigned[compStart1 + 1] = false
             for i = 1, nFlds do                                                  -- Accumulate rest of the fields.
                local cOff = args[3 * i]
                if notAssigned[cOff + 1] then
                   self:_assignOffsetOneFld(args[3 * i - 2], args[3 * i - 1], cOff)
                   notAssigned[cOff + 1] = false
                else
                   self:_accumulateOffsetOneFld(args[3 * i - 2], args[3 * i - 1], cOff)
                end
             end
          end or
          function(self, c1, fld1, ...)
             assert(false, "CartField:combineOffset: Combine only works on numeric fields")
          end,
      accumulateRange = isNumberType and
          function(self, c1, fld1, ...)
             local args = { ... } -- Package up rest of args as table.
             local nFlds = (#args - 1) / 2
             local rng = args[#args]
             self:_accumulateOneFldRange(c1, fld1, rng) -- Accumulate first field.
             for i = 1, nFlds do                -- Accumulate rest of the fields.
                self:_accumulateOneFldRange(args[2 * i - 1], args[2 * i], rng)
             end
          end or
          function(self, c1, fld1, ...)
             assert(false, "CartField:accumulate: Accumulate only works on numeric fields")
          end,
      accumulateOffsetRange = isNumberType and
          function(self, c1, fld1, compStart1, ...)
             local args = { ... } -- Package up rest of args as table.
             local nFlds = #args / 3
             local rng = args[#args]
             self:_accumulateOffsetOneFldRange(c1, fld1, compStart1, rng) -- Accumulate first field.
             for i = 1, nFlds do                                  -- Accumulate rest of the fields
                self:_accumulateOffsetOneFldRange(args[3 * i - 2], args[3 * i - 1], args[3 * i], rng)
             end
          end or
          function(self, c1, fld1, compStart1, ...)
             assert(false, "CartField:accumulateOffset: Accumulate only works on numeric fields")
          end,
      combineRange = isNumberType and
          function(self, c1, fld1, ...)
             local args = { ... } -- Package up rest of args as table.
             local nFlds = #args / 2
             local rng = args[#args]
             self:_assignRange(c1, fld1, rng) -- Assign first field.
             for i = 1, nFlds do             -- Accumulate rest of the fields.
                self:_accumulateOneFldRange(args[2 * i - 1], args[2 * i], rng)
             end
          end or
          function(self, c1, fld1, ...)
             assert(false, "CartField:combine: Combine only works on numeric fields")
          end,
      combineOffsetRange = isNumberType and
          function(self, c1, fld1, compStart1, ...)
             local args = { ... } -- Package up rest of args as table.
             local nFlds = #args / 3
             local rng = args[#args]
             local notAssigned = {}
             for i = 1, self:numComponents() do table.insert(notAssigned, true) end -- Boolean indicates if already assigned.
             self:_assignOffsetOneFldRange(c1, fld1, compStart1, rng)             -- Assign first field.
             notAssigned[compStart1 + 1] = false
             for i = 1, nFlds do                                                  -- Accumulate rest of the fields.
                local cOff = args[3 * i]
                if notAssigned[cOff + 1] then
                   self:_assignOffsetOneFldRange(args[3 * i - 2], args[3 * i - 1], cOff, rng)
                   notAssigned[cOff + 1] = false
                else
                   self:_accumulateOffsetOneFldRange(args[3 * i - 2], args[3 * i - 1], cOff, rng)
                end
             end
          end or
          function(self, c1, fld1, ...)
             assert(false, "CartField:combineOffset: Combine only works on numeric fields")
          end,
      scale = isNumberType and
          function(self, fact)
             self._zeroForOps:scale(fact)
          end or
          function(self, fact)
             assert(false, "CartField:scale: Scale only works on numeric fields")
          end,
      scaleByCell = function(self, factByCell)
         self._zeroForOps:scale_by_cell(factByCell._zeroForOps)
      end,
      shiftc = function(self, val, comp)
         self._zeroForOps:shiftc(val, comp)
      end,
      shiftcRange = function(self, val, comp, rng)
         self._zeroForOps:shiftc(val, comp, rng)
      end,
      defaultLayout = function(self)
         if defaultLayout == rowMajLayout then
            return "row-major"
         end
         return "col-major"
      end,
      layout = function(self)
         if self._layout == rowMajLayout then
            return "row-major"
         end
         return "col-major"
      end,
      lowerGhost = function(self)
         return self._lowerGhost
      end,
      upperGhost = function(self)
         return self._upperGhost
      end,
      lowerGhostVec = function(self)
         return self._lowerGhostVec
      end,
      upperGhostVec = function(self)
         return self._upperGhostVec
      end,
      localRange = function(self)
         return self._localRange
      end,
      localExtRange = function(self)  -- includes ghost cells
         return self._localExtRange
      end,
      localEdgeRange = function(self)
         return self._localEdgeRange
      end,
      localExtEdgeRange = function(self)
         return self._localExtEdgeRange
      end,
      globalRange = function(self)
         return self._globalRange
      end,
      globalExtRange = function(self)  -- includes ghost cells
         return self._globalExtRange
      end,
      localGhostRangeLower = function(self)
         return self._localGhostRngLo
      end,
      localGhostRangeUpper = function(self)
         return self._localGhostRngUp
      end,
      globalGhostRangeLower = function(self)
         return self._globalGhostRngLo
      end,
      globalGhostRangeUpper = function(self)
         return self._globalGhostRngUp
      end,
      localGlobalGhostRangeIntersectLower = function(self)
         return self._locGloGhostRngInterLo
      end,
      localGlobalGhostRangeIntersectUpper = function(self)
         return self._locGloGhostRngInterUp
      end,
      localSkinRangeLower = function(self)
         return self._localSkinRngLo
      end,
      localSkinRangeUpper = function(self)
         return self._localSkinRngUp
      end,
      globalSkinRangeLower = function(self)
         return self._globalSkinRngLo
      end,
      globalSkinRangeUpper = function(self)
         return self._globalSkinRngUp
      end,
      localGlobalSkinRangeIntersectLower = function(self)
         return self._locGloSkinRngInterLo
      end,
      localGlobalSkinRangeIntersectUpper = function(self)
         return self._locGloSkinRngInterUp
      end,
      localRangeIter = function(self)
         if self._layout == rowMajLayout then
            return self._localRange:rowMajorIter(self._localStartIdx, self._localNumBump)
         end
         return self._localRange:colMajorIter(self._localStartIdx, self._localNumBump)
      end,
      localExtRangeIter = function(self)  -- includes ghost cells
         local lext = self:localRange():extend(self:lowerGhost(), self:upperGhost())
         if self._layout == rowMajLayout then
            return lext:rowMajorIter(self._localExtStartIdx, self._localExtNumBump)
         end
         return lext:colMajorIter(self._localExtStartIdx, self._localExtNumBump)
      end,
      size = function(self)
         return self._size
      end,
      indexer = function(self)  -- linear indexer taking (i,j,...)
         return indexerMakerFuncs[self._layout](self:localExtRange())
      end,
      genIndexer = function(self)  -- linear indexer taking indices as a vector
         return genIndexerMakerFuncs[self._layout](self:localExtRange())
      end,
      get = function(self, k)           -- k is an integer returned by a linear indexer
         local loc = (k - 1) * self._numComponents -- (k-1) as k is 1-based index
         return fcompct(self._numComponents, self._data + loc)
      end,
      getDataPtrAt = function(self, k)  -- k is an integer returned by a linear indexer
         local loc = (k - 1) * self._numComponents -- (k-1) as k is 1-based index
         return self._data + loc
      end,
      write = function(self, fName)
	 local meta = ffi.new("struct gkyl_array_meta")

	 if self._metaData then
	    local metaDataPacked = mpack.pack(self._metaData)
	    meta.meta_sz = string.len(metaDataPacked)
	    meta.meta = ffi.new("char[?]", string.len(metaDataPacked)+1)
	    ffi.copy(meta.meta, metaDataPacked)
	 else
	    meta.meta_sz = 0
	 end

	 local fullNm = GKYL_OUT_PREFIX .. "_" .. fName -- Concatenate prefix.
         ffi.C.gkyl_grid_sub_array_write(self._grid._zero, self._localExtRange, meta, self._zero, fullNm)
      end,
      read = function(self, fName)  --> time-stamp, frame-number
	 local meta = ffi.new("struct gkyl_array_header_info")
	 ffi.C.gkyl_grid_sub_array_header_read(self._grid._zero, meta, fName)
	 self._metaData = 0
	 if (meta.meta_size > 0) then
	    self._metaData = mpack.unpack(ffi.string(meta.meta, meta.meta_size))
	    ffi.C.gkyl_array_header_info_release(meta)
	 end

	 ffi.C.gkyl_grid_sub_array_read(self._grid._zero, self._localExtRange, self._zero, fName)
         return 0.0, 0
      end,
      sync = function(self, syncPeriodicDirs_)
         local syncPeriodicDirs = xsys.pickBool(syncPeriodicDirs_, true)
         local mess = self._grid:getMessenger()

         -- MF 2023/03/29: quick fix for runs without decomposition/messenger.
         -- Would be better to create indirection pointing to functions defined
         -- at creation so this if-statement is not used every time.
         if mess then
            mess:syncCartField(self, mess:getConfComm())
            if self._syncPeriodicDirs and syncPeriodicDirs then
               mess:syncPeriodicCartField(self, mess:getConfComm())
            end
         else
            if self._syncPeriodicDirs and syncPeriodicDirs then
               self._field_periodic_sync(self, self._zeroForOps:data())
            end
         end
      end,
      -- This method is an alternative function for applying periodic boundary
      -- conditions when Mpi.Comm_size(nodeComm) = 1 and we do not need to call
      -- Send/Recv to copy the skin cell data into ghost cells.
      periodicCopy = function(self, syncPeriodicDirs_)
         local syncPeriodicDirs = xsys.pickBool(syncPeriodicDirs_, true)
         if self._syncPeriodicDirs and syncPeriodicDirs then
            self._field_periodic_copy(self)
         end
      end,
      periodicCopyInDir = function(self, dir)
         self._field_periodic_copy_indir(self, dir)
      end,
      setBasisId = function(self, basisId)
         self._basisId = basisId
      end,
      getBasisId = function(self)
         return self._basisId
      end,
      getMetaData = function(self)
         return self._metaData
      end,
      compatible = function(self, fld)
         return field_compatible(self, fld)
      end,
      checkRange = function(self, fld)
         return field_check_range(self, fld)
      end,
      reduce = isNumberType and
          function(self, opIn)
             self._zeroForOps:reduceRange(self.localReductionVal:data(), opIn, self._localRange)
             if self.useDevice then self.localReductionVal_h:copy(self.localReductionVal) end

             -- Input 'opIn' must be one of the binary operations in binOpFuncs.
             local grid = self._grid
             local localVal = {}
             Mpi.Allreduce(self.localReductionVal_h:data(), self.globalReductionVal:data(),
                self._numComponents, elctCommType, reduceOpsMPI[opIn], grid:commSet().host)

             --self.localReductionVal:copy(self.globalReductionVal)
             for k = 1, self._numComponents do localVal[k] = self.globalReductionVal:data()[k - 1] end
             return localVal
             --return self.localReductionVal:data()
          end or
          function(self, opIn)
             assert(false, "CartField:reduce: Reduce only works on numeric fields")
          end,
      copyRangeToBuffer = function(self, rgn, dataPointer)
         -- Copy the data in a range of this CartField to a buffer.
         self._zeroForOps:copy_to_buffer(dataPointer, rgn)
      end,
      copyRangeFromBuffer = function(self, rgn, dataPointer)
         -- Copy the data in a range of this CartField to a buffer.
         self._zeroForOps:copy_from_buffer(dataPointer, rgn)
      end,
      _field_periodic_sync = function(self, dataPtr)
         local comm = self._grid:commSet().host -- Communicator to use.
         if not Mpi.Is_comm_valid(comm) then
            return                                  -- No need to do anything if communicator is not valid
         end

         -- Immediately return if nothing to sync.
         if self._lowerGhost == 0 and self._upperGhost == 0 then return end

         local decomposedRange = self._grid:decomposedRange()
         assert((decomposedRange:numSubDomains() == 1) and (not self._syncCorners), "CartField: sync without messenger only allowed for serial runs.")
         return self:periodicCopy()
      end,
      _field_periodic_copy_indir = function(self, dir)
         -- First get region for skin cells for upper ghost region (the lower skin cells).
         local skinRgnLower = self._localSkinRngLo[dir]
         local skinRgnUpper = self._localSkinRngUp[dir]
         local ghostRgnLower = self._localGhostRngLo[dir]
         local ghostRgnUpper = self._localGhostRngUp[dir]

         local periodicBuffUpper = self._upperPeriodicBuff[dir]
         -- Copy skin cells into temporary buffer.
         self:copyRangeToBuffer(skinRgnUpper, periodicBuffUpper:data())
         -- Get region for looping over upper ghost cells and copy lower skin cells into upper ghost cells.
         self:copyRangeFromBuffer(ghostRgnLower, periodicBuffUpper:data())

         -- Now do the same, but for the skin cells for the lower ghost region (the upper skin cells).
         local periodicBuffLower = self._lowerPeriodicBuff[dir]
         self:copyRangeToBuffer(skinRgnLower, periodicBuffLower:data())
         self:copyRangeFromBuffer(ghostRgnUpper, periodicBuffLower:data())
      end,
      _field_periodic_copy = function(self)
         local grid = self._grid
         for dir = 1, self._ndim do
            if grid:isDirPeriodic(dir) then self:_field_periodic_copy_indir(dir) end
         end
      end,
   }

   return Field
end

return {
   new_field_ct = Field_meta_ctor,
   Field = Field_meta_ctor(typeof("double")),
}
