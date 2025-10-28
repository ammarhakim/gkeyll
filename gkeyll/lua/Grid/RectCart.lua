-- Gkyl ------------------------------------------------------------------------
--
-- Uniform Cartesian grids.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- System libraries.
local ffi = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- Gkyl libraries.
local DecompRegionCalc = require "Lib.CartDecomp"
local Lin = require "Lib.Linalg"
local lume = require "Lib.lume"
local Mpi = require "Comm.Mpi"
local Proto = require "Lib.Proto"
local Range = require "Lib.Range"

local ZeroRectCart = require "Grid.ZeroRectCart"

local rectCartSz = sizeof("struct gkyl_rect_grid")

-- Determine local domain index.
local function getSubDomIndex(comm)
   local idx = -1
   if Mpi.Is_comm_valid(comm) then
      idx = Mpi.Comm_rank(comm)+1 -- Sub-domains are 1-indexed.
   end
   return idx
end

-- RectCart --------------------------------------------------------------------
--
-- A uniform Cartesian grid
--------------------------------------------------------------------------------

-- Uniform Cartesian grid
local RectCart = Proto()

function RectCart:init(tbl)
   local cells = tbl.cells
   self._ndim = #cells
   local lo = tbl.lower and tbl.lower or {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
   local up = tbl.upper and tbl.upper or {1.0, 1.0, 1.0, 1.0, 1.0, 1.0}

   self._periodicDirs = tbl.periodicDirs and tbl.periodicDirs or {}
   self._isDirPeriodic = {false, false, false, false, false, false}
   for _, dir in ipairs(self._periodicDirs) do
      self._isDirPeriodic[dir] = true
   end

   self._lower = Lin.Vec(self._ndim)
   self._upper = Lin.Vec(self._ndim)
   self._dx = Lin.Vec(self._ndim)
   self._numCells = Lin.IntVec(self._ndim)
   self._currIdx = Lin.IntVec(self._ndim)

   self._vol = 1.0
   self._gridVol = 1.0
   for d = 1, #cells do
      self._lower[d], self._upper[d] = lo[d], up[d]
      self._numCells[d] = cells[d]
      self._dx[d] = (up[d]-lo[d])/cells[d]
      self._vol = self._vol*self._dx[d]
      self._gridVol = self._gridVol*(up[d]-lo[d])
   end

   -- create gkyl_rect_grid structure for use in g0
   self._zero = ffi.new("struct gkyl_rect_grid")
   ffiC.gkyl_rect_grid_init(self._zero, self._ndim, self._lower:data(), self._upper:data(), self._numCells:data())

   -- Compute global range.
   local l, u = {}, {}
   for d = 1, #cells do
      l[d], u[d] = tbl.rangeLower and tbl.rangeLower[d] or 1, tbl.rangeUpper and tbl.rangeUpper[d] or cells[d]
   end
   self._globalRange = Range.Range(l, u)   
   self._localRange = Range.Range(l, u)
   self._block = 1   -- Block number for use in parallel communications.
   
   self._messenger = tbl.messenger or nil  -- Object managing communications.
   self.decomp = tbl.decomposition and tbl.decomposition or nil  -- Decomposition.
   if self.decomp then
      assert(self.decomp:ndim() == self._ndim,
	     "Decomposition dimensions must be same as grid dimensions!")

      -- In parallel, we need to adjust local range. 
      self._commSet = self.decomp:commSet()
      self._decomposedRange = self.decomp:decompose(self._globalRange)
      local subDomIdx = getSubDomIndex(self._commSet.host)
      self._block = subDomIdx
      local localRange = self._decomposedRange:subDomain(subDomIdx)
      self._localRange:copy(localRange)
      self._cuts = {}
      for i = 1, self._ndim do 
         assert(self.decomp:cuts(i) <= self._numCells[i],
                "Cannot have more decomposition cuts than cells in any dimension!")
         self._cuts[i] = self.decomp:cuts(i) 
      end
   else
      -- Create a dummy decomp and use it to set the grid.
      local cuts = {}
      for i = 1, self._ndim do cuts[i] = 1 end
      local dec1 = DecompRegionCalc.CartProd { cuts = cuts }
      self._commSet = dec1:commSet()
      self._decomposedRange = dec1:decompose(self._globalRange)
      self._block = 1
      self._cuts = cuts
   end
end

-- Member functions.
function RectCart:id() return "uniform" end
function RectCart:commSet() return self._commSet end 
function RectCart:subGridId() return self._block end
function RectCart:subGridIdByDim(idx) 
   local invIdxr = self._decomposedRange:cutsInvIndexer()
   invIdxr(self._block, idx)
end
function RectCart:decomposedRange() return self._decomposedRange end
function RectCart:ndim() return self._ndim end
function RectCart:lower(dir) return self._lower[dir] end
function RectCart:mid(dir) return (self._lower[dir]+self._upper[dir])/2 end
function RectCart:upper(dir) return self._upper[dir] end
function RectCart:numCells(dir) return self._numCells[dir] end
function RectCart:totalNumCells()
   total = 1
   for d = 1, self._ndim do
      total = total*self._numCells[d] 
   end
   return total
end
function RectCart:localRange() return self._localRange end
function RectCart:globalRange() return self._globalRange end
function RectCart:isDirPeriodic(dir) return self._isDirPeriodic[dir] end
function RectCart:getPeriodicDirs() return self._periodicDirs end
function RectCart:setPeriodicDirs(periodicDirs)
   self._periodicDirs = periodicDirs
   self._isDirPeriodic = {false, false, false, false, false, false}
   for _, dir in ipairs(self._periodicDirs) do
      self._isDirPeriodic[dir] = true
   end
end
function RectCart:cuts(dir) return self._cuts[dir] end
function RectCart:setIndex(idxIn)
   if type(idxIn) == "cdata" then
      idxIn:copyInto(self._currIdx)
   else 
      for d = 1, self._ndim do self._currIdx[d] = idxIn[d] end
   end
end
function RectCart:dx(dir) return self._dx[dir] end
function RectCart:getDx(dxOut) 
   if type(dxOut) == "cdata" then
      self._dx:copyInto(dxOut)
   else 
      for d = 1, self._ndim do dxOut[d] = self._dx[d] end
   end
end
function RectCart:cellCenterInDir(d)
   return self:lower(d) + (self._currIdx[d]-0.5)*self:dx(d)
end
function RectCart:cellLowerInDir(d)
   return self:lower(d) + (self._currIdx[d]-1)*self:dx(d)
end
function RectCart:cellUpperInDir(d)
   return self:lower(d) + (self._currIdx[d])*self:dx(d)
end
function RectCart:cellCenter(xc)
   for d = 1, self._ndim do
      xc[d] = self:cellCenterInDir(d)
   end
end
function RectCart:cellVolume() return self._vol end
function RectCart:gridVolume() return self._gridVol end
function RectCart:findCell(point, cellIdx, pickLower, knownIdx)
   -- Return multidimensional index (cellIdx) of cell containing a point (n-dimensional coordinate).
   -- When multiple cells share that point, chooseLower=true returns the lower cell.
   -- Can specify a subset of known cell coordinates in knownIdx (e.g. knownIdx={nil,7} if you know cellIdx[2]=7).
   local pickLower = pickLower==nil and true or pickLower
   local knownIdx = knownIdx==nil and {} or knownIdx
   if #knownIdx == 0 then for d=1,self._ndim do knownIdx[d]=nil end end

   local searchDim, searchNum, dimTrans = {}, 0, {}
   for d=1,self._ndim do
      if knownIdx[d]==nil then
         searchNum = searchNum+1
         searchDim[searchNum] = d
         dimTrans[d] = searchNum
      else
         dimTrans[d] = nil
      end
   end

   local lessEq = function(a,b) return a<=b end
   local compPoints = function(pA, pB, compFunc)
      local isTrue = true
      for d = 1, searchNum do isTrue = isTrue and compFunc(pA[d],pB[d]) end
      return isTrue
   end
   local isInCell = function(pIn, iIn)
      local eps = 1.e-14
      local checkIdx = {}
      for d=1, self._ndim do checkIdx[d] = knownIdx[d]==nil and iIn[dimTrans[d]] or knownIdx[d] end
      self:setIndex(checkIdx)
      local inCell = true
      for d = 1, self._ndim do
         inCell = inCell and (self:cellLowerInDir(d)-eps<=pIn[d] and pIn[d]<=self:cellUpperInDir(d)+eps)
      end
      return inCell
   end

   -- Below we use a binary search. That is, if the i-th coordinate of the point in
   -- question is below(above) the ith-coordinate of the current mid-point, we search
   -- the lower(upper) half along that direction in the next iteration.
   local iStart, iEnd, iMid, iNew = {}, {}, {}, {}
   for d = 1, searchNum do
      iStart[d], iEnd[d] = self._globalRange:lower(searchDim[d]), self._globalRange:upper(searchDim[d])
      iMid[d], iNew[d] = 0, 0
   end
   while compPoints(iStart, iEnd, lessEq) do
      for d = 1, searchNum do iMid[d] = math.floor( (iStart[d]+iEnd[d])/2 ) end -- Calculate middle.
      if isInCell(point, iMid) then
         -- Check if neighboring cells also contains this point.
         local hcCells = {iMid}
         for d = 1, searchNum do
            local prevDimCells = #hcCells
            for pC = 1, prevDimCells do
               for pm = -1, 1, 2 do
                  local iNew = lume.clone(hcCells[pC])
                  iNew[d] = iNew[d]+pm
                  if (iNew[d] >= self._globalRange:lower(searchDim[d])) and   -- Do not consider ghost cells.
                     (iNew[d] <= self._globalRange:upper(searchDim[d])) then
                     table.insert(hcCells, iNew)
                  end
               end
            end
         end
         local cellsFound = {}
         for cI = 1, #hcCells do
            if isInCell(point, hcCells[cI]) then table.insert(cellsFound, hcCells[cI]) end
         end
         local idxOut = cellsFound[1]
         if pickLower then   -- Choose the lower-left most cell.
            for cI = 2, #cellsFound do
               local isHigher = false
               for d = 1, searchNum do isHigher = isHigher or (idxOut[d] > cellsFound[cI][d]) end
               if isHigher then idxOut = cellsFound[cI] end
            end
         else                -- Choose the upper-right most cell.
            for cI = 2, #cellsFound do
               local isLower = false
               for d = 1, searchNum do isLower = isLower or (idxOut[d] < cellsFound[cI][d]) end
               if isLower then idxOut = cellsFound[cI] end
            end
         end
         for d=1, self._ndim do cellIdx[d] = knownIdx[d]==nil and idxOut[dimTrans[d]] or knownIdx[d] end
         break
      else
         for d = 1, searchNum do 
            if point[searchDim[d]] < self:cellLowerInDir(searchDim[d]) then
               iEnd[d] = iMid[d]-1
            elseif point[searchDim[d]] > self:cellUpperInDir(searchDim[d]) then
               iStart[d] = iMid[d]+1
            end
         end
      end
   end
end

function RectCart:getNodalCoords() return nil end
function RectCart:write(basis, rankInComm) end

function RectCart:getMappings(dir)
   return nil
end

local _numMetricElems = { 1, 3, 6 }
-- size of vector to store metric elements
function RectCart:numMetricElems()
   return _numMetricElems[self:ndim()]
end

function RectCart:calcMetric(xc, gOut)
   if self._ndim == 1 then
      gOut[1] = 1.0 -- g_11
   elseif self._ndim == 2 then
      -- g_11, g_12, g_22
      gOut[1], gOut[2], gOut[3] = 1.0, 0.0, 1.0
   elseif self._ndim == 3 then
      -- g_11, g_12, g_13, g_22, g_23, g_33
      gOut[1], gOut[2], gOut[3], gOut[4], gOut[5], gOut[6] = 1.0, 0.0, 0.0, 1.0, 0.0, 1.0
   else
      assert(false, "RectCart:calcMetric does not support more than 3 dimensions!")
   end
end

function RectCart:calcContraMetric(xc, gOut)
   if self._ndim == 1 then
      gOut[1] = 1.0 -- g_11
   elseif self._ndim == 2 then
      -- g^11, g^12, g^22
      gOut[1], gOut[2], gOut[3] = 1.0, 0.0, 1.0
   elseif self._ndim == 3 then
      -- g^11, g^12, g^13, g^22, g^23, g^33
      gOut[1], gOut[2], gOut[3], gOut[4], gOut[5], gOut[6] = 1.0, 0.0, 0.0, 1.0, 0.0, 1.0
   else
      assert(false, "RectCart:calcContraMetric does not support more than 3 dimensions!")
   end
end

function RectCart:calcJacobian(xc)
   return 1.0
end

function RectCart:childGrid(keepDims)
   -- Collect the ingredients needed for a child grid: a grid with a subset of the
   -- directions of the original (parent) grid.
   -- Cannot output a childGrid itself because that would require a recursive RectCart updater.
   local childDim = #keepDims
   assert(childDim>0 and childDim<=self._ndim, "RectCart:childGrid cannot have dimensions < 1 or greater than parent grid.")

   local childLower, childUpper, childCells, childPeriodicDirs = {}, {}, {}, {}
   local dI, pIdx = 0, 0
   for _, d in ipairs(lume.sort(keepDims)) do
      dI = dI+1
      childLower[dI] = self:lower(d)
      childUpper[dI] = self:upper(d)
      childCells[dI] = self:numCells(d)
      if self:isDirPeriodic(d) then
         pIdx = pIdx+1
         childPeriodicDirs[pIdx] = dI
      end
   end

   local childDecomp = nil
   if self.decomp then
      local childComm, childWriteRank, childCuts = self.decomp:childDecomp(keepDims)
      childDecomp = DecompRegionCalc.CartProd {
         comm = childComm,       cuts = childCuts,
         writeRank = childWriteRank,
         __serTesting = true,
      }
   end

   local childGridIngredients = {
      lower = childLower,  periodicDirs = childPeriodicDirs,
      upper = childUpper,  decomposition = childDecomp,
      cells = childCells,  messenger = self._messenger,
   }
   return childGridIngredients
end

return RectCart
