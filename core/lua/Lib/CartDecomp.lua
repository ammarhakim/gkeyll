-- Gkyl ------------------------------------------------------------------------
--
-- Decomposition algorithms for rectangular grids.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- System libraries.
local ffi = require "ffi"

-- Gkyl libraries.
local Lin = require "Lib.Linalg"
local Mpi = require "Comm.Mpi"
local PrimeFactor = require "Lib.PrimeFactor"
local Proto = require "Lib.Proto"
local Range = require "Lib.Range"
local xsys = require "xsys"
local lume = require "Lib.lume"

-- Create constructor to store vector of Range objects.
local RangeVec = Lin.new_vec_ct(ffi.typeof("struct gkyl_range"))
-- Create constructor to store integer pairs.
local PairsVec = Lin.new_vec_ct(ffi.typeof("struct { int lower, upper; } "))

local wrapToRange = function(xIn,xMin,xMax)
   return (((xIn - xMin) % (xMax+1 - xMin)) + (xMax+1 - xMin)) % (xMax+1 - xMin) + xMin
end

-- DecomposedRange --------------------------------------------------------------
--
-- The decomposition object's decompose() method creates a
-- DecomposedRange object.
--------------------------------------------------------------------------------

local DecomposedRange = Proto()
-- Constructor to make object that stores decomposed ranges.
function DecomposedRange:init(decomp)
   local ones, upper = {}, {}
   for d = 1, decomp:ndim() do
      ones[d] = 1; upper[d] = decomp:cuts(d)
   end

   self._cutsRange = Range.Range(ones, upper)
   self._commSet = decomp:commSet() -- Set of communicators.
   self._domains = RangeVec(self._cutsRange:volume())

   -- For periodic directions store "skeleton" sub-domains, i.e. those
   -- that touch the domain boundary.
   self._periodicDomPairs = {} -- Table of size 'ndim'.
   for d = 1, decomp:ndim() do
      local sz = self._cutsRange:shorten(d):volume() -- Number of skeleton subdomains.
      self._periodicDomPairs[d] = PairsVec(sz)       -- Store as integer pair.
   end

   -- (we can do this here as skeleton sub-domain IDs only depends on
   -- cuts and not on domain to be decomposed).
   local cutsIdxr = Range.makeColMajorGenIndexer(self._cutsRange)
   -- Loop over each direction, adding boundary sub-regions to
   -- pair-list to allow communication in periodic directions.
   for dir = 1, decomp:ndim() do
      -- Loop over "shortened" range which are basically
      -- sub-domains that lie on the lower domain boundary.
      local shortRange = self._cutsRange:shorten(dir)
      local c = 1
      for idx in shortRange:colMajorIter() do
         local idxp = idx:copy()
         idxp[dir] = idxp[dir]+decomp:cuts(dir)-1 -- Upper index.
        
         self._periodicDomPairs[dir][c].lower = cutsIdxr(idx)  -- Lower sub-domain on boundary.
         self._periodicDomPairs[dir][c].upper = cutsIdxr(idxp) -- Upper sub-domain on boundary.
         c = c+1
      end
   end   

   -- Identify the corner neighbors of skeleton sub-domains.
   local dimRemain = {}  -- Remaining dimensions when a dimension is removed.
   for d1 = 1, decomp:ndim() do
      dimRemain[d1] = {}
      for d2 = 1, decomp:ndim() do dimRemain[d1][d2] = d2 end
      table.remove(dimRemain[d1],d1)
   end
   -- The idea below is the following. For each direction get the subdomain coordinate
   -- of the face neighbor. Then the corner neighbors are -1/+1 step in the other dimensions
   -- from this face neighbor. At the end we remove the face neighbor from the list. 
   self._cornerNeigh = {}
   local idxLo, idxUp = {}, {}
   for dir = 1, decomp:ndim() do
      local shortRange = self._cutsRange:shorten(dir)
      self._cornerNeigh[dir] = {}
      local c = 0
      for idx in shortRange:colMajorIter() do
         for d=1,decomp:ndim() do idxLo[d], idxUp[d] = idx[d], idx[d] end
         idxUp[dir] = idxUp[dir]+decomp:cuts(dir)-1

         c = c+1
         self._cornerNeigh[dir][c] = {}
         table.insert(self._cornerNeigh[dir][c],{lower=lume.clone(idxLo), upper=lume.clone(idxUp), dirs={dir}})
         for _, dr in ipairs(dimRemain[dir]) do
            local cDoms = lume.deepclone(self._cornerNeigh[dir][c])
            for _, cd in pairs(cDoms) do
               local newNeighM, newNeighP = lume.deepclone(cd), lume.deepclone(cd)
               newNeighM["upper"][dr] = wrapToRange(cd["upper"][dr]-1,self._cutsRange:lower(dr),
                                                                      self._cutsRange:upper(dr))
               table.insert(newNeighM["dirs"],-dr)
               newNeighP["upper"][dr] = wrapToRange(cd["upper"][dr]+1,self._cutsRange:lower(dr),
                                                                      self._cutsRange:upper(dr))
               table.insert(newNeighP["dirs"], dr) 
               table.insert(self._cornerNeigh[dir][c],newNeighM)
               table.insert(self._cornerNeigh[dir][c],newNeighP)
            end
         end
         table.remove(self._cornerNeigh[dir][c],1)

         -- Translate the subdomain coordinate to a rank ID.
         for _, w in pairs(self._cornerNeigh[dir][c]) do
            w["lower"] = cutsIdxr(w["lower"])
            w["upper"] = cutsIdxr(w["upper"])
         end
      end
   end   

end

-- Set callable methods.
function DecomposedRange:commSet() return self._commSet end
function DecomposedRange:cutsRange() return self._cutsRange end
function DecomposedRange:ndim() return self._cutsRange:ndim() end
function DecomposedRange:cuts(dir) return self._cutsRange:upper(dir) end
function DecomposedRange:numSubDomains() return self._cutsRange:volume() end
function DecomposedRange:subDomain(k) return self._domains[k] end
function DecomposedRange:boundarySubDomainIds(dir) return self._periodicDomPairs[dir] end
function DecomposedRange:boundarySubDomainCornerIds(dir) return self._cornerNeigh[dir] end
function DecomposedRange:cutsIndexer() return Range.makeColMajorGenIndexer(self._cutsRange) end
function DecomposedRange:cutsInvIndexer() return Range.makeColMajorInvIndexer(self._cutsRange) end

-- CartProdDecomp --------------------------------------------------------------
--
-- Decomposition of grid into smaller boxes, specified by number of
-- cuts in each direction.
--------------------------------------------------------------------------------

local CartProdDecomp = Proto()
-- Constructor to make new cart-prod decomp.
function CartProdDecomp:init(tbl)
   local ones = {}
   for d = 1, #tbl.cuts do ones[d] = 1 end

   self._cutsRange = Range.Range(ones, tbl.cuts)

   local comm_host, comm_dev, comm_default = Mpi.COMM_WORLD, nil, Mpi.COMM_WORLD
   -- Use a different communicator if one is specified.
   if tbl.comm then comm_host, comm_default = tbl.comm, tbl.comm end
   -- Specify host and device communicators.
   if tbl.comms then
      comm_host, comm_dev = tbl.comms.host, tbl.comms.device
      comm_default = GKYL_USE_GPU and comm_dev or comm_host
   end

   -- Denote specific ranks from which writes should happen (defaults to all ranks).
   local writeRank = tbl.writeRank or Mpi.Comm_rank(comm_host)

   -- Check if total number of domains specified by 'cuts' matches
   -- number of MPI ranks in comm.
   if not tbl.__serTesting then   -- Skip check if just testing.
      if Mpi.Is_comm_valid(comm_host) then
	 if Mpi.Comm_size(comm_host) ~= self._cutsRange:volume() then
	    assert(false,
		   string.format(
		      "CartProdDecomp: Number of sub-domains (%d) and comm size (%d) must match",
		      self._cutsRange:volume(), Mpi.Comm_size(comm_host)))
	 end
      end
   end

   -- Store various communicators in a table.
   self._commSet = {
      default = comm_default,  writeRank = writeRank,
      host = comm_host,     device = comm_dev,
   }
end

-- Set callable methods.
function CartProdDecomp:commSet() return self._commSet end
function CartProdDecomp:ndim() return self._cutsRange:ndim() end
function CartProdDecomp:cuts(dir) return self._cutsRange:upper(dir) end

function CartProdDecomp:decompose(range) -- Decompose range.
   local decompRgn = DecomposedRange(self)
   
   local shapes, starts = {}, {} -- To store shapes, start in each direction.
   -- Compute shape and start index in each direction.
   for dir = 1, range:ndim() do
      shapes[dir] = Lin.IntVec(self:cuts(dir))
      local baseShape = math.floor(range:shape(dir)/self:cuts(dir))
      local remCells = range:shape(dir) % self:cuts(dir) -- Extra cells left over.
      for c = 1, self:cuts(dir) do
	 shapes[dir][c] = baseShape + (remCells>0 and 1 or 0) -- Add extra cell if any remain.
	 remCells = remCells-1
      end

      starts[dir] = Lin.IntVec(self:cuts(dir))
      starts[dir][1] = range:lower(dir)
      for c = 2, self:cuts(dir) do
	 starts[dir][c] = starts[dir][c-1] + shapes[dir][c-1]
      end
   end

   local c = 1
   -- Add regions.
   for idx in self._cutsRange:colMajorIter() do
      local l, u = {}, {}
      for dir = 1, range:ndim() do
	 l[dir] = starts[dir][idx[dir]]
	 u[dir] = l[dir]+shapes[dir][idx[dir]]-1
      end
      decompRgn._domains[c] = Range.Range(l, u)
      c = c+1
   end

   return decompRgn
end

function CartProdDecomp:childDecomp(keepDir)
   -- Obtain the ingredients needed to create a Cartesian decomposition of a lower
   -- dimension, keeping the dimensions in 'keepDir' (keepDir should be a table of
   -- the directions to be kept, e.g. {1,3}). This method doesn't create the new decomposition
   -- because that would require CartProdDecomp creating a new instance of itself.
   local parentDim, childDim = self._cutsRange:ndim(), #keepDir
   assert(childDim>0, "CartProdDecomp: the table of dimensions passed (keepDir) must have at lease one element.")
   assert(childDim<=parentDim, "CartProdDecomp: the table of dimensions passed (keepDir) must have fewer or same number of dimensions as parent decomposition.")
   local isDirKept = {}
   for d=1,parentDim do isDirKept[d] = false end
   for _, d in ipairs(keepDir) do isDirKept[d] = true end

   local parentCuts, childCuts = {}, {}
   for d=1,parentDim do parentCuts[d] = self:cuts(d) or 1 end
   for d=1,childDim do childCuts[d] = self:cuts(keepDir[d]) or 1 end

   local parentRank = Mpi.Comm_rank(self:commSet().host)
   local childRank, childWriteRank = parentRank, parentRank
   if (parentDim > 1) and (childDim < parentDim) then
      -- The following assumes a colum-major order distribution of MPI processes.
      if parentDim == 2 then
         if keepDir[1] == 1 then 
            childRank = math.floor(parentRank/parentCuts[1])
         elseif keepDir[1] == 2 then 
            childRank = parentRank % parentCuts[1]
         end
      else  -- This assumes that there's no decomposition along directions>3. 
         if childDim == 1 then
            if keepDir[1] == 1 then 
               childRank = math.floor(parentRank/parentCuts[1])
            elseif keepDir[1] == 2 then 
               childRank = parentRank % parentCuts[1] + parentCuts[1]*math.floor(parentRank/(parentCuts[1]*parentCuts[2]))
            elseif keepDir[1] == 3 then 
               childRank = parentRank % (parentCuts[1]*parentCuts[2])
            end
         elseif childDim == 2 then
            if ((keepDir[1] == 1) and (keepDir[2] == 2)) or ((keepDir[1] == 2) and (keepDir[2] == 1)) then
               childRank = math.floor(parentRank/(parentCuts[1]*parentCuts[2]))
            elseif ((keepDir[1] == 2) and (keepDir[2] == 3)) or ((keepDir[1] == 3) and (keepDir[2] == 2)) then
               childRank = (parentRank % (parentCuts[1]*parentCuts[2])) % parentCuts[1]
            elseif ((keepDir[1] == 3) and (keepDir[2] == 1)) or ((keepDir[1] == 1) and (keepDir[2] == 3)) then
               childRank = math.floor(parentRank/parentCuts[1]) % parentCuts[2]
            end
         end
      end
--      -- (UNFINISHED) The following assumes a row-major order distribution of MPI processes.
--      if parentDim == 2 then
--         if keepDir[1] == 1 then 
--            childRank = parentRank % parentCuts[2]
--         elseif keepDir[1] == 2 then 
--            childRank = math.floor(parentRank/parentCuts[1])
--         end
--      elseif parentDim == 3 then
--         if childDim == 1 then
--            if keepDir[1] == 1 then 
--               childRank = parentRank % (parentCuts[2]*parentCuts[3])
--            elseif keepDir[1] == 2 then 
--               childRank = parentRank
--            elseif keepDir[1] == 3 then 
--               childRank = parentRank
--            end
--         elseif childDim == 2 then
--            if ((keepDir[1] == 1) and (keepDir[2] == 2)) or ((keepDir[1] == 2) and (keepDir[2] == 1)) then
--               childRank = parentRank % parentCuts[3]
--            elseif ((keepDir[1] == 2) and (keepDir[2] == 3)) or ((keepDir[1] == 3) and (keepDir[2] == 2)) then
--               childRank = math.floor(parentRank/(parentCuts[2]*parentCuts[3]))
--            elseif ((keepDir[1] == 3) and (keepDir[2] == 1)) or ((keepDir[1] == 1) and (keepDir[2] == 3)) then
--            end
--         end
--      end

      local childWriteRank = -1   -- MF: I think this should be Mpi.PROC_NULL
      -- For now assume only the ranks with the lowest childRank do IO. 
      if childRank == 0 then childWriteRank = parentRank end
   end
   local childComm = Mpi.Comm_split(self:commSet().host, childRank, parentRank)

   return childComm, childWriteRank, childCuts
end

--------------------------------------------------------------------------------
--
-- Set utility functions to construct a reasonable decomposition when
-- user has not specified one explicitly
--
--------------------------------------------------------------------------------

-- check if numSubDoms^(1/ndim) is an integer
local function isEvenDecomp(ndim, numSubDoms, dofs)

   -- This code is tricky as round off errors can mess the simple (and
   -- obvious) test for even decomposition, i.e. math.ceil(v) =
   -- math.floor(v) where v = numSubDoms^(1/ndim)
   
   local c = math.pow(numSubDoms, 1/ndim)
   local cdown = math.floor(c)
   local cup = math.ceil(c)

   local cuts = {}
   if cdown^ndim == numSubDoms then
      for i = 1, ndim do cuts[i] = dofs[i]==1 and 1 or cdown end
      return true, cuts
   end
   if cup^ndim == numSubDoms then
      for i = 1, ndim do cuts[i] = dofs[i]==1 and 1 or cup end
      return true, cuts      
   end
   return false, cuts
end

local function calcCuts(ndim, numSubDoms, dofs)

   -- first check if we can create an "even" decomposition
   local evFlag, evCuts = isEvenDecomp(ndim, numSubDoms, dofs)
   if evFlag then return evCuts end
   
   -- Compute prime factorization of numSubDoms and then use this to
   -- compute "cuts". This algorithm may not produce the best possible
   -- decomposition, but it is hard to tell what "best" means. For
   -- most reasonable cases (say numSubDoms is not a prime or small
   -- multiple of a prime) the decompositions are probably okay.
   
   local factors = PrimeFactor.all(numSubDoms)
   local div, rem = math.floor(#factors/ndim), #factors%ndim
   local cuts, idx = {}, 1
   for d = 1, ndim do
      cuts[d] = 1
      for i = 1, div do
	 cuts[d] = cuts[d]*factors[idx]
	 idx = idx+1
      end
      if rem > 0 then
	 cuts[d] = cuts[d]*factors[idx]
	 idx = idx+1
	 rem = rem-1
      end
   end
   return cuts
end

-- Construct a reasonable decomposition (cuts) given dimensions,
-- number of sub-domains and cells along each dimension. The returned
-- cuts table can be used in CartProdDecomp to decompose a given range
-- object.
local function makeCuts(ndim, numSubDoms, cells)
   -- This initial sort is required as calcCuts does not return cuts
   -- in any particular order. The sort ensures that the rearrangement
   -- of cuts (done in the second for-loop) then matches order of
   -- cells table
   local rawCuts = calcCuts(ndim, numSubDoms, cells); table.sort(rawCuts)

   -- The idea here is that we need cuts arranged in same order as
   -- cell sizes. This is needed as "rawCuts" is not returned in any
   -- particular order, and needs to be rearranged to give a
   -- reasonable decomposition
   
   local cellsWithIdx = {}
   for i = 1, #cells do
      cellsWithIdx[i] = { cells[i], i }
   end
   table.sort(cellsWithIdx, function(l,r) return l[1]<r[1] end)

   local cuts = {}
   for i = 1, #cells do
      cuts[i] = rawCuts[cellsWithIdx[i][2]]
   end
   
   return cuts
end

return {
   CartProd = CartProdDecomp,
   makeCuts = makeCuts,
}
