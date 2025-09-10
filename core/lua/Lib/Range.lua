-- Gkyl ------------------------------------------------------------------------
--
-- Range object
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- Gkyl libraries
local Lin = require "Lib.Linalg"

local _M = {}
local _P = {} -- private namespace

-- these are used to fetch consistent indexer and iterators
_M.rowMajor = 1
_M.colMajor = 2

-- Range ----------------------------------------------------------------------
--
-- A range object, representing a ndim integer index set. 
--------------------------------------------------------------------------------

ffi.cdef [[ 
/**
 * Range object, representing an N-dimensional integer index
 * set. Lower and upper limits are inclusive.
 */
struct gkyl_range {
  int _ndim; // number of dimension
  int _lower[7]; // lower bound
  int _upper[7]; // upper bound (inclusive)
  long _volume; // total volume of range
    
  // do not access directly
  uint32_t flags; // Flags for internal use
  int ilo[7]; // for use in inverse indexer
  long ac[8]; // coefficients for indexing
  long iac[8]; // for use in sub-range inverse indexer
  long linIdxZero; // linear index of {0,0,...}
  int nsplit, tid; // number of splits, split ID

  // FOR CUDA ONLY
  int nthreads, nblocks; // CUDA kernel launch specifiers for range-based ops
};

/**
 * Initialize new range object.
 *
 * @param rng Range object to initialize
 * @param ndim Dimension of range to create.
 * @param lower Lower indices of range
 * @param upper Upper indices of range
 */
void gkyl_range_init(struct gkyl_range *rng, int ndim,
  const int *lower, const int *upper);

/**
 * Return 1 if range is a sub-range.
 *
 * @param rng Range object
 * @return 1 if true, 0 otherwise
 */
int gkyl_range_is_sub_range(const struct gkyl_range *rng);

/**
 * Create a sub-range from a given range. The sub-range must be fully
 * contained in the parent range or else it will be truncated. The
 * sub-range and the parent range will returns the same linear index
 * for a given index.
 *
 * @param rng New range object to initialize
 * @param bigrng Parent range object 
 * @param sublower Lower indices of sub-range
 * @param subupper Upper indices of sub-range
 */
void gkyl_sub_range_init(struct gkyl_range *rng,
  const struct gkyl_range *bigrng, const int *sublower, const int *subupper);

/**
 * Return range which has some directions removed by setting the index
 * in those directions to fixed values. The "deflated" range has lower
 * dimension than the parent 'rng' object. The indexing into the
 * returned lower dimensional range gives the same index as the
 * corresponding location in the parent range (with the missing
 * indices set to 'locDir[dir]').
 *
 * @param srng Deflated range.
 * @param rng Range object to deflate
 * @param remDir 'ndim' Array of flags: 0 to keep direction, 1 to remove
 * @param loc Index to set removed direction.
 */
void gkyl_range_deflate(struct gkyl_range* srng,
  const struct gkyl_range* rng, const int *remDir, const int *locDir);
]]
local rTy = typeof("struct gkyl_range")
local rSz = sizeof(typeof("struct gkyl_range"))

-- generic iterator function creator: only difference between row- and
-- col-major order is the order in which the indices are incremented
local function make_range_iter(si, ei, incr)
   return function (iterState)
      if iterState.isEmpty then return nil end -- nothing to do for empty range

      if iterState.numBumps >= iterState.maxBumps then
	 return nil -- if bumped enough times aborting
      else
	 iterState.numBumps = iterState.numBumps+1
      end
      
      local idx = iterState.currIdx
      if iterState.isFirst then
	 iterState.isFirst = false
	 return idx
      end
      
      local range = iterState.range
      for dir = si, ei, incr do
	 idx[dir] = idx[dir]+1
	 if idx[dir] > range:upper(dir) then
	    idx[dir] = range:lower(dir)
	 else
	    return idx
	 end
      end
   end
end

-- compare a and b
local function cmpInt6(a, b)
   return a[0] == b[0] and a[1] == b[1] and a[2] == b[2] and a[3] == b[3] and a[4] == b[4] and a[5] == b[5]
end

-- Range object meta-type
local range_mt = {
   __new = function (self, lower, upper)
      -- lower and upper are tables of integers. Indices are inclusive
      local r = new(self)
      for d = 0, 5 do
      	 r._lower[d], r._upper[d] = 0, 0
      end
      local ndim = #lower
      r._ndim = ndim
      assert(r._ndim < 7, "Range.new: Only objects upto 6D are supported")
      for d = 1, #lower do
	 r._lower[d-1] = lower[d]
	 r._upper[d-1] = upper[d]
      end
      for d = 1, ndim do
	 -- adjust to give zero volume range if upper is less than lower
	 if r._upper[d-1] < r._lower[d-1] then
	    r._upper[d-1] = r._lower[d-1]-1
	 end
      end
      ffiC.gkyl_range_init(r, r._ndim, r._lower, r._upper)
      return r
   end,
   __eq = function (self, r)
      if self._ndim ~= r._ndim then return false end
      return cmpInt6(self._lower, r._lower) and cmpInt6(self._upper, r._upper)
   end,
   __index = {
      ndim = function (self)
	 return self._ndim
      end,
      lower = function (self, dir)
	 return self._lower[dir-1]
      end,
      upper = function (self, dir)
	 return self._upper[dir-1]
      end,
      copy = function (self, rng)
         ffi.copy(self, rng, rSz)
      end,
      lowerAsVec = function (self)
	 local v = Lin.IntVec(self:ndim())
	 for dir = 1, self:ndim() do  v[dir] = self._lower[dir-1] end
	 return v
      end,
      upperAsVec = function (self)
	 local v = Lin.IntVec(self:ndim())
	 for dir = 1, self:ndim() do v[dir] = self._upper[dir-1] end
	 return v
      end,
      shape = function (self, dir)
         local shapeOut
         if dir then
            shapeOut = math.max(0, self._upper[dir-1]-self._lower[dir-1]+1)
         else
            shapeOut = {}
            for dir = 1, self._ndim do shapeOut[dir] = math.max(0, self._upper[dir-1]-self._lower[dir-1]+1) end
         end
         return shapeOut
      end,
      volume = function (self)
	 local v = 1
	 for dir = 1, self._ndim do
	    v = v*math.max(0, self._upper[dir-1]-self._lower[dir-1]+1)
	 end
	 return v
      end,
      extend = function (self, lExt, uExt)
	 local r = new(rTy)
	 r._ndim = self:ndim()
	 for dir = 1, self:ndim() do
	    r._lower[dir-1], r._upper[dir-1] = self:lower(dir)-lExt, self:upper(dir)+uExt
	 end
         ffiC.gkyl_range_init(r, r._ndim, r._lower, r._upper)
	 return r
      end,
      extendDir = function (self, extDir, lExt, uExt)
	 local r = new(rTy)
	 r._ndim = self:ndim()
	 for dir = 1, self:ndim() do
	    r._lower[dir-1], r._upper[dir-1] = self:lower(dir), self:upper(dir)
	 end
	 r._lower[extDir-1], r._upper[extDir-1] = self:lower(extDir)-lExt, self:upper(extDir)+uExt
         ffiC.gkyl_range_init(r, r._ndim, r._lower, r._upper)
	 return r
      end,
      extendDirs = function (self, dirList, lExt, uExt)
	 local r = new(rTy)
	 r._ndim = self:ndim()
	 for dir = 1, self:ndim() do
	    r._lower[dir-1], r._upper[dir-1] = self:lower(dir), self:upper(dir)
	 end
	 for d = 1, #dirList do
	    local extDir = dirList[d]
	    r._lower[extDir-1], r._upper[extDir-1] = self:lower(extDir)-lExt, self:upper(extDir)+uExt
	 end
         ffiC.gkyl_range_init(r, r._ndim, r._lower, r._upper)
	 return r
      end,
      subRange = function(self, sublower, subupper)
	 local r = new(rTy)
         for d = 1, #sublower do
            r._lower[d-1] = sublower[d]
            r._upper[d-1] = subupper[d]
         end
         ffiC.gkyl_sub_range_init(r, self, r._lower, r._upper)
	 return r
      end,
      isSubRange = function(self)
         return ffiC.gkyl_range_is_sub_range(self)
      end, 
      selectFirst = function (self, n)
	 local r = new(rTy)
	 r._ndim = n
	 for dir = 1, n do
	    r._lower[dir-1], r._upper[dir-1] = self:lower(dir), self:upper(dir)
	 end
         ffiC.gkyl_range_init(r, r._ndim, r._lower, r._upper)
	 return r	 
      end,
      selectLast = function (self, n)
	 local r = new(rTy)
	 r._ndim = n
	 for dir = 1, n do
	    local i = self:ndim()-n+dir
	    r._lower[dir-1], r._upper[dir-1] = self:lower(i), self:upper(i)
	 end
         ffiC.gkyl_range_init(r, r._ndim, r._lower, r._upper)
	 return r	 
      end,      
      lowerSkin = function (self, dir, nGhost)
	 local r = new(rTy)
	 r._ndim = self:ndim()
	 for d = 1, self:ndim() do
	    r._lower[d-1], r._upper[d-1] = self:lower(d), self:upper(d)
	 end
	 r._upper[dir-1] = self:lower(dir)+nGhost-1
         ffiC.gkyl_sub_range_init(r, self, r._lower, r._upper)
	 return r
      end,
      upperSkin = function (self, dir, nGhost)
	 local r = new(rTy)
	 r._ndim = self:ndim()
	 for d = 1, self:ndim() do
	    r._lower[d-1], r._upper[d-1] = self:lower(d), self:upper(d)
	 end
	 r._lower[dir-1] = self:upper(dir)-nGhost+1
         ffiC.gkyl_sub_range_init(r, self, r._lower, r._upper)
	 return r
      end,
      lowerGhost = function (self, dir, nGhost)
	 local r = new(rTy)
	 r._ndim = self:ndim()
	 for d = 1, self:ndim() do
	    r._lower[d-1], r._upper[d-1] = self:lower(d), self:upper(d)
	 end
	 r._lower[dir-1], r._upper[dir-1] = self:lower(dir)-nGhost, self:lower(dir)-1
         ffiC.gkyl_range_init(r, r._ndim, r._lower, r._upper)
	 return r
      end,
      upperGhost = function (self, dir, nGhost)
	 local r = new(rTy)
	 r._ndim = self:ndim()
	 for d = 1, self:ndim() do
	    r._lower[d-1], r._upper[d-1] = self:lower(d), self:upper(d)
	 end
	 r._lower[dir-1], r._upper[dir-1] = self:upper(dir)+1, self:upper(dir)+nGhost
         ffiC.gkyl_range_init(r, r._ndim, r._lower, r._upper)
	 return r
      end,
      deflate = function (self, removeDir, locationInDir)
	 local r = new(rTy)
	 local remDir, locDir = Lin.IntVec(self:ndim()), Lin.IntVec(self:ndim())
         for d = 1, self:ndim() do
	    remDir[d], locDir[d] = removeDir[d], locationInDir[d]
         end
         ffiC.gkyl_range_deflate(r, self, remDir:data(), locDir:data())
         return r
      end,
      shorten = function (self, dir, len)
         len = len or 1
	 local r = new(rTy)
	 r._ndim = self:ndim()
	 for d = 1, self:ndim() do
	    r._lower[d-1], r._upper[d-1] = self:lower(d), self:upper(d)
	 end
	 r._upper[dir-1] = r._lower[dir-1]+len-1
         ffiC.gkyl_range_init(r, r._ndim, r._lower, r._upper)
	 return r
      end,
      shortenFromBelow = function (self, dir, len)
	 if len == nil then len = 1 end
	 local r = new(rTy)
	 r._ndim = self:ndim()
	 for d = 1, self:ndim() do
	    r._lower[d-1], r._upper[d-1] = self:lower(d), self:upper(d)
	 end
	 r._lower[dir-1] = r._upper[dir-1]-(len-1)
         ffiC.gkyl_range_init(r, r._ndim, r._lower, r._upper)
	 return r
      end,
      shift = function (self, offsets)
	 local r = new(rTy)
	 r._ndim = self:ndim()
	 for d = 1, self:ndim() do
	    r._lower[d-1], r._upper[d-1] = self:lower(d)+offsets[d], self:upper(d)+offsets[d]
	 end
         ffiC.gkyl_range_init(r, r._ndim, r._lower, r._upper)
	 return r
      end,
      shiftInDir = function (self, dir, offset)
	 local r = new(rTy)
	 r._ndim = self:ndim()
	 for d = 1, self:ndim() do
	    r._lower[d-1], r._upper[d-1] = self:lower(d), self:upper(d)
	 end	 
	 r._lower[dir-1], r._upper[dir-1] = self:lower(dir)+offset, self:upper(dir)+offset
         ffiC.gkyl_range_init(r, r._ndim, r._lower, r._upper)
	 return r
      end,
      intersect = function (self, rgn)
	 local lo, up = Lin.IntVec(self:ndim()), Lin.IntVec(self:ndim())
	 for d = 1, self:ndim() do
	    lo[d] = math.max(self:lower(d), rgn:lower(d))
	    up[d] = math.min(self:upper(d), rgn:upper(d))
	 end
	 return _M.Range(lo, up)
      end,
      isIntersectionEmpty = function (self, rgn)
	 for d = 1, self:ndim() do
	    if math.min(self:upper(d), rgn:upper(d)) < math.max(self:lower(d), rgn:lower(d)) then
	       return true
	    end
	 end
	 return false
      end,
      isDifferenceInDirEmpty = function (self, rgn, dir)
         -- Difference = relative complement = self\rgn. See the difference method below.
         -- The difference is empty if the ranges are the same.
         if self:lower(dir)==rgn:lower(dir) and self:upper(dir)==rgn:upper(dir) then
            return true
         else
            return false
	 end
      end,
      isDifferenceEmpty = function (self, rgn)
         -- Difference = relative complement = self\rgn. See the difference method below.
         -- The difference is empty if the ranges are the same.
         local isSame = true
	 for d = 1, self:ndim() do
            isSame = isSame and self:isDifferenceInDirEmpty(rgn, d)
	 end
	 return isSame
      end,
      difference = function (self, rgn)
         -- Compute the relative complement, or difference, self\rgn: elements in self but not in rgn.
         -- Or also seen as self - intersection(self,rgn).
         if self:isIntersectionEmpty(rgn) then return _M.Range(self:lowerAsVec(),self:upperAsVec()) end
         if self:isDifferenceEmpty(rgn) then return nil end
         local lo, up = Lin.IntVec(self:ndim()), Lin.IntVec(self:ndim())
         for d = 1, self:ndim() do
            local myLo, myUp = self:lower(d), self:upper(d)
            if self:isDifferenceInDirEmpty(rgn,d) then
               lo[d], up[d] = myLo, myUp
            else
               local interLo, interUp = math.max(myLo, rgn:lower(d)), math.min(myUp, rgn:upper(d))
               lo[d] = myLo<interLo and myLo or interUp+1
               up[d] = myLo<interLo and interLo-1 or myUp
            end
         end
         return _M.Range(lo, up)
      end,
      contains = function (self, idx)
	 for d = 1, self:ndim() do
	    if idx[d]<self:lower(d) or idx[d]>self:upper(d) then
	       return false
	    end
	 end
	 return true
      end,
      _iter = function (self, iter_func, idxStart, maxBumps)
	  -- package up iterator state into table
	 local iterState = {
	    isFirst  = true, numBumps = 0,
	    isEmpty  = self:volume() == 0 and true or false,
	    maxBumps = maxBumps,
	    range    = self
	 }
	 iterState.currIdx = Lin.IntVec(self:ndim())
	 for dir = 1, self:ndim() do
	    iterState.currIdx[dir] = idxStart[dir]
	 end

	 return iter_func, iterState
      end,
      colMajorIter = function (self, idxStartIn, maxBump)
	 local idxStart = idxStartIn
	 if idxStart == nil then
	    idxStart = Lin.IntVec(self:ndim())
	    for dir = 1, self:ndim() do
	       idxStart[dir] = self:lower(dir)
	    end
	 end
	 return self:_iter( make_range_iter(1, self:ndim(), 1), idxStart, maxBump and maxBump or self:volume() )
      end,
      rowMajorIter = function (self, idxStartIn, maxBump)
	 local idxStart = idxStartIn
	 if idxStart == nil then
	    idxStart = Lin.IntVec(self:ndim())
	    for dir = 1, self:ndim() do
	       idxStart[dir] = self:lower(dir)
	    end
	 end	 
	 return self:_iter( make_range_iter(self:ndim(), 1, -1), idxStart, maxBump and maxBump or self:volume() )
      end,
      iter = function (self, ordering, idxStartIn, maxBump)
	 if ordering == _M.rowMajor then
	    return self:rowMajorIter(idxStartIn, maxBump)
	 else
	    return self:colMajorIter(idxStartIn, maxBump)
	 end
      end,
      indexer = function (self, ordering)
	 return _M.makeIndexer(ordering, self)
      end,
      genIndexer = function (self, ordering)
	 return _M.makeGenIndexer(ordering, self)
      end,
   }
}
-- construct Range object, attaching meta-type to it
_M.Range = metatype(rTy, range_mt)

-- Indexers --------------------------------------------------------------------
--
-- Linear indexers, mapping n-dimensional index to a linear index
--------------------------------------------------------------------------------

-- The following set of functions takes explicit N-dimensional
-- (i,j,...) indices and map it to an integer.
local function getIndex1(ac, i1)
   return ac[0] + i1*ac[1]
end
local function getIndex2(ac, i1, i2)
   return ac[0]+i1*ac[1]+i2*ac[2]
end
local function getIndex3(ac, i1, i2, i3)
   return ac[0]+i1*ac[1]+i2*ac[2]+i3*ac[3]
end
local function getIndex4(ac, i1, i2, i3, i4)
   return ac[0]+i1*ac[1]+i2*ac[2]+i3*ac[3]+i4*ac[4]
end
local function getIndex5(ac, i1, i2, i3, i4, i5)
   return ac[0]+i1*ac[1]+i2*ac[2]+i3*ac[3]+i4*ac[4]+i5*ac[5]
end
local function getIndex6(ac, i1, i2, i3, i4, i5, i6)
   return ac[0]+i1*ac[1]+i2*ac[2]+i3*ac[3]+i4*ac[4]+i5*ac[5]+i6*ac[6]
end
-- package these up into a table
local indexerFunctions = {
   getIndex1, getIndex2, getIndex3, getIndex4, getIndex5, getIndex6
}

-- The following set of functions take a index with N-dim entries and
-- map it to an integer.
local function getGenIndex1(ac, idx)
   return ac[0] + idx[1]*ac[1]
end
local function getGenIndex2(ac, idx)
   return ac[0]+idx[1]*ac[1]+idx[2]*ac[2]
end
local function getGenIndex3(ac, idx)
   return ac[0]+idx[1]*ac[1]+idx[2]*ac[2]+idx[3]*ac[3]
end
local function getGenIndex4(ac, idx)
   return ac[0]+idx[1]*ac[1]+idx[2]*ac[2]+idx[3]*ac[3]+idx[4]*ac[4]
end
local function getGenIndex5(ac, idx)
   return ac[0]+idx[1]*ac[1]+idx[2]*ac[2]+idx[3]*ac[3]+idx[4]*ac[4]+idx[5]*ac[5]
end
local function getGenIndex6(ac, idx)
   return ac[0]+idx[1]*ac[1]+idx[2]*ac[2]+idx[3]*ac[3]+idx[4]*ac[4]+idx[5]*ac[5]+idx[6]*ac[6]
end
-- package these up into a table
local genIndexerFunctions = {
   getGenIndex1, getGenIndex2, getGenIndex3, getGenIndex4, getGenIndex5, getGenIndex6
}

-- create coefficients for row-major indexer  given "range" object
local function calcRowMajorIndexerCoeff(range)
   local ac = new("int[7]")
   local ndim = range:ndim()
   ac[ndim] = 1
   for i = ndim-1, 1, -1 do
      ac[i] = ac[i+1]*range:shape(i+1)
   end
   local start = 0
   for i = 1, ndim do
      start = start + ac[i]*range:lower(i)
   end
   ac[0] = 1-start
   return ac
end
_P.calcRowMajorIndexerCoeff = calcRowMajorIndexerCoeff

-- create coefficients for column-major indexer  given "range" object
local function calcColMajorIndexerCoeff(range)
   local ac = new("int[7]")
   local ndim = range:ndim()
   ac[1] = 1
   for i = 2, ndim do
      ac[i] = ac[i-1]*range:shape(i-1)
   end
   local start = 0
   for i = 1, ndim do
      start = start + ac[i]*range:lower(i)
   end
   ac[0] = 1-start
   return ac
end
_P.calcColMajorIndexerCoeff = calcColMajorIndexerCoeff

-- Following functions return an indexer function given a range
-- object. The returned function takes a (i,j,...) index.
function _M.makeRowMajorIndexer(range)
   local ac = calcRowMajorIndexerCoeff(range)
   local idxr = indexerFunctions[range:ndim()]
   return function (...)
      return idxr(ac, ...)
   end
end
function _M.makeColMajorIndexer(range)
   local ac = calcColMajorIndexerCoeff(range)
   local idxr = indexerFunctions[range:ndim()]
   return function (...)
      return idxr(ac, ...)
   end
end

-- Following functions return an indexer function given a range
-- object. The returned function takes a single idx parameter.
function _M.makeRowMajorGenIndexer(range)
   local ac = calcRowMajorIndexerCoeff(range)
   local idxr = genIndexerFunctions[range:ndim()]
   return function (idx)
      return idxr(ac, idx)
   end
end
function _M.makeColMajorGenIndexer(range)
   local ac = calcColMajorIndexerCoeff(range)
   local idxr = genIndexerFunctions[range:ndim()]
   return function (idx)
      return idxr(ac, idx)
   end
end

-- generic functions that dispatch on layout
function _M.makeIndexer(ordering, range)
   if ordering == _M.rowMajor then
      return _M.makeRowMajorIndexer(range)
   else
      return _M.makeColMajorIndexer(range)
   end
end
function _M.makeGenIndexer(ordering, range)
   if ordering == _M.rowMajor then
      return _M.makeRowMajorGenIndexer(range)
   else
      return _M.makeColMajorGenIndexer(range)
   end
end

-- InvIndexers ----------------------------------------------------------------
--
-- Inverse indexers, mapping a linear index to n-dimensional index
--------------------------------------------------------------------------------

ffi.cdef [[ 
  typedef struct { int quot, rem; } div_t;
  div_t div(int n, int d); 
]]
-- function to compute quotient and remainder from integer division
local function div(n, d)
   local r = ffiC.div(n, d) -- call std C library
   return r.quot, r.rem
end

-- inverse function is implemented as a template to unroll inner
-- loop. (We need to pass the 'div' function as the generated code
-- does not have access to the enclosing scope)
local getRowMajorInvIndexerTempl = xsys.template([[
return function (ac, range, loc, idx, div)
   local n = loc-1
|   for i = 1, NDIM do
      local q, r = div(n, ac[${i}])
      idx[${i}] = q+range:lower(${i})
      n = r
|   end
end
]])

-- package up into a table
local invRowIndexerFunctions = {
   loadstring(getRowMajorInvIndexerTempl { NDIM = 1 })(),
   loadstring(getRowMajorInvIndexerTempl { NDIM = 2 })(),
   loadstring(getRowMajorInvIndexerTempl { NDIM = 3 })(),
   loadstring(getRowMajorInvIndexerTempl { NDIM = 4 })(),
   loadstring(getRowMajorInvIndexerTempl { NDIM = 5 })(),
   loadstring(getRowMajorInvIndexerTempl { NDIM = 6 })()
}

-- compute index given linear location in range (NOT USED BUT PROVIDED
-- FOR REFERENCE TO UNDERSTAND THE TEMPLATE)
local function getRowMajorInvIndexer(ac, range, loc, idx)
   local n = loc-1
   for i = 1, range:ndim() do
      local q, r = div(n, ac[i])
      idx[i] = q+range:lower(i)
      n = r
   end
end

-----------

-- inverse function is implemented as a template to unroll inner loop.
-- (We need to pass the 'div' function as the generated code does not
-- have access to the enclosing scope)
local getColMajorInvIndexerTempl = xsys.template([[
return function (ac, range, loc, idx, div)
   local n = loc-1
|   for i = NDIM, 1, -1 do
      local q, r = div(n, ac[${i}])
      idx[${i}] = q+range:lower(${i})
      n = r
|   end
end
]])

-- package up into a table
local invColIndexerFunctions = {
   loadstring(getColMajorInvIndexerTempl { NDIM = 1 })(),
   loadstring(getColMajorInvIndexerTempl { NDIM = 2 })(),
   loadstring(getColMajorInvIndexerTempl { NDIM = 3 })(),
   loadstring(getColMajorInvIndexerTempl { NDIM = 4 })(),
   loadstring(getColMajorInvIndexerTempl { NDIM = 5 })(),
   loadstring(getColMajorInvIndexerTempl { NDIM = 6 })()
}

-- compute index given linear location in range (NOT USED BUT PROVIDED
-- FOR REFERENCE TO UNDERSTAND THE TEMPLATE)
local function getColMajorInvIndexer(ac, range, loc, idx)
   local n = loc-1
   for i = range:ndim(), 1, -1 do
      local q, r = div(n, ac[i])
      idx[i] = q+range:lower(i)
      n = r
   end
end

-- Following functions return inverse indexers: i.e. give a scalar (1
-- based) location, they set an (NDIM) index into the range
function _M.makeRowMajorInvIndexer(range)
   local ac = calcRowMajorIndexerCoeff(range)
   local invFunc = invRowIndexerFunctions[range:ndim()]
   return function (loc, idx)
      return invFunc(ac, range, loc, idx, div)
   end
end
function _M.makeColMajorInvIndexer(range)
   local ac = calcColMajorIndexerCoeff(range)
   local invFunc = invColIndexerFunctions[range:ndim()]
   return function (loc, idx)
      return invFunc(ac, range, loc, idx, div)
   end
end


return _M
