-- Gkyl ------------------------------------------------------------------------
--
-- Vectors, matrices and linear algebra
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- Vector ----------------------------------------------------------------------
--
-- A chunk of memory with length stored. One can also create vector of
-- structures (can not be VLS or VLA)
--------------------------------------------------------------------------------

-- Allocate vector with type "ct" and size "n"
local function vec_alloc(ct, n)
   local v = new(ct, n)
   v._n, v._p = n, v._a
   return v
end

local function vec_alloc_from_data(ct, n, data)
   local v = new(ct, 1) -- allocate a single element (this is never used)
   v._n, v._p = n, data
   return v
end

-- copy vector "x" into vector "y"
local function vec_memcpy(y, x)
  copy(y._p, x._p, sizeof(y:elemType())*y._n)
end

-- Vector meta-object creator
local function new_vec_ct(elct)
   local vec_mf = {
      elemType = function (self)	 
	 return elct
      end,
      data = function (self)
	 return self._p
      end,
      copy = function (self)
	 local v = vec_alloc(self, self._n)
	 vec_memcpy(v, self)
	 return v
      end,
      copyInto = function (self, vecOut)
	 vec_memcpy(vecOut, self)
      end,
      setValues = function (self, vals)
	 local n = math.min(#vals, #self)
	 for i = 1, n do
	    self._p[i-1] = vals[i]
	 end
	 return self
      end,
   }
   local vec_mt = {
      __new = function(ct, n, data)
	 return data and vec_alloc_from_data(ct, n, data) or vec_alloc(ct, n)
      end,
      __len = function (self)
	 return self._n
      end,
      __index = function (self, k)
	 if type(k) == "number" then
	    return self._p[k-1]
	 else
	    return vec_mf[k]
	 end
      end,
      __newindex = function (self, k, v)
	 self._p[k-1] = v
      end,
   }
   local ct = typeof("struct { int _n; $* _p; $ _a[?]; }", elct, elct)
   return metatype(ct, vec_mt)
end

-- Matrix ----------------------------------------------------------------------
--
-- A two-dimensional matrix object, stored in row major order
--------------------------------------------------------------------------------

-- Allocate matrix with type "ct" and size "n,m"
local function mat_alloc(ct, n, m)
   local v = new(ct, n*m)
   v._n, v._m, v._p = n, m, v._a
   return v
end

local function mat_alloc_from_data(ct, n, m, data)
   local v = new(ct, 1) -- allocate a single element (this is never used)
   v._n, v._m, v._p = n, m, data
   return v
end

-- copy matrix "x" into matrix "y"
local function mat_memcpy(y, x)
   copy(y._p, x._p, sizeof(y:elemType())*y._n*y._m)
end

-- Accessor object for values in a column
local function new_mat_row_ct(elct)
   local mat_row_mt = {
      __index = function (self, k)
	 return self._cdata[k-1]
      end,
      __newindex = function (self, k, v)
	 self._cdata[k-1] = v
      end,
   }
   return metatype(typeof("struct { int ncols; $* _cdata; }", elct), mat_row_mt)
end

-- Matrix meta-object creator
local function new_mat_ct(elct)
   local rowct = new_mat_row_ct(elct)
   local mat_mf = {
      elemType = function (self)
	 return elct
      end,
      data = function (self)
	 return self._p
      end,
      copy = function (self)
	 local m = mat_alloc(self, self._n, self._m)
	 mat_memcpy(m, self)
	 return m
      end,
      copyInto = function (self, matOut)
	 mat_memcpy(matOut, self)
      end,
      numRows = function (self)
	 return self._n
      end,
      numCols = function (self)
	 return self._m
      end,
      g = function(self, i, j)
	 return self._p[(i-1)*self._m + (j-1)]
      end,
      s = function(self, i, j, v)
	 self._p[(i-1)*self._m + (j-1)] = v
      end,
      getRow = function(self, i)
         return self._p+(i-1)*self._m
      end,
   }
   local mat_mt = {
      __new = function(ct, n, m, data)
	 return data and mat_alloc_from_data(ct, n, m, data) or mat_alloc(ct, n, m)
      end,
      __index = function (self, k)
	 if type(k) == "number" then
	    return rowct(self._m, self._p+(k-1)*self._m)
	 else
	    return mat_mf[k]
	 end
      end,
   }
   local ct = typeof("struct { int _n, _m; $* _p; $ _a[?]; }", elct, elct)
   return metatype(ct, mat_mt)
end

return {
   new_vec_ct = new_vec_ct,
   Vec = new_vec_ct(typeof("double")),
   UIntVec = new_vec_ct(typeof("unsigned")),
   IntVec = new_vec_ct(typeof("int")),
   BoolVec = new_vec_ct(typeof("bool")),
   FloatVec = new_vec_ct(typeof("float")),
   ComplexVec = new_vec_ct(typeof("double complex")),
   FloatComplexVec = new_vec_ct(typeof("float complex")),
   CharVec = new_vec_ct(typeof("char")),
   
   -- some specialized types
   UInt64Vec = new_vec_ct(typeof("uint64_t")),
   Int64Vec = new_vec_ct(typeof("int64_t")),
   UInt32Vec = new_vec_ct(typeof("uint32_t")),
   Int32ec = new_vec_ct(typeof("int32_t")),   

   new_mat_ct = new_mat_ct,
   Mat = new_mat_ct(typeof("double")),
   IntMat = new_mat_ct(typeof("int")),
   FloatMat = new_mat_ct(typeof("float")),
   ComplexMat = new_mat_ct(typeof("double complex")),
   FloatComplexMat = new_mat_ct(typeof("float complex")),
}
