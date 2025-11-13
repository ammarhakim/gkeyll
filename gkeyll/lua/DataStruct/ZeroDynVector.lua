-- Gkyl ------------------------------------------------------------------------
--
-- Lua interface to gkylzero's gkyl_dynvec structure
---
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local ffi = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"
local new, sizeof, typeof, metatype = xsys.from(ffi,
     "new, sizeof, typeof, metatype")

local Alloc = require "Lib.Alloc"
local cuda
require "Lib.ZeroUtil"

_M = {}

ffi.cdef [[
typedef struct gkyl_dynvec_tag* gkyl_dynvec;

struct gkyl_dynvec_tag {
  enum gkyl_elem_type type; // type of data stored in vector
  size_t elemsz, ncomp; // size of elements, number of 'components'

  size_t cloc; // current location to insert data into
  size_t csize; // current number of elements allocated

  size_t esznc; // elemsz*ncomp
  void *data; // pointer to data
  double *tm_mesh; // time stamps
  
  struct gkyl_ref_count ref_count;  
};

// Element type and number of components
struct gkyl_dynvec_etype_ncomp {
  enum gkyl_elem_type type; // type of data stored in vector
  size_t ncomp; // number of 'components'
};

/**
 * Create a new new dynvec. Delete using gkyl_dynvec_release method.
 * 
 * @param type Type of data in vector
 * @param ncomp Number of components
 * @return Newly allocated vector.
 */
gkyl_dynvec gkyl_dynvec_new(enum gkyl_elem_type type, size_t ncomp);

/**
 * Get element type stored in dynvec.
 *
 * @param vec Dynvec object
 * @return Element type
 */
int gkyl_dynvec_elem_type(gkyl_dynvec vec);

/**
 * Get number of components stored
 *
 * @param vec Dynvec object
 * @return Number of components
 */
int gkyl_dynvec_ncomp(gkyl_dynvec vec);

/**
 * Reserve @a rsize more elements so additional @a rsize append calls
 * do not require memory allocations.
 *
 * @param vec Vector to reserve data for
 * @param rsize Additional number of elements to reserve
 */
void gkyl_dynvec_reserve_more(gkyl_dynvec vec, size_t rsize);

/**
 * Append data to vector. You must ensure the data has the proper type
 * and ncomp number of elements.
 * 
 * @param vec Vector to append to
 * @param tm Time-stamp of data
 * @param data to append
 */
void gkyl_dynvec_append(gkyl_dynvec vec, double tm, const void *data);

/**
 * Get data at @a idx location. You must ensure the data has the
 * proper type and ncomp number of elements.
 * 
 * @param vec Vector
 * @param idx Index of data to fetch
 * @param data On return, data is copied in this buffer
 * @return If idx is not in range, false is returned
 */
bool gkyl_dynvec_get(const gkyl_dynvec vec, size_t idx, void *data);

/**
 * Get last appended data. You must ensure the data has the proper
 * type and ncomp number of elements.
 * 
 * @param vec Vector
 * @param data On return, data is copied in this buffer
 * @return If vector is empty, return false.
 */
bool gkyl_dynvec_getlast(const gkyl_dynvec vec, void *data);

/**
 * Get time-stamp for data at index idx.
 * 
 * @param vec Vector
 * @return Time-stamp of last data at idx
 */
double gkyl_dynvec_get_tm(const gkyl_dynvec vec, size_t idx);

/**
 * Get last appended data time-stamp.
 * 
 * @param vec Vector
 * @return Time-stamp of last appened data
 */
double gkyl_dynvec_getlast_tm(const gkyl_dynvec vec);

/**
 * Get size of dynvec.
 * 
 * @param vec Vector 
 * @return Number of elements in vector
 */
size_t gkyl_dynvec_size(const gkyl_dynvec vec);

/**
 * Get capacity of dynvec.
 *
 * @param vec Vector
 * @return Capacity (elements that can be stored witout reallocation)
 */
size_t gkyl_dynvec_capacity(const gkyl_dynvec vec);

/**
 * Get capacity of dynvec.
 * 
 * @param vec Vector 
 * @return Capacity of vector
 */
size_t gkyl_dynvec_capacity(const gkyl_dynvec vec);

/**
 * Clear contents of the vector
 *
 * @param vec Vector to clear
 */
void gkyl_dynvec_clear(gkyl_dynvec vec);

/**
 * Clear contents of the vector, but keep the last @a num inserted
 * elements, shifting them to the start of the vector. This is
 * typically useful when the vector has been written to file and the
 * data needs to be flushed, but the vector still is in use.
 *
 * If @a num is larger than the number of elements in the vector
 * then nothing is done.
 *
 * @param vec Vector to clear
 * @param num Number of final elements to keep.
 */
void gkyl_dynvec_clear_all_but(gkyl_dynvec vec, size_t num);

/**
 * Acquire a reference to the dynvec. Delete using gkyl_dynvec_release method.
 *
 * @param vec Vector to acquire reference from
 * @return Dynamic vector
 */
gkyl_dynvec gkyl_dynvec_acquire(const gkyl_dynvec vec);

/**
 * Write out dynvec to file. File is overwritten with new data, and
 * existing contents will lost.
 *
 * @param vec Vector to write
 * @param fname Name of output file.
 * @return 0 if succeeded.
 */
int gkyl_dynvec_write(const gkyl_dynvec vec, const char *fname);

/**
 * Write out dynvec to file. The dynvec is appened to the end of the
 * file if it already exists.
 *
 * @param vec Vector to write
 * @param fname Name of output file.
 * @return 0 if succeeded.
 */
int gkyl_dynvec_awrite(const gkyl_dynvec vec, const char *fname);

/**
 * Read number of components from the dynvec file.
 *
 * @param Name of input file
 * @return Element type and number of components
 */
struct gkyl_dynvec_etype_ncomp gkyl_dynvec_read_ncomp(const char *fname);

/**
 * Read dynvector from file, appending data to end of the
 * vector. Existing data in vector is retained, and read data is
 * appended.
 *
 * @param vec Vector to read into
 * @param fname Name of input file.
 */
bool gkyl_dynvec_read(gkyl_dynvec vec, const char *fname);

/**
 * Convert contents of dynvector to array. Time mesh is not copied to
 * the array but it returned as a seperate array. The input arrays
 * must be preallocated to be big enough to contain all the data.
 *
 * @param vec Dynvector to convert
 * @param tm_mesh On output, time-mesh of data
 * @param dyndata On output, data in dynamic array
 */
void gkyl_dynvec_to_array(const gkyl_dynvec vec, struct gkyl_array *tm_mesh,
  struct gkyl_array *dyndata);

/**
 * Release dynvec.
 *
 * @param vec Vector to release
 */
void gkyl_dynvec_release(gkyl_dynvec vec);

]]

local longSz = sizeof("long")

-- Various types for arrays of basic C-types
_M.int    = 'GKYL_INT'
_M.int64  = 'GKYL_INT64'
_M.float  = 'GKYL_FLOAT'
_M.double = 'GKYL_DOUBLE'
_M.long   = typeof('long')
_M.char   = typeof('char')

local function getArrayTypeCode(atype)
   if atype == typeof("int") then
      return 1
   elseif atype == typeof("int64") then
      return 2
   elseif atype == typeof("float") then
      return 3
   elseif atype == typeof("double") then
      return 4
   end
   return 42 -- user-defined type
end

local function getType(enum)
   if enum == 0 then
      return "int"
   elseif enum == 1 then
      return "int64"
   elseif enum == 2 then
      return "float"
   elseif enum == 3 then
      return "double"
   end
end

-- DynVector ctype
local DynVecCt = typeof("struct gkyl_dynvec_tag")

local dynvec_fn = {
   elemType = function(self)
      return ffiC.gkyl_dynvec_elem_type(self)
   end,
   numComponents = function(self)
      return ffiC.gkyl_dynvec_ncomp(self)
   end,
   appendData = function(self, tm, vals)
      local f = ffi.new(getType(self:elemType()) .. "[?]", self:numComponents())
      for i = 1, self:numComponents() do
	 f[i-1] = vals[i]
      end
      ffiC.gkyl_dynvec_append(self, tm, f)
   end,
   write = function(self, fName, tmStamp, frame)
      local fullNm = GKYL_OUT_PREFIX .. "_" .. fName
      return ffiC.gkyl_dynvec_write(self, fullNm)
   end,
}

local dynvec_mt = {
   __new = function (self, atype, ncomp)
      return ffi.gc(ffiC.gkyl_dynvec_new(atype, ncomp), ffiC.gkyl_dynvec_release)
   end,
   __index = dynvec_fn
}
local DynVecCtor = metatype(DynVecCt, dynvec_mt)

_M.DynVector = function (tbl)
   return DynVecCtor(_M.double, tbl.numComponents)
end

_M.DynVectorGen = function (atype, ncomp)
   return DynVecCtor(atype, ncomp)
end

_M.DynVectorFromFile = function(fName)
   local enc = ffi.C.gkyl_dynvec_read_ncomp(fName)
   local dv = _M.DynVectorGen(enc.type, enc.ncomp)
   ffi.C.gkyl_dynvec_read(dv,  fName)
   return dv
end

return _M
