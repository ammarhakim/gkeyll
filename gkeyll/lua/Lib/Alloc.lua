-- Gkyl ------------------------------------------------------------------------
--
-- Allocators for memory not managed by Lua garbage-collector: allows
-- allocation of blocks not bound by Lua GC limitations.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

ffi.cdef [[
  void* malloc(size_t size);
  void* calloc(size_t nitems, size_t size);
  void* realloc(void *ptr, size_t new_size);
  void free(void *ptr);
]]

local function malloc(sz)
   local d = ffiC.malloc(sz)
   assert(d, "Alloc.malloc: Unable to allocate memory!")
   return d
end
local function calloc(nitems, sz)
   local d = ffiC.calloc(nitems, sz)
   assert(d, "Alloc.calloc: Unable to allocate memory!")
   return d
end
local function realloc(dold, sz)
   local d = ffiC.realloc(dold, sz)
   assert(d, "Alloc.realloc: Unable to allocate memory!")
   return d
end
local function free(d)
   if d then
      ffiC.free(d); d = nil
   end
end

-- Alloc ----------------------------------------------------------------------
--
-- A dynamically sizable 1D array type
--------------------------------------------------------------------------------

local function Alloc_meta_ctor(elct)
   local elmSz = sizeof(elct) -- element size in bytes
   -- block size in bytes
   local blockSz = 2048*elmSz
   -- function to adjust allocation bytes to block-size: elmSz is the
   -- size of each element, num is number of such elements
   local function calcAdjustedSize(num)
      local adjBytes = (math.floor(elmSz*num/blockSz)+1)*blockSz
      local adjNum = math.floor(adjBytes/elmSz)
      return adjBytes, adjNum
   end

   -- copy function for non-numeric types: this is used in methods
   -- that set array values if the element type stored is not numeric
   local isNumberType = false
   local copyElemFunc = nil
   if ffi.istype(new(elct), new("double")) then
      isNumberType = true
   elseif ffi.istype(new(elct), new("float")) then
      isNumberType = true
   elseif ffi.istype(new(elct), new("int")) then
      isNumberType = true
   elseif ffi.istype(new(elct), new("complex")) then
      isNumberType = true
   else
      isNumberType = false
      copyElemFunc = function (dst, src) ffi.copy(dst, src, elmSz) end
   end   

   -- Use calloc to allocate 'num' elements of memory for type
   -- 'elct'. The type 'ct' is a struct holding size and pointer to
   -- allocated memory. This allocator rounds out the memory to "blockSz"
   -- bytes to prevent fragmentation of the C memory allocator.
   --
   -- NOTE: We use "calloc" and not malloc as calloc zeros out the
   -- memory. Uninitialized memory can cause random crashes in LuaJIT.
   local function alloc(ct, num)
      local adjBytes, adjNum = calcAdjustedSize(num)
      
      local v = new(ct)
      v._capacity = 0
      v._data = calloc(adjNum, elmSz)
      v._capacity = adjNum
      v._size = num
      return v
   end
   
   local alloc_funcs = {
      elemType = function (self)
	 return elct
      end,
      elemSize = function (self)
	 return sizeof(elct)
      end,
      data = function (self)
	 return self._data
      end,
      capacity = function (self)
	 return self._capacity
      end,
      size = function(self)
	 return self._size
      end,
      expand = function (self, nSize)
	 if nSize < self._capacity then
	    self._size = nSize
	 else
	    local adjBytes, adjNum = calcAdjustedSize(nSize)
	    self._data = realloc(self._data, adjBytes)
	    self._capacity = adjNum
	    self._size = nSize
	 end
      end,
      clear = function(self)
	 local adjBytes, adjNum = calcAdjustedSize(0) -- single block
	 self._data = realloc(self._data, adjBytes)
	 self._capacity = adjNum
	 self._size = 0
      end,
      fill = copyElemFunc and
	 function(self, v)
	    for i = 0, self._size-1 do
	       copyElemFunc(self._data[i], v)
	    end
	 end or
	 function (self, v)
	    for i = 0, self._size-1 do
	       self._data[i] = v
	    end
	 end,
      push = copyElemFunc and
	 function (self, v)
	    self:expand(self._size+1)
	    copyElemFunc(self._data[self._size-1], v)
	 end or
	 function (self, v)
	    self:expand(self._size+1)
	    self._data[self._size-1] = v
	 end,
      last = copyElemFunc and
	 function (self)
	    local e = new(elct)
	    copyElemFunc(e, self._data[self._size-1])
	    return e
	 end or
	 function (self)
	    return self._data[self._size-1]
	 end,
      assignLast = copyElemFunc and
	 function (self, v)
	    copyElemFunc(self._data[self._size-1], v)
	 end or
	 function (self, v)
	    self._data[self._size-1] = v
	 end,
      popLast = function(self)
	 assert(self._size > 0, "Can't pop from empty Alloc array")
	 local l = self:last()
	 self._size = self._size-1 -- just reduce size by one
	 return l
      end,
      delete = function (self)
	 if self._capacity == 0 then
	    return false
	 else
	    free(self._data)
	    self._capacity = 0
	 end
	 return true
      end,
   }
   local alloc_mt = {
      __new = function (ct, num)
	 if num then
	    return alloc(ct, num)
	 else
	    return alloc(ct, 0)
	 end
      end,
      __index = function (self, k)
	 -- for non-numeric types, a reference will be returned: this
	 -- is a little dangerous, but it allows for more "natural"
	 -- style of programming when working with C-structs
	 if type(k) == "number" then
	    return self._data[k-1]
	 else
	    return alloc_funcs[k]
	 end
      end,
      __newindex = copyElemFunc and
	 function (self, k, v)
	    copyElemFunc(self._data[k-1], v)
	 end or
	 function (self, k, v)
	    self._data[k-1] = v
	 end,
      __gc = function (self)
	 self:delete()
      end,
   }
   return metatype(typeof("struct { int _size; int _capacity; $* _data; }", elct), alloc_mt)
end

-- function to create an allocator for custom type
local function createAllocator(typeStr)
   return Alloc_meta_ctor(ffi.typeof(typeStr))
end

return {
   Alloc_meta_ctor = Alloc_meta_ctor,
   Double = createAllocator("double"),
   Float = createAllocator("float"),
   Int = createAllocator("int"),
   Complex = createAllocator("complex"),
   malloc = malloc,
   realloc = realloc,
   calloc = calloc,
   free = free,
   createAllocator = createAllocator,
}
