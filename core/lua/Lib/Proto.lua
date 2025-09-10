-- Gkyl ------------------------------------------------------------------------
--
-- Base object for prototype based object-system
-- 
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- call objects's init() (ctor) method if one is provided
local function newmember(proto, ...)
   local obj = setmetatable({}, proto)

   if obj.__gc then
      local prox = newproxy(true)
      getmetatable(prox).__gc = function() obj:__gc() end
      obj[prox] = true
   end

   if GKYL_USE_DEVICE then
      -- when we are using CUDA we may have a initDevice method to
      -- construct things on a device
      if obj.init then
	 if obj.initDevice then
	    return obj, obj:init(...), obj:initDevice(...)
	 else
	    return obj, obj:init(...)
	 end
      end
   else
      if obj.init then return obj, obj:init(...) end
   end
   return obj
end
-- metatable for Proto: redirect call to new()
local protometa = {  __call = function(self, ...) return self:new(...) end }

-- check if object is contructed from a Proto 'pr'
local function isa(self, pr)
   assert(pr, "isa expected a Proto object")
   if self.proto == pr then return true end
   if self.proto.super then return isa(self.proto.super, pr) end
   return false
end

-- ProtoTable-------------------------------------------------------------------
--
-- Provides object to store memembers/methods
--------------------------------------------------------------------------------

local ProtoTable = {}
ProtoTable.__index = ProtoTable -- so functions can be called
-- append all (key,value) pairs from each table passed to new() into
-- ProtoTable
function ProtoTable.new(...)
   local t = setmetatable({}, ProtoTable)
   for _, o in ipairs{...} do
      for k, v in pairs(o) do t[k] = v end
   end
   return t
end
-- make ProtoTable object callable, redirecting call to new()
setmetatable(ProtoTable, { __call = function(t, ...) return ProtoTable.new(...) end })

-- Proto -----------------------------------------------------------------------
--
-- Provides a primitive prototype object which can be used to define
-- objects in Gkyl
--------------------------------------------------------------------------------
local function Proto(...)
   local pr = ProtoTable(...) -- append all base vars/methods into this one
   pr.proto = pr

   -- set base object pointers
   local parents = {...}
   pr.super = parents[1]
   pr.supers = parents
   
   pr.__index = pr -- redirect calls to those defined by objects
   pr.new = newmember -- redirect new() to newmember
   pr.isa = isa -- for type checking
   
   -- won't fail unlike isa() method if supplied type is not a table
   pr.is = function(x)
      return type(x) == 'table' and x.isa and x:isa(pr)
   end

   setmetatable(pr, protometa)
   return pr
end

return Proto
