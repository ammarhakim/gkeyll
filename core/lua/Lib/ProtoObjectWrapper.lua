-- Gkyl ------------------------------------------------------------------------
--
-- Create a wrapper for Proto objects to allow appending another table
-- before the init() method is called.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local xsys = require "xsys"

return function (obj)
   return function (baseTbl)
      return function (auxTbl)
	 return obj( xsys.table.union(baseTbl, auxTbl) )
      end
   end
end
