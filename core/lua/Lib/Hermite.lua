-- Gkyl ------------------------------------------------------------------------
--
-- Normalized Hermite functions and their derivatives:
--
-- h(m,x) = exp(-x^2)*H(m,x)/sqrt(2^m*pi*m!)
--
-- where H(m,x) are Physicist's Hermite polynomials.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local Lin = require "Lib.Linalg"

_M = {}

function _M.H(m, x)
   -- Normalized Hermite function

   local hm2 = 1/math.sqrt(math.pi)*math.exp(-x^2)
   local hm1 = 1/math.sqrt(2*math.pi)*2*x*math.exp(-x^2)

   if m == 0 then return hm2 end
   if m == 1 then return hm1 end
   
   local hm = 0
   for n = 2, m do
      hm = math.sqrt(2/n)*x*hm1 - math.sqrt((n-1)/n)*hm2
      hm2, hm1 = hm1, hm
   end
   return hm
end

function _M.dH(m, x)
   -- Derivative of normalized Hermite function
   return -math.sqrt(2*(m+1))*_M.H(m,x)
end

function _M.Hallm(mmax, x)
   -- Normalized Hermite functions at a given x for all m<=mmax.
   --
   -- NOTE: Returned array is indexed starting 0

   local _data = ffi.new("double[?]", mmax+1)
   local hm = Lin.Vec(mmax+1, _data+1)

   if mmax > 0 then hm[0] = 1/math.sqrt(math.pi)*math.exp(-x^2) end
   if mmax > 1 then hm[1] = 1/math.sqrt(2*math.pi)*2*x*math.exp(-x^2) end

   for n = 2, mmax do
      hm[n] = math.sqrt(2/n)*x*hm[n-1] - math.sqrt((n-1)/n)*hm[n-2]
   end
   return hm
end

function _M.dHallm(mmax, x)
   -- Derivative of normalized Hermite function at a given x for all
   -- m<=mmax.
   -- 
   -- NOTE: Returned array is indexed starting 0

   local hm = _M.Hallm(mmax, x)
   for m = 0, mmax do
      hm[m] = -math.sqrt(2*(m+1))*hm[m]
   end
   return hm
end

return _M
