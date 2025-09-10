-- Automatic differentiation with recursive dual numbers.
-- This extends the functionality of sci/diff.lua to be able to calculate arbitrary order derivatives.
-- Given a function f(x1,x2,...), can compute derivative with respect to xi with
-- dn.deriv(f,i)(x1,x2,...)
-- Can compute higher order derivatives by calling deriv multiple times, e.g.
-- dn.deriv(dn.deriv(dn.deriv(f,i),i),j)(x1,x2,...)

local math = require "sci.math".generic
local xsys = require "xsys"

local
abs, acos, asin, atan, atan2, ceil, cos, cosh, deg, exp, floor, fmod, frexp,
huge, ldexp, log, log10, max, min, modf, pi, pow, rad, random, randomseed, sin, 
sinh, sqrt, tan, tanh,
round, step, sign,
phi, iphi, gamma, loggamma, logbeta, beta
= xsys.from(math, [[
abs, acos, asin, atan, atan2, ceil, cos, cosh, deg, exp, floor, fmod, frexp,
huge, ldexp, log, log10, max, min, modf, pi, pow, rad, random, randomseed, sin, 
sinh, sqrt, tan, tanh,
round, step, sign,
phi, iphi, gamma, loggamma, logbeta, beta
]])

local dn = {}
dn.__index = dn

setmetatable(dn, {
  __call = function(self, ...)
     return self.new(...)
  end,
})

function dn.new(v, a)
  local self = setmetatable( {}, dn)
  self._v = v
  self._a = a or 0
  return self
end
dn.__unm = function(x)
  return dn(-x._v, -x._a)
end
dn.__add = function(x, y) 
  if type(x) == "number" then x = dn(x) end
  if type(y) == "number" then y = dn(y) end
  return dn(x._v + y._v, x._a + y._a)
end
dn.__sub = function(x, y) 
  if type(x) == "number" then x = dn(x) end
  if type(y) == "number" then y = dn(y) end
  return dn(x._v - y._v, x._a - y._a)
end
dn.__mul = function(x, y)
  if type(x) == "number" then x = dn(x) end
  if type(y) == "number" then y = dn(y) end
  return dn(x._v*y._v, x._a*y._v + y._a*x._v)
end
dn.__div = function(x, y)
  if type(x) == "number" then x = dn(x) end
  if type(y) == "number" then y = dn(y) end
  return dn(x._v/y._v, (x._a*y._v - y._a*x._v)/y._v^2)
end
dn.__pow = function(x, y) -- Optimized version.
  if type(y) == "number" then
    return dn(x._v^y, y*x._v^(y-1)*x._a)
  elseif type(x) == "number" then
    return dn(x^y._v, x^y._v*log(x)*y._a)
  else
    return dn(x._v^y._v, x._v^y._v*(log(x._v)*y._a + y._v/x._v*x._a))
  end
end
-- unable to correctly overload comparison operations __eq, __lt, and __le
-- so we need to define actual functions for these operations
dn.eq = function(x, y) 
  if type(x) == "number" and type(y) == "number" then return x == y end
  if type(x) == "number" then x = dn(x) end
  if type(y) == "number" then y = dn(y) end
  return dn.eq(x._v, y._v)
end
dn.lt = function(x, y)
  if type(x) == "number" and type(y) == "number" then return x < y end
  if type(x) == "number" then x = dn(x) end
  if type(y) == "number" then y = dn(y) end
  return dn.lt(x._v, y._v)
end
dn.le = function(x, y)
  if type(x) == "number" and type(y) == "number" then return x <= y end
  if type(x) == "number" then x = dn(x) end
  if type(y) == "number" then y = dn(y) end
  return dn.le(x._v, y._v)
end
dn.__tostring = function(x)
  return tostring(x._v) -- Better to mimic behavior of numbers.
end

dn.copy = function(x)
  return dn(x)
end
dn.sin  = function(x) return dn(sin(x._v), x._a*cos(x._v)) end
dn.cos  = function(x) return dn(cos(x._v), -x._a*sin(x._v)) end
dn.tan  = function(x) return dn(tan(x._v),  x._a*(1 + tan(x._v)^2)) end
dn.asin = function(x) return dn(asin(x._v), x._a/sqrt(1 - x._v^2)) end
dn.acos = function(x) return dn(acos(x._v), -x._a/sqrt(1 - x._v^2)) end
dn.atan = function(x) return dn(atan(x._v), x._a/(1 + x._v^2)) end
dn.sinh = function(x) return dn(sinh(x._v), x._a*cosh(x._v)) end
dn.cosh = function(x) return dn(cosh(x._v), x._a*sinh(x._v)) end
dn.tanh = function(x) return dn(tanh(x._v), x._a*(1 - tanh(x._v)^2)) end
dn.exp  = function(x) return dn(exp(x._v),  x._a*exp(x._v)) end
dn.log  = function(x) return dn(log(x._v),  x._a/x._v) end
dn.sqrt = function(x) return dn(sqrt(x._v), x._a/(2*sqrt(x._v))) end
dn.abs  = function(x) return dn(abs(x._v),  x._a*sign(x._v)) end
-- Stick to dn type to improve type stability:
dn.floor = function(x) return dn(floor(x._v), 0) end
dn.ceil  = function(x) return dn(ceil(x._v),  0) end
-- Stick to dn type to improve type stability:
dn.round = function(x) return dn(round(x._v), 0) end
dn.step  = function(x) return dn(step(x._v),  0) end
dn.sign  = function(x) return dn(sign(x._v),  0) end

-- Function to compute a derivative of an arbitrary function.
-- Inputs:
--  f : a function f(x1, x2, ...) that returns one or more value(s)
--  di: the index of the argument to differentiate f with respect to
dn.deriv = function(f,di)
  return function(...)
    xn = {...}
    for i=1, #xn do
      -- determine which arg to differentiate w.r.t.
      if i == (di or 1) then 
        -- for the arg that is differentiated, substitute dn(x,1) for x
        xn[i] = dn(xn[i],1) 
      else
        -- for other xn, substitute dn(x,0) for x
        xn[i] = dn(xn[i],0)
      end
    end
    retvals = {f(unpack(xn))}
    for i=1, #retvals do
      if type(retvals[i])~="number" then
        retvals[i] = retvals[i]._a
      else
        retvals[i] = 0
      end
    end
    return unpack(retvals)
  end
end

-- same as deriv, but arguments of f are expected as a table
dn.derivt = function(f,di)
  return function(xn)
    local dxn = {}
    for i=1, #xn do
      -- determine which arg to differentiate w.r.t.
      if i == (di or 1) then 
        -- for the arg that is differentiated, substitute dn(x,1) for x
        dxn[i] = dn(xn[i],1) 
      else
        -- for other args, substitute dn(x,0) for x
        dxn[i] = dn(xn[i],0)
      end
    end
    retvals = {f(dxn)}
    for i=1, #retvals do
      if type(retvals[i])~="number" then
        retvals[i] = retvals[i]._a
      else
        retvals[i] = 0
      end
    end
    return unpack(retvals)
  end
end

local df = dn.deriv

return {
  dn = dn,
  df = dn.deriv,
  deriv = dn.deriv,
  derivt = dn.derivt,
  eq = dn.eq,
  lt = dn.lt,
  le = dn.le,
}
