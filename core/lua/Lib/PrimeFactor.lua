-- Gkyl ------------------------------------------------------------------------
--
-- Prime-factorization of integers
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local _M = {}

-- The algorithm used here is based on one presented by Knuth and
-- Pardo (1976). It is pretty fast for purposes of Gkeyll but could be
-- sped up somewhat if needed.
-- 
-- Knuth, D., & Pardo, L. T. (1976). Analysis of a Simple
-- Factrorization Algorithm. Theoretical Computer Science, 3, 321–348.

-- Largest prime-factor
function _M.largest(n)
   local m, d, p = n, 2, 1
    while d^2 <= m do
       if m%d == 0 then
	  p, m = d, m/d
       else
	  d = d+1
       end
    end
    return m
end

-- All prime-factors in increasing order
function _M.all(n)
   local m, d, p = n, 2, {}
    while d^2 <= m do
       if m%d == 0 then
	  table.insert(p, d)
	  m = m/d
       else
	  d = d+1
       end
    end
    table.insert(p, m)
    return p
end

return _M
