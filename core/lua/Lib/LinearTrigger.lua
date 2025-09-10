-- Gkyl ------------------------------------------------------------------------
--
-- Create function that triggers events when a equally spaced barrier
-- is crossed. "You can enter each barrier only once". Mr. Hyde
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- tStart: start point of trigger
-- tEnd: end point of trigger
-- nBarrier: Number of barriers
return function (tStart, tEnd, nBarrier)
   local tCurr, tEnd = tStart, tEnd
   local dt = (tEnd-tStart)/nBarrier
   local curr = 0

   return function (t)
      local status = false
      if t >= tCurr then
	 status = true
	 curr = curr + 1
	 tCurr = tCurr + dt
      end
      return status, curr
   end
end

-- G2 code
-- return function (tStart, tEnd, nBarrier)
--    local tFrame = (tEnd-tStart)/nBarrier
--    local nextTrigger = tStart
--    local eps = tFrame*1e-9 -- not sure what this should be
--    return function (t)
--       if t>=nextTrigger or math.abs(t-tEnd) < eps then
-- 	 local n = math.floor((t-tStart)/tFrame)+1
-- 	 nextTrigger = n*tFrame
-- 	 return true
--       end
--       return false
--    end
-- end
