-- Gkyl ------------------------------------------------------------------------
--
-- Function to generate a unique ID everytime it is called
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Time = require "Lib.Time"

local bitsize = 32  -- bitsize assumed for Lua VM. See randomseed function below.
local lua_version = tonumber(_VERSION:match("%d%.*%d*"))  -- grab Lua version used

-- Improved randomseed function (See
-- https://github.com/Tieske/uuid/blob/master/src/uuid.lua)
--
-- Lua 5.1 and 5.2 both truncate the seed given if it exceeds the
-- integer range. If this happens, the seed will be 0 or 1 and all
-- randomness will be gone (each application run will generate the
-- same sequence of random numbers in that case). This improved
-- version drops the most significant bits in those cases to get the
-- seed within the proper range again.
--
-- @param seed the random seed to set (integer from 0 - 2^32, negative values will be made positive)
-- @return the (potentially modified) seed used

local function randomseed(seed)
  seed = math.floor(math.abs(seed))
  if seed >= (2^bitsize) then
    -- integer overflow, so reduce to prevent a bad seed
    seed = seed - math.floor(seed / 2^bitsize) * (2^bitsize)
  end
  if lua_version < 5.2 then
    -- 5.1 uses (incorrect) signed int
    math.randomseed(seed - 2^(bitsize-1))
  else
    -- 5.2 uses (correct) unsigned int
    math.randomseed(seed)
  end
  return seed
end

-- seed the generator
randomseed( Time.time()*10000 )

-- UUID generator
return function ()
   local random = math.random
   local template ='xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'
   return string.gsub(
      template, '[xy]', function (c)
	 local v = (c == 'x') and random(0, 0xf) or random(8, 0xb)
	 return string.format('%x', v)
   end)
end
