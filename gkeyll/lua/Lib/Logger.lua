-- Gkyl ------------------------------------------------------------------------
--
-- Simple logger for use in sims
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Only load Mpi when the executable was built with MPI support.
-- Without MPI there is always exactly one process (rank 0), so we
-- provide a minimal stub that satisfies the same interface.
local Mpi
if GKYL_HAVE_MPI then
   Mpi = require "Comm.Mpi"
else
   Mpi = {
      COMM_WORLD  = nil,
      Comm_rank   = function(comm) return 0 end,
   }
end

local ffi  = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- write message to list of output streams
local function writeToFile(outStreams, msg)
   for _, s in ipairs(outStreams) do
      s:write(msg); s:flush()
   end
end

-- create logger
local function makeLogger(tbl)
   local writeRank = 0
   if tbl.writeRank then writeRank = tbl.writeRank end

   local fName
   if tbl.logFile then
      fName = GKYL_OUT_PREFIX .. "_" .. tbl.logFile .. ".log"
   else
      fName = GKYL_OUT_PREFIX .. "_" .. writeRank .. ".log"
   end

   local comm = tbl.comm and tbl.comm or Mpi.COMM_WORLD
   local rank = Mpi.Comm_rank(comm)

   local outStreams = {}
   if writeRank == rank then
      outStreams[1] = io.stdout
      if tbl.logToFile then
	 outStreams[2] = io.open(fName, "w")
      end
   end

   return function (msg)
      writeToFile(outStreams, msg)
   end
end

return makeLogger
