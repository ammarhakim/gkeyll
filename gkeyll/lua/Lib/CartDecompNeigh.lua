-- Gkyl ------------------------------------------------------------------------
--
-- Compute and store neighbor information for communication in parallel.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local Range = require "Lib.Range"
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"

local CartDecompNeigh = Proto()

-- Constructor to make new neighbor object.
function CartDecompNeigh:init(decomp)
   self._decomp = decomp -- Decomposition object.
   self._neighData = {}     -- Neighbor data.
end

-- Set callable methods.
function CartDecompNeigh:calcFaceCommNeigh(lowerGhost, upperGhost)
   self._neighData = {} -- clear out existing data
   local ndim = self._decomp:ndim()
   local numSubDomains = self._decomp:numSubDomains()
   for kme = 1, numSubDomains do
      local nlst = {} -- List of neighbors
      for d = 1, ndim  do
	 -- Expand sub-domain in direction `d`.
	 local expSubDom = self._decomp:subDomain(kme):extendDir(d, lowerGhost, upperGhost)
	 -- Loop over all other sub-domains and intersect.
	 for ku = 1, numSubDomains do
	    if ku == kme then goto continue end -- No self-intersections.
	    if not expSubDom:isIntersectionEmpty(self._decomp:subDomain(ku)) then
	       table.insert(nlst, ku)   -- Insert subDomain index into list of neighbors.
	    end
	    ::continue::
	 end
      end
      self._neighData[kme] = nlst
   end
end

function CartDecompNeigh:calcAllCommNeigh(lowerGhost, upperGhost)
   self._neighData = {} -- Clear out existing data.
   local ndim = self._decomp:ndim()
   local numSubDomains = self._decomp:numSubDomains()
   for kme = 1, numSubDomains do
      local nlst = {} -- List of neighbors
      -- Expand sub-domain.
      local expSubDom = self._decomp:subDomain(kme):extend(lowerGhost, upperGhost)
      -- Loop over all other sub-domains and intersect.
      for ku = 1, numSubDomains do
	 if ku == kme then goto continue end -- No self-intersections.
	 if not expSubDom:isIntersectionEmpty(self._decomp:subDomain(ku)) then
	    table.insert(nlst, ku) -- Insert subDomain index into list of neighbors.
	 end
	 ::continue::
      end
      self._neighData[kme] = nlst
   end
end

function CartDecompNeigh:neighborData(k)
   return self._neighData[k] and self._neighData[k] or {}
end

return CartDecompNeigh
