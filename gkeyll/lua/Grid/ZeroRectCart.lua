-- Gkyl ------------------------------------------------------------------------
--
-- Lua interface to gkylzero's gkyl_dynvec structure
---
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local ffi = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"
local new, sizeof, typeof, metatype = xsys.from(ffi,
     "new, sizeof, typeof, metatype")

local Alloc = require "Lib.Alloc"
local Lin = require "Lib.Linalg"
local Range = require "Lib.Range"

local cuda
require "Lib.ZeroUtil"

_M = {}

ffi.cdef [[

/**
 * Rectangular grid object.
 */
struct gkyl_rect_grid {
  int _ndim; // number of dimensions
  double _lower[7]; // lower-left corner
  double _upper[7]; // upper-right corner
  int _cells[7]; // number of cells    
  double _dx[7]; // cell spacing
  double _cellVolume; // cell volume
};

/**
 * Create new grid object.
 *
 * @param grid Grid object to initialize.
 * @param ndim Dimension of grid
 * @param lower Coordinates of lower-left corner of grid
 * @param upper Coordinates of upper-right corner of grid
 * @param cells Number of cells in each direction
 */
void gkyl_rect_grid_init(struct gkyl_rect_grid *grid, int ndim,
  const double *lower, const double *upper, const int *cells);

/**
 * Create range and extended ranges from grid and ghost-cell data. The
 * range is a sub-range of the extended range.
 *
 * @param grid Grid to compute ranges for
 * @param nghost Number of ghost-cells in each direction
 * @param ext_range On output, extended range spanning grid+ghost-cells
 * @param range On output, range spanning grid. Sub-range of ext_range.
 */
void gkyl_create_grid_ranges(const struct gkyl_rect_grid *grid,
  const int *nghost, struct gkyl_range *ext_range, struct gkyl_range *range);

/**
 * Compare grids
 *
 * @param grid1 Grid object to compare
 * @param grid2 Grid object to compare
 * @return true if the grids are the same, false otherwise
 */
bool gkyl_rect_grid_cmp(const struct gkyl_rect_grid *grid1, struct gkyl_rect_grid *grid2);
]]


-- Rect Grid ctype
local RectCartCt = typeof("struct gkyl_rect_grid")

local rect_cart_fn = {
   compare = function(self, grid)
      return ffi.C.gkyl_rect_grid_cmp(self, grid)
   end,
   -- returns range, ext_range
   createGridRanges = function(self, nghost)
      local ndim = self._ndim
      local vnghost = Lin.IntVec(ndim)
      for i = 1, ndim do vnghost[i] = nghost[i] end
      
      local range = ffi.new("struct gkyl_range")
      local ext_range = ffi.new("struct gkyl_range")
      ffi.C.gkyl_create_grid_ranges(self, vnghost:data(), ext_range, range)
      return range, ext_range
   end,
}

local rect_cart_mt = {
   __new = function (self, ndim, lower, upper, cells)
      local vlower = Lin.Vec(ndim)
      local vupper = Lin.Vec(ndim)
      local vcells = Lin.IntVec(ndim)
      local grid = ffi.new("struct gkyl_rect_grid")
      ffiC.gkyl_rect_grid_init(grid, ndim, vlower:data(), vupper:data(), vells:data())
      return grid
   end,
   __index = rect_cart_fn
}
local Rect_CartCtor = metatype(RectCartCt, rect_cart_mt)

-- Compare two grids
function _M.compare(grid1, grid2)
   return ffi.C.gkyl_rect_grid_cmp(grid1, grid2)
end

return _M
