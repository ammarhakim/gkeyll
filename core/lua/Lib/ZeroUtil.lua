local ffi = require "ffi"

ffi.cdef [[
enum gkyl_elem_type { GKYL_INT, GKYL_INT_64, GKYL_FLOAT, GKYL_DOUBLE, GKYL_USER };

/**
 * Object holding use count and pointer to destructor function.
 */
struct gkyl_ref_count {
  void (*free)(const struct gkyl_ref_count* );
  int count;
};

/**
 * Type of function to project.
 *
 * @param t Time to evaluate function
 * @param xn Coordinates for evaluation
 * @param fout Output vector of 'num_ret_vals'
 * @param ctx Context for function evaluation. Can be NULL
 */
typedef void (*evalf_t)(double t, const double *xn, double *fout, void *ctx);

/**
 * Return .gkyl file type. Returns -1 if file does not exist or is not
 * a gkyl file. file_types = 1 is a field, file_type 2 is a dynvector.
 *
 * @param fname File name
 * @param file type.
 */
int gkyl_get_gkyl_file_type(const char *fname);
]]

local _M = { }

function _M.gkylFileType(fname)
   local file_types = {
      [1] = "field",
      [2] = "dynvector",
      [3] = "multi-range-field",
      [4] = "block-topology",
      [5] = "multi-block-meta"
   }
   local ftype = ffi.C.gkyl_get_gkyl_file_type(fname)
   if ftype < 0 then return "not-gkyl" end
   return file_types[ftype]
end

return _M
