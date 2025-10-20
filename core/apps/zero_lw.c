#ifdef GKYL_HAVE_LUA

#include <gkyl_alloc.h>
#include <gkyl_lua_utils.h>
#include <gkyl_lw_priv.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_zero_lw.h>

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>

#include <string.h>

struct gkyl_tool_args *
gkyl_tool_args_new(lua_State *L)
{
  struct gkyl_tool_args *args = gkyl_malloc(sizeof(*args));

  with_lua_global(L, "GKYL_COMMANDS") {
    args->argc = glua_objlen(L);
    args->argv = 0;

    if (args->argc > 0) {
      args->argv = gkyl_malloc(args->argc*sizeof(char*));
    
      for (int i = 1; i <= glua_objlen(L); ++i)  {
        const char *av = glua_tbl_iget_string(L, i, "x");
        args->argv[i-1] = gkyl_malloc(1+strlen(av));
        strcpy(args->argv[i-1], av);
      }
    }
  }

  return args;
}

void
gkyl_tool_args_release(struct gkyl_tool_args* args)
{
  for (int i=0; i<args->argc; ++i)
    gkyl_free(args->argv[i]);
  gkyl_free(args->argv);
  gkyl_free(args);
}

/* *********************/
/* Rect decomp methods */
/* *********************/

// Metatable name for top-level Vlasov App
#define RECT_DECOMP_METATABLE_NM "GkeyllZero.Zero.RectDecomp"

// Lua userdata object for holding Vlasov app and run parameters
struct rect_decomp_lw {
  struct gkyl_rect_decomp *decomp; // Decomp object
};

// G0.RectDecomp.new { cells = { 100, 100}, cuts = { 2, 2 } }
static int
rect_decomp_lw_new(lua_State *L)
{
  struct rect_decomp_lw *rd_lw = gkyl_malloc(sizeof(*rd_lw));
  
  int ndim = 1;
  int cells[GKYL_MAX_DIM] = { 0 }, cuts[GKYL_MAX_DIM] = { 1 };
  
  with_lua_tbl_tbl(L, "cells") {
    ndim = glua_objlen(L);
    for (int d=0; d<ndim; ++d)
      cells[d] = glua_tbl_iget_integer(L, d+1, 0);
  }

  with_lua_tbl_tbl(L, "cuts") {
    for (int d=0; d<ndim; ++d)
      cuts[d] = glua_tbl_iget_integer(L, d+1, 0);
  }

  rd_lw->decomp = gkyl_rect_decomp_new_from_cuts_and_cells(ndim, cuts, cells);
  
  // create Lua userdata ...
  struct rect_decomp_lw **l_rd_lw = lua_newuserdata(L, sizeof(struct rect_decomp_lw*));
  *l_rd_lw = rd_lw; // ... point it to the rect decomp pointer

  // set metatable
  luaL_getmetatable(L, RECT_DECOMP_METATABLE_NM);
  lua_setmetatable(L, -2);  
  
  return 1;
}

// Clean up memory allocated for decomp
static int
rect_decomp_lw_gc(lua_State *L)
{
  struct rect_decomp_lw **l_rd_lw = GKYL_CHECK_UDATA(L, RECT_DECOMP_METATABLE_NM);
  struct rect_decomp_lw *rd_lw = *l_rd_lw;

  gkyl_rect_decomp_release(rd_lw->decomp);
  gkyl_free(*l_rd_lw);
  
  return 0;
}

// rect_decomp constructor
static struct luaL_Reg rect_decomp_ctor[] = {
  { "new",  rect_decomp_lw_new },
  { 0, 0 }
};

// rect_decomp methods
static struct luaL_Reg rect_decomp_funcs[] = {
  {0, 0}
};

static void
rect_decomp_openlibs(lua_State *L)
{
  do {
    luaL_newmetatable(L, RECT_DECOMP_METATABLE_NM);

    lua_pushstring(L, "__gc");
    lua_pushcfunction(L, rect_decomp_lw_gc);
    lua_settable(L, -3);

    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    luaL_register(L, NULL, rect_decomp_funcs);
    
    luaL_register(L, "G0.Zero.RectDecomp", rect_decomp_ctor);
    
  } while (0);
}

void
gkyl_zero_lw_openlibs(lua_State *L)
{
  // Push empty global table called "G0".
  lua_newtable(L);
  lua_setglobal(L, "G0");
  
  rect_decomp_openlibs(L);

  // Register types for species, field, and Poisson solver boundary conditions.
  gkyl_register_distribution_moment_types(L);
  gkyl_register_species_bc_types(L);
  gkyl_register_field_bc_types(L);
}

#endif
