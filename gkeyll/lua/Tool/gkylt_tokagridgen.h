#pragma once

#include <stdbool.h>
#include <gkyl_gk_block_geom.h>

// Forward declarations
struct gkyl_gk_block_geom;
struct gkyl_gk_block_geom_info_lite; 
// Moved from gkyl_toka_grid_gen.h to reduce headers
enum gkyl_toka_grid_gen_null_points {
    GKYL_TOKA_GRID_GEN_SINGLE_NULL, // Single Null  Configuration
    GKYL_TOKA_GRID_GEN_DOUBLE_NULL, // Double Null Configuration
};


// Function pointer type for plate functions
typedef void (*plate_func_t)(double s, double* RZ);

struct gkylt_tokagridgen_inp {
    enum gkyl_toka_grid_gen_null_points toka_type; //Tokamak configuration
  
    double wpsi[4]; // Width of psi for different regions
    int cell_count[6]; // number of cells in computational space
        
    // Plate function pointers - these will be set from Lua input
    plate_func_t plate_func_lower_outer;
    plate_func_t plate_func_upper_outer;
    plate_func_t plate_func_lower_inner;
    plate_func_t plate_func_upper_inner;
    
    // Region-specific parameters (conditional based on toka_type)
    // For double_null: lower_PF[6], outer_SOL[7], inner_SOL[6], core[5]
    // For single_null: lower_PF[8], SOL[7], core[8]
    double lower_PF[8]; // lower private flux region parameters (max size for both types)
    double outer_SOL[7]; // outer scrape-off layer parameters (DN only)
    double inner_SOL[6]; // inner scrape-off layer parameters (DN only)
    double SOL[7]; // scrape-off layer parameters (SN only)
    double core[7]; // core region parameters (max size for both types)
    
    const char *geqdsk_path; // name for file with psi(R,Z)
    const char *out_prefix; // output prefix
    const char *name; // application name
    bool use_half_domain; // flag to indicate if half domain should be used (DN only)
  };
  
struct gkyl_gk_block_geom* gkylt_tokagridgen(const struct gkylt_tokagridgen_inp *tginp);

// Lite accessor to export safe block geometry info to Lua
void gkyl_gk_block_geom_get_block_lite(const struct gkyl_gk_block_geom *bgeom,
  int bidx, struct gkyl_gk_block_geom_info_lite *out);

