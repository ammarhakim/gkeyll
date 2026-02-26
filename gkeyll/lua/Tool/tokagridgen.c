#include <gkyl_array_rio.h>
#include <gkylt_tokagridgen.h>


#include <stc/cstr.h>
#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <gkyl_mpi_comm.h>
#include <gkyl_null_comm.h>

// Lightweight POD struct to expose only required fields to Lua FFI
struct gkyl_gk_block_geom_info_lite {
  double lower[3], upper[3];
  int cells[3], cuts[3];

  int geometry_id;

  struct {
    char filepath[1024];
    int rz_poly_order;
    int flux_poly_order;
    int reflect; // bool as int for FFI stability
    int use_gpu; // bool as int for FFI stability
  } efit_info;

  struct {
    int ftype;
    int half_domain; // bool as int
    double rclose, rleft, rright, rmin, rmax;
    double zmin, zmax, zmin_left, zmin_right, zmax_left, zmax_right;
    int plate_spec; // bool as int
    int inexact_roots; // bool as int
    int use_cubics; // bool as int
    int use_hyperbolic_numbers; // bool as int
  } tok_grid_info;

  struct { int bid; int dir; int edge; } connections[3][2];
};

// Fills a lite struct for block bidx from the heavy structs
void gkyl_gk_block_geom_get_block_lite(const struct gkyl_gk_block_geom *bgeom,
  int bidx, struct gkyl_gk_block_geom_info_lite *out)
{
  const struct gkyl_gk_block_geom_info *info = gkyl_gk_block_geom_get_block(bgeom, bidx);

  // Copy basic arrays
  for (int d = 0; d < 3; ++d) {
    out->lower[d] = info->lower[d];
    out->upper[d] = info->upper[d];
    out->cells[d] = info->cells[d];
    out->cuts[d] = info->cuts[d];
  }

  // Geometry id
  out->geometry_id = info->geometry.geometry_id;

  // efit_info
  memset(out->efit_info.filepath, 0, sizeof(out->efit_info.filepath));
  strncpy(out->efit_info.filepath, info->geometry.efit_info.filepath,
    sizeof(out->efit_info.filepath)-1);
  out->efit_info.rz_poly_order = info->geometry.efit_info.rz_poly_order;
  out->efit_info.flux_poly_order = info->geometry.efit_info.flux_poly_order;
  out->efit_info.reflect = info->geometry.efit_info.reflect ? 1 : 0;
  out->efit_info.use_gpu = info->geometry.efit_info.use_gpu ? 1 : 0;

  // tok_grid_info
  out->tok_grid_info.ftype = info->geometry.tok_grid_info.ftype;
  out->tok_grid_info.half_domain = info->geometry.tok_grid_info.half_domain ? 1 : 0;
  out->tok_grid_info.rclose = info->geometry.tok_grid_info.rclose;
  out->tok_grid_info.rleft = info->geometry.tok_grid_info.rleft;
  out->tok_grid_info.rright = info->geometry.tok_grid_info.rright;
  out->tok_grid_info.rmin = info->geometry.tok_grid_info.rmin;
  out->tok_grid_info.rmax = info->geometry.tok_grid_info.rmax;
  out->tok_grid_info.zmin = info->geometry.tok_grid_info.zmin;
  out->tok_grid_info.zmax = info->geometry.tok_grid_info.zmax;
  out->tok_grid_info.zmin_left = info->geometry.tok_grid_info.zmin_left;
  out->tok_grid_info.zmin_right = info->geometry.tok_grid_info.zmin_right;
  out->tok_grid_info.zmax_left = info->geometry.tok_grid_info.zmax_left;
  out->tok_grid_info.zmax_right = info->geometry.tok_grid_info.zmax_right;
  out->tok_grid_info.plate_spec = info->geometry.tok_grid_info.plate_spec ? 1 : 0;
  out->tok_grid_info.inexact_roots = info->geometry.tok_grid_info.inexact_roots ? 1 : 0;
  out->tok_grid_info.use_cubics = info->geometry.tok_grid_info.use_cubics ? 1 : 0;
  out->tok_grid_info.use_hyperbolic_numbers = info->geometry.tok_grid_info.use_hyperbolic_numbers ? 1 : 0;

  // connections
  for (int d = 0; d < 3; ++d) {
    for (int e = 0; e < 2; ++e) {
      out->connections[d][e].bid = info->connections[d][e].bid;
      out->connections[d][e].dir = info->connections[d][e].dir;
      out->connections[d][e].edge = info->connections[d][e].edge;
    }
  }
}

struct gkylt_tokagridgen_output gkylt_tokagridgen(const struct gkylt_tokagridgen_inp *inp)
{
    // Plate functions are now passed as function pointers from Lua input

    const char *name = inp->name;
    int cell_count[6];
    double wpsi[4];
    for (int d=0; d<6; ++d) {
        cell_count[d] = inp->cell_count[d];
    }
    for (int d=0; d<4; ++d) {
        wpsi[d] = inp->wpsi[d];
    }
     bool reflect = (inp->toka_type == GKYL_TOKA_GRID_GEN_DOUBLE_NULL);
     bool use_half_domain = inp->use_half_domain;
     struct gkyl_efit_inp efit_inp = {
        .rz_poly_order = 2,
        .flux_poly_order = 1,
        .reflect = reflect,
        .use_gpu = false,
      };
    strncpy(efit_inp.filepath, inp->geqdsk_path, sizeof(efit_inp.filepath)-1);
    efit_inp.filepath[sizeof(efit_inp.filepath)-1] = '\0';   
    struct gkyl_efit *efit = gkyl_efit_new(&efit_inp);
    double psisep = efit->psisep;
    double psicenter = efit ->simag;    
    char psi_filename[256];
    snprintf(psi_filename, sizeof(psi_filename), "%s_psi.gkyl", name);
    gkyl_grid_sub_array_write(&efit->rzgrid, &efit->rzlocal, 0, efit->psizr, psi_filename);
    gkyl_efit_release(efit);

    enum gkyl_toka_grid_gen_null_points toka_type = inp->toka_type;
    
    // Declare individual variables for both configurations
    // Double null variables
    double rleft_PF, rright_PF, rmin_PF, rmax_PF, zmin_left_PF, zmin_right_PF;
    double rclose_outer_SOL, rleft_outer_SOL, rright_outer_SOL, rmin_outer_SOL, rmax_outer_SOL, zmin_outer_SOL, zmax_outer_SOL;
    double rleft_inner_SOL, rright_inner_SOL, rmin_inner_SOL, rmax_inner_SOL, zmin_inner_SOL, zmax_inner_SOL;
    double rclose_outer_core, rclose_inner_core, rleft_core, rright_core, rmin_core, rmax_core;
    
    // Single null variables
    double zmin_PF, zmax_PF;
    double rclose_SOL, rleft_SOL, rright_SOL, rmin_SOL, rmax_SOL, zmin_SOL, zmax_SOL;
    double rclose_core, zmin_core, zmax_core;
    
    // Initialize all variables to zero
    rleft_PF = rright_PF = rmin_PF = rmax_PF = zmin_left_PF = zmin_right_PF = 0.0;
    rclose_outer_SOL = rleft_outer_SOL = rright_outer_SOL = rmin_outer_SOL = rmax_outer_SOL = zmin_outer_SOL = zmax_outer_SOL = 0.0;
    rleft_inner_SOL = rright_inner_SOL = rmin_inner_SOL = rmax_inner_SOL = zmin_inner_SOL = zmax_inner_SOL = 0.0;
    rclose_core = rleft_core = rright_core = rmin_core = rmax_core = 0.0;
    zmin_PF = zmax_PF = 0.0;
    rclose_SOL = rleft_SOL = rright_SOL = rmin_SOL = rmax_SOL = zmin_SOL = zmax_SOL = 0.0;
    rclose_outer_core = rclose_inner_core = zmin_core = zmax_core = 0.0;
    
    if (toka_type == GKYL_TOKA_GRID_GEN_DOUBLE_NULL) {
        // Double null configuration: lower_PF[6], outer_SOL[7], inner_SOL[6], core[5]
        // Load values from inp struct into individual variables
        rleft_PF = inp->lower_PF[0];
        rright_PF = inp->lower_PF[1];
        rmin_PF = inp->lower_PF[2];
        rmax_PF = inp->lower_PF[3];
        zmin_left_PF = inp->lower_PF[4];
        zmin_right_PF = inp->lower_PF[5];
        
        rclose_outer_SOL = inp->outer_SOL[0];
        rleft_outer_SOL = inp->outer_SOL[1];
        rright_outer_SOL = inp->outer_SOL[2];
        rmin_outer_SOL = inp->outer_SOL[3];
        rmax_outer_SOL = inp->outer_SOL[4];
        zmin_outer_SOL = inp->outer_SOL[5];
        zmax_outer_SOL = inp->outer_SOL[6];
        
        rleft_inner_SOL = inp->inner_SOL[0];
        rright_inner_SOL = inp->inner_SOL[1];
        rmin_inner_SOL = inp->inner_SOL[2];
        rmax_inner_SOL = inp->inner_SOL[3];
        zmin_inner_SOL = inp->inner_SOL[4];
        zmax_inner_SOL = inp->inner_SOL[5];
        
        rclose_inner_core = inp->core[0];
        rclose_outer_core = inp ->core[1];
        rleft_core = inp->core[2];
        rright_core = inp->core[3];
        rmin_core = inp->core[4];
        rmax_core = inp->core[5];
    } else if (toka_type == GKYL_TOKA_GRID_GEN_SINGLE_NULL) {
        // Single null configuration: lower_PF[8], SOL[7], core[8]
        // Load values from inp struct into individual variables
        rleft_PF = inp->lower_PF[0];
        rright_PF = inp->lower_PF[1];
        rmin_PF = inp->lower_PF[2];
        rmax_PF = inp->lower_PF[3];
        zmin_PF = inp->lower_PF[4];
        zmax_PF = inp->lower_PF[5];
        zmin_left_PF = inp->lower_PF[6];
        zmin_right_PF = inp->lower_PF[7];
        
        rclose_SOL = inp->SOL[0];
        rleft_SOL = inp->SOL[1];
        rright_SOL = inp->SOL[2];
        rmin_SOL = inp->SOL[3];
        rmax_SOL = inp->SOL[4];
        zmin_SOL = inp->SOL[5];
        zmax_SOL = inp->SOL[6];
        
        rclose_core = inp->core[0];
        rleft_core = inp->core[1];
        rright_core = inp->core[2];
        rmin_core = inp->core[3];
        rmax_core = inp->core[4];
        zmin_core = inp->core[5];
        zmax_core = inp->core[6];
    }

  /*    // Diagnostic print statements for all individual variables
    printf("\n=== DIAGNOSTIC: Individual Variable Values ===\n");
    
    if (toka_type == GKYL_TOKA_GRID_GEN_DOUBLE_NULL) {
        printf("Double Null Configuration:\n");
        printf("  lower_PF variables:\n");
        printf("    rleft_PF = %g\n", rleft_PF);
        printf("    rright_PF = %g\n", rright_PF);
        printf("    rmin_PF = %g\n", rmin_PF);
        printf("    rmax_PF = %g\n", rmax_PF);
        printf("    zmin_left_PF = %g\n", zmin_left_PF);
        printf("    zmin_right_PF = %g\n", zmin_right_PF);
        
        printf("  outer_SOL variables:\n");
        printf("    rclose_outer_SOL = %g\n", rclose_outer_SOL);
        printf("    rleft_outer_SOL = %g\n", rleft_outer_SOL);
        printf("    rright_outer_SOL = %g\n", rright_outer_SOL);
        printf("    rmin_outer_SOL = %g\n", rmin_outer_SOL);
        printf("    rmax_outer_SOL = %g\n", rmax_outer_SOL);
        printf("    zmin_outer_SOL = %g\n", zmin_outer_SOL);
        printf("    zmax_outer_SOL = %g\n", zmax_outer_SOL);
        
        printf("  inner_SOL variables:\n");
        printf("    rleft_inner_SOL = %g\n", rleft_inner_SOL);
        printf("    rright_inner_SOL = %g\n", rright_inner_SOL);
        printf("    rmin_inner_SOL = %g\n", rmin_inner_SOL);
        printf("    rmax_inner_SOL = %g\n", rmax_inner_SOL);
        printf("    zmin_inner_SOL = %g\n", zmin_inner_SOL);
        printf("    zmax_inner_SOL = %g\n", zmax_inner_SOL);
        
        printf("  core variables:\n");
        printf("    rclose_inner_core = %g\n", rclose_inner_core);
        printf("    rclose_outer_core = %g\n", rclose_outer_core);
        printf("    rleft_core = %g\n", rleft_core);
        printf("    rright_core = %g\n", rright_core);
        printf("    rmin_core = %g\n", rmin_core);
        printf("    rmax_core = %g\n", rmax_core);
    } else if (toka_type == GKYL_TOKA_GRID_GEN_SINGLE_NULL) {
        printf("Single Null Configuration:\n");
        printf("  lower_PF variables:\n");
        printf("    rleft_PF = %g\n", rleft_PF);
        printf("    rright_PF = %g\n", rright_PF);
        printf("    rmin_PF = %g\n", rmin_PF);
        printf("    rmax_PF = %g\n", rmax_PF);
        printf("    zmin_PF = %g\n", zmin_PF);
        printf("    zmax_PF = %g\n", zmax_PF);
        printf("    zmin_left_PF = %g\n", zmin_left_PF);
        printf("    zmin_right_PF = %g\n", zmin_right_PF);
        
        printf("  SOL variables:\n");
        printf("    rclose_SOL = %g\n", rclose_SOL);
        printf("    rleft_SOL = %g\n", rleft_SOL);
        printf("    rright_SOL = %g\n", rright_SOL);
        printf("    rmin_SOL = %g\n", rmin_SOL);
        printf("    rmax_SOL = %g\n", rmax_SOL);
        printf("    zmin_SOL = %g\n", zmin_SOL);
        printf("    zmax_SOL = %g\n", zmax_SOL);
        
        printf("  core variables:\n");
        printf("    rclose_core = %g\n", rclose_core);
        printf("    rleft_core = %g\n", rleft_core);
        printf("    rright_core = %g\n", rright_core);
        printf("    rmin_core = %g\n", rmin_core);
        printf("    rmax_core = %g\n", rmax_core);
        printf("    zmin_core = %g\n", zmin_core);
        printf("    zmax_core = %g\n", zmax_core);
    }

    printf("Loaded values:\n");
    printf("    Cells: [%d, %d, %d, %d, %d, %d]\n", cell_count[0], cell_count[1], cell_count[2], cell_count[3], cell_count[4], cell_count[5]);
    printf("    Psi widths: [%g, %g, %g, %g]\n", wpsi[0], wpsi[1], wpsi[2], wpsi[3]);
    printf("    GEQDSK file: %s\n", inp->geqdsk_path);
    printf("    Separatrix psi: %g\n", psisep);
    printf("    Center psi: %g\n\n", psicenter);
    printf("    Use Half Domain: %s\n", use_half_domain ? "true" : "false");
    printf("=== END DIAGNOSTIC ===\n\n");  */

    double deltapsi;
    if (psisep > psicenter){
        deltapsi = psisep - psicenter; 
    } else {
        deltapsi = psicenter - psisep;
    }

    double wpf, wout, win, wcore, wsol;
    if (toka_type == GKYL_TOKA_GRID_GEN_DOUBLE_NULL) {
        wpf = wpsi[0] * deltapsi;
        wout = wpsi[1] * deltapsi;
        win = wpsi[2] * deltapsi;
        wcore = wpsi[3] * deltapsi;
    } else { // Single null case
        wpf = wpsi[0] * deltapsi;
        wsol = wpsi[1] * deltapsi;
        wcore = wpsi[2] * deltapsi;
    }
 
  
    // Initialize psi values for different regions
    double psi_lo_outer_sol, psi_up_outer_sol;
    double psi_lo_core, psi_up_core;
    double psi_lo_pf, psi_up_pf;
    double psi_lo_inner_sol, psi_up_inner_sol;
    double psi_lo_sol, psi_up_sol;
    if (toka_type == GKYL_TOKA_GRID_GEN_DOUBLE_NULL) {
      if (psisep > psicenter) {
          // For positive psi case (psisep > psicenter)
          psi_lo_outer_sol = psisep;
          psi_up_outer_sol = psisep + wout;
          
          psi_lo_core = psisep - wcore;
          psi_up_core = psisep;
          
          psi_lo_pf = psisep - wpf;
          psi_up_pf = psisep;
          
          psi_lo_inner_sol = psisep;
          psi_up_inner_sol = psisep + win;
      } else {
          // For negative psi case (psisep < psicenter)
          psi_lo_outer_sol = psisep - wout;
          psi_up_outer_sol = psisep;
          
          psi_lo_core = psisep ;
          psi_up_core = psisep + wcore;
          
          psi_lo_pf = psisep;
          psi_up_pf = psisep + wpf;
          
          psi_lo_inner_sol = psisep - win;
          psi_up_inner_sol = psisep;
      }
    } else{
        if (psisep > psicenter) {
          // For positive psi case (psisep > psicenter)
          psi_lo_sol = psisep;
          psi_up_sol = psisep + wsol;
          
          psi_lo_core = psisep - wcore;
          psi_up_core = psisep;
          
          psi_lo_pf = psisep - wpf;
          psi_up_pf = psisep;

      } else {
          // For negative psi case (psisep < psicenter)
          psi_lo_sol = psisep - wsol;
          psi_up_sol = psisep;
          
          psi_lo_core = psisep;
          psi_up_core = psisep + wcore;
          
          psi_lo_pf = psisep;
          psi_up_pf = psisep + wpf;          
      }
    }



    int npsi_pf = cell_count[0];
    int npsi_inner_sol = cell_count[1];
    int npsi_outer_sol = cell_count[1];
    int npsi_core = cell_count[2];


  
    double ntheta_lower_outer = cell_count[3];
    double ntheta_middle_outer = cell_count[4];
    double ntheta_middle_inner = cell_count[4];
    double ntheta_lower_inner  = cell_count[5];

  
    double Lz = (M_PI-1e-14)*2.0;
    double theta_lo = -Lz/2.0, theta_up = Lz/2.0;

    struct gkyl_gk_block_geom *bgeom = NULL;
    
    if (toka_type == GKYL_TOKA_GRID_GEN_DOUBLE_NULL) {     
      if (use_half_domain == true) {
        bgeom = gkyl_gk_block_geom_new(2, 8);
        if (psisep > psicenter){
          // For positive psi case (psisep > psicenter)
          // block 0. Lower outer PF region.
          gkyl_gk_block_geom_set_block(bgeom, 0, &(struct gkyl_gk_block_geom_info) {
              .lower = { psi_lo_pf, theta_lo},
              .upper = { psi_up_pf, theta_up},
              .cells = { npsi_pf, ntheta_lower_outer },
              .cuts = { 1, 1 },
              .geometry = {
                .world = {0.0},
                .geometry_id = GKYL_TOKAMAK,
                .efit_info = efit_inp,
                .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                  .ftype = GKYL_PF_LO_R,
                  .half_domain = true,
                  .rleft = rleft_PF,   // rleft_PF
                  .rright = rright_PF,  // rright_PF
                  .rmin = rmin_PF,    // rmin_PF
                  .rmax = rmax_PF,    // rmax_PF
                  .zmin_left = zmin_left_PF,  // zmin_left_PF
                  .zmin_right = zmin_right_PF, // zmin_right_PF
                  .plate_spec = true,
                  .plate_func_lower = inp->plate_func_lower_outer,
                  .plate_func_upper = inp->plate_func_lower_inner,
                }
              },
              
              .connections[0] = { // x-direction connections
                { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL},
                { .bid = 1, .dir = 0, .edge = GKYL_LOWER_POSITIVE}  // physical boundary
              },
              .connections[1] = { // z-direction connections
                { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
                { .bid = 5, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
              }
            }
          );
        
          // block 1. Lower outer SOL.
          gkyl_gk_block_geom_set_block(bgeom, 1, &(struct gkyl_gk_block_geom_info) {
              .lower = { psi_lo_outer_sol, theta_lo},
              .upper = { psi_up_outer_sol,  theta_up},
              .cells = { npsi_outer_sol, ntheta_lower_outer},
              .cuts = { 1, 1 },
              .geometry = {
                .world = {0.0},
                .geometry_id = GKYL_TOKAMAK,
                .efit_info = efit_inp,
                .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                  .ftype = GKYL_DN_SOL_OUT_LO,
                  .half_domain = true,
                  .rclose = rclose_outer_SOL, // rclose_outer_SOL
                  .rleft = rleft_outer_SOL,  // rleft_outer_SOL
                  .rright = rright_outer_SOL, // rright_outer_SOL
                  .rmin = rmin_outer_SOL,   // rmin_outer_SOL
                  .rmax = rmax_outer_SOL,   // rmax_outer_SOL
                  .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
                  .zmin = zmin_outer_SOL, // zmin_outer_SOL
                  .zmax = zmax_outer_SOL, // zmax_outer_SOL
                  .plate_spec = true,
                  .plate_func_lower = inp->plate_func_lower_outer,
                  .plate_func_upper = inp->plate_func_upper_outer,
                }
              },
              
              .connections[0] = { // x-direction connections
                { .bid = 0, .dir = 0, .edge = GKYL_UPPER_POSITIVE}, // physical boundary
                { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
              },
              .connections[1] = { // z-direction connections
                { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
                { .bid = 2, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
              }
            }
          );
        
          // block 2. Middle outer SOL.
          gkyl_gk_block_geom_set_block(bgeom, 2, &(struct gkyl_gk_block_geom_info) {
              .lower = { psi_lo_outer_sol, theta_lo },
              .upper = { psi_up_outer_sol, theta_up },
              .cells = { npsi_outer_sol, ntheta_middle_outer},
              .cuts = { 1, 1 },
              .geometry = {
                .world = {0.0},
                .geometry_id = GKYL_TOKAMAK,
                .efit_info = efit_inp,
                .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                  .ftype = GKYL_DN_SOL_OUT_MID,
                  .half_domain = true,
                  .rclose = rclose_outer_SOL, // rclose_outer_SOL
                  .rleft = rleft_outer_SOL,  // rleft_outer_SOL
                  .rright = rright_outer_SOL, // rright_outer_SOL
                  .rmin = rmin_outer_SOL,   // rmin_outer_SOL
                  .rmax = rmax_outer_SOL,   // rmax_outer_SOL
                  .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
                  .zmin = zmin_outer_SOL, // zmin_outer_SOL
                  .zmax = zmax_outer_SOL, // zmax_outer_SOL
                  .plate_spec = true,
                  .plate_func_lower = inp->plate_func_lower_outer,
                  .plate_func_upper = inp->plate_func_upper_outer,
                }
              },
              
              .connections[0] = { // x-direction connections
                { .bid = 6, .dir = 0, .edge = GKYL_UPPER_POSITIVE}, // physical boundary
                { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
              },
              .connections[1] = { // z-direction connections
                { .bid = 1, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
                { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL},
              }
            }
          );
          
          // block 7. Middle inner SOL.
          gkyl_gk_block_geom_set_block(bgeom, 3, &(struct gkyl_gk_block_geom_info) {
              .lower = { psi_lo_inner_sol, theta_lo },
              .upper = { psi_up_inner_sol, theta_up },
              .cells = { npsi_inner_sol, ntheta_middle_inner},
              .cuts = { 1, 1 },
              .geometry = {
                .world = {0.0},
                .geometry_id = GKYL_TOKAMAK,
                .efit_info = efit_inp,
                .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                  .ftype = GKYL_DN_SOL_IN_MID,
                  .half_domain = true,
                  .rleft = rleft_inner_SOL,  // rleft_inner_SOL
                  .rright = rright_inner_SOL, // rright_inner_SOL
                  .rmin = rmin_inner_SOL,   // rmin_inner_SOL
                  .rmax = rmax_inner_SOL,   // rmax_inner_SOL
                  .zmin = zmin_inner_SOL,   // zmin_inner_SOL
                  .zmax = zmax_inner_SOL,   // zmax_inner_SOL  
                  .plate_spec = true,
                  .plate_func_upper = inp->plate_func_upper_inner,
                  .plate_func_lower= inp->plate_func_lower_inner,
                }
              },
              
              .connections[0] = { // x-direction connections
                { .bid = 7, .dir = 0, .edge = GKYL_UPPER_POSITIVE}, // physical boundary
                { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}
              },
              .connections[1] = { // z-direction connections
                { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL},
                { .bid = 4, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
              }
            }
          );
        
          // block 8. Lower inner SOL.
          gkyl_gk_block_geom_set_block(bgeom, 4, &(struct gkyl_gk_block_geom_info) {
              .lower = { psi_lo_inner_sol, theta_lo },
              .upper = { psi_up_inner_sol, theta_up },
              .cells = { npsi_inner_sol, ntheta_lower_inner},
              .cuts = { 1, 1 },
              .geometry = {
                .world = {0.0},
                .geometry_id = GKYL_TOKAMAK,
                .efit_info = efit_inp,
                .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                  .ftype = GKYL_DN_SOL_IN_LO,
                  .half_domain = true,
                  .rleft = rleft_inner_SOL,  // rleft_inner_SOL
                  .rright = rright_inner_SOL, // rright_inner_SOL
                  .rmin = rmin_inner_SOL,   // rmin_inner_SOL
                  .rmax = rmax_inner_SOL,   // rmax_inner_SOL
                  .zmin = zmin_inner_SOL,   // zmin_inner_SOL
                  .zmax = zmax_inner_SOL,   // zmax_inner_SOL  
                  .plate_spec = true,
                  .plate_func_upper = inp->plate_func_upper_inner,
                  .plate_func_lower= inp->plate_func_lower_inner,
                }
              },
              
              .connections[0] = { // x-direction connections
                { .bid = 5, .dir = 0, .edge = GKYL_UPPER_POSITIVE}, // physical boundary
                { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}
              },
              .connections[1] = { // z-direction connections
                { .bid = 3, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
                { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}
              }
            }
          );
        
          // block 9. Lower inner PF region.
          gkyl_gk_block_geom_set_block(bgeom, 5, &(struct gkyl_gk_block_geom_info) {
              .lower = { psi_lo_pf, theta_lo},
              .upper = { psi_up_pf, theta_up},
              .cells = { npsi_pf, ntheta_lower_inner },
              .cuts = { 1, 1 },
              .geometry = {
                .world = {0.0},
                .geometry_id = GKYL_TOKAMAK,
                .efit_info = efit_inp,
                .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                  .ftype = GKYL_PF_LO_L,
                  .half_domain = true,
                  .rleft = rleft_PF,   // rleft_PF
                  .rright = rright_PF,  // rright_PF
                  .rmin = rmin_PF,    // rmin_PF
                  .rmax = rmax_PF,    // rmax_PF
                  .zmin_left = zmin_left_PF,  // zmin_left_PF
                  .zmin_right = zmin_right_PF, // zmin_right_PF
                  .plate_spec = true,
                  .plate_func_lower = inp->plate_func_lower_outer,
                  .plate_func_upper = inp->plate_func_lower_inner,
                }
              },
              
              .connections[0] = { // x-direction connections
                { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL},
                { .bid = 4, .dir = 0, .edge = GKYL_LOWER_POSITIVE}  // physical boundary
              },
              .connections[1] = { // z-direction connections
                { .bid = 0, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
                { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL} // physical boundary
              }
            }
          );
        
        
          // block 10. outer core.
          gkyl_gk_block_geom_set_block(bgeom, 6, &(struct gkyl_gk_block_geom_info) {
              .lower = { psi_lo_core, theta_lo},
              .upper = { psi_up_core,  theta_up},
              .cells = { npsi_core, ntheta_middle_outer},
              .cuts = { 1, 1 },
              .geometry = {
                .world = {0.0},
                .geometry_id = GKYL_TOKAMAK,
                .efit_info = efit_inp,
                .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                  .ftype = GKYL_CORE_R,
                  .half_domain = true,
                  .rclose = rclose_outer_core,    // rclose_core
                  .rleft = rleft_core,     // rleft_core
                  .rright = rright_core,    // rright_core
                  .rmin = rmin_core,      // rmin_core
                  .rmax = rmax_core,      // rmax_core
                  .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
                }
              },
              
              .connections[0] = { // x-direction connections
                { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
                { .bid = 2, .dir = 0, .edge = GKYL_LOWER_POSITIVE}, // physical boundary
              },
              .connections[1] = { // z-direction connections
                { .bid = 7, .dir = 1, .edge = GKYL_UPPER_POSITIVE}, // physical boundary
                { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL},
              }
            }
          );
        
          // block 11. Inner Core.
          gkyl_gk_block_geom_set_block(bgeom, 7, &(struct gkyl_gk_block_geom_info) {
              .lower = { psi_lo_core, theta_lo },
              .upper = { psi_up_core, theta_up },
              .cells = { npsi_core, ntheta_middle_inner},
              .cuts = { 1, 1 },
              .geometry = {
                .world = {0.0},
                .geometry_id = GKYL_TOKAMAK,
                .efit_info = efit_inp,
                .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                  .ftype = GKYL_CORE_L,
                  .half_domain = true,
                  .rclose = rclose_inner_core,    // rclose_core
                  .rleft = rleft_core,     // rleft_core
                  .rright = rright_core,    // rright_core
                  .rmin = rmin_core,      // rmin_core
                  .rmax = rmax_core,      // rmax_core
                  .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
                }
              },
              
              .connections[0] = { // x-direction connections
                { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
                { .bid = 3, .dir = 0, .edge = GKYL_LOWER_POSITIVE}, // physical boundary
              },
              .connections[1] = { // z-direction connections
                { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL},
                { .bid = 6, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
              }
            }
          );
        }

        else{ // For negative psi case (psisep < psicenter)
          // block 0. Lower outer PF region.
          gkyl_gk_block_geom_set_block(bgeom, 0, &(struct gkyl_gk_block_geom_info) {
              .lower = { psi_lo_pf, theta_lo},
              .upper = { psi_up_pf, theta_up},
              .cells = { npsi_pf, ntheta_lower_outer },
              .cuts = { 1, 1 },
              .geometry = {
                .world = {0.0},
                .geometry_id = GKYL_TOKAMAK,
                .efit_info = efit_inp,
                .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                  .ftype = GKYL_PF_LO_R,
                  .half_domain = true,
                  .rleft = rleft_PF,   // rleft_PF
                  .rright = rright_PF,  // rright_PF
                  .rmin = rmin_PF,    // rmin_PF
                  .rmax = rmax_PF,    // rmax_PF
                  .zmin_left = zmin_left_PF,  // zmin_left_PF
                  .zmin_right = zmin_right_PF, // zmin_right_PF
                  .plate_spec = true,
                  .plate_func_lower = inp->plate_func_lower_outer,
                  .plate_func_upper = inp->plate_func_lower_inner,
                }
              },
              
              .connections[0] = { // x-direction connections
                { .bid = 1, .dir = 0, .edge = GKYL_UPPER_POSITIVE },
                { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }  // physical boundary
              },
              .connections[1] = { // z-direction connections
                { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
                { .bid = 5, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
              }
            }
          );
        
          // block 1. Lower outer SOL.
          gkyl_gk_block_geom_set_block(bgeom, 1, &(struct gkyl_gk_block_geom_info) {
              .lower = { psi_lo_outer_sol, theta_lo},
              .upper = { psi_up_outer_sol,  theta_up},
              .cells = { npsi_outer_sol, ntheta_lower_outer},
              .cuts = { 1, 1 },
              .geometry = {
                .world = {0.0},
                .geometry_id = GKYL_TOKAMAK,
                .efit_info = efit_inp,
                .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                  .ftype = GKYL_DN_SOL_OUT_LO,
                  .half_domain = true,
                  .rclose = rclose_outer_SOL, // rclose_outer_SOL
                  .rleft = rleft_outer_SOL,  // rleft_outer_SOL
                  .rright = rright_outer_SOL, // rright_outer_SOL
                  .rmin = rmin_outer_SOL,   // rmin_outer_SOL
                  .rmax = rmax_outer_SOL,   // rmax_outer_SOL
                  .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
                  .zmin = zmin_outer_SOL, // zmin_outer_SOL
                  .zmax = zmax_outer_SOL, // zmax_outer_SOL
                  .plate_spec = true,
                  .plate_func_lower = inp->plate_func_lower_outer,
                  .plate_func_upper = inp->plate_func_upper_outer,
                }
              },
              
              .connections[0] = { // x-direction connections
                { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
                { .bid = 0, .dir = 0, .edge = GKYL_LOWER_POSITIVE }
              },
              .connections[1] = { // z-direction connections
                { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
                { .bid = 2, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
              }
            }
          );
        
          // block 2. Middle outer SOL.
          gkyl_gk_block_geom_set_block(bgeom, 2, &(struct gkyl_gk_block_geom_info) {
              .lower = { psi_lo_outer_sol, theta_lo },
              .upper = { psi_up_outer_sol, theta_up },
              .cells = { npsi_outer_sol, ntheta_middle_outer},
              .cuts = { 1, 1 },
              .geometry = {
                .world = {0.0},
                .geometry_id = GKYL_TOKAMAK,
                .efit_info = efit_inp,
                .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                  .ftype = GKYL_DN_SOL_OUT_MID,
                  .half_domain = true,
                  .rclose = rclose_outer_SOL, // rclose_outer_SOL
                  .rleft = rleft_outer_SOL,  // rleft_outer_SOL
                  .rright = rright_outer_SOL, // rright_outer_SOL
                  .rmin = rmin_outer_SOL,   // rmin_outer_SOL
                  .rmax = rmax_outer_SOL,   // rmax_outer_SOL
                  .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
                  .zmin = zmin_outer_SOL, // zmin_outer_SOL
                  .zmax = zmax_outer_SOL, // zmax_outer_SOL
                  .plate_spec = true,
                  .plate_func_lower = inp->plate_func_lower_outer,
                  .plate_func_upper = inp->plate_func_upper_outer,
                }
              },
              
              .connections[0] = { // x-direction connections
                { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
                { .bid = 6, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
              },
              .connections[1] = { // z-direction connections
                { .bid = 1, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
                { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL},
              }
            }
          );
          
          // block 7. Middle inner SOL.
          gkyl_gk_block_geom_set_block(bgeom, 3, &(struct gkyl_gk_block_geom_info) {
              .lower = { psi_lo_inner_sol, theta_lo },
              .upper = { psi_up_inner_sol, theta_up },
              .cells = { npsi_inner_sol, ntheta_middle_inner},
              .cuts = { 1, 1 },
              .geometry = {
                .world = {0.0},
                .geometry_id = GKYL_TOKAMAK,
                .efit_info = efit_inp,
                .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                  .ftype = GKYL_DN_SOL_IN_MID,
                  .half_domain = true,
                  .rleft = rleft_inner_SOL,  // rleft_inner_SOL
                  .rright = rright_inner_SOL, // rright_inner_SOL
                  .rmin = rmin_inner_SOL,   // rmin_inner_SOL
                  .rmax = rmax_inner_SOL,   // rmax_inner_SOL
                  .zmin = zmin_inner_SOL,   // zmin_inner_SOL
                  .zmax = zmax_inner_SOL,   // zmax_inner_SOL  
                  .plate_spec = true,
                  .plate_func_upper = inp->plate_func_upper_inner,
                  .plate_func_lower= inp->plate_func_lower_inner,
                }
              },
              
              .connections[0] = { // x-direction connections
                { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
                { .bid = 7, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
              },
              .connections[1] = { // z-direction connections
                { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL},
                { .bid = 4, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
              }
            }
          );
        
          // block 8. Lower inner SOL.
          gkyl_gk_block_geom_set_block(bgeom, 4, &(struct gkyl_gk_block_geom_info) {
              .lower = { psi_lo_inner_sol, theta_lo },
              .upper = { psi_up_inner_sol, theta_up },
              .cells = { npsi_inner_sol, ntheta_lower_inner},
              .cuts = { 1, 1 },
              .geometry = {
                .world = {0.0},
                .geometry_id = GKYL_TOKAMAK,
                .efit_info = efit_inp,
                .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                  .ftype = GKYL_DN_SOL_IN_LO,
                  .half_domain = true,
                  .rleft = rleft_inner_SOL,  // rleft_inner_SOL
                  .rright = rright_inner_SOL, // rright_inner_SOL
                  .rmin = rmin_inner_SOL,   // rmin_inner_SOL
                  .rmax = rmax_inner_SOL,   // rmax_inner_SOL
                  .zmin = zmin_inner_SOL,   // zmin_inner_SOL
                  .zmax = zmax_inner_SOL,   // zmax_inner_SOL  
                  .plate_spec = true,
                  .plate_func_upper = inp->plate_func_upper_inner,
                  .plate_func_lower= inp->plate_func_lower_inner,
                }
              },
              
              .connections[0] = { // x-direction connections
                { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
                { .bid = 5, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
              },
              .connections[1] = { // z-direction connections
                { .bid = 3, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
                { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}
              }
            }
          );
        
          // block 9. Lower inner PF region.
          gkyl_gk_block_geom_set_block(bgeom, 5, &(struct gkyl_gk_block_geom_info) {
              .lower = { psi_lo_pf, theta_lo},
              .upper = { psi_up_pf, theta_up},
              .cells = { npsi_pf, ntheta_lower_inner },
              .cuts = { 1, 1 },
              .geometry = {
                .world = {0.0},
                .geometry_id = GKYL_TOKAMAK,
                .efit_info = efit_inp,
                .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                  .ftype = GKYL_PF_LO_L,
                  .half_domain = true,
                  .rleft = rleft_PF,   // rleft_PF
                  .rright = rright_PF,  // rright_PF
                  .rmin = rmin_PF,    // rmin_PF
                  .rmax = rmax_PF,    // rmax_PF
                  .zmin_left = zmin_left_PF,  // zmin_left_PF
                  .zmin_right = zmin_right_PF, // zmin_right_PF
                  .plate_spec = true,
                  .plate_func_lower = inp->plate_func_lower_outer,
                  .plate_func_upper = inp->plate_func_lower_inner,
                }
              },
              
              .connections[0] = { // x-direction connections
                { .bid = 4, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
                { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
              },
              .connections[1] = { // z-direction connections
                { .bid = 0, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
                { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL} // physical boundary
              }
            }
          );
        
        
          // block 10. outer core.
          gkyl_gk_block_geom_set_block(bgeom, 6, &(struct gkyl_gk_block_geom_info) {
              .lower = { psi_lo_core, theta_lo},
              .upper = { psi_up_core,  theta_up},
              .cells = { npsi_core, ntheta_middle_outer},
              .cuts = { 1, 1 },
              .geometry = {
                .world = {0.0},
                .geometry_id = GKYL_TOKAMAK,
                .efit_info = efit_inp,
                .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                  .ftype = GKYL_CORE_R,
                  .half_domain = true,
                  .rclose = rclose_outer_core,    // rclose_core
                  .rleft = rleft_core,     // rleft_core
                  .rright = rright_core,    // rright_core
                  .rmin = rmin_core,      // rmin_core
                  .rmax = rmax_core,      // rmax_core
                  .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
                }
              },
              
              .connections[0] = { // x-direction connections
                { .bid = 2, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
                { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
              },
              .connections[1] = { // z-direction connections
                { .bid = 7, .dir = 1, .edge = GKYL_UPPER_POSITIVE}, // physical boundary
                { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL},
              }
            }
          );
        
          // block 11. Inner Core.
          gkyl_gk_block_geom_set_block(bgeom, 7, &(struct gkyl_gk_block_geom_info) {
              .lower = { psi_lo_core, theta_lo },
              .upper = { psi_up_core, theta_up },
              .cells = { npsi_core, ntheta_middle_inner},
              .cuts = { 1, 1 },
              .geometry = {
                .world = {0.0},
                .geometry_id = GKYL_TOKAMAK,
                .efit_info = efit_inp,
                .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                  .ftype = GKYL_CORE_L,
                  .half_domain = true,
                  .rclose = rclose_inner_core,    // rclose_core
                  .rleft = rleft_core,     // rleft_core
                  .rright = rright_core,    // rright_core
                  .rmin = rmin_core,      // rmin_core
                  .rmax = rmax_core,      // rmax_core
                  .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
                }
              },
              
              .connections[0] = { // x-direction connections
                { .bid = 3, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
                { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
              },
              .connections[1] = { // z-direction connections
                { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL},
                { .bid = 6, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
              }
            }
          );
        }
      }
      else { // For use_half_domain == false
        bgeom = gkyl_gk_block_geom_new(2, 12);
        if (psisep > psicenter){
          // For positive psi case (psisep > psicenter)
          // block 0. Lower outer PF region.
          gkyl_gk_block_geom_set_block(bgeom, 0, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_pf, theta_lo},
            .upper = { psi_up_pf, theta_up},
            .cells = { npsi_pf, ntheta_lower_outer },
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_PF_LO_R,
                .rright = rright_PF,
                .rleft = rleft_PF,
                .rmin = rmin_PF,
                .rmax = rmax_PF,
                .zmin_right = zmin_right_PF,
                .zmin_left = zmin_left_PF,
                .plate_spec = true,
                .plate_func_lower = inp->plate_func_lower_outer,
                .plate_func_upper = inp->plate_func_lower_inner,
              }
            },
            
            .connections[0] = { // x-direction connections
              { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL },
              { .bid = 1, .dir = 0, .edge = GKYL_LOWER_POSITIVE }  // physical boundary
            },
            .connections[1] = { // z-direction connections
              { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
              { .bid = 9, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
            }
          }
        );

        // block 1. Lower outer SOL.
        gkyl_gk_block_geom_set_block(bgeom, 1, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_outer_sol, theta_lo},
            .upper = { psi_up_outer_sol,  theta_up},
            .cells = { npsi_outer_sol, ntheta_lower_outer},
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              //.geometry_id = GKYL_GEOMETRY_FROMFILE,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_DN_SOL_OUT_LO,
                .rclose = rclose_outer_SOL,       // Closest R to region of interest
                .rright = rright_outer_SOL,       // Closest R to outboard SOL
                .rleft = rleft_outer_SOL,        // closest R to inboard SOL
                .rmin = rmin_outer_SOL,         // smallest R in machine
                .rmax = rmax_outer_SOL,         // largest R in machine
                .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
                .zmin = zmin_outer_SOL,
                .zmax = zmax_outer_SOL,
                .plate_spec = true,
                .plate_func_lower = inp->plate_func_lower_outer,
                .plate_func_upper = inp->plate_func_upper_outer,
              }
            },
            
            .connections[0] = { // x-direction connections
              { .bid = 0, .dir = 0, .edge = GKYL_UPPER_POSITIVE}, // physical boundary
              { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }
            },
            .connections[1] = { // z-direction connections
              { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
              { .bid = 2, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
            }
          }
        );

        // block 2. Middle outer SOL.
        gkyl_gk_block_geom_set_block(bgeom, 2, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_outer_sol, theta_lo },
            .upper = { psi_up_outer_sol, theta_up },
            .cells = { npsi_outer_sol, ntheta_middle_outer},
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              //.geometry_id = GKYL_GEOMETRY_FROMFILE,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_DN_SOL_OUT_MID,
                .rclose = rclose_outer_SOL,       // Closest R to region of interest
                .rright = rright_outer_SOL,       // Closest R to outboard SOL
                .rleft = rleft_outer_SOL,        // closest R to inboard SOL
                .rmin = rmin_outer_SOL,         // smallest R in machine
                .rmax = rmax_outer_SOL,         // largest R in machine
                .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
                .zmin = zmin_outer_SOL,
                .zmax = zmax_outer_SOL,
                .plate_spec = true,
                .plate_func_lower = inp->plate_func_lower_outer,
                .plate_func_upper = inp->plate_func_upper_outer,
              }
            },
            
            .connections[0] = { // x-direction connections
              { .bid = 10, .dir = 0, .edge = GKYL_UPPER_POSITIVE}, // physical boundary
              { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}
            },
            .connections[1] = { // z-direction connections
              { .bid = 1, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
              { .bid = 3, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
            }
          }
        );

        // block 3. Upper outer SOL.
        gkyl_gk_block_geom_set_block(bgeom, 3, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_outer_sol, theta_lo},
            .upper = { psi_up_outer_sol, theta_up},
            .cells = { npsi_outer_sol, ntheta_lower_outer},
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              //.geometry_id = GKYL_GEOMETRY_FROMFILE,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_DN_SOL_OUT_UP,
                .rclose = rclose_outer_SOL,       // Closest R to region of interest
                .rright = rright_outer_SOL,       // Closest R to outboard SOL
                .rleft = rleft_outer_SOL,        // closest R to inboard SOL
                .rmin = rmin_outer_SOL,         // smallest R in machine
                .rmax = rmax_outer_SOL,         // largest R in machine
                .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
                .zmin = zmin_outer_SOL,
                .zmax = zmax_outer_SOL,
                .plate_spec = true, 
                .plate_func_lower = inp->plate_func_lower_outer,
                .plate_func_upper = inp->plate_func_upper_outer,
              }
            },
            
            .connections[0] = { // x-direction connections
              { .bid = 4, .dir = 0, .edge = GKYL_UPPER_POSITIVE}, // physical boundary
              { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}
            },
            .connections[1] = { // z-direction connections
              { .bid = 2, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
              { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
            }
          }
        );

        // block 4. Upper outer PF region.
        gkyl_gk_block_geom_set_block(bgeom, 4, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_pf, theta_lo},
            .upper = { psi_up_pf, theta_up},
            .cells = { npsi_pf, ntheta_lower_outer},
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_PF_UP_R,
                .rright = rright_PF,
                .rleft = rleft_PF ,
                .rmin = rmin_PF,
                .rmax = rmax_PF,
                .zmax_right = -zmin_right_PF,
                .zmax_left = -zmin_left_PF,
                .plate_spec = true,
                .plate_func_lower = inp->plate_func_upper_inner,
                .plate_func_upper = inp->plate_func_upper_outer,
              }
            },
            
            .connections[0] = { // x-direction connections
              { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL},
              { .bid = 3, .dir = 0, .edge = GKYL_LOWER_POSITIVE}  // physical boundary
            },
            .connections[1] = { // z-direction connections
              { .bid = 5, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
              { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL} // physical boundary
            }
          }
        );

        // block 5. Upper inner PF region.
        gkyl_gk_block_geom_set_block(bgeom, 5, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_pf, theta_lo},
            .upper = { psi_up_pf, theta_up},
            .cells = { npsi_pf, ntheta_lower_inner},
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_PF_UP_L,
                .rright = rright_PF,
                .rleft = rleft_PF ,
                .rmin = rmin_PF,
                .rmax = rmax_PF,
                .zmax_right = -zmin_right_PF,
                .zmax_left = -zmin_left_PF,
                .plate_spec = true,
                .plate_func_lower = inp->plate_func_upper_inner,
                .plate_func_upper = inp->plate_func_upper_outer,
              }
            },
            
            .connections[0] = { // x-direction connections
              { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL},
              { .bid = 6, .dir = 0, .edge = GKYL_LOWER_POSITIVE}  // physical boundary
            },
            .connections[1] = { // z-direction connections
              { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL }, // physical boundary
              { .bid = 4, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
            }
          }
        );

        // block 6. Upper inner SOL.
        gkyl_gk_block_geom_set_block(bgeom, 6, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_inner_sol, theta_lo },
            .upper = { psi_up_inner_sol, theta_up },
            .cells = { npsi_inner_sol, ntheta_lower_inner},
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_DN_SOL_IN_UP,
                .rleft = rleft_inner_SOL,
                .rright= rright_inner_SOL,
                .rmin = rmin_inner_SOL,
                .rmax = rmax_inner_SOL,
                .zmin = zmin_inner_SOL,  
                .zmax = zmax_inner_SOL, 
                .plate_spec = true,
                .plate_func_upper = inp->plate_func_upper_inner,
                .plate_func_lower= inp->plate_func_lower_inner,
              }
            },
            
            .connections[0] = { // x-direction connections
              { .bid = 5, .dir = 0, .edge = GKYL_UPPER_POSITIVE}, // physical boundary
              { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}
            },
            .connections[1] = { // z-direction connections
              { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
              { .bid = 7, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
            }
          }
        );

        // block 7. Middle inner SOL.
        gkyl_gk_block_geom_set_block(bgeom, 7, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_inner_sol, theta_lo },
            .upper = { psi_up_inner_sol, theta_up },
            .cells = { npsi_inner_sol, ntheta_middle_inner},
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_DN_SOL_IN_MID,
                .rleft = rleft_inner_SOL,
                .rright= rright_inner_SOL,
                .rmin = rmin_inner_SOL,
                .rmax = rmax_inner_SOL,
                .zmin = zmin_inner_SOL,  
                .zmax = zmax_inner_SOL, 
                .plate_spec = true,
                .plate_func_upper = inp->plate_func_upper_inner,
                .plate_func_lower= inp->plate_func_lower_inner,
              }
            },
            
            .connections[0] = { // x-direction connections
              { .bid = 11, .dir = 0, .edge = GKYL_UPPER_POSITIVE}, // physical boundary
              { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}
            },
            .connections[1] = { // z-direction connections
              { .bid = 6, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
              { .bid = 8, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
            }
          }
        );

        // block 8. Lower inner SOL.
        gkyl_gk_block_geom_set_block(bgeom, 8, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_inner_sol, theta_lo },
            .upper = { psi_up_inner_sol, theta_up },
            .cells = { npsi_inner_sol, ntheta_lower_inner},
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_DN_SOL_IN_LO,
                .rleft = rleft_inner_SOL,
                .rright= rright_inner_SOL,
                .rmin = rmin_inner_SOL,
                .rmax = rmax_inner_SOL,
                .zmin = zmin_inner_SOL,  
                .zmax = zmax_inner_SOL, 
                .plate_spec = true,
                .plate_func_upper = inp->plate_func_upper_inner,
                .plate_func_lower= inp->plate_func_lower_inner,
              }
            },
            
            .connections[0] = { // x-direction connections
              { .bid = 9, .dir = 0, .edge = GKYL_LOWER_POSITIVE}, // physical boundary
              { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}
            },
            .connections[1] = { // z-direction connections
              { .bid = 7, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
              { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}
            }
          }
        );

        // block 9. Lower inner PF region.
        gkyl_gk_block_geom_set_block(bgeom, 9, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_pf, theta_lo},
            .upper = { psi_up_pf, theta_up},
            .cells = { npsi_pf, ntheta_lower_inner },
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_PF_LO_L,
                .rleft = rleft_PF,   // rleft_PF
                .rright = rright_PF,  // rright_PF
                .rmin = rmin_PF,    // rmin_PF
                .rmax = rmax_PF,    // rmax_PF
                .zmin_left = zmin_left_PF,  // zmin_left_PF
                .zmin_right = zmin_right_PF, // zmin_right_PF
                .plate_spec = true,
                .plate_func_lower = inp->plate_func_lower_outer,
                .plate_func_upper = inp->plate_func_lower_inner,
              }
            },
            
            .connections[0] = { // x-direction connections
              { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL},
              { .bid = 8, .dir = 0, .edge = GKYL_LOWER_POSITIVE}  // physical boundary
            },
            .connections[1] = { // z-direction connections
              { .bid = 0, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
              { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL} // physical boundary
            }
          }
        );


        // block 10. outer core.
        gkyl_gk_block_geom_set_block(bgeom, 10, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_core, theta_lo},
            .upper = { psi_up_core,  theta_up},
            .cells = { npsi_core, ntheta_middle_outer},
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_CORE_R,
                .rclose = rclose_outer_core,    // rclose_core
                .rleft = rleft_core,     // rleft_core
                .rright = rright_core,    // rright_core
                .rmin = rmin_core,      // rmin_core
                .rmax = rmax_core,      // rmax_core
                .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
              }
            },
            
            .connections[0] = { // x-direction connections
              { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL},
              { .bid = 2, .dir = 0, .edge = GKYL_LOWER_POSITIVE}  // physical boundary
            },
            .connections[1] = { // z-direction connections
              { .bid = 11, .dir = 1, .edge = GKYL_UPPER_POSITIVE}, // physical boundary
              { .bid = 11, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
            }
          }
        );

        // block 11. Inner Core.
        gkyl_gk_block_geom_set_block(bgeom, 11, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_core, theta_lo },
            .upper = { psi_up_core, theta_up },
            .cells = { npsi_core, ntheta_middle_inner},
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_CORE_L,
                .rclose = rclose_inner_core,    // rclose_core
                .rleft = rleft_core,     // rleft_core
                .rright = rright_core,    // rright_core
                .rmin = rmin_core,      // rmin_core
                .rmax = rmax_core,      // rmax_core
                .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
              }
            },
            
            .connections[0] = { // x-direction connections
              { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL},
              { .bid = 7, .dir = 0, .edge = GKYL_LOWER_POSITIVE}  // physical boundary
            },
            .connections[1] = { // z-direction connections
              { .bid = 10, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
              { .bid = 10, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
            }
          }
        );          

        }
        else{  // For negative psi case (psisep < psicenter)
          // block 0. Lower outer PF region.
          gkyl_gk_block_geom_set_block(bgeom, 0, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_pf, theta_lo},
            .upper = { psi_up_pf, theta_up},
            .cells = { npsi_pf, ntheta_lower_outer },
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_PF_LO_R,
                .rright = rright_PF,
                .rleft = rleft_PF,
                .rmin = rmin_PF,
                .rmax = rmax_PF,
                .zmin_right = zmin_right_PF,
                .zmin_left = zmin_left_PF,
                .plate_spec = true,
                .plate_func_lower = inp->plate_func_lower_outer,
                .plate_func_upper = inp->plate_func_lower_inner,
              }
            },
            
            .connections[0] = { // x-direction connections
              { .bid = 1, .dir = 0, .edge = GKYL_UPPER_POSITIVE },
              { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }  // physical boundary
            },
            .connections[1] = { // z-direction connections
              { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
              { .bid = 9, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
            }
          }
        );

        // block 1. Lower outer SOL.
        gkyl_gk_block_geom_set_block(bgeom, 1, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_outer_sol, theta_lo},
            .upper = { psi_up_outer_sol,  theta_up},
            .cells = { npsi_outer_sol, ntheta_lower_outer},
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              //.geometry_id = GKYL_GEOMETRY_FROMFILE,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_DN_SOL_OUT_LO,
                .rclose = rclose_outer_SOL,       // Closest R to region of interest
                .rright = rright_outer_SOL,       // Closest R to outboard SOL
                .rleft = rleft_outer_SOL,        // closest R to inboard SOL
                .rmin = rmin_outer_SOL,         // smallest R in machine
                .rmax = rmax_outer_SOL,         // largest R in machine
                .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
                .zmin = zmin_outer_SOL,
                .zmax = zmax_outer_SOL,
                .plate_spec = true,
                .plate_func_lower = inp->plate_func_lower_outer,
                .plate_func_upper = inp->plate_func_upper_outer,
              }
            },
            
            .connections[0] = { // x-direction connections
              { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
              { .bid = 0, .dir = 0, .edge = GKYL_LOWER_POSITIVE }
            },
            .connections[1] = { // z-direction connections
              { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
              { .bid = 2, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
            }
          }
        );

        // block 2. Middle outer SOL.
        gkyl_gk_block_geom_set_block(bgeom, 2, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_outer_sol, theta_lo },
            .upper = { psi_up_outer_sol, theta_up },
            .cells = { npsi_outer_sol, ntheta_middle_outer},
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              //.geometry_id = GKYL_GEOMETRY_FROMFILE,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_DN_SOL_OUT_MID,
                .rclose = rclose_outer_SOL,       // Closest R to region of interest
                .rright = rright_outer_SOL,       // Closest R to outboard SOL
                .rleft = rleft_outer_SOL,        // closest R to inboard SOL
                .rmin = rmin_outer_SOL,         // smallest R in machine
                .rmax = rmax_outer_SOL,         // largest R in machine
                .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
                .zmin = zmin_outer_SOL,
                .zmax = zmax_outer_SOL,
                .plate_spec = true,
                .plate_func_lower = inp->plate_func_lower_outer,
                .plate_func_upper = inp->plate_func_upper_outer,
              }
            },
            
            .connections[0] = { // x-direction connections
              { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
              { .bid = 10, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
            },
            .connections[1] = { // z-direction connections
              { .bid = 1, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
              { .bid = 3, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
            }
          }
        );

        // block 3. Upper outer SOL.
        gkyl_gk_block_geom_set_block(bgeom, 3, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_outer_sol, theta_lo},
            .upper = { psi_up_outer_sol, theta_up},
            .cells = { npsi_outer_sol, ntheta_lower_outer},
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              //.geometry_id = GKYL_GEOMETRY_FROMFILE,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_DN_SOL_OUT_UP,
                .rclose = rclose_outer_SOL,       // Closest R to region of interest
                .rright = rright_outer_SOL,       // Closest R to outboard SOL
                .rleft = rleft_outer_SOL,        // closest R to inboard SOL
                .rmin = rmin_outer_SOL,         // smallest R in machine
                .rmax = rmax_outer_SOL,         // largest R in machine
                .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
                .zmin = zmin_outer_SOL,
                .zmax = zmax_outer_SOL,
                .plate_spec = true, 
                .plate_func_lower = inp->plate_func_lower_outer,
                .plate_func_upper = inp->plate_func_upper_outer,
              }
            },
            
            .connections[0] = { // x-direction connections
              { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
              { .bid = 4, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
            },
            .connections[1] = { // z-direction connections
              { .bid = 2, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
              { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
            }
          }
        );

        // block 4. Upper outer PF region.
        gkyl_gk_block_geom_set_block(bgeom, 4, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_pf, theta_lo},
            .upper = { psi_up_pf, theta_up},
            .cells = { npsi_pf, ntheta_lower_outer},
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_PF_UP_R,
                .rright = rright_PF,
                .rleft = rleft_PF ,
                .rmin = rmin_PF,
                .rmax = rmax_PF,
                .zmax_right = -zmin_right_PF,
                .zmax_left = -zmin_left_PF,
                .plate_spec = true,
                .plate_func_lower = inp->plate_func_upper_inner,
                .plate_func_upper = inp->plate_func_upper_outer,
              }
            },
            
            .connections[0] = { // x-direction connections
              { .bid = 3, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
              { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
            },
            .connections[1] = { // z-direction connections
              { .bid = 5, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
              { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL} // physical boundary
            }
          }
        );

        // block 5. Upper inner PF region.
        gkyl_gk_block_geom_set_block(bgeom, 5, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_pf, theta_lo},
            .upper = { psi_up_pf, theta_up},
            .cells = { npsi_pf, ntheta_lower_inner},
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_PF_UP_L,
                .rright = rright_PF,
                .rleft = rleft_PF ,
                .rmin = rmin_PF,
                .rmax = rmax_PF,
                .zmax_right = -zmin_right_PF,
                .zmax_left = -zmin_left_PF,
                .plate_spec = true,
                .plate_func_lower = inp->plate_func_upper_inner,
                .plate_func_upper = inp->plate_func_upper_outer,
              }
            },
            
            .connections[0] = { // x-direction connections
              { .bid = 6, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
              { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
            },
            .connections[1] = { // z-direction connections
              { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL }, // physical boundary
              { .bid = 4, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
            }
          }
        );

        // block 6. Upper inner SOL.
        gkyl_gk_block_geom_set_block(bgeom, 6, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_inner_sol, theta_lo },
            .upper = { psi_up_inner_sol, theta_up },
            .cells = { npsi_inner_sol, ntheta_lower_inner},
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_DN_SOL_IN_UP,
                .rleft = rleft_inner_SOL,
                .rright= rright_inner_SOL,
                .rmin = rmin_inner_SOL,
                .rmax = rmax_inner_SOL,
                .zmin = zmin_inner_SOL,  
                .zmax = zmax_inner_SOL, 
                .plate_spec = true,
                .plate_func_upper = inp->plate_func_upper_inner,
                .plate_func_lower= inp->plate_func_lower_inner,
              }
            },
            
            .connections[0] = { // x-direction connections
              { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
              { .bid = 5, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
            },
            .connections[1] = { // z-direction connections
              { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // physical boundary
              { .bid = 7, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
            }
          }
        );

        // block 7. Middle inner SOL.
        gkyl_gk_block_geom_set_block(bgeom, 7, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_inner_sol, theta_lo },
            .upper = { psi_up_inner_sol, theta_up },
            .cells = { npsi_inner_sol, ntheta_middle_inner},
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_DN_SOL_IN_MID,
                .rleft = rleft_inner_SOL,
                .rright= rright_inner_SOL,
                .rmin = rmin_inner_SOL,
                .rmax = rmax_inner_SOL,
                .zmin = zmin_inner_SOL,  
                .zmax = zmax_inner_SOL, 
                .plate_spec = true,
                .plate_func_upper = inp->plate_func_upper_inner,
                .plate_func_lower= inp->plate_func_lower_inner,
              }
            },
            
            .connections[0] = { // x-direction connections
              { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
              { .bid = 11, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
            },
            .connections[1] = { // z-direction connections
              { .bid = 6, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
              { .bid = 8, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
            }
          }
        );

        // block 8. Lower inner SOL.
        gkyl_gk_block_geom_set_block(bgeom, 8, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_inner_sol, theta_lo },
            .upper = { psi_up_inner_sol, theta_up },
            .cells = { npsi_inner_sol, ntheta_lower_inner},
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_DN_SOL_IN_LO,
                .rleft = rleft_inner_SOL,
                .rright= rright_inner_SOL,
                .rmin = rmin_inner_SOL,
                .rmax = rmax_inner_SOL,
                .zmin = zmin_inner_SOL,  
                .zmax = zmax_inner_SOL, 
                .plate_spec = true,
                .plate_func_upper = inp->plate_func_upper_inner,
                .plate_func_lower= inp->plate_func_lower_inner,
              }
            },
            
            .connections[0] = { // x-direction connections
              { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}, // physical boundary
              { .bid = 9, .dir = 0, .edge = GKYL_LOWER_POSITIVE}
            },
            .connections[1] = { // z-direction connections
              { .bid = 7, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
              { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}
            }
          }
        );

        // block 9. Lower inner PF region.
        gkyl_gk_block_geom_set_block(bgeom, 9, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_pf, theta_lo},
            .upper = { psi_up_pf, theta_up},
            .cells = { npsi_pf, ntheta_lower_inner },
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_PF_LO_L,
                .rleft = rleft_PF,   // rleft_PF
                .rright = rright_PF,  // rright_PF
                .rmin = rmin_PF,    // rmin_PF
                .rmax = rmax_PF,    // rmax_PF
                .zmin_left = zmin_left_PF,  // zmin_left_PF
                .zmin_right = zmin_right_PF, // zmin_right_PF
                .plate_spec = true,
                .plate_func_lower = inp->plate_func_lower_outer,
                .plate_func_upper = inp->plate_func_lower_inner,
              }
            },
            
            .connections[0] = { // x-direction connections
              { .bid = 8, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
              { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
            },
            .connections[1] = { // z-direction connections
              { .bid = 0, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
              { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL} // physical boundary
            }
          }
        );


        // block 10. outer core.
        gkyl_gk_block_geom_set_block(bgeom, 10, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_core, theta_lo},
            .upper = { psi_up_core,  theta_up},
            .cells = { npsi_core, ntheta_middle_outer},
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_CORE_R,
                .rclose = rclose_outer_core,    // rclose_core
                .rleft = rleft_core,     // rleft_core
                .rright = rright_core,    // rright_core
                .rmin = rmin_core,      // rmin_core
                .rmax = rmax_core,      // rmax_core
                .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
              }
            },
            
            .connections[0] = { // x-direction connections
              { .bid = 2, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
              { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
            },
            .connections[1] = { // z-direction connections
              { .bid = 11, .dir = 1, .edge = GKYL_UPPER_POSITIVE}, // physical boundary
              { .bid = 11, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
            }
          }
        );

        // block 11. Inner Core.
        gkyl_gk_block_geom_set_block(bgeom, 11, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_core, theta_lo },
            .upper = { psi_up_core, theta_up },
            .cells = { npsi_core, ntheta_middle_inner},
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_CORE_L,
                .rclose = rclose_inner_core,    // rclose_core
                .rleft = rleft_core,     // rleft_core
                .rright = rright_core,    // rright_core
                .rmin = rmin_core,      // rmin_core
                .rmax = rmax_core,      // rmax_core
                .use_cubics = false, // Whether to use cubic representation of psi(R,Z) for field line tracing
              }
            },
            
            .connections[0] = { // x-direction connections
              { .bid = 7, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
              { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}  // physical boundary
            },
            .connections[1] = { // z-direction connections
              { .bid = 10, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
              { .bid = 10, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
            }
          }
        );


        }
      }
      }

    else if (toka_type == GKYL_TOKA_GRID_GEN_SINGLE_NULL) {       
      bgeom = gkyl_gk_block_geom_new(2, 6);
      if (psisep > psicenter){
        // For positive psi case (psisep > psicenter)
        // Block 0: outer private flux (PF) region.
        gkyl_gk_block_geom_set_block(bgeom, 0, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_pf, theta_lo},
            .upper = { psi_up_pf, theta_up},
            .cells = { npsi_pf, ntheta_lower_outer },
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_PF_LO_R,
                .rleft = rleft_PF,
                .rright = rright_PF,
                .rmin = rmin_PF,
                .rmax = rmax_PF,
                .zmin = zmin_PF,
                .zmax = zmax_PF,
                .zmin_left = zmin_left_PF,
                .zmax_right = zmin_right_PF,
                .plate_spec = true,
                .plate_func_lower = inp->plate_func_lower_outer,
                .plate_func_upper = inp->plate_func_lower_inner,
              }
            },
      
            .connections[0] = { // x-direction.
              { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL },  // Physical boundary.
              { .bid = 1, .dir = 0, .edge = GKYL_LOWER_POSITIVE },
            },
            .connections[1] = { // z-direction.
              { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // Physical boundary.
              { .bid = 4, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
            }
          }
        );
      
        // Block 1: lower outer SOL.
        gkyl_gk_block_geom_set_block(bgeom, 1, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_sol, theta_lo},
            .upper = { psi_up_sol,  theta_up},
            .cells = { npsi_outer_sol, ntheta_lower_outer},
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_LSN_SOL_LO,
                .rclose = rclose_SOL,
                .rleft = rleft_SOL,
                .rright = rright_SOL,
                .rmin = rmin_SOL,
                .rmax = rmax_SOL,
                .zmin = zmin_SOL,
                .zmax = zmax_SOL,
                .plate_spec = true,
                .plate_func_lower = inp->plate_func_lower_outer,
                .plate_func_upper = inp->plate_func_lower_inner,
              }
            },
            
            .connections[0] = { // x-direction.
              { .bid = 0, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
              { .bid = 1, .dir = 0, .edge = GKYL_PHYSICAL }, // Physical boundary.
            },
            .connections[1] = { // z-direction.
              { .bid = 1, .dir = 1, .edge = GKYL_PHYSICAL}, // Physical boundary.
              { .bid = 2, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
            }
          }
        );
      
        // Block 2: mid SOL.
        gkyl_gk_block_geom_set_block(bgeom, 2, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_sol, theta_lo },
            .upper = { psi_up_sol, theta_up },
            .cells = { npsi_outer_sol, ntheta_middle_outer},
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_LSN_SOL_MID,
                .rclose = rclose_SOL,
                .rleft = rleft_SOL,
                .rright = rright_SOL,
                .rmin = rmin_SOL,
                .rmax = rmax_SOL,
                .zmin = zmin_SOL,
                .zmax = zmax_SOL,
                .plate_spec = true,
                .plate_func_lower = inp->plate_func_lower_outer,
                .plate_func_upper = inp->plate_func_lower_inner,
              }
            },
            
            .connections[0] = { // x-direction.
              { .bid = 5, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
              { .bid = 2, .dir = 0, .edge = GKYL_PHYSICAL }, // Physical boundary.
            },
            .connections[1] = { // z-direction.
              { .bid = 1, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
              { .bid = 3, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
            }
          }
        );
      
        // Block 3: lower inner SOL.
        gkyl_gk_block_geom_set_block(bgeom, 3, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_sol, theta_lo },
            .upper = { psi_up_sol, theta_up },
            .cells = { npsi_inner_sol, ntheta_lower_inner},
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_LSN_SOL_UP,
                .rclose = rclose_SOL,
                .rleft = rleft_SOL,
                .rright = rright_SOL,
                .rmin = rmin_SOL,
                .rmax = rmax_SOL,
                .zmin = zmin_SOL,
                .zmax = zmax_SOL,
                .plate_spec = true,
                .plate_func_lower = inp->plate_func_lower_outer,
                .plate_func_upper = inp->plate_func_lower_inner,
              }
            },
            
            .connections[0] = { // x-direction.
              { .bid = 4, .dir = 0, .edge = GKYL_UPPER_POSITIVE},
              { .bid = 3, .dir = 0, .edge = GKYL_PHYSICAL }, // Physical boundary.
            },
            .connections[1] = { // z-direction.
              { .bid = 2, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
              { .bid = 3, .dir = 1, .edge = GKYL_PHYSICAL}, // Physical boundary.
            }
          }
        );
      
        // Block 4: inner private flux (PF) region.
        gkyl_gk_block_geom_set_block(bgeom, 4, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_pf, theta_lo},
            .upper = { psi_up_pf, theta_up},
            .cells = { npsi_pf, ntheta_lower_inner },
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_PF_LO_L,
                .rleft = rleft_PF,
                .rright = rright_PF,
                .rmin = rmin_PF,
                .rmax = rmax_PF,
                .zmin = zmin_PF,
                .zmax = zmax_PF,
                .zmin_left = zmin_left_PF,
                .zmax_right = zmin_right_PF,
                .plate_spec = true,
                .plate_func_lower = inp->plate_func_lower_outer,
                .plate_func_upper = inp->plate_func_lower_inner,
              }
            },
      
            .connections[0] = { // x-direction.
              { .bid = 4, .dir = 0, .edge = GKYL_PHYSICAL },  // Physical boundary.
              { .bid = 3, .dir = 0, .edge = GKYL_LOWER_POSITIVE },
            },
            .connections[1] = { // z-direction.
              { .bid = 0, .dir = 1, .edge = GKYL_UPPER_POSITIVE}, // Physical boundary.
              { .bid = 4, .dir = 1, .edge = GKYL_PHYSICAL},
            }
          }
        );
      
        // Block 5: core region.
        gkyl_gk_block_geom_set_block(bgeom, 5, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_core, theta_lo},
            .upper = { psi_up_core,  theta_up},
            .cells = { npsi_core, ntheta_middle_outer},
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_CORE,
                .rclose = rclose_core,
                .rleft = rleft_core,
                .rright = rright_core,
                .rmin = rmin_core,
                .rmax = rmax_core,
                .zmin = zmin_core,
                .zmax = zmax_core,
              }
            },
      
            .connections[0] = { // x-direction.
              { .bid = 5, .dir = 0, .edge = GKYL_PHYSICAL },  // Physical boundary.
              { .bid = 2, .dir = 0, .edge = GKYL_LOWER_POSITIVE },
            },
            .connections[1] = { // z-direction.
              { .bid = 5, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
              { .bid = 5, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
            }
          }
        );
      }
      
      else{
        // For negative psi case (psisep < psicenter)
        // Block 0: outer private flux (PF) region.
        gkyl_gk_block_geom_set_block(bgeom, 0, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_pf, theta_lo},
            .upper = { psi_up_pf, theta_up},
            .cells = { npsi_pf, ntheta_lower_outer },
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_PF_LO_R,
                .rleft = rleft_PF,   // rleft_PF
                .rright = rright_PF,  // rright_PF
                .rmin = rmin_PF,    // rmin_PF
                .rmax = rmax_PF,    // rmax_PF
                .zmin = zmin_PF,    // zmin_PF
                .zmax = zmax_PF,    // zmax_PF
                .zmin_left = zmin_left_PF,  // zmin_left_PF
                .zmax_right = zmin_right_PF, // zmin_right_PF
                .plate_spec = true,
                .plate_func_lower = inp->plate_func_lower_outer,
                .plate_func_upper = inp->plate_func_lower_inner,
              }
            },
      
            .connections[0] = { // x-direction.
              { .bid = 1, .dir = 0, .edge = GKYL_UPPER_POSITIVE },
              { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL },  // Physical boundary.
            },
            .connections[1] = { // z-direction.
              { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL}, // Physical boundary.
              { .bid = 4, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
            }
          }
        );
      
        // Block 1: lower outer SOL.
        gkyl_gk_block_geom_set_block(bgeom, 1, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_sol, theta_lo},
            .upper = { psi_up_sol,  theta_up},
            .cells = { npsi_outer_sol, ntheta_lower_outer},
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_LSN_SOL_LO,
                .rclose = rclose_SOL,     // rclose_SOL
                .rleft = rleft_SOL,      // rleft_SOL
                .rright = rright_SOL,     // rright_SOL
                .rmin = rmin_SOL,       // rmin_SOL
                .rmax = rmax_SOL,       // rmax_SOL
                .zmin = zmin_SOL,       // zmin_SOL
                .zmax = zmax_SOL,       // zmax_SOL
                .plate_spec = true,
                .plate_func_lower = inp->plate_func_lower_outer,
                .plate_func_upper = inp->plate_func_lower_inner,
              }
            },
            
            .connections[0] = { // x-direction.
              { .bid = 1, .dir = 0, .edge = GKYL_PHYSICAL }, // Physical boundary.
              { .bid = 0, .dir = 0, .edge = GKYL_LOWER_POSITIVE},
            },
            .connections[1] = { // z-direction.
              { .bid = 1, .dir = 1, .edge = GKYL_PHYSICAL}, // Physical boundary.
              { .bid = 2, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
            }
          }
        );
      
        // Block 2: mid SOL.
        gkyl_gk_block_geom_set_block(bgeom, 2, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_sol, theta_lo },
            .upper = { psi_up_sol, theta_up },
            .cells = { npsi_outer_sol, ntheta_middle_outer},
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_LSN_SOL_MID,
                .rclose = rclose_SOL,     // rclose_SOL
                .rleft = rleft_SOL,      // rleft_SOL
                .rright = rright_SOL,     // rright_SOL
                .rmin = rmin_SOL,       // rmin_SOL
                .rmax = rmax_SOL,       // rmax_SOL
                .zmin = zmin_SOL,       // zmin_SOL
                .zmax = zmax_SOL,       // zmax_SOL
                .plate_spec = true,
                .plate_func_lower = inp->plate_func_lower_outer,
                .plate_func_upper = inp->plate_func_lower_inner,
              }
            },
            
            .connections[0] = { // x-direction.
              { .bid = 2, .dir = 0, .edge = GKYL_PHYSICAL }, // Physical boundary.
              { .bid = 5, .dir = 0, .edge = GKYL_LOWER_POSITIVE},
            },
            .connections[1] = { // z-direction.
              { .bid = 1, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
              { .bid = 3, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
            }
          }
        );
      
        // Block 3: lower inner SOL.
        gkyl_gk_block_geom_set_block(bgeom, 3, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_sol, theta_lo },
            .upper = { psi_up_sol, theta_up },
            .cells = { npsi_inner_sol, ntheta_lower_inner},
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_LSN_SOL_UP,
                .rclose = rclose_SOL,     // rclose_SOL
                .rleft = rleft_SOL,      // rleft_SOL
                .rright = rright_SOL,     // rright_SOL
                .rmin = rmin_SOL,       // rmin_SOL
                .rmax = rmax_SOL,       // rmax_SOL
                .zmin = zmin_SOL,       // zmin_SOL
                .zmax = zmax_SOL,       // zmax_SOL
                .plate_spec = true,
                .plate_func_lower = inp->plate_func_lower_outer,
                .plate_func_upper = inp->plate_func_lower_inner,
              }
            },
            
            .connections[0] = { // x-direction.
              { .bid = 3, .dir = 0, .edge = GKYL_PHYSICAL }, // Physical boundary.
              { .bid = 4, .dir = 0, .edge = GKYL_LOWER_POSITIVE},
            },
            .connections[1] = { // z-direction.
              { .bid = 2, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
              { .bid = 3, .dir = 1, .edge = GKYL_PHYSICAL}, // Physical boundary.
            }
          }
        );
      
        // Block 4: inner private flux (PF) region.
        gkyl_gk_block_geom_set_block(bgeom, 4, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_pf, theta_lo},
            .upper = { psi_up_pf, theta_up},
            .cells = { npsi_pf, ntheta_lower_inner },
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_PF_LO_L,
                .rleft = rleft_PF,   // rleft_PF
                .rright = rright_PF,  // rright_PF
                .rmin = rmin_PF,    // rmin_PF
                .rmax = rmax_PF,    // rmax_PF
                .zmin = zmin_PF,    // zmin_PF
                .zmax = zmax_PF,    // zmax_PF
                .zmin_left = zmin_left_PF,  // zmin_left_PF
                .zmax_right = zmin_right_PF, // zmin_right_PF
                .plate_spec = true,
                .plate_func_lower = inp->plate_func_lower_outer,
                .plate_func_upper = inp->plate_func_lower_inner,
              }
            },
      
            .connections[0] = { // x-direction.
              { .bid = 3, .dir = 0, .edge = GKYL_UPPER_POSITIVE },
              { .bid = 4, .dir = 0, .edge = GKYL_PHYSICAL },  // Physical boundary.
            },
            .connections[1] = { // z-direction.
              { .bid = 0, .dir = 1, .edge = GKYL_UPPER_POSITIVE}, // Physical boundary.
              { .bid = 4, .dir = 1, .edge = GKYL_PHYSICAL},
            }
          }
        );
      
        // Block 5: core region.
        gkyl_gk_block_geom_set_block(bgeom, 5, &(struct gkyl_gk_block_geom_info) {
            .lower = { psi_lo_core, theta_lo},
            .upper = { psi_up_core,  theta_up},
            .cells = { npsi_core, ntheta_middle_outer},
            .cuts = { 1, 1 },
            .geometry = {
              .world = {0.0},
              .geometry_id = GKYL_TOKAMAK,
              .efit_info = efit_inp,
              .tok_grid_info = (struct gkyl_tok_geo_grid_inp) {
                .ftype = GKYL_CORE,
                .rclose = rclose_core,    // rclose_outer_core
                .rleft = rleft_core,     // rleft_core
                .rright = rright_core,    // rright_core
                .rmin = rmin_core,      // rmin_core
                .rmax = rmax_core,      // rmax_core
                .zmin = zmin_core,      // zmin_core
                .zmax = zmax_core,      // zmax_core
              }
            },
      
            .connections[0] = { // x-direction.
              { .bid = 2, .dir = 0, .edge = GKYL_UPPER_POSITIVE },
              { .bid = 5, .dir = 0, .edge = GKYL_PHYSICAL },  // Physical boundary.
            },
            .connections[1] = { // z-direction.
              { .bid = 5, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
              { .bid = 5, .dir = 1, .edge = GKYL_LOWER_POSITIVE},
            }
          }
        );
      }

    }
     printf("Block geometry creation %s\n", bgeom ? "successful" : "failed");

     // Construct communicator for use in app.
    struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(true, false, stderr);

    struct gkyl_gyrokinetic_multib app_inp = {
      .name = "",

      .cdim = 2, .vdim = 2,
      .poly_order = 1,
      .basis_type = GKYL_BASIS_MODAL_SERENDIPITY,
      .cfl_frac = 1.0,

      .gk_block_geom = bgeom,
      .comm = comm
    };

    // Set application name safely into fixed-size array
    snprintf(app_inp.name, sizeof(app_inp.name), "%s", name ? name : "tokagrid");

    // Create app object.
    struct gkyl_gyrokinetic_multib_app *app = gkyl_gyrokinetic_multib_app_new_geom(&app_inp);


    freeresources:
    // Free resources after simulation completion.
    gkyl_gyrokinetic_multib_app_release_geom(app);
    gkyl_gyrokinetic_comms_release(comm);  
 
    // Create settings struct with the values from the function execution
    struct gkylt_tokagridgen_settings settings = {
        .toka_type = inp->toka_type,
        .use_half_domain = inp->use_half_domain,
        .psisep = psisep  // psisep is already calculated in the function
    };
    
    // Create and return the result struct
    struct gkylt_tokagridgen_output output = {
        .bgeom = bgeom,
        .settings = settings
    };
    
    return output;

}

