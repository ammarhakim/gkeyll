#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_position_map.h>
#include <gkyl_util.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_gk_field_priv.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_comm_io.h>

#include <assert.h>
#include <float.h>
#include <time.h>

static void
gk_field_3x_write_twistshift(struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  // Write the discretized shift (for TS BCs) to file.
  int comm_rank, comm_size;
  gkyl_comm_get_rank(app->comm, &comm_rank);
  gkyl_comm_get_size(app->comm, &comm_size);

  const char *vars[] = {"x","y","z"};
  const char *edge[] = {"lower","upper"};
  const char *fmt = "%s-bc_%s%s_twistshift.gkyl";

  struct gk_species *gks = &app->species[0];
  for (int i = 0; i < 2*app->cdim; i++) {
    if (gks->info.bcs[i].type == GKYL_BC_GK_SPECIES_IWL ||
        gks->info.bcs[i].type == GKYL_BC_GK_SPECIES_TWISTSHIFT) {

      int dir = gks->info.bcs[i].dir;
      int edi = gks->info.bcs[i].edge;
      if (comm_rank == 0 && edi == GKYL_LOWER_EDGE) {
        struct gkyl_bc_twistshift *bc_ts = f->bc_ts_lo;
        
        struct gkyl_rect_grid shear_grid;
        struct gkyl_range shear_r;
        struct gkyl_basis shift_b;
        struct gkyl_array *shift_dg = gkyl_bc_twistshift_get_shift_objects(bc_ts, &shear_grid, &shear_r, &shift_b);

        // Twistshift updater stores the shift on a restricted range (the core) but a full
        // grid. Create a restricted grid for I/O.
        struct gkyl_rect_grid shear_grid_core;
        double lower[1], upper[1];
        int cells[] = {shear_r.volume};
        if (app->gk_geom->geqdsk_sign_convention == 0) {
          // x increases towards SOL.
          lower[0] = shear_grid.lower[0];
          upper[0] = shear_grid.lower[0] + shear_grid.dx[0]*cells[0];
          gkyl_rect_grid_init(&shear_grid_core, shear_grid.ndim, lower, upper, cells);
        }
        else {
          // x increases towards SOL.
          lower[0] = shear_grid.upper[0] - shear_grid.dx[0]*cells[0];
          upper[0] = shear_grid.upper[0];
          gkyl_rect_grid_init(&shear_grid_core, shear_grid.ndim, lower, upper, cells);
        }

        // Package metadata for shift file.
        struct gkyl_msgpack_map_elem io_meta_shift_dg[] = {
          { .key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = shift_b.poly_order },
          { .key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = shift_b.id }
        };
        int io_meta_shift_dg_len = sizeof(io_meta_shift_dg)/sizeof(io_meta_shift_dg[0]);
        int io_meta_shift_len[] = {app->io_meta_basic_len, io_meta_shift_dg_len};
        const struct gkyl_msgpack_map_elem* io_meta_shift[] = {app->io_meta_basic, io_meta_shift_dg};
        struct gkyl_msgpack_data *mt_shift = gkyl_msgpack_create_union(sizeof(io_meta_shift_len)/sizeof(int),
          io_meta_shift_len, io_meta_shift);

        int sz = gkyl_calc_strlen(fmt, app->name, vars[dir], edge[edi]);
        char fileNm[sz+1]; // ensures no buffer overflow
        sprintf(fileNm, fmt, app->name, vars[dir], edge[edi]);

        gkyl_grid_sub_array_write(&shear_grid_core, &shear_r, mt_shift, shift_dg, fileNm);

        gkyl_array_release(shift_dg);
        gkyl_msgpack_data_release(mt_shift);
      }
    }
  }
}

static void
gk_field_fem_projection_par_rho_iwl_2x(gkyl_gyrokinetic_app *app, struct gk_field *field,
  struct gkyl_array *arr_dg, struct gkyl_array *arr_fem)
{
  // Project a DG charge density onto the parallel FEM basis to make it
  // continuous along z using different BCs in the core and SOL.

  // Gather the DG array into a global (in z) array.
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, arr_dg, field->rho_c_global_dg);

  // Smooth the DG array.
  gkyl_fem_parproj_set_rhs(field->fem_parproj_rho_core, field->rho_c_global_dg, field->rho_c_global_dg);
  gkyl_fem_parproj_solve(field->fem_parproj_rho_core, field->phi_fem);
  gkyl_fem_parproj_set_rhs(field->fem_parproj_rho_sol, field->rho_c_global_dg, field->rho_c_global_dg);
  gkyl_fem_parproj_solve(field->fem_parproj_rho_sol, field->phi_fem);

  // Copy global, continuous FEM array to a local array.
  gkyl_array_copy_range_to_range(arr_fem, field->phi_fem, &app->local, &field->global_sub_range);
}

static void
gk_field_fem_projection_par_phi_iwl_2x(gkyl_gyrokinetic_app *app, struct gk_field *field,
  struct gkyl_array *arr_dg, struct gkyl_array *arr_fem)
{
  // Project a DG potential onto the parallel FEM basis to make it
  // continuous along z using different BCs in the core and SOL.

  // Gather the DG array into a global (in z) array.
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, arr_dg, field->rho_c_global_dg);

  // Smooth the the DG array.
  gkyl_fem_parproj_set_rhs(field->fem_parproj_phi_core, field->rho_c_global_dg, field->rho_c_global_dg);
  gkyl_fem_parproj_solve(field->fem_parproj_phi_core, field->phi_fem);
  gkyl_fem_parproj_set_rhs(field->fem_parproj_phi_sol, field->rho_c_global_dg, field->rho_c_global_dg);
  gkyl_fem_parproj_solve(field->fem_parproj_phi_sol, field->phi_fem);

  // Copy global, continuous FEM array to a local array.
  gkyl_array_copy_range_to_range(arr_fem, field->phi_fem, &app->local, &field->global_sub_range);
}

static void
gk_field_fem_projection_par_phi_iwl_3x(gkyl_gyrokinetic_app *app, struct gk_field *field,
  struct gkyl_array *arr_dg, struct gkyl_array *arr_fem)
{
  // Project a DG field onto the parallel FEM basis to make it
  // continuous along z (or to solve a Poisson equation in 1x),
  // using different BCs in the core and SOL.

  // Gather the DG array into a global (in z) array.
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, arr_dg, field->rho_c_global_dg);

  // Apply TS BC in the core lower parallel boundary, and
  // fill core upper parallel boundary ghost with skin boundary value.
  gkyl_array_copy_range_to_range(field->rho_c_global_dg, field->rho_c_global_dg,
    &app->global_lower_ghost_par_core, &app->global_upper_skin_par_core);
  gkyl_bc_twistshift_advance(field->bc_ts_lo, field->rho_c_global_dg, field->rho_c_global_dg);
  gkyl_bc_basic_gyrokinetic_advance(field->gfss_bc_op_core_up, field->bc_buffer, field->rho_c_global_dg);

  // Smooth the the DG array.
  gkyl_fem_parproj_set_rhs(field->fem_parproj_phi_core, field->rho_c_global_dg, field->rho_c_global_dg);
  gkyl_fem_parproj_solve(field->fem_parproj_phi_core, field->phi_fem);
  gkyl_fem_parproj_set_rhs(field->fem_parproj_phi_sol, field->rho_c_global_dg, field->rho_c_global_dg);
  gkyl_fem_parproj_solve(field->fem_parproj_phi_sol, field->phi_fem);

  // Copy global, continuous FEM array to a local array.
  gkyl_array_copy_range_to_range(arr_fem, field->phi_fem, &app->local, &field->global_sub_range);
}

static void
gk_field_2x3x_fill_fem_parproj_bias_lines(struct gkyl_gyrokinetic_app *app, struct gk_field *f, struct gkyl_poisson_bc *poisson_bcs)
{
  // Create a bias line list that includes the perpendicular BCs if
  // they are Dirichlet, and the bias lines from the input file.
  
  // Create a temporary biased line list.
  int par_dir = app->cdim-1; // Parallel direction index.
  int num_bias_line = 0;
  int bl_idx = 0;
  int num_bias_line_in = 0;
  if (f->info.bias_line_list)
    num_bias_line_in = f->info.bias_line_list->num_bias_line;

  size_t bl_sz = (2+num_bias_line_in) * sizeof(struct gkyl_poisson_bias_line);
  struct gkyl_poisson_bias_line *bias_lines_buff = gkyl_malloc(bl_sz);
  if ((app->gk_geom->geqdsk_sign_convention == 0) && (poisson_bcs->lo_type[0] == GKYL_POISSON_DIRICHLET)) {
    // psi increases towards SOL.
    struct gkyl_poisson_bias_line *bl;
    // Core x boundary at lower z boundary.
    bl = &(bias_lines_buff[bl_idx]);
    bl->perp_dirs[0] = 0;
    bl->perp_dirs[1] = par_dir;
    bl->perp_coords[0] = app->grid.lower[0];
    bl->perp_coords[1] = app->grid.lower[par_dir];
    bl->val = poisson_bcs->lo_value[0].v[0];
    num_bias_line++;
    bl_idx++;
    // Core x boundary at upper z boundary.
    bl = &(bias_lines_buff[bl_idx]);
    bl->perp_dirs[0] = 0;
    bl->perp_dirs[1] = par_dir;
    bl->perp_coords[0] = app->grid.lower[0];
    bl->perp_coords[1] = app->grid.upper[par_dir];
    bl->val = poisson_bcs->lo_value[0].v[0];
    num_bias_line++;
    bl_idx++;
  } 
  else if ((app->gk_geom->geqdsk_sign_convention != 0) && (poisson_bcs->up_type[0] == GKYL_POISSON_DIRICHLET)) {
    // psi increases towards core.
    struct gkyl_poisson_bias_line *bl;
    // Core x boundary at lower z boundary.
    bl = &(bias_lines_buff[bl_idx]);
    bl->perp_dirs[0] = 0;
    bl->perp_dirs[1] = par_dir;
    bl->perp_coords[0] = app->grid.upper[0];
    bl->perp_coords[1] = app->grid.lower[par_dir];
    bl->val = poisson_bcs->up_value[0].v[0];
    num_bias_line++;
    bl_idx++;
    // Core x boundary at upper z boundary.
    bl = &(bias_lines_buff[bl_idx]);
    bl->perp_dirs[0] = 0;
    bl->perp_dirs[1] = par_dir;
    bl->perp_coords[0] = app->grid.upper[0];
    bl->perp_coords[1] = app->grid.upper[par_dir];
    bl->val = poisson_bcs->up_value[0].v[0];
    num_bias_line++;
    bl_idx++;
  } 
  for (int i=0; i<num_bias_line_in; i++) {
    struct gkyl_poisson_bias_line *bl = &(bias_lines_buff[bl_idx]);
    struct gkyl_poisson_bias_line *bl_inp = &(f->info.bias_line_list->bl[i]);
    bl->perp_dirs[0]   = bl_inp->perp_dirs[0]  ;
    bl->perp_dirs[1]   = bl_inp->perp_dirs[1]  ;
    bl->perp_coords[0] = bl_inp->perp_coords[0];
    bl->perp_coords[1] = bl_inp->perp_coords[1];
    bl->val = bl_inp->val;
    num_bias_line++;
    bl_idx++;
  }
  
  // Copy temporary bias line list into app.
  bl_sz = num_bias_line * sizeof(struct gkyl_poisson_bias_line);
  f->fem_parproj_bias_line_list.num_bias_line = num_bias_line;
  f->fem_parproj_bias_line_list.bl = gkyl_malloc(bl_sz);
  memcpy(f->fem_parproj_bias_line_list.bl, bias_lines_buff, bl_sz);
  gkyl_free(bias_lines_buff);
}

static void
gk_field_2x3x_add_IWL_updaters(struct gkyl_gyrokinetic_app *app, struct gk_field *f, struct gkyl_poisson_bc *poisson_bcs)
{
  // Allocation ranges and updaters for IWL field solve.

  // Parallel smoother for the charge density.
  enum gkyl_fem_parproj_bc_type fem_parproj_bc_rho_core, fem_parproj_bc_rho_sol;
  enum gkyl_fem_parproj_bc_type fem_parproj_bc_phi_core, fem_parproj_bc_phi_sol;

  if (app->cdim == 2) {
    fem_parproj_bc_rho_core = GKYL_FEM_PARPROJ_PERIODIC;
    fem_parproj_bc_rho_sol  = GKYL_FEM_PARPROJ_NONE;
    fem_parproj_bc_phi_core = GKYL_FEM_PARPROJ_PERIODIC;
    fem_parproj_bc_phi_sol  = GKYL_FEM_PARPROJ_NONE;

    f->fem_projection_par_rho_func = gk_field_fem_projection_par_rho_iwl_2x;
    f->fem_projection_par_phi_func = gk_field_fem_projection_par_phi_iwl_2x;
  }
  else if (app->cdim == 3) {
    fem_parproj_bc_rho_core = GKYL_FEM_PARPROJ_DIRICHLET_GHOST;
    fem_parproj_bc_rho_sol  = GKYL_FEM_PARPROJ_NONE;
    fem_parproj_bc_phi_core = GKYL_FEM_PARPROJ_DIRICHLET_GHOST;
    fem_parproj_bc_phi_sol  = GKYL_FEM_PARPROJ_NONE;

    f->fem_projection_par_rho_func = gk_field_fem_projection_par;
    f->fem_projection_par_phi_func = gk_field_fem_projection_par_phi_iwl_3x;

    // Take the TS function from the parallel BC of the first species.
    int par_dir = app->cdim-1; // Parallel direction index.
    struct gk_species *gks = &app->species[0];
    const struct gkyl_gyrokinetic_bc *par_lower_bc;
    for (int i = 0; i < 2*app->cdim; i++) {
      if ( gks->info.bcs[i].dir == par_dir && gks->info.bcs[i].type == GKYL_BC_GK_SPECIES_IWL) {
        if (gks->info.bcs[i].edge == GKYL_LOWER_EDGE) {
          par_lower_bc = (const struct gkyl_gyrokinetic_bc *) &gks->info.bcs[i];
          break;
        }
      }
    }

    // TS BC updater for up to low TS for the lower edge. This sets ghost_L = T_LU(ghost_L).
    int ghost[] = {1, 1, 1};
    struct gkyl_bc_twistshift_inp T_LU_lo = {
      .bc_dir = par_dir,
      .shift_dir = 1, // y shift.
      .shear_dir = 0, // shift varies with x.
      .edge = GKYL_LOWER_EDGE,
      .cdim = app->cdim,
      .bcdir_ext_update_r = app->global_par_ext_core,
      .num_ghost = ghost, // one ghost per config direction
      .basis = app->basis,
      .grid = app->grid,
      .shift_func = par_lower_bc->aux_profile,
      .shift_func_ctx = par_lower_bc->aux_ctx,
      .use_gpu = app->use_gpu,
    };
    f->bc_ts_lo = gkyl_bc_twistshift_new(&T_LU_lo);

    long buff_sz = GKYL_MAX2(app->global_lower_ghost_par_sol.volume, app->global_lower_ghost_par_core.volume);
    f->bc_buffer = mkarr(app->use_gpu, app->basis.num_basis, buff_sz);

    f->gfss_bc_op_core_up = gkyl_bc_basic_gyrokinetic_new(par_dir, GKYL_UPPER_EDGE, GKYL_BC_GK_FIELD_BOUNDARY_VALUE,
      app->basis_on_dev, &app->global_upper_skin_par_core, &app->global_upper_ghost_par_core,
      app->basis.num_basis, app->cdim, app->use_gpu);

    // Write the discrete shift to file.
    gk_field_3x_write_twistshift(app, f);
  }

  // Parallel smoother for the charge density.
  f->fem_parproj_rho_core = gkyl_fem_parproj_new(&app->global_core, &app->grid, &app->basis,
    fem_parproj_bc_rho_core, 0, 0, 0, app->use_gpu);
  f->fem_parproj_rho_sol = gkyl_fem_parproj_new(&app->global_sol, &app->grid, &app->basis,
    fem_parproj_bc_rho_sol, 0, 0, 0, app->use_gpu);

  // Fill bias line list for fem_parproj_phi.
  gk_field_2x3x_fill_fem_parproj_bias_lines(app, f, poisson_bcs);
    
  // Parallel smoother for the potential.
  f->fem_parproj_phi_core = gkyl_fem_parproj_new(&app->global_core, &app->grid, &app->basis,
    fem_parproj_bc_phi_core, &f->fem_parproj_bias_line_list, 0, 0, app->use_gpu);
  f->fem_parproj_phi_sol = gkyl_fem_parproj_new(&app->global_sol, &app->grid, &app->basis,
    fem_parproj_bc_phi_sol, &f->fem_parproj_bias_line_list, 0, 0, app->use_gpu);
    
}

static void
gk_field_rhs_poisson_perp_2x3x(struct gkyl_gyrokinetic_app *app, struct gk_field *field)
{
  // Smooth the charge density along z.
  field->fem_projection_par_rho_func(app, field, field->rho_c, field->rho_c);

  // Solve the Poisson equation.
  gkyl_fem_poisson_perp_set_rhs(field->fem_poisson_perp, field->rho_c);
  gkyl_fem_poisson_perp_solve(field->fem_poisson_perp, field->phi_smooth);

  // Smooth the potential along z.
  field->fem_projection_par_phi_func(app, field, field->phi_smooth, field->phi_smooth);

  // Finish the Poisson solve with FLR effects.
  field->invert_flr(app, field, field->phi_smooth);
}

static void
gk_field_ohm_solve(struct gkyl_gyrokinetic_app *app, struct gk_field *field){
  struct timespec wst = gkyl_wall_clock();
  
  gkyl_fem_poisson_perp_update_lhs(field->fem_apardot_solver, field->lapWeightAmpere, field->dApartdtSlvr_kSq);
  gkyl_fem_poisson_perp_set_rhs(field->fem_apardot_solver, field->currentDensdot);
  gkyl_fem_poisson_perp_solve(field->fem_apardot_solver, field->apardot);

  field->invert_flr(app, field, field->apardot);

  app->stat.field_apar_solve_tm += gkyl_time_diff_now_sec(wst);
}

static void 
gk_field_ampere_solve_enabled(gkyl_gyrokinetic_app *app, struct gk_field *field){
  struct timespec wst = gkyl_wall_clock();

  gkyl_fem_poisson_perp_set_rhs(field->fem_apar_solver, field->currentDens);
  gkyl_fem_poisson_perp_solve(field->fem_apar_solver, field->apar);

  // Smooth along z.
  // gk_field_fem_projection_par(app, field, field->apar, field->apar_smooth_aux);
  // gkyl_array_copy_range(field->apar, field->apar_smooth_aux, &app->local_ext);

  field->invert_flr(app, field, field->apar);

  app->stat.field_apar_solve_tm += gkyl_time_diff_now_sec(wst);
}

static void 
gk_field_ampere_solve_none(gkyl_gyrokinetic_app *app, struct gk_field *field){
  // Do nothing.
}

static void
gk_field_em_rhs_enabled(gkyl_gyrokinetic_app *app, struct gk_field *field, const struct gkyl_array *f_in[],  struct gkyl_array *rhs_in[])
{
  gk_field_accumulate_current_dens_dot(app, field, rhs_in);
  gk_field_accumulate_ohms_kSq(app, field, f_in);
  gk_field_ohm_solve(app, field);
}

static void
gk_field_em_rhs_none(gkyl_gyrokinetic_app *app, struct gk_field *field, const struct gkyl_array *f_in[],  struct gkyl_array *rhs_in[])
{
  // Do nothing.
}

static void
gk_field_fem_release_2x3x(const gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  gkyl_array_release(f->rho_c);
  gkyl_array_release(f->rho_c_global_dg);
  gkyl_array_release(f->rho_c_global_smooth);
  gkyl_array_release(f->phi_fem);
  gkyl_array_release(f->phi_smooth);
  gkyl_array_release(f->apar);
  gkyl_array_release(f->apardot);
  gkyl_array_release(f->apar_curr);
  gkyl_array_release(f->apar1);
  gkyl_array_release(f->aparnew);

  if (f->is_em) {
    gkyl_array_release(f->apar_smooth_aux);
    gkyl_array_release(f->currentDens);
    gkyl_array_release(f->currentDensdot);
    gkyl_array_release(f->lapWeightAmpere);
    gkyl_array_release(f->dApartdtSlvr_kSq);
    gkyl_fem_parproj_release(f->fem_apar_parproj);
    gkyl_fem_poisson_perp_release(f->fem_apar_solver);
    gkyl_fem_poisson_perp_release(f->fem_apardot_solver);
    if (app->use_gpu) {
      gkyl_array_release(f->apar_host);
      gkyl_array_release(f->apardot_host);
    }
  }

  if (app->use_gpu) {
    gkyl_array_release(f->phi_host);
  }

  gkyl_array_release(f->epsilon);

  gkyl_fem_poisson_perp_release(f->fem_poisson_perp);
  if (f->is_dirichletvar) {
    gkyl_array_release(f->phi_bc);
  }
  
  gkyl_fem_parproj_release(f->fem_parproj);

  gkyl_array_integrate_release(f->calc_em_energy);

  // Release IWL updaters.
  if (f->gkfield_id == GKYL_GK_FIELD_ES_IWL || f->gkfield_id == GKYL_GK_FIELD_EM_IWL) {
    gkyl_free(f->fem_parproj_bias_line_list.bl);
    gkyl_fem_parproj_release(f->fem_parproj_rho_core);
    gkyl_fem_parproj_release(f->fem_parproj_phi_core);
    gkyl_fem_parproj_release(f->fem_parproj_rho_sol);
    gkyl_fem_parproj_release(f->fem_parproj_phi_sol);

    if (app->cdim == 3) {
      gkyl_bc_twistshift_release(f->bc_ts_lo);
      gkyl_bc_basic_gyrokinetic_release(f->gfss_bc_op_core_up);
      gkyl_array_release(f->bc_buffer);
    }
  }
  
  if (f->use_flr) {
    gk_field_flr_release(app, f);
  }
}

void
gk_field_fem_new_2x3x(struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  // Create global subrange we'll copy the field solver solution from (into local).
  gkyl_sub_range_intersect(&f->global_sub_range, &app->global, &app->local);

  // Allocate arrays for charge density.
  f->rho_c = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->rho_c_global_dg = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
  f->rho_c_global_smooth = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);

  // Allocate arrays for electrostatic potential.
  f->phi_fem = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
  f->phi_smooth = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->apar = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->apar_curr = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->apar1 = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->aparnew = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->apardot = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

  // Setup electromagnetic variables if needed.
  if (f->is_em) {
    f->apar_host = f->apar;
    f->apardot_host = f->apardot;
    if (app->use_gpu) {
      f->apar_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
      f->apardot_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
    }

    f->apar_smooth_aux = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    f->currentDens = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    f->currentDensdot = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    f->lapWeightAmpere = mkarr(app->use_gpu, (2*(app->cdim/3)+1)*app->basis.num_basis, app->local_ext.volume);
    f->dApartdtSlvr_kSq = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  }

  // Allocate phi_host for I/O.
  f->phi_host = f->phi_smooth;
  if (app->use_gpu) {
    f->phi_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
  }

  if (f->gkfield_id == GKYL_GK_FIELD_ADIABATIC) {
    f->accumulate_rhoc_func = gk_field_accumulate_rho_c_adiabatic;
  } else {
    f->accumulate_rhoc_func = gk_field_accumulate_rho_c_poisson;
  }

  double polarization_weight = 0.0;
  double polarization_bmag = f->info.polarization_bmag ? f->info.polarization_bmag : app->bmag_ref;
  // Linearized polarization density
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *s = &app->species[i];
    polarization_weight += s->info.polarization_density*s->info.mass/pow(polarization_bmag,2);
  }
  // Allocate array for the polarization weight times geometric coefficients.
  f->epsilon = mkarr(app->use_gpu, (2*(app->cdim/3)+1)*app->basis.num_basis, app->local_ext.volume);
  
  // Initialize the polarization weight.
  struct gkyl_array *Jgij[3] = {app->gk_geom->geo_int.gxxj, app->gk_geom->geo_int.gxyj, app->gk_geom->geo_int.gyyj};
  for (int i=0; i<app->cdim-2/app->cdim; i++) {
    gkyl_array_set_offset(f->epsilon, polarization_weight, Jgij[i], i*app->basis.num_basis);
  }

  if (f->is_em) {
    // Compute the weights in front of the laplacian operator (metric x 1/mu0)
    for (int i=0; i<app->cdim-2/app->cdim; i++) {
      gkyl_array_set_offset(f->lapWeightAmpere, 1.0/f->info.mu0, Jgij[i], i*app->basis.num_basis);
    }
    gkyl_array_set(f->apar_energy_fac, -0.5, f->lapWeightAmpere);
  }

  bool bc_is_np[GKYL_MAX_CDIM]; // Is the BC in this direction non-periodic?
  for (int d=0; d<app->cdim; ++d) bc_is_np[d] = true;
  for (int d=0; d<app->num_periodic_dir; ++d) {
    bc_is_np[app->periodic_dirs[d]] = false;
  }

  // Translate input file BCs into Poisson BCs.
  struct gkyl_poisson_bc poisson_bcs = { };
  for (int d=0; d<app->cdim-1; d++) {
    if (bc_is_np[d]) {
      struct gkyl_gyrokinetic_bc *bc_lo = gk_fetch_bc_with_dir_edge(f->info.poisson_bcs, 2*app->cdim, d, GKYL_LOWER_EDGE);
      if (bc_lo != 0) {
        poisson_bcs.lo_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(bc_lo->type);
        for (int i=0; i<3; i++) {
          poisson_bcs.lo_value[d].v[i] = bc_lo->value[i];
        }
      }

      struct gkyl_gyrokinetic_bc *bc_up = gk_fetch_bc_with_dir_edge(f->info.poisson_bcs, 2*app->cdim, d, GKYL_UPPER_EDGE);
      if (bc_up != 0) {
        poisson_bcs.up_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(bc_up->type);
        for (int i=0; i<3; i++) {
          poisson_bcs.up_value[d].v[i] = bc_up->value[i];
        }
      }
    } else {
      poisson_bcs.lo_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(GKYL_BC_GK_FIELD_PERIODIC);
      poisson_bcs.up_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(GKYL_BC_GK_FIELD_PERIODIC);
    }
  }

  // Initialize the Poisson solver.
  f->fem_poisson_perp = gkyl_fem_poisson_perp_new(&app->local, &app->grid, app->basis,
    &poisson_bcs, f->info.bias_line_list, f->epsilon, NULL, app->use_gpu);

  f->phi_bc = 0;
  f->is_dirichletvar = false;
  for (int i=0; i<2*app->cdim; i++) {
    f->is_dirichletvar = f->is_dirichletvar ||
                          (f->info.poisson_bcs[i].type == GKYL_BC_GK_FIELD_DIRICHLET_VARYING ||
                           f->info.poisson_bcs[i].type == GKYL_BC_GK_FIELD_DIRICHLET_VARYING);
  }

  if (f->is_dirichletvar) {
    // Project the spatially varying BC if the user specifies it.
    f->phi_bc = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
    struct gkyl_array *phi_bc_ho = mkarr(false, f->phi_bc->ncomp, f->phi_bc->size);

    for (int d=0; d<app->cdim; d++) {
      struct gkyl_gyrokinetic_bc *bc_lo = gk_fetch_bc_with_dir_edge(f->info.poisson_bcs, 2*app->cdim, d, GKYL_LOWER_EDGE);
      if (bc_lo != 0) {
        if (bc_lo->type == GKYL_BC_GK_FIELD_DIRICHLET_VARYING) {
          gkyl_eval_on_nodes *phibc_proj = gkyl_eval_on_nodes_new(&app->grid, &app->basis, 
            1, bc_lo->aux_profile, bc_lo->aux_ctx);
          gkyl_eval_on_nodes_advance(phibc_proj, 0.0, &app->local_lower_skin[d], phi_bc_ho);
          gkyl_eval_on_nodes_release(phibc_proj);
        }
      }
      struct gkyl_gyrokinetic_bc *bc_up = gk_fetch_bc_with_dir_edge(f->info.poisson_bcs, 2*app->cdim, d, GKYL_UPPER_EDGE);
      if (bc_up != 0) {
        if (bc_up->type == GKYL_BC_GK_FIELD_DIRICHLET_VARYING) {
          gkyl_eval_on_nodes *phibc_proj = gkyl_eval_on_nodes_new(&app->grid, &app->basis, 
            1, bc_up->aux_profile, bc_up->aux_ctx);
          gkyl_eval_on_nodes_advance(phibc_proj, 0.0, &app->local_lower_skin[d], phi_bc_ho);
          gkyl_eval_on_nodes_release(phibc_proj);
        }
      }
    }
    gkyl_array_copy(f->phi_bc, phi_bc_ho);
    gkyl_array_release(phi_bc_ho);
  }

  // Potential smoothing (in z) updater
  enum gkyl_fem_parproj_bc_type fem_parproj_bc = GKYL_FEM_PARPROJ_NONE;
  for (int d=0; d<app->num_periodic_dir; ++d)
    if (app->periodic_dirs[d] == app->cdim-1) fem_parproj_bc = GKYL_FEM_PARPROJ_PERIODIC;

  f->fem_parproj = gkyl_fem_parproj_new(&app->global, &app->grid, &app->basis,
    fem_parproj_bc, 0, 0, 0, app->use_gpu);

  f->fem_projection_par_rho_func = gk_field_fem_projection_par;
  f->fem_projection_par_phi_func = gk_field_fem_projection_par;

  // Setup EM solvers.
  f->ampere_solve = gk_field_ampere_solve_none;
  f->em_rhs_func = gk_field_em_rhs_none;
  if (f->is_em) {
    // Translate input file BCs into Ampere BCs.
    for (int d=0; d<app->cdim-1; d++) {
      if (bc_is_np[d]) {
        struct gkyl_gyrokinetic_bc *bc_lo = gk_fetch_bc_with_dir_edge(f->info.ampere_bcs, 2*app->cdim, d, GKYL_LOWER_EDGE);
        if (bc_lo != 0) {
          f->ampere_bcs.lo_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(bc_lo->type);
          for (int i=0; i<3; i++)
            f->ampere_bcs.lo_value[d].v[i] = bc_lo->value[i];
        }

        struct gkyl_gyrokinetic_bc *bc_up = gk_fetch_bc_with_dir_edge(f->info.ampere_bcs, 2*app->cdim, d, GKYL_UPPER_EDGE);
        if (bc_up != 0) {
          f->ampere_bcs.up_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(bc_up->type);
          for (int i=0; i<3; i++)
            f->ampere_bcs.up_value[d].v[i] = bc_up->value[i];
        }
      } else {
        f->ampere_bcs.lo_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(GKYL_BC_GK_FIELD_PERIODIC);
        f->ampere_bcs.up_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(GKYL_BC_GK_FIELD_PERIODIC);
      }
    }
    // FEM solver for Ampere equation (solved at IC only).
    f->fem_apar_solver = gkyl_fem_poisson_perp_new(&app->local, &app->grid, app->basis,
        &f->ampere_bcs, 0, f->lapWeightAmpere, NULL, app->use_gpu);
    // FEM solver for Ohm's law (evolves d(Apart)/dt).
    f->fem_apardot_solver = gkyl_fem_poisson_perp_new(&app->local, &app->grid, app->basis,
      &f->ampere_bcs, 0, f->lapWeightAmpere, f->dApartdtSlvr_kSq, app->use_gpu);

    // We have to advance the updater once to be able to update the kSq later.
    gkyl_fem_poisson_perp_set_rhs(f->fem_apardot_solver, f->currentDensdot);
    gkyl_fem_poisson_perp_solve(f->fem_apardot_solver, f->apardot);
    
    // FEM smoother for Aparallel.
    f->fem_parproj_ampere_bc = GKYL_FEM_PARPROJ_NONE;
    f->fem_apar_parproj = gkyl_fem_parproj_new(&app->global, &app->grid, 
       &app->basis, f->fem_parproj_ampere_bc, 0, 0, 0, app->use_gpu);

    f->ampere_solve = gk_field_ampere_solve_enabled;
    f->em_rhs_func = f->info.is_apar_static ? gk_field_em_rhs_none : gk_field_em_rhs_enabled;
  }

  // Updater for field energy calculation.
  gkyl_array_set(f->es_energy_fac, 0.5, f->epsilon);
  f->calc_em_energy = gkyl_array_integrate_new(&app->grid, &app->basis, 
    1, GKYL_ARRAY_INTEGRATE_OP_EPS_GRADPERP_SQ, app->use_gpu);

  // Create operator needed for FLR effects.
  f->use_flr = false;
  f->invert_flr = gk_field_invert_flr_none;
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *s = &app->species[i];
    if (s->info.flr.type) {
      f->use_flr = f->use_flr || s->info.flr.type;
    }
  }
  if (f->use_flr) {
    gk_field_flr_new(app, f);
  }

  // Twist-and-shift boundary condition for phi and skin surface from ghost to impose phi periodicity at z=-pi.
  if (f->gkfield_id == GKYL_GK_FIELD_ES_IWL || f->gkfield_id == GKYL_GK_FIELD_EM_IWL) {
    gk_field_2x3x_add_IWL_updaters(app, f, &poisson_bcs);
  }

  // Set the pointer to the function that computes phi.
  f->rhs_phi_func = gk_field_rhs_poisson_perp_2x3x;

  // Set pointer to function that releases memory.
  f->release_func = gk_field_fem_release_2x3x;
}
