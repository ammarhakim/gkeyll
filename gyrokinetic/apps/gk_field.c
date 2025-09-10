#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_bc_basic.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_position_map.h>
#include <gkyl_util.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_comm_io.h>

#include <assert.h>
#include <float.h>
#include <time.h>

static void
gk_field_invert_flr(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *phi)
{
  gkyl_deflated_fem_poisson_advance(field->flr_op, phi, phi, phi);
}

static void
gk_field_invert_flr_none(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *phi)
{
}

static void
eval_on_nodes_c2p_position_func(const double *xcomp, double *xphys, void *ctx)
{
  struct gkyl_position_map *gpm = ctx;
  gkyl_position_map_eval_mc2nu(gpm, xcomp, xphys);
}

static void
gk_field_calc_energy_dt_active(gkyl_gyrokinetic_app *app, const struct gk_field *field, double dt, double *energy_reduced)
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_array_integrate_advance(field->calc_em_energy, field->phi_smooth, 
    1.0/dt, field->es_energy_fac, &app->local, &app->local, energy_reduced);
  app->stat.phidot_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_field_calc_energy_dt_none(gkyl_gyrokinetic_app *app, const struct gk_field *field, double dt, double *energy_reduced)
{
}

void
gk_field_calc_energy_dt(gkyl_gyrokinetic_app *app, const struct gk_field *field, double dt, double *energy_reduced)
{
  field->calc_energy_dt_func(app, field, dt, energy_reduced);
}

static void
gk_field_add_TSBC_and_SSFG_updaters(struct gkyl_gyrokinetic_app *app, struct gk_field *f)
{
  // We take the first species to copy the function for the TS BC
  struct gk_species *gks = &app->species[0];
  // Get the z BC info from the first species in our app
  const struct gkyl_gyrokinetic_bcs *bcz = &gks->info.bcz;
  // Define the parallel direction index (handle 2x and 3x cases).
  int zdir = app->cdim - 1;

  // TSBC updaters
  int ghost[] = {1, 1, 1};
  if (app->cdim == 3) {
    //TS BC updater for up to low TS for the lower edge
    //this sets ghost_L = T_LU(ghost_L)
    struct gkyl_bc_twistshift_inp T_LU_lo = {
      .bc_dir = zdir,
      .shift_dir = 1, // y shift.
      .shear_dir = 0, // shift varies with x.
      .edge = GKYL_LOWER_EDGE,
      .cdim = app->cdim,
      .bcdir_ext_update_r = app->local_par_ext_core,
      .num_ghost = ghost, // one ghost per config direction
      .basis = app->basis,
      .grid = app->grid,
      .shift_func = bcz->lower.aux_profile,
      .shift_func_ctx = bcz->lower.aux_ctx,
      .use_gpu = app->use_gpu,
    };
    // Add the forward TS updater to f
    f->bc_T_LU_lo = gkyl_bc_twistshift_new(&T_LU_lo);
  }

  // Add the SSFG updater for lower and upper application.
  f->ssfg_lo = gkyl_skin_surf_from_ghost_new(zdir, GKYL_LOWER_EDGE,
    app->basis, &app->lower_skin_par_core, &app->lower_ghost_par_core, app->use_gpu);
}

static void
gk_field_enforce_zbc(const gkyl_gyrokinetic_app *app, const struct gk_field *field, struct gkyl_array *finout)
{
  // Apply the periodicity to fill the ghost cells
  int num_periodic_dir = 1; // we need only periodicity in z
  int zdir = app->cdim - 1;
  int periodic_dirs[] = {zdir};
  gkyl_comm_array_per_sync(app->comm, &app->local, &app->local_ext,
    num_periodic_dir, periodic_dirs, finout); 
  
  // // Update the lower z ghosts with twist-and-shift if we are in 3x2v
  if(app->cdim == 3) {
    gkyl_bc_twistshift_advance(field->bc_T_LU_lo, finout, finout);
  }

  // Synchronize the array between the MPI processes to erase inner ghosts modification (handle multi GPU case)
  gkyl_comm_array_sync(app->comm, &app->local, &app->local_ext, finout);

  // Force the lower skin surface value to match the ghost cell at the node position.
  gkyl_skin_surf_from_ghost_advance(field->ssfg_lo, finout);
}

static void
gk_field_enforce_zbc_none(const gkyl_gyrokinetic_app *app, const struct gk_field *field, struct gkyl_array *finout)
{
}

// initialize field object
struct gk_field* 
gk_field_new(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app)
{
  struct gk_field *f = gkyl_malloc(sizeof(struct gk_field));

  f->info = gk->field;

  f->gkfield_id = f->info.gkfield_id ? f->info.gkfield_id : GKYL_GK_FIELD_ES;

  f->calc_init_field = !f->info.zero_init_field;
  f->update_field = !f->info.is_static;
  // The combination update_field=true, calc_init_field=false is not allowed.
  assert(!(f->update_field && (!f->calc_init_field)));

  // allocate arrays for charge density
  f->rho_c = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->rho_c_global_dg = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
  f->rho_c_global_smooth = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);

  // allocate arrays for electrostatic potential
  // global phi (only used in 1x simulations)
  f->phi_fem = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
  // local phi (assuming domain decomposition is *only* in z right now)
  f->phi_smooth = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

  if (f->gkfield_id == GKYL_GK_FIELD_EM) {
    f->apar_fem = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    f->apardot_fem = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  }

  f->init_phi_pol = false;
  if (f->info.polarization_potential) {
    // Project the initial potential onto a p+1 tensor basis and compute the polarization
    // density to use use by species in calculating the initial ion density.
    f->init_phi_pol = true;
    struct gkyl_basis phi_pol_basis;
    gkyl_cart_modal_tensor(&phi_pol_basis, app->cdim, app->poly_order+1);

    f->phi_pol = mkarr(app->use_gpu, phi_pol_basis.num_basis, app->local_ext.volume);
    struct gkyl_array *phi_pol_ho = app->use_gpu? mkarr(false, f->phi_pol->ncomp, f->phi_pol->size)
                                                : gkyl_array_acquire(f->phi_pol);

    struct gkyl_eval_on_nodes *phi_pol_proj = gkyl_eval_on_nodes_inew( &(struct gkyl_eval_on_nodes_inp){
      .grid = &app->grid,
      .basis = &phi_pol_basis,
      .num_ret_vals = 1,
      .eval = f->info.polarization_potential,
      .ctx = f->info.polarization_potential_ctx,
      .c2p_func = eval_on_nodes_c2p_position_func,
      .c2p_func_ctx = app->position_map,
    });

    gkyl_eval_on_nodes_advance(phi_pol_proj, 0.0, &app->local, phi_pol_ho);
    gkyl_array_copy(f->phi_pol, phi_pol_ho);
    
    gkyl_eval_on_nodes_release(phi_pol_proj);
    gkyl_array_release(phi_pol_ho);
  }

  // Create global subrange we'll copy the field solver solution from (into local).
  int intersect = gkyl_sub_range_intersect(&f->global_sub_range, &app->global, &app->local);

  // detect if this process contains an edge in the z dimension by comparing the local and global indices
  // for applying bias plan at the extremal z values only.
  int ndim = app->grid.ndim;
  f->info.poisson_bcs.contains_lower_z_edge = f->global_sub_range.lower[ndim-1] == app->global.lower[ndim-1];
  f->info.poisson_bcs.contains_upper_z_edge = f->global_sub_range.upper[ndim-1] == app->global.upper[ndim-1];

  f->epsilon = 0;
  struct gkyl_array *epsilon_global = 0;
  f->kSq = 0;  // not currently used by fem_perp_poisson
  double polarization_weight = 0.0; 
  double es_energy_fac_1d_adiabatic = 0.0; 
  if (f->gkfield_id == GKYL_GK_FIELD_BOLTZMANN) {
    polarization_weight = 1.0; 
    f->ambi_pot = gkyl_ambi_bolt_potential_new(&app->grid, &app->basis, 
      f->info.electron_mass, f->info.electron_charge, f->info.electron_temp, app->use_gpu);
    // Sheath_vals contains both the density and potential sheath values.
    for (int j=0; j<app->cdim; ++j) {
      f->sheath_vals[2*j]   = mkarr(app->use_gpu, 2*app->basis.num_basis, app->local_ext.volume);
      f->sheath_vals[2*j+1] = mkarr(app->use_gpu, 2*app->basis.num_basis, app->local_ext.volume);
    }
  } else {

    // Allocate array for the polarization weight times geometric coefficients.
    f->epsilon = mkarr(app->use_gpu, (2*(app->cdim/3)+1)*app->basis.num_basis, app->local_ext.volume);

    double polarization_bmag = f->info.polarization_bmag ? f->info.polarization_bmag : app->bmag_ref;
    // Linearized polarization density
    for (int i=0; i<app->num_species; ++i) {
      struct gk_species *s = &app->species[i];
      polarization_weight += s->info.polarization_density*s->info.mass/pow(polarization_bmag,2);
    }
    if (app->cdim == 1) {
      // Need to set weight to kperpsq*polarizationWeight for use in potential smoothing.
      gkyl_array_copy(f->epsilon, app->gk_geom->jacobgeo);
      gkyl_array_scale(f->epsilon, polarization_weight);
      gkyl_array_scale(f->epsilon, f->info.kperpSq);

      if (f->gkfield_id == GKYL_GK_FIELD_ADIABATIC) {
        // Add the contribution from adiabatic electrons (in principle any
        // species can be adiabatic, which we can add suport for later).
        double n_s0 = f->info.electron_density;
        double q_s = f->info.electron_charge;
        double T_s = f->info.electron_temp;
        double quasineut_contr = q_s*n_s0*q_s/T_s;
        
        struct gkyl_array *epsilon_adiab = mkarr(app->use_gpu, f->epsilon->ncomp, f->epsilon->size);
        gkyl_array_copy(epsilon_adiab, app->gk_geom->jacobgeo);
        gkyl_array_scale(epsilon_adiab, quasineut_contr);
        gkyl_array_accumulate(f->epsilon, 1., epsilon_adiab);
        gkyl_array_release(epsilon_adiab);

        es_energy_fac_1d_adiabatic = 0.5*quasineut_contr; 
      }

      // Gather gather epsilon for (global) smoothing in z.
      epsilon_global = mkarr(app->use_gpu, f->epsilon->ncomp, app->global_ext.volume);
      gkyl_comm_array_allgather(app->comm, &app->local, &app->global, f->epsilon, epsilon_global);
    }
    else if (app->cdim > 1) {
      // Initialize the polarization weight.
      struct gkyl_array *Jgij[3] = {app->gk_geom->gxxj, app->gk_geom->gxyj, app->gk_geom->gyyj};
      for (int i=0; i<app->cdim-2/app->cdim; i++) {
        gkyl_array_set_offset(f->epsilon, polarization_weight, Jgij[i], i*app->basis.num_basis);
      }

      // Initialize the Poisson solver.
      f->deflated_fem_poisson = gkyl_deflated_fem_poisson_new(app->grid, app->basis_on_dev, app->basis,
        app->local, f->global_sub_range, f->epsilon, 0, f->info.poisson_bcs, f->info.bias_plane_list, app->use_gpu);
      f->fem_poisson = gkyl_fem_poisson_perp_new(&app->local, &app->grid, app->basis,
        &f->info.poisson_bcs, f->epsilon, NULL, app->use_gpu);

      f->phi_bc = 0;
      f->is_dirichletvar = false;
      for (int d=0; d<app->cdim; d++) f->is_dirichletvar = f->is_dirichletvar ||
        (f->info.poisson_bcs.lo_type[d] == GKYL_POISSON_DIRICHLET_VARYING ||
         f->info.poisson_bcs.up_type[d] == GKYL_POISSON_DIRICHLET_VARYING);
      if (f->is_dirichletvar) {
        // Project the spatially varying BC if the user specifies it.
        f->phi_bc = mkarr(app->use_gpu, app->basis.num_basis, app->global_ext.volume);
        struct gkyl_array *phi_bc_ho = mkarr(false, f->phi_bc->ncomp, f->phi_bc->size);

        gkyl_eval_on_nodes *phibc_proj = gkyl_eval_on_nodes_new(&app->grid, &app->basis, 
          1, f->info.poisson_bcs.bc_value_func, f->info.poisson_bcs.bc_value_func_ctx);
        gkyl_eval_on_nodes_advance(phibc_proj, 0.0, &app->global, phi_bc_ho);
        gkyl_array_copy(f->phi_bc, phi_bc_ho);

        gkyl_eval_on_nodes_release(phibc_proj);
        gkyl_array_release(phi_bc_ho);
      }
    }
  }

  // Potential smoothing (in z) updater
  if (f->gkfield_id == GKYL_GK_FIELD_ES_IWL) {
    enum gkyl_fem_parproj_bc_type fem_parproj_bc_core, fem_parproj_bc_sol;
    if (app->cdim == 2) {
      fem_parproj_bc_core = GKYL_FEM_PARPROJ_PERIODIC;
      fem_parproj_bc_sol = GKYL_FEM_PARPROJ_NONE;
    } else {
      fem_parproj_bc_core = GKYL_FEM_PARPROJ_NONE;
      fem_parproj_bc_sol = GKYL_FEM_PARPROJ_NONE;
    }

    f->fem_parproj_core = gkyl_fem_parproj_new(&app->global_core, &app->basis,
      fem_parproj_bc_core, 0, 0, app->use_gpu);
    f->fem_parproj_sol = gkyl_fem_parproj_new(&app->global_sol, &app->basis,
      fem_parproj_bc_sol, 0, 0, app->use_gpu);
  } 
  else {
    enum gkyl_fem_parproj_bc_type fem_parproj_bc = GKYL_FEM_PARPROJ_NONE;
    for (int d=0; d<app->num_periodic_dir; ++d)
      if (app->periodic_dirs[d] == app->cdim-1) fem_parproj_bc = GKYL_FEM_PARPROJ_PERIODIC;

    f->fem_parproj = gkyl_fem_parproj_new(&app->global, &app->basis,
      fem_parproj_bc, epsilon_global, 0, app->use_gpu);
  }

  if (epsilon_global)
    gkyl_array_release(epsilon_global);

  f->phi_host = f->phi_smooth;  
  if (app->use_gpu) {
    f->phi_host = mkarr(false, app->basis.num_basis, app->local_ext.volume);
    f->em_energy_red = gkyl_cu_malloc(sizeof(double[1]));
    f->em_energy_red_global = gkyl_cu_malloc(sizeof(double[1]));
  } else {
    f->em_energy_red = gkyl_malloc(sizeof(double[1]));
    f->em_energy_red_global = gkyl_malloc(sizeof(double[1]));
  }

  f->integ_energy = gkyl_dynvec_new(GKYL_DOUBLE, 1);
  f->is_first_energy_write_call = true;

  // Factors for ES energy. 
  f->es_energy_fac = mkarr(app->use_gpu, (2*(app->cdim/3)+1)*app->basis.num_basis, app->local_ext.volume);
  f->es_energy_fac_1d = 0.0;
  if (app->cdim==1) {
    if (f->gkfield_id == GKYL_GK_FIELD_BOLTZMANN)
      f->es_energy_fac_1d = polarization_weight;
    else
      f->es_energy_fac_1d = polarization_weight*f->info.kperpSq + es_energy_fac_1d_adiabatic;

    f->calc_em_energy = gkyl_array_integrate_new(&app->grid, &app->basis, 
      1, GKYL_ARRAY_INTEGRATE_OP_SQ, app->use_gpu);
  }
  else {
    if (f->gkfield_id != GKYL_GK_FIELD_BOLTZMANN)
      gkyl_array_set(f->es_energy_fac, 0.5, f->epsilon);

    f->calc_em_energy = gkyl_array_integrate_new(&app->grid, &app->basis, 
      1, GKYL_ARRAY_INTEGRATE_OP_EPS_GRADPERP_SQ, app->use_gpu);
  }

  f->calc_energy_dt_func = gk_field_calc_energy_dt_none;
  if (f->info.time_rate_diagnostics) {
    f->calc_energy_dt_func = gk_field_calc_energy_dt_active;
    if (app->use_gpu) {
      f->em_energy_red_new = gkyl_cu_malloc(sizeof(double[1]));
      f->em_energy_red_old = gkyl_cu_malloc(sizeof(double[1]));
      gkyl_cu_memset(f->em_energy_red_new, 0, sizeof(double[1]));
      gkyl_cu_memset(f->em_energy_red_old, 0, sizeof(double[1]));
    } else {
      f->em_energy_red_new = gkyl_malloc(sizeof(double[1]));
      f->em_energy_red_old = gkyl_malloc(sizeof(double[1]));
      memset(f->em_energy_red_new, 0, sizeof(double[1]));
      memset(f->em_energy_red_old, 0, sizeof(double[1]));
    }
    f->integ_energy_dot = gkyl_dynvec_new(GKYL_DOUBLE, 1);
    f->is_first_energy_dot_write_call = true;
  }

  // setup biased lower wall (same size as electrostatic potential), by default is 0.0
  f->phi_wall_lo = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->has_phi_wall_lo = false;
  f->phi_wall_lo_evolve = false;
  if (f->info.phi_wall_lo) {
    f->has_phi_wall_lo = true;
    if (f->info.phi_wall_lo_evolve)
      f->phi_wall_lo_evolve = f->info.phi_wall_lo_evolve;

    f->phi_wall_lo_host = f->phi_wall_lo;
    if (app->use_gpu) 
      f->phi_wall_lo_host = mkarr(false, f->phi_wall_lo->ncomp, f->phi_wall_lo->size);

    f->phi_wall_lo_proj = gkyl_eval_on_nodes_new(&app->grid, &app->basis, 
      1, f->info.phi_wall_lo, f->info.phi_wall_lo_ctx);

    // Compute phi_wall_lo at t = 0
    gkyl_eval_on_nodes_advance(f->phi_wall_lo_proj, 0.0, &app->local_ext, f->phi_wall_lo_host);
    if (app->use_gpu) // note: phi_wall_lo_host is same as phi_wall_lo when not on GPUs
      gkyl_array_copy(f->phi_wall_lo, f->phi_wall_lo_host);
  }

  // setup biased upper wall (same size as electrostatic potential), by default is 0.0
  f->phi_wall_up = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  f->has_phi_wall_up = false;
  f->phi_wall_up_evolve = false;
  if (f->info.phi_wall_up) {
    f->has_phi_wall_up = true;
    if (f->info.phi_wall_up_evolve)
      f->phi_wall_up_evolve = f->info.phi_wall_up_evolve;

    f->phi_wall_up_host = f->phi_wall_up;
    if (app->use_gpu) 
      f->phi_wall_up_host = mkarr(false, f->phi_wall_up->ncomp, f->phi_wall_up->size);

    f->phi_wall_up_proj = gkyl_eval_on_nodes_new(&app->grid, &app->basis,
      1, f->info.phi_wall_up, f->info.phi_wall_up_ctx);

    // Compute phi_wall_up at t = 0
    gkyl_eval_on_nodes_advance(f->phi_wall_up_proj, 0.0, &app->local_ext, f->phi_wall_up_host);
    if (app->use_gpu) // note: phi_wall_up_host is same as phi_wall_up when not on GPUs
      gkyl_array_copy(f->phi_wall_up, f->phi_wall_up_host);
  }

  // Create operator needed for FLR effects.
  f->use_flr = false;
  f->invert_flr = gk_field_invert_flr_none;
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *s = &app->species[i];
    if (s->info.flr.type)
      f->use_flr = f->use_flr || s->info.flr.type;
  }
  if (f->use_flr) {
    assert(app->cdim > 1);
    f->invert_flr = gk_field_invert_flr;

    double flr_weight = 0.0; 
    for (int i=0; i<app->num_species; ++i) {
      struct gk_species *s = &app->species[i];
      double gyroradius_bmag = s->info.flr.bmag ? s->info.flr.bmag : app->bmag_ref;
      flr_weight += s->info.flr.Tperp*s->info.mass/(pow(s->info.charge*gyroradius_bmag,2.0));
    }
    // Initialize the weight in the Laplacian operator.
    f->flr_rhoSq_sum = mkarr(app->use_gpu, (2*(app->cdim-1)-1)*app->basis.num_basis, app->local_ext.volume);
    gkyl_array_set_offset(f->flr_rhoSq_sum, flr_weight, app->gk_geom->gxxj, 0*app->basis.num_basis);
    if (app->cdim > 2) {
      gkyl_array_set_offset(f->flr_rhoSq_sum, flr_weight, app->gk_geom->gxyj, 1*app->basis.num_basis);
      gkyl_array_set_offset(f->flr_rhoSq_sum, flr_weight, app->gk_geom->gyyj, 2*app->basis.num_basis);
    }
    // Initialize the factor multiplying the field in the FLR operator.
    f->flr_kSq = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    gkyl_array_shiftc(f->flr_kSq, -pow(sqrt(2.0),app->cdim), 0); // Sets kSq=-1.

    // If domain is not periodic use Dirichlet BCs.
    struct gkyl_poisson_bc flr_bc;
    for (int d=0; d<app->cdim-1; d++) {
      flr_bc.lo_type[d] = f->info.poisson_bcs.lo_type[d];
      flr_bc.lo_value[d] = f->info.poisson_bcs.lo_value[d];
      if (flr_bc.lo_type[d] != GKYL_POISSON_PERIODIC)
        flr_bc.lo_type[d] = GKYL_POISSON_DIRICHLET_VARYING;

      flr_bc.up_type[d] = f->info.poisson_bcs.up_type[d];
      flr_bc.up_value[d] = f->info.poisson_bcs.up_value[d];
      if (flr_bc.up_type[d] != GKYL_POISSON_PERIODIC)
        flr_bc.up_type[d] = GKYL_POISSON_DIRICHLET_VARYING;
    }
    // Deflated Poisson solve is performed on range assuming decomposition is *only* in z.
    f->flr_op = gkyl_deflated_fem_poisson_new(app->grid, app->basis_on_dev, app->basis,
      app->local, app->local, f->flr_rhoSq_sum, f->flr_kSq, flr_bc, NULL, app->use_gpu);
  }

    // twist-and-shift boundary condition for phi and skin surface from ghost to impose phi periodicity at z=-pi
  if (f->gkfield_id == GKYL_GK_FIELD_ES_IWL){
    gk_field_add_TSBC_and_SSFG_updaters(app,f);
    f->enforce_zbc = gk_field_enforce_zbc;
  } else {
    f->enforce_zbc = gk_field_enforce_zbc_none;
  }
  
  return f;
}

void
gk_field_calc_phi_wall(gkyl_gyrokinetic_app *app, struct gk_field *field, double tm)
{
  if (field->has_phi_wall_lo && field->phi_wall_lo_evolve) {
    gkyl_eval_on_nodes_advance(field->phi_wall_lo_proj, tm, &app->local_ext, field->phi_wall_lo_host);
    if (app->use_gpu) // note: phi_wall_lo_host is same as phi_wall_lo when not on GPUs
      gkyl_array_copy(field->phi_wall_lo, field->phi_wall_lo_host);
  }
  if (field->has_phi_wall_up && field->phi_wall_up_evolve) {
    gkyl_eval_on_nodes_advance(field->phi_wall_up_proj, tm, &app->local_ext, field->phi_wall_up_host);
    if (app->use_gpu) // note: phi_wall_up_host is same as phi_wall_up when not on GPUs
      gkyl_array_copy(field->phi_wall_up, field->phi_wall_up_host);
  }
}

void
gk_field_accumulate_rho_c(gkyl_gyrokinetic_app *app, struct gk_field *field, 
  const struct gkyl_array *fin[])
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_array_clear(field->rho_c, 0.0);
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *s = &app->species[i];

    gk_species_moment_calc(&s->m0, s->local, app->local, fin[i]);
    if (field->gkfield_id == GKYL_GK_FIELD_BOLTZMANN) {
      // For Boltzmann electrons, we only need ion density (and the ion density
      // times the conf-space Jacobian), not charge density.
      // Rescale moment by inverse of Jacobian.
      gkyl_dg_div_op_range(s->m0.mem_geo, app->basis, 0, field->rho_c, 0, s->m0.marr, 0, 
        app->gk_geom->jacobgeo, &app->local);  

      // We also need the M0 flux of the boundary flux through the z
      // boundaries. Put it in the ghost cells of f and take its moment.
      gk_species_bflux_get_flux(&s->bflux, app->cdim-1, GKYL_LOWER_EDGE, s->f1, &s->lower_ghost[app->cdim-1]);
      gk_species_moment_calc(&s->m0, s->lower_ghost[app->cdim-1], app->lower_ghost[app->cdim-1], s->f1);

      gk_species_bflux_get_flux(&s->bflux, app->cdim-1, GKYL_UPPER_EDGE, s->f1, &s->upper_ghost[app->cdim-1]);
      gk_species_moment_calc(&s->m0, s->upper_ghost[app->cdim-1], app->upper_ghost[app->cdim-1], s->f1);
    } else {
      // Gyroaverage the density if needed.
      s->gyroaverage(app, s, s->m0.marr, s->m0_gyroavg);

      gkyl_array_accumulate_range(field->rho_c, s->info.charge, s->m0_gyroavg, &app->local);
      if (field->gkfield_id == GKYL_GK_FIELD_ADIABATIC) {
        // Add the background (electron) charge density.
        double n_s0 = field->info.electron_density;
        double q_s = field->info.electron_charge;
        gkyl_array_shiftc_range(field->rho_c, q_s*n_s0*sqrt(2.0), 0, &app->local);
      }
    }
  } 
  app->stat.field_phi_rhs_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_field_calc_ambi_pot_sheath_vals(gkyl_gyrokinetic_app *app, struct gk_field *field)
{
  // Note that the M0 moment of boundary fluxes along z should
  // be stored in the ghost cells of m0.marr at this point.
  int idx_par = app->cdim-1;
  int off = 2*idx_par;

  int comm_sz;
  gkyl_comm_get_size(app->comm, &comm_sz);

  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *s = &app->species[i];

    // Assumes symmetric sheath BCs for now only in 1D
    // NOTE: this relies on the accumulate_rho_c calling gk_species_moment_calc(s->m0)
    // to calculate the particle flux and place it in the ghost cells of s->m0.marr.
    gkyl_ambi_bolt_potential_sheath_calc(field->ambi_pot, GKYL_LOWER_EDGE, 
      &app->lower_skin[idx_par], &app->lower_ghost[idx_par], app->gk_geom->cmag, 
      app->gk_geom->jacobtot_inv, s->m0.marr, field->rho_c, s->m0.marr, field->sheath_vals[off]);
    gkyl_ambi_bolt_potential_sheath_calc(field->ambi_pot, GKYL_UPPER_EDGE, 
      &app->upper_skin[idx_par], &app->upper_ghost[idx_par], app->gk_geom->cmag,
      app->gk_geom->jacobtot_inv, s->m0.marr, field->rho_c, s->m0.marr, field->sheath_vals[off+1]);

    // Broadcast the sheath values from skin processes to other processes.
    gkyl_comm_array_bcast(app->comm, field->sheath_vals[off]  , field->sheath_vals[off], 0);
    gkyl_comm_array_bcast(app->comm, field->sheath_vals[off+1], field->sheath_vals[off+1], comm_sz-1);

    // Copy upper sheath values into lower ghost & add to lower sheath values for averaging.
    gkyl_array_copy_range_to_range(field->sheath_vals[off+1], field->sheath_vals[off+1],
      &app->lower_ghost[idx_par], &app->upper_ghost[idx_par]);
    gkyl_array_accumulate(field->sheath_vals[off], 1., field->sheath_vals[off+1]);
    gkyl_array_scale(field->sheath_vals[off], 0.5);
  } 
}

static void
gk_field_fem_projection_par(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *arr_dg, struct gkyl_array *arr_fem)
{
  // Project a DG field onto the parallel FEM basis to make it
  // continuous along z (or to solve a Poisson equation in 1x).

  // Gather the DG array into a global (in z) array.
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, arr_dg, field->rho_c_global_dg);

  // Smooth the the DG array.
  gkyl_fem_parproj_set_rhs(field->fem_parproj, field->rho_c_global_dg, field->rho_c_global_dg);
  gkyl_fem_parproj_solve(field->fem_parproj, field->phi_fem);

  // Copy global, continuous FEM array to a local array.
  gkyl_array_copy_range_to_range(arr_fem, field->phi_fem, &app->local, &field->global_sub_range);
}

void
gk_field_rhs(gkyl_gyrokinetic_app *app, struct gk_field *field)
{
  // Compute the electrostatic potential.
  struct timespec wst = gkyl_wall_clock();
  if (field->gkfield_id == GKYL_GK_FIELD_BOLTZMANN) { 
    // Compute sheath density n_i,s and potential phi_s = (Te/e)*ln(n_i,s*v_te/(sqrt(2*pi)*Gamma_i)).
    gk_field_calc_ambi_pot_sheath_vals(app, app->field);

    // Solve phi = phi_s + (Te/e)*ln(n_i/n_i,s).
    gkyl_ambi_bolt_potential_phi_calc(field->ambi_pot, &app->local, &app->local_ext,
      field->rho_c, field->sheath_vals[2*(app->cdim-1)], field->phi_smooth);

    // Smooth the potential along z.
    gk_field_fem_projection_par(app, field, field->phi_smooth, field->phi_smooth);
  }
  else {
    if (app->cdim == 1) {
      // Solve the Poisson equation in 1x with the parallel FEM projection.
      gk_field_fem_projection_par(app, field, field->rho_c, field->phi_smooth);
    }
    else if (app->cdim > 1) {
      // Gather charge density into global array.
      // Smooth the charge density. Input is rho_c_global_dg, globally smoothed in z,
      // and then output should be in *local* phi_smooth.
      if (field->gkfield_id == GKYL_GK_FIELD_ES_IWL) {
        gkyl_comm_array_allgather(app->comm, &app->local, &app->global, field->rho_c, field->rho_c_global_dg);
        gkyl_fem_parproj_set_rhs(field->fem_parproj_core, field->rho_c_global_dg, field->rho_c_global_dg);
        gkyl_fem_parproj_solve(field->fem_parproj_core, field->rho_c_global_smooth);
        gkyl_fem_parproj_set_rhs(field->fem_parproj_sol, field->rho_c_global_dg, field->rho_c_global_dg);
        gkyl_fem_parproj_solve(field->fem_parproj_sol, field->rho_c_global_smooth);
        gkyl_deflated_fem_poisson_advance(field->deflated_fem_poisson, field->rho_c_global_smooth,
          field->phi_bc, field->phi_smooth);
      }
      else {
//        // This workflow solves Poisson on planes (doesn't conserve energy).
//        gkyl_comm_array_allgather(app->comm, &app->local, &app->global, field->rho_c, field->rho_c_global_dg);
//        gkyl_fem_parproj_set_rhs(field->fem_parproj, field->rho_c_global_dg, field->rho_c_global_dg);
//        gkyl_fem_parproj_solve(field->fem_parproj, field->rho_c_global_smooth);
//        gkyl_deflated_fem_poisson_advance(field->deflated_fem_poisson, field->rho_c_global_smooth,
//          field->phi_bc, field->phi_smooth);

        // Smooth the charge density along z.
        gk_field_fem_projection_par(app, field, field->rho_c, field->rho_c);

        // Solve the Poisson equation.
        gkyl_fem_poisson_perp_set_rhs(field->fem_poisson, field->rho_c);
        gkyl_fem_poisson_perp_solve(field->fem_poisson, field->phi_smooth);

        // Smooth the potential along z.
        gk_field_fem_projection_par(app, field, field->phi_smooth, field->phi_smooth);
      }

      // Finish the Poisson solve with FLR effects.
      field->invert_flr(app, field, field->phi_smooth);

      // Enforce a BC of the field in the parallel direction.
      field->enforce_zbc(app, field, field->phi_smooth);

    }
  }
  app->stat.field_phi_solve_tm += gkyl_time_diff_now_sec(wst);
}

static struct gkyl_app_restart_status
header_from_file(gkyl_gyrokinetic_app *app, const char *fname)
{
  struct gkyl_app_restart_status rstat = { .io_status = 2 };
  
  FILE *fp = 0;
  with_file(fp, fname, "r") {
    struct gkyl_rect_grid grid;
    struct gkyl_array_header_info hdr;
    rstat.io_status = gkyl_grid_sub_array_header_read_fp(&grid, &hdr, fp);

    if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
      if (hdr.etype != GKYL_DOUBLE)
        rstat.io_status = GKYL_ARRAY_RIO_DATA_MISMATCH;
    }

    struct gyrokinetic_output_meta meta =
      gk_meta_from_mpack( &(struct gkyl_msgpack_data) {
          .meta = hdr.meta,
          .meta_sz = hdr.meta_size
        }
      );

    rstat.frame = meta.frame;
    rstat.stime = meta.stime;

    gkyl_grid_sub_array_header_release(&hdr);
  }
  
  return rstat;
}

void
gk_field_file_import_init(struct gkyl_gyrokinetic_app *app, struct gkyl_gyrokinetic_ic_import inp)
{
  // Import initial condition from a file.
  struct gkyl_app_restart_status rstat = header_from_file(app, inp.file_name);

  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    struct gkyl_app_restart_status rstat;
    rstat.io_status = gkyl_comm_array_read(app->comm, &app->grid, &app->local, app->field->phi_host, inp.file_name);
    gkyl_array_copy(app->field->phi_smooth, app->field->phi_host);
  }
  else
    assert(false);
}

void
gk_field_project_init(struct gkyl_gyrokinetic_app *app)
{
  // Project the initial field.
  struct gkyl_eval_on_nodes *phi_proj = gkyl_eval_on_nodes_new(&app->grid, &app->basis,
    1, app->field->info.init_field_profile, app->field->info.init_field_profile_ctx);
  gkyl_eval_on_nodes_advance(phi_proj, 0.0, &app->local, app->field->phi_host);
  gkyl_eval_on_nodes_release(phi_proj);
  gkyl_array_copy(app->field->phi_smooth, app->field->phi_host);
}

void
gk_field_calc_energy(gkyl_gyrokinetic_app *app, double tm, const struct gk_field *field)
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_array_integrate_advance(field->calc_em_energy, field->phi_smooth, 
    1.0, field->es_energy_fac, &app->local, &app->local, field->em_energy_red);

  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, field->em_energy_red, field->em_energy_red_global);

  double energy_global[1] = { 0.0 };
  if (app->use_gpu)
    gkyl_cu_memcpy(energy_global, field->em_energy_red_global, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
  else
    energy_global[0] = field->em_energy_red_global[0];
  
  if (app->cdim == 1)
    energy_global[0] *= field->es_energy_fac_1d; 

  gkyl_dynvec_append(field->integ_energy, tm, energy_global);

  if (field->info.time_rate_diagnostics) {
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, field->em_energy_red_old, field->em_energy_red_global);
    double energy_dot_global_old[1] = { 0.0 };
    if (app->use_gpu)
      gkyl_cu_memcpy(energy_dot_global_old, field->em_energy_red_global, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
    else
      energy_dot_global_old[0] = field->em_energy_red_global[0];
    if (app->cdim == 1)
      energy_dot_global_old[0] *= field->es_energy_fac_1d; 

    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 1, field->em_energy_red_new, field->em_energy_red_global);
    double energy_dot_global_new[1] = { 0.0 };
    if (app->use_gpu)
      gkyl_cu_memcpy(energy_dot_global_new, field->em_energy_red_global, sizeof(double[1]), GKYL_CU_MEMCPY_D2H);
    else
      energy_dot_global_new[0] = field->em_energy_red_global[0];
    if (app->cdim == 1)
      energy_dot_global_new[0] *= field->es_energy_fac_1d; 

    double energy_dot_global[1] = { 0.0 };
    energy_dot_global[0] = energy_dot_global_new[0] - energy_dot_global_old[0];

    gkyl_dynvec_append(field->integ_energy_dot, tm, energy_dot_global);
  }
  app->stat.field_diag_calc_tm += gkyl_time_diff_now_sec(wst);
}

// release resources for field
void
gk_field_release(const gkyl_gyrokinetic_app* app, struct gk_field *f)
{
  gkyl_array_release(f->rho_c);
  gkyl_array_release(f->rho_c_global_dg);
  gkyl_array_release(f->rho_c_global_smooth);
  gkyl_array_release(f->phi_fem);
  gkyl_array_release(f->phi_smooth);

  if (f->gkfield_id == GKYL_GK_FIELD_EM) {
    gkyl_array_release(f->apar_fem);
    gkyl_array_release(f->apardot_fem);
  }

  if (f->init_phi_pol) {
    gkyl_array_release(f->phi_pol);
  }

  gkyl_array_release(f->es_energy_fac);
  if (f->gkfield_id == GKYL_GK_FIELD_BOLTZMANN) {
    gkyl_ambi_bolt_potential_release(f->ambi_pot);
    for (int i=0; i<2*app->cdim; ++i) 
      gkyl_array_release(f->sheath_vals[i]);
  } 
  else {
    gkyl_array_release(f->epsilon);
    if (app->cdim > 1) {
      gkyl_deflated_fem_poisson_release(f->deflated_fem_poisson);
      gkyl_fem_poisson_perp_release(f->fem_poisson);
      if (f->is_dirichletvar) {
        gkyl_array_release(f->phi_bc);
      }
    }
  }

  if (app->cdim > 1 && f->gkfield_id == GKYL_GK_FIELD_ES_IWL) {
    gkyl_fem_parproj_release(f->fem_parproj_core);
    gkyl_fem_parproj_release(f->fem_parproj_sol);
  } else {
    gkyl_fem_parproj_release(f->fem_parproj);
  }

  // Release TS BS and SSFG updater
  if (f->gkfield_id == GKYL_GK_FIELD_ES_IWL) {
    if(app->cdim == 3) {
      gkyl_bc_twistshift_release(f->bc_T_LU_lo);
    }
    gkyl_skin_surf_from_ghost_release(f->ssfg_lo);
  }

  gkyl_dynvec_release(f->integ_energy);
  gkyl_array_integrate_release(f->calc_em_energy);
  if (app->use_gpu) {
    gkyl_array_release(f->phi_host);
    gkyl_cu_free(f->em_energy_red);
    gkyl_cu_free(f->em_energy_red_global);
  } else {
    gkyl_free(f->em_energy_red);
    gkyl_free(f->em_energy_red_global);
  }

  if (f->info.time_rate_diagnostics) {
    gkyl_dynvec_release(f->integ_energy_dot);
    if (app->use_gpu) {
      gkyl_cu_free(f->em_energy_red_new);
      gkyl_cu_free(f->em_energy_red_old);
    } else {
      gkyl_free(f->em_energy_red_new);
      gkyl_free(f->em_energy_red_old);
    }
  }

  gkyl_array_release(f->phi_wall_lo);
  gkyl_array_release(f->phi_wall_up);
  if (f->has_phi_wall_lo) {
    gkyl_eval_on_nodes_release(f->phi_wall_lo_proj);
    if (app->use_gpu) 
      gkyl_array_release(f->phi_wall_lo_host);
  }
  if (f->has_phi_wall_up) {
    gkyl_eval_on_nodes_release(f->phi_wall_up_proj);
    if (app->use_gpu) 
      gkyl_array_release(f->phi_wall_up_host);
  }

  if (f->use_flr) {
    gkyl_array_release(f->flr_rhoSq_sum);
    gkyl_array_release(f->flr_kSq);
    gkyl_deflated_fem_poisson_release(f->flr_op);
  }

  gkyl_free(f);
}

