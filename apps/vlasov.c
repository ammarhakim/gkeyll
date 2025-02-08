#include <stdarg.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dflt.h>
#include <gkyl_dynvec.h>
#include <gkyl_null_comm.h>

#include <gkyl_vlasov_priv.h>

gkyl_vlasov_app*
gkyl_vlasov_app_new(struct gkyl_vm *vm)
{
  disable_denorm_float();

  assert(vm->num_species <= GKYL_MAX_SPECIES);

  gkyl_vlasov_app *app = gkyl_malloc(sizeof(gkyl_vlasov_app));

  int cdim = app->cdim = vm->cdim;
  int vdim = app->vdim = vm->vdim;
  int pdim = cdim+vdim;
  int poly_order = app->poly_order = vm->poly_order;
  int ns = app->num_species = vm->num_species;
  int nsf = app->num_fluid_species = vm->num_fluid_species;

  double cfl_frac = vm->cfl_frac == 0 ? 1.0 : vm->cfl_frac;
  app->cfl = cfl_frac;

#ifdef GKYL_HAVE_CUDA
  app->use_gpu = vm->use_gpu;
#else
  app->use_gpu = false; // can't use GPUs if we don't have them!
#endif

  app->num_periodic_dir = vm->num_periodic_dir;
  for (int d=0; d<cdim; ++d)
    app->periodic_dirs[d] = vm->periodic_dirs[d];

  strcpy(app->name, vm->name);
  app->tcurr = 0.0; // reset on init

  if (app->use_gpu) {
    // allocate device basis if we are using GPUs
    app->basis_on_dev.basis = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    app->basis_on_dev.confBasis = gkyl_cu_malloc(sizeof(struct gkyl_basis));
  }
  else {
    app->basis_on_dev.basis = &app->basis;
    app->basis_on_dev.confBasis = &app->confBasis;
  }

  // basis functions
  switch (vm->basis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      gkyl_cart_modal_serendip(&app->confBasis, cdim, poly_order);
      if (vdim > 0) { 
        if (poly_order == 1) {
          /* Force hybrid basis (p=2 in velocity space). */
          gkyl_cart_modal_hybrid(&app->basis, cdim, vdim);
          gkyl_cart_modal_serendip(&app->velBasis, vdim, 2);
        }
        else {
          gkyl_cart_modal_serendip(&app->basis, pdim, poly_order);
          gkyl_cart_modal_serendip(&app->velBasis, vdim, poly_order);
        }
      }
      if (app->use_gpu) {
        gkyl_cart_modal_serendip_cu_dev(app->basis_on_dev.confBasis, cdim, poly_order);
        if (vdim > 0) {
          if (poly_order == 1) {
            /* Force hybrid basis (p=2 in velocity space). */
            gkyl_cart_modal_hybrid_cu_dev(app->basis_on_dev.basis, cdim, vdim); 
          }
          else {
            gkyl_cart_modal_serendip_cu_dev(app->basis_on_dev.basis, pdim, poly_order);
          }
        }
      }
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      gkyl_cart_modal_tensor(&app->confBasis, cdim, poly_order);
      if (vdim > 0) {
        gkyl_cart_modal_tensor(&app->basis, pdim, poly_order);
        gkyl_cart_modal_tensor(&app->velBasis, vdim, poly_order);
      }
      if (app->use_gpu) {
        gkyl_cart_modal_tensor_cu_dev(app->basis_on_dev.confBasis, cdim, poly_order);
        if (vdim > 0) {
          gkyl_cart_modal_tensor_cu_dev(app->basis_on_dev.basis, pdim, poly_order);
        }
      }
      break;

    default:
      assert(false);
      break;
  }

  gkyl_rect_grid_init(&app->grid, cdim, vm->lower, vm->upper, vm->cells);

  int ghost[] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&app->grid, ghost, &app->global_ext, &app->global);

  if (vm->has_low_inp) {
    // create local and local_ext from user-supplied local range
    gkyl_create_ranges(&vm->low_inp.local_range, ghost, &app->local_ext, &app->local);
    
    if (vm->low_inp.comm)
      app->comm = gkyl_comm_acquire(vm->low_inp.comm);
    else {
      int cuts[3] = { 1, 1, 1 };
      struct gkyl_rect_decomp *rect_decomp =
        gkyl_rect_decomp_new_from_cuts(cdim, cuts, &app->global);
      
      app->comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
          .decomp = rect_decomp,
          .use_gpu = app->use_gpu
        }
      );

      gkyl_rect_decomp_release(rect_decomp);
    }
  }
  else {
    // global and local ranges are same, and so just copy
    memcpy(&app->local, &app->global, sizeof(struct gkyl_range));
    memcpy(&app->local_ext, &app->global_ext, sizeof(struct gkyl_range));

    int cuts[3] = { 1, 1, 1 };
    struct gkyl_rect_decomp *rect_decomp =
      gkyl_rect_decomp_new_from_cuts(cdim, cuts, &app->global);
    
    app->comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .decomp = rect_decomp,
        .use_gpu = app->use_gpu
      }
    );
    
    gkyl_rect_decomp_release(rect_decomp);
  }
  // local skin and ghost ranges for configuration space fields
  for (int dir=0; dir<cdim; ++dir) {
    gkyl_skin_ghost_ranges(&app->lower_skin[dir], &app->lower_ghost[dir], dir, GKYL_LOWER_EDGE, &app->local_ext, ghost); 
    gkyl_skin_ghost_ranges(&app->upper_skin[dir], &app->upper_ghost[dir], dir, GKYL_UPPER_EDGE, &app->local_ext, ghost);
  }

  // Configuration space geometry initialization
  // Note: *only* uses a p=1 DG representation of the geometry (JJ: 11/24/23)
  app->c2p_ctx = app->mapc2p = 0;  
  app->has_mapc2p = vm->mapc2p ? true : false;

  if (app->has_mapc2p) {
    // initialize computational to physical space mapping
    app->c2p_ctx = vm->c2p_ctx;
    app->mapc2p = vm->mapc2p;

    // we project mapc2p on p=1 basis functions
    struct gkyl_basis basis;
    gkyl_cart_modal_tensor(&basis, cdim, 1);

    // initialize DG field representing mapping
    struct gkyl_array *c2p = mkarr(false, cdim*basis.num_basis, app->local_ext.volume);
    gkyl_eval_on_nodes *ev_c2p = gkyl_eval_on_nodes_new(&app->grid, &basis, cdim, vm->mapc2p, vm->c2p_ctx);
    gkyl_eval_on_nodes_advance(ev_c2p, 0.0, &app->local_ext, c2p);

    // write DG projection of mapc2p to file
    cstr fileNm = cstr_from_fmt("%s-mapc2p.gkyl", app->name);
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, c2p, fileNm.str);
    cstr_drop(&fileNm);

    gkyl_array_release(c2p);
    gkyl_eval_on_nodes_release(ev_c2p);
  }

  // create geometry object
  app->geom = gkyl_wave_geom_new(&app->grid, &app->local_ext,
    app->mapc2p, app->c2p_ctx, app->use_gpu);

  app->has_field = !vm->skip_field; // note inversion of truth value
  if (app->has_field)
    app->field = vm_field_new(vm, app);

  // allocate space to store species objects
  app->species = ns>0 ? gkyl_malloc(sizeof(struct vm_species[ns])) : 0;
  for (int i=0; i<ns; ++i)
    app->species[i] = (struct vm_species) { };

  // set info for each species: this needs to be done here as we need
  // to access species name from vm_species_init
  for (int i=0; i<ns; ++i)
    app->species[i].info = vm->species[i];

  // allocate space to store fluid species objects
  app->fluid_species = nsf>0 ? gkyl_malloc(sizeof(struct vm_fluid_species[nsf])) : 0;

  // set info for each fluid species: this needs to be done here as we
  // need to access species name from vm_fluid_species_init
  for (int i=0; i<nsf; ++i)
    app->fluid_species[i].info = vm->fluid_species[i];

  // initialize each species
  for (int i=0; i<ns; ++i) 
    vm_species_init(vm, app, &app->species[i]);

  // initialize each species cross-species terms: this has to be done here
  // as need pointers to colliding species' collision objects
  // allocated in the previous step
  for (int i=0; i<ns; ++i)
    if (app->species[i].collision_id == GKYL_LBO_COLLISIONS
      && app->species[i].lbo.num_cross_collisions) {
      vm_species_lbo_cross_init(app, &app->species[i], &app->species[i].lbo);
    }

  // initialize each species source terms: this has to be done here
  // as they may initialize a bflux updater for their source species
  for (int i=0; i<ns; ++i)
    if (app->species[i].source_id)
      vm_species_source_init(app, &app->species[i], &app->species[i].src);

  // initialize each fluid species
  // Fluid species must be initialized after kinetic species, as some fluid species couple
  // to kinetic species and pointers are allocated by the kinetic species objects
  for (int i=0; i<nsf; ++i)
    vm_fluid_species_init(vm, app, &app->fluid_species[i]);

  for (int i=0; i<nsf; ++i)
    if (app->fluid_species[i].source_id)
      vm_fluid_species_source_init(app, &app->fluid_species[i], &app->fluid_species[i].src);

  // Check if there are both any fluid species and an EM field. 
  // If there are, initialize the implicit fluid-EM coupling solver.
  app->has_fluid_em_coupling = false;
  if (nsf > 0 && app->has_field) {
    app->has_fluid_em_coupling = true;
    app->fl_em = vm_fluid_em_coupling_init(app);
  }

  // Use implicit BGK collisions if specified
  app->has_implicit_coll_scheme = false;
  for (int i=0; i<ns; ++i){
    if (vm->species[i].collisions.has_implicit_coll_scheme){
      app->has_implicit_coll_scheme = true;
    }
  }

  // Set the appropriate update function for taking a single time step
  // If we have implicit fluid-EM coupling or implicit BGK collisions, 
  // we perform a first-order operator split and treat those terms implicitly.
  // Otherwise, we default to an SSP-RK3 method. 
  if (app->has_implicit_coll_scheme || app->has_fluid_em_coupling) {
    app->update_func = vlasov_update_op_split;
  }
  else {
    app->update_func = vlasov_update_ssp_rk3;
  }

  // initialize stat object
  app->stat = (struct gkyl_vlasov_stat) {
    .use_gpu = app->use_gpu,
    .stage_2_dt_diff = { DBL_MAX, 0.0 },
    .stage_3_dt_diff = { DBL_MAX, 0.0 },
  };

  return app;
}

struct vm_species *
vm_find_species(const gkyl_vlasov_app *app, const char *nm)
{
  for (int i=0; i<app->num_species; ++i)
    if (strcmp(nm, app->species[i].info.name) == 0)
      return &app->species[i];
  return 0;
}

int
vm_find_species_idx(const gkyl_vlasov_app *app, const char *nm)
{
  for (int i=0; i<app->num_species; ++i)
    if (strcmp(nm, app->species[i].info.name) == 0)
      return i;
  return -1;
}

struct vm_fluid_species *
vm_find_fluid_species(const gkyl_vlasov_app *app, const char *nm)
{
  for (int i=0; i<app->num_fluid_species; ++i)
    if (strcmp(nm, app->fluid_species[i].info.name) == 0)
      return &app->fluid_species[i];
  return 0;
}

int
vm_find_fluid_species_idx(const gkyl_vlasov_app *app, const char *nm)
{
  for (int i=0; i<app->num_fluid_species; ++i)
    if (strcmp(nm, app->fluid_species[i].info.name) == 0)
      return i;
  return -1;
}

void
gkyl_vlasov_app_apply_ic(gkyl_vlasov_app* app, double t0)
{
  app->tcurr = t0;
  if (app->has_field) 
    gkyl_vlasov_app_apply_ic_field(app, t0);
  for (int i=0; i<app->num_species; ++i)
    gkyl_vlasov_app_apply_ic_species(app, i, t0);
  for (int i=0; i<app->num_fluid_species; ++i)
    gkyl_vlasov_app_apply_ic_fluid_species(app, i, t0);
}

void
gkyl_vlasov_app_apply_ic_field(gkyl_vlasov_app* app, double t0)
{
  app->tcurr = t0;

  struct timespec wtm = gkyl_wall_clock();
  vm_field_apply_ic(app, app->field, t0);
  app->stat.init_field_tm += gkyl_time_diff_now_sec(wtm);

  vm_field_apply_bc(app, app->field, app->field->em);
}

void
gkyl_vlasov_app_apply_ic_species(gkyl_vlasov_app* app, int sidx, double t0)
{
  assert(sidx < app->num_species);

  app->tcurr = t0;
  struct timespec wtm = gkyl_wall_clock();
  vm_species_apply_ic(app, &app->species[sidx], t0);
  app->stat.init_species_tm += gkyl_time_diff_now_sec(wtm);

  vm_species_apply_bc(app, &app->species[sidx], app->species[sidx].f);
}

void
gkyl_vlasov_app_apply_ic_fluid_species(gkyl_vlasov_app* app, int sidx, double t0)
{
  assert(sidx < app->num_fluid_species);

  app->tcurr = t0;
  struct timespec wtm = gkyl_wall_clock();
  vm_fluid_species_apply_ic(app, &app->fluid_species[sidx], t0);
  app->stat.init_fluid_species_tm += gkyl_time_diff_now_sec(wtm);

  vm_fluid_species_apply_bc(app, &app->fluid_species[sidx], app->fluid_species[sidx].fluid);
}

void
gkyl_vlasov_app_calc_mom(gkyl_vlasov_app* app)
{
  for (int i=0; i<app->num_species; ++i) {
    struct vm_species *vm_s = &app->species[i];

    for (int m=0; m<vm_s->info.num_diag_moments; ++m) {
      struct timespec wst = gkyl_wall_clock();
      vm_species_moment_calc(&vm_s->moms[m], vm_s->local, app->local, vm_s->f);
      app->stat.mom_tm += gkyl_time_diff_now_sec(wst);
      app->stat.nmom += 1;
    }
  }
}

void
gkyl_vlasov_app_calc_integrated_mom(gkyl_vlasov_app* app, double tm)
{
  int vdim = app->vdim;
  double avals[2+vdim], avals_global[2+vdim];

  struct timespec wst = gkyl_wall_clock();

  for (int i=0; i<app->num_species; ++i) {
    struct vm_species *vm_s = &app->species[i];

    struct timespec wst = gkyl_wall_clock();

    vm_species_moment_calc(&vm_s->integ_moms, vm_s->local, app->local, vm_s->f);
    // reduce to compute sum over whole domain, append to diagnostics
    if (app->use_gpu) {
      gkyl_array_reduce_range(vm_s->red_integ_diag, vm_s->integ_moms.marr, GKYL_SUM, &app->local);
      gkyl_cu_memcpy(avals, vm_s->red_integ_diag, sizeof(double[2+vdim]), GKYL_CU_MEMCPY_D2H);
    }
    else {
      gkyl_array_reduce_range(avals, vm_s->integ_moms.marr_host, GKYL_SUM, &app->local);
    }
    gkyl_comm_all_reduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 2+vdim, avals, avals_global);
    gkyl_dynvec_append(vm_s->integ_diag, tm, avals_global);

    if (vm_s->source_id) {
      vm_species_moment_calc(&vm_s->src.integ_moms, vm_s->local, app->local, vm_s->src.source); 
      // reduce to compute sum over whole domain, append to diagnostics
      if (app->use_gpu) {
        gkyl_array_reduce_range(vm_s->src.red_integ_diag, vm_s->src.integ_moms.marr, GKYL_SUM, &app->local);
        gkyl_cu_memcpy(avals, vm_s->src.red_integ_diag, sizeof(double[2+vdim]), GKYL_CU_MEMCPY_D2H);
      }
      else {
        gkyl_array_reduce_range(avals, vm_s->integ_moms.marr_host, GKYL_SUM, &app->local);
      }
      gkyl_comm_all_reduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 2+vdim, avals, avals_global);
      gkyl_dynvec_append(vm_s->src.integ_diag, tm, avals_global);
    }

    app->stat.mom_tm += gkyl_time_diff_now_sec(wst);
    app->stat.nmom += 1;
  }

  double avals_fluid[6], avals_fluid_global[6];
  for (int i=0; i<app->num_fluid_species; ++i) {
    struct vm_fluid_species *f = &app->fluid_species[i];

    gkyl_array_clear(f->integ_mom, 0.0);

    vm_fluid_species_prim_vars(app, f, f->fluid);
    gkyl_dg_calc_fluid_integrated_vars(f->calc_fluid_vars, &app->local, 
      f->fluid, f->u, f->p, f->integ_mom);
    gkyl_array_scale_range(f->integ_mom, app->grid.cellVolume, &app->local);
    if (app->use_gpu) {
      gkyl_array_reduce_range(f->red_integ_diag, f->integ_mom, GKYL_SUM, &app->local);
      gkyl_cu_memcpy(avals_fluid, f->red_integ_diag, sizeof(double[6]), GKYL_CU_MEMCPY_D2H);
    }
    else { 
      gkyl_array_reduce_range(avals_fluid, f->integ_mom, GKYL_SUM, &app->local);
    }

    gkyl_comm_all_reduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 5, avals_fluid, avals_fluid_global);
    gkyl_dynvec_append(f->integ_diag, tm, avals_fluid_global);
  }

  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.ndiag += 1;
}

void
gkyl_vlasov_app_calc_integrated_L2_f(gkyl_vlasov_app* app, double tm)
{
  struct timespec wst = gkyl_wall_clock();
  for (int i=0; i<app->num_species; ++i) {
    struct vm_species *vm_s = &app->species[i];
    vm_species_calc_L2(app, tm, vm_s);
  }
  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.ndiag += 1;
}

void
gkyl_vlasov_app_calc_field_energy(gkyl_vlasov_app* app, double tm)
{
  struct timespec wst = gkyl_wall_clock();
  vm_field_calc_energy(app, tm, app->field);
  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.ndiag += 1;
}

void
gkyl_vlasov_app_write(gkyl_vlasov_app* app, double tm, int frame)
{
  app->stat.nio += 1;
  struct timespec wtm = gkyl_wall_clock();
  
  if (app->has_field)
    gkyl_vlasov_app_write_field(app, tm, frame);
  for (int i=0; i<app->num_species; ++i) {
    gkyl_vlasov_app_write_species(app, i, tm, frame);
    if(app->species[i].info.output_f_lte) {
      gkyl_vlasov_app_write_species_lte(app, i, tm, frame);
    }
    if ((app->species[i].collision_id == GKYL_FPO_COLLISIONS) && (app->species[i].fpo.write_diagnostics)) {
      gkyl_vlasov_app_write_species_fpo(app, i, tm, frame);
    }
  }
  for (int i=0; i<app->num_fluid_species; ++i) {
    gkyl_vlasov_app_write_fluid_species(app, i, tm, frame);
  }

  app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
}

void
gkyl_vlasov_app_write_field(gkyl_vlasov_app* app, double tm, int frame)
{
  const char *fmt = "%s-field_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, frame);

  if (app->use_gpu) {
    // copy data from device to host before writing it out
    gkyl_array_copy(app->field->em_host, app->field->em);
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, app->field->em_host, fileNm);
  }
  else {
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, app->field->em, fileNm);
  }
}

void
gkyl_vlasov_app_write_species(gkyl_vlasov_app* app, int sidx, double tm, int frame)
{
  struct vm_species *vm_s = &app->species[sidx];

  const char *fmt = "%s-%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, vm_s->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, vm_s->info.name, frame);

  if (app->use_gpu) {
    // copy data from device to host before writing it out
    gkyl_array_copy(vm_s->f_host, vm_s->f);
  }
  gkyl_comm_array_write(vm_s->comm, &vm_s->grid, &vm_s->local, 
    vm_s->f_host, fileNm);  

  if (vm_s->source_id) {
    if (vm_s->src.write_source) {
      // Write out the source distribution function
      const char *fmt_source = "%s-%s_source_%d.gkyl";
      int sz_source = gkyl_calc_strlen(fmt_source, app->name, vm_s->info.name, frame);
      char fileNm_source[sz_source+1]; // ensures no buffer overflow
      snprintf(fileNm_source, sizeof fileNm_source, fmt_source, app->name, vm_s->info.name, frame);
      // copy data from device to host before writing it out
      if (app->use_gpu) {
        gkyl_array_copy(vm_s->src.source_host, vm_s->src.source);
      }

      gkyl_comm_array_write(vm_s->comm, &vm_s->grid, &vm_s->local, 
        vm_s->src.source_host, fileNm); 
    }
  }  
}

void
gkyl_vlasov_app_write_species_lte(gkyl_vlasov_app* app, int sidx, double tm, int frame)
{
  const char *fmt = "%s-%s_%d_lte.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, app->species[sidx].info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, app->species[sidx].info.name, frame);

  if (app->species[sidx].info.output_f_lte) {
    vm_species_lte(app, &app->species[sidx], &app->species[sidx].lte, app->species[sidx].f);
  }
  
  if (app->use_gpu) {
    // copy data from device to host before writing it out
    gkyl_array_copy(app->species[sidx].f_host, app->species[sidx].lte.f_lte);
    gkyl_comm_array_write(app->species[sidx].comm, &app->species[sidx].grid, &app->species[sidx].local,
      app->species[sidx].f_host, fileNm);
  }
  else {
    gkyl_comm_array_write(app->species[sidx].comm, &app->species[sidx].grid, &app->species[sidx].local,
      app->species[sidx].lte.f_lte, fileNm);
  }
}

void
gkyl_vlasov_app_write_species_fpo(gkyl_vlasov_app* app, int sidx, double tm, int frame)
{
  const char *fmt_h = "%s-%s_H_%d.gkyl";
  const char *fmt_g = "%s-%s_G_%d.gkyl";
  const char *fmt_drag = "%s-%s_drag_coeff_%d.gkyl";
  const char *fmt_diff = "%s-%s_diff_coeff_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt_drag, app->name, app->species[sidx].info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow

  if (app->species[sidx].collision_id == GKYL_FPO_COLLISIONS) {
    vm_species_fpo_drag_diff_coeffs(app, &app->species[sidx], 
      &app->species[sidx].fpo, app->species[sidx].f);
  }
  
  if (app->use_gpu) {
    // copy data from device to host before writing it out
    gkyl_array_copy(app->species[sidx].fpo.h_host, app->species[sidx].fpo.h);
    gkyl_array_copy(app->species[sidx].fpo.g_host, app->species[sidx].fpo.g);
    gkyl_array_copy(app->species[sidx].fpo.drag_coeff_host, app->species[sidx].fpo.drag_coeff);
    gkyl_array_copy(app->species[sidx].fpo.diff_coeff_host, app->species[sidx].fpo.diff_coeff);

    // Write H
    snprintf(fileNm, sizeof fileNm, fmt_h, app->name, app->species[sidx].info.name, frame);
    gkyl_comm_array_write(app->species[sidx].comm, &app->species[sidx].grid, 
      &app->species[sidx].local, app->species[sidx].fpo.h_host, fileNm);

    // Write G
    snprintf(fileNm, sizeof fileNm, fmt_g, app->name, app->species[sidx].info.name, frame);
    gkyl_comm_array_write(app->species[sidx].comm, &app->species[sidx].grid, 
      &app->species[sidx].local, app->species[sidx].fpo.g_host, fileNm);

    // Write drag coefficient
    snprintf(fileNm, sizeof fileNm, fmt_drag, app->name, app->species[sidx].info.name, frame);
    gkyl_comm_array_write(app->species[sidx].comm, &app->species[sidx].grid, 
      &app->species[sidx].local, app->species[sidx].fpo.drag_coeff_host, fileNm);

    // Write diffusion tensor
    snprintf(fileNm, sizeof fileNm, fmt_diff, app->name, app->species[sidx].info.name, frame);
    gkyl_comm_array_write(app->species[sidx].comm, &app->species[sidx].grid, 
      &app->species[sidx].local, app->species[sidx].fpo.diff_coeff_host, fileNm);
  }
  else {
    snprintf(fileNm, sizeof fileNm, fmt_h, app->name, app->species[sidx].info.name, frame);
    gkyl_comm_array_write(app->species[sidx].comm, &app->species[sidx].grid, 
      &app->species[sidx].local, app->species[sidx].fpo.h, fileNm);

    snprintf(fileNm, sizeof fileNm, fmt_g, app->name, app->species[sidx].info.name, frame);
    gkyl_comm_array_write(app->species[sidx].comm, &app->species[sidx].grid, 
      &app->species[sidx].local, app->species[sidx].fpo.g, fileNm);

    snprintf(fileNm, sizeof fileNm, fmt_drag, app->name, app->species[sidx].info.name, frame);
    gkyl_comm_array_write(app->species[sidx].comm, &app->species[sidx].grid, 
      &app->species[sidx].local, app->species[sidx].fpo.drag_coeff, fileNm);

    snprintf(fileNm, sizeof fileNm, fmt_diff, app->name, app->species[sidx].info.name, frame);
    gkyl_comm_array_write(app->species[sidx].comm, &app->species[sidx].grid, 
      &app->species[sidx].local, app->species[sidx].fpo.diff_coeff, fileNm);
  }
}

void
gkyl_vlasov_app_write_fluid_species(gkyl_vlasov_app* app, int sidx, double tm, int frame)
{
  const char *fmt = "%s-%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, app->fluid_species[sidx].info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, app->fluid_species[sidx].info.name, frame);
  // copy data from device to host before writing it out
  if (app->use_gpu) 
    gkyl_array_copy(app->fluid_species[sidx].fluid_host, app->fluid_species[sidx].fluid);

  gkyl_comm_array_write(app->comm, &app->grid, &app->local,
    app->fluid_species[sidx].fluid_host, fileNm);
}

void
gkyl_vlasov_app_write_mom(gkyl_vlasov_app* app, double tm, int frame)
{
  for (int i=0; i<app->num_species; ++i) {
    struct vm_species *vm_s = &app->species[i];

    for (int m=0; m<vm_s->info.num_diag_moments; ++m) {

      const char *fmt = "%s-%s_%s_%d.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, vm_s->info.name,
        vm_s->info.diag_moments[m], frame);
      char fileNm[sz+1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, vm_s->info.name,
        vm_s->info.diag_moments[m], frame);

      if (app->use_gpu) {
        gkyl_array_copy(vm_s->moms[m].marr_host, vm_s->moms[m].marr);
      }
      gkyl_comm_array_write(app->comm, &app->grid, &app->local, vm_s->moms[m].marr_host, fileNm);

      if (vm_s->source_id) {
        if (vm_s->src.write_source) {
          const char *fmt_source = "%s-%s_source_%s_%d.gkyl";
          int sz_source = gkyl_calc_strlen(fmt, app->name, vm_s->info.name,
            vm_s->info.diag_moments[m], frame);
          char fileNm_source[sz_source+1]; // ensures no buffer overflow
          snprintf(fileNm_source, sizeof fileNm_source, fmt_source, app->name, vm_s->info.name,
            vm_s->info.diag_moments[m], frame);

          if (app->use_gpu) {
            gkyl_array_copy(vm_s->src.moms[m].marr_host, vm_s->src.moms[m].marr);
          }
          gkyl_comm_array_write(app->comm, &app->grid, &app->local, vm_s->src.moms[m].marr_host, fileNm_source); 
        }
      }      
    }
  }
}

void
gkyl_vlasov_app_write_integrated_mom(gkyl_vlasov_app *app)
{
  for (int i=0; i<app->num_species; ++i) {
    struct vm_species *vm_s = &app->species[i];

    int rank;
    gkyl_comm_get_rank(app->comm, &rank);
    if (rank == 0) {
      // write out integrated diagnostic moments
      const char *fmt = "%s-%s-%s.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, vm_s->info.name,
        "imom");
      char fileNm[sz+1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, vm_s->info.name,
        "imom");

      if (vm_s->is_first_integ_write_call) {
        gkyl_dynvec_write(vm_s->integ_diag, fileNm);
        vm_s->is_first_integ_write_call = false;
      }
      else {
        gkyl_dynvec_awrite(vm_s->integ_diag, fileNm);
      }

      if (vm_s->source_id) {
        if (vm_s->src.write_source) { 
          // write out integrated diagnostic moments from sources
          const char *fmt_source = "%s-%s-source-%s.gkyl";
          int sz_source = gkyl_calc_strlen(fmt, app->name, vm_s->info.name,
            "imom");
          char fileNm_source[sz_source+1]; // ensures no buffer overflow
          snprintf(fileNm_source, sizeof fileNm_source, fmt_source, app->name, vm_s->info.name,
            "imom");

          if (vm_s->src.is_first_integ_write_call) {
            gkyl_dynvec_write(vm_s->src.integ_diag, fileNm_source);
            vm_s->src.is_first_integ_write_call = false;
          }
          else {
            gkyl_dynvec_awrite(vm_s->src.integ_diag, fileNm_source);
          }          
        }
      }      
    }
    gkyl_dynvec_clear(vm_s->integ_diag);
    if (vm_s->source_id) {
      if (vm_s->src.write_source) { 
        gkyl_dynvec_clear(vm_s->src.integ_diag);
      }
    }
  }
}

void
gkyl_vlasov_app_write_integrated_L2_f(gkyl_vlasov_app* app)
{
  for (int i=0; i<app->num_species; ++i) {
    struct vm_species *vm_s = &app->species[i];

    int rank;
    gkyl_comm_get_rank(app->comm, &rank);
    if (rank == 0) {
      // write out integrated L^2
      const char *fmt = "%s-%s-%s.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, vm_s->info.name,
        "L2");
      char fileNm[sz+1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, vm_s->info.name,
        "L2");

      if (vm_s->is_first_integ_L2_write_call) {
        // write to a new file (this ensure previous output is removed)
        gkyl_dynvec_write(vm_s->integ_L2_f, fileNm);
        vm_s->is_first_integ_L2_write_call = false;
      }
      else {
        // append to existing file
        gkyl_dynvec_awrite(vm_s->integ_L2_f, fileNm);
      }
    }
    gkyl_dynvec_clear(vm_s->integ_L2_f);
  }
}

void
gkyl_vlasov_app_write_field_energy(gkyl_vlasov_app* app)
{
  // write out diagnostic moments
  const char *fmt = "%s-field-energy.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name);

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);

  if (rank == 0) {
    if (app->field->is_first_energy_write_call) {
      // write to a new file (this ensure previous output is removed)
      gkyl_dynvec_write(app->field->integ_energy, fileNm);
      app->field->is_first_energy_write_call = false;
    }
    else {
      // append to existing file
      gkyl_dynvec_awrite(app->field->integ_energy, fileNm);
    }
  }
  gkyl_dynvec_clear(app->field->integ_energy);
}

void
gkyl_vlasov_app_write_lte_corr_status(gkyl_vlasov_app* app)
{
  for (int i=0; i<app->num_species; ++i) {
    struct vm_species *vm_s = &app->species[i];

    if (vm_s->collision_id == GKYL_BGK_COLLISIONS) {

      int rank;
      gkyl_comm_get_rank(app->comm, &rank);
      if (rank == 0) {
         // write out correction statistics
        const char *fmt = "%s-%s-%s.gkyl";
        int sz = gkyl_calc_strlen(fmt, app->name, vm_s->info.name,
          "corr-lte-stat");
        char fileNm[sz+1]; // ensures no buffer overflow
        snprintf(fileNm, sizeof fileNm, fmt, app->name, vm_s->info.name,
          "corr-lte-stat"); 
          
        if (vm_s->bgk.lte.is_first_corr_status_write_call) {
          // write to a new file (this ensure previous output is removed)
          gkyl_dynvec_write(vm_s->bgk.lte.corr_stat, fileNm);
          vm_s->bgk.lte.is_first_corr_status_write_call = false;
        }
        else {
          // append to existing file
          gkyl_dynvec_awrite(vm_s->bgk.lte.corr_stat, fileNm);
        }
      }
      gkyl_dynvec_clear(vm_s->bgk.lte.corr_stat);
    }
  } 
}

struct gkyl_update_status
gkyl_vlasov_update(gkyl_vlasov_app* app, double dt)
{
  app->stat.nup += 1;

  struct timespec wst = gkyl_wall_clock();
  struct gkyl_update_status status = app->update_func(app, dt);
  app->tcurr += status.dt_actual;

  app->stat.total_tm += gkyl_time_diff_now_sec(wst);

  // Check for any CUDA errors during time step
  if (app->use_gpu)
    checkCuda(cudaGetLastError());
  return status;
}

struct gkyl_vlasov_stat
gkyl_vlasov_app_stat(gkyl_vlasov_app* app)
{
  vm_species_coll_tm(app);
  vm_species_bgk_niter(app);
  vm_species_tm(app);
  vm_species_rad_tm(app);
  return app->stat;
}

void
gkyl_vlasov_app_species_ktm_rhs(gkyl_vlasov_app* app, int update_vol_term)
{
  for (int i=0; i<app->num_species; ++i) {

    struct vm_species *species = &app->species[i];

    const struct gkyl_array *fin = species->f;
    struct gkyl_array *rhs = species->f1;

    gkyl_array_clear(rhs, 0.0);
    gkyl_dg_updater_vlasov_advance(species->slvr, &species->local, 
      fin, species->cflrate, rhs); 
  }
}

static void
range_stat_write(gkyl_vlasov_app* app, const char *nm, const struct gkyl_range *r, FILE *fp)
{
  gkyl_vlasov_app_cout(app, fp, " %s_cells : [ ", nm);
  for (int i=0; i<r->ndim; ++i)
    gkyl_vlasov_app_cout(app, fp, " %d, ", gkyl_range_shape(r, i));
  gkyl_vlasov_app_cout(app, fp, " ],\n");
}

// ensure stats across processors are made consistent
static void
comm_reduce_app_stat(const gkyl_vlasov_app* app,
  const struct gkyl_vlasov_stat *local, struct gkyl_vlasov_stat *global)
{
  int comm_sz;
  gkyl_comm_get_size(app->comm, &comm_sz);
  if (comm_sz == 1) {
    memcpy(global, local, sizeof(struct gkyl_vlasov_stat));
    return;
  }

  global->use_gpu = local->use_gpu;

  enum { NUP, NFEULER, NSTAGE_2_FAIL, NSTAGE_3_FAIL, L_END };
  int64_t l_red[] = {
    [NUP] = local->nup,
    [NFEULER] = local->nfeuler,
    [NSTAGE_2_FAIL] = local->nstage_2_fail,
    [NSTAGE_3_FAIL] = local->nstage_3_fail
  };

  int64_t l_red_global[L_END];
  gkyl_comm_all_reduce(app->comm, GKYL_INT_64, GKYL_MAX, L_END, l_red, l_red_global);

  global->nup = l_red_global[NUP];
  global->nfeuler = l_red_global[NFEULER];
  global->nstage_2_fail = l_red_global[NSTAGE_2_FAIL];
  global->nstage_3_fail = l_red_global[NSTAGE_3_FAIL];  

  int64_t l_red_bgk_corr[app->num_species];
  for (int s=0; s<app->num_species; ++s) {
    l_red_bgk_corr[s] = local->niter_self_bgk_corr[s];
  }

  int64_t l_red_global_bgk_corr[app->num_species];
  gkyl_comm_all_reduce(app->comm, GKYL_INT_64, GKYL_MAX, app->num_species, 
    l_red_bgk_corr, l_red_global_bgk_corr);

  for (int s=0; s<app->num_species; ++s) {
    global->niter_self_bgk_corr[s] = l_red_bgk_corr[s];
  }

  enum {
    TOTAL_TM, RK3_TM, FL_EM_TM, 
    INIT_SPECIES_TM, INIT_FLUID_SPECIES_TM, INIT_FIELD_TM, 
    SPECIES_RHS_TM, FLUID_SPECIES_RHS_TM, FLUID_SPECIES_VARS_TM, 
    SPECIES_COLL_MOM_TM, SPECIES_COL_TM, SPECIES_RAD_TM, SPECIES_LTE_TM, 
    FIELD_RHS_TM, CURRENT_TM,
    SPECIES_OMEGA_CFL_TM, FIELD_OMEGA_CFL_TM, MOM_TM, DIAG_TM, IO_TM,
    SPECIES_BC_TM, FLUID_SPECIES_BC_TM, FIELD_BC_TM,
    D_END
  };

  double d_red[D_END] = {
    [TOTAL_TM] = local->total_tm,
    [RK3_TM] = local->rk3_tm,
    [FL_EM_TM] = local->fl_em_tm,
    [INIT_SPECIES_TM] = local->init_species_tm,
    [INIT_FLUID_SPECIES_TM] = local->init_fluid_species_tm,
    [INIT_FIELD_TM] = local->field_rhs_tm,
    [SPECIES_RHS_TM] = local->species_rhs_tm,
    [FLUID_SPECIES_RHS_TM] = local->fluid_species_rhs_tm,
    [FLUID_SPECIES_VARS_TM] = local->fluid_species_vars_tm,
    [SPECIES_COLL_MOM_TM] = local->species_coll_mom_tm,
    [SPECIES_COL_TM] = local->species_coll_tm,
    [SPECIES_RAD_TM] = local->species_rad_tm,
    [SPECIES_LTE_TM] = local->species_lte_tm,
    [FIELD_RHS_TM] = local->field_rhs_tm,
    [CURRENT_TM] = local->current_tm,
    [SPECIES_OMEGA_CFL_TM] = local->species_omega_cfl_tm,
    [FIELD_OMEGA_CFL_TM] = local->field_omega_cfl_tm,
    [MOM_TM] = local->mom_tm,
    [DIAG_TM] = local->diag_tm,
    [IO_TM] = local->io_tm,
    [SPECIES_BC_TM] = local->species_bc_tm,
    [FLUID_SPECIES_BC_TM] = local->fluid_species_bc_tm,
    [FIELD_BC_TM] = local->field_bc_tm
  };

  double d_red_global[D_END];
  gkyl_comm_all_reduce(app->comm, GKYL_DOUBLE, GKYL_MAX, D_END, d_red, d_red_global);
  
  global->total_tm = d_red_global[TOTAL_TM];
  global->rk3_tm = d_red_global[RK3_TM];
  global->fl_em_tm = d_red_global[FL_EM_TM];
  global->init_species_tm = d_red_global[INIT_SPECIES_TM];
  global->init_fluid_species_tm = d_red_global[INIT_FLUID_SPECIES_TM];
  global->field_rhs_tm = d_red_global[INIT_FIELD_TM];
  global->species_rhs_tm = d_red_global[SPECIES_RHS_TM];
  global->fluid_species_rhs_tm = d_red_global[FLUID_SPECIES_RHS_TM];
  global->fluid_species_vars_tm = d_red_global[FLUID_SPECIES_VARS_TM];
  global->species_coll_mom_tm = d_red_global[SPECIES_COLL_MOM_TM];
  global->species_coll_tm = d_red_global[SPECIES_COL_TM];
  global->species_rad_tm = d_red_global[SPECIES_RAD_TM];
  global->species_lte_tm = d_red_global[SPECIES_LTE_TM];
  global->field_rhs_tm = d_red_global[FIELD_RHS_TM];
  global->current_tm = d_red_global[CURRENT_TM];
  global->species_omega_cfl_tm = d_red_global[SPECIES_OMEGA_CFL_TM];
  global->field_omega_cfl_tm = d_red_global[FIELD_OMEGA_CFL_TM];
  global->mom_tm = d_red_global[MOM_TM];
  global->diag_tm = d_red_global[DIAG_TM];
  global->io_tm = d_red_global[IO_TM];
  global->species_bc_tm = d_red_global[SPECIES_BC_TM];
  global->fluid_species_bc_tm = d_red_global[FLUID_SPECIES_BC_TM];
  global->field_bc_tm = d_red_global[FIELD_BC_TM];

  // misc data needing reduction

  gkyl_comm_all_reduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 2, local->stage_2_dt_diff,
    global->stage_2_dt_diff);
  gkyl_comm_all_reduce(app->comm, GKYL_DOUBLE, GKYL_MAX, 2, local->stage_3_dt_diff,
    global->stage_3_dt_diff);

  gkyl_comm_all_reduce(app->comm, GKYL_DOUBLE, GKYL_MAX, GKYL_MAX_SPECIES, local->species_lbo_coll_drag_tm,
    global->species_lbo_coll_drag_tm);
  gkyl_comm_all_reduce(app->comm, GKYL_DOUBLE, GKYL_MAX, GKYL_MAX_SPECIES, local->species_lbo_coll_diff_tm,
    global->species_lbo_coll_diff_tm);
}

void
gkyl_vlasov_app_stat_write(gkyl_vlasov_app* app)
{
  const char *fmt = "%s-%s";
  int sz = gkyl_calc_strlen(fmt, app->name, "stat.json");
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, "stat.json");

  int num_ranks;
  gkyl_comm_get_size(app->comm, &num_ranks);

  char buff[70];
  time_t t = time(NULL);
  struct tm curr_tm = *localtime(&t);

  vm_species_coll_tm(app);
  vm_species_bgk_niter(app);
  vm_species_tm(app);
  vm_species_rad_tm(app);

  struct gkyl_vlasov_stat stat = { };
  comm_reduce_app_stat(app, &app->stat, &stat);
  
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  // append to existing file so we have a history of different runs
  FILE *fp = 0;
  if (rank == 0) fp = fopen(fileNm, "a");

  gkyl_vlasov_app_cout(app, fp, "{\n");

  if (strftime(buff, sizeof buff, "%c", &curr_tm))
    gkyl_vlasov_app_cout(app, fp, " date : %s,\n", buff);

  gkyl_vlasov_app_cout(app, fp, " use_gpu : %d,\n", stat.use_gpu);
  gkyl_vlasov_app_cout(app, fp, " num_ranks : %d,\n", num_ranks); 
  
  for (int s=0; s<app->num_species; ++s)
    range_stat_write(app, app->species[s].info.name, &app->species[s].global, fp);
  
  gkyl_vlasov_app_cout(app, fp, " nup : %ld,\n", stat.nup);
  gkyl_vlasov_app_cout(app, fp, " nfeuler : %ld,\n", stat.nfeuler);
  gkyl_vlasov_app_cout(app, fp, " nstage_2_fail : %ld,\n", stat.nstage_2_fail);
  gkyl_vlasov_app_cout(app, fp, " nstage_3_fail : %ld,\n", stat.nstage_3_fail);

  gkyl_vlasov_app_cout(app, fp, " stage_2_dt_diff : [ %lg, %lg ],\n",
    stat.stage_2_dt_diff[0], stat.stage_2_dt_diff[1]);
  gkyl_vlasov_app_cout(app, fp, " stage_3_dt_diff : [ %lg, %lg ],\n",
    stat.stage_3_dt_diff[0], stat.stage_3_dt_diff[1]);

  gkyl_vlasov_app_cout(app, fp, " total_tm : %lg,\n", stat.total_tm);
  gkyl_vlasov_app_cout(app, fp, " rk3_tm : %lg,\n", stat.rk3_tm);
  gkyl_vlasov_app_cout(app, fp, " fluid_em_coupling_tm : %lg,\n", stat.fl_em_tm);
  gkyl_vlasov_app_cout(app, fp, " init_species_tm : %lg,\n", stat.init_species_tm);
  if (app->has_field)
    gkyl_vlasov_app_cout(app, fp, " init_field_tm : %lg,\n", stat.init_field_tm);
  
  gkyl_vlasov_app_cout(app, fp, " species_rhs_tm : %lg,\n", stat.species_rhs_tm);

  for (int s=0; s<app->num_species; ++s) {
    gkyl_vlasov_app_cout(app, fp, " species_coll_drag_tm[%d] : %lg,\n", s,
      stat.species_lbo_coll_drag_tm[s]);
    gkyl_vlasov_app_cout(app, fp, " species_coll_diff_tm[%d] : %lg,\n", s,
      stat.species_lbo_coll_diff_tm[s]);
    gkyl_vlasov_app_cout(app, fp, " niter_self_bgk_corr[%d] : %ld,\n", s, 
      stat.niter_self_bgk_corr[s]);
  }

  gkyl_vlasov_app_cout(app, fp, " species_coll_mom_tm : %lg,\n", stat.species_coll_mom_tm);
  gkyl_vlasov_app_cout(app, fp, " species_coll_tm : %lg,\n", stat.species_coll_tm);
  gkyl_vlasov_app_cout(app, fp, " species_rad_tm : %lg,\n", stat.species_rad_tm);

  gkyl_vlasov_app_cout(app, fp, " species_lte_tm : %lg,\n", stat.species_lte_tm);

  gkyl_vlasov_app_cout(app, fp, " species_bc_tm : %lg,\n", stat.species_bc_tm);
  
  gkyl_vlasov_app_cout(app, fp, " fluid_species_rhs_tm : %lg,\n", stat.fluid_species_rhs_tm);

  gkyl_vlasov_app_cout(app, fp, " fluid_species_bc_tm : %lg,\n", stat.fluid_species_bc_tm);

  if (app->has_field) {
    gkyl_vlasov_app_cout(app, fp, " field_rhs_tm : %lg,\n", stat.field_rhs_tm);
    gkyl_vlasov_app_cout(app, fp, " field_bc_tm : %lg,\n", stat.field_bc_tm);
    
    gkyl_vlasov_app_cout(app, fp, " current_tm : %lg,\n", stat.current_tm);
  }

  gkyl_vlasov_app_cout(app, fp, " nmom : %ld,\n", stat.nmom);
  gkyl_vlasov_app_cout(app, fp, " mom_tm : %lg\n", stat.mom_tm);

  gkyl_vlasov_app_cout(app, fp, " ndiag : %ld,\n", stat.ndiag);
  gkyl_vlasov_app_cout(app, fp, " diag_tm : %lg\n", stat.diag_tm);
  
  gkyl_vlasov_app_cout(app, fp, " nspecies_omega_cfl : %ld,\n", stat.nspecies_omega_cfl);
  gkyl_vlasov_app_cout(app, fp, " species_omega_cfl_tm : %lg\n", stat.species_omega_cfl_tm);

  gkyl_vlasov_app_cout(app, fp, " nfield_omega_cfl : %ld,\n", stat.nfield_omega_cfl);
  gkyl_vlasov_app_cout(app, fp, " field_omega_cfl_tm : %lg\n", stat.field_omega_cfl_tm);

  gkyl_vlasov_app_cout(app, fp, " nio : %ld,\n", stat.nio);
  gkyl_vlasov_app_cout(app, fp, " io_tm : %lg\n", stat.io_tm);
  
  gkyl_vlasov_app_cout(app, fp, "}\n");

  if (rank == 0)
    fclose(fp);  

}

// private function to handle variable argument list for printing
static void
v_vlasov_app_cout(const gkyl_vlasov_app* app, FILE *fp, const char *fmt, va_list argp)
{
  int rank, r = 0;
  gkyl_comm_get_rank(app->comm, &rank);
  if ((rank == 0) && fp)
    vfprintf(fp, fmt, argp);
}

void
gkyl_vlasov_app_cout(const gkyl_vlasov_app* app, FILE *fp, const char *fmt, ...)
{
  va_list argp;
  va_start(argp, fmt);
  v_vlasov_app_cout(app, fp, fmt, argp);
  va_end(argp);
}

void
gkyl_vlasov_app_release(gkyl_vlasov_app* app)
{
  for (int i=0; i<app->num_species; ++i)
    vm_species_release(app, &app->species[i]);
  for (int i=0; i<app->num_fluid_species; ++i)
    vm_fluid_species_release(app, &app->fluid_species[i]);
  if (app->num_species > 0)
    gkyl_free(app->species);
  if (app->num_fluid_species > 0)
    gkyl_free(app->fluid_species);
  if (app->has_field)
    vm_field_release(app, app->field);
  if (app->has_fluid_em_coupling)
    vm_fluid_em_coupling_release(app, app->fl_em);

  gkyl_comm_release(app->comm);

  gkyl_wave_geom_release(app->geom);

  if (app->use_gpu) {
    gkyl_cu_free(app->basis_on_dev.basis);
    gkyl_cu_free(app->basis_on_dev.confBasis);
  }

  gkyl_free(app);
}
