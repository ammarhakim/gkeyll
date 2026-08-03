#include <gkyl_gk_block_geom.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_elem_type_priv.h>
#include <gkyl_efit.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_gyrokinetic_multib_priv.h>
#include <gkyl_multib_conn.h>

#include <mpack.h>
#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

// Compute total number of ranges specified by cuts.
static inline int
calc_cuts(int ndim, const int *cuts)
{
  int tc = 1;
  for (int d=0; d<ndim; ++d) tc *= cuts[d];
  return tc;
}

// Simple linear search to check if val occurs in lst.
static bool
has_int(int n, int val, const int *lst)
{
  for (int i=0; i<n; ++i)
    if (val == lst[i])
      return true;
  return false;
}

// Compute total and maximum number of cuts.
static void
calc_tot_and_max_cuts(const struct gkyl_gk_block_geom *gk_block_geom, int tot_max[2])
{
  int ndim = gkyl_gk_block_geom_ndim(gk_block_geom);
  int num_blocks = gkyl_gk_block_geom_num_blocks(gk_block_geom);

  int max_cuts = 0, tot_cuts = 0;
  for (int i=0; i<num_blocks; ++i) {
    const struct gkyl_gk_block_geom_info *bgi = gkyl_gk_block_geom_get_block(gk_block_geom, i);
    int ncuts = calc_cuts(ndim, bgi->cuts);
    max_cuts = ncuts > max_cuts ? ncuts : max_cuts;
    tot_cuts += ncuts;
  }
  tot_max[0] = tot_cuts;
  tot_max[1] = max_cuts;
}

// Construct single-block App geometry for given block ID.
static struct gkyl_gyrokinetic_app *
singleb_app_new_geom(const struct gkyl_gyrokinetic_multib *mbinp, int bid,
  const struct gkyl_gyrokinetic_multib_app *mbapp)
{
  const char *tok_geo_trace = getenv("GKYL_TOK_GEO_TRACE");
  bool trace_tok_geo = tok_geo_trace && tok_geo_trace[0] != '\0' && tok_geo_trace[0] != '0';
  if (trace_tok_geo) {
    char block_id[32];
    snprintf(block_id, sizeof block_id, "%d", bid);
    setenv("GKYL_TOK_GEO_TRACE_BLOCK", block_id, 1);
    fprintf(stderr, "GKYL_TOK_GEO_TRACE block_start bid=%d\n", bid);
    fflush(stderr);
  }

  // For kinetic simulations, block dimension defined configuration-space dimensionality.
  int cdim = gkyl_gk_block_geom_ndim(mbapp->gk_block_geom);
  int num_blocks = gkyl_gk_block_geom_num_blocks(mbapp->gk_block_geom);

  const struct gkyl_gk_block_geom_info *bgi =
    gkyl_gk_block_geom_get_block(mbapp->gk_block_geom, bid);

  // Construct top-level single-block input on the heap; this struct is large
  // enough to overflow the default stack in geometry-only multiblock tests.
  struct gkyl_gk *app_inp = gkyl_calloc(1, sizeof(*app_inp));

  strcpy(app_inp->name, mbinp->name);
  if (num_blocks > 1) {
    cstr app_name = cstr_from_fmt("%s_b%d", mbinp->name, bid);
    strcpy(app_inp->name, app_name.str);
    cstr_drop(&app_name);
  }

  // Set the configuration-space extents, cells, and geometry.
  app_inp->cdim = cdim;
  for (int i=0; i<cdim; ++i) {
    app_inp->lower[i] = bgi->lower[i];
    app_inp->upper[i] = bgi->upper[i];
    app_inp->cells[i] = bgi->cells[i];
  }

  // Set z dir grid extents based on tokamak global normalization
  if (bgi->geometry.geometry_id == GKYL_GEOMETRY_TOKAMAK || bgi->geometry.geometry_id == GKYL_GEOMETRY_FROMFILE) {
    gkyl_gk_geometry_tok_set_grid_extents(bgi->geometry.efit_info, bgi->geometry.tok_grid_info, &app_inp->lower[cdim-1], &app_inp->upper[cdim-1]);
    gkyl_gk_block_geom_reset_block_extents(mbapp->gk_block_geom, bid, app_inp->lower, app_inp->upper);
  }

  app_inp->geometry = bgi->geometry;
  // This constructor only builds geometry; species and neutral inputs are not
  // populated in the single-block geometry app.
  app_inp->num_species = 0;
  app_inp->num_neut_species = 0;

  app_inp->poly_order = mbinp->poly_order;
  app_inp->basis_type = mbinp->basis_type;
  app_inp->cfl_frac = mbinp->cfl_frac;
  app_inp->cfl_frac_omegaH = mbinp->cfl_frac_omegaH;

  struct gkyl_comm *comm = mbapp->block_comms[bid];

  struct gkyl_app_parallelism_inp parallel_inp = {};
  parallel_inp.use_gpu = mbinp->use_gpu;
  for (int d=0; d<cdim; ++d) parallel_inp.cuts[d] = bgi->cuts[d];
  parallel_inp.comm = comm;
  // Copy parallelism input into app input.
  memcpy(&app_inp->parallelism, &parallel_inp, sizeof(struct gkyl_app_parallelism_inp));

  app_inp->num_periodic_dir = mbinp->num_periodic_dir;
  for(int i = 0; i < mbinp->cdim; i++)
    app_inp->periodic_dirs[i] = mbinp->periodic_dirs[i];

  app_inp->metadata.num_attributes = mbapp->io_meta_basic_len;
  app_inp->metadata.attributes = mbapp->io_meta_basic;

  struct gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new_geom(app_inp);
  gkyl_free(app_inp);
  if (trace_tok_geo) {
    fprintf(stderr, "GKYL_TOK_GEO_TRACE block_done bid=%d\n", bid);
    fflush(stderr);
  }
  return app;
}

// Construct single-block App solver for given block ID.
static void
singleb_app_new_solver(const struct gkyl_gyrokinetic_multib *mbinp, int bid,
  const struct gkyl_gyrokinetic_multib_app *mbapp, struct gkyl_gyrokinetic_app *app)
{
  // For kinetic simulations, block dimension defined configuration-space dimensionality.
  int cdim = gkyl_gk_block_geom_ndim(mbapp->gk_block_geom);
  int num_blocks = gkyl_gk_block_geom_num_blocks(mbapp->gk_block_geom);

  const struct gkyl_gk_block_geom_info *bgi =
    gkyl_gk_block_geom_get_block(mbapp->gk_block_geom, bid);

  // Construct top-level single-block input on the heap; this struct is large
  // enough to overflow the default stack in multiblock solver construction.
  struct gkyl_gk *app_inp = gkyl_calloc(1, sizeof(*app_inp));

  // Set the configuration-space extents, cells.
  app_inp->cdim = cdim;
  for (int i=0; i<cdim; ++i) {
    app_inp->lower[i] = bgi->lower[i];
    app_inp->upper[i] = bgi->upper[i];
    app_inp->cells[i] = bgi->cells[i];
  }

  int num_species = app_inp->num_species = mbinp->num_species;
  int num_neut_species = app_inp->num_neut_species = mbinp->num_neut_species; 

  app_inp->poly_order = mbinp->poly_order;
  app_inp->basis_type = mbinp->basis_type;
  app_inp->cfl_frac = mbinp->cfl_frac; 
  app_inp->cfl_frac_omegaH = mbinp->cfl_frac_omegaH; 
  app_inp->eirene = mbinp->eirene;

  for (int i=0; i<num_species; ++i) {
    const struct gkyl_gyrokinetic_multib_species *sp = &mbinp->species[i];
    
    struct gkyl_gyrokinetic_species species_inp = { };
    strcpy(species_inp.name, sp->name);

    species_inp.charge = sp->charge;
    species_inp.mass = sp->mass;

    // Velocity-space information
    int vdim = species_inp.vdim = sp->vdim;
    for (int v=0; v<vdim; ++v) {
      species_inp.lower[v] = sp->lower[v];
      species_inp.upper[v] = sp->upper[v];
      species_inp.cells[v] = sp->cells[v];      
    }
    species_inp.mapc2p = sp->mapc2p;

    // Species physics modules.
    species_inp.collisionless = sp->collisionless;
    species_inp.collisions = sp->collisions;
    species_inp.anomalous_diffusion = sp->anomalous_diffusion;
    species_inp.radiation = sp->radiation;
    species_inp.react = sp->react;
    species_inp.react_neut = sp->react_neut; 

    // Species diagnostics
    species_inp.num_diag_moments = sp->num_diag_moments;
    for (int n=0; n<species_inp.num_diag_moments; ++n) {
      species_inp.diag_moments[n] = sp->diag_moments[n];
    }
    species_inp.num_integrated_diag_moments = sp->num_integrated_diag_moments;
    for (int n=0; n<species_inp.num_integrated_diag_moments; ++n) {
      species_inp.integrated_diag_moments[n] = sp->integrated_diag_moments[n];
    }
    species_inp.time_rate_diagnostics = sp->time_rate_diagnostics;
    species_inp.boundary_flux_diagnostics = sp->boundary_flux_diagnostics;

    // Choose proper block-specific species input.
    const struct gkyl_gyrokinetic_multib_species_pb *sp_pb = &sp->blocks[0];
    if (!sp->duplicate_across_blocks) {
      for (int i=0; i<num_blocks; ++i) {
        if (bid == sp->blocks[i].block_id) {
          sp_pb = &sp->blocks[i];
          break;
        }
      }
    }
    species_inp.projection = sp_pb->projection;
    species_inp.source = sp_pb->source;
    species_inp.polarization_density = sp_pb->polarization_density;  

    // By default, skip BCs altogether.
    for (int i=0; i<2*GKYL_MAX_CDIM; i++) {
      species_inp.bcs[i].type = GKYL_BC_GK_SKIP;
      species_inp.bcs[i].type = GKYL_BC_GK_SKIP;
    }

    int pardir = cdim-1;
    int num_below = gkyl_multib_conn_get_num_connected(mbapp->block_topo,
      bid, pardir, 0, GKYL_CONN_BELOW);
    int num_above = gkyl_multib_conn_get_num_connected(mbapp->block_topo,
      bid, pardir, 0, GKYL_CONN_ABOVE);

    int bc_count_sp[num_blocks];
    for (int i=0; i<num_blocks; i++)
      bc_count_sp[i] = 0;

    if (cdim ==3) {
      if (num_below > 0) {
        species_inp.bcs[bc_count_sp[bid]].dir = pardir;
        species_inp.bcs[bc_count_sp[bid]].edge = GKYL_LOWER_EDGE;
        species_inp.bcs[bc_count_sp[bid]].type = GKYL_BC_GK_SPECIES_TWISTSHIFT;
        bc_count_sp[bid] += 1;
      }

      if (num_above > 0) {
        species_inp.bcs[bc_count_sp[bid]].dir = pardir;
        species_inp.bcs[bc_count_sp[bid]].edge = GKYL_UPPER_EDGE;
        species_inp.bcs[bc_count_sp[bid]].type = GKYL_BC_GK_SPECIES_TWISTSHIFT;
        bc_count_sp[bid] += 1;
      }
    }

    // Set species physical BCs.
    for (int i=0; i<sp->num_physical_bcs; ++i) {
      if (bid == sp->bcs[i].bidx) {
        species_inp.bcs[bc_count_sp[bid]].dir = sp->bcs[i].dir;
        species_inp.bcs[bc_count_sp[bid]].edge = sp->bcs[i].edge;
        species_inp.bcs[bc_count_sp[bid]].type = sp->bcs[i].type;
        species_inp.bcs[bc_count_sp[bid]].aux_profile = sp->bcs[i].aux_profile;
        species_inp.bcs[bc_count_sp[bid]].aux_ctx = sp->bcs[i].aux_ctx;
        species_inp.bcs[bc_count_sp[bid]].projection = sp->bcs[i].projection;
        for (int k=0; k<3; ++k)
          species_inp.bcs[bc_count_sp[bid]].value[k] = sp->bcs[i].value[k];

        bc_count_sp[bid] += 1;
      }
    }

    // Copy species input into app input.
    memcpy(&app_inp->species[i], &species_inp, sizeof(struct gkyl_gyrokinetic_species));
  }

  for (int i=0; i<num_neut_species; ++i) {
    const struct gkyl_gyrokinetic_multib_neut_species *nsp = &mbinp->neut_species[i];
    
    struct gkyl_gyrokinetic_neut_species neut_species_inp = { };
    strcpy(neut_species_inp.name, nsp->name);

    neut_species_inp.mass = nsp->mass; 
    neut_species_inp.is_static = nsp->is_static; 

    // Velocity space information (neutrals are 3V)
    int vdim = neut_species_inp.vdim = nsp->vdim;
    for (int v=0; v<vdim; ++v) {
      neut_species_inp.lower[v] = nsp->lower[v];
      neut_species_inp.upper[v] = nsp->upper[v];
      neut_species_inp.cells[v] = nsp->cells[v];
    }
    neut_species_inp.mapc2p = nsp->mapc2p;

    // Neutral species physics modules
    neut_species_inp.react_neut = nsp->react_neut;

    // Neutral species diagnostics
    neut_species_inp.num_diag_moments = nsp->num_diag_moments;
    for (int n=0; n<neut_species_inp.num_diag_moments; ++n) {
      neut_species_inp.diag_moments[n] = nsp->diag_moments[n];
    }

    // choose proper block-specific species input
    const struct gkyl_gyrokinetic_multib_neut_species_pb *nsp_pb = &nsp->blocks[0];
    if (!nsp->duplicate_across_blocks) {
      for (int i=0; i<num_blocks; ++i) {
        if (bid == nsp->blocks[i].block_id) {
          nsp_pb = &nsp->blocks[i];
          break;
        }
      }
    }
    neut_species_inp.projection = nsp_pb->projection;
    neut_species_inp.source = nsp_pb->source;

    // By default, skip BCs altogether.
    for (int i=0; i<2*GKYL_MAX_CDIM; i++) {
      neut_species_inp.bcs[i].type = GKYL_BC_GK_SKIP;
      neut_species_inp.bcs[i].type = GKYL_BC_GK_SKIP;
    }

    // Set species physical BCs.
    int bc_count_nsp[num_blocks];
    for (int i=0; i<num_blocks; i++)
      bc_count_nsp[i] = 0;

    for (int i=0; i<nsp->num_physical_bcs; ++i) {
      if (bid == nsp->bcs[i].bidx) {
        neut_species_inp.bcs[bc_count_nsp[bid]].dir = nsp->bcs[i].dir;
        neut_species_inp.bcs[bc_count_nsp[bid]].edge = nsp->bcs[i].edge;
        neut_species_inp.bcs[bc_count_nsp[bid]].type = nsp->bcs[i].type;
        neut_species_inp.bcs[bc_count_nsp[bid]].aux_profile = nsp->bcs[i].aux_profile;
        neut_species_inp.bcs[bc_count_nsp[bid]].aux_ctx = nsp->bcs[i].aux_ctx;
        neut_species_inp.bcs[bc_count_nsp[bid]].projection = nsp->bcs[i].projection;
        for (int k=0; k<3; ++k)
          neut_species_inp.bcs[bc_count_nsp[bid]].value[k] = nsp->bcs[i].value[k];

        bc_count_nsp[bid] += 1;
      }
    }

    // Copy neutral species input into app input.
    memcpy(&app_inp->neut_species[i], &neut_species_inp, sizeof(struct gkyl_gyrokinetic_neut_species));
  } 

  // Initialize the single-block field solver (only used for num_blocks=1).
  const struct gkyl_gyrokinetic_multib_field *fld = &mbinp->field;
  struct gkyl_gyrokinetic_field field_inp = { };
  field_inp.gkfield_id = fld->gkfield_id;
  field_inp.kperpSq = fld->kperpSq; 
  field_inp.time_rate_diagnostics = fld->time_rate_diagnostics; 
  
  // Adiabatic electron inputs.
  field_inp.electron_mass = fld->electron_mass;
  field_inp.electron_charge = fld->electron_charge;
  field_inp.electron_density = fld->electron_density; 
  field_inp.electron_temp = fld->electron_temp; 
  
  // BCs.
  for (int d=0; d<cdim; d++) {
    for (int i=0; i<2; i++) {
      // Set it to Dirichlet first, reset below. This avoids problems in
      // creating the single block field solve, which may not be used.
      field_inp.poisson_bcs[2*d+i].dir = d;
      field_inp.poisson_bcs[2*d+i].edge = i==0? GKYL_LOWER_EDGE : GKYL_UPPER_EDGE;
      field_inp.poisson_bcs[2*d+i].type = GKYL_BC_GK_FIELD_DIRICHLET;
      for (int k=0; k<3; ++k)
        field_inp.poisson_bcs[i].value[k] = -1.0e3;
    }
  }
  
  for (int i=0; i<fld->num_physical_bcs; i++) { 
    if (bid == fld->bcs[i].bidx) {
      struct gkyl_gyrokinetic_bc *bc_curr = gk_fetch_bc_with_dir_edge(field_inp.poisson_bcs, 2*cdim, fld->bcs[i].dir, fld->bcs[i].edge);
      bc_curr->type = fld->bcs[i].type;
      bc_curr->aux_profile = fld->bcs[i].aux_profile;
      bc_curr->aux_ctx = fld->bcs[i].aux_ctx;
      for (int k=0; k<3; ++k)
        bc_curr->value[k] = fld->bcs[i].value[k];
    }
  }
  
  const struct gkyl_gyrokinetic_multib_field_pb *fld_pb = &fld->blocks[0];
  // Choose proper block-specific field input.
  if (!fld->duplicate_across_blocks) {
    for (int i=0; i<num_blocks; ++i) {
      if (bid == fld->blocks[i].block_id) {
        const struct gkyl_gyrokinetic_multib_field_pb *fld_pb = &fld->blocks[i];
        break;
      }
    }
  }
  
  if (!fld->duplicate_across_blocks) {
    for (int i=0; i<num_blocks; ++i) {
      if (bid == fld->blocks[i].block_id) {
        fld_pb = &fld->blocks[i];
        break;
      }
    }
  }
  
  field_inp.polarization_bmag = fld_pb->polarization_bmag ? fld_pb->polarization_bmag : mbapp->bmag_ref;
  field_inp.kperpSq = fld_pb->kperpSq;
  field_inp.time_rate_diagnostics = fld_pb->time_rate_diagnostics; 
  
  field_inp.phi_wall_lo_ctx = fld_pb->phi_wall_lo_ctx; 
  field_inp.phi_wall_lo = fld_pb->phi_wall_lo; 
  field_inp.phi_wall_lo_evolve = fld_pb->phi_wall_lo_evolve; 
  
  field_inp.phi_wall_up_ctx = fld_pb->phi_wall_up_ctx; 
  field_inp.phi_wall_up = fld_pb->phi_wall_up; 
  field_inp.phi_wall_up_evolve = fld_pb->phi_wall_up_evolve;   
  

  // Copy field input into app input.
  memcpy(&app_inp->field, &field_inp, sizeof(struct gkyl_gyrokinetic_field));  

  gkyl_gyrokinetic_app_new_solver(app_inp, app);
  gkyl_free(app_inp);
}

gkyl_gyrokinetic_multib_app*
gkyl_gyrokinetic_multib_app_new_geom(const struct gkyl_gyrokinetic_multib *mbinp)
{
  int my_rank, num_ranks;
  gkyl_comm_get_rank(mbinp->comm, &my_rank);
  gkyl_comm_get_size(mbinp->comm, &num_ranks);

  int tot_max[2];
  calc_tot_and_max_cuts(mbinp->gk_block_geom, tot_max);
  if ((num_ranks > tot_max[0]) || (num_ranks < tot_max[1])) {
    fprintf(stderr, "\nSpecified %d total cuts but provided %d processes, \
and the maximum number of cuts in a block is %d\n\n", tot_max[0], num_ranks, tot_max[1]);
    return 0;
  }

  struct gkyl_gyrokinetic_multib_app *mbapp = gkyl_malloc(sizeof(*mbapp));

  strcpy(mbapp->name, mbinp->name);
  mbapp->comm = gkyl_comm_acquire(mbinp->comm);  
  mbapp->use_gpu = mbinp->use_gpu;
  
  mbapp->gk_block_geom = gkyl_gk_block_geom_acquire(mbinp->gk_block_geom);
  mbapp->block_topo = gkyl_gk_block_geom_topo(mbinp->gk_block_geom);
  
  int cdim = gkyl_gk_block_geom_ndim(mbapp->gk_block_geom);
  int num_blocks = gkyl_gk_block_geom_num_blocks(mbapp->gk_block_geom);

  // Construct round-robin decomposition.
  int *branks = gkyl_malloc(sizeof(int[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    const struct gkyl_gk_block_geom_info *bgi = gkyl_gk_block_geom_get_block(mbapp->gk_block_geom, i);
    branks[i] = calc_cuts(cdim, bgi->cuts);
  }
  mbapp->round_robin = gkyl_rrobin_decomp_new(num_ranks, num_blocks, branks);

  int num_local_blocks = 0;
  mbapp->local_blocks = gkyl_malloc(sizeof(int[num_blocks]));

  int lidx = 0;
  int *rank_list = gkyl_malloc(sizeof(int[num_ranks])); // This is larger than needed.

  mbapp->decomp = gkyl_malloc(num_blocks*sizeof(struct gkyl_rect_decomp*));
  
  // Construct list of block communicators: there are as many
  // communicators as blocks. Not all communicators are valid on each
  // rank. The total number of valid communicators is
  // num_local_blocks.
  mbapp->block_comms = gkyl_malloc(num_blocks*sizeof(struct gkyl_comm *));
  for (int i=0; i<num_blocks; ++i) {
    gkyl_rrobin_decomp_getranks(mbapp->round_robin, i, rank_list);

    bool is_my_rank_in_decomp = has_int(branks[i], my_rank, rank_list);

    if (is_my_rank_in_decomp) {
      mbapp->local_blocks[lidx++] = i;
      num_local_blocks += 1;      
    }

    const struct gkyl_gk_block_geom_info *bgi = gkyl_gk_block_geom_get_block(mbapp->gk_block_geom, i);
    struct gkyl_range block_global_range;
    gkyl_create_global_range(cdim, bgi->cells, &block_global_range);

    mbapp->decomp[i] = gkyl_rect_decomp_new_from_cuts(
      cdim, bgi->cuts, &block_global_range);

    bool status;
    mbapp->block_comms[i] = gkyl_comm_create_comm_from_ranks(mbinp->comm,
      branks[i], rank_list, mbapp->decomp[i], &status);
  }

  mbapp->num_local_blocks = num_local_blocks;  

  // Create multiblock metadata to pass to each block, and combine with user's metadata.
  const char *fmt_btopo = "%s-block_topo.gkyl";
  int sz = gkyl_calc_strlen(fmt_btopo, mbapp->name);
  char fileNm_btopo[sz+1]; // ensures no buffer overflow
  snprintf(fileNm_btopo, sizeof fileNm_btopo, fmt_btopo, mbapp->name);

  // Basic metadata for I/O (including metadata optional from user).
  const char* build_id = GIT_COMMIT_ID;
  const char* build_date = GKYL_BUILD_DATE;
  struct gkyl_msgpack_map_elem io_meta_default[] = {
    { .key = "changeset", .elem_type = GKYL_MP_STRING, .cval = (char *)build_id },
    { .key = "builddate", .elem_type = GKYL_MP_STRING, .cval = (char *)build_date },
    { .key = "is_multib", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 1 },
    { .key = "topo_file", .elem_type = GKYL_MP_STRING, .cval = fileNm_btopo },
  };
  const struct gkyl_msgpack_map_elem *io_meta_union[] = {io_meta_default, mbinp->metadata.attributes};
  int io_meta_union_len[] = {sizeof(io_meta_default)/sizeof(io_meta_default[0]), mbinp->metadata.num_attributes};

  mbapp->io_meta_basic = gkyl_msgpack_map_elem_union(sizeof(io_meta_union)/sizeof(io_meta_union[0]),
    io_meta_union_len, io_meta_union, &mbapp->io_meta_basic_len);

  printf("Rank %d handles %d Apps\n", my_rank, num_local_blocks);
  for (int i=0; i<num_local_blocks; ++i)
    printf("  Rank %d handles block %d\n", my_rank, mbapp->local_blocks[i]);

  mbapp->num_species = 0;
  mbapp->num_neut_species = 0;
  mbapp->update_field = 0;
  mbapp->singleb_apps = 0;

  if (num_local_blocks > 0) {
    mbapp->num_species = mbinp->num_species;
    mbapp->num_neut_species = mbinp->num_neut_species;
    mbapp->update_field = !mbinp->skip_field; // Note inversion of truth value (default: update field).


    mbapp->singleb_apps = gkyl_malloc(num_local_blocks*sizeof(struct gkyl_gyrokinetic_app*));
  }

  for (int i=0; i<mbinp->num_species; ++i)
    strcpy(mbapp->species_name[i], mbinp->species[i].name);

  for (int i=0; i<mbinp->num_neut_species; ++i)
    strcpy(mbapp->neut_species_name[i], mbinp->neut_species[i].name);  

  // Write the block topo file.
  gkyl_gyrokinetic_multib_app_write_topo(mbapp);

  // Create single-block grids and geometries.
  for (int i=0; i<num_local_blocks; ++i) {
    mbapp->singleb_apps[i] = singleb_app_new_geom(mbinp, mbapp->local_blocks[i], mbapp);
  }

  // Create connections needed for conf-space syncs.
  int ghost[] = { 1, 1, 1 };
  mbapp->mbcc_sync_conf = gkyl_malloc(sizeof(struct gkyl_mbcc_sr));
  mbapp->mbcc_sync_conf->send = gkyl_malloc(mbapp->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
  mbapp->mbcc_sync_conf->recv = gkyl_malloc(mbapp->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));

  for (int bI=0; bI<num_local_blocks; ++bI) {
    int bid = mbapp->local_blocks[bI];

    gkyl_rrobin_decomp_getranks(mbapp->round_robin, bid, rank_list);
    int brank = -1;
    for (int i=0; i<branks[bid]; ++i)
      if (rank_list[i] == my_rank) brank = i;

    mbapp->mbcc_sync_conf->recv[bI] = gkyl_multib_comm_conn_new_recv(bid, brank,
      ghost, &mbapp->block_topo->conn[bid], mbapp->decomp);
    mbapp->mbcc_sync_conf->send[bI] = gkyl_multib_comm_conn_new_send(bid, brank,
      ghost, &mbapp->block_topo->conn[bid], mbapp->decomp);

    struct gkyl_multib_comm_conn *mbcc_s = mbapp->mbcc_sync_conf->send[bI], *mbcc_r = mbapp->mbcc_sync_conf->recv[bI];
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];

    for (int ns=0; ns<mbcc_s->num_comm_conn; ++ns) {
      struct gkyl_comm_conn *ccs = &mbcc_s->comm_conn[ns];
      int rankIdx = ccs->rank;
      gkyl_rrobin_decomp_getranks(mbapp->round_robin, ccs->block_id, rank_list);
      ccs->rank = rank_list[rankIdx];
      gkyl_sub_range_init(&ccs->range, &sbapp->local_ext, ccs->range.lower, ccs->range.upper);
    }
    for (int nr=0; nr<mbcc_r->num_comm_conn; ++nr) {
      struct gkyl_comm_conn *ccr = &mbcc_r->comm_conn[nr];
      int rankIdx = ccr->rank;
      gkyl_rrobin_decomp_getranks(mbapp->round_robin, ccr->block_id, rank_list);
      ccr->rank = rank_list[rankIdx];
      gkyl_sub_range_init(&ccr->range, &sbapp->local_ext, ccr->range.lower, ccr->range.upper);
    }

    gkyl_multib_comm_conn_sort(mbcc_r);
    gkyl_multib_comm_conn_sort(mbcc_s);
  }

  // Sync the conf-space volume Jacobian needed for syncing quantities that include a
  // jacobgeo factor in them.
  struct gkyl_array *jacs_vol[mbapp->num_local_blocks];
  for (int b=0; b<mbapp->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
    jacs_vol[b] = sbapp->gk_geom->geo_int.jacobgeo_ghost;
  }
  gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbapp->num_local_blocks, mbapp->local_blocks,
    mbapp->mbcc_sync_conf->send, mbapp->mbcc_sync_conf->recv, jacs_vol, jacs_vol);

  // Sync the surface conf-space Jacobian, compute its reciprocal, and store its
  // product with the Jacobian of this block in the ghost cell.
  for (int d = 0; d<cdim; d++) {
    struct gkyl_array *jacs[mbapp->num_local_blocks];
    for (int b=0; b<mbapp->num_local_blocks; ++b) {
      struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
      struct gk_geom_surf geo_surf = sbapp->gk_geom->geo_surf[d];
      jacs[b] = geo_surf.jacobgeo_ratio;
      gkyl_array_copy_range(jacs[b], geo_surf.jacobgeo, &sbapp->local_lower_skin[d]);
      gkyl_array_copy_range_to_range(jacs[b], geo_surf.jacobgeo, &sbapp->local_upper_skin[d], &sbapp->local_upper_ghost[d]);
    }
    gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbapp->num_local_blocks, mbapp->local_blocks,
      mbapp->mbcc_sync_conf->send, mbapp->mbcc_sync_conf->recv, jacs, jacs);

    for (int b=0; b<mbapp->num_local_blocks; ++b) {
      struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
      struct gk_geom_surf geo_surf = sbapp->gk_geom->geo_surf[d];
      struct gkyl_array *jacgeo = mkarr(mbapp->use_gpu, geo_surf.jacobgeo_ratio->ncomp, geo_surf.jacobgeo_ratio->size);

      gkyl_array_set_range(jacgeo, 1.0, geo_surf.jacobgeo_ratio, &sbapp->local_lower_ghost[d]);
      gkyl_array_set_range(jacgeo, 1.0, geo_surf.jacobgeo_ratio, &sbapp->local_upper_ghost[d]);
      gkyl_dg_inv_op_range(&sbapp->gk_geom->surf_basis, 0, geo_surf.jacobgeo_ratio, 0, jacgeo, &sbapp->local_lower_ghost[d]);
      gkyl_dg_inv_op_range(&sbapp->gk_geom->surf_basis, 0, geo_surf.jacobgeo_ratio, 0, jacgeo, &sbapp->local_upper_ghost[d]);
      gkyl_array_copy_range_to_range(jacgeo, geo_surf.jacobgeo, &sbapp->local_lower_ghost[d], &sbapp->local_lower_skin[d]);
      gkyl_array_copy_range(jacgeo, geo_surf.jacobgeo, &sbapp->local_upper_ghost[d]);
      gkyl_dg_mul_op_range(&sbapp->gk_geom->surf_basis, 0, geo_surf.jacobgeo_ratio,
        0, jacgeo, 0, geo_surf.jacobgeo_ratio, &sbapp->local_lower_ghost[d]);
      gkyl_dg_mul_op_range(&sbapp->gk_geom->surf_basis, 0, geo_surf.jacobgeo_ratio,
        0, jacgeo, 0, geo_surf.jacobgeo_ratio, &sbapp->local_upper_ghost[d]);
      gkyl_array_clear_range(geo_surf.jacobgeo_ratio, 0.0, &sbapp->local);
      gkyl_array_shiftc_range(geo_surf.jacobgeo_ratio, pow(sqrt(2.0),sbapp->cdim), 0, &sbapp->local);

      gkyl_array_release(jacgeo);
    }
  }

  const struct gkyl_gk_block_geom_info *bgi0 = gkyl_gk_block_geom_get_block(mbapp->gk_block_geom, 0);

  if (cdim > 1 && bgi0->geometry.geometry_id == GKYL_GEOMETRY_TOKAMAK) {
    for (int d = 0; d<cdim; d++) {
      struct gkyl_array *deltats[mbapp->num_local_blocks];
      for (int b=0; b<mbapp->num_local_blocks; ++b) {
        struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
        struct gk_geom_surf geo_surf = sbapp->gk_geom->geo_surf[d];
        deltats[b] = geo_surf.deltats;
        gkyl_array_copy_range_to_range(deltats[b], deltats[b], &sbapp->local_upper_skin[d], &sbapp->local_upper_ghost[d]);
      }
      gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbapp->num_local_blocks, mbapp->local_blocks,
        mbapp->mbcc_sync_conf->send, mbapp->mbcc_sync_conf->recv, deltats, deltats);
    }
  }

  return mbapp;
}

gkyl_gyrokinetic_multib_app* gkyl_gyrokinetic_multib_app_new(const struct gkyl_gyrokinetic_multib *mbinp)
{
  int my_rank, num_ranks;
  gkyl_comm_get_rank(mbinp->comm, &my_rank);
  gkyl_comm_get_size(mbinp->comm, &num_ranks);

  int tot_max[2];
  calc_tot_and_max_cuts(mbinp->gk_block_geom, tot_max);
  if ((num_ranks > tot_max[0]) || (num_ranks < tot_max[1])) {
    fprintf(stderr, "\nSpecified %d total cuts but provided %d processes, \
and the maximum number of cuts in a block is %d\n\n", tot_max[0], num_ranks, tot_max[1]);
    return 0;
  }

  struct gkyl_gyrokinetic_multib_app *mbapp = gkyl_malloc(sizeof(*mbapp));

  strcpy(mbapp->name, mbinp->name);
  mbapp->comm = gkyl_comm_acquire(mbinp->comm);  
  mbapp->use_gpu = mbinp->use_gpu;
  
  mbapp->gk_block_geom = gkyl_gk_block_geom_acquire(mbinp->gk_block_geom);
  mbapp->block_topo = gkyl_gk_block_geom_topo(mbinp->gk_block_geom);
  
  int cdim = gkyl_gk_block_geom_ndim(mbapp->gk_block_geom);
  int num_blocks = gkyl_gk_block_geom_num_blocks(mbapp->gk_block_geom);

  // Construct round-robin decomposition.
  int *branks = gkyl_malloc(sizeof(int[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    const struct gkyl_gk_block_geom_info *bgi = gkyl_gk_block_geom_get_block(mbapp->gk_block_geom, i);
    branks[i] = calc_cuts(cdim, bgi->cuts);
  }
  mbapp->round_robin = gkyl_rrobin_decomp_new(num_ranks, num_blocks, branks);

  int num_local_blocks = 0;
  mbapp->local_blocks = gkyl_malloc(sizeof(int[num_blocks]));

  int lidx = 0;
  int *rank_list = gkyl_malloc(sizeof(int[num_ranks])); // This is larger than needed.

  mbapp->decomp = gkyl_malloc(num_blocks*sizeof(struct gkyl_rect_decomp*));
  
  // Construct list of block communicators: there are as many
  // communicators as blocks. Not all communicators are valid on each
  // rank. The total number of valid communicators is
  // num_local_blocks.
  mbapp->block_comms = gkyl_malloc(num_blocks*sizeof(struct gkyl_comm *));
  for (int i=0; i<num_blocks; ++i) {
    gkyl_rrobin_decomp_getranks(mbapp->round_robin, i, rank_list);

    bool is_my_rank_in_decomp = has_int(branks[i], my_rank, rank_list);

    if (is_my_rank_in_decomp) {
      mbapp->local_blocks[lidx++] = i;
      num_local_blocks += 1;      
    }

    const struct gkyl_gk_block_geom_info *bgi = gkyl_gk_block_geom_get_block(mbapp->gk_block_geom, i);
    struct gkyl_range block_global_range;
    gkyl_create_global_range(cdim, bgi->cells, &block_global_range);

    mbapp->decomp[i] = gkyl_rect_decomp_new_from_cuts(
      cdim, bgi->cuts, &block_global_range);

    bool status;
    mbapp->block_comms[i] = gkyl_comm_create_comm_from_ranks(mbinp->comm,
      branks[i], rank_list, mbapp->decomp[i], &status);
  }

  mbapp->num_local_blocks = num_local_blocks;  

  // Create multiblock metadata to pass to each block, and combine with user's metadata.
  const char *fmt_btopo = "%s-block_topo.gkyl";
  int sz = gkyl_calc_strlen(fmt_btopo, mbapp->name);
  char fileNm_btopo[sz+1]; // ensures no buffer overflow
  snprintf(fileNm_btopo, sizeof fileNm_btopo, fmt_btopo, mbapp->name);

  // Basic metadata for I/O (including metadata optional from user).
  const char* build_id = GIT_COMMIT_ID;
  const char* build_date = GKYL_BUILD_DATE;
  struct gkyl_msgpack_map_elem io_meta_default[] = {
    { .key = "changeset", .elem_type = GKYL_MP_STRING, .cval = (char *)build_id },
    { .key = "builddate", .elem_type = GKYL_MP_STRING, .cval = (char *)build_date },
    { .key = "is_multib", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 1 },
    { .key = "topo_file", .elem_type = GKYL_MP_STRING, .cval = fileNm_btopo },
  };
  const struct gkyl_msgpack_map_elem *io_meta_union[] = {io_meta_default, mbinp->metadata.attributes};
  int io_meta_union_len[] = {sizeof(io_meta_default)/sizeof(io_meta_default[0]), mbinp->metadata.num_attributes};

  mbapp->io_meta_basic = gkyl_msgpack_map_elem_union(sizeof(io_meta_union)/sizeof(io_meta_union[0]),
    io_meta_union_len, io_meta_union, &mbapp->io_meta_basic_len);

  // Write the block topo file.
  gkyl_gyrokinetic_multib_app_write_topo(mbapp);

  printf("Rank %d handles %d Apps\n", my_rank, num_local_blocks);
  for (int i=0; i<num_local_blocks; ++i)
    printf("  Rank %d handles block %d\n", my_rank, mbapp->local_blocks[i]);

  mbapp->num_species = 0;
  mbapp->num_neut_species = 0;
  mbapp->update_field = 0;
  mbapp->singleb_apps = 0;

  if (num_local_blocks > 0) {
    mbapp->num_species = mbinp->num_species;
    mbapp->num_neut_species = mbinp->num_neut_species;
    mbapp->update_field = !mbinp->skip_field; // Note inversion of truth value (default: update field).

    mbapp->singleb_apps = gkyl_malloc(num_local_blocks*sizeof(struct gkyl_gyrokinetic_app*));
  }

  for (int i=0; i<mbinp->num_species; ++i)
    strcpy(mbapp->species_name[i], mbinp->species[i].name);

  for (int i=0; i<mbinp->num_neut_species; ++i)
    strcpy(mbapp->neut_species_name[i], mbinp->neut_species[i].name);  

  // Create single-block grids and geometries.
  for (int i=0; i<num_local_blocks; ++i)
    mbapp->singleb_apps[i] = singleb_app_new_geom(mbinp, mbapp->local_blocks[i], mbapp);

  // Set bmag_ref.
  double bmag_min_local = DBL_MAX;
  double bmag_min_global;
  for (int i=0; i<num_local_blocks; ++i) {
    double bmag_min = gkyl_gk_geometry_reduce_bmag(mbapp->singleb_apps[i]->gk_geom, GKYL_MIN);
    bmag_min_local = GKYL_MIN2(bmag_min_local, bmag_min);
  }
  gkyl_comm_allreduce_host(mbapp->comm, GKYL_DOUBLE, GKYL_MIN, 1, &bmag_min_local, &bmag_min_global);

  double bmag_max_local = -DBL_MAX;
  double bmag_max_global;
  for (int i=0; i<num_local_blocks; ++i) {
    double bmag_max = gkyl_gk_geometry_reduce_bmag(mbapp->singleb_apps[i]->gk_geom, GKYL_MAX);
    bmag_max_local = GKYL_MAX2(bmag_max_local, bmag_max);
  }
  gkyl_comm_allreduce_host(mbapp->comm, GKYL_DOUBLE, GKYL_MAX, 1, &bmag_max_local, &bmag_max_global);

  mbapp->bmag_ref = (bmag_max_global + bmag_min_global)/2.0;

  // Create connections needed for conf-space syncs.
  int ghost[] = { 1, 1, 1 };
  mbapp->mbcc_sync_conf = gkyl_malloc(sizeof(struct gkyl_mbcc_sr));
  mbapp->mbcc_sync_conf->send = gkyl_malloc(mbapp->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
  mbapp->mbcc_sync_conf->recv = gkyl_malloc(mbapp->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));

  for (int bI=0; bI<num_local_blocks; ++bI) {
    int bid = mbapp->local_blocks[bI];

    gkyl_rrobin_decomp_getranks(mbapp->round_robin, bid, rank_list);
    int brank = -1;
    for (int i=0; i<branks[bid]; ++i)
      if (rank_list[i] == my_rank) brank = i;

    mbapp->mbcc_sync_conf->recv[bI] = gkyl_multib_comm_conn_new_recv(bid, brank, 
      ghost, &mbapp->block_topo->conn[bid], mbapp->decomp);
    mbapp->mbcc_sync_conf->send[bI] = gkyl_multib_comm_conn_new_send(bid, brank,
      ghost, &mbapp->block_topo->conn[bid], mbapp->decomp);

    struct gkyl_multib_comm_conn *mbcc_s = mbapp->mbcc_sync_conf->send[bI], *mbcc_r = mbapp->mbcc_sync_conf->recv[bI];
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];

    for (int ns=0; ns<mbcc_s->num_comm_conn; ++ns) {
      // Translate the "rank" in gkyl_multib_comm_conn (right now it is a rank index).
      struct gkyl_comm_conn *ccs = &mbcc_s->comm_conn[ns];
      int rankIdx = ccs->rank;
      gkyl_rrobin_decomp_getranks(mbapp->round_robin, ccs->block_id, rank_list);
      ccs->rank = rank_list[rankIdx];
      // Make range a sub range.
      gkyl_sub_range_init(&ccs->range, &sbapp->local_ext, ccs->range.lower, ccs->range.upper);
    }
    for (int nr=0; nr<mbcc_r->num_comm_conn; ++nr) {
      // Translate the "rank" in gkyl_multib_comm_conn (right now it is a rank index).
      struct gkyl_comm_conn *ccr = &mbcc_r->comm_conn[nr];
      int rankIdx = ccr->rank;
      gkyl_rrobin_decomp_getranks(mbapp->round_robin, ccr->block_id, rank_list);
      ccr->rank = rank_list[rankIdx];
      // Make range a sub range.
      gkyl_sub_range_init(&ccr->range, &sbapp->local_ext, ccr->range.lower, ccr->range.upper);
    }

    // Sort connections according to rank and block ID (needed by NCCL).
    gkyl_multib_comm_conn_sort(mbcc_r);
    gkyl_multib_comm_conn_sort(mbcc_s);
  }

  // Sync the conf-space volume Jacobian needed for syncing quantities that include a
  // jacobgeo factor in them.
  struct gkyl_array *jacs_vol[mbapp->num_local_blocks];
  for (int b=0; b<mbapp->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
    jacs_vol[b] = sbapp->gk_geom->geo_int.jacobgeo_ghost;
  }
  // Sync across blocks.
  gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbapp->num_local_blocks, mbapp->local_blocks,
    mbapp->mbcc_sync_conf->send, mbapp->mbcc_sync_conf->recv, jacs_vol, jacs_vol);

  // Sync the surface conf-space Jacobian, compute its reciprocal, and 
  // store its product with the Jacobian of this block (in the ghost cell).
  for (int d = 0; d<cdim; d++) {
    // Sync jacobgeo.
    struct gkyl_array *jacs[mbapp->num_local_blocks];
    for (int b=0; b<mbapp->num_local_blocks; ++b) {
      struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
      struct gk_geom_surf geo_surf = sbapp->gk_geom->geo_surf[d];
      jacs[b] = geo_surf.jacobgeo_ratio;
      gkyl_array_copy_range(jacs[b], geo_surf.jacobgeo, &sbapp->local_lower_skin[d]);
      gkyl_array_copy_range_to_range(jacs[b], geo_surf.jacobgeo, &sbapp->local_upper_skin[d], &sbapp->local_upper_ghost[d]);
    }
    gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbapp->num_local_blocks, mbapp->local_blocks,
      mbapp->mbcc_sync_conf->send, mbapp->mbcc_sync_conf->recv, jacs, jacs);

    for (int b=0; b<mbapp->num_local_blocks; ++b) {
      struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
      struct gk_geom_surf geo_surf = sbapp->gk_geom->geo_surf[d];
      struct gkyl_array *jacgeo = mkarr(mbapp->use_gpu, geo_surf.jacobgeo_ratio->ncomp, geo_surf.jacobgeo_ratio->size);

      // Compute 1/jacobgeo in ghost cells.
      gkyl_array_set_range(jacgeo, 1.0, geo_surf.jacobgeo_ratio, &sbapp->local_lower_ghost[d]);
      gkyl_array_set_range(jacgeo, 1.0, geo_surf.jacobgeo_ratio, &sbapp->local_upper_ghost[d]);
      gkyl_dg_inv_op_range(&sbapp->gk_geom->surf_basis, 0, geo_surf.jacobgeo_ratio, 0, jacgeo, &sbapp->local_lower_ghost[d]);
      gkyl_dg_inv_op_range(&sbapp->gk_geom->surf_basis, 0, geo_surf.jacobgeo_ratio, 0, jacgeo, &sbapp->local_upper_ghost[d]);
      // Multiply by the Jacobian of this block.
      gkyl_array_copy_range_to_range(jacgeo, geo_surf.jacobgeo, &sbapp->local_lower_ghost[d], &sbapp->local_lower_skin[d]);
      gkyl_array_copy_range(jacgeo, geo_surf.jacobgeo, &sbapp->local_upper_ghost[d]);
      gkyl_dg_mul_op_range(&sbapp->gk_geom->surf_basis, 0, geo_surf.jacobgeo_ratio,
        0, jacgeo, 0, geo_surf.jacobgeo_ratio, &sbapp->local_lower_ghost[d]);
      gkyl_dg_mul_op_range(&sbapp->gk_geom->surf_basis, 0, geo_surf.jacobgeo_ratio,
        0, jacgeo, 0, geo_surf.jacobgeo_ratio, &sbapp->local_upper_ghost[d]);
      // Set the ratio to 1 in the interior (shouldn't be in use).
      gkyl_array_clear_range(geo_surf.jacobgeo_ratio, 0.0, &sbapp->local);
      gkyl_array_shiftc_range(geo_surf.jacobgeo_ratio, pow(sqrt(2.0),sbapp->cdim), 0, &sbapp->local);

      gkyl_array_release(jacgeo);
    }

  }

  const struct gkyl_gk_block_geom_info *bgi0 = gkyl_gk_block_geom_get_block(mbapp->gk_block_geom, 0);

  if (cdim > 1 && bgi0->geometry.geometry_id == GKYL_GEOMETRY_TOKAMAK) {
    // Sync the surface deltats. Need to copy the upper ghost into the upper skin before syncing.
    for (int d = 0; d<cdim; d++) {
      struct gkyl_array *deltats[mbapp->num_local_blocks];
      for (int b=0; b<mbapp->num_local_blocks; ++b) {
        struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
        struct gk_geom_surf geo_surf = sbapp->gk_geom->geo_surf[d];
        deltats[b] = geo_surf.deltats;
        gkyl_array_copy_range_to_range(deltats[b], deltats[b], &sbapp->local_upper_skin[d], &sbapp->local_upper_ghost[d]);
      }
      gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbapp->num_local_blocks, mbapp->local_blocks,
        mbapp->mbcc_sync_conf->send, mbapp->mbcc_sync_conf->recv, deltats, deltats);
    }
    // Accumulate the appropriate shift
    for (int b=0; b<mbapp->num_local_blocks; ++b) {
      int par_dir = cdim-1;
      struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
      struct gkyl_array *delta_ts = sbapp->gk_geom->geo_surf[par_dir].deltats;
      struct gkyl_array *buffer = mkarr(sbapp->use_gpu, delta_ts->ncomp, delta_ts->size);

      gkyl_array_copy_range_to_range(buffer, delta_ts, &sbapp->local_upper_skin[par_dir], &sbapp->local_upper_ghost[par_dir]);
      gkyl_array_accumulate_range(delta_ts, -1.0, buffer, &sbapp->local_upper_skin[par_dir]);

      gkyl_array_copy_range_to_range(buffer, delta_ts, &sbapp->local_lower_skin[par_dir], &sbapp->local_lower_ghost[par_dir]);
      gkyl_array_accumulate_range(delta_ts, -1.0, buffer, &sbapp->local_lower_skin[par_dir]);

      gkyl_array_release(buffer);

      // Deflate delta_ts.
      gyrokinetic_deflate_delta_ts(sbapp, delta_ts);

      // Write the twistshift shift.
      gyrokinetic_app_write_ts_shift(sbapp);
    }
  }


  // Create the rest of the single-block solvers.
  for (int i=0; i<num_local_blocks; ++i)
    singleb_app_new_solver(mbinp, mbapp->local_blocks[i], mbapp, mbapp->singleb_apps[i]);

  // Create the MB field app.
  mbapp->field = gk_multib_field_new(mbinp, mbapp);


  // Create connections needed for syncing charged species phase-space quantities.
  mbapp->mbcc_sync_charged = gkyl_malloc(mbapp->num_species * sizeof(struct gkyl_mbcc_sr));
  for (int i=0; i<mbinp->num_species; ++i) {
    mbapp->mbcc_sync_charged[i].send = gkyl_malloc(mbapp->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
    mbapp->mbcc_sync_charged[i].recv = gkyl_malloc(mbapp->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
  }
  for (int bI=0; bI<num_local_blocks; ++bI) {
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];
    struct gkyl_multib_comm_conn *mbcc_cs = mbapp->mbcc_sync_conf->send[bI], *mbcc_cr = mbapp->mbcc_sync_conf->recv[bI];
    for (int i=0; i<mbinp->num_species; ++i) {
      mbapp->mbcc_sync_charged[i].send[bI] = gkyl_multib_comm_conn_new(mbcc_cs->num_comm_conn, mbcc_cs->comm_conn);
      mbapp->mbcc_sync_charged[i].recv[bI] = gkyl_multib_comm_conn_new(mbcc_cr->num_comm_conn, mbcc_cr->comm_conn);
      struct gkyl_multib_comm_conn *mbcc_ps = mbapp->mbcc_sync_charged[i].send[bI],
                                   *mbcc_pr = mbapp->mbcc_sync_charged[i].recv[bI];
      // Extend ranges to include velocity space, and make them a sub range.
      struct gk_species *gks = &sbapp->species[i];
      for (int ns=0; ns<mbcc_cs->num_comm_conn; ++ns) {
        struct gkyl_comm_conn *ccs_conf = &mbcc_cs->comm_conn[ns];
        struct gkyl_comm_conn *ccs_phase = &mbcc_ps->comm_conn[ns];
        struct gkyl_range phase_r;
        gkyl_range_ten_prod(&phase_r, &ccs_conf->range, &gks->local_vel);
        gkyl_sub_range_init(&ccs_phase->range, &gks->local_ext, phase_r.lower, phase_r.upper);
      }
      for (int nr=0; nr<mbcc_cr->num_comm_conn; ++nr) {
        struct gkyl_comm_conn *ccs_conf = &mbcc_cr->comm_conn[nr];
        struct gkyl_comm_conn *ccs_phase = &mbcc_pr->comm_conn[nr];
        struct gkyl_range phase_r;
        gkyl_range_ten_prod(&phase_r, &ccs_conf->range, &gks->local_vel);
        gkyl_sub_range_init(&ccs_phase->range, &gks->local_ext, phase_r.lower, phase_r.upper);
      }
    }
  }

  // Create connections needed for syncing neutral species phase-space quantities.
  mbapp->mbcc_sync_neut = gkyl_malloc(mbapp->num_neut_species * sizeof(struct gkyl_mbcc_sr));
  for (int i=0; i<mbinp->num_neut_species; ++i) {
    mbapp->mbcc_sync_neut[i].send = gkyl_malloc(mbapp->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
    mbapp->mbcc_sync_neut[i].recv = gkyl_malloc(mbapp->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
  }
  for (int bI=0; bI<num_local_blocks; ++bI) {
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];
    struct gkyl_multib_comm_conn *mbcc_cs = mbapp->mbcc_sync_conf->send[bI], *mbcc_cr = mbapp->mbcc_sync_conf->recv[bI];
    for (int i=0; i<mbinp->num_neut_species; ++i) {
      mbapp->mbcc_sync_neut[i].send[bI] = gkyl_multib_comm_conn_new(mbcc_cs->num_comm_conn, mbcc_cs->comm_conn);
      mbapp->mbcc_sync_neut[i].recv[bI] = gkyl_multib_comm_conn_new(mbcc_cr->num_comm_conn, mbcc_cr->comm_conn);
      struct gkyl_multib_comm_conn *mbcc_ps = mbapp->mbcc_sync_neut[i].send[bI],
                                   *mbcc_pr = mbapp->mbcc_sync_neut[i].recv[bI];
      // Extend ranges to include velocity space, and make them a sub range.
      struct gk_neut_species *gkns = &sbapp->neut_species[i];
      for (int ns=0; ns<mbcc_cs->num_comm_conn; ++ns) {
        struct gkyl_comm_conn *ccs_conf = &mbcc_cs->comm_conn[ns];
        struct gkyl_comm_conn *ccs_phase = &mbcc_ps->comm_conn[ns];
        struct gkyl_range phase_r;
        gkyl_range_ten_prod(&phase_r, &ccs_conf->range, &gkns->local_vel);
        gkyl_sub_range_init(&ccs_phase->range, &gkns->local_ext, phase_r.lower, phase_r.upper);
      }
      for (int nr=0; nr<mbcc_cr->num_comm_conn; ++nr) {
        struct gkyl_comm_conn *ccs_conf = &mbcc_cr->comm_conn[nr];
        struct gkyl_comm_conn *ccs_phase = &mbcc_pr->comm_conn[nr];
        struct gkyl_range phase_r;
        gkyl_range_ten_prod(&phase_r, &ccs_conf->range, &gkns->local_vel);
        gkyl_sub_range_init(&ccs_phase->range, &gkns->local_ext, phase_r.lower, phase_r.upper);
      }
    }
  }

  // Sync the effective diffusivity of the anomalous diffusion operator.
  // Assume they either all have anomalous diffusion or none of them do.
  bool any_anomalous_diff = false;
  for (int i=0; i<mbapp->num_species; ++i) {
    bool has_anomalous_diff = true;
    for (int b=0; b<mbapp->num_local_blocks; ++b) {
      struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
      has_anomalous_diff = has_anomalous_diff && sbapp->species[i].anom_diff.anom_diff_id;
    }

    if (has_anomalous_diff) {
      any_anomalous_diff = true;
      // Divide diffD by dz before the transfer
      for (int b=0; b<mbapp->num_local_blocks; ++b) {
        struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
        int d  = 0;
        gkyl_array_scale_range(sbapp->species[i].anom_diff.diffD, sbapp->grid.dx[sbapp->cdim-1], &sbapp->global);
      }

      // Sync
      struct gkyl_array *gkad_nu[mbapp->num_local_blocks];
      for (int b=0; b<mbapp->num_local_blocks; ++b) {
        struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
        gkad_nu[b] = sbapp->species[i].anom_diff.diffD;
      }
      gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbapp->num_local_blocks, mbapp->local_blocks,
        mbapp->mbcc_sync_conf->send, mbapp->mbcc_sync_conf->recv, gkad_nu, gkad_nu);

      // Multiply by diffD dz after the transfer to achieve rescaling
      for (int b=0; b<mbapp->num_local_blocks; ++b) {
        struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
        int d  = 0;
        gkyl_array_scale_range(sbapp->species[i].anom_diff.diffD, 1.0/sbapp->grid.dx[sbapp->cdim-1], &sbapp->global_ext);
      }

    }
  }

  if (any_anomalous_diff) {
    // Sync the interior conf-space reciprocal Jacobian.
    struct gkyl_array *jacs_inv_vol[mbapp->num_local_blocks];
    for (int b=0; b<mbapp->num_local_blocks; ++b) {
      struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
      jacs_inv_vol[b] = sbapp->gk_geom->geo_int.jacobgeo_inv;
    }
    gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbapp->num_local_blocks, mbapp->local_blocks,
      mbapp->mbcc_sync_conf->send, mbapp->mbcc_sync_conf->recv, jacs_inv_vol, jacs_inv_vol);
  }

  mbapp->stat = (struct gkyl_gyrokinetic_stat) {};

  mbapp->dts = gkyl_dynvec_new(GKYL_DOUBLE, 1); // Dynvector to store time steps.
  mbapp->is_first_dt_write_call = true;

  gkyl_free(rank_list);
  gkyl_free(branks);

  return mbapp;
}

void
gyrokinetic_multib_calc_field(struct gkyl_gyrokinetic_multib_app* app, double tcurr,
  const struct gkyl_array *fin[], struct gkyl_array **bflux[])
{
  struct timespec wtm = gkyl_wall_clock();
  // Compute fields.
  if (app->update_field) {
    // Solve the field equation.
    gk_multib_field_rhs(app, app->field, fin, bflux);
  }
  app->stat.field_tm += gkyl_time_diff_now_sec(wtm);
}

static void
gyrokinetic_multib_apply_bc(struct gkyl_gyrokinetic_multib_app* app, double tcurr,
  struct gkyl_array *distf[], struct gkyl_array *distf_neut[])
{
  // Apply boundary conditions in each block (including intrablock sync).
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    int li_charged = b * app->num_species;
    int li_neut = b * app->num_neut_species;
    for (int i=0; i<app->num_species; ++i) {
      gk_species_apply_bc(sbapp, &sbapp->species[i], distf[li_charged+i]);
    }
    for (int i=0; i<app->num_neut_species; ++i) {
      gk_neut_species_apply_bc(sbapp, &sbapp->neut_species[i], distf_neut[li_neut+i]);
    }
  }

  // Sync blocks.
  struct timespec wst = gkyl_wall_clock();
  int cdim = gkyl_gk_block_geom_ndim(app->gk_block_geom);
  for (int i=0; i<app->num_species; ++i) {
    // Sync charged species.
    struct gkyl_array *fs[app->num_local_blocks];
    for (int b=0; b<app->num_local_blocks; ++b) {
      int li_charged = b * app->num_species;
      fs[b] = distf[li_charged+i];
    }
    gkyl_multib_comm_conn_array_transfer(app->comm, app->num_local_blocks, app->local_blocks,
      app->mbcc_sync_charged[i].send, app->mbcc_sync_charged[i].recv, fs, fs);
  }
  app->stat.species_bc_tm += gkyl_time_diff_now_sec(wst);

  struct timespec wst_neut = gkyl_wall_clock();
  struct gkyl_gyrokinetic_app *sbapp0 = app->singleb_apps[0];
  for (int i=0; i<app->num_neut_species; ++i) {
    // Sync neutral species.
    if (!sbapp0->neut_species[i].info.is_static) {
      struct gkyl_array *fs[app->num_local_blocks];
      for (int b=0; b<app->num_local_blocks; ++b) {
        int li_neut = b * app->num_neut_species;
        fs[b] = distf_neut[li_neut+i];
      }
      gkyl_multib_comm_conn_array_transfer(app->comm, app->num_local_blocks, app->local_blocks,
        app->mbcc_sync_neut[i].send, app->mbcc_sync_neut[i].recv, fs, fs);
    }
  }
  app->stat.neut_species_bc_tm += gkyl_time_diff_now_sec(wst_neut);

}

void
gyrokinetic_multib_calc_field_and_apply_bc(struct gkyl_gyrokinetic_multib_app* app, double tcurr,
  struct gkyl_array *distf[], struct gkyl_array **bflux[], struct gkyl_array *distf_neut[])
{
  // Compute fields and apply BCs.

  // Compute the field.
  // MF 2024/09/27/: Need the cast here for consistency. Fixing
  // this may require removing 'const' from a lot of places.
  gyrokinetic_multib_calc_field(app, tcurr, (const struct gkyl_array **) distf, bflux);

  // Apply boundary conditions.
  struct timespec wst = gkyl_wall_clock();
  gyrokinetic_multib_apply_bc(app, tcurr, distf, distf_neut);
  app->stat.bc_tm += gkyl_time_diff_now_sec(wst);
}

void
gkyl_gyrokinetic_multib_app_apply_ic(gkyl_gyrokinetic_multib_app* app, double t0)
{
  app->tcurr = t0;
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_multib_app_apply_ic_species(app, i, t0);
  }
  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_multib_app_apply_ic_neut_species(app, i, t0);
  }  

  // Compute the fields and apply BCs.
  struct gkyl_array *distf[app->num_species * app->num_local_blocks];
  struct gkyl_array **bflux[app->num_species * app->num_local_blocks];
  struct gkyl_array *distf_neut[app->num_neut_species * app->num_local_blocks];
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    int li_charged = b * app->num_species;
    int li_neut = b * app->num_neut_species;
    for (int i=0; i<app->num_species; ++i) {
      distf[li_charged+i] = sbapp->species[i].f;
      bflux[li_charged+i] = sbapp->species[i].bflux.f1;
    }
    for (int i=0; i<app->num_neut_species; ++i) {
      distf_neut[li_neut+i] = sbapp->neut_species[i].f;
    }
  }
  gyrokinetic_multib_calc_field_and_apply_bc(app, t0, distf, bflux, distf_neut);
}

void
gkyl_gyrokinetic_multib_app_apply_ic_species(gkyl_gyrokinetic_multib_app* app, int sidx, double t0)
{
  app->tcurr = t0;
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_apply_ic_species(app->singleb_apps[b], sidx, t0);
  }
  gkyl_comm_barrier(app->comm);
}

void
gkyl_gyrokinetic_multib_app_apply_ic_neut_species(gkyl_gyrokinetic_multib_app* app, int sidx, double t0)
{
  app->tcurr = t0;
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_apply_ic_neut_species(app->singleb_apps[b], sidx, t0);
  }
  gkyl_comm_barrier(app->comm);
}


struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_file_field(gkyl_gyrokinetic_multib_app *app, const char *fname)
{
 // TO DO
  return (struct gkyl_app_restart_status) { };

}

struct gkyl_app_restart_status 
gkyl_gyrokinetic_multib_app_from_file_species(gkyl_gyrokinetic_multib_app *app, int sidx,
  const char *fname)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}

struct gkyl_app_restart_status 
gkyl_gyrokinetic_multib_app_from_file_neut_species(gkyl_gyrokinetic_multib_app *app, int sidx,
  const char *fname)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_read_from_frame(gkyl_gyrokinetic_multib_app *app, int frame)
{
  struct gkyl_app_restart_status rstat;
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    for (int i=0; i<app->num_neut_species; i++) {
      int neut_frame = frame;
      if (sbapp->neut_species[i].info.is_static) {
        neut_frame = 0;
      }
      rstat = gkyl_gyrokinetic_app_from_frame_neut_species(sbapp, i, neut_frame);
    }
    for (int i=0; i<app->num_species; i++) {
      rstat = gkyl_gyrokinetic_app_from_frame_species(sbapp, i, frame);
    }
  }
  
  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    // Compute the fields and apply BCs.
    struct gkyl_array *distf[app->num_species * app->num_local_blocks];
    struct gkyl_array **bflux[app->num_species * app->num_local_blocks];
    struct gkyl_array *distf_neut[app->num_neut_species * app->num_local_blocks];
    for (int b=0; b<app->num_local_blocks; ++b) {
      struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
      int li_charged = b * app->num_species;
      int li_neut = b * app->num_neut_species;
      for (int i=0; i<app->num_species; ++i) {
        distf[li_charged+i] = sbapp->species[i].f;
        bflux[li_charged+i] = sbapp->species[i].bflux.f;
      }
      for (int i=0; i<app->num_neut_species; ++i) {
        distf_neut[li_neut+i] = sbapp->neut_species[i].f;
      }
    }
//    if (app->update_field && app->field->gkfield_id == GKYL_GK_FIELD_BOLTZMANN) {
//      for (int i=0; i<app->num_species; ++i) {
//        struct gk_species *s = &app->species[i];
//
//        // Compute advection speeds so we can compute the initial boundary flux.
//        gkyl_gk_collisionless_flux_alpha_surf(s->calc_gk_vars, 
//          &app->local, &s->local, &s->local_ext, app->field->phi_smooth,
//          s->alpha_surf, s->sgn_alpha_surf, s->const_sgn_alpha);
//
//        // Compute and store (in the ghost cell of of out) the boundary fluxes.
//        // NOTE: this overwrites ghost cells that may be used for sourcing.
//        gk_species_bflux_rhs(app, s, &s->bflux, distf[i], distf[i]);
//      }
//    }
    gyrokinetic_multib_calc_field_and_apply_bc(app, rstat.stime, distf, bflux, distf_neut);
  }

  struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[0];
  sbapp->field->is_first_energy_write_call = false; // Append to existing diagnostic.
                                                    //
  app->is_first_dt_write_call = false;

  return rstat;
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_frame_field(gkyl_gyrokinetic_multib_app *app, int frame)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_frame_species(gkyl_gyrokinetic_multib_app *app, int sidx, int frame)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_frame_neut_species(gkyl_gyrokinetic_multib_app *app, int sidx, int frame)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}

// private function to handle variable argument list for printing
static void
v_gyrokinetic_multib_app_cout(const gkyl_gyrokinetic_multib_app* app, FILE *fp, const char *fmt, va_list argp)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if ((rank == 0) && fp)
    vfprintf(fp, fmt, argp);
}

void
gkyl_gyrokinetic_multib_app_cout(const gkyl_gyrokinetic_multib_app* app, FILE *fp, const char *fmt, ...)
{
  va_list argp;
  va_start(argp, fmt);
  v_gyrokinetic_multib_app_cout(app, fp, fmt, argp);
  va_end(argp);
}

void
gkyl_gyrokinetic_multib_app_write_topo(const gkyl_gyrokinetic_multib_app* app)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (0 == rank) {
    cstr file_name = cstr_from_fmt("%s-block_topo.gkyl", app->name);
    gkyl_block_topo_write(app->block_topo, file_name.str);
    cstr_drop(&file_name);
  }
}

static struct gkyl_gk_geometry_inp
gyrokinetic_multib_geometry_inp_from_block(const struct gkyl_gk_block_geom_info *bgi,
  struct gkyl_comm *comm)
{
  struct gkyl_gk_geometry_inp geometry_inp = {
    .geometry_id = bgi->geometry.geometry_id,
    .c2p_ctx = bgi->geometry.c2p_ctx,
    .mapc2p = bgi->geometry.mapc2p,
    .bfield_ctx = bgi->geometry.bfield_ctx,
    .bfield_func = bgi->geometry.bfield_func,
    .efit_info = bgi->geometry.efit_info,
    .tok_grid_info = bgi->geometry.tok_grid_info,
    .mirror_grid_info = bgi->geometry.mirror_grid_info,
    .position_map = 0,
    .comm = comm,
    .world = {
      bgi->geometry.world[0],
      bgi->geometry.world[1],
      bgi->geometry.world[2],
    },
    .has_LCFS = bgi->geometry.has_LCFS,
    .x_LCFS = bgi->geometry.x_LCFS,
  };
  strncpy(geometry_inp.geometry_path, bgi->geometry.geometry_path,
    sizeof(geometry_inp.geometry_path)-1);
  geometry_inp.geometry_path[sizeof(geometry_inp.geometry_path)-1] = '\0';
  return geometry_inp;
}

static bool
jacobgeo_ratio_diag_coeff_is_set(double val)
{
  union { double d; uint64_t u; } bits = { .d = val };
  union { double d; uint64_t u; } huge = { .d = 1.0e300 };
  uint64_t abs_bits = bits.u & UINT64_C(0x7fffffffffffffff);
  return abs_bits < huge.u;
}

static void
gyrokinetic_multib_app_write_jacobgeo_ratio_diag(gkyl_gyrokinetic_multib_app *app)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);

  cstr file_name = rank == 0
    ? cstr_from_fmt("%s-jacobgeo_ratio_diagnostics.csv", app->name)
    : cstr_from_fmt("%s-jacobgeo_ratio_diagnostics_rank%d.csv", app->name, rank);
  FILE *fp = fopen(file_name.str, "w");
  cstr_drop(&file_name);
  if (!fp)
    return;

  fprintf(fp, "app,rank,block,dir,edge,num_cells,raw_coeff0_min,raw_coeff0_max,raw_coeff0_mean,ratio_min,ratio_max,ratio_mean\n");
  int cdim = gkyl_gk_block_geom_ndim(app->gk_block_geom);
  const double surf_modal_one = pow(sqrt(2.0), cdim-1);

  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    int bid = app->local_blocks[b];
    for (int d=0; d<cdim; ++d) {
      for (int e=0; e<2; ++e) {
        const char *edge = e == 0 ? "lower_ghost" : "upper_ghost";
        const struct gkyl_range *range = e == 0
          ? &sbapp->local_lower_ghost[d] : &sbapp->local_upper_ghost[d];

        double raw_min = DBL_MAX, raw_max = -DBL_MAX, raw_sum = 0.0;
        double ratio_min = DBL_MAX, ratio_max = -DBL_MAX, ratio_sum = 0.0;
        long count = 0;

        struct gkyl_range_iter iter;
        gkyl_range_iter_init(&iter, range);
        while (gkyl_range_iter_next(&iter)) {
          long loc = gkyl_range_idx(&sbapp->local_ext, iter.idx);
          const double *val = gkyl_array_cfetch(sbapp->gk_geom->geo_surf[d].jacobgeo_ratio, loc);
          double raw = val[0];
          if (!jacobgeo_ratio_diag_coeff_is_set(raw))
            continue;
          double ratio = raw/surf_modal_one;
          raw_min = GKYL_MIN2(raw_min, raw);
          raw_max = GKYL_MAX2(raw_max, raw);
          raw_sum += raw;
          ratio_min = GKYL_MIN2(ratio_min, ratio);
          ratio_max = GKYL_MAX2(ratio_max, ratio);
          ratio_sum += ratio;
          count += 1;
        }

        if (count > 0) {
          fprintf(fp, "%s,%d,%d,%d,%s,%ld,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e\n",
            app->name, rank, bid, d, edge, count,
            raw_min, raw_max, raw_sum/count,
            ratio_min, ratio_max, ratio_sum/count);
        }
      }
    }
  }

  fclose(fp);
}

static void
gyrokinetic_multib_app_write_flux_weight_rescale_diag(gkyl_gyrokinetic_multib_app *app)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);

  cstr file_name = rank == 0
    ? cstr_from_fmt("%s-flux_weight_rescale_diagnostics.csv", app->name)
    : cstr_from_fmt("%s-flux_weight_rescale_diagnostics_rank%d.csv", app->name, rank);
  FILE *fp = fopen(file_name.str, "w");
  cstr_drop(&file_name);
  if (!fp)
    return;

  fprintf(fp, "app,rank,block,dir,edge,quantity,num_cells,coeff0_raw_min,coeff0_raw_max,coeff0_raw_mean,coeff0_ratio_min,coeff0_ratio_max,coeff0_ratio_mean,coeff0_ratio_scaled_min,coeff0_ratio_scaled_max,coeff0_ratio_scaled_mean\n");
  int cdim = gkyl_gk_block_geom_ndim(app->gk_block_geom);
  const double surf_modal_one = pow(sqrt(2.0), cdim-1);

  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    int bid = app->local_blocks[b];
    for (int d=0; d<cdim; ++d) {
      const struct gk_geom_surf *geo_surf = &sbapp->gk_geom->geo_surf[d];
      struct {
        const char *name;
        const struct gkyl_array *array;
      } quantities[] = {
        { "jacobtot_inv", geo_surf->jacobtot_inv },
        { "B3", geo_surf->B3 },
        { "lenr", geo_surf->lenr },
        { "cmag", geo_surf->cmag },
      };

      for (int e=0; e<2; ++e) {
        const char *edge = e == 0 ? "lower_ghost" : "upper_ghost";
        const struct gkyl_range *range = e == 0
          ? &sbapp->local_lower_ghost[d] : &sbapp->local_upper_ghost[d];

        for (int q=0; q<sizeof(quantities)/sizeof(quantities[0]); ++q) {
          int surf_ncomp = geo_surf->jacobgeo_ratio->ncomp;
          int side_offset = quantities[q].array->ncomp == 2*surf_ncomp
            ? (e == 0 ? 0 : surf_ncomp) : 0;
          double raw_min = DBL_MAX, raw_max = -DBL_MAX, raw_sum = 0.0;
          double ratio_min = DBL_MAX, ratio_max = -DBL_MAX, ratio_sum = 0.0;
          double scaled_min = DBL_MAX, scaled_max = -DBL_MAX, scaled_sum = 0.0;
          long count = 0;

          struct gkyl_range_iter iter;
          gkyl_range_iter_init(&iter, range);
          while (gkyl_range_iter_next(&iter)) {
            int skin_idx[GKYL_MAX_DIM];
            gkyl_copy_int_arr(cdim, iter.idx, skin_idx);
            skin_idx[d] += e == 0 ? 1 : -1;
            long skin_loc = gkyl_range_idx(&sbapp->local_ext, skin_idx);
            long ghost_loc = gkyl_range_idx(&sbapp->local_ext, iter.idx);
            const double *raw_val = gkyl_array_cfetch(quantities[q].array, skin_loc);
            const double *ratio_val = gkyl_array_cfetch(geo_surf->jacobgeo_ratio, ghost_loc);
            double raw = raw_val[side_offset];
            double ratio_raw = ratio_val[0];
            if (!jacobgeo_ratio_diag_coeff_is_set(raw) || !jacobgeo_ratio_diag_coeff_is_set(ratio_raw))
              continue;

            double ratio = ratio_raw/surf_modal_one;
            double scaled = raw*ratio;
            raw_min = GKYL_MIN2(raw_min, raw);
            raw_max = GKYL_MAX2(raw_max, raw);
            raw_sum += raw;
            ratio_min = GKYL_MIN2(ratio_min, ratio);
            ratio_max = GKYL_MAX2(ratio_max, ratio);
            ratio_sum += ratio;
            scaled_min = GKYL_MIN2(scaled_min, scaled);
            scaled_max = GKYL_MAX2(scaled_max, scaled);
            scaled_sum += scaled;
            count += 1;
          }

          if (count > 0) {
            fprintf(fp, "%s,%d,%d,%d,%s,%s,%ld,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e\n",
              app->name, rank, bid, d, edge, quantities[q].name, count,
              raw_min, raw_max, raw_sum/count,
              ratio_min, ratio_max, ratio_sum/count,
              scaled_min, scaled_max, scaled_sum/count);
          }
        }
      }
    }
  }

  fclose(fp);
}

static void
gyrokinetic_multib_app_write_flux_weight_modal_product_diag(gkyl_gyrokinetic_multib_app *app)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);

  cstr file_name = rank == 0
    ? cstr_from_fmt("%s-flux_weight_modal_product_diagnostics.csv", app->name)
    : cstr_from_fmt("%s-flux_weight_modal_product_diagnostics_rank%d.csv", app->name, rank);
  FILE *fp = fopen(file_name.str, "w");
  cstr_drop(&file_name);
  if (!fp)
    return;

  fprintf(fp, "app,rank,block,dir,edge,quantity,cell_ord,idx0,idx1,idx2,coeff,side_offset,raw_side_coeff,ratio_coeff,product_coeff\n");
  int cdim = gkyl_gk_block_geom_ndim(app->gk_block_geom);

  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    int bid = app->local_blocks[b];
    for (int d=0; d<cdim; ++d) {
      const struct gk_geom_surf *geo_surf = &sbapp->gk_geom->geo_surf[d];
      int surf_ncomp = geo_surf->jacobgeo_ratio->ncomp;
      struct {
        const char *name;
        const struct gkyl_array *array;
      } quantities[] = {
        { "jacobtot_inv", geo_surf->jacobtot_inv },
        { "B3", geo_surf->B3 },
        { "lenr", geo_surf->lenr },
        { "cmag", geo_surf->cmag },
      };

      for (int e=0; e<2; ++e) {
        const char *edge = e == 0 ? "lower_ghost" : "upper_ghost";
        const struct gkyl_range *range = e == 0
          ? &sbapp->local_lower_ghost[d] : &sbapp->local_upper_ghost[d];

        for (int q=0; q<sizeof(quantities)/sizeof(quantities[0]); ++q) {
          int side_offset = quantities[q].array->ncomp == 2*surf_ncomp
            ? (e == 0 ? 0 : surf_ncomp) : 0;
          struct gkyl_array *raw_side = mkarr(app->use_gpu, surf_ncomp, quantities[q].array->size);
          struct gkyl_array *product = mkarr(app->use_gpu, surf_ncomp, quantities[q].array->size);
          gkyl_array_clear(raw_side, 0.0);
          gkyl_array_clear(product, 0.0);

          struct gkyl_range_iter iter;
          gkyl_range_iter_init(&iter, range);
          while (gkyl_range_iter_next(&iter)) {
            int skin_idx[GKYL_MAX_DIM];
            gkyl_copy_int_arr(cdim, iter.idx, skin_idx);
            skin_idx[d] += e == 0 ? 1 : -1;
            long skin_loc = gkyl_range_idx(&sbapp->local_ext, skin_idx);
            long ghost_loc = gkyl_range_idx(&sbapp->local_ext, iter.idx);
            const double *raw_val = gkyl_array_cfetch(quantities[q].array, skin_loc);
            double *raw_side_val = gkyl_array_fetch(raw_side, ghost_loc);
            for (int k=0; k<surf_ncomp; ++k)
              raw_side_val[k] = raw_val[side_offset+k];
          }

          gkyl_dg_mul_op_range(&sbapp->gk_geom->surf_basis, 0, product,
            0, raw_side, 0, geo_surf->jacobgeo_ratio, range);

          long cell_ord = 0;
          gkyl_range_iter_init(&iter, range);
          while (gkyl_range_iter_next(&iter)) {
            long ghost_loc = gkyl_range_idx(&sbapp->local_ext, iter.idx);
            const double *raw_side_val = gkyl_array_cfetch(raw_side, ghost_loc);
            const double *ratio_val = gkyl_array_cfetch(geo_surf->jacobgeo_ratio, ghost_loc);
            const double *product_val = gkyl_array_cfetch(product, ghost_loc);
            for (int k=0; k<surf_ncomp; ++k) {
              if (!jacobgeo_ratio_diag_coeff_is_set(raw_side_val[k]) ||
                  !jacobgeo_ratio_diag_coeff_is_set(ratio_val[k]) ||
                  !jacobgeo_ratio_diag_coeff_is_set(product_val[k]))
                continue;
              fprintf(fp, "%s,%d,%d,%d,%s,%s,%ld,%d,%d,%d,%d,%d,%.17e,%.17e,%.17e\n",
                app->name, rank, bid, d, edge, quantities[q].name, cell_ord,
                iter.idx[0], iter.idx[1], iter.idx[2], k, side_offset,
                raw_side_val[k], ratio_val[k], product_val[k]);
            }
            cell_ord += 1;
          }

          gkyl_array_release(raw_side);
          gkyl_array_release(product);
        }
      }
    }
  }

  fclose(fp);
}

static const char*
edge_name_from_index(int edge)
{
  return edge == 0 ? "lower" : "upper";
}

static const char*
oriented_edge_name(enum gkyl_oriented_edge edge)
{
  switch (edge) {
  case GKYL_LOWER_POSITIVE: return "lower_positive";
  case GKYL_LOWER_NEGATIVE: return "lower_negative";
  case GKYL_UPPER_POSITIVE: return "upper_positive";
  case GKYL_UPPER_NEGATIVE: return "upper_negative";
  case GKYL_PHYSICAL: return "physical";
  }
  return "unknown";
}

static int
oriented_edge_index(enum gkyl_oriented_edge edge)
{
  switch (edge) {
  case GKYL_LOWER_POSITIVE:
  case GKYL_LOWER_NEGATIVE:
    return 0;
  case GKYL_UPPER_POSITIVE:
  case GKYL_UPPER_NEGATIVE:
    return 1;
  case GKYL_PHYSICAL:
    return -1;
  }
  return -1;
}

static const char*
orientation_name(enum gkyl_oriented_edge edge)
{
  switch (edge) {
  case GKYL_LOWER_POSITIVE:
  case GKYL_UPPER_POSITIVE:
    return "positive";
  case GKYL_LOWER_NEGATIVE:
  case GKYL_UPPER_NEGATIVE:
    return "negative";
  case GKYL_PHYSICAL:
    return "physical";
  }
  return "unknown";
}

enum pointwise_geometry_field {
  POINTWISE_VALID,
  POINTWISE_R,
  POINTWISE_Z,
  POINTWISE_PHI,
  POINTWISE_PSI_REQUESTED,
  POINTWISE_PSI_EVALUATED,
  POINTWISE_PSI_RESIDUAL,
  POINTWISE_EPSI_R,
  POINTWISE_EPSI_PHI,
  POINTWISE_EPSI_Z,
  POINTWISE_ETHETA_R,
  POINTWISE_ETHETA_PHI,
  POINTWISE_ETHETA_Z,
  POINTWISE_JACOBIAN_SIGNED,
  POINTWISE_JACOBIAN_ABS,
  POINTWISE_NUM_FIELDS
};

static int
pointwise_int_pow(int base, int exponent)
{
  int result = 1;
  for (int i=0; i<exponent; ++i)
    result *= base;
  return result;
}

static struct gkyl_array*
pointwise_host_copy(const struct gkyl_array *src)
{
  struct gkyl_array *dest = gkyl_array_new(src->type, src->ncomp, src->size);
  gkyl_array_copy(dest, src);
  return dest;
}

// Evaluate the same equilibrium representation selected by tok_geo.c.  In
// particular, do not use the cubic interpolant when the map was constructed
// from the quadratic/DG representation.
static double
pointwise_eval_psi(const struct gkyl_efit *efit, bool use_cubics,
  double R, double Z)
{
  if (!efit || !jacobgeo_ratio_diag_coeff_is_set(R)
      || !jacobgeo_ratio_diag_coeff_is_set(Z))
    return 0.0;

  if (use_cubics) {
    double xn[2] = { R, Z }, psi = 0.0;
    efit->evf->eval_cubic(0.0, xn, &psi, efit->evf->ctx);
    return psi;
  }

  int idx[2];
  idx[0] = GKYL_MIN2(efit->rzlocal.upper[0], GKYL_MAX2(efit->rzlocal.lower[0],
    efit->rzlocal.lower[0]
      +(int) floor((R-efit->rzgrid.lower[0])/efit->rzgrid.dx[0])));
  idx[1] = GKYL_MIN2(efit->rzlocal.upper[1], GKYL_MAX2(efit->rzlocal.lower[1],
    efit->rzlocal.lower[1]
      +(int) floor((Z-efit->rzgrid.lower[1])/efit->rzgrid.dx[1])));
  long loc = gkyl_range_idx(&efit->rzlocal, idx);
  const double *coeffs = gkyl_array_cfetch(efit->psizr, loc);
  double xc[2], eta[2];
  gkyl_rect_grid_cell_center(&efit->rzgrid, idx, xc);
  eta[0] = (R-xc[0])/(0.5*efit->rzgrid.dx[0]);
  eta[1] = (Z-xc[1])/(0.5*efit->rzgrid.dx[1]);
  return efit->rzbasis.eval_expand(eta, coeffs);
}

// Recover the radial logical coordinate used by gkyl_tok_geo_calc_surface.
// The current surface-node construction is reliable for p1, which is also
// the polynomial order used by the multiblock X-point cases diagnosed here.
static double
pointwise_requested_psi(const struct gkyl_gyrokinetic_app *sbapp,
  int dir, const int *surf_idx)
{
  int poly_order = sbapp->basis.poly_order;
  if (poly_order != 1)
    return 0.0;

  double dx = sbapp->grid.dx[0];
  double logical_lower = sbapp->grid.lower[0]
    +(sbapp->local.lower[0]-sbapp->global.lower[0])*dx;
  double logical_psi;
  if (dir == 0) {
    int node = surf_idx[0]-sbapp->gk_geom->nrange_surf[dir].lower[0];
    logical_psi = logical_lower+node*dx;
  }
  else {
    int node = surf_idx[0]-sbapp->gk_geom->nrange_surf[dir].lower[0];
    int cell = node/2, quad_node = node%2;
    double eta = quad_node == 0 ? -1.0/sqrt(3.0) : 1.0/sqrt(3.0);
    logical_psi = logical_lower+(cell+0.5*(1.0+eta))*dx;
  }

  double requested_psi = 0.0;
  sbapp->position_map->maps[0](0.0, &logical_psi, &requested_psi,
    sbapp->position_map->ctxs[0]);
  return requested_psi;
}

static void
pointwise_fill_geometry_record(const struct gkyl_gyrokinetic_app *sbapp,
  const struct gkyl_efit *efit, bool use_cubics, int dir,
  const int *surf_idx, const struct gkyl_array *mc2p_nodal_fd,
  const struct gkyl_array *dxdz_nodal,
  const struct gkyl_array *jacobgeo_signed_nodal,
  const struct gkyl_array *jacobgeo_nodal, double *record)
{
  long loc = gkyl_range_idx(&sbapp->gk_geom->nrange_surf[dir], surf_idx);
  const double *position = gkyl_array_cfetch(mc2p_nodal_fd, loc);
  const double *tangents = gkyl_array_cfetch(dxdz_nodal, loc);
  const double *jacobian_signed = gkyl_array_cfetch(jacobgeo_signed_nodal, loc);
  const double *jacobian_abs = gkyl_array_cfetch(jacobgeo_nodal, loc);

  double R = position[0], Z = position[1], phi = position[2];
  double cos_phi = cos(phi), sin_phi = sin(phi);
  double requested_psi = pointwise_requested_psi(sbapp, dir, surf_idx);
  double evaluated_psi = pointwise_eval_psi(efit, use_cubics, R, Z);

  record[POINTWISE_R] = R;
  record[POINTWISE_Z] = Z;
  record[POINTWISE_PHI] = phi;
  record[POINTWISE_PSI_REQUESTED] = requested_psi;
  record[POINTWISE_PSI_EVALUATED] = evaluated_psi;
  record[POINTWISE_PSI_RESIDUAL] = evaluated_psi-requested_psi;

  // dxdz_nodal stores Cartesian e_psi in components 0:2 and Cartesian
  // e_theta in components 6:8.  Rotate those physical vectors to the
  // cylindrical basis at the stored phi before writing them.
  record[POINTWISE_EPSI_R] = cos_phi*tangents[0]+sin_phi*tangents[1];
  record[POINTWISE_EPSI_PHI] = -sin_phi*tangents[0]+cos_phi*tangents[1];
  record[POINTWISE_EPSI_Z] = tangents[2];
  record[POINTWISE_ETHETA_R] = cos_phi*tangents[6]+sin_phi*tangents[7];
  record[POINTWISE_ETHETA_PHI] = -sin_phi*tangents[6]+cos_phi*tangents[7];
  record[POINTWISE_ETHETA_Z] = tangents[8];
  record[POINTWISE_JACOBIAN_SIGNED] = jacobian_signed[0];
  record[POINTWISE_JACOBIAN_ABS] = jacobian_abs[0];

  bool valid = sbapp->basis.poly_order == 1;
  for (int f=POINTWISE_R; f<POINTWISE_NUM_FIELDS; ++f)
    valid = valid && jacobgeo_ratio_diag_coeff_is_set(record[f]);
  record[POINTWISE_VALID] = valid ? 1.0 : 0.0;
}

static void
pointwise_fprint_geometry_record(FILE *fp, const double *record)
{
  for (int i=0; i<POINTWISE_NUM_FIELDS; ++i)
    fprintf(fp, ",%.17e", record[i]);
}

static void
pointwise_write_geometry_header(FILE *fp)
{
  fprintf(fp, "app,rank,block,dir,edge,partner_block,partner_dir,partner_edge,orientation,comparison_status,interface_cell,node_ordinal,interface_node_index,cell_idx0,cell_idx1,cell_idx2,surface_idx0,surface_idx1,surface_idx2,quad_idx0,quad_idx1,quad_idx2");
  const char *prefixes[] = { "local", "partner" };
  const char *fields[] = {
    "valid", "R", "Z", "phi", "psi_requested", "psi_evaluated", "psi_residual",
    "e_psi_R", "e_psi_phi", "e_psi_Z", "e_theta_R", "e_theta_phi",
    "e_theta_Z", "jacobian_signed", "jacobian_abs"
  };
  for (int p=0; p<2; ++p)
    for (int f=0; f<POINTWISE_NUM_FIELDS; ++f)
      fprintf(fp, ",%s_%s", prefixes[p], fields[f]);
  fputc('\n', fp);
}

static void
gyrokinetic_multib_app_write_interface_pointwise_diag(
  gkyl_gyrokinetic_multib_app *app)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);

  cstr file_name = rank == 0
    ? cstr_from_fmt("%s-interface_pointwise_geometry_diagnostics.csv", app->name)
    : cstr_from_fmt("%s-interface_pointwise_geometry_diagnostics_rank%d.csv", app->name, rank);
  FILE *fp = fopen(file_name.str, "w");
  if (!fp) {
    int saved_errno = errno;
    fprintf(stderr, "Unable to open pointwise interface geometry diagnostic '%s': %s\n",
      file_name.str, strerror(saved_errno));
    cstr_drop(&file_name);
    return;
  }
  cstr_drop(&file_name);
  pointwise_write_geometry_header(fp);

  int cdim = gkyl_gk_block_geom_ndim(app->gk_block_geom);
  if (cdim < 2 || cdim > 3) {
    fprintf(stderr, "Pointwise interface geometry diagnostic supports cdim=2 or 3; got cdim=%d\n", cdim);
    fclose(fp);
    return;
  }

  struct gkyl_efit *efit[app->num_local_blocks];
  struct gkyl_efit_inp efit_inp[app->num_local_blocks];
  bool use_cubics[app->num_local_blocks];
  for (int b=0; b<app->num_local_blocks; ++b) {
    int bid = app->local_blocks[b];
    const struct gkyl_gk_block_geom_info *bgi =
      gkyl_gk_block_geom_get_block(app->gk_block_geom, bid);
    efit[b] = 0;
    use_cubics[b] = false;
    if (bgi->geometry.geometry_id == GKYL_GEOMETRY_TOKAMAK) {
      efit_inp[b] = bgi->geometry.efit_info;
      efit_inp[b].use_gpu = false;
      efit[b] = gkyl_efit_new(&efit_inp[b]);
      use_cubics[b] = bgi->geometry.tok_grid_info.use_cubics;
    }
  }

  for (int d=0; d<cdim; ++d) {
    int poly_order = app->singleb_apps[0]->basis.poly_order;
    int num_quad = poly_order+1;
    int nodes_per_face = pointwise_int_pow(num_quad, cdim-1);
    int packed_ncomp = 2*nodes_per_face*POINTWISE_NUM_FIELDS;
    struct gkyl_array *side[app->num_local_blocks];
    struct gkyl_array *side_host[app->num_local_blocks];

    for (int b=0; b<app->num_local_blocks; ++b) {
      struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
      side[b] = mkarr(app->use_gpu, packed_ncomp, sbapp->local_ext.volume);
      side_host[b] = app->use_gpu
        ? mkarr(false, packed_ncomp, sbapp->local_ext.volume) : side[b];
      gkyl_array_clear(side_host[b], 0.0);

      if (!efit[b]) {
        if (app->use_gpu)
          gkyl_array_copy(side[b], side_host[b]);
        continue;
      }

      const struct gk_geom_surf *geo_surf = &sbapp->gk_geom->geo_surf[d];
      struct gkyl_array *mc2p_nodal_fd = pointwise_host_copy(geo_surf->mc2p_nodal_fd);
      struct gkyl_array *dxdz_nodal = pointwise_host_copy(geo_surf->dxdz_nodal);
      struct gkyl_array *jacobgeo_signed_nodal =
        pointwise_host_copy(geo_surf->jacobgeo_signed_nodal);
      struct gkyl_array *jacobgeo_nodal = pointwise_host_copy(geo_surf->jacobgeo_nodal);

      for (int e=0; e<2; ++e) {
        bool on_block_edge = e == 0
          ? sbapp->local.lower[d] == sbapp->global.lower[d]
          : sbapp->local.upper[d] == sbapp->global.upper[d];
        if (!on_block_edge)
          continue;

        const struct gkyl_range *skin_range = e == 0
          ? &sbapp->local_lower_skin[d] : &sbapp->local_upper_skin[d];
        struct gkyl_range_iter iter;
        gkyl_range_iter_init(&iter, skin_range);
        while (gkyl_range_iter_next(&iter)) {
          long cell_loc = gkyl_range_idx(&sbapp->local_ext, iter.idx);
          double *packed = gkyl_array_fetch(side_host[b], cell_loc);

          for (int n=0; n<nodes_per_face; ++n) {
            int remainder = n;
            int surf_idx[GKYL_MAX_DIM] = { 0 };
            for (int td=0; td<cdim; ++td) {
              if (td == d) {
                surf_idx[td] = e == 0
                  ? sbapp->gk_geom->nrange_surf[d].lower[td]
                  : sbapp->gk_geom->nrange_surf[d].upper[td];
              }
              else {
                int q = remainder%num_quad;
                remainder /= num_quad;
                surf_idx[td] = sbapp->gk_geom->nrange_surf[d].lower[td]
                  +(iter.idx[td]-sbapp->local.lower[td])*num_quad+q;
              }
            }
            pointwise_fill_geometry_record(sbapp, efit[b], use_cubics[b], d,
              surf_idx, mc2p_nodal_fd, dxdz_nodal,
              jacobgeo_signed_nodal, jacobgeo_nodal,
              &packed[(e*nodes_per_face+n)*POINTWISE_NUM_FIELDS]);
          }
        }
      }

      gkyl_array_release(mc2p_nodal_fd);
      gkyl_array_release(dxdz_nodal);
      gkyl_array_release(jacobgeo_signed_nodal);
      gkyl_array_release(jacobgeo_nodal);
      if (app->use_gpu)
        gkyl_array_copy(side[b], side_host[b]);
    }

    int transfer_status = gkyl_multib_comm_conn_array_transfer(app->comm,
      app->num_local_blocks,
      app->local_blocks, app->mbcc_sync_conf->send, app->mbcc_sync_conf->recv,
      side, side);
    if (transfer_status != 0)
      fprintf(stderr, "Pointwise interface geometry transfer failed in direction %d with status %d\n",
        d, transfer_status);
    if (app->use_gpu)
      for (int b=0; b<app->num_local_blocks; ++b)
        gkyl_array_copy(side_host[b], side[b]);

    for (int b=0; b<app->num_local_blocks; ++b) {
      if (!efit[b])
        continue;

      struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
      int bid = app->local_blocks[b];
      for (int e=0; e<2; ++e) {
        const struct gkyl_target_edge *partner =
          &app->block_topo->conn[bid].connections[d][e];
        if (partner->edge == GKYL_PHYSICAL)
          continue;

        bool on_block_edge = e == 0
          ? sbapp->local.lower[d] == sbapp->global.lower[d]
          : sbapp->local.upper[d] == sbapp->global.upper[d];
        if (!on_block_edge)
          continue;

        int partner_edge = oriented_edge_index(partner->edge);
        bool same_direction = partner->dir == d;
        const char *comparison_status = sbapp->basis.poly_order != 1
          ? "unsupported_poly_order"
          : (same_direction ? "raw_unoriented" : "unsupported_cross_direction");
        const struct gkyl_range *ghost_range = e == 0
          ? &sbapp->local_lower_ghost[d] : &sbapp->local_upper_ghost[d];

        struct gkyl_range_iter iter;
        gkyl_range_iter_init(&iter, ghost_range);
        while (gkyl_range_iter_next(&iter)) {
          int skin_idx[GKYL_MAX_DIM] = { 0 };
          gkyl_copy_int_arr(cdim, iter.idx, skin_idx);
          skin_idx[d] += e == 0 ? 1 : -1;
          long skin_loc = gkyl_range_idx(&sbapp->local_ext, skin_idx);
          long ghost_loc = gkyl_range_idx(&sbapp->local_ext, iter.idx);
          const double *local_packed = gkyl_array_cfetch(side_host[b], skin_loc);
          const double *partner_packed = gkyl_array_cfetch(side_host[b], ghost_loc);

          long interface_cell = 0, cell_stride = 1;
          for (int td=0; td<cdim; ++td) {
            if (td == d)
              continue;
            interface_cell += (skin_idx[td]-sbapp->global.lower[td])*cell_stride;
            cell_stride *= gkyl_range_shape(&sbapp->global, td);
          }

          for (int n=0; n<nodes_per_face; ++n) {
            int remainder = n;
            int surf_idx[GKYL_MAX_DIM] = { -1, -1, -1 };
            int quad_idx[GKYL_MAX_DIM] = { -1, -1, -1 };
            long interface_node = 0, node_stride = 1;
            for (int td=0; td<cdim; ++td) {
              if (td == d) {
                surf_idx[td] = e == 0
                  ? sbapp->gk_geom->nrange_surf[d].lower[td]
                  : sbapp->gk_geom->nrange_surf[d].upper[td];
              }
              else {
                quad_idx[td] = remainder%num_quad;
                remainder /= num_quad;
                surf_idx[td] = sbapp->gk_geom->nrange_surf[d].lower[td]
                  +(skin_idx[td]-sbapp->local.lower[td])*num_quad+quad_idx[td];
                long global_node = (skin_idx[td]-sbapp->global.lower[td])*num_quad
                  +quad_idx[td];
                interface_node += global_node*node_stride;
                node_stride *= gkyl_range_shape(&sbapp->global, td)*num_quad;
              }
            }

            const double *local_record =
              &local_packed[(e*nodes_per_face+n)*POINTWISE_NUM_FIELDS];
            double missing_partner[POINTWISE_NUM_FIELDS];
            for (int f=0; f<POINTWISE_NUM_FIELDS; ++f)
              missing_partner[f] = 0.0;
            const double *partner_record = same_direction && partner_edge >= 0
              ? &partner_packed[(partner_edge*nodes_per_face+n)*POINTWISE_NUM_FIELDS]
              : missing_partner;

            fprintf(fp, "%s,%d,%d,%d,%s,%d,%d,%s,%s,%s,%ld,%d,%ld,%d,%d,%d,%d,%d,%d,%d,%d,%d",
              app->name, rank, bid, d, edge_name_from_index(e), partner->bid,
              partner->dir, oriented_edge_name(partner->edge),
              orientation_name(partner->edge), comparison_status,
              interface_cell, n, interface_node,
              cdim > 0 ? skin_idx[0] : -1,
              cdim > 1 ? skin_idx[1] : -1,
              cdim > 2 ? skin_idx[2] : -1,
              surf_idx[0], surf_idx[1], surf_idx[2],
              quad_idx[0], quad_idx[1], quad_idx[2]);
            pointwise_fprint_geometry_record(fp, local_record);
            pointwise_fprint_geometry_record(fp, partner_record);
            fputc('\n', fp);
          }
        }
      }
    }

    for (int b=0; b<app->num_local_blocks; ++b) {
      if (app->use_gpu)
        gkyl_array_release(side_host[b]);
      gkyl_array_release(side[b]);
    }
  }

  for (int b=0; b<app->num_local_blocks; ++b)
    if (efit[b])
      gkyl_efit_release(efit[b]);
  fclose(fp);
}

static const struct gkyl_array*
gyrokinetic_multib_surf_quantity(const struct gk_geom_surf *geo_surf, const char *name)
{
  if (strcmp(name, "jacobgeo") == 0) return geo_surf->jacobgeo;
  if (strcmp(name, "jacobgeo_signed") == 0)
    return geo_surf->jacobgeo_signed;
  if (strcmp(name, "bmag") == 0) return geo_surf->bmag;
  if (strcmp(name, "jacobtot_inv") == 0) return geo_surf->jacobtot_inv;
  if (strcmp(name, "B3") == 0) return geo_surf->B3;
  if (strcmp(name, "cmag") == 0) return geo_surf->cmag;
  if (strcmp(name, "lenr") == 0) return geo_surf->lenr;
  return 0;
}

static void
gyrokinetic_multib_copy_surf_quantity_skin(gkyl_gyrokinetic_multib_app *app,
  int dir, const char *quantity, struct gkyl_array **side)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    const struct gk_geom_surf *geo_surf = &sbapp->gk_geom->geo_surf[dir];
    const struct gkyl_array *arr = gyrokinetic_multib_surf_quantity(geo_surf, quantity);

    gkyl_array_clear(side[b], 0.0);
    gkyl_array_copy_range(side[b], arr, &sbapp->local_lower_skin[dir]);
    gkyl_array_copy_range_to_range(side[b], arr,
      &sbapp->local_upper_skin[dir], &sbapp->local_upper_ghost[dir]);
  }
}

static void
gyrokinetic_multib_app_write_interface_partner_diag(gkyl_gyrokinetic_multib_app *app)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);

  cstr file_name = rank == 0
    ? cstr_from_fmt("%s-interface_partner_diagnostics.csv", app->name)
    : cstr_from_fmt("%s-interface_partner_diagnostics_rank%d.csv", app->name, rank);
  FILE *fp = fopen(file_name.str, "w");
  cstr_drop(&file_name);
  if (!fp)
    return;

  fprintf(fp, "app,rank,block,dir,edge,partner_block,partner_dir,partner_edge,quantity,num_cells,num_coeff,local_min,local_max,local_mean,partner_min,partner_max,partner_mean,diff_min,diff_max,diff_l1,diff_l2,diff_maxabs,rel_l1,rel_l2\n");

  int cdim = gkyl_gk_block_geom_ndim(app->gk_block_geom);
  const char *quantities[] = { "jacobgeo", "jacobgeo_signed", "bmag", "jacobtot_inv", "B3", "cmag", "lenr" };

  for (int d=0; d<cdim; ++d) {
    for (int q=0; q<sizeof(quantities)/sizeof(quantities[0]); ++q) {
      struct gkyl_array *side[app->num_local_blocks];
      for (int b=0; b<app->num_local_blocks; ++b) {
        struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
        const struct gkyl_array *arr =
          gyrokinetic_multib_surf_quantity(&sbapp->gk_geom->geo_surf[d], quantities[q]);
        side[b] = mkarr(app->use_gpu, arr->ncomp, arr->size);
      }

      gyrokinetic_multib_copy_surf_quantity_skin(app, d, quantities[q], side);
      gkyl_multib_comm_conn_array_transfer(app->comm, app->num_local_blocks, app->local_blocks,
        app->mbcc_sync_conf->send, app->mbcc_sync_conf->recv, side, side);

      for (int b=0; b<app->num_local_blocks; ++b) {
        struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
        int bid = app->local_blocks[b];
        int ncomp = side[b]->ncomp;

        for (int e=0; e<2; ++e) {
          const struct gkyl_target_edge *partner =
            &app->block_topo->conn[bid].connections[d][e];
          if (partner->edge == GKYL_PHYSICAL)
            continue;

          const struct gkyl_range *ghost_range = e == 0
            ? &sbapp->local_lower_ghost[d] : &sbapp->local_upper_ghost[d];

          double local_min = DBL_MAX, local_max = -DBL_MAX, local_sum = 0.0;
          double partner_min = DBL_MAX, partner_max = -DBL_MAX, partner_sum = 0.0;
          double diff_min = DBL_MAX, diff_max = -DBL_MAX, diff_l1 = 0.0, diff_l2 = 0.0, diff_maxabs = 0.0;
          double partner_l1 = 0.0, partner_l2 = 0.0;
          long count = 0, cell_count = 0;

          struct gkyl_range_iter iter;
          gkyl_range_iter_init(&iter, ghost_range);
          while (gkyl_range_iter_next(&iter)) {
            int skin_idx[GKYL_MAX_DIM];
            gkyl_copy_int_arr(cdim, iter.idx, skin_idx);
            skin_idx[d] += e == 0 ? 1 : -1;

            long skin_loc = gkyl_range_idx(&sbapp->local_ext, skin_idx);
            long ghost_loc = gkyl_range_idx(&sbapp->local_ext, iter.idx);
            const double *local_val = gkyl_array_cfetch(side[b], skin_loc);
            const double *partner_val = gkyl_array_cfetch(side[b], ghost_loc);

            for (int k=0; k<ncomp; ++k) {
              double lv = local_val[k], pv = partner_val[k];
              if (!jacobgeo_ratio_diag_coeff_is_set(lv) || !jacobgeo_ratio_diag_coeff_is_set(pv))
                continue;

              double diff = lv-pv;
              local_min = GKYL_MIN2(local_min, lv);
              local_max = GKYL_MAX2(local_max, lv);
              local_sum += lv;
              partner_min = GKYL_MIN2(partner_min, pv);
              partner_max = GKYL_MAX2(partner_max, pv);
              partner_sum += pv;
              diff_min = GKYL_MIN2(diff_min, diff);
              diff_max = GKYL_MAX2(diff_max, diff);
              diff_l1 += fabs(diff);
              diff_l2 += diff*diff;
              diff_maxabs = GKYL_MAX2(diff_maxabs, fabs(diff));
              partner_l1 += fabs(pv);
              partner_l2 += pv*pv;
              count += 1;
            }
            cell_count += 1;
          }

          if (count > 0) {
            double rel_l1 = partner_l1 > 0.0 ? diff_l1/partner_l1 : 0.0;
            double rel_l2 = partner_l2 > 0.0 ? sqrt(diff_l2/partner_l2) : 0.0;
            fprintf(fp, "%s,%d,%d,%d,%s,%d,%d,%s,%s,%ld,%d,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e\n",
              app->name, rank, bid, d, edge_name_from_index(e),
              partner->bid, partner->dir, oriented_edge_name(partner->edge), quantities[q],
              cell_count, ncomp, local_min, local_max, local_sum/count,
              partner_min, partner_max, partner_sum/count,
              diff_min, diff_max, diff_l1, sqrt(diff_l2), diff_maxabs,
              rel_l1, rel_l2);
          }
        }
      }

      for (int b=0; b<app->num_local_blocks; ++b)
        gkyl_array_release(side[b]);
    }
  }

  fclose(fp);
}

void
gkyl_gyrokinetic_multib_app_write_geometry(gkyl_gyrokinetic_multib_app *app)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    int bid = app->local_blocks[b];
    const struct gkyl_gk_block_geom_info *bgi =
      gkyl_gk_block_geom_get_block(app->gk_block_geom, bid);
    struct gkyl_gk_geometry_inp geometry_inp =
      gyrokinetic_multib_geometry_inp_from_block(bgi, app->block_comms[bid]);
    gkyl_gyrokinetic_app_write_geometry(app->singleb_apps[b], &geometry_inp);
  }

  gyrokinetic_multib_app_write_jacobgeo_ratio_diag(app);
  gyrokinetic_multib_app_write_flux_weight_rescale_diag(app);
  gyrokinetic_multib_app_write_flux_weight_modal_product_diag(app);
  gyrokinetic_multib_app_write_interface_partner_diag(app);
  gyrokinetic_multib_app_write_interface_pointwise_diag(app);
}

void
gkyl_gyrokinetic_multib_app_write_field(gkyl_gyrokinetic_multib_app *app, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_field(app->singleb_apps[b], tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_field_energy(gkyl_gyrokinetic_multib_app* app, double tm)
{
  if (app->update_field) {
    for (int b=0; b<app->num_local_blocks; ++b) {
      gkyl_gyrokinetic_app_calc_field_energy(app->singleb_apps[b], tm);
    }
  }
}

void
gkyl_gyrokinetic_multib_app_write_field_energy(gkyl_gyrokinetic_multib_app* app)
{
  if (app->update_field) {
    for (int b=0; b<app->num_local_blocks; ++b) {
      gkyl_gyrokinetic_app_write_field_energy(app->singleb_apps[b]);
    }
  }
}

void
gkyl_gyrokinetic_multib_app_write_eirene(gkyl_gyrokinetic_multib_app *app, double tm, int frame)
{

  for (int b=0; b<app->num_local_blocks; ++b)
    gkyl_gyrokinetic_app_write_eirene_diagnostics(app->singleb_apps[b], tm, frame);

  struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[0];
  cstr fileNm = cstr_from_fmt("%snew_data_flag", sbapp->eirene->info.output_data_path);
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (0 == rank) {
    FILE *fp = fopen(fileNm.str, "w");
    if (fp == NULL)
        return;
    fprintf(fp, "%d\n", frame);
    fclose(fp);
  }
  cstr_drop(&fileNm);
}

void
gkyl_gyrokinetic_multib_app_calc_eirene_integrated_diagnostics(gkyl_gyrokinetic_multib_app *app, double tm)
{
  for (int b=0; b<app->num_local_blocks; ++b)
    gkyl_gyrokinetic_app_calc_eirene_integrated_diagnostics(app->singleb_apps[b], tm);
}

void
gkyl_gyrokinetic_multib_app_write_eirene_integrated_diagnostics(gkyl_gyrokinetic_multib_app *app)
{
  for (int b=0; b<app->num_local_blocks; ++b)
    gkyl_gyrokinetic_app_write_eirene_integrated_diagnostics(app->singleb_apps[b]);
}

//
// ............. Species outputs ............... //
// 
void
gkyl_gyrokinetic_multib_app_write_species(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_neut_species(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_mom(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species_mom(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_species_integrated_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_calc_species_integrated_mom(app->singleb_apps[b], sidx, tm);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_neut_species_integrated_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_calc_neut_species_integrated_mom(app->singleb_apps[b], sidx, tm);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_species_boundary_flux_integrated_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_calc_species_boundary_flux_integrated_mom(app->singleb_apps[b], sidx, tm);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_integrated_mom(app->singleb_apps[b], sidx);
  }
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species_integrated_mom(app->singleb_apps[b], sidx);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_boundary_flux_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_boundary_flux_integrated_mom(app->singleb_apps[b], sidx);
  }
}

//
// ............. Source outputs ............... //
// 
void
gkyl_gyrokinetic_multib_app_write_species_source(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_source(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_source(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species_source(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_source_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_source_mom(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_source_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species_source_mom(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_species_source_integrated_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_calc_species_source_integrated_mom(app->singleb_apps[b], sidx, tm);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_neut_species_source_integrated_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_calc_neut_species_source_integrated_mom(app->singleb_apps[b], sidx, tm);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_source_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_source_integrated_mom(app->singleb_apps[b], sidx);
  }
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_source_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species_source_integrated_mom(app->singleb_apps[b], sidx);
  }
}

//
// ............. BGK Source outputs ............... //
// 
void
gkyl_gyrokinetic_multib_app_calc_species_source_bgk_integrated_diagnostics(gkyl_gyrokinetic_multib_app* app, int sidx, double tm)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_calc_species_source_bgk_integrated_diagnostics(app->singleb_apps[b], sidx, tm);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_source_bgk_integrated_diagnostics(gkyl_gyrokinetic_multib_app *app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_source_bgk_integrated_diagnostics(app->singleb_apps[b], sidx);
  }
}

//
// ............. LTE outputs ............... //
// 
void
gkyl_gyrokinetic_multib_app_write_species_lte_max_corr_status(gkyl_gyrokinetic_multib_app* app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_lte_max_corr_status(app->singleb_apps[b], sidx);
  }
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_lte_max_corr_status(gkyl_gyrokinetic_multib_app* app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species_lte_max_corr_status(app->singleb_apps[b], sidx);
  }
}

//
// ............. Collision outputs ............... //
// 
void
gkyl_gyrokinetic_multib_app_write_species_lbo_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_lbo_mom(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_bgk_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_bgk_mom(app->singleb_apps[b], sidx, tm, frame);
  }
}

//
// ............. Radiation outputs ............... //
// 
void
gkyl_gyrokinetic_multib_app_write_species_rad_drag(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_rad_drag(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_rad_emissivity(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_rad_emissivity(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_species_rad_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx, double tm)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_calc_species_rad_integrated_mom(app->singleb_apps[b], sidx, tm);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_rad_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_rad_integrated_mom(app->singleb_apps[b], sidx);
  }
}

//
// ............. Functions that group several outputs for a single species ............... //
//
void
gkyl_gyrokinetic_multib_app_write_species_phase(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  gkyl_gyrokinetic_multib_app_write_species(app, sidx, tm, frame);

  gkyl_gyrokinetic_multib_app_write_species_source(app, sidx, tm, frame);

  gkyl_gyrokinetic_multib_app_write_species_rad_drag(app, sidx, tm, frame);
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_phase(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  gkyl_gyrokinetic_multib_app_write_neut_species(app, sidx, tm, frame);

  gkyl_gyrokinetic_multib_app_write_neut_species_source(app, sidx, tm, frame);
}

void
gkyl_gyrokinetic_multib_app_write_species_conf(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_conf(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_conf(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species_conf(app->singleb_apps[b], sidx, tm, frame);
  }
}

//
// ............. Functions that group several species outputs ............... //
// 
void
gkyl_gyrokinetic_multib_app_write_mom(gkyl_gyrokinetic_multib_app* app, double tm, int frame)
{
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_species_mom(app, i, tm, frame);
    gkyl_gyrokinetic_multib_app_write_species_source_mom(app, i, tm, frame);
    gkyl_gyrokinetic_multib_app_write_species_lbo_mom(app, i, tm, frame);
    gkyl_gyrokinetic_multib_app_write_species_bgk_mom(app, i, tm, frame);
    gkyl_gyrokinetic_multib_app_write_species_rad_emissivity(app, i, tm, frame);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_neut_species_mom(app, i, tm, frame);
    gkyl_gyrokinetic_multib_app_write_neut_species_source_mom(app, i, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_integrated_mom(gkyl_gyrokinetic_multib_app* app, double tm)
{
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_multib_app_calc_species_integrated_mom(app, i, tm);
    gkyl_gyrokinetic_multib_app_calc_species_source_integrated_mom(app, i, tm);
    gkyl_gyrokinetic_multib_app_calc_species_source_bgk_integrated_diagnostics(app, i, tm);
    gkyl_gyrokinetic_multib_app_calc_species_rad_integrated_mom(app, i, tm);
    gkyl_gyrokinetic_multib_app_calc_species_boundary_flux_integrated_mom(app, i, tm);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_multib_app_calc_neut_species_integrated_mom(app, i, tm);
    gkyl_gyrokinetic_multib_app_calc_neut_species_source_integrated_mom(app, i, tm);
  }

  gkyl_gyrokinetic_multib_app_calc_eirene_integrated_diagnostics(app, tm);
}

void
gkyl_gyrokinetic_multib_app_write_integrated_mom(gkyl_gyrokinetic_multib_app *app)
{
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_species_integrated_mom(app, i);
    gkyl_gyrokinetic_multib_app_write_species_source_integrated_mom(app, i);
    gkyl_gyrokinetic_multib_app_write_species_source_bgk_integrated_diagnostics(app, i);
    gkyl_gyrokinetic_multib_app_write_species_lte_max_corr_status(app, i);
    gkyl_gyrokinetic_multib_app_write_species_rad_integrated_mom(app, i);
    gkyl_gyrokinetic_multib_app_write_species_boundary_flux_integrated_mom(app, i);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_neut_species_integrated_mom(app, i);
    gkyl_gyrokinetic_multib_app_write_neut_species_source_integrated_mom(app, i);
    gkyl_gyrokinetic_multib_app_write_neut_species_lte_max_corr_status(app, i);
  }

  gkyl_gyrokinetic_multib_app_write_eirene_integrated_diagnostics(app);
}

void
gkyl_gyrokinetic_multib_app_write_conf(gkyl_gyrokinetic_multib_app* app, double tm, int frame)
{
  gkyl_gyrokinetic_multib_app_write_field(app, tm, frame);

  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_species_conf(app, i, tm, frame);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_neut_species_conf(app, i, tm, frame);
  }

  gkyl_gyrokinetic_multib_app_write_eirene(app, tm, frame);
}

void
gkyl_gyrokinetic_multib_app_write_phase(gkyl_gyrokinetic_multib_app* app, double tm, int frame)
{
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_species_phase(app, i, tm, frame);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_neut_species_phase(app, i, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write(gkyl_gyrokinetic_multib_app* app, double tm, int frame)
{
  gkyl_gyrokinetic_multib_app_write_phase(app, tm, frame);

  gkyl_gyrokinetic_multib_app_write_conf(app, tm, frame);
}

//
// ............. End of write functions ............... //
// 

struct gkyl_update_status
gkyl_gyrokinetic_multib_update(gkyl_gyrokinetic_multib_app* app, double dt)
{
  app->stat.nup += 1;
  struct timespec wst = gkyl_wall_clock();

  struct gkyl_update_status status = gyrokinetic_multib_update_ssp_rk3(app, dt);
  app->tcurr += status.dt_actual;

  app->stat.time_loop_tm += gkyl_time_diff_now_sec(wst);

  // Check for any CUDA errors during time step
  if (app->use_gpu)
    checkCuda(cudaGetLastError());

  return status;
}

struct gkyl_gyrokinetic_stat
gkyl_gyrokinetic_multib_app_stat(gkyl_gyrokinetic_multib_app* app)
{
  for (int i=0; i<app->num_species; ++i) {
    app->stat.n_iter_corr[i] = 0;
    app->stat.num_corr[i] = 0;
  }
  for (int i=0; i<app->num_neut_species; ++i) {
    app->stat.neut_n_iter_corr[i] = 0;
    app->stat.neut_num_corr[i] = 0;
  }

  for (int b=0; b<app->num_local_blocks; ++b) {
    // Add time spent on various operations for each local block.
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    struct gkyl_gyrokinetic_stat sb_stat = gkyl_gyrokinetic_app_stat(sbapp);

    app->stat.field_phi_rhs_tm += sb_stat.field_phi_rhs_tm;

    app->stat.species_lte_tm += sb_stat.species_lte_tm;
    app->stat.species_coll_mom_tm += sb_stat.species_coll_mom_tm;
    app->stat.species_coll_tm += sb_stat.species_coll_tm;
    app->stat.species_rad_mom_tm += sb_stat.species_rad_mom_tm;
    app->stat.species_rad_tm += sb_stat.species_rad_tm;
    app->stat.species_react_mom_tm += sb_stat.species_react_mom_tm;
    app->stat.species_react_tm += sb_stat.species_react_tm;
    app->stat.species_omega_cfl_tm += sb_stat.species_omega_cfl_tm;

    app->stat.neut_species_lte_tm += sb_stat.neut_species_lte_tm;
    app->stat.neut_species_coll_mom_tm += sb_stat.neut_species_coll_mom_tm;
    app->stat.neut_species_coll_tm += sb_stat.neut_species_coll_tm;
    app->stat.neut_species_react_mom_tm += sb_stat.neut_species_react_mom_tm;
    app->stat.neut_species_react_tm += sb_stat.neut_species_react_tm;
    app->stat.neut_species_omega_cfl_tm += sb_stat.neut_species_omega_cfl_tm;

    app->stat.species_bc_tm += sb_stat.species_bc_tm;
    app->stat.neut_species_bc_tm += sb_stat.neut_species_bc_tm;

    app->stat.species_io_tm += sb_stat.species_io_tm;
    app->stat.species_diag_calc_tm += sb_stat.species_diag_calc_tm;
    app->stat.species_diag_io_tm += sb_stat.species_diag_io_tm;

    app->stat.neut_species_io_tm += sb_stat.neut_species_io_tm;
    app->stat.neut_species_diag_calc_tm += sb_stat.neut_species_diag_calc_tm;
    app->stat.neut_species_diag_io_tm += sb_stat.neut_species_diag_io_tm;

    app->stat.n_species_omega_cfl += sb_stat.n_species_omega_cfl;
    app->stat.n_mom += sb_stat.n_mom;    
    app->stat.n_diag += sb_stat.n_diag;
    app->stat.n_io += sb_stat.n_io;
    app->stat.n_diag_io += sb_stat.n_diag_io;

    app->stat.n_neut_species_omega_cfl += sb_stat.n_neut_species_omega_cfl;
    app->stat.n_neut_mom += sb_stat.n_neut_mom;  
    app->stat.n_neut_diag += sb_stat.n_neut_diag;
    app->stat.n_neut_io += sb_stat.n_neut_io;
    app->stat.n_neut_diag_io += sb_stat.n_neut_diag_io;

    for (int i=0; i<app->num_species; ++i) {
      app->stat.n_iter_corr[i] += sb_stat.n_iter_corr[i];
      app->stat.num_corr[i] += sb_stat.num_corr[i];
    }
    for (int i=0; i<app->num_neut_species; ++i) {
      app->stat.neut_n_iter_corr[i] += sb_stat.neut_n_iter_corr[i];
      app->stat.neut_num_corr[i] += sb_stat.neut_num_corr[i];
    }

    // Timers not yet computed in app directly.
    app->stat.time_rate_diags_tm += sb_stat.fdot_tm + sb_stat.phidot_tm;
    app->stat.pos_shift_tm += sb_stat.species_pos_shift_tm + sb_stat.neut_species_pos_shift_tm + sb_stat.pos_shift_quasineut_tm;
    app->stat.io_tm += sb_stat.species_io_tm + sb_stat.species_diag_calc_tm + sb_stat.species_diag_io_tm + sb_stat.neut_species_io_tm
      + sb_stat.neut_species_diag_calc_tm + sb_stat.neut_species_diag_io_tm + sb_stat.field_io_tm + sb_stat.field_diag_calc_tm
      + sb_stat.field_diag_io_tm + sb_stat.app_io_tm;
    
    // Additions of several timers.
    app->stat.fwd_euler_sum_tm += sb_stat.species_coll_mom_tm + sb_stat.species_react_mom_tm + sb_stat.neut_species_coll_mom_tm
      + sb_stat.neut_species_react_mom_tm + sb_stat.species_rad_mom_tm + sb_stat.species_gyroavg_tm + sb_stat.species_collisionless_tm
      + sb_stat.species_coll_tm + sb_stat.species_diffusion_tm + sb_stat.species_rad_tm + sb_stat.species_react_tm
      + sb_stat.species_bflux_calc_tm+sb_stat.species_bflux_moms_tm + sb_stat.species_omega_cfl_tm + sb_stat.species_src_tm
      + sb_stat.neut_species_collisionless_tm + sb_stat.neut_species_coll_tm + sb_stat.neut_species_react_tm
      + sb_stat.neut_species_omega_cfl_tm + sb_stat.neut_species_src_tm + sb_stat.dfdt_dt_reduce_tm;

    app->stat.time_rate_diags_sum_tm += sb_stat.fdot_tm + sb_stat.phidot_tm;

    app->stat.pos_shift_sum_tm += sb_stat.species_pos_shift_tm + sb_stat.neut_species_pos_shift_tm + sb_stat.pos_shift_quasineut_tm;

    app->stat.io_sum_tm += sb_stat.species_io_tm + sb_stat.species_diag_calc_tm + sb_stat.species_diag_io_tm + sb_stat.neut_species_io_tm
    + sb_stat.neut_species_diag_calc_tm + sb_stat.neut_species_diag_io_tm + sb_stat.field_io_tm + sb_stat.field_diag_calc_tm
    + sb_stat.field_diag_io_tm + sb_stat.app_io_tm;

  }

  app->stat.fwd_euler_sum_tm += app->stat.fwd_euler_step_f_tm + app->stat.dfdt_dt_reduce_tm;
  app->stat.bc_sum_tm = app->stat.species_bc_tm + app->stat.neut_species_bc_tm;
  app->stat.field_sum_tm = app->stat.field_phi_rhs_tm  + app->stat.field_phi_solve_tm;
  app->stat.time_stepper_sum_tm = app->stat.fwd_euler_tm  + app->stat.field_tm + app->stat.bc_tm
    + app->stat.time_rate_diags_tm + app->stat.pos_shift_tm + app->stat.time_stepper_arithmetic_tm;

  return app->stat;
}

static inline
double
ratio_to_percent(double num, double den, double alt)
{
  return den > 1e-12 ? 100.*num/den : alt;
}

void
gkyl_gyrokinetic_multib_app_print_timings(gkyl_gyrokinetic_multib_app* app, FILE *iostream)
{
  struct gkyl_gyrokinetic_stat *stat = &app->stat;

  double bflux_tm = stat->species_bflux_calc_tm+stat->species_bflux_moms_tm;

  gkyl_gyrokinetic_multib_app_cout(app, iostream, "Timing:\n");
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "  - Time loop:                         %.4e sec.\n", stat->time_loop_tm);
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Forward Euler:                   %.4e sec. / %4.2f %%.\n", stat->fwd_euler_tm, ratio_to_percent(stat->fwd_euler_tm,stat->time_loop_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Collision moments (charged):   %.4e sec. / %4.2f %%.\n", stat->species_coll_mom_tm          , ratio_to_percent(stat->species_coll_mom_tm          ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Reaction moments (charged):    %.4e sec. / %4.2f %%.\n", stat->species_react_mom_tm         , ratio_to_percent(stat->species_react_mom_tm         ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Collision moments (neutral):   %.4e sec. / %4.2f %%.\n", stat->neut_species_coll_mom_tm     , ratio_to_percent(stat->neut_species_coll_mom_tm     ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Reaction moments (neutral):    %.4e sec. / %4.2f %%.\n", stat->neut_species_react_mom_tm    , ratio_to_percent(stat->neut_species_react_mom_tm    ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Radiation moments:             %.4e sec. / %4.2f %%.\n", stat->species_rad_mom_tm           , ratio_to_percent(stat->species_rad_mom_tm           ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Species gyroaverage:           %.4e sec. / %4.2f %%.\n", stat->species_gyroavg_tm           , ratio_to_percent(stat->species_gyroavg_tm           ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Collisionless terms (charged): %.4e sec. / %4.2f %%.\n", stat->species_collisionless_tm     , ratio_to_percent(stat->species_collisionless_tm     ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Collision terms (charged):     %.4e sec. / %4.2f %%.\n", stat->species_coll_tm              , ratio_to_percent(stat->species_coll_tm              ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Diffusion (charged):           %.4e sec. / %4.2f %%.\n", stat->species_diffusion_tm         , ratio_to_percent(stat->species_diffusion_tm         ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Radiation terms:               %.4e sec. / %4.2f %%.\n", stat->species_rad_tm               , ratio_to_percent(stat->species_rad_tm               ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Reaction terms (charged):      %.4e sec. / %4.2f %%.\n", stat->species_react_tm             , ratio_to_percent(stat->species_react_tm             ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Boundary fluxes (charged):     %.4e sec. / %4.2f %%.\n", bflux_tm                           , ratio_to_percent(bflux_tm                           ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ omega_cfl (charged):           %.4e sec. / %4.2f %%.\n", stat->species_omega_cfl_tm         , ratio_to_percent(stat->species_omega_cfl_tm         ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Sources (charged):             %.4e sec. / %4.2f %%.\n", stat->species_src_tm               , ratio_to_percent(stat->species_src_tm               ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Collisionless terms (neutral): %.4e sec. / %4.2f %%.\n", stat->neut_species_collisionless_tm, ratio_to_percent(stat->neut_species_collisionless_tm,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Collision terms (neutral):     %.4e sec. / %4.2f %%.\n", stat->neut_species_coll_tm         , ratio_to_percent(stat->neut_species_coll_tm         ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Reaction terms (neutral):      %.4e sec. / %4.2f %%.\n", stat->neut_species_react_tm        , ratio_to_percent(stat->neut_species_react_tm        ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ omega_cfl (neutral):           %.4e sec. / %4.2f %%.\n", stat->neut_species_omega_cfl_tm    , ratio_to_percent(stat->neut_species_omega_cfl_tm    ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Sources (neutral):             %.4e sec. / %4.2f %%.\n", stat->neut_species_src_tm          , ratio_to_percent(stat->neut_species_src_tm          ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Time step reduction:           %.4e sec. / %4.2f %%.\n", stat->dfdt_dt_reduce_tm            , ratio_to_percent(stat->dfdt_dt_reduce_tm            ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Step f:                        %.4e sec. / %4.2f %%.\n", stat->fwd_euler_step_f_tm          , ratio_to_percent(stat->fwd_euler_step_f_tm          ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Accounted for:                 %4.2f %%.\n", ratio_to_percent(stat->fwd_euler_sum_tm, stat->fwd_euler_tm, 100.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Field solves:                    %.4e sec. / %4.2f %%.\n", stat->field_tm, ratio_to_percent(stat->field_tm,stat->time_loop_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Phi eqn RHS:                   %.4e sec. / %4.2f %%.\n", stat->field_phi_rhs_tm  , ratio_to_percent(stat->field_phi_rhs_tm  ,stat->field_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Phi eqn solve:                 %.4e sec. / %4.2f %%.\n", stat->field_phi_solve_tm, ratio_to_percent(stat->field_phi_solve_tm,stat->field_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Accounted for:                 %4.2f %%.\n", ratio_to_percent(stat->field_sum_tm, stat->field_tm, 100.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Boundary conditions::            %.4e sec. / %4.2f %%.\n", stat->bc_tm, ratio_to_percent(stat->bc_tm,stat->time_loop_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Species (charged):             %.4e sec. / %4.2f %%.\n", stat->species_bc_tm     , ratio_to_percent(stat->species_bc_tm     ,stat->bc_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Species (neutral):             %.4e sec. / %4.2f %%.\n", stat->neut_species_bc_tm, ratio_to_percent(stat->neut_species_bc_tm,stat->bc_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Accounted for:                 %4.2f %%.\n", ratio_to_percent(stat->bc_sum_tm, stat->bc_tm, 100.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Time rate diagnostics:           %.4e sec. / %4.2f %%.\n", stat->time_rate_diags_tm, ratio_to_percent(stat->time_rate_diags_tm,stat->time_loop_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Charged species:               %.4e sec. / %4.2f %%.\n", stat->fdot_tm  , ratio_to_percent(stat->fdot_tm  ,stat->time_rate_diags_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Phi:                           %.4e sec. / %4.2f %%.\n", stat->phidot_tm, ratio_to_percent(stat->phidot_tm,stat->time_rate_diags_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Accounted for:                 %4.2f %%.\n", ratio_to_percent(stat->time_rate_diags_sum_tm, stat->time_rate_diags_tm, 100.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Positivity shift:                %.4e sec. / %4.2f %%.\n", stat->pos_shift_tm, ratio_to_percent(stat->pos_shift_tm,stat->time_loop_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Species (charged):             %.4e sec. / %4.2f %%.\n", stat->species_pos_shift_tm     , ratio_to_percent(stat->species_pos_shift_tm     ,stat->pos_shift_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Species (neutral):             %.4e sec. / %4.2f %%.\n", stat->neut_species_pos_shift_tm, ratio_to_percent(stat->neut_species_pos_shift_tm,stat->pos_shift_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Quasineutrality:               %.4e sec. / %4.2f %%.\n", stat->pos_shift_quasineut_tm   , ratio_to_percent(stat->pos_shift_quasineut_tm   ,stat->pos_shift_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Accounted for:                 %4.2f %%.\n", ratio_to_percent(stat->pos_shift_sum_tm, stat->pos_shift_tm, 100.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Time stepper arithmetic:         %.4e sec. / %4.2f %%.\n", stat->time_stepper_arithmetic_tm, ratio_to_percent(stat->time_stepper_arithmetic_tm,stat->time_loop_tm, 0.0));

  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Accounted for:                   %4.2f %%.\n", ratio_to_percent(stat->time_stepper_sum_tm, stat->time_loop_tm, 100.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "  - I/O:                               %.4e sec.\n", stat->io_tm);
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * f write (charged):               %.4e sec. / %4.2f %%.\n", stat->species_io_tm            , ratio_to_percent(stat->species_io_tm            , stat->io_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Species diag calc (charged):     %.4e sec. / %4.2f %%.\n", stat->species_diag_calc_tm     , ratio_to_percent(stat->species_diag_calc_tm     , stat->io_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Species diag write (charged):    %.4e sec. / %4.2f %%.\n", stat->species_diag_io_tm       , ratio_to_percent(stat->species_diag_io_tm       , stat->io_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * f write (neutral):               %.4e sec. / %4.2f %%.\n", stat->neut_species_io_tm       , ratio_to_percent(stat->neut_species_io_tm       , stat->io_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Species diag calc (neutral):     %.4e sec. / %4.2f %%.\n", stat->neut_species_diag_calc_tm, ratio_to_percent(stat->neut_species_diag_calc_tm, stat->io_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Species diag write (neutral):    %.4e sec. / %4.2f %%.\n", stat->neut_species_diag_io_tm  , ratio_to_percent(stat->neut_species_diag_io_tm  , stat->io_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Field write:                     %.4e sec. / %4.2f %%.\n", stat->field_io_tm              , ratio_to_percent(stat->field_io_tm              , stat->io_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Field diag calc:                 %.4e sec. / %4.2f %%.\n", stat->field_diag_calc_tm       , ratio_to_percent(stat->field_diag_calc_tm       , stat->io_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Field diag write:                %.4e sec. / %4.2f %%.\n", stat->field_diag_io_tm         , ratio_to_percent(stat->field_diag_io_tm         , stat->io_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Common write:                    %.4e sec. / %4.2f %%.\n", stat->app_io_tm                , ratio_to_percent(stat->app_io_tm                , stat->io_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Accounted for:                   %4.2f %%.\n", ratio_to_percent(stat->io_sum_tm, stat->io_tm, 100.0));
}

void
gkyl_gyrokinetic_multib_app_stat_write(gkyl_gyrokinetic_multib_app* app)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    gkyl_gyrokinetic_app_stat_write(sbapp);
  }
}

void
gkyl_gyrokinetic_multib_app_write_dt(gkyl_gyrokinetic_multib_app* app)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);

  if (rank == 0) {
    // Write integrated diagnostic moments.
    const char *fmt = "%s-%s.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, "dt");
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, "dt");

    struct timespec wtm = gkyl_wall_clock();
    if (app->is_first_dt_write_call) {
      struct gkyl_msgpack_map_elem io_meta_phi[] = {
        { .key = "Description", .elem_type = GKYL_MP_STRING, .cval = "Time step size." }
      };
      int io_meta_len[] = {app->io_meta_basic_len, 1};
      const struct gkyl_msgpack_map_elem* io_meta[] = {app->io_meta_basic, io_meta_phi};
      struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

      gkyl_dynvec_write_wmeta(app->dts, fileNm, mt);
      app->is_first_dt_write_call = false;
      gkyl_msgpack_data_release(mt);
    }
    else {
      gkyl_dynvec_awrite(app->dts, fileNm);
    }
    app->stat.app_io_tm += gkyl_time_diff_now_sec(wtm);
    app->stat.n_diag_io += 1;
  }
  gkyl_dynvec_clear(app->dts);
}

void
gkyl_gyrokinetic_multib_app_save_dt(gkyl_gyrokinetic_multib_app* app, double tm, double dt)
{
  gkyl_dynvec_append(app->dts, tm, &dt);
}

void
gkyl_gyrokinetic_multib_app_release_geom(gkyl_gyrokinetic_multib_app* mbapp)
{
  if (mbapp->singleb_apps) {
    for (int i=0; i<mbapp->num_local_blocks; ++i)
      gkyl_gyrokinetic_app_release_geom(mbapp->singleb_apps[i]);
    gkyl_free(mbapp->singleb_apps);
  }  

  for (int bI=0; bI<mbapp->num_local_blocks; ++bI) {
    gkyl_multib_comm_conn_release(mbapp->mbcc_sync_conf->send[bI]);
    gkyl_multib_comm_conn_release(mbapp->mbcc_sync_conf->recv[bI]);
  }
  gkyl_free(mbapp->mbcc_sync_conf->send);
  gkyl_free(mbapp->mbcc_sync_conf->recv);
  gkyl_free(mbapp->mbcc_sync_conf);

  int num_blocks = gkyl_gk_block_geom_num_blocks(mbapp->gk_block_geom);

  for (int i=0; i<num_blocks; ++i)
    gkyl_comm_release(mbapp->block_comms[i]);
  gkyl_free(mbapp->block_comms);

  for (int i=0; i<num_blocks; ++i)
    gkyl_rect_decomp_release(mbapp->decomp[i]);
  gkyl_free(mbapp->decomp);

  gkyl_free(mbapp->local_blocks);    

  gkyl_rrobin_decomp_release(mbapp->round_robin);
  
  gkyl_gk_block_geom_release(mbapp->gk_block_geom);
  gkyl_block_topo_release(mbapp->block_topo);
  
  gkyl_comm_release(mbapp->comm);

  gkyl_msgpack_map_elem_release(mbapp->io_meta_basic_len, mbapp->io_meta_basic);

  gkyl_free(mbapp);
}

void
gkyl_gyrokinetic_multib_app_release(gkyl_gyrokinetic_multib_app* mbapp)
{
  for (int i=0; i<mbapp->num_neut_species; ++i) {
    for (int bI=0; bI<mbapp->num_local_blocks; ++bI) {
      gkyl_multib_comm_conn_release(mbapp->mbcc_sync_neut[i].send[bI]);
      gkyl_multib_comm_conn_release(mbapp->mbcc_sync_neut[i].recv[bI]);
    }
    gkyl_free(mbapp->mbcc_sync_neut[i].send);
    gkyl_free(mbapp->mbcc_sync_neut[i].recv);
  }
  gkyl_free(mbapp->mbcc_sync_neut);

  for (int i=0; i<mbapp->num_species; ++i) {
    for (int bI=0; bI<mbapp->num_local_blocks; ++bI) {
      gkyl_multib_comm_conn_release(mbapp->mbcc_sync_charged[i].send[bI]);
      gkyl_multib_comm_conn_release(mbapp->mbcc_sync_charged[i].recv[bI]);
    }
    gkyl_free(mbapp->mbcc_sync_charged[i].send);
    gkyl_free(mbapp->mbcc_sync_charged[i].recv);
  }
  gkyl_free(mbapp->mbcc_sync_charged);

  for (int bI=0; bI<mbapp->num_local_blocks; ++bI) {
    gkyl_multib_comm_conn_release(mbapp->mbcc_sync_conf->send[bI]);
    gkyl_multib_comm_conn_release(mbapp->mbcc_sync_conf->recv[bI]);
  }
  gkyl_free(mbapp->mbcc_sync_conf->send);
  gkyl_free(mbapp->mbcc_sync_conf->recv);
  gkyl_free(mbapp->mbcc_sync_conf);

  if (mbapp->singleb_apps) {
    for (int i=0; i<mbapp->num_local_blocks; ++i)
      gkyl_gyrokinetic_app_release(mbapp->singleb_apps[i]);
    gkyl_free(mbapp->singleb_apps);
  }  

  gk_multib_field_release(mbapp->field);

  int num_blocks = gkyl_gk_block_geom_num_blocks(mbapp->gk_block_geom);

  for (int i=0; i<num_blocks; ++i)
    gkyl_comm_release(mbapp->block_comms[i]);
  gkyl_free(mbapp->block_comms);

  for (int i=0; i<num_blocks; ++i)
    gkyl_rect_decomp_release(mbapp->decomp[i]);
  gkyl_free(mbapp->decomp);

  gkyl_free(mbapp->local_blocks);    

  gkyl_rrobin_decomp_release(mbapp->round_robin);
  
  gkyl_gk_block_geom_release(mbapp->gk_block_geom);
  gkyl_block_topo_release(mbapp->block_topo);
  
  gkyl_comm_release(mbapp->comm);

  gkyl_dynvec_release(mbapp->dts);

  gkyl_free(mbapp);
}
