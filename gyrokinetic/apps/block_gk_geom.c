#include <gkyl_gk_block_geom.h>
#include <gkyl_alloc.h>
#include <gkyl_tok_geo.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Geometry info for all blocks in simulation
struct gkyl_gk_block_geom {
  int ndim; // dimension
  int num_blocks; // total number of blocks
  struct gkyl_gk_block_geom_info *blocks; // info for each block
  struct gkyl_block_topo *btopo; // topology of blocks
  
  struct gkyl_ref_count ref_count;
};

static void
gk_block_geom_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_gk_block_geom *bgeom = container_of(ref, struct gkyl_gk_block_geom, ref_count);
  gkyl_free(bgeom->blocks);
  gkyl_block_topo_release(bgeom->btopo);
  gkyl_free(bgeom);
}

struct gkyl_gk_block_geom*
gkyl_gk_block_geom_new(int ndim, int nblocks)
{
  struct gkyl_gk_block_geom *bgeom = gkyl_malloc(sizeof(struct gkyl_gk_block_geom));
  bgeom->ndim = ndim;
  bgeom->num_blocks = nblocks;
  bgeom->blocks = gkyl_calloc(sizeof(struct gkyl_gk_block_geom_info), nblocks);

  bgeom->btopo = gkyl_block_topo_new(ndim, nblocks);

  bgeom->ref_count = gkyl_ref_count_init(gk_block_geom_free);

  return bgeom;
}

int
gkyl_gk_block_geom_ndim(const struct gkyl_gk_block_geom *bgeom)
{
  return bgeom->ndim;
}

int
gkyl_gk_block_geom_num_blocks(const struct gkyl_gk_block_geom *bgeom)
{
  return bgeom->num_blocks;
}

void
gkyl_gk_block_geom_set_block(struct gkyl_gk_block_geom *bgeom, int bidx,
  const struct gkyl_gk_block_geom_info *info)
{
  memcpy(&bgeom->blocks[bidx], info, sizeof(struct gkyl_gk_block_geom_info));
  
  for (int d=0; d<bgeom->ndim; ++d)
    bgeom->blocks[bidx].cuts[d] = info->cuts[d] > 0 ? info->cuts[d] : 1;
  
  // set topology information
  for (int i=0; i<bgeom->ndim; ++i)
    for (int e=0; e<2; ++e)
      bgeom->btopo->conn[bidx].connections[i][e] = info->connections[i][e];  
}

void
gkyl_gk_block_geom_reset_block_extents(struct gkyl_gk_block_geom *bgeom, int bidx, double *lower, double *upper)
{
  struct gkyl_gk_block_geom_info *bgi = &bgeom->blocks[bidx];
  for (int i = 0; i < bgeom->ndim; ++i) {
    bgi->lower[i] = lower[i];
    bgi->upper[i] = upper[i];
  }
}

void
gkyl_gk_block_geom_apply_xpt_seam_selection(struct gkyl_gk_block_geom *bgeom,
  int bidx, double coefficient, double bound)
{
  struct gkyl_gk_block_geom_info *bgi = &bgeom->blocks[bidx];
  struct gkyl_tok_geo_grid_inp *inp = &bgi->geometry.tok_grid_info;
  inp->relaxed_xpt_seam = true;
  inp->relaxed_xpt_seam_sweep = true;
  inp->relaxed_xpt_seam_delta_s_coeff = coefficient;
  inp->relaxed_xpt_seam_delta_s_bound = bound;
}

const struct gkyl_gk_block_geom_info*
gkyl_gk_block_geom_get_block(const struct gkyl_gk_block_geom *bgeom, int bidx)
{
  return &bgeom->blocks[bidx];
}

// Report interfaces whose two blocks disagree about taking the extended
// construction.
//
// The extended path reparameterizes a block's separatrix row -- the theta
// ladder, and the |grad psi| poloidal measure -- while the legacy path leaves
// it on plain arc length.  Two blocks that SHARE that row and disagree about
// taking it therefore trace the same curve and still place different nodes
// along it.  Measured on asdex/tcv, which declare straight_xpt_ray on 3 of 6
// blocks: the two sides agree on arc length to 7e-05 relative and on both
// endpoints to 1e-14 m, and differ only in the interior, by 35-66 mm -- which
// is 0.6-0.8 cells at x1 and DOUBLES at every refinement, because it is a fixed
// offset rather than a convergence error.
//
// Every interface that fails is such a mixed pair, and every uniform pair
// passes; step (12 of 12 declared) and the NSTX-U half-domain set (all blocks
// declared) have no mixed pair and no seam failure.  So this is the invariant
// worth stating, and it is a property of the DECLARATION -- checkable here,
// on any device, at any resolution, without constructing anything.
//
// Note the participation predicate deliberately comes from
// gkyl_tok_geo_uses_extended_construction() rather than being restated: it also
// depends on half_domain and on an environment override, and a second copy
// would drift.
//
// This REPORTS by default and does not change the verdict, because existing
// multiblock cases carry mixed interfaces today and would begin failing at
// setup.  Set GKYL_TOK_STRICT_SEAM_PARTICIPATION=1 to make it an error.
static int
gk_block_geom_check_seam_participation(const struct gkyl_gk_block_geom *bgeom)
{
  const char *strict_env = getenv("GKYL_TOK_STRICT_SEAM_PARTICIPATION");
  bool strict = strict_env && strict_env[0] != '\0' && strict_env[0] != '0';
  const char *diag_env = getenv("GKYL_TOK_SEAM_PARTICIPATION_DIAG");
  bool diag = diag_env && diag_env[0] != '\0' && diag_env[0] != '0';
  int nmixed = 0, nchecked = 0;

  // "No mixed interface" and "the check never ran" look identical from the
  // outside, and that cost a debugging cycle when this was first wired up.
  // The diagnostic distinguishes them by naming every block and every edge it
  // actually examined.
  if (diag) {
    fprintf(stderr, "TOK_SEAM_PARTICIPATION_DIAG num_blocks=%d ndim=%d\n",
      bgeom->num_blocks, bgeom->ndim);
    for (int i=0; i<bgeom->num_blocks; ++i) {
      const struct gkyl_gk_block_geom_info *bi = &bgeom->blocks[i];
      fprintf(stderr,
        "TOK_SEAM_PARTICIPATION_DIAG block=%d geometry_id=%d ftype=%d "
        "straight_xpt_ray=%d half_domain=%d extended=%d\n",
        i, (int) bi->geometry.geometry_id,
        (int) bi->geometry.tok_grid_info.ftype,
        (int) bi->geometry.tok_grid_info.straight_xpt_ray,
        (int) bi->geometry.tok_grid_info.half_domain,
        (int) gkyl_tok_geo_uses_extended_construction(&bi->geometry.tok_grid_info));
      for (int d=0; d<bgeom->ndim; ++d)
        for (int e=0; e<2; ++e) {
          const struct gkyl_target_edge *te = &bgeom->btopo->conn[i].connections[d][e];
          fprintf(stderr,
            "TOK_SEAM_PARTICIPATION_DIAG   block=%d dir=%d edge=%s -> bid=%d edgekind=%d\n",
            i, d, e ? "upper" : "lower", te->bid, (int) te->edge);
        }
    }
  }

  for (int i=0; i<bgeom->num_blocks; ++i) {
    const struct gkyl_gk_block_geom_info *bi = &bgeom->blocks[i];
    if (bi->geometry.geometry_id != GKYL_GEOMETRY_TOKAMAK)
      continue;

    for (int d=0; d<bgeom->ndim; ++d) {
      for (int e=0; e<2; ++e) {
        const struct gkyl_target_edge *te = &bgeom->btopo->conn[i].connections[d][e];
        if (te->edge == GKYL_PHYSICAL)
          continue;
        int j = te->bid;
        // Visit each interface once, and skip a block joined to itself: a
        // self-periodic block cannot disagree with itself.
        if (j <= i)
          continue;
        if (j < 0 || j >= bgeom->num_blocks)
          continue; // malformed; the topology check reports this itself

        const struct gkyl_gk_block_geom_info *bj = &bgeom->blocks[j];
        if (bj->geometry.geometry_id != GKYL_GEOMETRY_TOKAMAK)
          continue;

        nchecked += 1;
        bool ext_i = gkyl_tok_geo_uses_extended_construction(&bi->geometry.tok_grid_info);
        bool ext_j = gkyl_tok_geo_uses_extended_construction(&bj->geometry.tok_grid_info);
        if (ext_i == ext_j)
          continue;

        if (nmixed == 0)
          fprintf(stderr,
            "TOK_SEAM_PARTICIPATION mixed extended construction across a declared "
            "interface: the two blocks reparameterize their shared row "
            "differently, so its interior nodes will not line up and the gap "
            "doubles at every refinement. Make participation uniform across the "
            "seam (declare straight_xpt_ray on both, or neither).\n");

        fprintf(stderr,
          "TOK_SEAM_PARTICIPATION block=%d ftype=%d extended=%d <-> block=%d "
          "ftype=%d extended=%d dir=%d edge=%s\n",
          i, (int) bi->geometry.tok_grid_info.ftype, (int) ext_i,
          j, (int) bj->geometry.tok_grid_info.ftype, (int) ext_j,
          d, e ? "upper" : "lower");
        nmixed += 1;
      }
    }
  }

  if (diag)
    fprintf(stderr,
      "TOK_SEAM_PARTICIPATION_DIAG interfaces_examined=%d mixed=%d\n",
      nchecked, nmixed);

  if (nmixed > 0 && strict) {
    fprintf(stderr,
      "TOK_SEAM_PARTICIPATION %d mixed interface(s); failing because "
      "GKYL_TOK_STRICT_SEAM_PARTICIPATION is set\n", nmixed);
    return 0;
  }
  return 1;
}

int
gkyl_gk_block_geom_check_consistency(const struct gkyl_gk_block_geom *bgeom)
{
  // MORE TESTS ARE NEEDED HERE
  //
  // Run BOTH checks and combine, rather than short-circuiting on the topology.
  // They answer independent questions -- "are the connections well formed" and
  // "do connected blocks agree about how they parameterize the row they share"
  // -- and letting the first gate the second means a declaration error stays
  // invisible in exactly the cases that have a topology problem too.
  int topo_ok = gkyl_block_topo_check_consistency(bgeom->btopo);
  int seam_ok = gk_block_geom_check_seam_participation(bgeom);
  return topo_ok && seam_ok;
}

struct gkyl_gk_block_geom *
gkyl_gk_block_geom_acquire(const struct gkyl_gk_block_geom* bgeom)
{
  gkyl_ref_count_inc(&bgeom->ref_count);
  return (struct gkyl_gk_block_geom*) bgeom;
}

struct gkyl_block_topo*
gkyl_gk_block_geom_topo(const struct gkyl_gk_block_geom *bgeom)
{
  return gkyl_block_topo_acquire(bgeom->btopo);
}

void
gkyl_gk_block_geom_release(struct gkyl_gk_block_geom* bgeom)
{
  gkyl_ref_count_dec(&bgeom->ref_count);
}
