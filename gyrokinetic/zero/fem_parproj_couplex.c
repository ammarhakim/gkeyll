#include <gkyl_fem_parproj_couplex.h>
#include <gkyl_fem_parproj_couplex_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <assert.h>

struct gkyl_fem_parproj_couplex*
gkyl_fem_parproj_couplex_new(const struct gkyl_range *solve_range,
  const struct gkyl_basis *basis, enum gkyl_fem_parproj_bc_type bctype,
  const struct gkyl_array *weight_left, const struct gkyl_array *weight_right, bool use_gpu)
{
  struct gkyl_fem_parproj_couplex *up = gkyl_malloc(sizeof(struct gkyl_fem_parproj_couplex));

  up->kernels = gkyl_malloc(sizeof(struct gkyl_fem_parproj_couplex_kernels));
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    up->kernels_cu = gkyl_cu_malloc(sizeof(struct gkyl_fem_parproj_couplex_kernels));
#endif

  up->solve_range = solve_range;
  up->ndim = solve_range->ndim;
  up->num_basis  = basis->num_basis;
  up->basis_type = basis->b_type;
  up->poly_order = basis->poly_order;
  up->pardir = up->ndim-1; // Assume parallel direction is always the last.
  up->isperiodic = bctype == GKYL_FEM_PARPROJ_PERIODIC;
  up->isdirichlet = bctype == GKYL_FEM_PARPROJ_DIRICHLET;
  up->use_gpu = use_gpu;

  up->has_weight_rhs = false;
  if (weight_right) {
    up->has_weight_rhs = true;
    up->weight_rhs = gkyl_array_acquire(weight_right);
  }

  bool has_weight_lhs = false;
  struct gkyl_array *weight_left_ho;
  if (weight_left) {
    has_weight_lhs = true;
    weight_left_ho = use_gpu? gkyl_array_new(GKYL_DOUBLE, weight_left->ncomp, weight_left->size)
                            : gkyl_array_acquire(weight_left);
    gkyl_array_copy(weight_left_ho, weight_left);
  }

  up->globalidx = gkyl_malloc(sizeof(long[up->num_basis]));

  for (int d=0; d<up->ndim; d++)
    up->num_cells[d] = up->solve_range->upper[d] - up->solve_range->lower[d] + 1;

  // Range of cells spanned by FEM basis.
  struct gkyl_range fem_range_sub;
  int sublower[GKYL_MAX_CDIM], subupper[GKYL_MAX_CDIM];
  for (int d=0; d<up->ndim; d++) {
    sublower[d] = up->solve_range->lower[d];
    subupper[d] = up->solve_range->upper[d];
  }
  if (up->ndim == 3) {
    subupper[1] = up->solve_range->lower[1];
  }
  gkyl_sub_range_init(&fem_range_sub, up->solve_range, sublower, subupper);
  
  // Range of perpendicular cells (y)
  struct gkyl_range perp_range_sub;
  for (int d=0; d<up->ndim; d++) {
    sublower[d] = up->solve_range->lower[d];
    subupper[d] = up->solve_range->lower[d];
  }
  if (up->ndim == 3) {
    subupper[1] = up->solve_range->upper[1];
  }
  gkyl_sub_range_init(&perp_range_sub, up->solve_range, sublower, subupper);

  // 2D range of cells spanned by FEM basis.
  sublower[0] = up->solve_range->lower[0];
  subupper[0] = up->solve_range->upper[0];
  for (int d=up->pardir; d<up->pardir+1; d++) {
    sublower[1] = up->solve_range->lower[d];
    subupper[1] = up->solve_range->upper[d];
  }
  gkyl_range_init(&up->fem_range, GKYL_MIN2(2,up->ndim), sublower, subupper);

  // 1D range of perpendicular cells (y).
  sublower[0] = perp_range_sub.lower[0];
  subupper[0] = perp_range_sub.upper[0];
  if (up->ndim == 3) {
    sublower[0] = perp_range_sub.lower[1];
    subupper[0] = perp_range_sub.upper[1];
  }
  gkyl_range_init(&up->perp_range, 1, sublower, subupper);

  up->numnodes_local = up->num_basis;
  up->numnodes_global = gkyl_fem_parproj_couplex_global_num_nodes(basis, up->isperiodic, up->num_cells);
  
  up->brhs = gkyl_array_new(GKYL_DOUBLE, 1, up->numnodes_global*perp_range_sub.volume); // Global right side vector.

  fem_parproj_couplex_choose_kernels(basis, has_weight_lhs, up->has_weight_rhs, up->isperiodic, up->isdirichlet, up->kernels);

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    fem_parproj_couplex_choose_kernels_cu(basis, up->has_weight_rhs, up->isperiodic, up->isdirichlet, up->kernels_cu);
#endif

  // We support two cases:
  //  a) No weight, or weight is a single number so we can divide the RHS by it.
  //     Then we solve Ax=B where A is the discrete FEM projection operator,
  //     and B is a matrix with a column for each perpendicular cell.
  //  b) There's a spatially dependent weight. Then we solve A_i x_i=B_i
  //     where there's a different A_i for each perp cell and B_i is a single
  //     column matrix.
  struct gkyl_range prob_range;
  int nrhs;
  if (up->ndim == 1) {
    nrhs = 1;
    gkyl_range_init(&prob_range, 1, &((int){1}), &((int){1}));
  }
  else {
    if (has_weight_lhs) {
      nrhs = 1;
      gkyl_range_init(&prob_range, up->perp_range.ndim, up->perp_range.lower, up->perp_range.upper);
    }
    else {
      nrhs = perp_range_sub.volume;
      gkyl_range_init(&prob_range, 1, &((int){1}), &((int){1}));
    }
  }

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    up->prob_cu = gkyl_culinsolver_prob_new(prob_range.volume, up->numnodes_global, up->numnodes_global, nrhs);
  else
    up->prob = gkyl_superlu_prob_new(prob_range.volume, up->numnodes_global, up->numnodes_global, nrhs);
#else
  up->prob = gkyl_superlu_prob_new(prob_range.volume, up->numnodes_global, up->numnodes_global, nrhs);
#endif

  // Assign non-zero elements in A.
  struct gkyl_mat_triples **tri = gkyl_malloc(prob_range.volume*sizeof(struct gkyl_mat_triples *));
  for (size_t i=0; i<prob_range.volume; i++) {
    tri[i] = gkyl_mat_triples_new(up->numnodes_global, up->numnodes_global);
#ifdef GKYL_HAVE_CUDA
    if (up->use_gpu) gkyl_mat_triples_set_rowmaj_order(tri[i]);
#endif
  }

//  printf(" ---- ---- --- ---- -\n");
  int idx0[GKYL_MAX_CDIM], idx1[GKYL_MAX_CDIM];
  gkyl_range_iter_init(&up->perp_iter, &prob_range);
  while (gkyl_range_iter_next(&up->perp_iter)) {
    long linidx_perp = gkyl_range_idx(&prob_range, up->perp_iter.idx);

    gkyl_range_iter_init(&up->fem_iter, &up->fem_range);
    while (gkyl_range_iter_next(&up->fem_iter)) {
      long linidx_fem = gkyl_range_idx(&up->fem_range, up->fem_iter.idx);

      idx1[0] = up->fem_iter.idx[0];
      for (size_t d=1; d<up->pardir; d++) idx1[d] = up->perp_iter.idx[0];
      idx1[up->pardir] = up->fem_iter.idx[GKYL_MIN2(up->pardir,1)];

      for (size_t d=0; d<up->ndim; d++) idx0[d] = idx1[d]-up->solve_range->lower[d];

      const double *wgt_p = NULL;
      if (has_weight_lhs) {
        long linidx = gkyl_range_idx(up->solve_range, idx1);
        wgt_p = gkyl_array_cfetch(weight_left_ho, linidx);
      }

      int keri = idx_to_inup_ker(up->fem_range.ndim, up->fem_range.upper, up->fem_iter.idx);
      up->kernels->l2g[keri](up->num_cells, idx0, up->globalidx);
//      if (up->ndim == 2)
//        printf("idx1=%d,%d | keri l2g=%d | globalidx=%ld,%ld,%ld,%ld",idx1[0],idx1[1],keri,up->globalidx[0],up->globalidx[1],up->globalidx[2],up->globalidx[3]);
//      else if (up->ndim == 3)
//        printf("idx1=%d,%d,%d | keri l2g=%d | globalidx=%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld",idx1[0],idx1[1],idx1[2],keri,
//          up->globalidx[0],up->globalidx[1],up->globalidx[2],up->globalidx[3],
//          up->globalidx[4],up->globalidx[5],up->globalidx[6],up->globalidx[7]
//          );

      // Apply the wgt*phi*basis stencil.
      keri = idx_to_inloup_ker(up->fem_range.ndim, up->fem_range.lower, up->fem_range.upper, up->fem_iter.idx);
//      printf(" keri lhs=%d | linidx_perp=%ld\n",keri,linidx_perp);
      up->kernels->lhsker[keri](wgt_p, up->globalidx, tri[linidx_perp]);
    }
  }

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    gkyl_culinsolver_amat_from_triples(up->prob_cu, tri);
  else
    gkyl_superlu_amat_from_triples(up->prob, tri);
#else
  gkyl_superlu_amat_from_triples(up->prob, tri);
#endif

  for (size_t i=0; i<prob_range.volume; i++)
    gkyl_mat_triples_release(tri[i]);
  gkyl_free(tri);

  if (weight_left)
    gkyl_array_release(weight_left_ho);

  return up;
}

void
gkyl_fem_parproj_couplex_set_rhs(struct gkyl_fem_parproj_couplex* up,
  const struct gkyl_array *rhsin, const struct gkyl_array *phibc)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    assert(gkyl_array_is_cu_dev(rhsin));
    if (phibc)
      assert(gkyl_array_is_cu_dev(phibc));

    gkyl_fem_parproj_couplex_set_rhs_cu(up, rhsin, phibc);
    return;
  }
#endif

  gkyl_array_clear(up->brhs, 0.0);
  double *brhs_p = gkyl_array_fetch(up->brhs, 0);

//  printf(" ---- ---- --- ---- -\n");
  int idx0[GKYL_MAX_CDIM], idx1[GKYL_MAX_CDIM];
  gkyl_range_iter_init(&up->perp_iter, &up->perp_range);
  while (gkyl_range_iter_next(&up->perp_iter)) {
    long linidx_perp = gkyl_range_idx(&up->perp_range, up->perp_iter.idx);

    gkyl_range_iter_init(&up->fem_iter, &up->fem_range);
    while (gkyl_range_iter_next(&up->fem_iter)) {
      long linidx_fem = gkyl_range_idx(&up->fem_range, up->fem_iter.idx);

      idx1[0] = up->fem_iter.idx[0];
      for (size_t d=1; d<up->pardir; d++) idx1[d] = up->perp_iter.idx[0];
      idx1[up->pardir] = up->fem_iter.idx[GKYL_MIN2(up->pardir,1)];

      for (size_t d=0; d<up->ndim; d++) idx0[d] = idx1[d]-up->solve_range->lower[d];

      int keri = idx_to_inup_ker(up->fem_range.ndim, up->fem_range.upper, up->fem_iter.idx);
      up->kernels->l2g[keri](up->num_cells, idx0, up->globalidx);
//      if (up->ndim == 2)
//        printf("idx1=%d,%d | keri l2g=%d | globalidx=%ld,%ld,%ld,%ld",idx1[0],idx1[1],keri,up->globalidx[0],up->globalidx[1],up->globalidx[2],up->globalidx[3]);
//      else if (up->ndim == 3)
//        printf("idx1=%d,%d,%d | keri l2g=%d | globalidx=%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld",idx1[0],idx1[1],idx1[2],keri,
//          up->globalidx[0],up->globalidx[1],up->globalidx[2],up->globalidx[3],
//          up->globalidx[4],up->globalidx[5],up->globalidx[6],up->globalidx[7]
//          );

      long linidx = gkyl_range_idx(up->solve_range, idx1);

      const double *wgt_p = up->has_weight_rhs? gkyl_array_cfetch(up->weight_rhs, linidx) : NULL;
      const double *phibc_p = up->isdirichlet? gkyl_array_cfetch(phibc, linidx) : NULL;
      const double *rhsin_p = gkyl_array_cfetch(rhsin, linidx);

      long perpProbOff = linidx_perp*up->numnodes_global;

      keri = idx_to_inloup_ker(up->fem_range.ndim, up->fem_range.lower, up->fem_range.upper, up->fem_iter.idx);
//      printf(" keri src=%d | linidx_perp=%ld\n",keri,linidx_perp);
      up->kernels->srcker[keri](wgt_p, rhsin_p, phibc_p, perpProbOff, up->globalidx, brhs_p);
    }
  }

  gkyl_superlu_brhs_from_array(up->prob, brhs_p);
}

void
gkyl_fem_parproj_couplex_solve(struct gkyl_fem_parproj_couplex* up, struct gkyl_array *phiout) {
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    assert(gkyl_array_is_cu_dev(phiout));
    gkyl_fem_parproj_couplex_solve_cu(up, phiout);
    return;
  }
#endif

  gkyl_superlu_solve(up->prob);

//  printf(" ---- ---- --- ---- -\n");
  int idx0[GKYL_MAX_CDIM], idx1[GKYL_MAX_CDIM];
  gkyl_range_iter_init(&up->perp_iter, &up->perp_range);
  while (gkyl_range_iter_next(&up->perp_iter)) {
    long linidx_perp = gkyl_range_idx(&up->perp_range, up->perp_iter.idx);

    gkyl_range_iter_init(&up->fem_iter, &up->fem_range);
    while (gkyl_range_iter_next(&up->fem_iter)) {
      long linidx_fem = gkyl_range_idx(&up->fem_range, up->fem_iter.idx);

      idx1[0] = up->fem_iter.idx[0];
      for (size_t d=1; d<up->pardir; d++) idx1[d] = up->perp_iter.idx[0];
      idx1[up->pardir] = up->fem_iter.idx[GKYL_MIN2(up->pardir,1)];

      for (size_t d=0; d<up->ndim; d++) idx0[d] = idx1[d]-up->solve_range->lower[d];

      int keri = idx_to_inup_ker(up->fem_range.ndim, up->fem_range.upper, up->fem_iter.idx);
      up->kernels->l2g[keri](up->num_cells, idx0, up->globalidx);

//      printf("idx1=%d,%d | keri l2g=%d | globalidx=%ld,%ld,%ld,%ld",idx1[0],idx1[1],keri,up->globalidx[0],up->globalidx[1],up->globalidx[2],up->globalidx[3]);
      long linidx = gkyl_range_idx(up->solve_range, idx1);

      double *phiout_p = gkyl_array_fetch(phiout, linidx);

      long perpProbOff = linidx_perp*up->numnodes_global;

//      printf(" keri sol=%d | linidx_perp=%ld\n",keri,linidx_perp);
      up->kernels->solker(gkyl_superlu_get_rhs_ptr(up->prob, 0), perpProbOff, up->globalidx, phiout_p);
    }
  }
}

void gkyl_fem_parproj_couplex_release(struct gkyl_fem_parproj_couplex *up)
{
  if (up->has_weight_rhs) {
    gkyl_array_release(up->weight_rhs);
  }
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_cu_free(up->kernels_cu);
    gkyl_culinsolver_prob_release(up->prob_cu);
  } else {
    gkyl_superlu_prob_release(up->prob);
  }
#else
  gkyl_superlu_prob_release(up->prob);
#endif
  gkyl_array_release(up->brhs);
  gkyl_free(up->globalidx);
  gkyl_free(up->kernels);
  gkyl_free(up);
}
