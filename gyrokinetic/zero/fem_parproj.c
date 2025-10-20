#include <gkyl_fem_parproj.h>
#include <gkyl_fem_parproj_priv.h>

struct gkyl_fem_parproj*
gkyl_fem_parproj_new(const struct gkyl_range *solve_range,
  const struct gkyl_basis *basis, enum gkyl_fem_parproj_bc_type bctype,
  const struct gkyl_array *weight_left, const struct gkyl_array *weight_right, bool use_gpu)
{
  struct gkyl_fem_parproj *up = gkyl_malloc(sizeof(struct gkyl_fem_parproj));

  up->kernels = gkyl_malloc(sizeof(struct gkyl_fem_parproj_kernels));
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    up->kernels_cu = gkyl_cu_malloc(sizeof(struct gkyl_fem_parproj_kernels));
  else
    up->kernels_cu = up->kernels;
#else
  up->kernels_cu = up->kernels;
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

  // Range of parallel cells, as a sub-range of up->solve_range.
  struct gkyl_range par_range;
  int sublower[GKYL_MAX_CDIM], subupper[GKYL_MAX_CDIM];
  for (int d=0; d<up->ndim; d++) {
    sublower[d] = up->solve_range->lower[d];
    subupper[d] = up->solve_range->lower[d];
  }
  subupper[up->pardir] = up->solve_range->upper[up->pardir];
  gkyl_sub_range_init(&par_range, up->solve_range, sublower, subupper);
  up->parnum_cells = par_range.volume;

  // Range of perpendicular cells.
  struct gkyl_range perp_range;
  gkyl_range_shorten_from_above(&perp_range, up->solve_range, up->pardir, 1);

  // 1D range of parallel cells.
  int lower1d[] = {par_range.lower[up->pardir]}, upper1d[] = {par_range.upper[up->pardir]};
  gkyl_range_init(&up->par_range1d, 1, lower1d, upper1d);
  // 2D range of perpendicular cells.
  gkyl_range_init(&up->perp_range2d, up->ndim==3 ? 2 : 1, perp_range.lower, perp_range.upper);

  // Compute the number of local and global nodes.
  up->numnodes_local = up->num_basis;
  up->numnodes_global = gkyl_fem_parproj_global_num_nodes(basis, up->isperiodic, par_range.volume);

  up->brhs = gkyl_array_new(GKYL_DOUBLE, 1, up->numnodes_global*perp_range.volume); // Global right side vector.

  // Select local-to-global mapping kernel:
  fem_parproj_choose_local2global_kernel(basis, up->isperiodic, up->kernels->l2g);

  // Select weighted LHS kernel (not always used):
  fem_parproj_choose_lhs_kernel(basis, up->isdirichlet, has_weight_lhs, up->kernels->lhsker);

  // Select RHS source kernel:
  fem_parproj_choose_srcstencil_kernel(basis, up->isdirichlet, up->has_weight_rhs, up->kernels->srcker);

  // Select kernel that fetches the solution:
  up->kernels->solker = fem_parproj_choose_solstencil_kernel(basis);

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    fem_parproj_choose_kernels_cu(basis, up->has_weight_rhs, up->isperiodic, up->isdirichlet, up->kernels_cu);
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
      gkyl_range_init(&prob_range, up->perp_range2d.ndim, up->perp_range2d.lower, up->perp_range2d.upper);
    }
    else {
      nrhs = perp_range.volume;
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

  int idx1[GKYL_MAX_CDIM];
  gkyl_range_iter_init(&up->perp_iter2d, &prob_range);
  while (gkyl_range_iter_next(&up->perp_iter2d)) {
    long perpidx = gkyl_range_idx(&prob_range, up->perp_iter2d.idx);

    gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d);
    while (gkyl_range_iter_next(&up->par_iter1d)) {
      long paridx = gkyl_range_idx(&up->par_range1d, up->par_iter1d.idx);

      const double *wgt_p = NULL;
      if (has_weight_lhs) {
        for (size_t d=0; d<up->pardir; d++) idx1[d] = up->perp_iter2d.idx[d];
        idx1[up->pardir] = up->par_iter1d.idx[0];
        long linidx = gkyl_range_idx(up->solve_range, idx1);
        wgt_p = gkyl_array_cfetch(weight_left_ho, linidx);
      }

      int keri = up->par_iter1d.idx[0] == up->parnum_cells? 1 : 0;
      up->kernels->l2g[keri](up->parnum_cells, paridx, up->globalidx);

      // Apply the wgt*phi*basis stencil.
      keri = idx_to_inloup_ker(up->parnum_cells, up->par_iter1d.idx[0]);
      up->kernels->lhsker[keri](wgt_p, up->globalidx, tri[perpidx]);
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
gkyl_fem_parproj_set_rhs(struct gkyl_fem_parproj* up, const struct gkyl_array *rhsin, const struct gkyl_array *phibc)
{

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    assert(gkyl_array_is_cu_dev(rhsin));
    if (phibc)
      assert(gkyl_array_is_cu_dev(phibc));

    gkyl_fem_parproj_set_rhs_cu(up, rhsin, phibc);
    return;
  }
#endif

  gkyl_array_clear(up->brhs, 0.0);
  double *brhs_p = gkyl_array_fetch(up->brhs, 0);

  int idx1[GKYL_MAX_CDIM];
  gkyl_range_iter_init(&up->perp_iter2d, &up->perp_range2d);
  while (gkyl_range_iter_next(&up->perp_iter2d)) {
    long perpidx = gkyl_range_idx(&up->perp_range2d, up->perp_iter2d.idx);

    gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d);
    while (gkyl_range_iter_next(&up->par_iter1d)) {
      long paridx = gkyl_range_idx(&up->par_range1d, up->par_iter1d.idx);

      for (size_t d=0; d<up->pardir; d++) idx1[d] = up->perp_iter2d.idx[d];
      idx1[up->pardir] = up->par_iter1d.idx[0];
      long linidx = gkyl_range_idx(up->solve_range, idx1);

      const double *wgt_p = up->has_weight_rhs? gkyl_array_cfetch(up->weight_rhs, linidx) : NULL;
      const double *phibc_p = up->isdirichlet? gkyl_array_cfetch(phibc, linidx) : NULL;
      const double *rhsin_p = gkyl_array_cfetch(rhsin, linidx);

      long perpProbOff = perpidx*up->numnodes_global;

      int keri = up->par_iter1d.idx[0] == up->parnum_cells? 1 : 0;
      up->kernels->l2g[keri](up->parnum_cells, paridx, up->globalidx);

      keri = idx_to_inloup_ker(up->parnum_cells, up->par_iter1d.idx[0]);
      up->kernels->srcker[keri](wgt_p, rhsin_p, phibc_p, perpProbOff, up->globalidx, brhs_p);
    }
  }

  gkyl_superlu_brhs_from_array(up->prob, brhs_p);

}

void
gkyl_fem_parproj_solve(struct gkyl_fem_parproj* up, struct gkyl_array *phiout) {
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    assert(gkyl_array_is_cu_dev(phiout));
    gkyl_fem_parproj_solve_cu(up, phiout);
    return;
  }
#endif

  gkyl_superlu_solve(up->prob);

  int idx1[GKYL_MAX_CDIM];
  gkyl_range_iter_init(&up->perp_iter2d, &up->perp_range2d);
  while (gkyl_range_iter_next(&up->perp_iter2d)) {
    long perpidx = gkyl_range_idx(&up->perp_range2d, up->perp_iter2d.idx);

    gkyl_range_iter_init(&up->par_iter1d, &up->par_range1d);
    while (gkyl_range_iter_next(&up->par_iter1d)) {
      long paridx = gkyl_range_idx(&up->par_range1d, up->par_iter1d.idx);

      for (size_t d=0; d<up->pardir; d++) idx1[d] = up->perp_iter2d.idx[d];
      idx1[up->pardir] = up->par_iter1d.idx[0];
      long linidx = gkyl_range_idx(up->solve_range, idx1);

      double *phiout_p = gkyl_array_fetch(phiout, linidx);

      long perpProbOff = perpidx*up->numnodes_global;

      int keri = up->par_iter1d.idx[0] == up->parnum_cells? 1 : 0;
      up->kernels->l2g[keri](up->parnum_cells, paridx, up->globalidx);

      up->kernels->solker(gkyl_superlu_get_rhs_ptr(up->prob, 0), perpProbOff, up->globalidx, phiout_p);
    }
  }

}

void gkyl_fem_parproj_release(struct gkyl_fem_parproj *up)
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
