#include <gkyl_alloc.h>
#include <gkyl_mat_triples.h>
#include <gkyl_superlu_ops.h>

#include <stdbool.h>

struct gkyl_superlu_prob {
  SuperMatrix **A, **B; // matrices in A_j x_j = B_j problems.
  SuperMatrix **L, **U; // L and U factors in LU decomposition.
  double *rhs; // right-hand side entries. 
  int *perm_c; // column permutation vector (re-used for each problem).
  int **perm_r; // row permutations from partial pivoting.
  int mrow, ncol; // A is a mrow x ncol matrix.
  int nnz; // number of non-zero entries in A.
  int nprob; // number of problems to solve (j=0,1,...,nprob-1).
  int nrhs; // number of columns of the RHS (B is an mrow x nrhs matrix).

  double **nzvals; // non-zero matrix elements.
  int **rowinds;   // row index of entries in nzval.
  int **colptrs;   // 1st entry of each column as index in nzval.

  int info, permc_spec;
  superlu_options_t options;
  SuperLUStat_t stat;
  trans_t trans;
  GlobalLU_t Glu;
  bool assigned_rhs;
};

gkyl_superlu_prob*
gkyl_superlu_prob_new(int nprob, int mrow, int ncol, int nrhs)
{
  assert((nprob==1) || (nrhs==1));

  struct gkyl_superlu_prob *prob = gkyl_malloc(sizeof(*prob));

  prob->nprob = nprob;
  prob->mrow = mrow;
  prob->ncol = ncol;
  prob->nrhs = nrhs;

  prob->A = gkyl_malloc(prob->nprob*sizeof(SuperMatrix *));
  prob->B = gkyl_malloc(prob->nprob*sizeof(SuperMatrix *));
  prob->L = gkyl_malloc(prob->nprob*sizeof(SuperMatrix *));
  prob->U = gkyl_malloc(prob->nprob*sizeof(SuperMatrix *));
  for (size_t k=0; k<prob->nprob; k++) {
    prob->A[k] = gkyl_malloc(sizeof(SuperMatrix));
    prob->B[k] = gkyl_malloc(sizeof(SuperMatrix));
    prob->L[k] = gkyl_malloc(sizeof(SuperMatrix));
    prob->U[k] = gkyl_malloc(sizeof(SuperMatrix));
  }

  prob->rhs = doubleMalloc(mrow*GKYL_MAX2(nprob,nrhs));
  prob->perm_c = intMalloc(ncol);
  prob->perm_r = gkyl_malloc(prob->nprob*sizeof(int *));
  for (size_t k=0; k<prob->nprob; k++)
    prob->perm_r[k] = intMalloc(mrow);

  // Set the default input options.
  set_default_options(&prob->options);
  prob->options.ColPerm = NATURAL;

  // Initialize the statistics variables.
  StatInit(&prob->stat);

  prob->trans = NOTRANS;

  prob->assigned_rhs = false;

  return prob;
}

void
gkyl_superlu_amat_from_triples(struct gkyl_superlu_prob *prob, struct gkyl_mat_triples **tri)
{
  prob->nnz = gkyl_mat_triples_size(tri[0]);
  for (size_t k=0; k<prob->nprob; k++) {
    assert(gkyl_mat_triples_size(tri[k]) == prob->nnz);  // No. of nonzeros must be the same for every problem.
    assert(gkyl_mat_triples_is_colmaj(tri[k]));  // triples must be in colmaj order for superlu.
  }

  // Allocate some memory needed in superlu. NOTE: this memory is
  // deleted when Destroy_CompCol_Matrix is called, and so we do not
  // need to do it ourselves.
  prob->nzvals = gkyl_malloc(prob->nprob*sizeof(double *));
  prob->rowinds = gkyl_malloc(prob->nprob*sizeof(int *));
  prob->colptrs = gkyl_malloc(prob->nprob*sizeof(int *));
  for (size_t k=0; k<prob->nprob; k++) {
    prob->nzvals[k] = doubleMalloc(prob->nnz); // non-zero matrix elements.
    prob->rowinds[k] = intMalloc(prob->nnz); // row index of entries in nzval.
    prob->colptrs[k] = intMalloc(prob->ncol+1); // 1st entry of each column as index in nzval.
  }

  bool *colptr_assigned = gkyl_malloc(prob->ncol*sizeof(bool));
  // Sorted (column-major order) keys (linear indices to flattened matrix).
  for (size_t k=0; k<prob->nprob; k++) {

    double *nzval = prob->nzvals[k];
    int *rowind = prob->rowinds[k];
    int *colptr = prob->colptrs[k];

    for (size_t i=0; i<prob->ncol; i++)
      colptr_assigned[i] = false;

    gkyl_mat_triples_iter *iter = gkyl_mat_triples_iter_new(tri[k]);
    for (size_t i=0; i<prob->nnz; ++i) {
      gkyl_mat_triples_iter_next(iter); // bump iterator.
      struct gkyl_mtriple mt = gkyl_mat_triples_iter_at(iter);
      size_t idx[2] = { mt.row, mt.col };
      
      nzval[i] = mt.val;
      rowind[i] = idx[0];
      if (!colptr_assigned[idx[1]]) {
        colptr[idx[1]] = i;
        colptr_assigned[idx[1]] = true;
      }
    }
    colptr[prob->ncol] = prob->nnz;
    gkyl_mat_triples_iter_release(iter);

    // Create matrix A. See SuperLU manual for definitions.
    dCreate_CompCol_Matrix(prob->A[k], prob->mrow, prob->ncol, prob->nnz,
      nzval, rowind, colptr, SLU_NC, SLU_D, SLU_GE);
  }
  
  gkyl_free(colptr_assigned);

  prob->options.Fact = DOFACT; // Haven't computed LU decomp yet.
}

void
gkyl_superlu_print_amat(struct gkyl_superlu_prob *prob)
{
  char strA[5];
  for (size_t k=0; k<prob->nprob; k++) {
    snprintf(strA, 5, "A%zu", k); // puts string into buffer
    dPrint_CompCol_Matrix(strA, prob->A[k]);
  }
}

void
gkyl_superlu_ludecomp(struct gkyl_superlu_prob *prob)
{
  /*
  *   Get column permutation vector perm_c[], according to permc_spec:
  * = 0: natural ordering
  * = 1: minimum degree on structure of A’*A
  * = 2: minimum degree on structure of A’+A
  * = 3: approximate minimum degree for unsymmetric matrices
  */
  int permc_spec = 0; 
  get_perm_c(permc_spec, prob->A[0], prob->perm_c);

  int *etree; // Column elimination tree.
  if ( !(etree = intMalloc(prob->ncol)) ) ABORT("superlu_ops: Malloc fails for etree[].");
  SuperMatrix AC; // permutation matrix time A.
  sp_preorder(&prob->options, prob->A[0], prob->perm_c, etree, &AC);

  int panel_size = sp_ienv(1);
  int relax = sp_ienv(2);
  dgstrf(&prob->options, &AC, relax, panel_size, etree, NULL, 0, prob->perm_c,
    prob->perm_r[0], prob->L[0], prob->U[0], &prob->Glu, &prob->stat, &prob->info);

  prob->options.Fact = prob->nprob==1? FACTORED : SamePattern; // LU decomp done.

  for (size_t k=1; k<prob->nprob; k++) {
    dgstrf(&prob->options, &AC, relax, panel_size, etree, NULL, 0, prob->perm_c,
      prob->perm_r[k], prob->L[k], prob->U[k], &prob->Glu, &prob->stat, &prob->info);
  }

  SUPERLU_FREE(etree);
  Destroy_CompCol_Permuted(&AC);
}

void
gkyl_superlu_brhs_from_triples(struct gkyl_superlu_prob *prob, struct gkyl_mat_triples *tri)
{
  long nnz_rhs = gkyl_mat_triples_size(tri);  // number of non-zero entries in RHS matrix B
  
  // sorted (column-major order) keys (linear indices to flattened matrix)
  gkyl_mat_triples_iter *iter = gkyl_mat_triples_iter_new(tri);
  for (size_t i=0; i<nnz_rhs; i++) {
    gkyl_mat_triples_iter_next(iter); // bump iterator
    struct gkyl_mtriple mt = gkyl_mat_triples_iter_at(iter);    
    prob->rhs[i] = mt.val;
  }
  gkyl_mat_triples_iter_release(iter);
  
  // Create RHS matrix B. See SuperLU manual for definitions.
  for (size_t k=0; k<prob->nprob; k++)
    dCreate_Dense_Matrix(prob->B[k], prob->mrow, prob->nrhs, &prob->rhs[k*prob->mrow], prob->mrow,
      SLU_DN, SLU_D, SLU_GE);

  prob->assigned_rhs = true;
}

void
gkyl_superlu_brhs_from_array(struct gkyl_superlu_prob *prob, const double *bin)
{
  for (size_t i=0; i<prob->mrow*GKYL_MAX2(prob->nprob,prob->nrhs); i++)
    prob->rhs[i] = bin[i];
  
  // Create RHS matrix B. See SuperLU manual for definitions.
  for (size_t k=0; k<prob->nprob; k++)
    dCreate_Dense_Matrix(prob->B[k], prob->mrow, prob->nrhs, &prob->rhs[k*prob->mrow], prob->mrow,
      SLU_DN, SLU_D, SLU_GE);

  prob->assigned_rhs = true;
}

void
gkyl_superlu_solve(struct gkyl_superlu_prob *prob)
{
  if (prob->options.Fact==FACTORED) {
    for (size_t k=0; k<prob->nprob; k++)
      dgstrs(prob->trans, prob->L[k], prob->U[k], prob->perm_c, prob->perm_r[k], prob->B[k], &prob->stat, &prob->info);
  } else {
    dgssv(&prob->options, prob->A[0], prob->perm_c, prob->perm_r[0], prob->L[0], prob->U[0],
      prob->B[0], &prob->stat, &prob->info);

//  MF 2023/05/25: My intention was to re-use the column permutation vector perm_c by
//  using SamePattern. But for some reason it errors out.
//    prob->options.Fact = prob->nprob==1? FACTORED : SamePattern; // LU decomp done.

    for (size_t k=1; k<prob->nprob; k++)
      dgssv(&prob->options, prob->A[k], prob->perm_c, prob->perm_r[k], prob->L[k], prob->U[k],
        prob->B[k], &prob->stat, &prob->info);

    prob->options.Fact = FACTORED; // LU decomp done.
  }

  for (size_t k=0; k<prob->nprob; k++) {
    if (prob->assigned_rhs)
      Destroy_SuperMatrix_Store(prob->B[k]);
  }
}

double
gkyl_superlu_get_rhs_ij(struct gkyl_superlu_prob *prob, long ielement, long jprob)
{
  return prob->rhs[jprob*prob->mrow+ielement];
}


double
gkyl_superlu_get_rhs_lin(struct gkyl_superlu_prob *prob, long loc)
{
  return prob->rhs[loc];
}

double*
gkyl_superlu_get_rhs_ptr(struct gkyl_superlu_prob *prob, long loc)
{
  return &prob->rhs[loc];
}

void
gkyl_superlu_prob_release(struct gkyl_superlu_prob *prob)
{
  SUPERLU_FREE (prob->rhs);
  SUPERLU_FREE (prob->perm_c);
  
  for (size_t k=0; k<prob->nprob; k++) {
    SUPERLU_FREE (prob->perm_r[k]);
    Destroy_CompCol_Matrix(prob->A[k]);
    gkyl_free(prob->A[k]);

    if (prob->options.Fact==FACTORED) {
      Destroy_SuperNode_Matrix(prob->L[k]);
      Destroy_CompCol_Matrix(prob->U[k]);
    }

    gkyl_free(prob->B[k]);
    gkyl_free(prob->L[k]);
    gkyl_free(prob->U[k]);
  }
  StatFree(&prob->stat);
  gkyl_free(prob->A);
  gkyl_free(prob->B);
  gkyl_free(prob->L);
  gkyl_free(prob->U);
  gkyl_free(prob->nzvals);
  gkyl_free(prob->rowinds);
  gkyl_free(prob->colptrs);
  gkyl_free(prob->perm_r);
  gkyl_free(prob);
}
