#include <string.h>
#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_vlasov_triad_geom.h>
#include <gkyl_vlasov_triad_geom_priv.h>

// This function needs to:
// 1. Project the basis at Gauss-Legendre points in a cell
// 2. Compute h_ij
// 3. Compute sqrt(det(h_ij))
// 4. Compute h_ij_inv
// 5. Compute conf_poisson tensor at nodal configuration space points
// 6. Convert nodal to modal representations of h_ij, h_ij_inv, sqrt(det(h_ij))
// (LATER FOR GR ONLY) Use this information to contruct the continious hamil

void
gkyl_vlasov_triad_geom_new(const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange, const struct gkyl_basis cbasis, 
  const struct gkyl_rect_grid *pgrid, const struct gkyl_range *prange, const struct gkyl_basis pbasis, 
  struct gkyl_vlasov_triad_geom_inp inp_basis_vectors, struct gkyl_array *cov_tangent_basis,
  struct gkyl_array *triad_basis, struct gkyl_array *h_ij,
  struct gkyl_array *h_ij_inv, struct gkyl_array *det_h, struct gkyl_array *conf_poisson_tensor)
{

  // INCOMPLETE
  // 1. Need to generalize to multiple dimensions (cdim dimension)
    // See: ctest_dg_basis_ops.c (test_cubic_evalf_2d) for details
  // 2. Need to save the evalf results nodally
  // 3. Need to construct remaining physical quantites with some nodal to modal conv.
    // >? perhaps the way I verified the hamiltonian is build nodally is useful?
  // 4. Should remove iteration over cdim, it should go over the entire conf grid.

  int cdim = cgrid->ndim;
  int pdim = pgrid->ndim;
  int vdim = pdim - cdim;
  struct gkyl_array *cov_tangent_basis_nodal;
  struct gkyl_array *triad_basis_nodal;
  struct gkyl_array *triad_basis_gradient_nodal;
  struct gkyl_array *h_ij_nodal;
  struct gkyl_array *h_ij_inv_nodal;
  struct gkyl_array *det_h_nodal;
  struct gkyl_array *conf_poisson_tensor_nodal;

  // Configuration space nodal grid construction from cgrid
  const double *clower = cgrid->lower, *cupper = cgrid->upper;
  const double *plower = pgrid->lower, *pupper = pgrid->upper;
  const int *cells = cgrid->cells;

  // Nodal grid used for constructing physical coordinates in configuration space
  double nc_lower[GKYL_MAX_CDIM];
  double nc_upper[GKYL_MAX_CDIM];
  int nc_cells[GKYL_MAX_CDIM];
  int nghost[GKYL_MAX_CDIM];
  for (int i=0; i<cdim; ++i) {
    nc_lower[i] = clower[i] - 0.5*cgrid->dx[i];
    nc_upper[i] = cupper[i] + 0.5*cgrid->dx[i];
    nc_cells[i] = cells[i] + 1; 
    nghost[i] = 1; // One ghost cell in conf space
  }

  struct gkyl_rect_grid nc_grid;
  struct gkyl_range nc_local, nc_local_ext;
  gkyl_rect_grid_init(&nc_grid, cdim, nc_lower, nc_upper, nc_cells);
  gkyl_create_grid_ranges(&nc_grid, nghost, &nc_local_ext, &nc_local);

  cov_tangent_basis_nodal = gkyl_array_new(GKYL_DOUBLE, vdim*vdim, nc_local.volume);
  triad_basis_nodal = gkyl_array_new(GKYL_DOUBLE, vdim*vdim, nc_local.volume);
  h_ij_nodal = gkyl_array_new(GKYL_DOUBLE, vdim*(vdim+1)/2, nc_local.volume);
  h_ij_inv_nodal = gkyl_array_new(GKYL_DOUBLE, vdim*(vdim+1)/2, nc_local.volume);
  det_h_nodal = gkyl_array_new(GKYL_DOUBLE, 1, nc_local.volume);
  conf_poisson_tensor_nodal = gkyl_array_new(GKYL_DOUBLE, vdim*(vdim+cdim)*((vdim+cdim)+1)/2, nc_local.volume);
  triad_basis_gradient_nodal = gkyl_array_new(GKYL_DOUBLE, vdim*vdim*cdim, nc_local.volume);
  double xn[GKYL_MAX_DIM];

  // Choose functions
  metric_t compute_h_ij_from_cov_tangent_basis = choose_metric_kern(vdim); 
  metric_inv_t compute_h_ij_inv_from_h_ij = choose_metric_inv_kern(vdim); 
  metric_det_t compute_det_h_ij = choose_metric_det_kern(vdim); 
  conf_poisson_tensor_t compute_triad_poisson_tensor_ij = choose_conf_poisson_tensor_kern(cdim, vdim); 

  // initialize nodal values in configuration space
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &nc_local);
  while (gkyl_range_iter_next(&iter)) {
    long nidx = gkyl_range_idx(&nc_local, iter.idx);
    
    gkyl_rect_grid_ll_node(cgrid, iter.idx, xn);

    double *pn_cov_tangent_basis = gkyl_array_fetch(cov_tangent_basis_nodal, nidx);
    double *pn_triad_basis = gkyl_array_fetch(triad_basis_nodal, nidx);
    double *pn_h_ij = gkyl_array_fetch(h_ij_nodal, nidx);
    double *pn_h_ij_inv = gkyl_array_fetch(h_ij_inv_nodal, nidx);
    double *pn_det_h = gkyl_array_fetch(det_h_nodal, nidx);
    double *pn_conf_poisson_tensor = gkyl_array_fetch(conf_poisson_tensor_nodal, nidx);
    double *pn_triad_basis_gradient = gkyl_array_fetch(triad_basis_gradient_nodal, nidx);

    // Evaluate the functions for basis and their gradients at a nodal point xn
    inp_basis_vectors.eval_cov_tangent_basis(0.0, xn, pn_cov_tangent_basis, inp_basis_vectors.ctx);
    inp_basis_vectors.eval_triad_basis(0.0, xn, pn_triad_basis, inp_basis_vectors.ctx);
    inp_basis_vectors.eval_triad_basis_gradient(0.0, xn, pn_triad_basis_gradient, inp_basis_vectors.ctx);

    // Compute the metric at nodal points
    compute_h_ij_from_cov_tangent_basis(pn_cov_tangent_basis, pn_h_ij);

    // Compute the metric inverse at nodal points
    compute_h_ij_inv_from_h_ij(pn_h_ij, pn_h_ij_inv);

    // Compute the sqrt(det(h_ij)) at nodal points
    compute_det_h_ij(pn_h_ij, pn_det_h);

    // Compute the Poisson Tensor in configruation space components at nodal points
    compute_triad_poisson_tensor_ij(pn_h_ij_inv, pn_triad_basis, pn_cov_tangent_basis,
      pn_triad_basis_gradient, pn_conf_poisson_tensor);
  }

  // Convert nodal to modal using the computed nodal quantities
  gkyl_eval_on_nodes *h_ij_proj;
  gkyl_eval_on_nodes *h_ij_inv_proj;
  gkyl_eval_on_nodes *det_h_proj;
  gkyl_eval_on_nodes *conf_poisson_tensor_proj;

  h_ij_proj = gkyl_eval_on_nodes_new(cgrid, &cbasis, vdim*(vdim+1)/2, NULL, NULL);
  h_ij_inv_proj = gkyl_eval_on_nodes_new(cgrid, &cbasis, vdim*(vdim+1)/2, NULL, NULL);
  det_h_proj = gkyl_eval_on_nodes_new(cgrid, &cbasis, 1, NULL, NULL);
  conf_poisson_tensor_proj = gkyl_eval_on_nodes_new(cgrid, &cbasis, vdim*(vdim+cdim)*((vdim+cdim)+1)/2, NULL, NULL);

  gkyl_range_iter_init(&iter, crange);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(crange, iter.idx);
    gkyl_eval_on_nodes_nod2mod(h_ij_proj, h_ij_nodal, gkyl_array_fetch(h_ij, lidx));
    gkyl_eval_on_nodes_nod2mod(h_ij_inv_proj, h_ij_inv_nodal, gkyl_array_fetch(h_ij_inv, lidx));
    gkyl_eval_on_nodes_nod2mod(det_h_proj, det_h_nodal, gkyl_array_fetch(det_h, lidx));
    gkyl_eval_on_nodes_nod2mod(conf_poisson_tensor_proj, conf_poisson_tensor_nodal, gkyl_array_fetch(conf_poisson_tensor, lidx));
  }

  gkyl_eval_on_nodes_release(h_ij_proj);
  gkyl_eval_on_nodes_release(h_ij_inv_proj);
  gkyl_eval_on_nodes_release(det_h_proj);
  gkyl_eval_on_nodes_release(conf_poisson_tensor_proj);

  // free temporary memory
  gkyl_array_release(cov_tangent_basis_nodal);
  gkyl_array_release(triad_basis_nodal);
  gkyl_array_release(triad_basis_gradient_nodal);
  gkyl_array_release(h_ij_nodal);
  gkyl_array_release(h_ij_inv_nodal);
  gkyl_array_release(det_h_nodal);
  gkyl_array_release(conf_poisson_tensor_nodal);
}