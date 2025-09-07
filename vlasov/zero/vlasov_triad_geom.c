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


static inline void
log_to_comp(int ndim, const double *eta,
  const double * GKYL_RESTRICT dx, const double * GKYL_RESTRICT xc,
  double* GKYL_RESTRICT xout)
{
  for (int d=0; d<ndim; ++d) xout[d] = 0.5*dx[d]*eta[d]+xc[d];
}

// Identity comp to phys coord mapping, for when user doesn't provide a map.
static inline void
c2p_identity(const double *xcomp, double *xphys, int ndim)
{
  for (int d=0; d<ndim; d++) xphys[d] = xcomp[d];
}


void
gkyl_vlasov_triad_geom_new(const struct gkyl_rect_grid *cgrid, const struct gkyl_range *crange, const struct gkyl_basis cbasis, 
  const struct gkyl_rect_grid *pgrid, const struct gkyl_range *prange, const struct gkyl_basis pbasis, 
  struct gkyl_vlasov_triad_geom_inp inp_basis_vectors, struct gkyl_array *h_ij,
  struct gkyl_array *h_ij_inv, struct gkyl_array *det_h, struct gkyl_array *conf_poisson_tensor)
{

  int num_pt_indices[3] = { 1 , 6, 18 };
  int cdim = cgrid->ndim;
  int pdim = pgrid->ndim;
  int vdim = pdim - cdim;

  // Choose functions
  metric_t compute_h_ij_from_cov_tangent_basis = choose_metric_kern(vdim); 
  metric_inv_t compute_h_ij_inv_from_h_ij = choose_metric_inv_kern(vdim); 
  metric_det_t compute_det_h_ij = choose_metric_det_kern(vdim); 
  conf_poisson_tensor_t compute_triad_poisson_tensor_ij = choose_conf_poisson_tensor_kern(vdim); 

  struct gkyl_array *cov_tangent_basis_at_nodes;
  struct gkyl_array *triad_basis_at_nodes;
  struct gkyl_array *triad_basis_gradient_at_nodes;
  struct gkyl_array *h_ij_at_nodes;
  struct gkyl_array *h_ij_inv_at_nodes;
  struct gkyl_array *det_h_at_nodes;
  struct gkyl_array *conf_poisson_tensor_at_nodes;

  // Convert nodal to modal using the computed nodal quantities
  gkyl_eval_on_nodes *h_ij_proj;
  gkyl_eval_on_nodes *h_ij_inv_proj;
  gkyl_eval_on_nodes *det_h_proj;
  gkyl_eval_on_nodes *conf_poisson_tensor_proj;

  h_ij_proj = gkyl_eval_on_nodes_new(cgrid, &cbasis, vdim*(vdim+1)/2, NULL, NULL);
  h_ij_inv_proj = gkyl_eval_on_nodes_new(cgrid, &cbasis, vdim*(vdim+1)/2, NULL, NULL);
  det_h_proj = gkyl_eval_on_nodes_new(cgrid, &cbasis, 1, NULL, NULL);
  conf_poisson_tensor_proj = gkyl_eval_on_nodes_new(cgrid, &cbasis, num_pt_indices[vdim-1], NULL, NULL);


  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];
  int num_basis = cbasis.num_basis;

  // Eval functions
  cov_tangent_basis_at_nodes = gkyl_array_new(GKYL_DOUBLE, vdim*vdim, num_basis);
  triad_basis_at_nodes = gkyl_array_new(GKYL_DOUBLE, vdim*vdim, num_basis);
  triad_basis_gradient_at_nodes = gkyl_array_new(GKYL_DOUBLE, vdim*vdim*vdim, num_basis);

  // Derived quantities from eval functions
  h_ij_at_nodes = gkyl_array_new(GKYL_DOUBLE, vdim*(vdim+1)/2, num_basis);
  h_ij_inv_at_nodes = gkyl_array_new(GKYL_DOUBLE, vdim*(vdim+1)/2, num_basis);
  det_h_at_nodes = gkyl_array_new(GKYL_DOUBLE, 1, num_basis);
  conf_poisson_tensor_at_nodes = gkyl_array_new(GKYL_DOUBLE, num_pt_indices[vdim-1], num_basis);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, crange);
  
  while (gkyl_range_iter_next(&iter)) {
    gkyl_rect_grid_cell_center(cgrid, iter.idx, xc);

    for (int i=0; i<num_basis; ++i) {

      // Grab the local node list
      struct gkyl_array *nodes = gkyl_array_new(GKYL_DOUBLE, cgrid->ndim, cbasis.num_basis);
      cbasis.node_list(gkyl_array_fetch(nodes, 0));
      log_to_comp(cgrid->ndim, gkyl_array_cfetch(nodes, i),
        cgrid->dx, xc, xmu);
      c2p_identity(xmu, xmu, cgrid->ndim);

      // Evaluate the functions for basis and their gradients at a nodal point xmu
      inp_basis_vectors.eval_cov_tangent_basis(0.0, xmu, gkyl_array_fetch(cov_tangent_basis_at_nodes, i), inp_basis_vectors.eval_cov_tangent_basis_ctx);
      inp_basis_vectors.eval_triad_basis(0.0, xmu, gkyl_array_fetch(triad_basis_at_nodes, i), inp_basis_vectors.eval_triad_basis_ctx);
      inp_basis_vectors.eval_triad_basis_gradient(0.0, xmu, gkyl_array_fetch(triad_basis_gradient_at_nodes, i), inp_basis_vectors.eval_triad_basis_gradient_ctx);

      // Fill the remaining nodal quantities:
      double *pn_cov_tangent_basis = gkyl_array_fetch(cov_tangent_basis_at_nodes, i);
      double *pn_triad_basis = gkyl_array_fetch(triad_basis_at_nodes, i);
      double *pn_h_ij = gkyl_array_fetch(h_ij_at_nodes, i);
      double *pn_h_ij_inv = gkyl_array_fetch(h_ij_inv_at_nodes, i);
      double *pn_det_h = gkyl_array_fetch(det_h_at_nodes, i);
      double *pn_conf_poisson_tensor = gkyl_array_fetch(conf_poisson_tensor_at_nodes, i);
      double *pn_triad_basis_gradient = gkyl_array_fetch(triad_basis_gradient_at_nodes, i);

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

    long lidx = gkyl_range_idx(crange, iter.idx);
    gkyl_eval_on_nodes_nod2mod(h_ij_proj, h_ij_at_nodes, gkyl_array_fetch(h_ij, lidx));
    gkyl_eval_on_nodes_nod2mod(h_ij_inv_proj, h_ij_inv_at_nodes, gkyl_array_fetch(h_ij_inv, lidx));
    gkyl_eval_on_nodes_nod2mod(det_h_proj, det_h_at_nodes, gkyl_array_fetch(det_h, lidx));
    gkyl_eval_on_nodes_nod2mod(conf_poisson_tensor_proj, conf_poisson_tensor_at_nodes, gkyl_array_fetch(conf_poisson_tensor, lidx));
  }


  gkyl_eval_on_nodes_release(h_ij_proj);
  gkyl_eval_on_nodes_release(h_ij_inv_proj);
  gkyl_eval_on_nodes_release(det_h_proj);
  gkyl_eval_on_nodes_release(conf_poisson_tensor_proj);

  // free temporary memory
  gkyl_array_release(cov_tangent_basis_at_nodes);
  gkyl_array_release(triad_basis_at_nodes);
  gkyl_array_release(triad_basis_gradient_at_nodes);
  gkyl_array_release(h_ij_at_nodes);
  gkyl_array_release(h_ij_inv_at_nodes);
  gkyl_array_release(det_h_at_nodes);
  gkyl_array_release(conf_poisson_tensor_at_nodes);
}