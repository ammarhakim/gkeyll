#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_math.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_mirror_grid_gen.h>

struct gkyl_mirror_geo_dg {
  struct gkyl_array *mapc2p; // Map from computational to physical space
  struct gkyl_array *mc2nu_pos; // Map from computational space to field aligned space

  struct gkyl_array *tang; // Cartesian components of tangent vectors stored in order e_1, e_2, e_3
  struct gkyl_array *dual; // Cartesian components of dual vectors stored in order e^1, e^2, e^3
  struct gkyl_array *dualmag; // Norms of the dual vectors: sqrt(e^i.e^i)
  struct gkyl_array *normals; // Cartesian components of normal vectors in order n^1, n^2, n^3

  struct gkyl_array *Jc; // Configuration space Jacobian J
  struct gkyl_array *Jc_inv; // Inverse of configuration space Jacobian 1/J
  struct gkyl_array *JB; // Phase space Jacobian = JB
  struct gkyl_array *JB_inv; // Inverse of phase space Jacobian 1/(JB)

  struct gkyl_array *metric_covar; // Metric coefficients g_{ij} stored in order g_11, g12, g_13, g_22, g_23, g_33
  struct gkyl_array *metric_covar_neut; // Metric coefficients g_{ij} for neutrals stored in order g_11, g12, g_13, g_22, g_23, g_33
  struct gkyl_array *metric_contr; // Metric coefficients g^{ij} stored in order g^{11}, g^{12}, g^{13}, g^{22}, g^{23}, g^{33}
  struct gkyl_array *metric_contr_neut; // Metric coefficients g^{ij} for neutrals stored in order g^{11}, g^{12}, g^{13}, g^{22}, g^{23}, g^{33}

  struct gkyl_array *gxxj; // g^{xx} * J for Poisson solve
  struct gkyl_array *gxyj; // g^{xy} * J for Poisson solve
  struct gkyl_array *gyyj; // g^{yy} * J for Poisson solve
  struct gkyl_array *gxzj; // g^{xz} * J for Poisson solve if z derivatives are kept
  
  struct gkyl_array *Bmag; // B magnitude
  struct gkyl_array *Bmag_inv; // Inverse of B: 1/B
  struct gkyl_array *Bmag_inv_sq; // Inverse of B squared: 1/B^2
  struct gkyl_array *b_covar; // Covariant components of magnetic field vector b_1, b_2, b_3
  struct gkyl_array *b_cart; // Cartesian components of magnetic field vector b_X, b_Y, b_Z
  
  struct gkyl_array *C; // C = JB/sqrt(g_33)
  struct gkyl_array *eps2; // 1 component. eps2 = Jg^33 - J/g_33 for Poisson solve if z derivatives are kept
};  

struct gkyl_mirror_geo_dg_inp {
  struct gkyl_mirror_geo_gen *mirror_geo; // Generated mirror geometry
  
  struct gkyl_range range;
  struct gkyl_range range_ext;
  struct gkyl_basis basis; // Basis for the geometry
  struct gkyl_rect_grid *comp_grid; // Computational space grid (psi, phi, z)
};

/**
 * Convert nodal geometry to DG geometry.
 *
 * @param inp Input parameters to the geometry generator
 * @return newly create mirror geometry object
 */
struct gkyl_mirror_geo_dg *gkyl_mirror_geo_dg_inew(const struct gkyl_mirror_geo_dg_inp *inp);

/**
 * Release the mirror geometry object.
 *
 * @param geo Geometry object to release
 */
void gkyl_mirror_geo_dg_release(struct gkyl_mirror_geo_dg *geo);
