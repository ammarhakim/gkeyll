#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_math.h>
#include <gkyl_mirror_grid_gen.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// This file aims to compute all the geometric quantities needed in
// the simulation at nodal values].
struct __attribute__((__packed__)) gkyl_mirror_geo_gen_geom {
  double rza_coord[3]; // Cylindrical coordinate radius, axial length, angle
  double psi; // Psi value at this node

  struct gkyl_vec3 tang[3]; // tangent vectors, e_i cartesian components
  struct gkyl_vec3 dual[3]; // dual vectors, e^i cartesian
  struct gkyl_vec3 normal[3]; // Vector pointing normal vectors n^1, n^2, n^3 cartesian
  double dualmag[3]; // Magnitude of dual vector

  // The metric tensors are stored in order g_11, g_12, g_13, g_22, g_23, g_33
  double metric_covar[6]; // g_ij covariant metric tensor cartesian
  double metric_contr[6]; // g^ij contravariant metric tensor cartesian
  double metric_covar_neut[6]; // g_ij covariant metric tensor for neutrals
  double metric_contr_neut[6]; // g^ij contravariant metric tensor for eutrals

  double Bmag; // Magnitude of magnetic field
  double Bmag_inv; // Inverse of magnetic field
  double Bmag_inv_sq;  // Inverse of magnetic field squared
  struct gkyl_vec3 b_covar; // Covariant components of magnetic field unit vector
  struct gkyl_vec3 b_cart; // Cartesian components of magnetic field unit vector

  double Jc; // Jacobian = e_1*(e_2 X e_3)  = 1/e^1*(e^2 X e^3)
  double Jc_inv; // Jacobian inverse
  double JB; // Jacobian times magnetic field
  double JB_inv; // Jacobian times magnetic field inverse

  double C; // JB/sqrt(g_33)
  double eps2; // epsilon^2 for the poisson solve Jg^33 - J/g_33
};

struct gkyl_mirror_geo_gen {
  struct gkyl_array *nodes_geom; // Geometric quantities at nodes:
  // this is an array of gkyl_mirror_geo_gen_geom objects
};  

struct gkyl_mirror_geo_gen_inp {
  const struct gkyl_rect_grid *comp_grid; // Computational space grid (psi, phi, z)
  struct gkyl_mirror_grid_gen *mirror_grid; // Generated mirror grid
  struct gkyl_range range; // Range of the computational grid
  struct gkyl_basis basis; // Basis for the geometry
};

/**
 * Create new geometry generator for mirror geometries. Return a NULL
 * pointer if the geometry generator failed.
 *
 * @param inp Input parameters to the grid generator
 * @return newly create mirror geometry object
 */
struct gkyl_mirror_geo_gen *gkyl_mirror_geo_gen_inew(const struct gkyl_mirror_geo_gen_inp *inp);

/**
 * Release the mirror geometry object.
 *
 * @param geo Geometry object to release
 */
void gkyl_mirror_geo_gen_release(struct gkyl_mirror_geo_gen *geo);
