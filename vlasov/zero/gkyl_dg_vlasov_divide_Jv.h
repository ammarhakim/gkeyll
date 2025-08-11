#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>

/**
 * Divide out the velocity-space Jacobian from Jf to obtain f for use in the DG update. 
 * This update is performed at Gauss-Legendre quadrature points because, by construction,
 * we know Jv at these points and our Jv is static in time, so we can obtain the weakly
 * equivalent f from Jf in our DG expansion (since we can always transform between the
 * Gauss-Legendre nodal basis and our modal basis). 
 * 
 * @param conf_basis Configuration-space basis. 
 * @param phase_basis Phase-space basis. 
 * @param vel_range Velocity-space range. 
 * @param phase_range Phase-space range. 
 * @param jacob_vel_gauss Full (up to 3V) velocity-space Jacobian at Gauss-Legendre quadrature points. 
 * @param Jf Input array Jf. 
 * @param f_no_J Output array f with velocity-space Jacobian divided out.  
 * @param use_gpu bool to determine if on GPU
 */
void gkyl_dg_vlasov_divide_Jv(const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  const struct gkyl_range *vel_range, const struct gkyl_range *phase_range, 
  const struct gkyl_array *jacob_vel_gauss, const struct gkyl_array *Jf, 
  struct gkyl_array *f_no_J, bool use_gpu);

/**
 * Divide out the velocity-space Jacobian from Jf to obtain f on
 * NV-GPU. See new() method for documentation.
 */  
void gkyl_dg_vlasov_divide_Jv_cu(const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  const struct gkyl_range *vel_range, const struct gkyl_range *phase_range, 
  const struct gkyl_array *jacob_vel_gauss, const struct gkyl_array *Jf, 
  struct gkyl_array *f_no_J);  

/**
 * Multiply the velocity-space Jacobian by f to obtain Jf use in the DG update. 
 * This update is performed at Gauss-Legendre quadrature points because, by construction,
 * we know Jv at these points and our Jv is static in time, so we can obtain the weakly
 * equivalent f from Jf in our DG expansion (since we can always transform between the
 * Gauss-Legendre nodal basis and our modal basis). 
 * 
 * @param conf_basis Configuration-space basis. 
 * @param phase_basis Phase-space basis. 
 * @param vel_range Velocity-space range. 
 * @param phase_range Phase-space range. 
 * @param jacob_vel_gauss Full (up to 3V) velocity-space Jacobian at Gauss-Legendre quadrature points. 
 * @param f_no_J Input array f with velocity-space Jacobian divided out.  
 * @param Jf Output array Jf. 
 * @param use_gpu bool to determine if on GPU
 */
void gkyl_dg_vlasov_rescale_Jv(const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  const struct gkyl_range *vel_range, const struct gkyl_range *phase_range, 
  const struct gkyl_array *jacob_vel_gauss, const struct gkyl_array *f_no_J, 
  struct gkyl_array *Jf, bool use_gpu);

/**
 * Multiply the velocity-space Jacobian by f to obtain Jf on
 * NV-GPU. See new() method for documentation.
 */  
void gkyl_dg_vlasov_rescale_Jv_cu(const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  const struct gkyl_range *vel_range, const struct gkyl_range *phase_range, 
  const struct gkyl_array *jacob_vel_gauss, const struct gkyl_array *f_no_J, 
  struct gkyl_array *Jf);    