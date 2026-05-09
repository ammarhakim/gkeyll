#include <math.h>
#include <gkyl_dg_gr_maxwell_kernels.h>

// Characteristic-based one-sided outflow flux at a +/- y boundary face for
// 2D serendipity p=1. See outflow_flux_x_2x_ser_p1.c for the rationale and
// the relation to A_dQ_y_calc.
GKYL_CU_DH double
outflow_flux_y_2x_ser_p1(const double *dxv, int edge,
  const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal,
  const double *det_h_nodal, const double *field_con_skin,
  double* GKYL_RESTRICT conf_flux_surf)
{
  double dx11 = 2.0/dxv[1];

  const double *shift_nodal_x = &shift_nodal[0];
  const double *shift_nodal_y = &shift_nodal[2];
  const double *shift_nodal_z = &shift_nodal[4];
  const double *h_xx_nodal = &h_ij_nodal[0];
  const double *h_xy_nodal = &h_ij_nodal[2];
  const double *h_xz_nodal = &h_ij_nodal[4];
  const double *h_yy_nodal = &h_ij_nodal[6];
  const double *h_yz_nodal = &h_ij_nodal[8];
  const double *h_zz_nodal = &h_ij_nodal[10];

  const double *JDx_s = &field_con_skin[0];
  const double *JDy_s = &field_con_skin[4];
  const double *JDz_s = &field_con_skin[8];
  const double *JBx_s = &field_con_skin[12];
  const double *JBy_s = &field_con_skin[16];
  const double *JBz_s = &field_con_skin[20];

  // Evaluate U_skin at the 2 face quadrature nodes.
  // Lower-y face: skin at right of face, evaluate at eta = -1
  //   (matches JD*_con_r_n in dg_gr_maxwell_alpha_quad_y_2x_ser_p1.c).
  // Upper-y face: skin at left of face, evaluate at eta = +1
  //   (matches JD*_con_l_n in the same kernel).
  double JDx_n[2], JDy_n[2], JDz_n[2], JBx_n[2], JBy_n[2], JBz_n[2];
  if (edge == 0) {
    JDx_n[0] = 0.8660254037844386*JDx_s[3]-0.8660254037844386*JDx_s[2]-0.5*JDx_s[1]+0.5*JDx_s[0];
    JDy_n[0] = 0.8660254037844386*JDy_s[3]-0.8660254037844386*JDy_s[2]-0.5*JDy_s[1]+0.5*JDy_s[0];
    JDz_n[0] = 0.8660254037844386*JDz_s[3]-0.8660254037844386*JDz_s[2]-0.5*JDz_s[1]+0.5*JDz_s[0];
    JBx_n[0] = 0.8660254037844386*JBx_s[3]-0.8660254037844386*JBx_s[2]-0.5*JBx_s[1]+0.5*JBx_s[0];
    JBy_n[0] = 0.8660254037844386*JBy_s[3]-0.8660254037844386*JBy_s[2]-0.5*JBy_s[1]+0.5*JBy_s[0];
    JBz_n[0] = 0.8660254037844386*JBz_s[3]-0.8660254037844386*JBz_s[2]-0.5*JBz_s[1]+0.5*JBz_s[0];
    JDx_n[1] = -(0.8660254037844386*JDx_s[3])-0.8660254037844386*JDx_s[2]+0.5*JDx_s[1]+0.5*JDx_s[0];
    JDy_n[1] = -(0.8660254037844386*JDy_s[3])-0.8660254037844386*JDy_s[2]+0.5*JDy_s[1]+0.5*JDy_s[0];
    JDz_n[1] = -(0.8660254037844386*JDz_s[3])-0.8660254037844386*JDz_s[2]+0.5*JDz_s[1]+0.5*JDz_s[0];
    JBx_n[1] = -(0.8660254037844386*JBx_s[3])-0.8660254037844386*JBx_s[2]+0.5*JBx_s[1]+0.5*JBx_s[0];
    JBy_n[1] = -(0.8660254037844386*JBy_s[3])-0.8660254037844386*JBy_s[2]+0.5*JBy_s[1]+0.5*JBy_s[0];
    JBz_n[1] = -(0.8660254037844386*JBz_s[3])-0.8660254037844386*JBz_s[2]+0.5*JBz_s[1]+0.5*JBz_s[0];
  } else {
    JDx_n[0] = -(0.8660254037844386*JDx_s[3])+0.8660254037844386*JDx_s[2]-0.5*JDx_s[1]+0.5*JDx_s[0];
    JDy_n[0] = -(0.8660254037844386*JDy_s[3])+0.8660254037844386*JDy_s[2]-0.5*JDy_s[1]+0.5*JDy_s[0];
    JDz_n[0] = -(0.8660254037844386*JDz_s[3])+0.8660254037844386*JDz_s[2]-0.5*JDz_s[1]+0.5*JDz_s[0];
    JBx_n[0] = -(0.8660254037844386*JBx_s[3])+0.8660254037844386*JBx_s[2]-0.5*JBx_s[1]+0.5*JBx_s[0];
    JBy_n[0] = -(0.8660254037844386*JBy_s[3])+0.8660254037844386*JBy_s[2]-0.5*JBy_s[1]+0.5*JBy_s[0];
    JBz_n[0] = -(0.8660254037844386*JBz_s[3])+0.8660254037844386*JBz_s[2]-0.5*JBz_s[1]+0.5*JBz_s[0];
    JDx_n[1] = 0.8660254037844386*JDx_s[3]+0.8660254037844386*JDx_s[2]+0.5*JDx_s[1]+0.5*JDx_s[0];
    JDy_n[1] = 0.8660254037844386*JDy_s[3]+0.8660254037844386*JDy_s[2]+0.5*JDy_s[1]+0.5*JDy_s[0];
    JDz_n[1] = 0.8660254037844386*JDz_s[3]+0.8660254037844386*JDz_s[2]+0.5*JDz_s[1]+0.5*JDz_s[0];
    JBx_n[1] = 0.8660254037844386*JBx_s[3]+0.8660254037844386*JBx_s[2]+0.5*JBx_s[1]+0.5*JBx_s[0];
    JBy_n[1] = 0.8660254037844386*JBy_s[3]+0.8660254037844386*JBy_s[2]+0.5*JBy_s[1]+0.5*JBy_s[0];
    JBz_n[1] = 0.8660254037844386*JBz_s[3]+0.8660254037844386*JBz_s[2]+0.5*JBz_s[1]+0.5*JBz_s[0];
  }

  double alpha_max = 0.0;
  for (int i = 0; i < 2; ++i) {
    double sx = shift_nodal_x[i], sy = shift_nodal_y[i], sz = shift_nodal_z[i];
    double hxx = h_xx_nodal[i], hxy = h_xy_nodal[i], hxz = h_xz_nodal[i];
    double hyy = h_yy_nodal[i], hyz = h_yz_nodal[i], hzz = h_zz_nodal[i];
    double alpha = lapse_nodal[i];
    double Jc = det_h_nodal[i];

    // Eigenvalues of the y-direction flux Jacobian at this face quad node.
    // Same expressions as dg_gr_maxwell_alpha_quad_y_2x_ser_p1.c lines 244-245.
    double sqrt_term = sqrt(hxx*hzz - hxz*hxz)/Jc;
    double l1 = -sy;
    double l2 = -sy + alpha*sqrt_term;
    double l3 = -sy - alpha*sqrt_term;

    double lambda[6];
    if (edge == 0) {
      // Lower-y boundary: outgoing = -y direction = NEGATIVE eigenvalues.
      lambda[0] = (l1 < 0.0) ? l1 : 0.0;
      lambda[1] = (l1 < 0.0) ? l1 : 0.0;
      lambda[2] = (l2 < 0.0) ? l2 : 0.0;
      lambda[3] = (l2 < 0.0) ? l2 : 0.0;
      lambda[4] = (l3 < 0.0) ? l3 : 0.0;
      lambda[5] = (l3 < 0.0) ? l3 : 0.0;
    } else {
      // Upper-y boundary: outgoing = +y direction = POSITIVE eigenvalues.
      lambda[0] = (l1 > 0.0) ? l1 : 0.0;
      lambda[1] = (l1 > 0.0) ? l1 : 0.0;
      lambda[2] = (l2 > 0.0) ? l2 : 0.0;
      lambda[3] = (l2 > 0.0) ? l2 : 0.0;
      lambda[4] = (l3 > 0.0) ? l3 : 0.0;
      lambda[5] = (l3 > 0.0) ? l3 : 0.0;
    }

    double U_n[6] = { JDx_n[i], JDy_n[i], JDz_n[i], JBx_n[i], JBy_n[i], JBz_n[i] };
    double F_n[6] = { 0.0 };
    A_dQ_y_calc(alpha, sx, sy, sz, hxx, hxy, hxz, hyy, hyz, hzz, Jc, lambda, U_n, F_n);

    conf_flux_surf[0  + i] = F_n[0];
    conf_flux_surf[2  + i] = F_n[1];
    conf_flux_surf[4  + i] = F_n[2];
    conf_flux_surf[6  + i] = F_n[3];
    conf_flux_surf[8  + i] = F_n[4];
    conf_flux_surf[10 + i] = F_n[5];

    double m = fabs(l1);
    if (fabs(l2) > m) m = fabs(l2);
    if (fabs(l3) > m) m = fabs(l3);
    if (m > alpha_max) alpha_max = m;
  }

  conf_flux_surf[12] = 0.0;
  conf_flux_surf[13] = 0.0;
  conf_flux_surf[14] = 0.0;
  conf_flux_surf[15] = 0.0;

  return 1.5*dx11*alpha_max;
}
