#include <math.h>
#include <gkyl_dg_gr_maxwell_kernels.h>

// Curved-norm LLF flux for the y-direction face, 2D serendipity p=1.
// See lax_flux_x_2x_ser_p1_curved.c for the math derivation; this is the
// y-direction analog with the appropriate eta-direction face quadrature
// evaluation and dir=1 conf_flux_surf offset.
GKYL_CU_DH double
lax_flux_y_2x_ser_p1_curved(const double *dxv, const int theta_pole,
      const double *h_ij_nodal, const double *J_c,
      const double *flux_l, const double *flux_r, const double *max_alpha_quad,
      const double *field_con_l, const double *field_con_r,
      double* GKYL_RESTRICT conf_flux_surf)
{
  double dx11 = 2.0/dxv[1];

  const double *h_xx = &h_ij_nodal[0];
  const double *h_xy = &h_ij_nodal[2];
  const double *h_xz = &h_ij_nodal[4];
  const double *h_yy = &h_ij_nodal[6];
  const double *h_yz = &h_ij_nodal[8];
  const double *h_zz = &h_ij_nodal[10];

  double *out = &conf_flux_surf[16]; // dir=1 slice begins at offset 16
  double alpha_max = 0.0;

  // Step 1: compute jumps dQ_n[component][quad_node] (y-direction face quadrature
  // evaluation, matching lax_flux_y_2x_ser_p1.c's formulas).
  double dQ_n[6][2];
  for (int comp = 0; comp < 6; ++comp) {
    const double *em_l = &field_con_l[comp*4];
    const double *em_r = &field_con_r[comp*4];

    double Q_l_n0 = -(0.8660254037844386*em_l[3])+0.8660254037844386*em_l[2]-0.5*em_l[1]+0.5*em_l[0];
    double Q_r_n0 =  (0.8660254037844386*em_r[3])-0.8660254037844386*em_r[2]-0.5*em_r[1]+0.5*em_r[0];
    dQ_n[comp][0] = Q_r_n0 - Q_l_n0;

    double Q_l_n1 = (0.8660254037844386*em_l[3])+0.8660254037844386*em_l[2]+0.5*em_l[1]+0.5*em_l[0];
    double Q_r_n1 = -(0.8660254037844386*em_r[3])-0.8660254037844386*em_r[2]+0.5*em_r[1]+0.5*em_r[0];
    dQ_n[comp][1] = Q_r_n1 - Q_l_n1;
  }

  // Step 2: metric-weighted dissipation per quad node.
  for (int n = 0; n < 2; ++n) {
    alpha_max = fmax(alpha_max, fabs(max_alpha_quad[n]));

    double diss[6];
    diss[0] = h_xx[n]*dQ_n[0][n] + h_xy[n]*dQ_n[1][n] + h_xz[n]*dQ_n[2][n];
    diss[1] = h_xy[n]*dQ_n[0][n] + h_yy[n]*dQ_n[1][n] + h_yz[n]*dQ_n[2][n];
    diss[2] = h_xz[n]*dQ_n[0][n] + h_yz[n]*dQ_n[1][n] + h_zz[n]*dQ_n[2][n];
    diss[3] = h_xx[n]*dQ_n[3][n] + h_xy[n]*dQ_n[4][n] + h_xz[n]*dQ_n[5][n];
    diss[4] = h_xy[n]*dQ_n[3][n] + h_yy[n]*dQ_n[4][n] + h_yz[n]*dQ_n[5][n];
    diss[5] = h_xz[n]*dQ_n[3][n] + h_yz[n]*dQ_n[4][n] + h_zz[n]*dQ_n[5][n];

    double a = fabs(max_alpha_quad[n]);
    for (int comp = 0; comp < 6; ++comp) {
      const double *flux_l_q = &flux_l[comp*2];
      const double *flux_r_q = &flux_r[comp*2];
      if (theta_pole) {
        out[n + comp*2] = 0.5*(flux_r_q[n] + flux_l_q[n]);
      } else {
        out[n + comp*2] = 0.5*((flux_r_q[n] + flux_l_q[n]) - a*diss[comp]);
      }
    }
  }

  return 1.5*dx11*alpha_max;
}
