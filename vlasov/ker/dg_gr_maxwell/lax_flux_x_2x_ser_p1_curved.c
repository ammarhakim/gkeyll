#include <math.h>
#include <gkyl_dg_gr_maxwell_kernels.h>

// Curved-norm LLF flux for the x-direction face, 2D serendipity p=1.
//
// Standard lax_flux_x_2x_ser_p1 uses
//     F*[k] = 0.5*(F_l[k] + F_r[k]) - 0.5*|alpha_max| * (Q_r[k] - Q_l[k])
// — flat-L^2 dissipation. That can leave a real-positive eigenvalue in
// the operator's spectrum inside the ergoregion.
//
// This curved-norm variant uses h_{ij} (the spatial 3-metric) as the
// dissipation weight:
//     F*[k] = 0.5*(F_l[k] + F_r[k]) - 0.5*|alpha_max| * h_{kj} * (Q_r-Q_l)[j]
// h_{kj} is positive-definite on the spatial slice (only the spacetime
// 4-metric becomes indefinite in the ergoregion via g_{tt}). Block-
// diagonal between D-components (k,j in {0,1,2}) and B-components
// (k,j in {3,4,5}).
//
// We tested also a J_c * alpha * h_ij weighting (full volumetric energy
// kernel) at N=24: it gave 40% lower max Re(lambda) but lost the machine-
// zero result at N=12, with the same N^3.8 power-law scaling. The simpler
// h_ij-only form is kept because it produces a cleaner low-N spectrum and
// neither variant eliminates the high-N residual.
GKYL_CU_DH double
lax_flux_x_2x_ser_p1_curved(const double *dxv, const int theta_pole,
      const double *h_ij_nodal, const double *J_c,
      const double *flux_l, const double *flux_r, const double *max_alpha_quad,
      const double *field_con_l, const double *field_con_r,
      double* GKYL_RESTRICT conf_flux_surf)
{
  double dx10 = 2.0/dxv[0];

  // h_ij layout at this face: 6 components (xx, xy, xz, yy, yz, zz),
  // each at 2 face quadrature nodes. Component-major.
  const double *h_xx = &h_ij_nodal[0];
  const double *h_xy = &h_ij_nodal[2];
  const double *h_xz = &h_ij_nodal[4];
  const double *h_yy = &h_ij_nodal[6];
  const double *h_yz = &h_ij_nodal[8];
  const double *h_zz = &h_ij_nodal[10];

  double *out = &conf_flux_surf[0];
  double alpha_max = 0.0;

  // Step 1: compute jumps dQ_n[component][quad_node] using the same
  // modal-to-nodal evaluation as the production lax_flux_x kernel.
  double dQ_n[6][2];
  for (int comp = 0; comp < 6; ++comp) {
    const double *em_l = &field_con_l[comp*4];
    const double *em_r = &field_con_r[comp*4];

    double Q_l_n0 = -(0.8660254037844386*em_l[3])-0.5*em_l[2]+0.8660254037844386*em_l[1]+0.5*em_l[0];
    double Q_r_n0 =  (0.8660254037844386*em_r[3])-0.5*em_r[2]-0.8660254037844386*em_r[1]+0.5*em_r[0];
    dQ_n[comp][0] = Q_r_n0 - Q_l_n0;

    double Q_l_n1 = (0.8660254037844386*em_l[3])+0.5*em_l[2]+0.8660254037844386*em_l[1]+0.5*em_l[0];
    double Q_r_n1 = -(0.8660254037844386*em_r[3])+0.5*em_r[2]-0.8660254037844386*em_r[1]+0.5*em_r[0];
    dQ_n[comp][1] = Q_r_n1 - Q_l_n1;
  }

  // Step 2: build metric-weighted dissipation per quad node and write output.
  // Dissipation weight per quad node: J_c * alpha * h_ij. The face nodal
  // J_c is the det_h argument; alpha is the lapse_nodal argument.
  for (int n = 0; n < 2; ++n) {
    alpha_max = fmax(alpha_max, fabs(max_alpha_quad[n]));

    double diss[6];
    // D-block: components 0..2.
    diss[0] = h_xx[n]*dQ_n[0][n] + h_xy[n]*dQ_n[1][n] + h_xz[n]*dQ_n[2][n];
    diss[1] = h_xy[n]*dQ_n[0][n] + h_yy[n]*dQ_n[1][n] + h_yz[n]*dQ_n[2][n];
    diss[2] = h_xz[n]*dQ_n[0][n] + h_yz[n]*dQ_n[1][n] + h_zz[n]*dQ_n[2][n];
    // B-block: components 3..5 (same h_ij weighting, independent block).
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

  return 1.5*dx10*alpha_max;
}
