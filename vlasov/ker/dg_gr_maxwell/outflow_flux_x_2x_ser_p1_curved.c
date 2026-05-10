#include <math.h>
#include <gkyl_dg_gr_maxwell_kernels.h>

// Curved-norm outflow flux at an x-direction boundary face, 2D ser p=1.
//
// The standard outflow kernel returns F_BC = A^outgoing . U_skin, which is
// the upwind characteristic flux at a boundary with no incoming wave. By
// construction it has no explicit dissipative term, since A^outgoing is
// already a one-sided projection. dgeev shows however that pairing this
// outflow with curved-norm LLF in the interior leaves a small residual
// real-positive eigenvalue (~3e-3 at N=12); the contribution comes from
// the boundary face's flux not damping in the curved-space energy norm.
//
// This curved-norm variant adds a metric-weighted dissipation term on top
// of the standard outflow output:
//
//   F_curved[k] = (A^outgoing . U_skin)[k] - 0.5 * |alpha_max| * h_kj * U_skin^j
//
// The dissipation matrix h_kj is the spatial 3-metric, block-diagonal
// between D-block (k,j in {0..2}) and B-block (k,j in {3..5}), and is
// positive-definite on the spatial slice. The 0.5 factor matches the
// standard LLF dissipation magnitude. |alpha_max| is the maximum
// characteristic speed at the face (max |l1|, |l2|, |l3|).
//
// Effect: a positive-definite damping of U_skin in the curved-space norm
// at the boundary face. Outgoing characteristics get a small extra damp
// (slightly over-dissipative compared to a true non-reflecting BC), but
// the trade-off is structural curved-norm stability.
GKYL_CU_DH double
outflow_flux_x_2x_ser_p1_curved(const double *dxv, int edge,
  const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal,
  const double *det_h_nodal, const double *field_con_skin,
  double* GKYL_RESTRICT conf_flux_surf)
{
  double dx10 = 2.0/dxv[0];

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

  // Same modal-to-nodal evaluation as outflow_flux_x_2x_ser_p1.c.
  double JDx_n[2], JDy_n[2], JDz_n[2], JBx_n[2], JBy_n[2], JBz_n[2];
  if (edge == 0) {
    JDx_n[0] = 0.8660254037844386*JDx_s[3]-0.5*JDx_s[2]-0.8660254037844386*JDx_s[1]+0.5*JDx_s[0];
    JDy_n[0] = 0.8660254037844386*JDy_s[3]-0.5*JDy_s[2]-0.8660254037844386*JDy_s[1]+0.5*JDy_s[0];
    JDz_n[0] = 0.8660254037844386*JDz_s[3]-0.5*JDz_s[2]-0.8660254037844386*JDz_s[1]+0.5*JDz_s[0];
    JBx_n[0] = 0.8660254037844386*JBx_s[3]-0.5*JBx_s[2]-0.8660254037844386*JBx_s[1]+0.5*JBx_s[0];
    JBy_n[0] = 0.8660254037844386*JBy_s[3]-0.5*JBy_s[2]-0.8660254037844386*JBy_s[1]+0.5*JBy_s[0];
    JBz_n[0] = 0.8660254037844386*JBz_s[3]-0.5*JBz_s[2]-0.8660254037844386*JBz_s[1]+0.5*JBz_s[0];
    JDx_n[1] = -(0.8660254037844386*JDx_s[3])+0.5*JDx_s[2]-0.8660254037844386*JDx_s[1]+0.5*JDx_s[0];
    JDy_n[1] = -(0.8660254037844386*JDy_s[3])+0.5*JDy_s[2]-0.8660254037844386*JDy_s[1]+0.5*JDy_s[0];
    JDz_n[1] = -(0.8660254037844386*JDz_s[3])+0.5*JDz_s[2]-0.8660254037844386*JDz_s[1]+0.5*JDz_s[0];
    JBx_n[1] = -(0.8660254037844386*JBx_s[3])+0.5*JBx_s[2]-0.8660254037844386*JBx_s[1]+0.5*JBx_s[0];
    JBy_n[1] = -(0.8660254037844386*JBy_s[3])+0.5*JBy_s[2]-0.8660254037844386*JBy_s[1]+0.5*JBy_s[0];
    JBz_n[1] = -(0.8660254037844386*JBz_s[3])+0.5*JBz_s[2]-0.8660254037844386*JBz_s[1]+0.5*JBz_s[0];
  } else {
    JDx_n[0] = -(0.8660254037844386*JDx_s[3])-0.5*JDx_s[2]+0.8660254037844386*JDx_s[1]+0.5*JDx_s[0];
    JDy_n[0] = -(0.8660254037844386*JDy_s[3])-0.5*JDy_s[2]+0.8660254037844386*JDy_s[1]+0.5*JDy_s[0];
    JDz_n[0] = -(0.8660254037844386*JDz_s[3])-0.5*JDz_s[2]+0.8660254037844386*JDz_s[1]+0.5*JDz_s[0];
    JBx_n[0] = -(0.8660254037844386*JBx_s[3])-0.5*JBx_s[2]+0.8660254037844386*JBx_s[1]+0.5*JBx_s[0];
    JBy_n[0] = -(0.8660254037844386*JBy_s[3])-0.5*JBy_s[2]+0.8660254037844386*JBy_s[1]+0.5*JBy_s[0];
    JBz_n[0] = -(0.8660254037844386*JBz_s[3])-0.5*JBz_s[2]+0.8660254037844386*JBz_s[1]+0.5*JBz_s[0];
    JDx_n[1] = 0.8660254037844386*JDx_s[3]+0.5*JDx_s[2]+0.8660254037844386*JDx_s[1]+0.5*JDx_s[0];
    JDy_n[1] = 0.8660254037844386*JDy_s[3]+0.5*JDy_s[2]+0.8660254037844386*JDy_s[1]+0.5*JDy_s[0];
    JDz_n[1] = 0.8660254037844386*JDz_s[3]+0.5*JDz_s[2]+0.8660254037844386*JDz_s[1]+0.5*JDz_s[0];
    JBx_n[1] = 0.8660254037844386*JBx_s[3]+0.5*JBx_s[2]+0.8660254037844386*JBx_s[1]+0.5*JBx_s[0];
    JBy_n[1] = 0.8660254037844386*JBy_s[3]+0.5*JBy_s[2]+0.8660254037844386*JBy_s[1]+0.5*JBy_s[0];
    JBz_n[1] = 0.8660254037844386*JBz_s[3]+0.5*JBz_s[2]+0.8660254037844386*JBz_s[1]+0.5*JBz_s[0];
  }

  double alpha_max = 0.0;
  for (int i = 0; i < 2; ++i) {
    double sx = shift_nodal_x[i], sy = shift_nodal_y[i], sz = shift_nodal_z[i];
    double hxx = h_xx_nodal[i], hxy = h_xy_nodal[i], hxz = h_xz_nodal[i];
    double hyy = h_yy_nodal[i], hyz = h_yz_nodal[i], hzz = h_zz_nodal[i];
    double alpha = lapse_nodal[i];
    double Jc = det_h_nodal[i];

    double sqrt_term = sqrt(hyy*hzz - hyz*hyz)/Jc;
    double l1 = -sx;
    double l2 = -sx + alpha*sqrt_term;
    double l3 = -sx - alpha*sqrt_term;

    double lambda[6];
    if (edge == 0) {
      lambda[0] = (l1 < 0.0) ? l1 : 0.0;
      lambda[1] = (l1 < 0.0) ? l1 : 0.0;
      lambda[2] = (l2 < 0.0) ? l2 : 0.0;
      lambda[3] = (l2 < 0.0) ? l2 : 0.0;
      lambda[4] = (l3 < 0.0) ? l3 : 0.0;
      lambda[5] = (l3 < 0.0) ? l3 : 0.0;
    } else {
      lambda[0] = (l1 > 0.0) ? l1 : 0.0;
      lambda[1] = (l1 > 0.0) ? l1 : 0.0;
      lambda[2] = (l2 > 0.0) ? l2 : 0.0;
      lambda[3] = (l2 > 0.0) ? l2 : 0.0;
      lambda[4] = (l3 > 0.0) ? l3 : 0.0;
      lambda[5] = (l3 > 0.0) ? l3 : 0.0;
    }

    // Standard outflow contribution: F_outgoing = A^outgoing . U_skin.
    double U_n[6] = { JDx_n[i], JDy_n[i], JDz_n[i], JBx_n[i], JBy_n[i], JBz_n[i] };
    double F_n[6] = { 0.0 };
    A_dQ_x_calc(alpha, sx, sy, sz, hxx, hxy, hxz, hyy, hyz, hzz, Jc, lambda, U_n, F_n);

    // Curved-norm dissipation correction. Weight is J_c * alpha * h_ij —
    // the kernel of the volumetric curved-space energy form,
    //   E = (1/2) Integral [J_c * alpha * h_ij * (D^iD^j + B^iB^j)] d^3x.
    //
    // The sign matches LLF's
    //   F_diss = -0.5*|alpha_max| * (J_c*alpha*h_ij) * (U_R - U_L)^j
    // with the implicit ghost = 0 convention:
    //   lower edge (U_L = ghost = 0, U_R = U_skin):
    //     F_diss = -0.5 * |alpha_max| * J_c*alpha*h * U_skin
    //   upper edge (U_L = U_skin, U_R = ghost = 0):
    //     F_diss = +0.5 * |alpha_max| * J_c*alpha*h * U_skin
    // The sign flip at the upper edge is essential: the boundary face
    // contributes +F to skin's RHS at the lower edge and -F at the upper
    // edge, so the dissipation has to reverse to stay damping in both
    // cases. Without the sign flip, the upper boundary's "dissipation"
    // is destabilizing and blows up the spectrum.
    double a_face = fabs(l1);
    if (fabs(l2) > a_face) a_face = fabs(l2);
    if (fabs(l3) > a_face) a_face = fabs(l3);
    double sign = (edge == 0) ? -1.0 : +1.0;

    double damp_D[3];
    damp_D[0] = hxx*U_n[0] + hxy*U_n[1] + hxz*U_n[2];
    damp_D[1] = hxy*U_n[0] + hyy*U_n[1] + hyz*U_n[2];
    damp_D[2] = hxz*U_n[0] + hyz*U_n[1] + hzz*U_n[2];
    double damp_B[3];
    damp_B[0] = hxx*U_n[3] + hxy*U_n[4] + hxz*U_n[5];
    damp_B[1] = hxy*U_n[3] + hyy*U_n[4] + hyz*U_n[5];
    damp_B[2] = hxz*U_n[3] + hyz*U_n[4] + hzz*U_n[5];

    F_n[0] += sign*0.5*a_face*damp_D[0];
    F_n[1] += sign*0.5*a_face*damp_D[1];
    F_n[2] += sign*0.5*a_face*damp_D[2];
    F_n[3] += sign*0.5*a_face*damp_B[0];
    F_n[4] += sign*0.5*a_face*damp_B[1];
    F_n[5] += sign*0.5*a_face*damp_B[2];

    conf_flux_surf[0  + i] = F_n[0];
    conf_flux_surf[2  + i] = F_n[1];
    conf_flux_surf[4  + i] = F_n[2];
    conf_flux_surf[6  + i] = F_n[3];
    conf_flux_surf[8  + i] = F_n[4];
    conf_flux_surf[10 + i] = F_n[5];

    if (a_face > alpha_max) alpha_max = a_face;
  }

  conf_flux_surf[12] = 0.0;
  conf_flux_surf[13] = 0.0;
  conf_flux_surf[14] = 0.0;
  conf_flux_surf[15] = 0.0;

  return 1.5*dx10*alpha_max;
}
