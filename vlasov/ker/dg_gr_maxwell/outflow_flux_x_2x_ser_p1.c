#include <math.h>
#include <stdio.h>
#include <gkyl_dg_gr_maxwell_kernels.h>

// Characteristic-based one-sided outflow flux at a +/- x boundary face for
// 2D serendipity p=1.
//
// At a domain boundary the standard LLF/Roe flux mixes the skin and ghost
// states: any reflective ghost (fixed-function around a non-zero background,
// copy, etc.) leaks information back into the skin via either the dissipative
// term (LLF) or the ghost-side characteristic projection (Roe). For a clean
// non-reflecting boundary we want
//
//   F_BC = A^outgoing . U_skin
//
// where A^outgoing is the flux Jacobian projected onto the outgoing
// characteristics (positive eigenvalues at the upper edge, negative at the
// lower edge). This is computed wave-by-wave using A_dQ_x_calc, which is
// linear in its dQ argument; passing U_skin in place of dQ gives A . U_skin
// for whatever eigenvalue spectrum is supplied.
//
// `edge` is GKYL_LOWER_EDGE (=0) for the lower-x face (ghost-on-left of skin)
// or GKYL_UPPER_EDGE (=1) for the upper-x face (ghost-on-right of skin).
GKYL_CU_DH double
outflow_flux_x_2x_ser_p1(const double *dxv, int edge,
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

  // Evaluate U_skin at the 2 face quadrature nodes.
  // Lower-x face: skin is at the right of the face, evaluate at xi = -1
  //   (matches JD*_con_r_n in dg_gr_maxwell_alpha_quad_x_2x_ser_p1.c).
  // Upper-x face: skin is at the left of the face, evaluate at xi = +1
  //   (matches JD*_con_l_n in the same kernel).
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

  // [INSTRUMENTATION] Print eigenvalues, U_skin, and resulting flux on the
  // first invocation per (edge) so we can confirm what the kernel actually
  // sees from a real simulation step. Bounded to 2 calls (one per edge).
  static int dbg_n_calls[2] = { 0, 0 };
  int do_dbg = (edge >= 0 && edge < 2 && dbg_n_calls[edge] == 0);

  double alpha_max = 0.0;
  for (int i = 0; i < 2; ++i) {
    double sx = shift_nodal_x[i], sy = shift_nodal_y[i], sz = shift_nodal_z[i];
    double hxx = h_xx_nodal[i], hxy = h_xy_nodal[i], hxz = h_xz_nodal[i];
    double hyy = h_yy_nodal[i], hyz = h_yz_nodal[i], hzz = h_zz_nodal[i];
    double alpha = lapse_nodal[i];
    double Jc = det_h_nodal[i];

    // Eigenvalues of the x-direction flux Jacobian at this quad node.
    // lambda_1 is the (multiplicity-2) advection at -beta^x; lambda_2/3 are
    // the fast/slow EM-wave pair. Same expressions as in
    // dg_gr_maxwell_alpha_quad_x_2x_ser_p1.c, lines 244-245.
    double sqrt_term = sqrt(hyy*hzz - hyz*hyz)/Jc;
    double l1 = -sx;
    double l2 = -sx + alpha*sqrt_term;
    double l3 = -sx - alpha*sqrt_term;

    if (do_dbg) {
      fprintf(stderr, "[outflow_flux_x ser_p1] FIRST CALL edge=%d quad_node=%d:\n"
        "    metric: alpha=%1.6e  shift=(%1.6e, %1.6e, %1.6e)  Jc=%1.6e\n"
        "    h_ij = [xx=%1.6e xy=%1.6e xz=%1.6e yy=%1.6e yz=%1.6e zz=%1.6e]\n"
        "    eigenvalues: l1=%+1.6e  l2=%+1.6e  l3=%+1.6e\n"
        "    U_skin@node = (JD: %+1.6e %+1.6e %+1.6e | JB: %+1.6e %+1.6e %+1.6e)\n",
        edge, i, alpha, sx, sy, sz, Jc,
        hxx, hxy, hxz, hyy, hyz, hzz,
        l1, l2, l3,
        JDx_n[i], JDy_n[i], JDz_n[i], JBx_n[i], JBy_n[i], JBz_n[i]);
    }

    // Build outgoing-only eigenvalue vector for A_dQ_x_calc.
    // Indexing in A_dQ_x_calc's lambda[6]: [D-l1, B-l1, D-l2, B-l2, D-l3, B-l3].
    double lambda[6];
    if (edge == 0) {
      // Lower-x boundary: outgoing waves leave the domain in the -x direction,
      // i.e., have NEGATIVE characteristic speeds. Keep only those.
      lambda[0] = (l1 < 0.0) ? l1 : 0.0;
      lambda[1] = (l1 < 0.0) ? l1 : 0.0;
      lambda[2] = (l2 < 0.0) ? l2 : 0.0;
      lambda[3] = (l2 < 0.0) ? l2 : 0.0;
      lambda[4] = (l3 < 0.0) ? l3 : 0.0;
      lambda[5] = (l3 < 0.0) ? l3 : 0.0;
    } else {
      // Upper-x boundary: outgoing = +x direction = POSITIVE eigenvalues.
      lambda[0] = (l1 > 0.0) ? l1 : 0.0;
      lambda[1] = (l1 > 0.0) ? l1 : 0.0;
      lambda[2] = (l2 > 0.0) ? l2 : 0.0;
      lambda[3] = (l2 > 0.0) ? l2 : 0.0;
      lambda[4] = (l3 > 0.0) ? l3 : 0.0;
      lambda[5] = (l3 > 0.0) ? l3 : 0.0;
    }

    // F_outflow = A^outgoing . U_skin, computed wave-by-wave by reusing
    // A_dQ_x_calc's machinery (linear in its dQ argument).
    double U_n[6] = { JDx_n[i], JDy_n[i], JDz_n[i], JBx_n[i], JBy_n[i], JBz_n[i] };
    double F_n[6] = { 0.0 };
    A_dQ_x_calc(alpha, sx, sy, sz, hxx, hxy, hxz, hyy, hyz, hzz, Jc, lambda, U_n, F_n);

    if (do_dbg) {
      fprintf(stderr, "    lambda_outgoing = (%+1.6e %+1.6e %+1.6e %+1.6e %+1.6e %+1.6e)\n"
        "    F_outflow@node  = (%+1.6e %+1.6e %+1.6e | %+1.6e %+1.6e %+1.6e)\n",
        lambda[0], lambda[1], lambda[2], lambda[3], lambda[4], lambda[5],
        F_n[0], F_n[1], F_n[2], F_n[3], F_n[4], F_n[5]);
    }

    // Write into the conf_flux_surf slot (component-major: 6 components,
    // 2 quad nodes), matching the layout produced by lax_flux_x_2x_ser_p1.
    conf_flux_surf[0  + i] = F_n[0];
    conf_flux_surf[2  + i] = F_n[1];
    conf_flux_surf[4  + i] = F_n[2];
    conf_flux_surf[6  + i] = F_n[3];
    conf_flux_surf[8  + i] = F_n[4];
    conf_flux_surf[10 + i] = F_n[5];

    // Track max |lambda| for cflrate.
    double m = fabs(l1);
    if (fabs(l2) > m) m = fabs(l2);
    if (fabs(l3) > m) m = fabs(l3);
    if (m > alpha_max) alpha_max = m;
  }
  if (do_dbg) dbg_n_calls[edge]++;

  // Cleaning slots (phi, psi) carry no outflow contribution: leave at zero.
  // Caller ensures the array is zeroed before the per-direction kernel runs;
  // we explicitly zero here too so that downstream surfx kernels see a clean
  // slot regardless of caller convention.
  conf_flux_surf[12] = 0.0;
  conf_flux_surf[13] = 0.0;
  conf_flux_surf[14] = 0.0;
  conf_flux_surf[15] = 0.0;

  // Same cflrate convention as lax_flux_x_2x_ser_p1: 1.5 * dx10 * |alpha_max|.
  return 1.5*dx10*alpha_max;
}
