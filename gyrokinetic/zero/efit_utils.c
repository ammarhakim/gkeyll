#include <gkyl_efit_priv.h>
#include <float.h>
#include <math.h>

static bool newton_raphson_seed(struct gkyl_efit *up, const double *coeffs,
  const double *xinit, double *xsol, bool cubics);

static inline double
gkyl_sq(double x)
{
  return x*x;
}

static inline double
gkyl_dist2(double R0, double Z0, double R1, double Z1)
{
  return gkyl_sq(R0-R1) + gkyl_sq(Z0-Z1);
}

static inline double
gkyl_clamp(double x, double lo, double hi)
{
  return x < lo ? lo : (x > hi ? hi : x);
}

static double
eval_laplacian_expand_2d_tensor_p2(int dir, const double *z, const double *f )
{
  const double z0 = z[0];
  const double z1 = z[1];
  if (dir == 0)
    return 11.25*f[8]*z1*z1+5.809475019311125*f[6]*z1-3.75*f[8]+3.354101966249685*f[4];
  if (dir == 1)
    return 11.25*f[8]*z0*z0+5.809475019311125*f[7]*z0-3.75*f[8]+3.354101966249685*f[5];

  return 0.0; // can't happen, suppresses warning

}

static double
eval_mixedpartial_expand_2d_tensor_p2(const double *z, const double *f )
{
  const double z0 = z[0];
  const double z1 = z[1];
  return 22.5*f[8]*z0*z1+5.809475019311125*f[7]*z1+5.809475019311125*f[6]*z0+1.5*f[3]; 

}

static void print_result(int n, double x[], double dx[], double errx, double errf, int niter)
{
  double x0 = x[0];
  double y0 = x[1];
  if (x0 >= -1 && x0 <= 1 && y0 >= -1 && y0 <= 1 ) {
    printf("x = ");
    for(int i=0; i<n; i++) printf(" %g ",  x[i]);
    printf("\n");
    printf("dx = ");
    for(int i=0; i<n; i++) printf("= %g ",  dx[i]);
    printf("\n");
    printf("errx = %g\n", errx);
    printf("errf = %g\n", errf);
  }

}

bool 
newton_raphson(struct gkyl_efit *up, const double *coeffs, double *xsol, bool cubics)
{
  return newton_raphson_seed(up, coeffs, (double[2]){0.0, 0.0}, xsol, cubics);
}

static bool
newton_raphson_seed(struct gkyl_efit *up, const double *coeffs, const double *xinit, double *xsol, bool cubics)
{
  int n = 2;
  double x[2] = {xinit[0], xinit[1]};
  double dx[2] = {0.0,0.0};
  double fjac[2][2];
  double fjac_inv[2][2];
  double fvec[2];
  double p[2];
  int ntrial = 100;
  double errx = 0.0;
  double errf = 0.0;
  for (int i=0;i<n;i++) xsol[i] = 0.0;

  for (int niter = 0; niter < ntrial; niter++) {
    if (cubics) {
      for(int i=0; i<n; i++) fvec[i] = up->rzbasis_cubic.eval_grad_expand(i,x,coeffs);
      fjac[0][0] = up->evf->eval_cubic_laplacian(0,x,coeffs);
      fjac[0][1] = up->evf->eval_cubic_mixedpartial(x,coeffs);
      fjac[1][0] = up->evf->eval_cubic_mixedpartial(x,coeffs);
      fjac[1][1] = up->evf->eval_cubic_laplacian(1,x,coeffs);
    }
    else {
      for(int i=0; i<n; i++) fvec[i] = up->rzbasis.eval_grad_expand(i,x,coeffs);
      fjac[0][0] = eval_laplacian_expand_2d_tensor_p2(0,x,coeffs);
      fjac[0][1] = eval_mixedpartial_expand_2d_tensor_p2(x,coeffs);
      fjac[1][0] = eval_mixedpartial_expand_2d_tensor_p2(x,coeffs);
      fjac[1][1] = eval_laplacian_expand_2d_tensor_p2(1,x,coeffs);
    }
    errf = 0.0;
    for (int i=0;i<n;i++) errf += fvec[i]*fvec[i];
    errf = sqrt(errf);
    if (errf <= 1e-18) {
      for (int i=0;i<n;i++) xsol[i] = x[i];
      return true;
    }
    for (int i=0;i<n;i++) p[i] = -fvec[i];

    double det = fjac[0][0]*fjac[1][1] - fjac[0][1]*fjac[1][0];
    fjac_inv[0][0] = fjac[1][1]/det;
    fjac_inv[1][1] = fjac[0][0]/det;
    fjac_inv[0][1] = -fjac[0][1]/det;
    fjac_inv[1][0] = -fjac[1][0]/det;

    dx[0] = fjac_inv[0][0]*p[0] + fjac_inv[0][1]*p[1];
    dx[1] = fjac_inv[1][0]*p[0] + fjac_inv[1][1]*p[1];

    errx = 0.0;
    for (int i=0;i<n;i++) {
      errx += fabs(dx[i]);
      x[i] += dx[i];
    }

    if (errx<=1e-18) {
      for (int i=0;i<n;i++) xsol[i] = x[i];
      return true;
    }

  }
  return false;

}

static bool
find_quadratic_xpt_near_target(gkyl_efit *up, double Rtar, double Ztar,
  double *Rout, double *Zout, double *psiout, double *metric_out)
{
  // Search a small stencil of P2 cells around the cubic X-point location and
  // run Newton seeded near the target point in reference coordinates.
  // This avoids selecting a different (spurious) critical point elsewhere.
  const double eps_ref = 1e-10;

  double point[2] = {Rtar, Ztar};
  int known[2] = {-1, -1};
  int base_idx[2] = {0, 0};
  gkyl_rect_grid_find_cell(&up->rzgrid, point, true, known, base_idx);

  bool found = false;
  double best_m = DBL_MAX;
  double best_R = 0.0, best_Z = 0.0, best_psi = DBL_MAX;

  // Stencil size: keep small but robust to boundary placement.
  for (int dR = -2; dR <= 2; ++dR) {
    for (int dZ = -2; dZ <= 2; ++dZ) {
      int idx[2] = { base_idx[0] + dR, base_idx[1] + dZ };
      if (!gkyl_range_contains_idx(&up->rzlocal, idx))
        continue;

      const double* psi = gkyl_array_cfetch(up->psizr, gkyl_range_idx(&up->rzlocal, idx));

      double xc[2];
      gkyl_rect_grid_cell_center(&up->rzgrid, idx, xc);

      double xinit[2];
      xinit[0] = 2.0*(Rtar - xc[0]) / up->rzgrid.dx[0];
      xinit[1] = 2.0*(Ztar - xc[1]) / up->rzgrid.dx[1];

      // If the target lies on/near a cell boundary, nudging the seed inward helps
      // avoid rejecting roots due to tiny out-of-range roundoff.
      xinit[0] = gkyl_clamp(xinit[0], -0.9, 0.9);
      xinit[1] = gkyl_clamp(xinit[1], -0.9, 0.9);

      double xsol[2];
      bool status = newton_raphson_seed(up, psi, xinit, xsol, false);
      if (!status)
        continue;

      double xr = xsol[0], zr = xsol[1];
      if (!(xr >= -1.0-eps_ref && xr <= 1.0+eps_ref && zr >= -1.0-eps_ref && zr <= 1.0+eps_ref))
        continue;

      xr = gkyl_clamp(xr, -1.0, 1.0);
      zr = gkyl_clamp(zr, -1.0, 1.0);
      double xrz[2] = {xr, zr};

      double R0 = up->rzgrid.dx[0]*xr/2.0 + xc[0];
      double Z0 = up->rzgrid.dx[1]*zr/2.0 + xc[1];
      double psi0 = up->rzbasis.eval_expand(xrz, psi);

      double m = gkyl_dist2(R0, Z0, Rtar, Ztar);
      if (m < best_m) {
        best_m = m;
        best_R = R0;
        best_Z = Z0;
        best_psi = psi0;
        found = true;
      }
    }
  }

  if (found) {
    *Rout = best_R;
    *Zout = best_Z;
    *psiout = best_psi;
    *metric_out = best_m;
  }
  return found;
}

int
find_xpts(gkyl_efit* up, double *Rxpt, double *Zxpt)
{
  // Use cubic X-point as a robust reference location.
  // This fixes cases where reflect=true (half-domain workflows) causes the quadratic
  // selector to lock onto a different in-bounds critical point.
  double Rcubic[2] = {0.0, 0.0};
  double Zcubic[2] = {0.0, 0.0};
  int num_cubic = find_xpts_cubic(up, Rcubic, Zcubic);
  bool have_cubic = num_cubic > 0;

  double Rtarget_top = 0.0, Ztarget_top = 0.0;
  double Rtarget_bot = 0.0, Ztarget_bot = 0.0;
  if (have_cubic) {
    double Rref = Rcubic[0];
    double Zref = Zcubic[0];
    double absZ = fabs(Zref);
    Rtarget_top = Rref;
    Ztarget_top = absZ;
    Rtarget_bot = Rref;
    Ztarget_bot = -absZ;
  }

  bool found_any = false;
  bool found_top = false, found_bot = false;
  double Rbest = 0.0, Zbest = 0.0, psibest = DBL_MAX;
  double metric_best = DBL_MAX;

  double Rtop = 0.0, Ztop = 0.0, psitop = DBL_MAX, metric_top = DBL_MAX;
  double Rbot = 0.0, Zbot = 0.0, psibot = DBL_MAX, metric_bot = DBL_MAX;

  // If we have a cubic reference, first try a targeted stencil search in the P2 field
  // around that location. This is much more reliable than a global scan when the P2
  // field contains multiple critical points.
  if (have_cubic) {
    if (up->reflect) {
      found_top = find_quadratic_xpt_near_target(up, Rtarget_top, Ztarget_top, &Rtop, &Ztop, &psitop, &metric_top);
      found_bot = find_quadratic_xpt_near_target(up, Rtarget_bot, Ztarget_bot, &Rbot, &Zbot, &psibot, &metric_bot);
      found_any = found_top || found_bot;
    }
    else {
      found_any = find_quadratic_xpt_near_target(up, Rtarget_top, Ztarget_top, &Rbest, &Zbest, &psibest, &metric_best);
    }
  }

  // Fallback: global scan (previous behavior), but scored against the cubic reference if available.
  if (!found_any) {

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->rzlocal);
  while (gkyl_range_iter_next(&iter)) {
    // Always scan the full domain. For reflect=true, the DG coefficients were already
    // reflected across Z=0 in efit.c, so both halves are available for selection.
    const double* psi = gkyl_array_cfetch(up->psizr, gkyl_range_idx(&up->rzlocal, iter.idx));
    double xsol[2];
    bool status = newton_raphson(up, psi, xsol, false);
    double x0 = xsol[0];
    double y0 = xsol[1];
    if (!(x0 >= -1 && x0 <= 1 && y0 >= -1 && y0 <= 1 && status))
      continue;

    double psi0 = up->rzbasis.eval_expand(xsol, psi);

    double xc[2];
    gkyl_rect_grid_cell_center(&up->rzgrid, iter.idx, xc);
    double R0 = up->rzgrid.dx[0]*x0/2.0 + xc[0];
    double Z0 = up->rzgrid.dx[1]*y0/2.0 + xc[1];

    double metric = DBL_MAX;
    if (have_cubic) {
      if (up->reflect) {
        if (Z0 >= 0)
          metric = gkyl_dist2(R0, Z0, Rtarget_top, Ztarget_top);
        else
          metric = gkyl_dist2(R0, Z0, Rtarget_bot, Ztarget_bot);
      }
      else {
        metric = gkyl_dist2(R0, Z0, Rtarget_top, Ztarget_top);
      }
    }
    else {
      metric = fabs(psi0 - up->sibry);
    }

    found_any = true;

    if (!up->reflect) {
      if (metric < metric_best) {
        metric_best = metric;
        Rbest = R0;
        Zbest = Z0;
        psibest = psi0;
      }
    }
    else {
      if (Z0 >= 0) {
        if (metric < metric_top) {
          metric_top = metric;
          Rtop = R0;
          Ztop = Z0;
          psitop = psi0;
          found_top = true;
        }
      }
      else {
        if (metric < metric_bot) {
          metric_bot = metric;
          Rbot = R0;
          Zbot = Z0;
          psibot = psi0;
          found_bot = true;
        }
      }
    }
  }

  } // end fallback global scan

  int num_xpts = 0;
  if (found_any) {
    if (up->reflect) {
      // Choose the closer of (top,bottom) to the cubic target, then mirror it.
      // This preserves symmetry while allowing selection from either half.
      bool have_choice = found_top || found_bot;
      if (have_choice) {
        double Rsel = 0.0, Zsel = 0.0, psisel = DBL_MAX;
        double msel = DBL_MAX;

        if (found_top) {
          Rsel = Rtop; Zsel = Ztop; psisel = psitop; msel = metric_top;
        }
        if (found_bot && (!found_top || metric_bot < msel)) {
          Rsel = Rbot; Zsel = Zbot; psisel = psibot; msel = metric_bot;
        }

        double absZ = fabs(Zsel);

        // If the quadratic-selected point is not very close to the cubic one,
        // prefer the cubic X-point for geometry construction. This is especially
        // important when reflect=true modifies the P2 field enough to introduce a
        // different in-bounds critical point.
        if (have_cubic) {
          double tol2 = gkyl_sq(0.25*up->rzgrid.dx[0]) + gkyl_sq(0.25*up->rzgrid.dx[1]);
          double d2 = gkyl_dist2(Rsel, absZ, Rtarget_top, Ztarget_top);
          if (d2 > tol2) {
            Rsel = Rtarget_top;
            absZ = Ztarget_top;
            psisel = up->psisep_cubic;
          }
        }

        num_xpts = 2;
        Rxpt[0] = Rsel;
        Rxpt[1] = Rsel;
        Zxpt[0] = -absZ;
        Zxpt[1] = absZ;
        up->psisep = psisel;
      }
    }
    else {
      if (have_cubic) {
        double tol2 = gkyl_sq(0.25*up->rzgrid.dx[0]) + gkyl_sq(0.25*up->rzgrid.dx[1]);
        double d2 = gkyl_dist2(Rbest, Zbest, Rtarget_top, Ztarget_top);
        if (d2 > tol2) {
          Rbest = Rtarget_top;
          Zbest = Ztarget_top;
          psibest = up->psisep_cubic;
        }
      }
      num_xpts = 1;
      Rxpt[0] = Rbest;
      Zxpt[0] = Zbest;
      up->psisep = psibest;
    }
  }
  return num_xpts;
}

int
find_xpts_cubic(gkyl_efit* up, double *Rxpt, double *Zxpt)
{
  bool found_xpt = false;
  double Rsep, Zsep;
  double psisep = DBL_MAX;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->rzlocal_cubic);
  while (gkyl_range_iter_next(&iter)) {
    if ((iter.idx[1] < gkyl_range_shape(&up->rzlocal_cubic,1)/2 + 1) || (!up->reflect)) {
      const double* psi = gkyl_array_cfetch(up->psizr_cubic, gkyl_range_idx(&up->rzlocal_cubic, iter.idx));
      double xsol[2];
      bool status = newton_raphson(up, psi, xsol, true);
      double x0 = xsol[0];
      double y0 = xsol[1];
      double psi0 = up->rzbasis_cubic.eval_expand(xsol, psi);
      if (x0 >= -1 && x0 <= 1 && y0 >= -1 && y0 <= 1 && status) {
        found_xpt = true;
        double xc[2];
        gkyl_rect_grid_cell_center(&up->rzgrid_cubic, iter.idx, xc);
        double R0 = up->rzgrid_cubic.dx[0]*x0/2.0 + xc[0];
        double Z0 = up->rzgrid_cubic.dx[1]*y0/2.0 + xc[1];
        if (fabs(psi0 - up->sibry) <= fabs(psisep - up->sibry)) {
          Rsep = R0;
          Zsep = Z0;
          psisep = psi0;
        }
      }
    }
  }

  int num_xpts = 0;
  if (found_xpt) {
    if (up->reflect) {
      num_xpts = 2;
      Rxpt[0] = Rsep;
      Rxpt[1] = Rsep;
      Zxpt[0] = Zsep;
      Zxpt[1] = -Zsep;
      up->psisep_cubic = psisep;
    }
    else {
      num_xpts = 1;
      Rxpt[0] = Rsep;
      Zxpt[0] = Zsep;
      up->psisep_cubic = psisep;
    }
  }
  return num_xpts;
}

