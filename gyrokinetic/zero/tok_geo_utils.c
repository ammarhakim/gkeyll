#include <gkyl_tok_geo_priv.h>
#include <assert.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

// Helper functions for finding turning points when necessary

void find_lower_turning_point(struct gkyl_tok_geo *geo, double psi_curr,
  double zup, double *zmin, double tolerance);

static double
tok_eval_psi_rz(const struct gkyl_tok_geo *geo, double R, double Z)
{
  if (geo->use_cubics) {
    double xn[2] = { R, Z }, psi;
    geo->efit->evf->eval_cubic(0.0, xn, &psi, geo->efit->evf->ctx);
    return psi;
  }

  int rzidx[2];
  rzidx[0] = fmin(geo->rzlocal.lower[0]
      + (int) floor((R-geo->rzgrid.lower[0])/geo->rzgrid.dx[0]),
    geo->rzlocal.upper[0]);
  rzidx[1] = fmin(geo->rzlocal.lower[1]
      + (int) floor((Z-geo->rzgrid.lower[1])/geo->rzgrid.dx[1]),
    geo->rzlocal.upper[1]);
  rzidx[0] = fmax(rzidx[0], geo->rzlocal.lower[0]);
  rzidx[1] = fmax(rzidx[1], geo->rzlocal.lower[1]);

  long loc = gkyl_range_idx(&geo->rzlocal, rzidx);
  const double *coeffs = gkyl_array_cfetch(geo->psiRZ, loc);
  double xc[2], xy[2];
  gkyl_rect_grid_cell_center(&geo->rzgrid, rzidx, xc);
  xy[0] = (R-xc[0])/(0.5*geo->rzgrid.dx[0]);
  xy[1] = (Z-xc[1])/(0.5*geo->rzgrid.dx[1]);
  return geo->rzbasis.eval_expand(xy, coeffs);
}

// Find the unique intersection of psi=psi0 with a shaped plate.  The search
// uses only the fixed target surface, so the cached X-point ray is identical
// in corner, interior, and surface geometry passes.
static bool
tok_plate_intersection(const struct gkyl_tok_geo *geo, plate_func plate,
  double psi0, double *rplate, double *zplate)
{
  if (!plate)
    return false;
  const int nsamp = 512;
  const double flux_tol = 1e-10*fmax(1.0, fabs(psi0));
  int nroot = 0;
  double last_root = -DBL_MAX, root_s = 0.0;
  double rz0[2], rz1[2];
  plate(0.0, rz0);
  double f0 = tok_eval_psi_rz(geo, rz0[0], rz0[1])-psi0;
  if (isfinite(f0) && fabs(f0) <= flux_tol) {
    root_s = 0.0; last_root = 0.0; nroot = 1;
  }
  for (int i=1; i<=nsamp; ++i) {
    double s1 = i/(double) nsamp;
    plate(s1, rz1);
    double f1 = tok_eval_psi_rz(geo, rz1[0], rz1[1])-psi0;
    bool endpoint = isfinite(f1) && fabs(f1) <= flux_tol;
    bool bracket = isfinite(f0) && isfinite(f1) && f0*f1 < 0.0;
    if (endpoint || bracket) {
      double sr = s1;
      if (!endpoint) {
        double lo = (i-1)/(double) nsamp, hi = s1, flo = f0;
        for (int k=0; k<64; ++k) {
          double mid = 0.5*(lo+hi), rzm[2];
          plate(mid, rzm);
          double fm = tok_eval_psi_rz(geo, rzm[0], rzm[1])-psi0;
          if (!isfinite(fm))
            return false;
          if (flo*fm <= 0.0)
            hi = mid;
          else { lo = mid; flo = fm; }
        }
        sr = 0.5*(lo+hi);
      }
      if (fabs(sr-last_root) > 2e-10) {
        root_s = sr; last_root = sr; ++nroot;
      }
    }
    rz0[0] = rz1[0]; rz0[1] = rz1[1]; f0 = f1;
  }
  if (nroot != 1)
    return false;
  double rz[2];
  plate(root_s, rz);
  *rplate = rz[0]; *zplate = rz[1];
  return isfinite(*rplate) && isfinite(*zplate);
}

enum tok_xpt_sector {
  TOK_XPT_CORE,
  TOK_XPT_PF,
  TOK_XPT_SOL_OUT,
  TOK_XPT_SOL_IN,
  TOK_XPT_UNSUPPORTED,
};

static enum tok_xpt_sector
tok_xpt_sector_from_ftype(enum gkyl_tok_geo_type ftype)
{
  switch (ftype) {
    case GKYL_GEOMETRY_TOKAMAK_CORE_R:
    case GKYL_GEOMETRY_TOKAMAK_CORE_L:
      return TOK_XPT_CORE;
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_R:
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_L:
      return TOK_XPT_PF;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID:
      return TOK_XPT_SOL_OUT;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO:
      return TOK_XPT_SOL_IN;
    default:
      return TOK_XPT_UNSUPPORTED;
  }
}

static bool
tok_xpt_ray_enabled(const struct gkyl_tok_geo_grid_inp *inp)
{
  enum tok_xpt_sector sector = tok_xpt_sector_from_ftype(inp->ftype);
  if (sector == TOK_XPT_UNSUPPORTED)
    return false;
  if (inp->straight_xpt_ray)
    return true;
  return inp->straight_core_xpt_ray && sector == TOK_XPT_CORE;
}

static bool
tok_nearest_root_at_z(const struct gkyl_tok_geo_grid_inp *inp,
  const struct gkyl_tok_geo *geo,
  enum tok_xpt_sector sector, double psi0, double z, double rx, double zx,
  double *rbest, double *d2best)
{
  double R[4] = { 0 }, dRdZ[4] = { 0 }, dR[4] = { 0 }, dZ[4] = { 0 };
  int nr = gkyl_tok_geo_R_psiZ(geo, psi0, z, 4, R, dRdZ, dR, dZ);
  bool found = false;
  *d2best = DBL_MAX;

  // An R threshold around the X point lets the inboard and outboard tests
  // overlap whenever the two roots approach the saddle.  Select the root
  // belonging to the requested SOL branch by its established branch seed
  // instead.  Core and PF searches remain global because their nearest point
  // can lie on either side of a turning point.
  int nlo = 0, nhi = nr;
  int selected = -1;
  if (sector == TOK_XPT_SOL_OUT || sector == TOK_XPT_SOL_IN) {
    double reference = sector == TOK_XPT_SOL_OUT ? inp->rright : inp->rleft;
    double best_reference_distance = DBL_MAX;
    for (int n=0; n<nr; ++n) {
      double distance = fabs(R[n]-reference);
      if (distance < best_reference_distance) {
        best_reference_distance = distance;
        selected = n;
      }
    }
    nlo = selected;
    nhi = selected+1;
  }
  for (int n=0; n<nr; ++n) {
    if (n < nlo || n >= nhi)
      continue;
    double d2 = SQ(R[n]-rx)+SQ(z-zx);
    if (d2 < *d2best) {
      *d2best = d2;
      *rbest = R[n];
      found = true;
    }
  }
  return found;
}

static bool
tok_nearest_point_on_surface(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx, enum tok_xpt_sector sector, double psi0,
  double rx, double zx, double *rnear, double *znear)
{
  const struct gkyl_tok_geo *geo = arc_ctx->geo;
  double zlo, zhi;
  if (sector == TOK_XPT_PF) {
    // These input bounds are fixed for the block.  Do not use the current
    // surface's plate intersections here: corner, interior, and surface
    // calculations begin at different radial quadrature points, but must all
    // initialize exactly the same ray on xpt_ray_psi0.
    zlo = fmax(geo->rzgrid.lower[1],
      fmin(inp->zmin_left, inp->zmin_right));
    zhi = zx;
  }
  else if (sector == TOK_XPT_SOL_OUT || sector == TOK_XPT_SOL_IN) {
    // A lower-SOL half surface runs from its physical divertor plate to the
    // magnetic-axis midplane.  The nearest point to the X point is allowed to
    // lie below the X point; restricting this search to X->midplane selected a
    // visibly non-nearest ray in shot 203730.
    double rplate = 0.0, zplate = 0.0;
    if (!tok_plate_intersection(geo, geo->plate_func_lower, psi0,
        &rplate, &zplate))
      return false;
    zlo = fmin(zplate, geo->zmaxis);
    zhi = fmax(zplate, geo->zmaxis);
  }
  else {
    zlo = fmin(zx, geo->zmaxis);
    zhi = fmax(zx, geo->zmaxis);
  }
  if (!(zhi > zlo))
    return false;

  // This search is intentionally global in Z. It is run only once per block,
  // and avoids making the selected ray depend on a branch-local Newton guess.
  const int nsamp = 4*geo->rzgrid.cells[1]+1;
  double best_d2 = DBL_MAX, best_r = 0.0, best_z = 0.0;
  int best_i = -1;
  for (int i=0; i<nsamp; ++i) {
    double z = zlo+(zhi-zlo)*i/(nsamp-1.0), r, d2;
    if (tok_nearest_root_at_z(inp, geo, sector, psi0, z, rx, zx, &r, &d2)
        && d2 < best_d2) {
      best_d2 = d2; best_r = r; best_z = z; best_i = i;
    }
  }
  if (best_i < 0)
    return false;

  // Refine the sampled minimum without differentiating through DG-cell
  // interfaces or assuming that only one R root exists at a given Z.
  double dzs = (zhi-zlo)/(nsamp-1.0);
  double a = fmax(zlo, best_z-dzs), b = fmin(zhi, best_z+dzs);
  const double gr = 0.6180339887498948482;
  for (int iter=0; iter<64; ++iter) {
    double ztrial[2] = { b-gr*(b-a), a+gr*(b-a) };
    double dtrial[2] = { DBL_MAX, DBL_MAX }, rtrial[2] = { best_r, best_r };
    for (int k=0; k<2; ++k)
      tok_nearest_root_at_z(inp, geo, sector, psi0, ztrial[k], rx, zx,
        &rtrial[k], &dtrial[k]);
    if (dtrial[0] < best_d2) {
      best_d2 = dtrial[0]; best_r = rtrial[0]; best_z = ztrial[0];
    }
    if (dtrial[1] < best_d2) {
      best_d2 = dtrial[1]; best_r = rtrial[1]; best_z = ztrial[1];
    }
    if (dtrial[0] <= dtrial[1])
      b = ztrial[1];
    else
      a = ztrial[0];
  }

  *rnear = best_r;
  *znear = best_z;
  return isfinite(best_r) && isfinite(best_z);
}

bool
tok_xpt_classify_branch_at_point(const struct gkyl_tok_geo_grid_inp *inp,
  const struct arc_length_ctx *arc_ctx, double Rpoint, double Zpoint,
  double psi, bool *resolved, bool *on_right)
{
  double R[4] = { 0 }, dRdZ[4] = { 0 }, dR[4] = { 0 }, dZ[4] = { 0 };
  int nr = gkyl_tok_geo_R_psiZ(arc_ctx->geo, psi, Zpoint,
    4, R, dRdZ, dR, dZ);
  *resolved = false;
  if (nr < 2)
    return true;

  double rr = choose_closest(inp->rright, R, R, nr);
  double rl = choose_closest(inp->rleft, R, R, nr);
  double scale = fmax(1.0, fmax(fabs(rr), fabs(rl)));
  if (fabs(rr-rl) <= 1e-8*scale)
    return true;

  double error_right = fabs(Rpoint-rr), error_left = fabs(Rpoint-rl);
  *on_right = error_right <= error_left;
  *resolved = true;
  return true;
}

static bool
tok_xpt_initialize_branch(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx, enum tok_xpt_sector sector,
  double rx, double zx)
{
  if (sector == TOK_XPT_SOL_OUT || sector == TOK_XPT_SOL_IN) {
    arc_ctx->xpt_ray_on_right = sector == TOK_XPT_SOL_OUT;
    arc_ctx->xpt_ray_branch_valid = true;
    return true;
  }

  // A nearest core/PF endpoint is commonly a turning point, where the two R
  // roots coalesce and cannot label the branch.  Probe the already-fixed ray
  // away from both endpoints, using psi evaluated at the probe itself.  This
  // makes branch identity independent of which corner/Gauss surface happened
  // to initialize this arc-length context first.
  const double probes[] = { 0.5, 0.25, 0.75, 0.125, 0.875 };
  bool found = false, selected_right = false;
  for (unsigned i=0; i<sizeof(probes)/sizeof(probes[0]); ++i) {
    double s = probes[i];
    double R = rx+s*(arc_ctx->xpt_ray_r0-rx);
    double Z = zx+s*(arc_ctx->xpt_ray_z0-zx);
    double psi = tok_eval_psi_rz(arc_ctx->geo, R, Z);
    bool resolved = false, on_right = false;
    if (!tok_xpt_classify_branch_at_point(inp, arc_ctx, R, Z, psi,
        &resolved, &on_right))
      return false;
    if (!resolved)
      continue;
    if (found && on_right != selected_right) {
      fprintf(stderr,
        "TOK_XPT_RAY_DIAG reason=ray_branch_crossing ftype=%d probe_s=%.17g\n",
        inp->ftype, s);
      return false;
    }
    found = true;
    selected_right = on_right;
  }
  if (!found) {
    fprintf(stderr,
      "TOK_XPT_RAY_DIAG reason=unresolved_ray_branch ftype=%d endpoint=(%.17g,%.17g)\n",
      inp->ftype, arc_ctx->xpt_ray_r0, arc_ctx->xpt_ray_z0);
    return false;
  }
  arc_ctx->xpt_ray_on_right = selected_right;
  arc_ctx->xpt_ray_branch_valid = true;
  return true;
}

static bool
tok_xpt_ray_anchor(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx, double psi_curr)
{
  const struct gkyl_tok_geo *geo = arc_ctx->geo;
  enum tok_xpt_sector sector = tok_xpt_sector_from_ftype(inp->ftype);
  arc_ctx->xpt_anchor_valid = false;
  double rx = geo->use_cubics ? geo->efit->Rxpt_cubic[0] : geo->efit->Rxpt[0];
  double zx = geo->use_cubics ? geo->efit->Zxpt_cubic[0] : geo->efit->Zxpt[0];

  if (!arc_ctx->xpt_ray_initialized) {
    if (!tok_nearest_point_on_surface(inp, arc_ctx, sector,
        arc_ctx->xpt_ray_psi0, rx, zx,
        &arc_ctx->xpt_ray_r0, &arc_ctx->xpt_ray_z0)) {
      fprintf(stderr, "TOK_XPT_RAY_DIAG reason=no_far_surface_point ftype=%d target_psi=%.17g\n",
        inp->ftype, arc_ctx->xpt_ray_psi0);
      return false;
    }

    double psi_endpoint = tok_eval_psi_rz(geo, arc_ctx->xpt_ray_r0,
      arc_ctx->xpt_ray_z0);
    double flux_scale = fmax(1.0,
      fmax(fabs(arc_ctx->xpt_ray_psi0), fabs(geo->psisep)));
    if (!isfinite(psi_endpoint) ||
        fabs(psi_endpoint-arc_ctx->xpt_ray_psi0) > 1e-9*flux_scale) {
      fprintf(stderr,
        "TOK_XPT_RAY_DIAG reason=far_surface_flux_residual ftype=%d residual=%.17g requested_psi=%.17g evaluated_psi=%.17g endpoint=(%.17g,%.17g)\n",
        inp->ftype, psi_endpoint-arc_ctx->xpt_ray_psi0,
        arc_ctx->xpt_ray_psi0, psi_endpoint, arc_ctx->xpt_ray_r0,
        arc_ctx->xpt_ray_z0);
      return false;
    }
    double delta = arc_ctx->xpt_ray_psi0-geo->psisep;
    if (!isfinite(delta) || fabs(delta) <= 256.0*DBL_EPSILON*flux_scale) {
      fprintf(stderr, "TOK_XPT_RAY_DIAG reason=degenerate_flux_span ftype=%d psi_x=%.17g psi_far=%.17g\n",
        inp->ftype, geo->psisep, arc_ctx->xpt_ray_psi0);
      return false;
    }

    // The quadratic DG representation can have a small local reversal at a
    // cell boundary even when the underlying equilibrium is smooth.  A global
    // monotonicity test therefore rejects otherwise unambiguous rays.  Check
    // finiteness here; uniqueness is tested below for every actual geometry
    // and finite-difference target before an anchor is accepted.
    for (int i=1; i<=256; ++i) {
      double s = i/256.0;
      double R = rx+s*(arc_ctx->xpt_ray_r0-rx);
      double Z = zx+s*(arc_ctx->xpt_ray_z0-zx);
      double q = (tok_eval_psi_rz(geo, R, Z)-geo->psisep)/delta;
      if (!isfinite(q)) {
        fprintf(stderr,
          "TOK_XPT_RAY_DIAG reason=nonfinite_initial_ray ftype=%d sample=%d endpoint=(%.17g,%.17g)\n",
          inp->ftype, i, arc_ctx->xpt_ray_r0, arc_ctx->xpt_ray_z0);
        return false;
      }
    }
    if (!tok_xpt_initialize_branch(inp, arc_ctx, sector, rx, zx))
      return false;
    arc_ctx->xpt_ray_initialized = true;
  }

  double psi1 = arc_ctx->xpt_ray_psi0;
  double delta = psi1-geo->psisep;
  double scale = fmax(1.0, fmax(fabs(psi_curr), fabs(geo->psisep)));
  double tol = 128.0*DBL_EPSILON*scale;
  if (fabs(psi_curr-geo->psisep) <= tol) {
    arc_ctx->xpt_anchor_r = rx;
    arc_ctx->xpt_anchor_z = zx;
    arc_ctx->xpt_anchor_valid = true;
    return true;
  }
  if (fabs(psi_curr-psi1) <= tol) {
    arc_ctx->xpt_anchor_r = arc_ctx->xpt_ray_r0;
    arc_ctx->xpt_anchor_z = arc_ctx->xpt_ray_z0;
    arc_ctx->xpt_anchor_valid = true;
    return true;
  }

  double dr = arc_ctx->xpt_ray_r0-rx, dz = arc_ctx->xpt_ray_z0-zx;
  double smax = DBL_MAX;
  if (dr > 0.0)
    smax = fmin(smax, (geo->rzgrid.upper[0]-rx)/dr);
  else if (dr < 0.0)
    smax = fmin(smax, (geo->rzgrid.lower[0]-rx)/dr);
  if (dz > 0.0)
    smax = fmin(smax, (geo->rzgrid.upper[1]-zx)/dz);
  else if (dz < 0.0)
    smax = fmin(smax, (geo->rzgrid.lower[1]-zx)/dz);
  smax *= 1.0-64.0*DBL_EPSILON;
  if (!isfinite(smax) || smax <= 0.0) {
    fprintf(stderr, "TOK_XPT_RAY_DIAG reason=invalid_ray_extent ftype=%d smax=%.17g\n",
      inp->ftype, smax);
    return false;
  }

  double qtarget = (psi_curr-geo->psisep)/delta;
  if (qtarget < -1e-8 || qtarget > 1.0+1e-8) {
    fprintf(stderr,
      "TOK_XPT_RAY_DIAG reason=target_outside_ray_span ftype=%d qtarget=%.17g psi=%.17g ray_psi=(%.17g,%.17g)\n",
      inp->ftype, qtarget, psi_curr, geo->psisep, psi1);
    return false;
  }
  double search_max = fmin(1.0, smax);
  double qprev = (tok_eval_psi_rz(geo, rx, zx)-geo->psisep)/delta;
  double qpeak = qprev;
  double sprev = 0.0, slo = -1.0, sup = -1.0;
  int upward_crossings = 0, downward_crossings = 0;
  double crossing_tol = 1e-13*fmax(1.0, fabs(qtarget));
  const int nsamp = 512;
  for (int i=1; i<=nsamp; ++i) {
    double s = search_max*i/(double) nsamp;
    double R = rx+s*dr, Z = zx+s*dz;
    double q = (tok_eval_psi_rz(geo, R, Z)-geo->psisep)/delta;
    if (!isfinite(q)) {
      fprintf(stderr,
        "TOK_XPT_RAY_DIAG reason=nonfinite_anchor_search ftype=%d sample=%d qtarget=%.17g\n",
        inp->ftype, i, qtarget);
      return false;
    }
    bool upward = qprev < qtarget-crossing_tol &&
      q >= qtarget-crossing_tol;
    bool downward = qprev > qtarget+crossing_tol &&
      q <= qtarget+crossing_tol;
    if (upward) {
      upward_crossings++;
      if (slo < 0.0) {
        slo = sprev;
        sup = s;
      }
    }
    if (downward)
      downward_crossings++;
    qprev = q;
    qpeak = fmax(qpeak, q);
    sprev = s;
  }
  if (slo < 0.0) {
    fprintf(stderr,
      "TOK_XPT_RAY_DIAG reason=no_anchor_bracket ftype=%d qtarget=%.17g qpeak=%.17g qend=%.17g smax=%.17g\n",
      inp->ftype, qtarget, qpeak, qprev, search_max);
    return false;
  }
  if (upward_crossings != 1) {
    fprintf(stderr,
      "TOK_XPT_RAY_DIAG reason=ambiguous_anchor_crossing ftype=%d qtarget=%.17g upward_crossings=%d downward_crossings=%d first_bracket=(%.17g,%.17g)\n",
      inp->ftype, qtarget, upward_crossings, downward_crossings, slo, sup);
    return false;
  }

  for (int n=0; n<80; ++n) {
    double smid = 0.5*(slo+sup);
    double R = rx+smid*dr, Z = zx+smid*dz;
    double qmid = (tok_eval_psi_rz(geo, R, Z)-geo->psisep)/delta;
    if (!isfinite(qmid)) {
      fprintf(stderr,
        "TOK_XPT_RAY_DIAG reason=nonfinite_anchor_refinement ftype=%d qtarget=%.17g s=%.17g\n",
        inp->ftype, qtarget, smid);
      return false;
    }
    if (qmid-qtarget >= 0.0)
      sup = smid;
    else
      slo = smid;
  }
  double s = 0.5*(slo+sup);
  arc_ctx->xpt_anchor_r = rx+s*dr;
  arc_ctx->xpt_anchor_z = zx+s*dz;
  double residual = tok_eval_psi_rz(geo, arc_ctx->xpt_anchor_r,
    arc_ctx->xpt_anchor_z)-psi_curr;
  if (!isfinite(residual) || fabs(residual) > 1e-9*scale) {
    fprintf(stderr, "TOK_XPT_RAY_DIAG reason=anchor_flux_residual ftype=%d residual=%.17g anchor=(%.17g,%.17g)\n",
      inp->ftype, residual, arc_ctx->xpt_anchor_r,
      arc_ctx->xpt_anchor_z);
    return false;
  }
  arc_ctx->xpt_anchor_valid = true;
  return true;
}

static bool
tok_xpt_anchor_on_right(const struct gkyl_tok_geo_grid_inp *inp,
  const struct arc_length_ctx *arc_ctx)
{
  if (!arc_ctx->xpt_ray_branch_valid) {
    fprintf(stderr,
      "TOK_XPT_RAY_DIAG reason=missing_cached_branch ftype=%d\n",
      inp->ftype);
    abort();
  }
  bool resolved = false, on_right = false;
  if (!tok_xpt_classify_branch_at_point(inp, arc_ctx,
      arc_ctx->xpt_anchor_r, arc_ctx->xpt_anchor_z, arc_ctx->psi,
      &resolved, &on_right))
    abort();
  if (resolved && on_right != arc_ctx->xpt_ray_on_right) {
    fprintf(stderr,
      "TOK_XPT_RAY_DIAG reason=ray_branch_crossing ftype=%d psi=%.17g anchor=(%.17g,%.17g) cached_right=%d resolved_right=%d\n",
      inp->ftype, arc_ctx->psi, arc_ctx->xpt_anchor_r,
      arc_ctx->xpt_anchor_z, arc_ctx->xpt_ray_on_right, on_right);
    abort();
  }
  return arc_ctx->xpt_ray_on_right;
}

static void
tok_xpt_map_interval(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx, double lo, double hi)
{
  double dtheta = inp->cgrid.upper[2]-inp->cgrid.lower[2];
  if (!isfinite(lo) || !isfinite(hi) || !(hi > lo) || !(dtheta > 0.0)) {
    fprintf(stderr,
      "TOK_XPT_RAY invalid arc interval ftype=%d psi=%.17g lo=%.17g hi=%.17g dtheta=%.17g\n",
      inp->ftype, arc_ctx->psi, lo, hi, dtheta);
    abort();
  }
  arc_ctx->xpt_map_arc_lo = lo;
  arc_ctx->xpt_map_arc_hi = hi;
  arc_ctx->xpt_map_darc_dtheta = (hi-lo)/dtheta;
  arc_ctx->xpt_map_valid = true;
}

static void
tok_configure_xpt_map(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx)
{
  arc_ctx->xpt_map_valid = false;
  arc_ctx->xpt_anchor_valid = false;
  if (!tok_xpt_ray_enabled(inp))
    return;
  if (!inp->half_domain) {
    fprintf(stderr,
      "TOK_XPT_RAY ftype=%d currently requires half_domain=true\n", inp->ftype);
    abort();
  }
  if (!tok_xpt_ray_anchor(inp, arc_ctx, arc_ctx->psi)) {
    fprintf(stderr,
      "TOK_XPT_RAY failed ftype=%d psi=%.17g target_psi=%.17g\n",
      inp->ftype, arc_ctx->psi, arc_ctx->xpt_ray_psi0);
    abort();
  }

  enum tok_xpt_sector sector = tok_xpt_sector_from_ftype(inp->ftype);
  double cut = 0.0, lo = 0.0, hi = 0.0;
  bool at_sep = fabs(arc_ctx->psi-arc_ctx->geo->psisep)
    <= 128.0*DBL_EPSILON*fmax(1.0, fabs(arc_ctx->geo->psisep));

  if (sector == TOK_XPT_CORE) {
    if (at_sep)
      cut = 0.0;
    else if (tok_xpt_anchor_on_right(inp, arc_ctx))
      cut = integrate_psi_contour_memo(arc_ctx->geo, arc_ctx->psi,
        arc_ctx->zmin, arc_ctx->xpt_anchor_z, inp->rright,
        true, false, arc_ctx->arc_memo_right);
    else
      cut = arc_ctx->arcL_right
        + integrate_psi_contour_memo(arc_ctx->geo, arc_ctx->psi,
          arc_ctx->xpt_anchor_z, arc_ctx->zmax, inp->rleft,
          true, false, arc_ctx->arc_memo_left);
    double out_mid = integrate_psi_contour_memo(arc_ctx->geo, arc_ctx->psi,
      arc_ctx->zmin, arc_ctx->geo->zmaxis, inp->rright,
      true, false, arc_ctx->arc_memo_right);
    double in_mid = arc_ctx->arcL_right
      + integrate_psi_contour_memo(arc_ctx->geo, arc_ctx->psi,
        arc_ctx->geo->zmaxis, arc_ctx->zmax, inp->rleft,
        true, false, arc_ctx->arc_memo_left);
    if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R) {
      lo = cut; hi = out_mid;
    }
    else {
      lo = in_mid; hi = cut;
    }
    while (hi <= lo)
      hi += arc_ctx->arcL_tot;
  }
  else if (sector == TOK_XPT_PF) {
    if (at_sep)
      cut = arc_ctx->arcL_right;
    else if (tok_xpt_anchor_on_right(inp, arc_ctx))
      cut = integrate_psi_contour_memo(arc_ctx->geo, arc_ctx->psi,
        arc_ctx->zmin_right, arc_ctx->xpt_anchor_z, inp->rright,
        true, false, arc_ctx->arc_memo_right);
    else
      cut = arc_ctx->arcL_right
        + integrate_psi_contour_memo(arc_ctx->geo, arc_ctx->psi,
          arc_ctx->xpt_anchor_z, arc_ctx->zmax, inp->rleft,
          true, false, arc_ctx->arc_memo_left);
    if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_R) {
      lo = 0.0; hi = cut;
    }
    else {
      lo = cut; hi = arc_ctx->arcL_tot;
    }
  }
  else if (sector == TOK_XPT_SOL_OUT) {
    cut = at_sep
      ? integrate_psi_contour_memo(arc_ctx->geo, arc_ctx->psi,
          arc_ctx->zmin, arc_ctx->xpt_anchor_z, inp->rright,
          true, false, arc_ctx->arc_memo)
      : integrate_psi_contour_memo(arc_ctx->geo, arc_ctx->psi,
          arc_ctx->zmin, arc_ctx->xpt_anchor_z, inp->rright,
          true, false, arc_ctx->arc_memo);
    double mid = integrate_psi_contour_memo(arc_ctx->geo, arc_ctx->psi,
      arc_ctx->zmin, arc_ctx->geo->zmaxis, inp->rright,
      true, false, arc_ctx->arc_memo);
    if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO) {
      lo = 0.0; hi = cut;
    }
    else {
      lo = cut; hi = mid;
    }
  }
  else if (sector == TOK_XPT_SOL_IN) {
    cut = integrate_psi_contour_memo(arc_ctx->geo, arc_ctx->psi,
      arc_ctx->xpt_anchor_z, arc_ctx->zmax, inp->rleft,
      true, false, arc_ctx->arc_memo);
    double mid = integrate_psi_contour_memo(arc_ctx->geo, arc_ctx->psi,
      arc_ctx->geo->zmaxis, arc_ctx->zmax, inp->rleft,
      true, false, arc_ctx->arc_memo);
    if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID) {
      lo = mid; hi = cut;
    }
    else {
      lo = cut; hi = arc_ctx->arcL_tot;
    }
  }
  tok_xpt_map_interval(inp, arc_ctx, lo, hi);
}

void
tok_prepare_ordered_map(struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx, double psi_curr)
{
  arc_ctx->psi = psi_curr;
  arc_ctx->xpt_map_valid = false;
  arc_ctx->xpt_anchor_valid = false;
  // Preserve the old coordinate algebra in callers while the ordered map
  // consumes theta directly; no physical contour length is needed here.
  arc_ctx->arcL_tot = 2.0*M_PI;
  // Full-domain single- and double-null blocks use topology-aware boundary
  // traces built lazily in tok_geo.c.  They can have two distinct X-point
  // rays (or a closed core seam), so the scalar lower-X-point anchor used by
  // the original half-domain path is intentionally bypassed here.
  if (inp->straight_xpt_ray && !inp->half_domain)
    return;
  if (!tok_xpt_ray_enabled(inp) ||
      !tok_xpt_ray_anchor(inp, arc_ctx, psi_curr)) {
    fprintf(stderr,
      "TOK_ORDERED_MAP failed to initialize ray ftype=%d psi=%.17g target_psi=%.17g\n",
      inp->ftype, psi_curr, arc_ctx->xpt_ray_psi0);
    abort();
  }
}


// This function will set zmax to be the upper turning point location
void find_upper_turning_point(struct gkyl_tok_geo *geo, double psi_curr, double zlo, double *zmax, double tolerance)
{
    double tol = tolerance ? tolerance : 1e-12 ;
    //Find the turning points
    double zlo_last;
    double zup=*zmax;
    zlo_last = zlo;
    double R[4], dRdZ[4];
    double dR[4], dZ[4];
    double Rup[4], dRdZup[4];
    double dRup[4], dZup[4];
    while(true){
      int nlo = R_psiZ(geo, psi_curr, zlo, 4, R, dRdZ, dR, dZ);
      int nup = R_psiZ(geo, psi_curr, zup, 4, Rup, dRdZup, dRup, dZup);
      //printf("nlo, nup = %d %d; zlo, zup = %g %g\n", nlo, nup, zlo, zup);
      if (nup > 0) { // This is for the PF_LO regions. Does not seem to break core_L or core_R
                  // However I need to think thos through more. I think it is ok only when xpt
                  // is known quite precisely.
        *zmax = zup;
        break;
      }
      if (nlo>=1){
        if(fabs(zlo-zup)<tol){
          *zmax = zlo;
          break;
        }
        zlo_last = zlo;
        zlo = (zlo+zup)/2;
      }
      if(nlo==0){
        zup = zlo;
        zlo = zlo_last;
      }
    }
}

// This function will set zmin to be the upper turning point location
void find_lower_turning_point(struct gkyl_tok_geo *geo, double psi_curr, double zup, double *zmin, double tolerance)
{
    double tol = tolerance ? tolerance : 1e-12 ;
    int nup = 0;
    double zlo=*zmin;
    double zup_last = zup;
    double R[4], dRdZ[4];
    double dR[4], dZ[4];
    double Rlo[4], dRdZlo[4];
    double dRlo[4], dZlo[4];
    while(true){
      int nup = R_psiZ(geo, psi_curr, zup, 4, R, dRdZ, dR, dZ);
      int nlo = R_psiZ(geo, psi_curr, zlo, 4, Rlo, dRdZlo, dRlo, dZlo);
      //printf("psi = %g;lo, nup = %d %d; zlo, zup = %g %g\n", psi_curr, nlo, nup, zlo, zup);
      if (nlo > 0) {
        *zmin = zlo;
        break;
      }
      if(nup>=1){
        if(fabs(zlo-zup)<tol){
          *zmin = zup;
          break;
        }
        zup_last = zup;
        zup = (zlo+zup)/2;
      }
      if(nup==0){
        zlo = zup;
        zup = zup_last;
      }
    }

}

// This function will set zmin to be the upper lower point location
void find_lower_turning_point_pf_up(struct gkyl_tok_geo *geo, double psi_curr, double zup, double *zmin)
{
    int nup = 0;
    double zlo=*zmin;
    double zup_last = zup;
    double R[4], dRdZ[4];
    double dR[4], dZ[4];
    double Rlo[4], dRdZlo[4];
    double dRlo[4], dZlo[4];
    while(true){
      int nup = R_psiZ(geo, psi_curr, zup, 4, R, dRdZ, dR, dZ);
      int nlo = R_psiZ(geo, psi_curr, zlo, 4, Rlo, dRlo, dRlo, dZlo);
      //if(nlo==1){
      //  if (Rlo[0] < geo->rleft)
      //    nlo=0;
      //}
      if (nlo > 0) {
        *zmin = zlo;
        break;
      }
      if(nup>=2){
        if(fabs(zlo-zup)<1e-12){
          *zmin = zup;
          break;
        }
        zup_last = zup;
        zup = (zlo+zup)/2;
      }
      if(nup==1){
        zlo = zup;
        zup = zup_last;
      }
      if(nup==0){
        zlo = zup;
        zup = zup_last;
      }
    }
}

// This function will set zmax to be the upper turning point location
void find_upper_turning_point_pf_lo(struct gkyl_tok_geo *geo, double psi_curr, double zlo, double *zmax)
{
    //Find the turning points
    double zlo_last;
    double zup=*zmax;
    zlo_last = zlo;
    double R[4], dRdZ[4];
    double dR[4], dZ[4];
    double Rup[4], dRdZup[4];
    double dRup[4], dZup[4];
    while(true){
      int nlo = R_psiZ(geo, psi_curr, zlo, 4, R, dRdZ, dR, dZ);
      int nup = R_psiZ(geo, psi_curr, zup, 4, Rup, dRdZup, dRup, dZup);
      if(nup==1){
        if (Rup[0] < geo->rleft)
          nup=0;
      }
      if (nup > 0) { // This is for the PF_LO regions. Does not seem to break core_L or core_R
                  // However I need to think thos through more. I think it is ok only when xpt
                  // is known quite precisely.
        *zmax = zup;
        break;
      }
      if (nlo>=2){
        if(fabs(zlo-zup)<1e-12){
          *zmax = zlo;
          break;
        }
        zlo_last = zlo;
        zlo = (zlo+zup)/2;
      }
      if(nlo==1){
        zup = zlo;
        zlo = zlo_last;
      }
      if(nlo==0){
        zup = zlo;
        zlo = zlo_last;
      }
    }
}




// Sets zmax if plate is specified
void set_upper_plate(struct gkyl_tok_geo *geo, struct arc_length_ctx* arc_ctx, struct plate_ctx* pctx, double psi_curr)
{
      double rzplate[2];
      pctx->psi_curr = psi_curr;
      pctx->lower=false;
      double a = 0;
      double b = 1;
      double fa = tok_plate_psi_func(a, pctx);
      double fb = tok_plate_psi_func(b, pctx);
      struct gkyl_qr_res res = gkyl_ridders(tok_plate_psi_func, pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smax = res.res;
      geo->plate_func_upper(smax, rzplate);
      arc_ctx->zmax = rzplate[1];
}

// Sets zmin if plate is specified
void set_lower_plate(struct gkyl_tok_geo *geo, struct arc_length_ctx* arc_ctx, struct plate_ctx* pctx, double psi_curr)
{
      double rzplate[2];
      pctx->psi_curr = psi_curr;
      pctx->lower=true;
      double a = 0;
      double b = 1;
      double fa = tok_plate_psi_func(a, pctx);
      double fb = tok_plate_psi_func(b, pctx);
      struct gkyl_qr_res res = gkyl_ridders(tok_plate_psi_func, pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smin = res.res;
      geo->plate_func_lower(smin, rzplate);
      arc_ctx->zmin = rzplate[1];
}

// Sets zmax if plate is specified
void set_upper_iwl_plate(struct gkyl_tok_geo *geo, struct arc_length_ctx* arc_ctx, struct plate_ctx* pctx, double psi_curr)
{
      double rzplate[2];
      pctx->psi_curr = psi_curr;
      pctx->lower=false;
      double a = 0;
      double b = 1;
      double fa = tok_plate_psi_func(a, pctx);
      double fb = tok_plate_psi_func(b, pctx);
      struct gkyl_qr_res res = gkyl_ridders(tok_plate_psi_func, pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smax = res.res;
      geo->plate_func_upper(smax, rzplate);
      arc_ctx->zmax_iwl_plate = rzplate[1];
      geo->rmin = rzplate[0];
}

// Sets zmin if plate is specified
void set_lower_iwl_plate(struct gkyl_tok_geo *geo, struct arc_length_ctx* arc_ctx, struct plate_ctx* pctx, double psi_curr)
{
      double rzplate[2];
      pctx->psi_curr = psi_curr;
      pctx->lower=true;
      double a = 0;
      double b = 1;
      double fa = tok_plate_psi_func(a, pctx);
      double fb = tok_plate_psi_func(b, pctx);
      struct gkyl_qr_res res = gkyl_ridders(tok_plate_psi_func, pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smin = res.res;
      geo->plate_func_lower(smin, rzplate);
      arc_ctx->zmin_iwl_plate = rzplate[1];
      geo->rmin = rzplate[0];
}

void 
tok_geo_set_extent(struct gkyl_tok_geo_grid_inp* inp, struct gkyl_tok_geo *geo, double *theta_lo, double *theta_up)
{
  geo->rleft = inp->rleft;
  geo->rright = inp->rright;

  geo->inexact_roots = inp->inexact_roots;

  geo->rmax = inp->rmax;
  geo->rmin = inp->rmin;
  int nzcells;
  if(geo->use_cubics)
    nzcells = geo->rzgrid_cubic.cells[1];
  else
    nzcells = geo->rzgrid.cells[1];
  double *arc_memo = gkyl_malloc(sizeof(double[nzcells]));
  double *arc_memo_left = gkyl_malloc(sizeof(double[nzcells]));
  double *arc_memo_right = gkyl_malloc(sizeof(double[nzcells]));

  struct arc_length_ctx arc_ctx = {
    .geo = geo,
    .arc_memo = arc_memo,
    .arc_memo_right = arc_memo_right,
    .arc_memo_left = arc_memo_left,
    .ftype = inp->ftype,
    .zmaxis = geo->zmaxis
  };
  struct plate_ctx pctx = {
    .geo = geo
  };

  double del = 1.0e-14;

  if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT || inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO || inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID || inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP) {
    // Immediately set rclose
    arc_ctx.rclose = inp->rright;
    // Set zmin and zmax either fixed or with plate
    if (geo->plate_spec){
      set_upper_plate(geo, &arc_ctx, &pctx, geo->psisep);
      set_lower_plate(geo, &arc_ctx, &pctx, geo->psisep);
    }
    else{
      arc_ctx.zmin = inp->zmin;
      arc_ctx.zmax = inp->zmax;
    }
    double zxpt_up = geo->use_cubics ? geo->efit->Zxpt_cubic[1] : geo->efit->Zxpt[1];
    double zxpt_lo = geo->use_cubics ? geo->efit->Zxpt_cubic[0] : geo->efit->Zxpt[0];
    // Set the arc length
    double arcL_tot = integrate_psi_contour_memo(geo, geo->psisep, arc_ctx.zmin, arc_ctx.zmax, arc_ctx.rclose, false, false, arc_memo);
    double arcL_lo = integrate_psi_contour_memo(geo, geo->psisep, arc_ctx.zmin, zxpt_lo, arc_ctx.rclose, false, false, arc_memo);
    double arcL_mid = integrate_psi_contour_memo(geo, geo->psisep, zxpt_lo, zxpt_up, arc_ctx.rclose, false, false, arc_memo);
    double arcL_up = integrate_psi_contour_memo(geo, geo->psisep, zxpt_up, arc_ctx.zmax, arc_ctx.rclose, false, false, arc_memo);
    if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT) {
      *theta_lo = -M_PI+del;
      *theta_up = M_PI-del;
    }
    else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO) {
      *theta_lo = -M_PI+del;
      *theta_up = -M_PI+del + arcL_lo/arcL_tot*2.0*M_PI;
    }
    else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID) {
      *theta_lo = -M_PI+del + arcL_lo/arcL_tot*2.0*M_PI;
      *theta_up = inp->half_domain ? 0.0 : M_PI-del - arcL_up/arcL_tot*2.0*M_PI;
    }
    else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP) {
      *theta_lo = M_PI-del - arcL_up/arcL_tot*2.0*M_PI;
      *theta_up = M_PI-del;
    }
  }

  else if(inp->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN || inp->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO || inp->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID || inp->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP){
    // Immediately set rclose
    arc_ctx.rclose = inp->rleft;
    // Set zmin and zmax either fixed or with plate
    if (geo->plate_spec){
      set_upper_plate(geo, &arc_ctx, &pctx, geo->psisep);
      set_lower_plate(geo, &arc_ctx, &pctx, geo->psisep);
    }
    else{
      arc_ctx.zmin = inp->zmin;
      arc_ctx.zmax = inp->zmax;
    }
    double zxpt_up = geo->use_cubics ? geo->efit->Zxpt_cubic[1] : geo->efit->Zxpt[1];
    double zxpt_lo = geo->use_cubics ? geo->efit->Zxpt_cubic[0] : geo->efit->Zxpt[0];
    // Set the arc Length
    double arcL_tot = integrate_psi_contour_memo(geo, geo->psisep, arc_ctx.zmin, arc_ctx.zmax, arc_ctx.rclose, false, false, arc_memo);
    double arcL_lo = integrate_psi_contour_memo(geo, geo->psisep, arc_ctx.zmin, zxpt_lo, arc_ctx.rclose, false, false, arc_memo);
    double arcL_mid = integrate_psi_contour_memo(geo, geo->psisep, zxpt_lo, zxpt_up, arc_ctx.rclose, false, false, arc_memo);
    double arcL_up = integrate_psi_contour_memo(geo, geo->psisep, zxpt_up, arc_ctx.zmax, arc_ctx.rclose, false, false, arc_memo);
    if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN) {
      *theta_lo = -M_PI+del;
      *theta_up = M_PI-del;
    }
    else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP) {
      *theta_lo = -M_PI+del;
      *theta_up = -M_PI+del + arcL_lo/arcL_tot*2.0*M_PI;
    }
    else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID) {
      *theta_lo = inp->half_domain ? 0.0 : -M_PI+del + arcL_lo/arcL_tot*2.0*M_PI;
      *theta_up = M_PI-del - arcL_up/arcL_tot*2.0*M_PI;
    }
    else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO) {
      *theta_lo = M_PI-del - arcL_up/arcL_tot*2.0*M_PI;
      *theta_up = M_PI-del;
    }
  }
  else if(inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE || inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R || inp->ftype ==  GKYL_GEOMETRY_TOKAMAK_CORE_L){
    // Immediately set rleft and rright. Will need both
    arc_ctx.rright = inp->rright;
    arc_ctx.rleft = inp->rleft;

    double zxpt_up = geo->use_cubics ? geo->efit->Zxpt_cubic[1] : geo->efit->Zxpt[1];
    double zxpt_lo = geo->use_cubics ? geo->efit->Zxpt_cubic[0] : geo->efit->Zxpt[0];
    arc_ctx.zmax = inp->zmax ? inp->zmax : zxpt_up; // Initial guess.
                                                  // zmax is specified for single null full core
    double zlo = geo->zmaxis;
    find_upper_turning_point(geo, geo->psisep, zlo, &arc_ctx.zmax, 0);
    arc_ctx.zmin = zxpt_lo; // Initial guess
    double zup = geo->zmaxis;
    find_lower_turning_point(geo, geo->psisep, zup, &arc_ctx.zmin, 0);
    // Done finding turning points
    arc_ctx.right = true;
    double arcL_r = integrate_psi_contour_memo(geo, geo->psisep, arc_ctx.zmin, arc_ctx.zmax, arc_ctx.rright,
      false, false, arc_memo_right);
    arc_ctx.right = false;
    double arcL_l = integrate_psi_contour_memo(geo, geo->psisep, arc_ctx.zmin, arc_ctx.zmax, arc_ctx.rleft,
      false, false, arc_memo_left);
    double arcL_tot = arcL_l + arcL_r;

    if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE) {
      *theta_lo = -M_PI+del;
      *theta_up = M_PI-del;
    }
    else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R) {
      *theta_lo = -M_PI+del;
      *theta_up = inp->half_domain? -M_PI+del + arcL_r/arcL_tot*2.0*M_PI/2.0 : -M_PI+del + arcL_r/arcL_tot*2.0*M_PI;
    }
    else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L) {
      *theta_lo = inp->half_domain ? M_PI-del - arcL_l/arcL_tot*2.0*M_PI/2.0 : M_PI-del - arcL_l/arcL_tot*2.0*M_PI;
      *theta_up = M_PI-del;
    }
  }
  else if(inp->ftype==GKYL_GEOMETRY_TOKAMAK_LSN_SOL || inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO || inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID || inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP){
    // Immediately set rleft and rright. Will need both
    arc_ctx.rright = inp->rright;
    arc_ctx.rleft = inp->rleft;

    //Find the  upper turning point
    arc_ctx.zmax = inp->zmax; // Initial guess
    double zlo = fmax(inp->zmin_left, inp->zmin_right);
    find_upper_turning_point(geo, geo->psisep, zlo, &arc_ctx.zmax, 0);

    // Set zmin left and zmin right wither with plate or fixed
    // This one can't be used with the general func for setting upper and lower plates because it uses zmin left and zmin right
    if (geo->plate_spec){
      double rzplate[2];
      pctx.psi_curr = geo->psisep;
      pctx.lower=false;
      double a = 0;
      double b = 1;
      double fa = tok_plate_psi_func(a, &pctx);
      double fb = tok_plate_psi_func(b, &pctx);
      struct gkyl_qr_res res = gkyl_ridders(tok_plate_psi_func, &pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smax = res.res;
      geo->plate_func_upper(smax, rzplate);
      arc_ctx.zmin_left = rzplate[1];

      pctx.lower=true;
      a = 0;
      b = 1;
      fa = tok_plate_psi_func(a, &pctx);
      fb = tok_plate_psi_func(b, &pctx);
      res = gkyl_ridders(tok_plate_psi_func, &pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smin = res.res;
      geo->plate_func_lower(smin, rzplate);
      arc_ctx.zmin_right = rzplate[1];
    }
    else{
      arc_ctx.zmin_left = inp->zmin_left;
      arc_ctx.zmin_right = inp->zmin_right;
    }

    // Done finding turning points
    double zxpt = geo->use_cubics ? geo->efit->Zxpt_cubic[0] : geo->efit->Zxpt[0];

    arc_ctx.right = true;
    double arcL_mid_r = integrate_psi_contour_memo(geo, geo->psisep, zxpt, arc_ctx.zmax, arc_ctx.rright,
      false, false, arc_memo_right);
    double arcL_lo = integrate_psi_contour_memo(geo, geo->psisep, arc_ctx.zmin_right, zxpt, arc_ctx.rright,
      false, false, arc_memo_right);
    arc_ctx.right = false;
    double arcL_mid_l = integrate_psi_contour_memo(geo, geo->psisep, zxpt, arc_ctx.zmax, arc_ctx.rleft,
      false, false, arc_memo_right);
    double arcL_up = integrate_psi_contour_memo(geo, geo->psisep, arc_ctx.zmin_left, zxpt, arc_ctx.rleft,
      false, false, arc_memo_left);
    double arcL_tot = arcL_lo + arcL_mid_l + arcL_mid_r + arcL_up;

    if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL) {
      *theta_lo = -M_PI+del;
      *theta_up = M_PI-del;
    }
    else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO) {
      *theta_lo = -M_PI+del;
      *theta_up = -M_PI+del + arcL_lo/arcL_tot*2.0*M_PI;
    }
    else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID) {
      *theta_lo = -M_PI+del + arcL_lo/arcL_tot*2.0*M_PI;
      double theta_pi = inp->straight_xpt_ray && !inp->half_domain
        ? M_PI-del : M_PI+del;
      *theta_up = theta_pi - arcL_up/arcL_tot*2.0*M_PI;
    }
    else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP) {
      double theta_pi = inp->straight_xpt_ray && !inp->half_domain
        ? M_PI-del : M_PI+del;
      *theta_lo = theta_pi - arcL_up/arcL_tot*2.0*M_PI;
      *theta_up = M_PI-del;
    }
  }

  else if(inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_R || inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_L){
    arc_ctx.rright = inp->rright;
    arc_ctx.rleft = inp->rleft;

    //Find the  upper turning point to set zmax
    double zxpt_lo = geo->use_cubics ? geo->efit->Zxpt_cubic[0] : geo->efit->Zxpt[0];
    arc_ctx.zmax = zxpt_lo; // Initial guess
    double zlo = fmin(inp->zmin_left, inp->zmin_right);
    find_upper_turning_point(geo, geo->psisep, zlo, &arc_ctx.zmax, 1e-15);

    // Set zmin left and zmin right wither with plate or fixed
    // This one can't be used with the general func for setting upper and lower plates because it uses zmin left and zmin right
    if (geo->plate_spec){
      double rzplate[2];
      pctx.psi_curr = geo->psisep;
      pctx.lower=false;
      double a = 0;
      double b = 1;
      double fa = tok_plate_psi_func(a, &pctx);
      double fb = tok_plate_psi_func(b, &pctx);
      struct gkyl_qr_res res = gkyl_ridders(tok_plate_psi_func, &pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smax = res.res;
      geo->plate_func_upper(smax, rzplate);
      arc_ctx.zmin_left = rzplate[1];


      pctx.lower=true;
      a = 0;
      b = 1;
      fa = tok_plate_psi_func(a, &pctx);
      fb = tok_plate_psi_func(b, &pctx);
      res = gkyl_ridders(tok_plate_psi_func, &pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smin = res.res;
      geo->plate_func_lower(smin, rzplate);
      arc_ctx.zmin_right = rzplate[1];
    }
    else{
      arc_ctx.zmin_left = inp->zmin_left;
      arc_ctx.zmin_right = inp->zmin_right;
    }

    // Set arc length
    arc_ctx.rclose = inp->rright;
    arc_ctx.right = true;
    double arcL_r = integrate_psi_contour_memo(geo, geo->psisep, arc_ctx.zmin_right, arc_ctx.zmax, arc_ctx.rright,
      false, false, arc_memo_right);

    // Immediately set rclose
    arc_ctx.rclose = inp->rleft;
    arc_ctx.right = false;
    double arcL_l = integrate_psi_contour_memo(geo, geo->psisep, arc_ctx.zmin_left, arc_ctx.zmax, arc_ctx.rleft,
      false, false, arc_memo_left);
    double arcL_tot = arcL_l + arcL_r;

    if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_R) {
      *theta_lo = -M_PI+del;
      *theta_up = -M_PI+del + arcL_r/arcL_tot*2.0*M_PI;
    }
    else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_L) {
      *theta_lo = M_PI-del - arcL_l/arcL_tot*2.0*M_PI;
      *theta_up = M_PI-del;
    }
  }

  else if(inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_UP_L || inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_UP_R){
    arc_ctx.rright = inp->rright;
    arc_ctx.rleft = inp->rleft;
    //Find the lower turning point to set zmin
    double zxpt_up = geo->use_cubics ? geo->efit->Zxpt_cubic[1] : geo->efit->Zxpt[1];
    arc_ctx.zmin = zxpt_up; // Initial guess
    double zup = fmax(inp->zmax_left,  inp->zmax_right);
    find_lower_turning_point(geo, geo->psisep, zup, &arc_ctx.zmin, 1e-15);
    // Done finding turning point

    // Set zmax left and zmax right wither with plate or fixed
    // This one can't be used with the general func for setting upper and lower plates because it uses zmax left and zmax right
    if (geo->plate_spec){
      double rzplate[2];
      pctx.psi_curr = geo->psisep;
      pctx.lower=false;
      double a = 0;
      double b = 1;
      double fa = tok_plate_psi_func(a, &pctx);
      double fb = tok_plate_psi_func(b, &pctx);
      struct gkyl_qr_res res = gkyl_ridders(tok_plate_psi_func, &pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smax = res.res;
      geo->plate_func_upper(smax, rzplate);
      arc_ctx.zmax_right= rzplate[1];

      pctx.lower=true;
      a = 0;
      b = 1;
      fa = tok_plate_psi_func(a, &pctx);
      fb = tok_plate_psi_func(b, &pctx);
      res = gkyl_ridders(tok_plate_psi_func, &pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smin = res.res;
      geo->plate_func_lower(smin, rzplate);
      arc_ctx.zmax_left= rzplate[1];
    }
    else{
      arc_ctx.zmax_left = inp->zmax_left;
      arc_ctx.zmax_right = inp->zmax_right;
    }

    // Immediately set rclose
    arc_ctx.rclose = inp->rleft;
    arc_ctx.right = false;
    double arcL_l = integrate_psi_contour_memo(geo, geo->psisep, arc_ctx.zmin, arc_ctx.zmax_left, arc_ctx.rleft,
      false, false, arc_memo_left);

    // Immediately set rclose
    arc_ctx.rclose = inp->rright;
    arc_ctx.right = true;
    double arcL_r = integrate_psi_contour_memo(geo, geo->psisep, arc_ctx.zmin, arc_ctx.zmax_right, arc_ctx.rright,
      false, false, arc_memo_right);
    double arcL_tot = arcL_r + arcL_l;

    if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_UP_L) {
      *theta_lo = -M_PI+del;
      *theta_up = -M_PI+del + arcL_l/arcL_tot*2.0*M_PI;
    }
    else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_UP_R) {
      *theta_lo = M_PI-del - arcL_r/arcL_tot*2.0*M_PI;
      *theta_up = M_PI-del;
    }

  }

  gkyl_free(arc_memo);
  gkyl_free(arc_memo_left);
  gkyl_free(arc_memo_right);

}

void
tok_find_endpoints(struct gkyl_tok_geo_grid_inp* inp, struct gkyl_tok_geo *geo, struct arc_length_ctx* arc_ctx, struct plate_ctx* pctx, double psi_curr, double alpha_curr, double* arc_memo, double* arc_memo_left, double* arc_memo_right){
  enum { PH_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates


  // Set psicurr no matter what
  arc_ctx->psi = psi_curr;
  arc_ctx->xpt_map_valid = false;
  arc_ctx->xpt_anchor_valid = false;

  if(inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE || inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R || inp->ftype ==  GKYL_GEOMETRY_TOKAMAK_CORE_L){
    // Immediately set rleft and rright. Will need both
    arc_ctx->rright = inp->rright;
    arc_ctx->rleft = inp->rleft;

    double zxpt_up = geo->use_cubics ? geo->efit->Zxpt_cubic[1] : geo->efit->Zxpt[1];
    double zxpt_lo = geo->use_cubics ? geo->efit->Zxpt_cubic[0] : geo->efit->Zxpt[0];
    arc_ctx->zmax = inp->zmax ? inp->zmax : zxpt_up; // Initial guess.
                                                  // zmax is specified for single null full core
    double zlo = geo->zmaxis;
    find_upper_turning_point(geo, psi_curr, zlo, &arc_ctx->zmax, 0);
    arc_ctx->zmin = zxpt_lo; // Initial guess
    double zup = geo->zmaxis;
    find_lower_turning_point(geo, psi_curr, zup, &arc_ctx->zmin, 0);
    if (fabs(psi_curr-geo->psisep) <=
        128.0*DBL_EPSILON*fmax(1.0, fabs(geo->psisep)))
      arc_ctx->zmin = zxpt_lo;
    // Done finding turning points
    arc_ctx->arcL_right = integrate_psi_contour_memo(geo, psi_curr, arc_ctx->zmin, arc_ctx->zmax, arc_ctx->rright,
      true, true, arc_memo_right);
    arc_ctx->right = false;
    arc_ctx->arcL_left = integrate_psi_contour_memo(geo, psi_curr, arc_ctx->zmin, arc_ctx->zmax, arc_ctx->rleft,
      true, true, arc_memo_left);
    arc_ctx->arcL_tot = arc_ctx->arcL_left + arc_ctx->arcL_right;


    // Adjust the starting point (theta=0) so that the core blocks
    // are up-down symmetric
    if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R) {
      double theta_extent = inp->half_domain ?  2.0*(inp->cgrid.upper[TH_IDX] - inp->cgrid.lower[TH_IDX]) : inp->cgrid.upper[TH_IDX] - inp->cgrid.lower[TH_IDX];
      double arcL_extent = theta_extent/(2.0*M_PI)*arc_ctx->arcL_tot;
      double extra_arcL = arcL_extent - arc_ctx->arcL_right;
      arc_ctx->arcL_start = extra_arcL/2.0;
    }
    else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L){
      double theta_extent = inp->half_domain ? 2.0*M_PI - 2.0*(inp->cgrid.upper[TH_IDX] - inp->cgrid.lower[TH_IDX]) : 2.0*M_PI - (inp->cgrid.upper[TH_IDX] - inp->cgrid.lower[TH_IDX]);
      double arcL_extent = theta_extent/(2.0*M_PI)*arc_ctx->arcL_tot;
      double extra_arcL = arcL_extent - arc_ctx->arcL_right;
      arc_ctx->arcL_start = extra_arcL/2.0;
    }
    arc_ctx->right = true;
    arc_ctx->phi_right = 0.0;
    arc_ctx->rclose = arc_ctx->rright;
    arc_ctx->phi_right = phi_func(alpha_curr, arc_ctx->zmax, arc_ctx) - alpha_curr;

    if(inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L) {
      arc_ctx->right = false;
      arc_ctx->rclose = inp->rleft;
    }

  }

  else if(inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_R || inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_L){
    arc_ctx->rright = inp->rright;
    arc_ctx->rleft = inp->rleft;

    //Find the  upper turning point to set zmax
    double zxpt_lo = geo->use_cubics ? geo->efit->Zxpt_cubic[0] : geo->efit->Zxpt[0];
    arc_ctx->zmax = zxpt_lo; // Initial guess
    double zlo = fmin(inp->zmin_left, inp->zmin_right);
    find_upper_turning_point(geo, psi_curr, zlo, &arc_ctx->zmax, 1e-15);
    if (fabs(psi_curr-geo->psisep) <=
        128.0*DBL_EPSILON*fmax(1.0, fabs(geo->psisep)))
      arc_ctx->zmax = zxpt_lo;

    // Set zmin left and zmin right wither with plate or fixed
    // This one can't be used with the general func for setting upper and lower plates because it uses zmin left and zmin right
    if (geo->plate_spec){
      double rzplate[2];
      pctx->psi_curr = psi_curr;
      pctx->lower=false;
      double a = 0;
      double b = 1;
      double fa = tok_plate_psi_func(a, pctx);
      double fb = tok_plate_psi_func(b, pctx);
      struct gkyl_qr_res res = gkyl_ridders(tok_plate_psi_func, pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smax = res.res;
      geo->plate_func_upper(smax, rzplate);
      arc_ctx->zmin_left = rzplate[1];

      pctx->lower=true;
      a = 0;
      b = 1;
      fa = tok_plate_psi_func(a, pctx);
      fb = tok_plate_psi_func(b, pctx);
      res = gkyl_ridders(tok_plate_psi_func, pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smin = res.res;
      geo->plate_func_lower(smin, rzplate);
      arc_ctx->zmin_right = rzplate[1];
    }
    else{
      arc_ctx->zmin_left = inp->zmin_left;
      arc_ctx->zmin_right = inp->zmin_right;
    }

    // Immediately set rclose
    arc_ctx->rclose = inp->rright;
    arc_ctx->right = true;
    // Set arc length
    arc_ctx->arcL_right = integrate_psi_contour_memo(geo, psi_curr, arc_ctx->zmin_right, arc_ctx->zmax, arc_ctx->rright,
      true, true, arc_memo_right);
    // Immediately set rclose
    arc_ctx->right = false;
    arc_ctx->rclose = inp->rleft;
    double arcL_l = integrate_psi_contour_memo(geo, psi_curr, arc_ctx->zmin_left, arc_ctx->zmax, arc_ctx->rleft,
      true, true, arc_memo_left);
    arc_ctx->arcL_tot = arcL_l + arc_ctx->arcL_right;

    if(inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_R) {
      arc_ctx->right = true;
      arc_ctx->rclose = inp->rright;
    }
    else if(inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_L) {
      arc_ctx->right = false;
      arc_ctx->rclose = inp->rleft;
    }

  }

  else if(inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_UP_L || inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_UP_R){
    arc_ctx->rright = inp->rright;
    arc_ctx->rleft = inp->rleft;
    //Find the lower turning point to set zmin
    double zxpt_up = geo->use_cubics ? geo->efit->Zxpt_cubic[1] : geo->efit->Zxpt[1];
    arc_ctx->zmin = zxpt_up; // Initial guess
    double zup = fmax(inp->zmax_left,  inp->zmax_right);
    find_lower_turning_point(geo, psi_curr, zup, &arc_ctx->zmin, 1e-15);
    // Done finding turning point

    // Set zmax left and zmax right wither with plate or fixed
    // This one can't be used with the general func for setting upper and lower plates because it uses zmax left and zmax right
    if (geo->plate_spec){
      double rzplate[2];
      pctx->psi_curr = psi_curr;
      pctx->lower=false;
      double a = 0;
      double b = 1;
      double fa = tok_plate_psi_func(a, pctx);
      double fb = tok_plate_psi_func(b, pctx);
      struct gkyl_qr_res res = gkyl_ridders(tok_plate_psi_func, pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smax = res.res;
      geo->plate_func_upper(smax, rzplate);
      arc_ctx->zmax_right= rzplate[1];

      pctx->lower=true;
      a = 0;
      b = 1;
      fa = tok_plate_psi_func(a, pctx);
      fb = tok_plate_psi_func(b, pctx);
      res = gkyl_ridders(tok_plate_psi_func, pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smin = res.res;
      geo->plate_func_lower(smin, rzplate);
      arc_ctx->zmax_left= rzplate[1];
    }
    else{
      arc_ctx->zmax_left = inp->zmax_left;
      arc_ctx->zmax_right = inp->zmax_right;
    }


    // Immediately set rclose
    arc_ctx->rclose = inp->rleft;
    arc_ctx->right = false;
    arc_ctx->arcL_left = integrate_psi_contour_memo(geo, psi_curr, arc_ctx->zmin, arc_ctx->zmax_left, arc_ctx->rleft,
      true, true, arc_memo_left);

    // Immediately set rclose
    arc_ctx->rclose = inp->rright;
    arc_ctx->right = true;
    double arcL_r = integrate_psi_contour_memo(geo, psi_curr, arc_ctx->zmin, arc_ctx->zmax_right, arc_ctx->rright,
      true, true, arc_memo_right);
    arc_ctx->arcL_tot = arcL_r + arc_ctx->arcL_left;

    if(inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_UP_R) {
      arc_ctx->right = true;
      arc_ctx->rclose = inp->rright;
    }
    else if(inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_UP_L) {
      arc_ctx->right = false;
      arc_ctx->rclose = inp->rleft;
    }
  }

  else if(inp->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT || inp->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO || inp->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID || inp->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP){
    // Immediately set rclose
    arc_ctx->rclose = inp->rright;
    // Set zmin and zmax either fixed or with plate
    if (geo->plate_spec){
      set_upper_plate(geo, arc_ctx, pctx, arc_ctx->psi);
      set_lower_plate(geo, arc_ctx, pctx, arc_ctx->psi);
    }
    else{
      arc_ctx->zmin = inp->zmin;
      arc_ctx->zmax = inp->zmax;
    }
    // Set the arc length
    arc_ctx->arcL_tot = integrate_psi_contour_memo(geo, psi_curr, arc_ctx->zmin, arc_ctx->zmax, arc_ctx->rclose, true, true, arc_memo);
  }

  else if(inp->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN || inp->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO || inp->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID || inp->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP){
    // Immediately set rclose
    arc_ctx->rclose = inp->rleft;
    // Set zmin and zmax either fixed or with plate
    if (geo->plate_spec){
      set_upper_plate(geo, arc_ctx, pctx, arc_ctx->psi);
      set_lower_plate(geo, arc_ctx, pctx, arc_ctx->psi);
    }
    else{
      arc_ctx->zmin = inp->zmin;
      arc_ctx->zmax = inp->zmax;
    }
    // Set the arc Length
    arc_ctx->arcL_tot = integrate_psi_contour_memo(geo, psi_curr, arc_ctx->zmin, arc_ctx->zmax, arc_ctx->rclose, true, true, arc_memo);
  }

  else if(inp->ftype==GKYL_GEOMETRY_TOKAMAK_LSN_SOL || inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO || inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID || inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP){
    // Immediately set rleft and rright. Will need both
    arc_ctx->rright = inp->rright;
    arc_ctx->rleft = inp->rleft;
    //Find the  upper turning point
    arc_ctx->zmax = inp->zmax; // Initial guess
    double zlo = fmax(inp->zmin_left, inp->zmin_right);
    find_upper_turning_point(geo, psi_curr, zlo, &arc_ctx->zmax, 0);

    // Set zmin left and zmin right wither with plate or fixed
    // This one can't be used with the general func for setting upper and lower plates because it uses zmin left and zmin right
    if (geo->plate_spec){
      double rzplate[2];
      pctx->psi_curr = psi_curr;
      pctx->lower=false;
      double a = 0;
      double b = 1;
      double fa = tok_plate_psi_func(a, pctx);
      double fb = tok_plate_psi_func(b, pctx);
      struct gkyl_qr_res res = gkyl_ridders(tok_plate_psi_func, pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smax = res.res;
      geo->plate_func_upper(smax, rzplate);
      arc_ctx->zmin_left = rzplate[1];

      pctx->lower=true;
      a = 0;
      b = 1;
      fa = tok_plate_psi_func(a, pctx);
      fb = tok_plate_psi_func(b, pctx);
      res = gkyl_ridders(tok_plate_psi_func, pctx,
        a, b, fa, fb, geo->root_param.max_iter, 1e-10);
      double smin = res.res;
      geo->plate_func_lower(smin, rzplate);
      arc_ctx->zmin_right = rzplate[1];
    }
    else{
      arc_ctx->zmin_left = inp->zmin_left;
      arc_ctx->zmin_right = inp->zmin_right;
    }

    // Done finding turning point
    arc_ctx->arcL_right = integrate_psi_contour_memo(geo, psi_curr, arc_ctx->zmin_right, arc_ctx->zmax, arc_ctx->rright,
      true, true, arc_memo_right);
    arc_ctx->right = false;
    double arcL_l = integrate_psi_contour_memo(geo, psi_curr, arc_ctx->zmin_left, arc_ctx->zmax, arc_ctx->rleft,
      true, true, arc_memo_left);
    arc_ctx->arcL_tot = arcL_l + arc_ctx->arcL_right;

    arc_ctx->right = true;
    arc_ctx->phi_right = 0.0;
    arc_ctx->rclose = arc_ctx->rright;
    arc_ctx->psi = psi_curr;
    arc_ctx->zmin = arc_ctx->zmin_right;
    arc_ctx->phi_right = phi_func(alpha_curr, arc_ctx->zmax, arc_ctx) - alpha_curr;
  }

  else if(inp->ftype == GKYL_GEOMETRY_TOKAMAK_IWL){
    // Immediately set rleft and rright. Will need both
    arc_ctx->rright = inp->rright;
    arc_ctx->rleft = inp->rleft;

    double R_lcfs[4], dRdZ_lcfs[4];
    double dR_lcfs[4], dZ_lcfs[4];
    int nr_lcfs = gkyl_tok_geo_R_psiZ(geo, geo->efit->sibry, geo->efit->zmaxis, 4, R_lcfs, dRdZ_lcfs, dR_lcfs, dZ_lcfs);
    if (nr_lcfs <= 0) {
      fprintf(stderr,
        "TOK_GEOMETRY IWL failed to find LCFS root at psi=%.17g Z=%.17g\n",
        geo->efit->sibry, geo->efit->zmaxis);
      abort();
    }
    double r_lcfs = nr_lcfs == 1 ? R_lcfs[0] : choose_closest(arc_ctx->rleft, R_lcfs, R_lcfs, nr_lcfs);
    double rz_lcfs[2];

    double rmin_old = geo->rmin;
    if (geo->plate_spec) {
      geo->plate_func_upper(0.0, rz_lcfs);
      if(fabs(rz_lcfs[0] - r_lcfs) > 1e-6) {
        fprintf(stderr, "The upper plate function has an error. It must return (R(s=0),Z(s=0)) = (%1.16f, %1.16f). \n", R_lcfs[0], geo->efit->zmaxis);
        assert(false);
      }
      geo->plate_func_lower(0.0, rz_lcfs);
      if(fabs(rz_lcfs[0] - r_lcfs) > 1e-6) {
        fprintf(stderr, "The lower plate function has an error. It must return (R(s=0),Z(s=0)) = (%1.16f, %1.16f). \n", R_lcfs[0], geo->efit->zmaxis);
        assert(false);
      }
  
      if ( (arc_ctx->psi >= geo->efit->sibry && geo->efit->sibry >= geo->efit->simag) ||
           (arc_ctx->psi <= geo->efit->sibry && geo->efit->sibry <= geo->efit->simag) ) {
        set_upper_iwl_plate(geo, arc_ctx, pctx, arc_ctx->psi);
        set_lower_iwl_plate(geo, arc_ctx, pctx, arc_ctx->psi);
      }
    }
    else {
      arc_ctx->zmin_iwl_plate = geo->zmaxis;
      arc_ctx->zmax_iwl_plate = geo->zmaxis;
    }

    arc_ctx->zmax = inp->zmax; // Initial guess.

    double zlo = geo->zmaxis;
    find_upper_turning_point(geo, psi_curr, zlo, &arc_ctx->zmax, 0);
    arc_ctx->zmin = inp->zmin; // Initial guess
    double zup = geo->zmaxis;
    find_lower_turning_point(geo, psi_curr, zup, &arc_ctx->zmin, 0);
    // Done finding turning points
    arc_ctx->zmin_iwl = arc_ctx->zmin;
    arc_ctx->zmax_iwl = arc_ctx->zmax;
    geo->rmin = rmin_old;

    arc_ctx->right = true;
    arc_ctx->arcL_q1 = integrate_psi_contour_memo(geo, psi_curr, geo->zmaxis, arc_ctx->zmax, arc_ctx->rright, false, false, arc_memo);
    arc_ctx->right = false;
    arc_ctx->arcL_q2 = integrate_psi_contour_memo(geo, psi_curr, arc_ctx->zmin_iwl_plate, arc_ctx->zmax, arc_ctx->rleft, false, false, arc_memo);
    arc_ctx->right = false;
    arc_ctx->arcL_q3 = integrate_psi_contour_memo(geo, psi_curr, arc_ctx->zmin, arc_ctx->zmax_iwl_plate, arc_ctx->rleft, false, false, arc_memo);
    arc_ctx->right = true;
    arc_ctx->arcL_q4 = integrate_psi_contour_memo(geo, psi_curr, arc_ctx->zmin, geo->zmaxis, arc_ctx->rright, false, false, arc_memo);

    arc_ctx->arcL_tot = arc_ctx->arcL_q1 + arc_ctx->arcL_q2 + arc_ctx->arcL_q3 + arc_ctx->arcL_q4;

    arc_ctx->right = true;
    arc_ctx->rclose = arc_ctx->rright;
    arc_ctx->phi_right = phi_func(alpha_curr, arc_ctx->zmax, arc_ctx) - alpha_curr;
  }

  tok_configure_xpt_map(inp, arc_ctx);
}



void
tok_set_ridders(struct gkyl_tok_geo_grid_inp* inp, struct arc_length_ctx* arc_ctx, double psi_curr, double arcL_curr,double* rclose, double *ridders_min, double* ridders_max){


  if(inp->ftype==GKYL_GEOMETRY_TOKAMAK_CORE){
    if(arcL_curr <= arc_ctx->arcL_right){
      *rclose = arc_ctx->rright;
      arc_ctx->right = true;
      *ridders_min = -arcL_curr;
      *ridders_max = arc_ctx->arcL_tot-arcL_curr;
    }
    else{
      *rclose = arc_ctx->rleft;
      arc_ctx->right = false;
      *ridders_min = arc_ctx->arcL_tot - arcL_curr;
      *ridders_max = -arcL_curr + arc_ctx->arcL_right;
    }
  }

  if(inp->ftype==GKYL_GEOMETRY_TOKAMAK_CORE_R || inp->ftype==GKYL_GEOMETRY_TOKAMAK_CORE_L){
    if (arc_ctx->xpt_map_valid) {
      while (arcL_curr < 0.0)
        arcL_curr += arc_ctx->arcL_tot;
      while (arcL_curr >= arc_ctx->arcL_tot)
        arcL_curr -= arc_ctx->arcL_tot;
      arc_ctx->pre = false;
      if (arcL_curr <= arc_ctx->arcL_right) {
        *rclose = arc_ctx->rright;
        arc_ctx->right = true;
        *ridders_min = -arcL_curr;
        *ridders_max = arc_ctx->arcL_tot-arcL_curr;
      }
      else {
        *rclose = arc_ctx->rleft;
        arc_ctx->right = false;
        *ridders_min = arc_ctx->arcL_tot-arcL_curr;
        *ridders_max = arc_ctx->arcL_right-arcL_curr;
      }
    }
    else if(arcL_curr <= arc_ctx->arcL_start){
      *rclose = arc_ctx->rleft;
      arc_ctx->right = false;
      arc_ctx->pre= true;
      *ridders_min = arc_ctx->arcL_start-arcL_curr;
      *ridders_max = arc_ctx->arcL_start - arc_ctx->arcL_left - arcL_curr;
    }
    else if(arcL_curr <= arc_ctx->arcL_right + arc_ctx->arcL_start){
      *rclose = arc_ctx->rright;
      arc_ctx->right = true;
      arc_ctx->pre= false;
      *ridders_min = arc_ctx->arcL_start -arcL_curr;
      *ridders_max = arc_ctx->arcL_right + arc_ctx->arcL_start - arcL_curr;
    }
    else{
      *rclose = arc_ctx->rleft;
      arc_ctx->right = false;
      arc_ctx->pre= false;
      *ridders_min = arc_ctx->arcL_tot + arc_ctx->arcL_start - arcL_curr;
      *ridders_max = -arcL_curr + arc_ctx->arcL_right+arc_ctx->arcL_start;
    }
  }


  else if(inp->ftype==GKYL_GEOMETRY_TOKAMAK_PF_LO_R || inp->ftype==GKYL_GEOMETRY_TOKAMAK_PF_LO_L){
    if(arcL_curr <= arc_ctx->arcL_right){
      *rclose = arc_ctx->rright;
      arc_ctx->right = true;
      *ridders_min = -arcL_curr;
      *ridders_max = arc_ctx->arcL_tot-arcL_curr;
      arc_ctx->zmin = arc_ctx->zmin_right;
    }
    else{
      *rclose = arc_ctx->rleft;
      arc_ctx->right = false;
      *ridders_min = arc_ctx->arcL_tot - arcL_curr;
      *ridders_max = -arcL_curr + arc_ctx->arcL_right;
      arc_ctx->zmin = arc_ctx->zmin_left;
    }
  }

  else if(inp->ftype==GKYL_GEOMETRY_TOKAMAK_PF_UP_R || inp->ftype==GKYL_GEOMETRY_TOKAMAK_PF_UP_L){
    if(arcL_curr <= arc_ctx->arcL_left){
      *rclose = arc_ctx->rleft;
      arc_ctx->right = false;
      *ridders_min = -arcL_curr + arc_ctx->arcL_left;
      *ridders_max = -arcL_curr;
      arc_ctx->zmax = arc_ctx->zmax_left;
    }
    else{
      *rclose = arc_ctx->rright;
      arc_ctx->right = true;
      *ridders_min = arc_ctx->arcL_left - arcL_curr;
      *ridders_max = -arcL_curr + arc_ctx->arcL_tot;
      arc_ctx->zmax = arc_ctx->zmax_right;
    }
  }

  else if( (arc_ctx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT) || (arc_ctx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT) || (arc_ctx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO) || (arc_ctx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID) || (arc_ctx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP) ){
    *ridders_min = -arcL_curr;
    *ridders_max = arc_ctx->arcL_tot-arcL_curr;
    *rclose = arc_ctx->rclose;
  }
  else if( (arc_ctx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN) || (arc_ctx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN) || (arc_ctx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO) || (arc_ctx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID) || (arc_ctx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP) ){
    *ridders_min = arc_ctx->arcL_tot-arcL_curr;
    *ridders_max = -arcL_curr;
    *rclose = arc_ctx->rclose;
  }
  else if(arc_ctx->ftype==GKYL_GEOMETRY_TOKAMAK_LSN_SOL || arc_ctx->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO || arc_ctx->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID || arc_ctx->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP){
    if(arcL_curr <= arc_ctx->arcL_right){
      *rclose = arc_ctx->rright;
      arc_ctx->right = true;
      *ridders_min = -arcL_curr;
      *ridders_max = arc_ctx->arcL_tot-arcL_curr;
      arc_ctx->zmin = arc_ctx->zmin_right;
    }
    else{
      *rclose = arc_ctx->rleft;
      arc_ctx->right = false;
      *ridders_min = arc_ctx->arcL_tot - arcL_curr;
      *ridders_max = -arcL_curr + arc_ctx->arcL_right;
      arc_ctx->zmin = arc_ctx->zmin_left;
    }
  }

  else if (arc_ctx->ftype==GKYL_GEOMETRY_TOKAMAK_IWL) {
    arc_ctx->q1 = false;
    arc_ctx->q2 = false;
    arc_ctx->q3 = false;
    arc_ctx->q4 = false;
    if(arcL_curr <= arc_ctx->arcL_q3) {
      arc_ctx->q3 = true;
      *rclose = arc_ctx->rleft;
      arc_ctx->right = false;
      *ridders_min = arc_ctx->arcL_q3 - arcL_curr;
      *ridders_max = -arcL_curr;
      arc_ctx->zmin = arc_ctx->zmin_iwl;
      arc_ctx->zmax = arc_ctx->zmax_iwl_plate;
    }

    else if(arcL_curr <= ( arc_ctx->arcL_q3 + arc_ctx->arcL_q4) ) {
      arc_ctx->q4 = true;
      *rclose = arc_ctx->rright;
      arc_ctx->right = true;
      *ridders_min = arc_ctx->arcL_q3 - arcL_curr;
      *ridders_max = arc_ctx->arcL_q3 + arc_ctx->arcL_q4 - arcL_curr;
      arc_ctx->zmin = arc_ctx->zmin_iwl;
      arc_ctx->zmax = arc_ctx->geo->zmaxis;
    }

    else if(arcL_curr <= (arc_ctx->arcL_q3 + arc_ctx->arcL_q4 + arc_ctx->arcL_q1) ) {
      arc_ctx->q1 = true;
      *rclose = arc_ctx->rright;
      arc_ctx->right = true;
      *ridders_min = arc_ctx->arcL_q3 + arc_ctx->arcL_q4 - arcL_curr;
      *ridders_max = arc_ctx->arcL_q3 + arc_ctx->arcL_q4 + arc_ctx->arcL_q1 - arcL_curr;
      arc_ctx->zmin = arc_ctx->geo->zmaxis;
      arc_ctx->zmax = arc_ctx->zmax_iwl;
    }

    else {
      arc_ctx->q2 = true;
      *rclose = arc_ctx->rleft;
      arc_ctx->right = false;
      *ridders_min = arc_ctx->arcL_q3 + arc_ctx->arcL_q4 + arc_ctx->arcL_q1 + arc_ctx->arcL_q2 - arcL_curr;
      *ridders_max = arc_ctx->arcL_q3 + arc_ctx->arcL_q4 + arc_ctx->arcL_q1 - arcL_curr;
      arc_ctx->zmin = arc_ctx->zmin_iwl_plate;
      arc_ctx->zmax = arc_ctx->zmax_iwl;
    }
  }

  arc_ctx->arcL = arcL_curr;
  arc_ctx->rclose = *rclose; // This would be unnecessary for all double null block cases. Only needed for SN and full core
}
