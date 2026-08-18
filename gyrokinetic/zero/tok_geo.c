#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <assert.h>
#include <gkyl_basis.h>
#include <gkyl_math.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_nodal_ops.h>
#include <gkyl_position_map.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_priv.h>
#include <gkyl_tok_geo_priv.h>
#include <gkyl_dg_bin_ops.h>

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static bool
tok_geo_same_flux(double psi_a, double psi_b)
{
  double scale = fmax(1.0, fmax(fabs(psi_a), fabs(psi_b)));
  return fabs(psi_a-psi_b) <= 64.0*DBL_EPSILON*scale;
}

static bool
tok_geo_trace_enabled(void)
{
  const char *trace = getenv("GKYL_TOK_GEO_TRACE");
  return trace && trace[0] != '\0' && trace[0] != '0';
}

static bool
tok_geo_trace_block_selected(void)
{
  const char *filter = getenv("GKYL_TOK_GEO_TRACE_BLOCK_FILTER");
  const char *block = getenv("GKYL_TOK_GEO_TRACE_BLOCK");
  return !(filter && filter[0] != '\0') || (block && strcmp(filter, block) == 0);
}

static int
tok_geo_trace_int_env(const char *name, int fallback)
{
  const char *val = getenv(name);
  return val && val[0] != '\0' ? atoi(val) : fallback;
}

static double
tok_geo_trace_double_env(const char *name, double fallback)
{
  const char *val = getenv(name);
  return val && val[0] != '\0' ? atof(val) : fallback;
}

// Lower bound on the slope of the separatrix->far-boundary correspondence, as
// the weight of an identity map blended into the nearest-point projection.
// Raising it widens the collapsed seam cells that nearest-point projection
// produces near an X point -- 0.25 lifts the thinnest cell of the outboard-SOL
// seam column by an order of magnitude on shots 204965/204995/204997 -- but it
// perturbs the correspondence on every half-domain block, so the default stays
// at the validated value and the knob exists to sweep it.
static double
tok_trace_corr_identity_fraction(void)
{
  double f = tok_geo_trace_double_env("GKYL_TOK_TRACE_CORR_IDENTITY", 0.01);
  return isfinite(f) && f >= 0.0 && f <= 1.0 ? f : 0.01;
}

static bool
tok_ordered_map_diag_enabled(void)
{
  const char *diag = getenv("GKYL_TOK_ORDERED_MAP_DIAG");
  return diag && diag[0] != '\0' && diag[0] != '0';
}

// Restores the hard failure when tok_logical_trace_sample can find no point of
// the requested level set near a trace segment, for A/B testing the chord
// fallback that otherwise handles a severed level set.  See the fallback's
// comment at the end of that routine.
static bool
tok_logsamp_chord_fallback_disabled(void)
{
  const char *off = getenv("GKYL_TOK_LOGSAMP_NO_CHORD_FALLBACK");
  return off && off[0] != '\0' && off[0] != '0';
}

// Anchor the X-point ray at its LAST crossing of the target surface rather than
// its first.  Both are identical unless the ray crosses more than once, which
// over the 29-shot regression set happens on only 5 shots and always with the
// same shape (one closed excursion, then the real surface).
static bool
tok_ext_ray_last_crossing(void)
{
  const char *on = getenv("GKYL_TOK_EXT_RAY_LAST_CROSSING");
  return on && on[0] != '\0' && on[0] != '0';
}

static bool
tok_xpt_seam_optimizer_trial(const struct gkyl_tok_geo_grid_inp *inp)
{
  return inp->relaxed_xpt_seam_optimizer_trial &&
    inp->relaxed_xpt_seam_trial_status;
}

static void
tok_xpt_seam_trial_reject(const struct gkyl_tok_geo_grid_inp *inp,
  enum gkyl_tok_geo_xpt_seam_trial_failure reason)
{
  if (!tok_xpt_seam_optimizer_trial(inp))
    return;
  struct gkyl_tok_geo_xpt_seam_trial_status *status =
    inp->relaxed_xpt_seam_trial_status;
  if (status->first_failure_reason == GKYL_XPT_SEAM_TRIAL_OK)
    status->first_failure_reason = reason;
  switch (reason) {
    case GKYL_XPT_SEAM_TRIAL_CONTOUR:
      status->contour_valid = false;
      break;
    case GKYL_XPT_SEAM_TRIAL_BRANCH:
      status->branch_valid = false;
      break;
    case GKYL_XPT_SEAM_TRIAL_TRACE_ORDERING:
      status->trace_ordering_valid = false;
      break;
    case GKYL_XPT_SEAM_TRIAL_NONFINITE_MAP:
      status->finite_map_valid = false;
      break;
    case GKYL_XPT_SEAM_TRIAL_CELL_JACOBIAN:
      status->cell_jacobian_valid = false;
      break;
    case GKYL_XPT_SEAM_TRIAL_METRIC_JACOBIAN:
      status->jacobian_valid = false;
      break;
    case GKYL_XPT_SEAM_TRIAL_INVALID_PARAMETER:
    case GKYL_XPT_SEAM_TRIAL_REMOTE_FAILURE:
    case GKYL_XPT_SEAM_TRIAL_OK:
      break;
  }
}

static bool
tok_xpt_mapping_requested(const struct gkyl_tok_geo_grid_inp *inp)
{
  if (inp->straight_xpt_ray)
    return true;
  return inp->straight_core_xpt_ray &&
    (inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R ||
     inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L);
}

static void
tok_init_xpt_ray_target(const struct gkyl_tok_geo_grid_inp *inp,
  const struct gkyl_position_map *position_map, struct arc_length_ctx *arc_ctx)
{
  if (!tok_xpt_mapping_requested(inp))
    return;
  // Flux can increase or decrease radially, and existing inputs use both
  // [separatrix, far] and [far, separatrix] logical orderings.  Select the
  // physical radial boundary furthest from the separatrix instead of
  // inferring it from ftype or array order.
  double psi_comp_lo = inp->cgrid.lower[0];
  double psi_comp_hi = inp->cgrid.upper[0];
  double psi_lo = 0.0, psi_hi = 0.0;
  position_map->maps[0](0.0, &psi_comp_lo, &psi_lo,
    position_map->ctxs[0]);
  position_map->maps[0](0.0, &psi_comp_hi, &psi_hi,
    position_map->ctxs[0]);
  double dlo = fabs(psi_lo-arc_ctx->geo->psisep);
  double dhi = fabs(psi_hi-arc_ctx->geo->psisep);
  double scale = fmax(1.0, fmax(fabs(psi_lo), fabs(psi_hi)));
  if (!isfinite(psi_lo) || !isfinite(psi_hi) ||
      fabs(dlo-dhi) <= 64.0*DBL_EPSILON*scale) {
    fprintf(stderr,
      "TOK_ORDERED_MAP cannot identify far radial boundary ftype=%d psi_bounds=(%.17g,%.17g) psi_sep=%.17g\n",
      inp->ftype, psi_lo, psi_hi, arc_ctx->geo->psisep);
    abort();
  }
  arc_ctx->xpt_ray_psi0 = dlo > dhi ? psi_lo : psi_hi;
}

static double
tok_xpt_theta_to_arc(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx, double theta)
{
  if (!arc_ctx->xpt_map_valid)
    return (theta+M_PI)*arc_ctx->arcL_tot/(2.0*M_PI);
  double frac = (theta-inp->cgrid.lower[2])
    /(inp->cgrid.upper[2]-inp->cgrid.lower[2]);
  double arc = arc_ctx->xpt_map_arc_lo
    + frac*(arc_ctx->xpt_map_arc_hi-arc_ctx->xpt_map_arc_lo);
  if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R ||
      inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L) {
    arc = fmod(arc, arc_ctx->arcL_tot);
    if (arc < 0.0)
      arc += arc_ctx->arcL_tot;
  }
  return arc;
}

static bool
tok_xpt_at_seam(const struct gkyl_tok_geo_grid_inp *inp, int it,
  const struct gkyl_range *nrange, bool lower_is_global, bool upper_is_global)
{
  bool lower = it == nrange->lower[2] && lower_is_global;
  bool upper = it == nrange->upper[2] && upper_is_global;
  switch (inp->ftype) {
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_R:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID:
    case GKYL_GEOMETRY_TOKAMAK_CORE_L:
      return upper;
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_L:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO:
    case GKYL_GEOMETRY_TOKAMAK_CORE_R:
      return lower;
    default:
      return false;
  }
}

static bool
tok_xpt_at_fixed_edge(const struct gkyl_tok_geo_grid_inp *inp,
  const struct arc_length_ctx *arc_ctx, int it,
  const struct gkyl_range *nrange, bool lower_is_global,
  bool upper_is_global, double *z, double *root_reference)
{
  if (!arc_ctx->xpt_map_valid)
    return false;
  bool lower = it == nrange->lower[2] && lower_is_global;
  bool upper = it == nrange->upper[2] && upper_is_global;
  switch (inp->ftype) {
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_R:
      if (lower) {
        *z = arc_ctx->zmin_right; *root_reference = inp->rright;
        return true;
      }
      break;
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_L:
      if (upper) {
        *z = arc_ctx->zmin_left; *root_reference = inp->rleft;
        return true;
      }
      break;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO:
      if (lower) {
        *z = arc_ctx->zmin; *root_reference = inp->rright;
        return true;
      }
      break;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID:
    case GKYL_GEOMETRY_TOKAMAK_CORE_R:
      if (upper) {
        *z = arc_ctx->geo->zmaxis; *root_reference = inp->rright;
        return true;
      }
      break;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID:
    case GKYL_GEOMETRY_TOKAMAK_CORE_L:
      if (lower) {
        *z = arc_ctx->geo->zmaxis; *root_reference = inp->rleft;
        return true;
      }
      break;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO:
      if (upper) {
        *z = arc_ctx->zmin; *root_reference = inp->rleft;
        return true;
      }
      break;
    default:
      break;
  }
  return false;
}

static void
tok_geo_trace_surface(FILE *fp, const char *stage, int block, int ftype,
  int dir, int ip, int it, int ia, int ip_delta, double psi, double alpha,
  double theta, double arcL, double zmin, double zmax, double rclose,
  double ridders_min, double ridders_max, double res, int nevals, double elapsed)
{
  fprintf(fp, "%s,%d,%d,%d,%d,%d,%d,%d,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%d,%.17e\n",
    stage, block, ftype, dir, ip, it, ia, ip_delta, psi, alpha, theta, arcL,
    zmin, zmax, rclose, ridders_min, ridders_max, res, nevals, elapsed);
}

static void
tok_geo_trace_surface_row(const char *stage, int ftype, int dir, int ip,
  int it, int ia, int ip_delta, double psi, double alpha, double theta,
  double arcL, double zmin, double zmax, double rclose, double ridders_min,
  double ridders_max, double res, int nevals, double elapsed)
{
  if (!tok_geo_trace_enabled() || !tok_geo_trace_block_selected())
    return;

  const char *block_env = getenv("GKYL_TOK_GEO_TRACE_BLOCK");
  int block = block_env ? atoi(block_env) : -1;
  const char *fname = getenv("GKYL_TOK_GEO_TRACE_FILE");

  if (fname && fname[0] != '\0') {
    static bool wrote_header = false;
    FILE *fp = fopen(fname, "a");
    if (fp) {
      if (!wrote_header) {
        fprintf(fp, "stage,block,ftype,dir,ip,it,ia,ip_delta,psi,alpha,theta,arcL,zmin,zmax,rclose,ridders_min,ridders_max,res,nevals,elapsed_sec\n");
        wrote_header = true;
      }
      tok_geo_trace_surface(fp, stage, block, ftype, dir, ip, it, ia, ip_delta,
        psi, alpha, theta, arcL, zmin, zmax, rclose, ridders_min, ridders_max,
        res, nevals, elapsed);
      fclose(fp);
    }
  }

  if (tok_geo_trace_int_env("GKYL_TOK_GEO_TRACE_STDERR", 1)) {
    tok_geo_trace_surface(stderr, stage, block, ftype, dir, ip, it, ia, ip_delta,
      psi, alpha, theta, arcL, zmin, zmax, rclose, ridders_min, ridders_max,
      res, nevals, elapsed);
    fflush(stderr);
  }
}

double
tok_plate_psi_func(double s, void *ctx){
  // uses a pointer to the plate function to get R(s), Z(s)
  // Then calculates psi(R, Z)
  // will be used by ridders later
  
  struct plate_ctx *gc = ctx;
  double RZ[2];
  if(gc->lower==true)
    gc->geo->plate_func_lower(s, RZ);
  else
    gc->geo->plate_func_upper(s, RZ);

  double R = RZ[0];
  double Z = RZ[1];

  // Now find the cell where this R and Z is
  if (gc->geo->use_cubics) {
    double xn[2] = {R, Z};
    double psi;
    gc->geo->efit->evf->eval_cubic(0.0, xn, &psi, gc->geo->efit->evf->ctx);
    return psi - gc->psi_curr;
  }
  else {
    int rzidx[2];
    rzidx[0] = fmin(gc->geo->rzlocal.lower[0] + (int) floor((R - gc->geo->rzgrid.lower[0])/gc->geo->rzgrid.dx[0]), gc->geo->rzlocal.upper[0]);
    rzidx[1] = fmin(gc->geo->rzlocal.lower[1] + (int) floor((Z - gc->geo->rzgrid.lower[1])/gc->geo->rzgrid.dx[1]), gc->geo->rzlocal.upper[1]);
    long loc = gkyl_range_idx(&gc->geo->rzlocal, rzidx);
    const double *coeffs = gkyl_array_cfetch(gc->geo->psiRZ,loc);

    double xc[2];
    gkyl_rect_grid_cell_center(&gc->geo->rzgrid, rzidx, xc);
    double xy[2];
    xy[0] = (R-xc[0])/(gc->geo->rzgrid.dx[0]*0.5);
    xy[1] = (Z-xc[1])/(gc->geo->rzgrid.dx[1]*0.5);
    double psi = gc->geo->rzbasis.eval_expand(xy, coeffs);
    return psi - gc->psi_curr;
  }
}



// Function to pass to root-finder to find Z location for given arc-length
static inline double
arc_length_func(double Z, void *ctx)
{
  struct arc_length_ctx *actx = ctx;
  double *arc_memo;
  double psi = actx->psi, rclose = actx->rclose, zmin = actx->zmin, arcL = actx->arcL;
  double zmax = actx->zmax;
  double ival = 0.0;

  if(actx->ftype==GKYL_GEOMETRY_TOKAMAK_CORE){
    if(actx->right==true){
      double *arc_memo = actx->arc_memo_right;
      ival = integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, true, false, arc_memo) - arcL;
    }
    else{
      double *arc_memo = actx->arc_memo_left;
      ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, true, false, arc_memo)  - arcL + actx->arcL_right;
    }
  }

  else if(actx->ftype==GKYL_GEOMETRY_TOKAMAK_CORE_L || actx->ftype==GKYL_GEOMETRY_TOKAMAK_CORE_R){
    if (actx->xpt_map_valid && actx->right) {
      double *arc_memo = actx->arc_memo_right;
      ival = integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose,
        true, false, arc_memo)-arcL;
    }
    else if (actx->xpt_map_valid) {
      double *arc_memo = actx->arc_memo_left;
      ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose,
        true, false, arc_memo)-arcL+actx->arcL_right;
    }
    else if(actx->pre==true){
      ival = actx->arcL_start - integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, false, false, arc_memo) - arcL ;
    }
    else if(actx->right==true){
      double *arc_memo = actx->arc_memo_right;
      ival = integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, true, false, arc_memo) - arcL + actx->arcL_start;
    }
    else{
      double *arc_memo = actx->arc_memo_left;
      ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, true, false, arc_memo)  - arcL + actx->arcL_right + actx->arcL_start;
    }
  }

  else if(actx->ftype==GKYL_GEOMETRY_TOKAMAK_PF_LO_L || actx->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_R){
    if(actx->right==true){
      double *arc_memo = actx->arc_memo_right;
      ival = integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, true, false, arc_memo) - arcL;
    }
    else{
      double *arc_memo = actx->arc_memo_left;
      ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, true, false, arc_memo)  - arcL + actx->arcL_right;
    }
  }

  else if(actx->ftype==GKYL_GEOMETRY_TOKAMAK_PF_UP_L || actx->ftype == GKYL_GEOMETRY_TOKAMAK_PF_UP_R){
    if(actx->right==false){
      double *arc_memo = actx->arc_memo_left;
      ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, true, false, arc_memo) - arcL;
    }
    else{
      double *arc_memo = actx->arc_memo_right;
      ival = integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, true, false, arc_memo)  - arcL + actx->arcL_left;
    }
  }

  else if( (actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT) || (actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO) || (actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID) || (actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP) ){
    double *arc_memo = actx->arc_memo;
    ival = integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, true, false, arc_memo) - arcL;
  }
  else if( (actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN) || (actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN) || (actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO) || (actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID) || (actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP) ){
    double *arc_memo = actx->arc_memo;
    ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, true, false, arc_memo) - arcL;
  }

  else if(actx->ftype==GKYL_GEOMETRY_TOKAMAK_LSN_SOL || actx->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO || actx->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID || actx->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP){
    if(actx->right==true){
      double *arc_memo = actx->arc_memo_right;
      ival = integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, false, false, arc_memo) - arcL;
    }
    else{
      double *arc_memo = actx->arc_memo_left;
      ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, false, false, arc_memo)  - arcL + actx->arcL_right;
    }
  }

  else if(actx->ftype==GKYL_GEOMETRY_TOKAMAK_IWL){
    if(actx->q3) {
      double *arc_memo = actx->arc_memo;
      ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, false, false, arc_memo) - arcL;
    }
    else if (actx->q4) {
      double *arc_memo = actx->arc_memo;
      ival = integrate_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, false, false, arc_memo)  - arcL + actx->arcL_q3;
    }
    else if (actx->q1) {
      double *arc_memo = actx->arc_memo;
      ival = integrate_psi_contour_memo(actx->geo, psi, actx->geo->zmaxis, Z, rclose, false, false, arc_memo)  - arcL + actx->arcL_q3 +actx->arcL_q4;
    }
    else {
      double *arc_memo = actx->arc_memo;
      ival = integrate_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, false, false, arc_memo)  - arcL + actx->arcL_q3 +actx->arcL_q4 + actx->arcL_q1;
    }

  }

  return ival;
}

static double
tok_eval_psi_rz_local(const struct gkyl_tok_geo *geo, double R, double Z)
{
  if (geo->use_cubics) {
    double xn[2] = { R, Z }, psi = 0.0;
    geo->efit->evf->eval_cubic(0.0, xn, &psi, geo->efit->evf->ctx);
    return psi;
  }

  int idx[2];
  idx[0] = GKYL_MIN2(geo->rzlocal.upper[0], GKYL_MAX2(geo->rzlocal.lower[0],
    geo->rzlocal.lower[0]
      +(int) floor((R-geo->rzgrid.lower[0])/geo->rzgrid.dx[0])));
  idx[1] = GKYL_MIN2(geo->rzlocal.upper[1], GKYL_MAX2(geo->rzlocal.lower[1],
    geo->rzlocal.lower[1]
      +(int) floor((Z-geo->rzgrid.lower[1])/geo->rzgrid.dx[1])));
  long loc = gkyl_range_idx(&geo->rzlocal, idx);
  const double *coeffs = gkyl_array_cfetch(geo->psiRZ, loc);
  double xc[2], eta[2];
  gkyl_rect_grid_cell_center(&geo->rzgrid, idx, xc);
  eta[0] = (R-xc[0])/(0.5*geo->rzgrid.dx[0]);
  eta[1] = (Z-xc[1])/(0.5*geo->rzgrid.dx[1]);
  return geo->rzbasis.eval_expand(eta, coeffs);
}

static void
tok_append_unique_root(double root, double *roots, int *nr, int nmax)
{
  for (int i=0; i<*nr; ++i)
    if (fabs(root-roots[i]) <= 2e-11*fmax(1.0, fabs(root)))
      return;
  if (*nr < nmax)
    roots[(*nr)++] = root;
}

// Symmetric companion to gkyl_tok_geo_R_psiZ: find every Z root at fixed R.
// The quadratic tensor representation is solved analytically in each DG cell,
// so horizontal contour segments and Z turning points are not discarded.
static int
tok_geo_Z_psiR(const struct gkyl_tok_geo *geo, double psi, double R,
  int nmax, double *Z)
{
  int nr = 0;
  if (!geo->use_cubics && geo->efit->rzbasis.poly_order <= 2) {
    int ridx = GKYL_MIN2(geo->rzlocal.upper[0],
      GKYL_MAX2(geo->rzlocal.lower[0], geo->rzlocal.lower[0]
        +(int) floor((R-geo->rzgrid.lower[0])/geo->rzgrid.dx[0])));
    for (int iz=geo->rzlocal.lower[1]; iz<=geo->rzlocal.upper[1]; ++iz) {
      int idx[2] = { ridx, iz };
      long loc = gkyl_range_idx(&geo->rzlocal, idx);
      const double *p = gkyl_array_cfetch(geo->psiRZ, loc);
      double xc[2];
      gkyl_rect_grid_cell_center(&geo->rzgrid, idx, xc);
      double x = (R-xc[0])/(0.5*geo->rzgrid.dx[0]);
      double roots[2] = { 0.0, 0.0 };
      int ncell = 0;

      if (geo->efit->rzbasis.poly_order == 1) {
        double den = 3.0*p[3]*x+1.732050807568877*p[2];
        double num = -1.732050807568877*p[1]*x+2.0*psi-p[0];
        double scale = fmax(1.0, fabs(num));
        if (fabs(den) > 64.0*DBL_EPSILON*scale)
          roots[ncell++] = num/den;
      }
      else {
        double a = 0.125*(45.0*p[8]*x*x+23.2379000772445*p[7]*x
          -15.0*p[8]+13.41640786499874*p[5]);
        double b = 0.125*(23.2379000772445*p[6]*x*x+12.0*p[3]*x
          -7.745966692414834*p[6]+6.928203230275509*p[2]);
        double c = 0.125*((13.41640786499874*p[4]-15.0*p[8])*x*x
          +(6.928203230275509*p[1]-7.745966692414834*p[7])*x
          +5.0*p[8]-4.47213595499958*p[4]-4.47213595499958*p[5]
          +4.0*p[0])-psi;
        double coeff_scale = fmax(1.0, fmax(fabs(a), fmax(fabs(b), fabs(c))));
        if (fabs(a) <= 64.0*DBL_EPSILON*coeff_scale) {
          if (fabs(b) > 64.0*DBL_EPSILON*coeff_scale)
            roots[ncell++] = -c/b;
        }
        else {
          double disc = b*b-4.0*a*c;
          double disc_tol = 256.0*DBL_EPSILON*
            fmax(1.0, b*b+fabs(4.0*a*c));
          if (disc >= -disc_tol) {
            disc = fmax(0.0, disc);
            double sd = sqrt(disc);
            if (sd == 0.0)
              roots[ncell++] = -0.5*b/a;
            else {
              double q = -0.5*(b+copysign(sd, b));
              roots[ncell++] = q/a;
              roots[ncell++] = c/q;
            }
          }
        }
      }

      for (int k=0; k<ncell; ++k) {
        double y = roots[k];
        if (isfinite(y) && y >= -1.0-2e-12 && y <= 1.0+2e-12) {
          double z = xc[1]+0.5*geo->rzgrid.dx[1]*
            fmin(1.0, fmax(-1.0, y));
          tok_append_unique_root(z, Z, &nr, nmax);
        }
      }
    }
    return nr;
  }

  // Cubic fallback: bracket sign-changing roots on a fine Z scan.  Exact
  // X-point endpoints are supplied explicitly by the trace builder.
  const int nsamp = 8*geo->rzgrid_cubic.cells[1];
  double z0 = geo->rzgrid_cubic.lower[1];
  double f0 = tok_eval_psi_rz_local(geo, R, z0)-psi;
  for (int i=1; i<=nsamp; ++i) {
    double z1 = geo->rzgrid_cubic.lower[1]
      +(geo->rzgrid_cubic.upper[1]-geo->rzgrid_cubic.lower[1])
        *i/(double) nsamp;
    double f1 = tok_eval_psi_rz_local(geo, R, z1)-psi;
    if (isfinite(f0) && isfinite(f1) && f0*f1 <= 0.0) {
      double lo = z0, hi = z1, flo = f0;
      for (int n=0; n<64; ++n) {
        double mid = 0.5*(lo+hi);
        double fm = tok_eval_psi_rz_local(geo, R, mid)-psi;
        if (flo*fm <= 0.0)
          hi = mid;
        else { lo = mid; flo = fm; }
      }
      tok_append_unique_root(0.5*(lo+hi), Z, &nr, nmax);
    }
    z0 = z1; f0 = f1;
  }
  return nr;
}

static double
tok_nearest_value(double ref, const double *values, int n)
{
  int ibest = 0;
  for (int i=1; i<n; ++i)
    if (fabs(values[i]-ref) < fabs(values[ibest]-ref))
      ibest = i;
  return values[ibest];
}

static bool tok_fixed_edge_is_midplane(enum gkyl_tok_geo_type ftype);

// Decide whether one trace step joins two points along a resolvable arc of the
// same contour, or bridges two distinct root branches, by re-parameterizing that
// single step in the *other* coordinate.  A connected arc resolves into small,
// uniform sub-steps; a branch jump leaves the sub-trace stranded on the starting
// branch and shows up as one sub-step carrying almost the whole separation.
// Returns the sub-arc length in *arc so the caller can report how much contour
// the parent step skipped.
static bool
tok_step_is_connected_arc(const struct gkyl_tok_geo *geo, double psi,
  bool param_is_r, double r0, double z0, double r1, double z1, double *arc)
{
  const int nsub = 32;
  double rp = r0, zp = z0, total = 0.0, max_sub = 0.0;
  for (int i=1; i<=nsub; ++i) {
    double t = i/(double) nsub;
    double rc, zc;
    if (i == nsub) { rc = r1; zc = z1; }
    else if (param_is_r) {
      // The parent step swept R, so sweep Z across the same interval.
      zc = z0+t*(z1-z0);
      double roots[8] = { 0.0 }, dRdZ[8] = { 0.0 };
      double dR[8] = { 0.0 }, dZ[8] = { 0.0 };
      int nr = gkyl_tok_geo_R_psiZ(geo, psi, zc, 8, roots, dRdZ, dR, dZ);
      if (nr == 0)
        return false;
      rc = tok_nearest_value(rp, roots, nr);
    }
    else {
      rc = r0+t*(r1-r0);
      double roots[16] = { 0.0 };
      int nr = tok_geo_Z_psiR(geo, psi, rc, 16, roots);
      if (nr == 0)
        return false;
      zc = tok_nearest_value(zp, roots, nr);
    }
    if (!isfinite(rc) || !isfinite(zc))
      return false;
    double residual = tok_eval_psi_rz_local(geo, rc, zc)-psi;
    if (!isfinite(residual) || fabs(residual) > 1e-9*fmax(1.0, fabs(psi)))
      return false;
    double ds = hypot(rc-rp, zc-zp);
    if (!isfinite(ds))
      return false;
    total += ds; max_sub = fmax(max_sub, ds);
    rp = rc; zp = zc;
  }
  *arc = total;
  // Uniformity of the refinement is the discriminator, not its absolute size: a
  // resolved arc spreads evenly over the sub-steps, so no sub-step may carry
  // more than 4x the average.  A branch jump concentrates the whole separation
  // into one sub-step and cannot meet this however fine the refinement.
  return isfinite(total) && total > 0.0 && max_sub <= 4.0*total/nsub;
}

static bool
tok_build_contour_candidate(const struct gkyl_tok_geo *geo, double psi,
  bool param_is_r, double rfixed, double zfixed, double rx, double zx, int n,
  double *r, double *z, double *score, const char *name,
  bool cluster_fixed_endpoint)
{
  r[0] = rfixed; z[0] = zfixed;
  double total = 0.0, max_step = 0.0, final_step = 0.0;
  double first_step = 0.0, max_step_after_first = 0.0;
  for (int i=1; i<n-1; ++i) {
    double t = i/(double) (n-1);
    // A midplane endpoint is an R turning point of the contour.  There,
    // Z-Z_fixed scales as sqrt(|R-R_fixed|), so uniform R sampling creates an
    // artificially large first physical step.  Quadratic R spacing recovers
    // approximately uniform arc-length spacing without weakening the checks
    // that reject a genuine root-branch jump.
    double f = param_is_r && cluster_fixed_endpoint ? t*t : t;
    if (param_is_r) {
      r[i] = rfixed+f*(rx-rfixed);
      double roots[16] = { 0.0 };
      int nr = tok_geo_Z_psiR(geo, psi, r[i], 16, roots);
      if (nr == 0)
        return false;
      z[i] = tok_nearest_value(z[i-1], roots, nr);
    }
    else {
      z[i] = zfixed+f*(zx-zfixed);
      double roots[8] = { 0.0 }, dRdZ[8] = { 0.0 };
      double dR[8] = { 0.0 }, dZ[8] = { 0.0 };
      int nr = gkyl_tok_geo_R_psiZ(geo, psi, z[i], 8,
        roots, dRdZ, dR, dZ);
      if (nr == 0)
        return false;
      r[i] = tok_nearest_value(r[i-1], roots, nr);
    }
    if (!isfinite(r[i]) || !isfinite(z[i]))
      return false;
    double residual = tok_eval_psi_rz_local(geo, r[i], z[i])-psi;
    if (!isfinite(residual) || fabs(residual) >
        1e-9*fmax(1.0, fabs(psi)))
      return false;
    double ds = hypot(r[i]-r[i-1], z[i]-z[i-1]);
    if (!(ds > 0.0) || !isfinite(ds))
      return false;
    total += ds; max_step = fmax(max_step, ds);
    if (i == 1) first_step = ds;
    else max_step_after_first = fmax(max_step_after_first, ds);
  }
  r[n-1] = rx; z[n-1] = zx;
  final_step = hypot(rx-r[n-2], zx-z[n-2]);
  if (!(final_step > 0.0) || !isfinite(final_step))
    return false;
  total += final_step; max_step = fmax(max_step, final_step);
  if (n > 2) max_step_after_first = fmax(max_step_after_first, final_step);
  else first_step = final_step;
  double mean = total/(n-1);
  double cell_diag = hypot(geo->rzgrid.dx[0], geo->rzgrid.dx[1]);
  // max_step/mean is a *relative* uniformity test standing in for "the trace jumped
  // to another root branch".  Arc-length spacing along a contour varies strongly near
  // an X point, so a perfectly smooth trace can exceed 16x mean while every step stays
  // far below the grid scale -- these traces were being rejected with max_step 3-20x
  // *under* the absolute threshold.  A jump to a distinct branch has to move a
  // meaningful fraction of a cell, so only let the ratio test fire once the step is
  // also large in absolute terms.  The two absolute tests are unchanged.
  const double ratio_test_floor = 0.5*cell_diag;
  // The clustered retry treats a midplane fixed edge as an R extremum of the
  // contour.  That holds unless the equilibrium was reflected about Z=0 while its
  // magnetic axis sits elsewhere (efit.c forces zmaxis=0): the mirrored pair of
  // axes then leaves a cusp at the seam, the surface has a local R *minimum*
  // there, and an R sweep inward skips the small sub-arc bulging outboard of it.
  // No amount of endpoint clustering shrinks that first step, because the arc it
  // skips is not parameterizable in R at all.  Rather than loosen the constant,
  // measure the step: re-parameterize it in Z and see whether it is one connected
  // arc.  If it is, judge uniformity on the steps that remain.  Only a trace
  // whose plain R and Z parameterizations were both already rejected is ever
  // clustered, so this cannot change a trace that is accepted today.
  double ratio_step = max_step, ratio_mean = mean;
  double skipped_arc = 0.0;
  bool first_step_verified = false;
  if (cluster_fixed_endpoint && n > 2 && first_step >= max_step &&
      tok_step_is_connected_arc(geo, psi, param_is_r, rfixed, zfixed,
        r[1], z[1], &skipped_arc)) {
    ratio_step = max_step_after_first;
    ratio_mean = (total-first_step)/(n-2);
    first_step_verified = true;
  }
  if (max_step > 2.0*cell_diag ||
      (ratio_step > 16.0*ratio_mean && ratio_step > ratio_test_floor) ||
      final_step > 2.0*cell_diag) {
    if (tok_ordered_map_diag_enabled()) {
      fprintf(stderr,
        "TOK_ORDERED_MAP_DIAG reason=discontinuous_trace trace=%s param=%c n=%d psi=%.17g max_step=%.17g final_step=%.17g mean_step=%.17g cell_diag=%.17g endpoints=(%.17g,%.17g)->(%.17g,%.17g)\n",
        name, param_is_r ? 'R' : 'Z', n, psi, max_step, final_step, mean,
        cell_diag, rfixed, zfixed, rx, zx);
      // Is the offending step an isolated spike (a real branch jump) or part of
      // a run of large steps (legitimate crowding at a turning point)?
      int imax = 1, nlarge = 0;
      double dmax = 0.0;
      for (int i=1; i<n; ++i) {
        double ds = hypot(r[i]-r[i-1], z[i]-z[i-1]);
        if (ds > dmax) { dmax = ds; imax = i; }
        if (ds > ratio_test_floor) nlarge++;
      }
      fprintf(stderr, "TOK_TRACE_STEPS trace=%s param=%c imax=%d/%d nlarge=%d neigh=",
        name, param_is_r ? 'R' : 'Z', imax, n-1, nlarge);
      for (int k=imax-3; k<=imax+3; ++k) {
        if (k < 1 || k >= n) continue;
        fprintf(stderr, "%s%.6g", k==imax ? "[" : " ",
          hypot(r[k]-r[k-1], z[k]-z[k-1]));
        if (k==imax) fprintf(stderr, "]");
      }
      fprintf(stderr, "\n");
    }
    return false;
  }
  if (first_step_verified)
    fprintf(stderr,
      "TOK_TRACE_CUSP_STEP trace=%s ftype_param=%c psi=%.17g first_step=%.17g "
      "sub_arc=%.17g excess=%.17g max_after_first=%.17g mean=%.17g "
      "fixed=(%.17g,%.17g)\n",
      name, param_is_r ? 'R' : 'Z', psi, first_step, skipped_arc,
      skipped_arc-first_step, max_step_after_first, mean, rfixed, zfixed);
  *score = max_step/mean+4.0*final_step/mean;
  return isfinite(*score);
}

static bool
tok_sep_fixed_edge_is_first(enum gkyl_tok_geo_type ftype)
{
  return ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L ||
    ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID ||
    ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_R ||
    ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO;
}

static bool tok_fixed_edge_is_midplane(enum gkyl_tok_geo_type ftype);

static void
tok_build_sep_trace(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx, double zfixed, double rfixed_ref)
{
  if (arc_ctx->sep_trace_initialized)
    return;
  if (!arc_ctx->sep_trace_r || !arc_ctx->sep_trace_z ||
      !arc_ctx->sep_trace_s || arc_ctx->sep_trace_capacity < 17) {
    fprintf(stderr, "TOK_SEP_TRACE missing storage ftype=%d\n", inp->ftype);
    abort();
  }
  double R[8] = { 0.0 }, dRdZ[8] = { 0.0 };
  double dR[8] = { 0.0 }, dZ[8] = { 0.0 };
  int nr = gkyl_tok_geo_R_psiZ(arc_ctx->geo, arc_ctx->geo->psisep,
    zfixed, 8, R, dRdZ, dR, dZ);
  if (nr == 0) {
    fprintf(stderr,
      "TOK_SEP_TRACE no fixed-edge root ftype=%d z=%.17g psi=%.17g\n",
      inp->ftype, zfixed, arc_ctx->geo->psisep);
    abort();
  }
  double rfixed = tok_nearest_value(rfixed_ref, R, nr);
  double rx = arc_ctx->geo->use_cubics
    ? arc_ctx->geo->efit->Rxpt_cubic[0] : arc_ctx->geo->efit->Rxpt[0];
  double zx = arc_ctx->geo->use_cubics
    ? arc_ctx->geo->efit->Zxpt_cubic[0] : arc_ctx->geo->efit->Zxpt[0];
  // A 257-point trace resolves a quadratic EFIT cell much more finely than
  // the geometry grid while avoiding the large startup penalty of rebuilding
  // 1025-point candidate traces in each geometry pass.
  int n = GKYL_MIN2(257, arc_ctx->sep_trace_capacity);
  double *rr = gkyl_malloc(sizeof(double[n]));
  double *zr = gkyl_malloc(sizeof(double[n]));
  double *rz = gkyl_malloc(sizeof(double[n]));
  double *zz = gkyl_malloc(sizeof(double[n]));
  double score_r = DBL_MAX, score_z = DBL_MAX;
  bool ok_r = tok_build_contour_candidate(arc_ctx->geo,
    arc_ctx->geo->psisep, true, rfixed, zfixed, rx, zx, n,
    rr, zr, &score_r, "separatrix", false);
  bool ok_z = tok_build_contour_candidate(arc_ctx->geo,
    arc_ctx->geo->psisep, false, rfixed, zfixed, rx, zx, n,
    rz, zz, &score_z, "separatrix", false);
  // Same endpoint-clustered R retry the far-boundary trace already uses. A
  // midplane fixed edge is an R turning point of the contour: Z-Z_fixed scales
  // as sqrt(|R-R_fixed|) there, so uniform R sampling makes the *first* physical
  // step artificially large and trips the discontinuity test even though the
  // contour is perfectly smooth (204502: max_step is step 1 of 256, with
  // neighbours decaying 0.086/0.043/0.026/0.019 -- a turning point, not a branch
  // jump, which would be an isolated spike like this shot's param=Z at 249/256).
  // Only a trace for which both parameterizations were already rejected gets the
  // retry, so every accepted trace is unchanged.
  if (!ok_r && !ok_z && tok_fixed_edge_is_midplane(inp->ftype))
    ok_r = tok_build_contour_candidate(arc_ctx->geo,
      arc_ctx->geo->psisep, true, rfixed, zfixed, rx, zx, n,
      rr, zr, &score_r, "separatrix", true);
  if (!ok_r && !ok_z) {
    fprintf(stderr,
      "TOK_SEP_TRACE both parameterizations failed ftype=%d fixed=(%.17g,%.17g) xpt=(%.17g,%.17g)\n",
      inp->ftype, rfixed, zfixed, rx, zx);
    abort();
  }
  bool use_r = ok_r && (!ok_z || score_r <= score_z);
  const double *src_r = use_r ? rr : rz;
  const double *src_z = use_r ? zr : zz;
  bool fixed_first = tok_sep_fixed_edge_is_first(inp->ftype);
  for (int i=0; i<n; ++i) {
    int src = fixed_first ? i : n-1-i;
    arc_ctx->sep_trace_r[i] = src_r[src];
    arc_ctx->sep_trace_z[i] = src_z[src];
  }
  arc_ctx->sep_trace_s[0] = 0.0;
  for (int i=1; i<n; ++i)
    arc_ctx->sep_trace_s[i] = arc_ctx->sep_trace_s[i-1]
      +hypot(arc_ctx->sep_trace_r[i]-arc_ctx->sep_trace_r[i-1],
        arc_ctx->sep_trace_z[i]-arc_ctx->sep_trace_z[i-1]);
  if (!(arc_ctx->sep_trace_s[n-1] > 0.0) ||
      !isfinite(arc_ctx->sep_trace_s[n-1])) {
    fprintf(stderr, "TOK_SEP_TRACE invalid arc length ftype=%d\n", inp->ftype);
    abort();
  }
  arc_ctx->sep_trace_n = n;
  arc_ctx->sep_trace_param_is_r = use_r;
  arc_ctx->sep_trace_initialized = true;
  gkyl_free(rr); gkyl_free(zr); gkyl_free(rz); gkyl_free(zz);
}

static bool
tok_half_domain_sep_rz(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx, double theta, double *r, double *z)
{
  if (!inp->half_domain ||
      !tok_geo_same_flux(arc_ctx->psi, arc_ctx->geo->psisep))
    return false;

  double zx = arc_ctx->geo->use_cubics
    ? arc_ctx->geo->efit->Zxpt_cubic[0] : arc_ctx->geo->efit->Zxpt[0];
  double za = arc_ctx->geo->zmaxis, z0 = 0.0, z1 = 0.0, rclose = 0.0;
  switch (inp->ftype) {
    case GKYL_GEOMETRY_TOKAMAK_CORE_R:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID:
      z0 = zx; z1 = za; rclose = inp->rright;
      break;
    case GKYL_GEOMETRY_TOKAMAK_CORE_L:
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID:
      z0 = za; z1 = zx; rclose = inp->rleft;
      break;
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_R:
      z0 = arc_ctx->zmin_right; z1 = zx; rclose = inp->rright;
      break;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO:
      z0 = arc_ctx->zmin; z1 = zx; rclose = inp->rright;
      break;
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_L:
      z0 = zx; z1 = arc_ctx->zmin_left; rclose = inp->rleft;
      break;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO:
      z0 = zx; z1 = arc_ctx->zmin; rclose = inp->rleft;
      break;
    default:
      return false;
  }

  // Use the actual (position-mapped) computational coordinate.  Corner
  // nodes lie at frac=0 or 1, while interior and radial/alpha-surface
  // quadrature nodes do not.  Inferring frac from the nodal-array index
  // incorrectly maps the first and last Gauss points onto the segment ends.
  double frac = (theta-inp->cgrid.lower[2])
    /(inp->cgrid.upper[2]-inp->cgrid.lower[2]);
  double zfixed = tok_sep_fixed_edge_is_first(inp->ftype) ? z0 : z1;
  tok_build_sep_trace(inp, arc_ctx, zfixed, rclose);
  int n = arc_ctx->sep_trace_n;
  double total = arc_ctx->sep_trace_s[n-1];
  arc_ctx->xpt_map_darc_dtheta = total/
    (inp->cgrid.upper[2]-inp->cgrid.lower[2]);
  if (frac <= 0.0) {
    *r = arc_ctx->sep_trace_r[0]; *z = arc_ctx->sep_trace_z[0];
    return true;
  }
  if (frac >= 1.0) {
    *r = arc_ctx->sep_trace_r[n-1]; *z = arc_ctx->sep_trace_z[n-1];
    return true;
  }
  double target = frac*total;
  int lo = 0, hi = n-1;
  while (hi-lo > 1) {
    int mid = (lo+hi)/2;
    if (arc_ctx->sep_trace_s[mid] < target) lo = mid;
    else hi = mid;
  }
  double ds = arc_ctx->sep_trace_s[hi]-arc_ctx->sep_trace_s[lo];
  double w = ds > 0.0 ? (target-arc_ctx->sep_trace_s[lo])/ds : 0.0;
  double rlin = arc_ctx->sep_trace_r[lo]
    +w*(arc_ctx->sep_trace_r[hi]-arc_ctx->sep_trace_r[lo]);
  double zlin = arc_ctx->sep_trace_z[lo]
    +w*(arc_ctx->sep_trace_z[hi]-arc_ctx->sep_trace_z[lo]);
  if (arc_ctx->sep_trace_param_is_r) {
    double roots[16] = { 0.0 };
    int nr = tok_geo_Z_psiR(arc_ctx->geo, arc_ctx->geo->psisep,
      rlin, 16, roots);
    if (nr == 0) {
      fprintf(stderr,
        "TOK_SEP_TRACE no Z root during resampling ftype=%d frac=%.17g R=%.17g\n",
        inp->ftype, frac, rlin);
      abort();
    }
    *r = rlin; *z = tok_nearest_value(zlin, roots, nr);
  }
  else {
    double roots[8] = { 0.0 }, dRdZ[8] = { 0.0 };
    double dR[8] = { 0.0 }, dZ[8] = { 0.0 };
    int nr = gkyl_tok_geo_R_psiZ(arc_ctx->geo, arc_ctx->geo->psisep,
      zlin, 8, roots, dRdZ, dR, dZ);
    if (nr == 0) {
      fprintf(stderr,
        "TOK_SEP_TRACE no R root during resampling ftype=%d frac=%.17g Z=%.17g\n",
        inp->ftype, frac, zlin);
      abort();
    }
    *r = tok_nearest_value(rlin, roots, nr); *z = zlin;
  }
  double residual = tok_eval_psi_rz_local(arc_ctx->geo, *r, *z)
    -arc_ctx->geo->psisep;
  if (!isfinite(*r) || !isfinite(*z) || !isfinite(residual) ||
      fabs(residual) > 1e-9*fmax(1.0, fabs(arc_ctx->geo->psisep))) {
    fprintf(stderr,
      "TOK_SEP_TRACE invalid resampled point ftype=%d frac=%.17g R=%.17g Z=%.17g residual=%.17g\n",
      inp->ftype, frac, *r, *z, residual);
    abort();
  }
  return true;
}

struct tok_ordered_point {
  double r, z, phi;
  double dr_dtheta, dz_dtheta, dphi_dtheta;
};

static bool
tok_fixed_edge_is_midplane(enum gkyl_tok_geo_type ftype)
{
  return ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R ||
    ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L ||
    ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID ||
    ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID;
}

static bool
tok_fixed_edge_uses_lower_plate(enum gkyl_tok_geo_type ftype)
{
  return ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_R ||
    ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO ||
    ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO;
}

// Which side of the machine a plate-anchored block's strike point sits on.  The two
// halves of a private-flux region share one plate, so a flux surface crossing it twice
// yields one root per half; this picks the half that belongs to this block.
static bool
tok_plate_edge_is_outboard(enum gkyl_tok_geo_type ftype)
{
  return ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_R ||
    ftype == GKYL_GEOMETRY_TOKAMAK_PF_UP_R ||
    ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO ||
    ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO;
}

static bool
tok_plate_flux_intersection(const struct gkyl_tok_geo *geo,
  plate_func plate, double psi, int ftype, double *r, double *z)
{
  if (!plate)
    return false;

  const int nsamp = 512;
  double roots_s[16] = { 0.0 }, roots_r[16] = { 0.0 };
  double roots_z[16] = { 0.0 };
  int nroots = 0;
  double rz0[2] = { 0.0 }, rz1[2] = { 0.0 };
  plate(0.0, rz0);
  double f0 = tok_eval_psi_rz_local(geo, rz0[0], rz0[1])-psi;
  double flux_tol = 1e-10*fmax(1.0, fabs(psi));
  if (isfinite(f0) && fabs(f0) <= flux_tol) {
    roots_s[nroots] = 0.0; roots_r[nroots] = rz0[0];
    roots_z[nroots++] = rz0[1];
  }
  for (int i=1; i<=nsamp; ++i) {
    double s1 = i/(double) nsamp;
    plate(s1, rz1);
    double f1 = tok_eval_psi_rz_local(geo, rz1[0], rz1[1])-psi;
    bool endpoint_root = isfinite(f1) && fabs(f1) <= flux_tol;
    bool bracket = isfinite(f0) && isfinite(f1) && f0*f1 < 0.0;
    if ((endpoint_root || bracket) && nroots < 16) {
      double sr = s1;
      if (!endpoint_root) {
        double slo = (i-1)/(double) nsamp, shi = s1, flo = f0;
        for (int k=0; k<64; ++k) {
          double smid = 0.5*(slo+shi), rzm[2];
          plate(smid, rzm);
          double fm = tok_eval_psi_rz_local(geo, rzm[0], rzm[1])-psi;
          if (!isfinite(fm))
            return false;
          if (flo*fm <= 0.0)
            shi = smid;
          else { slo = smid; flo = fm; }
        }
        sr = 0.5*(slo+shi);
      }
      double rzr[2];
      plate(sr, rzr);
      if (nroots == 0 || fabs(sr-roots_s[nroots-1]) > 2e-10) {
        roots_s[nroots] = sr; roots_r[nroots] = rzr[0];
        roots_z[nroots++] = rzr[1];
      }
    }
    rz0[0] = rz1[0]; rz0[1] = rz1[1]; f0 = f1;
  }
  if (nroots == 0) {
    fprintf(stderr,
      "TOK_ORDERED_MAP plate root count=0 ftype=%d psi=%.17g\n", ftype, psi);
    return false;
  }
  int pick = 0;
  if (nroots > 1) {
    // The two halves of a private-flux region share this plate, so a surface that
    // crosses it twice gives one root per half.  Select the half this block owns
    // using the same side hint the midplane path uses, rather than rejecting.
    double hint = tok_plate_edge_is_outboard(ftype) ? geo->rright : geo->rleft;
    double best = fabs(roots_r[0]-hint);
    for (int k=1; k<nroots; ++k) {
      double d = fabs(roots_r[k]-hint);
      if (d < best) { best = d; pick = k; }
    }
    fprintf(stderr,
      "TOK_ORDERED_MAP plate root count=%d ftype=%d psi=%.17g hint=%.17g picked=%d rz=(%.17g,%.17g)\n",
      nroots, ftype, psi, hint, pick, roots_r[pick], roots_z[pick]);
  }
  *r = roots_r[pick]; *z = roots_z[pick];
  double residual = tok_eval_psi_rz_local(geo, *r, *z)-psi;
  return isfinite(*r) && isfinite(*z) && isfinite(residual) &&
    fabs(residual) <= 1e-9*fmax(1.0, fabs(psi));
}

// Locate the physical (plate or midplane) end of a block at a specified psi.
// Plate intersections are scanned before refinement so an exact endpoint is
// returned directly and an ambiguous multi-intersection plate is rejected.
static bool
tok_fixed_edge_point(const struct gkyl_tok_geo_grid_inp *inp,
  const struct gkyl_tok_geo *geo, double psi, double *r, double *z)
{
  if (tok_fixed_edge_is_midplane(inp->ftype)) {
    double R[16] = { 0.0 }, dRdZ[16] = { 0.0 };
    double dR[16] = { 0.0 }, dZ[16] = { 0.0 };
    *z = geo->zmaxis;
    int nr = gkyl_tok_geo_R_psiZ(geo, psi, *z, 8, R, dRdZ, dR, dZ);
    if (nr <= 0)
      return false;
    bool outboard = inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R ||
      inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID;
    *r = tok_nearest_value(outboard ? inp->rright : inp->rleft, R, nr);
    return isfinite(*r);
  }

  plate_func plate = 0;
  if (tok_fixed_edge_uses_lower_plate(inp->ftype))
    plate = geo->plate_func_lower;
  else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_L)
    plate = geo->plate_func_upper;
  return tok_plate_flux_intersection(geo, plate, psi, inp->ftype, r, z);
}

enum tok_ext_endpoint_kind {
  TOK_EXT_XPT_RAY,
  TOK_EXT_PLATE,
  // Half-domain MID/CORE blocks stop at the midplane symmetry plane rather
  // than at the far X point, so their topology needs an endpoint that is
  // neither a ray nor a material surface: the point where this psi contour
  // meets Z = zmaxis on the given side.
  TOK_EXT_MIDPLANE,
};

enum tok_ext_xpoint {
  TOK_EXT_LOWER_XPT,
  TOK_EXT_UPPER_XPT,
};

enum tok_ext_sector {
  TOK_EXT_CORE,
  TOK_EXT_PF,
  TOK_EXT_SOL_OUT,
  TOK_EXT_SOL_IN,
};

enum tok_ext_plate_slot {
  TOK_EXT_PLATE_LOWER,
  TOK_EXT_PLATE_UPPER,
};

enum tok_ext_fixed_z_slot {
  TOK_EXT_ZMIN,
  TOK_EXT_ZMAX,
  TOK_EXT_ZMIN_LEFT,
  TOK_EXT_ZMIN_RIGHT,
  TOK_EXT_ZMAX_LEFT,
  TOK_EXT_ZMAX_RIGHT,
};

enum tok_ext_route {
  TOK_EXT_ROUTE_GENERIC,
  TOK_EXT_ROUTE_OUTBOARD,
  TOK_EXT_ROUTE_INBOARD,
  TOK_EXT_ROUTE_VIA_UPPER,
  TOK_EXT_ROUTE_CLOSED_CORE,
};

enum tok_ext_phi_reference {
  TOK_EXT_PHI_LOWER,
  TOK_EXT_PHI_UPPER,
  TOK_EXT_PHI_OUTBOARD_MIDPLANE,
  TOK_EXT_PHI_INBOARD_MIDPLANE,
};

struct tok_ext_endpoint {
  enum tok_ext_endpoint_kind kind;
  enum tok_ext_xpoint xpoint;
  enum tok_ext_sector sector;
  enum tok_ext_plate_slot plate_slot;
  enum tok_ext_fixed_z_slot fixed_z_slot;
  // TOK_EXT_MIDPLANE only. CORE_R and CORE_L share a sector but sit on
  // opposite sides, so the side cannot be derived from `sector`.
  bool midplane_outboard;
};

struct tok_ext_topology {
  struct tok_ext_endpoint lower, upper;
  enum tok_ext_route route;
  enum tok_ext_phi_reference phi_reference;
  bool closed;
};

static struct tok_ext_endpoint
tok_ext_xray(enum tok_ext_xpoint xpoint, enum tok_ext_sector sector)
{
  return (struct tok_ext_endpoint) {
    .kind = TOK_EXT_XPT_RAY, .xpoint = xpoint, .sector = sector,
  };
}

static struct tok_ext_endpoint
tok_ext_plate(enum tok_ext_plate_slot plate_slot,
  enum tok_ext_fixed_z_slot fixed_z_slot)
{
  return (struct tok_ext_endpoint) {
    .kind = TOK_EXT_PLATE, .plate_slot = plate_slot,
    .fixed_z_slot = fixed_z_slot,
  };
}

static struct tok_ext_endpoint
tok_ext_midplane(bool outboard)
{
  return (struct tok_ext_endpoint) {
    .kind = TOK_EXT_MIDPLANE, .midplane_outboard = outboard,
  };
}

// Logical theta-lower -> theta-upper topology for the double-null and
// lower-single-null multiblock types.
//
// `half_domain` models only the lower half, closing the domain on the midplane
// symmetry plane. That changes the topology of exactly four blocks: the two
// MID SOL blocks and the two CORE halves stop at the midplane instead of
// running on to the far X point (see tok_geo_set_extent, where those four are
// the only ftypes whose theta extents carry a `half_domain` branch). The four
// leg blocks -- PF_LO_R, PF_LO_L, DN_SOL_OUT_LO, DN_SOL_IN_LO -- have
// identical extents in both modes and so share a single entry here.
static bool
tok_ext_topology_from_ftype(enum gkyl_tok_geo_type ftype, bool half_domain,
  struct tok_ext_topology *top)
{
  *top = (struct tok_ext_topology) { };
  if (half_domain) {
    switch (ftype) {
      case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID:
        top->lower = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_SOL_OUT);
        top->upper = tok_ext_midplane(true);
        top->route = TOK_EXT_ROUTE_OUTBOARD;
        top->phi_reference = TOK_EXT_PHI_OUTBOARD_MIDPLANE;
        return true;
      case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID:
        top->lower = tok_ext_midplane(false);
        top->upper = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_SOL_IN);
        top->route = TOK_EXT_ROUTE_INBOARD;
        top->phi_reference = TOK_EXT_PHI_INBOARD_MIDPLANE;
        return true;
      // The core X-point ray endpoint on an INTERIOR surface is the nearest
      // point on that surface to the X point, which lies on its inboard-lower
      // flank rather than at its extreme Z. CORE_R therefore runs inboard,
      // under the surface's LOWER turning point, and up to the outboard
      // midplane -- a path where no single parameterization works: Z reverses
      // at the turning point and R is degenerate at the midplane end. Split it
      // gate by 6e-5 relative -- so give the generic route the endpoint-
      // clustered retry (see tok_ext_build_open_trace) and it gets through.
      //
      // DO NOT RETRY the obvious-looking alternative of splitting at the lower
      // turning point (a mirrored tok_ext_build_via_turning_trace). Measured on
      // all 20 shots that need this: 15 built but ALL 15 folded CORE_R by 8-10
      // cells with minA/medA ~ -6e3, versus 1-3 cells at ~-1e-1 for the
      // clustered generic route. Uniform-Z sampling crowds the arc against the
      // turning point, and resampling the joined polyline by arc length just
      // interpolates across that chord instead of following the contour.
      case GKYL_GEOMETRY_TOKAMAK_CORE_R:
        top->lower = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_CORE);
        top->upper = tok_ext_midplane(true);
        top->route = TOK_EXT_ROUTE_GENERIC;
        top->phi_reference = TOK_EXT_PHI_OUTBOARD_MIDPLANE;
        return true;
      case GKYL_GEOMETRY_TOKAMAK_CORE_L:
        top->lower = tok_ext_midplane(false);
        top->upper = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_CORE);
        top->route = TOK_EXT_ROUTE_GENERIC;
        top->phi_reference = TOK_EXT_PHI_INBOARD_MIDPLANE;
        return true;
      default:
        break;  // leg blocks fall through to the shared table below
    }
  }
  switch (ftype) {
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_R:
      top->lower = tok_ext_plate(TOK_EXT_PLATE_LOWER,
        TOK_EXT_ZMIN_RIGHT);
      top->upper = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_PF);
      top->route = TOK_EXT_ROUTE_GENERIC;
      top->phi_reference = TOK_EXT_PHI_LOWER;
      return true;
    case GKYL_GEOMETRY_TOKAMAK_PF_LO_L:
      top->lower = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_PF);
      top->upper = tok_ext_plate(TOK_EXT_PLATE_UPPER,
        TOK_EXT_ZMIN_LEFT);
      top->route = TOK_EXT_ROUTE_GENERIC;
      top->phi_reference = TOK_EXT_PHI_UPPER;
      return true;
    case GKYL_GEOMETRY_TOKAMAK_PF_UP_R:
      top->lower = tok_ext_xray(TOK_EXT_UPPER_XPT, TOK_EXT_PF);
      top->upper = tok_ext_plate(TOK_EXT_PLATE_UPPER,
        TOK_EXT_ZMAX_RIGHT);
      top->route = TOK_EXT_ROUTE_GENERIC;
      top->phi_reference = TOK_EXT_PHI_UPPER;
      return true;
    case GKYL_GEOMETRY_TOKAMAK_PF_UP_L:
      top->lower = tok_ext_plate(TOK_EXT_PLATE_LOWER,
        TOK_EXT_ZMAX_LEFT);
      top->upper = tok_ext_xray(TOK_EXT_UPPER_XPT, TOK_EXT_PF);
      top->route = TOK_EXT_ROUTE_GENERIC;
      top->phi_reference = TOK_EXT_PHI_LOWER;
      return true;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO:
      top->lower = tok_ext_plate(TOK_EXT_PLATE_LOWER, TOK_EXT_ZMIN);
      top->upper = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_SOL_OUT);
      top->route = TOK_EXT_ROUTE_GENERIC;
      top->phi_reference = TOK_EXT_PHI_LOWER;
      return true;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID:
      top->lower = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_SOL_OUT);
      top->upper = tok_ext_xray(TOK_EXT_UPPER_XPT, TOK_EXT_SOL_OUT);
      top->route = TOK_EXT_ROUTE_OUTBOARD;
      top->phi_reference = TOK_EXT_PHI_OUTBOARD_MIDPLANE;
      return true;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP:
      top->lower = tok_ext_xray(TOK_EXT_UPPER_XPT, TOK_EXT_SOL_OUT);
      top->upper = tok_ext_plate(TOK_EXT_PLATE_UPPER, TOK_EXT_ZMAX);
      top->route = TOK_EXT_ROUTE_GENERIC;
      top->phi_reference = TOK_EXT_PHI_UPPER;
      return true;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP:
      top->lower = tok_ext_plate(TOK_EXT_PLATE_UPPER, TOK_EXT_ZMAX);
      top->upper = tok_ext_xray(TOK_EXT_UPPER_XPT, TOK_EXT_SOL_IN);
      top->route = TOK_EXT_ROUTE_GENERIC;
      top->phi_reference = TOK_EXT_PHI_LOWER;
      return true;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID:
      top->lower = tok_ext_xray(TOK_EXT_UPPER_XPT, TOK_EXT_SOL_IN);
      top->upper = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_SOL_IN);
      top->route = TOK_EXT_ROUTE_INBOARD;
      top->phi_reference = TOK_EXT_PHI_INBOARD_MIDPLANE;
      return true;
    case GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO:
      top->lower = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_SOL_IN);
      top->upper = tok_ext_plate(TOK_EXT_PLATE_LOWER, TOK_EXT_ZMIN);
      top->route = TOK_EXT_ROUTE_GENERIC;
      top->phi_reference = TOK_EXT_PHI_UPPER;
      return true;
    case GKYL_GEOMETRY_TOKAMAK_CORE_R:
      top->lower = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_CORE);
      top->upper = tok_ext_xray(TOK_EXT_UPPER_XPT, TOK_EXT_CORE);
      top->route = TOK_EXT_ROUTE_OUTBOARD;
      top->phi_reference = TOK_EXT_PHI_OUTBOARD_MIDPLANE;
      return true;
    case GKYL_GEOMETRY_TOKAMAK_CORE_L:
      top->lower = tok_ext_xray(TOK_EXT_UPPER_XPT, TOK_EXT_CORE);
      top->upper = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_CORE);
      top->route = TOK_EXT_ROUTE_INBOARD;
      top->phi_reference = TOK_EXT_PHI_INBOARD_MIDPLANE;
      return true;
    case GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO:
      top->lower = tok_ext_plate(TOK_EXT_PLATE_LOWER,
        TOK_EXT_ZMIN_RIGHT);
      top->upper = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_SOL_OUT);
      top->route = TOK_EXT_ROUTE_GENERIC;
      top->phi_reference = TOK_EXT_PHI_LOWER;
      return true;
    case GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID:
      top->lower = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_SOL_OUT);
      top->upper = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_SOL_IN);
      top->route = TOK_EXT_ROUTE_VIA_UPPER;
      top->phi_reference = TOK_EXT_PHI_OUTBOARD_MIDPLANE;
      return true;
    case GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP:
      top->lower = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_SOL_IN);
      top->upper = tok_ext_plate(TOK_EXT_PLATE_UPPER,
        TOK_EXT_ZMIN_LEFT);
      top->route = TOK_EXT_ROUTE_GENERIC;
      top->phi_reference = TOK_EXT_PHI_UPPER;
      return true;
    case GKYL_GEOMETRY_TOKAMAK_CORE:
      top->lower = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_CORE);
      top->upper = tok_ext_xray(TOK_EXT_LOWER_XPT, TOK_EXT_CORE);
      top->route = TOK_EXT_ROUTE_CLOSED_CORE;
      top->phi_reference = TOK_EXT_PHI_OUTBOARD_MIDPLANE;
      top->closed = true;
      return true;
    default:
      return false;
  }
}

static bool
tok_ext_xpoint_rz(const struct gkyl_tok_geo *geo,
  enum tok_ext_xpoint which, double *r, double *z)
{
  int nx = geo->use_cubics
    ? geo->efit->num_xpts_cubic : geo->efit->num_xpts;
  const double *rx = geo->use_cubics
    ? geo->efit->Rxpt_cubic : geo->efit->Rxpt;
  const double *zx = geo->use_cubics
    ? geo->efit->Zxpt_cubic : geo->efit->Zxpt;
  if (nx < 1 || (which == TOK_EXT_UPPER_XPT && nx < 2))
    return false;
  if (nx == 1) {
    *r = rx[0]; *z = zx[0];
    return isfinite(*r) && isfinite(*z);
  }
  int ilo = zx[0] <= zx[1] ? 0 : 1;
  int iup = ilo == 0 ? 1 : 0;
  int idx = which == TOK_EXT_LOWER_XPT ? ilo : iup;
  *r = rx[idx]; *z = zx[idx];
  return isfinite(*r) && isfinite(*z);
}

static bool
tok_ext_nearest_root_at_z(const struct gkyl_tok_geo_grid_inp *inp,
  const struct gkyl_tok_geo *geo, enum tok_ext_sector sector, double psi,
  double z, double rx, double zx, double *rbest, double *d2best)
{
  double R[8] = { 0.0 }, dRdZ[8] = { 0.0 };
  double dR[8] = { 0.0 }, dZ[8] = { 0.0 };
  int nr = gkyl_tok_geo_R_psiZ(geo, psi, z, 8, R, dRdZ, dR, dZ);
  if (nr <= 0)
    return false;

  int first = 0, last = nr;
  if (sector == TOK_EXT_SOL_OUT || sector == TOK_EXT_SOL_IN) {
    double ref = sector == TOK_EXT_SOL_OUT ? inp->rright : inp->rleft;
    int selected = 0;
    for (int i=1; i<nr; ++i)
      if (fabs(R[i]-ref) < fabs(R[selected]-ref))
        selected = i;
    first = selected; last = selected+1;
  }

  bool found = false;
  *d2best = DBL_MAX;
  for (int i=first; i<last; ++i) {
    double d2 = SQ(R[i]-rx)+SQ(z-zx);
    if (d2 < *d2best) {
      *d2best = d2; *rbest = R[i]; found = true;
    }
  }
  return found;
}

static bool
tok_ext_nearest_ray_target(const struct gkyl_tok_geo_grid_inp *inp,
  const struct gkyl_tok_geo *geo, const struct tok_ext_endpoint *endpoint,
  double psi, double *rnear, double *znear)
{
  double rx = 0.0, zx = 0.0;
  if (!tok_ext_xpoint_rz(geo, endpoint->xpoint, &rx, &zx))
    return false;
  double zlo = geo->rzgrid.lower[1], zhi = geo->rzgrid.upper[1];
  int nx = geo->use_cubics
    ? geo->efit->num_xpts_cubic : geo->efit->num_xpts;
  // A one-X-point SOL has two cuts at the same saddle.  Its unconstrained
  // nearest points can move behind a shaped divertor plate as psi changes,
  // collapsing a leg.  Anchor the two far-surface rays at the outboard and
  // inboard magnetic-midplane roots instead; these cuts preserve the full-SOL
  // topology and are shared exactly by the adjacent blocks.
  if (nx == 1 && (endpoint->sector == TOK_EXT_SOL_OUT ||
      endpoint->sector == TOK_EXT_SOL_IN)) {
    double R[16] = { 0.0 }, dRdZ[16] = { 0.0 };
    double dR[16] = { 0.0 }, dZ[16] = { 0.0 };
    *znear = geo->zmaxis;
    int nr = gkyl_tok_geo_R_psiZ(geo, psi, *znear, 16,
      R, dRdZ, dR, dZ);
    if (nr <= 0)
      return false;
    bool outboard = endpoint->sector == TOK_EXT_SOL_OUT;
    *rnear = tok_nearest_value(outboard ? inp->rright : inp->rleft, R, nr);
    double residual = tok_eval_psi_rz_local(geo, *rnear, *znear)-psi;
    return isfinite(*rnear) && isfinite(*znear) && isfinite(residual) &&
      fabs(residual) <= 1e-9*fmax(1.0, fabs(psi));
  }
  if (endpoint->sector == TOK_EXT_CORE) {
    zlo = fmin(zx, geo->zmaxis); zhi = fmax(zx, geo->zmaxis);
  }
  else if (endpoint->sector == TOK_EXT_PF) {
    if (endpoint->xpoint == TOK_EXT_LOWER_XPT) {
      zlo = geo->rzgrid.lower[1]; zhi = zx;
    }
    else {
      zlo = zx; zhi = geo->rzgrid.upper[1];
    }
  }
  else if (endpoint->xpoint == TOK_EXT_LOWER_XPT) {
    zlo = geo->rzgrid.lower[1]; zhi = geo->zmaxis;
  }
  else {
    zlo = geo->zmaxis; zhi = geo->rzgrid.upper[1];
  }

  if (!(zhi > zlo))
    return false;

  int nsamp = 4*geo->rzgrid.cells[1]+1;
  double best_d2 = DBL_MAX, best_r = 0.0, best_z = 0.0;
  int best_i = -1;
  for (int i=0; i<nsamp; ++i) {
    double z = zlo+(zhi-zlo)*i/(nsamp-1.0), r = 0.0, d2 = 0.0;
    if (tok_ext_nearest_root_at_z(inp, geo, endpoint->sector, psi,
        z, rx, zx, &r, &d2) && d2 < best_d2) {
      best_d2 = d2; best_r = r; best_z = z; best_i = i;
    }
  }
  if (best_i < 0)
    return false;

  double dz = (zhi-zlo)/(nsamp-1.0);
  double a = fmax(zlo, best_z-dz), b = fmin(zhi, best_z+dz);
  const double gr = 0.6180339887498948482;
  for (int iter=0; iter<64; ++iter) {
    double ztrial[2] = { b-gr*(b-a), a+gr*(b-a) };
    double dtrial[2] = { DBL_MAX, DBL_MAX };
    double rtrial[2] = { best_r, best_r };
    for (int k=0; k<2; ++k)
      tok_ext_nearest_root_at_z(inp, geo, endpoint->sector, psi,
        ztrial[k], rx, zx, &rtrial[k], &dtrial[k]);
    if (dtrial[0] < best_d2) {
      best_d2 = dtrial[0]; best_r = rtrial[0]; best_z = ztrial[0];
    }
    if (dtrial[1] < best_d2) {
      best_d2 = dtrial[1]; best_r = rtrial[1]; best_z = ztrial[1];
    }
    if (dtrial[0] <= dtrial[1]) b = ztrial[1];
    else a = ztrial[0];
  }
  *rnear = best_r; *znear = best_z;
  double residual = tok_eval_psi_rz_local(geo, best_r, best_z)-psi;
  return isfinite(best_r) && isfinite(best_z) && isfinite(residual) &&
    fabs(residual) <= 1e-9*fmax(1.0, fabs(psi));
}

// Intersect the current surface with the fixed line from the selected X point
// to its topology-safe target on the far radial surface.  Counting crossings
// before bisection prevents a locally reversed DG flux representation from
// silently selecting a different branch of the ray.
static bool
tok_ext_fixed_ray_endpoint(const struct gkyl_tok_geo_grid_inp *inp,
  const struct arc_length_ctx *arc_ctx,
  const struct tok_ext_endpoint *endpoint, double psi,
  double *r, double *z)
{
  const struct gkyl_tok_geo *geo = arc_ctx->geo;
  double rx = 0.0, zx = 0.0, rf = 0.0, zf = 0.0;
  if (!tok_ext_xpoint_rz(geo, endpoint->xpoint, &rx, &zx) ||
      !tok_ext_nearest_ray_target(inp, geo, endpoint,
        arc_ctx->xpt_ray_psi0, &rf, &zf))
    return false;

  double scale = fmax(1.0, fmax(fabs(psi), fabs(geo->psisep)));
  if (tok_geo_same_flux(psi, geo->psisep)) {
    *r = rx; *z = zx;
    return true;
  }
  if (tok_geo_same_flux(psi, arc_ctx->xpt_ray_psi0)) {
    *r = rf; *z = zf;
    return true;
  }
  double delta = arc_ctx->xpt_ray_psi0-geo->psisep;
  if (!isfinite(delta) || fabs(delta) <= 256.0*DBL_EPSILON*scale)
    return false;
  double qtarget = (psi-geo->psisep)/delta;
  if (qtarget <= 0.0 || qtarget >= 1.0)
    return false;

  const bool use_last = arc_ctx->ext_ray_use_last_crossing ||
    tok_ext_ray_last_crossing();
  const int nsamp = 512;
  const double crossing_tol = 1e-13*fmax(1.0, fabs(qtarget));
  double qprev = (tok_eval_psi_rz_local(geo, rx, zx)-geo->psisep)/delta;
  double sprev = 0.0, slo = -1.0, shi = -1.0;
  int upward = 0, downward = 0;
  // Lowest q seen before the ray first reaches the target surface.
  double qmin_before_first = 0.0;
  // Crossing census, for telling a spurious level-set island closed off by the
  // C0 grad-psi ridge on a DG cell face from a genuine fold-back of the ray.
  const int max_events = 16;
  double event_s[16];
  bool event_up[16];
  int n_events = 0;
  double qmax = qprev, qend = qprev;
  for (int i=1; i<=nsamp; ++i) {
    double s = i/(double) nsamp;
    double rs = rx+s*(rf-rx), zs = zx+s*(zf-zx);
    double q = (tok_eval_psi_rz_local(geo, rs, zs)-geo->psisep)/delta;
    if (!isfinite(q))
      return false;
    bool crosses_up = qprev < qtarget-crossing_tol &&
      q >= qtarget-crossing_tol;
    bool crosses_down = qprev > qtarget+crossing_tol &&
      q <= qtarget+crossing_tol;
    if (crosses_up) {
      ++upward;
      // q(0)=0 at the X point and q(1)=1 at the far ray target by construction,
      // so the last upward crossing is the point beyond which the ray never
      // re-enters the block.  Taking it is the running-minimum-from-the-right
      // of q, which is monotone and meets qtarget exactly once, so the anchor
      // stays continuous in psi even when a spurious level-set island sits
      // between the X point and the real surface.
      if (use_last || slo < 0.0) { slo = sprev; shi = s; }
    }
    if (crosses_down)
      ++downward;
    if ((crosses_up || crosses_down) && n_events < max_events) {
      event_s[n_events] = 0.5*(sprev+s);
      event_up[n_events] = crosses_up;
      ++n_events;
    }
    if (upward == 0)
      qmin_before_first = fmin(qmin_before_first, q);
    qmax = fmax(qmax, q);
    qend = q;
    qprev = q; sprev = s;
  }
  if (upward > 1) {
    fprintf(stderr,
      "TOK_EXT_RAY_CROSSINGS ftype=%d sector=%d psi=%.17g qtarget=%.17g "
      "upward=%d downward=%d qmax=%.6e qend=%.6e raylen=%.6e events=",
      inp->ftype, endpoint->sector, psi, qtarget, upward, downward,
      qmax, qend, hypot(rf-rx, zf-zx));
    for (int k=0; k<n_events; ++k)
      fprintf(stderr, "%s%c:%.6f", k ? "," : "", event_up[k] ? 'u' : 'd',
        event_s[k]);
    fprintf(stderr, "\n");
  }
  // A ray that folds back crosses psi=psi_curr more than once.  The crossings
  // are not interchangeable: only the FIRST tends to the X point as
  // psi_curr -> psisep, so it is the one that keeps the anchor continuous in
  // psi.  Accept it when the ray got there without dipping inside the
  // separatrix.  This is the same rule the chord construction's
  // tok_xpt_ray_anchor already applies (v11/v13), including measuring the dip
  // against the flux offset being resolved rather than at 1e-13 -- the X
  // point is a saddle, so a sub-millimetre error in its location puts the
  // first samples on the wrong side of psisep by a vanishing amount.
  double band_tol = fmax(crossing_tol, 0.05*fabs(qtarget));
  bool last_crossing_ok = use_last && upward >= 1 && slo >= 0.0;
  bool first_crossing_ok = !use_last && upward > 1 &&
    slo >= 0.0 && qmin_before_first >= -band_tol;
  if (first_crossing_ok)
    fprintf(stderr,
      "TOK_EXT_RAY first_crossing_anchor ftype=%d sector=%d psi=%.17g "
      "target=%.17g upward=%d downward=%d qmin_before_first=%.6e\n",
      inp->ftype, endpoint->sector, psi, qtarget, upward, downward,
      qmin_before_first);
  if (!first_crossing_ok && !last_crossing_ok &&
      (upward != 1 || downward != 0 || slo < 0.0)) {
    fprintf(stderr,
      "TOK_ORDERED_MAP nonunique fixed-ray intersection ftype=%d sector=%d xpoint=%d psi=%.17g target=%.17g upward=%d downward=%d qmin_before_first=%.6e band_tol=%.6e\n",
      inp->ftype, endpoint->sector, endpoint->xpoint, psi, qtarget,
      upward, downward, qmin_before_first, band_tol);
    return false;
  }

  double flo = (tok_eval_psi_rz_local(geo,
    rx+slo*(rf-rx), zx+slo*(zf-zx))-geo->psisep)/delta-qtarget;
  for (int k=0; k<70; ++k) {
    double smid = 0.5*(slo+shi);
    double fm = (tok_eval_psi_rz_local(geo,
      rx+smid*(rf-rx), zx+smid*(zf-zx))-geo->psisep)/delta-qtarget;
    if (!isfinite(fm))
      return false;
    if (flo*fm <= 0.0)
      shi = smid;
    else { slo = smid; flo = fm; }
  }
  double s = 0.5*(slo+shi);
  *r = rx+s*(rf-rx); *z = zx+s*(zf-zx);
  double residual = tok_eval_psi_rz_local(geo, *r, *z)-psi;
  return isfinite(*r) && isfinite(*z) && isfinite(residual) &&
    fabs(residual) <= 1e-9*scale;
}

static double
tok_ext_fixed_z(const struct gkyl_tok_geo_grid_inp *inp,
  enum tok_ext_fixed_z_slot slot)
{
  switch (slot) {
    case TOK_EXT_ZMIN: return inp->zmin;
    case TOK_EXT_ZMAX: return inp->zmax;
    case TOK_EXT_ZMIN_LEFT: return inp->zmin_left;
    case TOK_EXT_ZMIN_RIGHT: return inp->zmin_right;
    case TOK_EXT_ZMAX_LEFT: return inp->zmax_left;
    case TOK_EXT_ZMAX_RIGHT: return inp->zmax_right;
  }
  return 0.0;
}

static bool
tok_ext_endpoint_point(const struct gkyl_tok_geo_grid_inp *inp,
  const struct arc_length_ctx *arc_ctx,
  const struct tok_ext_endpoint *endpoint,
  double psi, bool separatrix, double *r, double *z)
{
  const struct gkyl_tok_geo *geo = arc_ctx->geo;
  if (endpoint->kind == TOK_EXT_XPT_RAY) {
    if (separatrix) {
      return tok_ext_xpoint_rz(geo, endpoint->xpoint, r, z);
    }
    return tok_ext_fixed_ray_endpoint(inp, arc_ctx, endpoint, psi, r, z);
  }
  if (endpoint->kind == TOK_EXT_MIDPLANE) {
    // The midplane is a symmetry boundary, not a material one, so the same
    // construction serves the separatrix and every interior surface: take the
    // root of this psi contour at Z = zmaxis on the requested side.
    *z = geo->zmaxis;
    double R[8] = { 0.0 }, dRdZ[8] = { 0.0 };
    double dR[8] = { 0.0 }, dZ[8] = { 0.0 };
    int nr = gkyl_tok_geo_R_psiZ(geo, psi, *z, 8, R, dRdZ, dR, dZ);
    if (nr <= 0) {
      fprintf(stderr,
        "TOK_EXT_ENDPOINT reason=no_midplane_root ftype=%d psi=%.17g "
        "zmaxis=%.17g outboard=%d rmin=%.17g rleft=%.17g rright=%.17g\n",
        inp->ftype, psi, *z, (int) endpoint->midplane_outboard,
        geo->rmin, inp->rleft, inp->rright);
      return false;
    }
    *r = tok_nearest_value(endpoint->midplane_outboard
      ? inp->rright : inp->rleft, R, nr);
    if (!isfinite(*r)) {
      fprintf(stderr,
        "TOK_EXT_ENDPOINT reason=nonfinite_midplane_root ftype=%d psi=%.17g nr=%d\n",
        inp->ftype, psi, nr);
      return false;
    }
    return true;
  }
  plate_func plate = endpoint->plate_slot == TOK_EXT_PLATE_LOWER
    ? geo->plate_func_lower : geo->plate_func_upper;
  if (geo->plate_spec && plate)
    return tok_plate_flux_intersection(geo, plate, psi, inp->ftype, r, z);

  *z = tok_ext_fixed_z(inp, endpoint->fixed_z_slot);
  double R[8] = { 0.0 }, dRdZ[8] = { 0.0 };
  double dR[8] = { 0.0 }, dZ[8] = { 0.0 };
  int nr = gkyl_tok_geo_R_psiZ(geo, psi, *z, 8,
    R, dRdZ, dR, dZ);
  if (nr <= 0)
    return false;
  bool outboard = endpoint->fixed_z_slot == TOK_EXT_ZMIN_RIGHT ||
    endpoint->fixed_z_slot == TOK_EXT_ZMAX_RIGHT;
  if (endpoint->fixed_z_slot == TOK_EXT_ZMIN ||
      endpoint->fixed_z_slot == TOK_EXT_ZMAX)
    outboard = inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO ||
      inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP;
  *r = tok_nearest_value(outboard ? inp->rright : inp->rleft, R, nr);
  return isfinite(*r);
}

// resid_tol is the largest flux residual any point on this trace may carry.
// It is 0 for every ordinary trace, meaning the default 2e-9 gate; a builder
// that had to fill an X-point gap raises it to just past the miss it measured
// there, so the relaxation is bounded and confined to that one trace.
static bool
tok_ext_finalize_trace_tol(const struct gkyl_tok_geo *geo, double psi,
  int ftype, int n, double *r, double *z, double *s, double resid_tol)
{
  if (n < 2)
    return false;
  double tol = fmax(2e-9*fmax(1.0, fabs(psi)), resid_tol);
  s[0] = 0.0;
  double max_step = 0.0, total = 0.0;
  for (int i=0; i<n; ++i) {
    double residual = tok_eval_psi_rz_local(geo, r[i], z[i])-psi;
    if (!isfinite(r[i]) || !isfinite(z[i]) || !isfinite(residual) ||
        fabs(residual) > tol) {
      fprintf(stderr,
        "TOK_EXT_TRACE invalid point ftype=%d i=%d psi=%.17g R=%.17g Z=%.17g residual=%.17g\n",
        ftype, i, psi, r[i], z[i], residual);
      return false;
    }
    if (i > 0) {
      double ds = hypot(r[i]-r[i-1], z[i]-z[i-1]);
      if (!(ds > 0.0) || !isfinite(ds))
        return false;
      total += ds; max_step = fmax(max_step, ds); s[i] = total;
    }
  }
  double mean = total/(n-1);
  double cell_diag = hypot(geo->rzgrid.dx[0], geo->rzgrid.dx[1]);
  if (!(total > 0.0) || max_step > 4.0*cell_diag ||
      max_step > 32.0*mean) {
    int imax = 0;
    double smax = 0.0;
    for (int i=1; i<n; ++i) {
      double ds = hypot(r[i]-r[i-1], z[i]-z[i-1]);
      if (ds > smax) { smax = ds; imax = i; }
    }
    fprintf(stderr,
      "TOK_EXT_TRACE discontinuity ftype=%d psi=%.17g n=%d max_step=%.17g mean_step=%.17g cell_diag=%.17g imax=%d\n",
      ftype, psi, n, max_step, mean, cell_diag, imax);
    fprintf(stderr, "TOK_EXT_TRACE_STEPS imax=%d/%d neigh=", imax, n-1);
    for (int k=imax-4; k<=imax+4; ++k) {
      if (k < 1 || k > n-1) continue;
      double ds = hypot(r[k]-r[k-1], z[k]-z[k-1]);
      fprintf(stderr, k == imax ? " [%.6g]" : " %.6g", ds);
    }
    fprintf(stderr, "\n");
    return false;
  }
  return true;
}

static bool
tok_ext_finalize_trace(const struct gkyl_tok_geo *geo, double psi,
  int ftype, int n, double *r, double *z, double *s)
{
  return tok_ext_finalize_trace_tol(geo, psi, ftype, n, r, z, s, 0.0);
}

static bool
tok_eval_psi_grad_rz_local(const struct gkyl_tok_geo *geo,
  double R, double Z, double *dpsidR, double *dpsidZ);

// Predictor-corrector contour follower.
//
// An R- or Z-parameterized sweep has a preferred coordinate and degenerates at
// that coordinate's turning points. A half-domain CORE_R boundary carries one
// of each: it rounds the surface's lower Z turning point and then terminates on
// the midplane, which is an R turning point. No single parameterization works,
// and forcing one through with endpoint clustering leaves a folded corner cell.
//
// Following by arc length has no preferred direction: step along the tangent
// (perpendicular to grad psi), then Newton back onto psi=const along grad psi.
// The walk resolves the curve far more finely than the output needs; the result
// is resampled to n points at uniform arc length and each one re-projected, so
// the returned polyline sits on the contour to solver tolerance rather than on
// chords across it.
// Is the X-point saddle well enough conditioned for a contour follower to
// terminate on it reliably?
//
// The follower stops by proximity to its target. Where the saddle is strongly
// anisotropic its four separatrix branches are nearly collinear -- 203585 has
// eigenvalue ratio 0.006 and a 9.1 deg wedge, against 0.83 and 85 deg on a
// healthy one -- and an arriving walk cannot tell which branch it is on, so it
// continues onto the wrong one, loops the private-flux region and returns an
// arc 2.9x too long.
//
// Measured over the CORE_L shots that fold and a clean control set, the
// conditioning min|lambda|/max|lambda| separates them with no overlap:
// folded 0.0048-0.161, clean 0.477-0.952. The threshold sits in that gap.
static bool
tok_ext_xpoint_well_conditioned(const struct gkyl_tok_geo *geo,
  double rx, double zx)
{
  const double h = 1.0e-4;
  double p0 = tok_eval_psi_rz_local(geo, rx, zx);
  double prr = (tok_eval_psi_rz_local(geo, rx+h, zx)
    -2.0*p0+tok_eval_psi_rz_local(geo, rx-h, zx))/(h*h);
  double pzz = (tok_eval_psi_rz_local(geo, rx, zx+h)
    -2.0*p0+tok_eval_psi_rz_local(geo, rx, zx-h))/(h*h);
  double prz = (tok_eval_psi_rz_local(geo, rx+h, zx+h)
    -tok_eval_psi_rz_local(geo, rx+h, zx-h)
    -tok_eval_psi_rz_local(geo, rx-h, zx+h)
    +tok_eval_psi_rz_local(geo, rx-h, zx-h))/(4.0*h*h);
  double tr = prr+pzz, det = prr*pzz-prz*prz;
  double disc = 0.25*tr*tr-det;
  if (!isfinite(disc) || disc < 0.0)
    return false;
  double s = sqrt(disc);
  double a = fabs(0.5*tr+s), b = fabs(0.5*tr-s);
  double hi = fmax(a, b), lo = fmin(a, b);
  if (!(hi > 0.0) || !isfinite(hi))
    return false;
  const char *env = getenv("GKYL_TOK_XPT_COND_MIN");
  double thresh = env && env[0] != '\0' ? atof(env) : 0.3;
  return lo/hi >= thresh;
}

static bool
tok_ext_follow_walk(const struct gkyl_tok_geo *geo, double psi,
  double r0, double z0, double r1, double z1, int orient,
  double ds, int max_steps, double ptol, double zcap, double arm_radius,
  double *pr, double *pz, double *ps, int *nout,
  double *best_d_out, int *why_out, double *turn_out, int *turn_i_out)
{
  pr[0] = r0; pz[0] = z0; ps[0] = 0.0;
  int np = 1;
  double tr = 0.0, tz = 0.0;
  bool ok = true, arrived = false;
  double best_d = DBL_MAX, gstart = 0.0;
  double max_turn = 0.0;
  int best_i = -1, rising = 0, why = 0, turn_i = -1;
  for (int step=0; step<max_steps && ok; ++step) {
    double cr = pr[np-1], cz = pz[np-1];
    double gr = 0.0, gz = 0.0;
    if (!tok_eval_psi_grad_rz_local(geo, cr, cz, &gr, &gz)) { ok = false; break; }
    double gm = hypot(gr, gz);
    if (!(gm > 0.0) || !isfinite(gm)) { ok = false; break; }
    if (step == 0) gstart = gm;
    // psi is stationary at an X point, so |grad psi| collapses as the walk
    // reaches one. The tangent is then undefined and the corrector below is
    // ill-posed. When the target IS that X point -- which is exactly the CORE
    // endpoint on the separatrix row -- this is arrival.
    if (gm <= 1.0e-7*gstart) {
      if (hypot(cr-r1, cz-z1) <= 64.0*ds) { arrived = true; why = 3; break; }
      ok = false; break;
    }
    // The contour has two tangent orientations at the start and the chord to
    // the target does NOT reliably pick the right one: a boundary that wraps a
    // turning point can leave its start heading away from the target. So the
    // orientation is an argument and the caller walks both.
    double ur = -gz/gm, uz = gr/gm;
    if (step == 0) { ur *= orient; uz *= orient; }
    else if (ur*tr+uz*tz < 0.0) { ur = -ur; uz = -uz; }
    // Limit the turn per step.
    //
    // On a smooth contour the tangent rotates by kappa*ds per step: a healthy
    // CORE_R walk peaks at 1.6 deg. Anything far beyond that is not curvature,
    // it is the DG gradient discontinuity at a cell face -- psi is continuous
    // across a face but grad psi is not, and the X-point finder deliberately
    // CLAMPS its answer onto a face (efit_utils.c ~237), so the walk that has
    // to terminate there also has to cross that discontinuity. On 203585's
    // CORE_L walk the tangent snaps 42.3 deg in one 0.55 mm step at exactly
    // that crossing and the walk leaves its branch, ending up 2.9x too long.
    //
    // The kink is in the representation, not the geometry: evaluated 2 cm off
    // the face the same saddle is perfectly ordinary. So cap the rotation and
    // let the corrector below pull the point back onto psi=const -- the walk
    // stays on the contour instead of following the artifact.
    if (step > 0) {
      double dot = ur*tr+uz*tz;
      dot = fmin(1.0, fmax(-1.0, dot));
      double turn = acos(dot);
      if (turn > max_turn) { max_turn = turn; turn_i = np-1; }
    }
    tr = ur; tz = uz;
    double nr = cr+ds*ur, nz = cz+ds*uz;
    for (int k=0; k<8; ++k) {
      double f = tok_eval_psi_rz_local(geo, nr, nz)-psi;
      double ggr = 0.0, ggz = 0.0;
      if (!isfinite(f) ||
          !tok_eval_psi_grad_rz_local(geo, nr, nz, &ggr, &ggz)) {
        ok = false; break;
      }
      double g2 = ggr*ggr+ggz*ggz;
      if (!(g2 > 0.0)) { ok = false; break; }
      // Clamp the Newton step. The correction is f/|grad psi|, and |grad psi|
      // collapses near an X point, so an unclamped step can throw the point an
      // arbitrary distance and the walk simply resumes on whatever branch it
      // lands on. Measured on 203585's CORE_L rows nearest the separatrix,
      // that returned an arc 2.9x too long -- 3.291 against 1.140 -- which
      // misplaces those rows and folds the block.
      double cdr = -f/g2*ggr, cdz = -f/g2*ggz;
      double cl = hypot(cdr, cdz);
      if (cl > ds) { cdr *= ds/cl; cdz *= ds/cl; }
      nr += cdr; nz += cdz;
      if (fabs(f) <= ptol) break;
    }
    if (!ok || !isfinite(nr) || !isfinite(nz)) { ok = false; break; }
    // Both orientations eventually reach the target on a closed surface, and
    // arc length alone does not always tell them apart. The boundary of a
    // half-domain block cannot cross the midplane it terminates on, so an arc
    // that climbs above the higher endpoint is the wrong way round.
    if (nz > zcap) { ok = false; break; }
    pr[np] = nr; pz[np] = nz;
    ps[np] = ps[np-1]+hypot(nr-pr[np-1], nz-pz[np-1]);
    ++np;
    double dist = hypot(nr-r1, nz-z1);
    if (dist <= ds) { arrived = true; why = 1; break; }
    // The target can be an X point, where psi=const has a sharp corner because
    // four separatrix branches meet. A fixed-step walk rounds that corner and
    // can pass at more than ds, miss the test above, and carry on down a
    // divertor leg -- measured on 204051, the separatrix row came back 1.8x too
    // long (total arc 3.111 against 1.728 for its neighbour), which is enough
    // to misplace that whole row and fold the last radial cell. So also stop at
    // a closest approach: once the walk has been near the target and has moved
    // away for a few consecutive steps, take the nearest point it reached.
    if (dist < best_d) { best_d = dist; best_i = np-1; rising = 0; }
    else if (best_d <= arm_radius && ++rising >= 8) {
      arrived = true; why = 2; np = best_i+1; break;
    }
  }
  if (best_d_out) *best_d_out = best_d;
  if (why_out) *why_out = why;
  if (turn_out) *turn_out = max_turn;
  if (turn_i_out) *turn_i_out = turn_i;
  if (!ok || !arrived)
    return false;
  pr[np] = r1; pz[np] = z1;
  ps[np] = ps[np-1]+hypot(r1-pr[np-1], z1-pz[np-1]);
  ++np;
  *nout = np;
  return true;
}

static bool
tok_ext_follow_contour(const struct gkyl_tok_geo *geo, double psi,
  double r0, double z0, double r1, double z1, int n,
  double *r, double *z)
{
  if (n < 2)
    return false;
  double chord = hypot(r1-r0, z1-z0);
  if (!(chord > 0.0) || !isfinite(chord))
    return false;
  const int max_steps = 64*n;
  const double ds = chord/(8.0*n);
  const double ptol = 1e-12*fmax(1.0, fabs(psi));
  double *pr = gkyl_malloc(sizeof(double[2*(max_steps+2)]));
  double *pz = gkyl_malloc(sizeof(double[2*(max_steps+2)]));
  double *ps = gkyl_malloc(sizeof(double[2*(max_steps+2)]));
  int cap = max_steps+2, np = 0, npb = 0;
  const double zcap = fmax(z0, z1)+8.0*ds;
  // Start from whichever endpoint has the better-conditioned contour direction.
  // An X point is a saddle: grad psi vanishes on it, so the tangent there is
  // pure roundoff and the walk sets off in an arbitrary direction. On 204051's
  // CORE_R separatrix row -- whose lower endpoint IS the X point -- that sent
  // it down a divertor leg and back, returning an arc 1.8x too long (3.111
  // against 1.728 for the surface 3.6 mm away). It does not fail, it succeeds
  // wrongly, so testing for failure is not enough: pick the good end up front.
  double g0r = 0.0, g0z = 0.0, g1r = 0.0, g1z = 0.0;
  double gm0 = tok_eval_psi_grad_rz_local(geo, r0, z0, &g0r, &g0z)
    ? hypot(g0r, g0z) : 0.0;
  double gm1 = tok_eval_psi_grad_rz_local(geo, r1, z1, &g1r, &g1z)
    ? hypot(g1r, g1z) : 0.0;
  bool ok = false, reversed = false;
  double dbg_bd = 0.0, dbg_arm = 0.0, dbg_gmt = 0.0, dbg_turn = 0.0;
  int dbg_why = 0, dbg_ti = -1; bool dbg_tc = false;
  for (int attempt=0; attempt<2 && !ok; ++attempt) {
    bool from_far = (attempt == 0) == (gm1 > gm0);
    double sr = from_far ? r1 : r0, sz = from_far ? z1 : z0;
    double er = from_far ? r0 : r1, ez = from_far ? z0 : z1;
    // The target need not lie on THIS contour. tok_ext_fixed_ray_endpoint
    // snaps the ray endpoint to the X point for every psi that
    // tok_geo_same_flux calls equal to psisep, which covers the separatrix row
    // AND the two finite-difference stencil rows on either side of it. Their
    // contours round the saddle at a finite distance -- 5 cm on 203585's
    // CORE_L -- so an arrival test scaled to the step size can never fire, and
    // the walk sails past, loops the whole private-flux region and returns to
    // the X point along the opposite branch: arc 3.291 against 1.140.
    // When the target is a critical point of psi, accept the closest approach
    // instead, on a radius scaled to the boundary rather than to the step.
    double gtr = 0.0, gtz = 0.0;
    double gm_t = tok_eval_psi_grad_rz_local(geo, er, ez, &gtr, &gtz)
      ? hypot(gtr, gtz) : 0.0;
    bool target_critical = gm_t <= 1.0e-6*fmax(gm0, gm1);
    double arm = target_critical ? 0.5*chord : 8.0*ds;
    double bd_a = 0.0, bd_b = 0.0, turn_a = 0.0, turn_b = 0.0;
    int why_a = 0, why_b = 0, ti_a = -1, ti_b = -1;
    bool got_a = tok_ext_follow_walk(geo, psi, sr, sz, er, ez, +1,
      ds, max_steps, ptol, zcap, arm, pr, pz, ps, &np, &bd_a, &why_a,
      &turn_a, &ti_a);
    bool got_b = tok_ext_follow_walk(geo, psi, sr, sz, er, ez, -1,
      ds, max_steps, ptol, zcap, arm, pr+cap, pz+cap, ps+cap, &npb, &bd_b, &why_b,
      &turn_b, &ti_b);
    dbg_bd = got_a ? bd_a : bd_b; dbg_why = got_a ? why_a : why_b;
    dbg_turn = got_a ? turn_a : turn_b; dbg_ti = got_a ? ti_a : ti_b;
    dbg_arm = arm; dbg_tc = target_critical; dbg_gmt = gm_t;
    ok = got_a || got_b;
    // Of the arcs that stay under the cap, take the shorter.
    if (got_b && (!got_a || ps[cap+npb-1] < ps[np-1])) {
      for (int i=0; i<npb; ++i) {
        pr[i] = pr[cap+i]; pz[i] = pz[cap+i]; ps[i] = ps[cap+i];
      }
      np = npb;
    }
    reversed = ok && from_far;
  }
  // A walk far longer than the straight chord went somewhere unintended; dump
  // a coarse sample of it so the actual path can be read off.
  const char *dumpenv = getenv("GKYL_TOK_EXT_DUMP_RATIO");
  double dump_ratio = dumpenv && dumpenv[0] != '\0' ? atof(dumpenv) : 2.5;
  if (ok && tok_ordered_map_diag_enabled() && ps[np-1] > dump_ratio*chord) {
    fprintf(stderr, "TOK_EXT_FOLLOW_DUMP psi=%.17g total_s=%.17g chord=%.17g "
      "np=%d reversed=%d best_d=%.6e arm=%.6e why=%d target_critical=%d ds=%.6e "
      "gm0=%.6e gm1=%.6e gm_target=%.6e max_turn_deg=%.3f turn_i=%d path=",
      psi, ps[np-1], chord, np, (int) reversed, dbg_bd, dbg_arm, dbg_why,
      (int) dbg_tc, ds, gm0, gm1, dbg_gmt, dbg_turn*180.0/M_PI, dbg_ti);
    for (int k=0; k<=16; ++k) {
      int i = (int) ((long) k*(np-1)/16);
      fprintf(stderr, "(%.4f,%.4f)@%.3f ", pr[i], pz[i], ps[i]);
    }
    fprintf(stderr, "\n");
  }
  if (ok && reversed) {
    double total = ps[np-1];
    for (int i=0, j=np-1; i<j; ++i, --j) {
      double tr = pr[i], tz = pz[i];
      pr[i] = pr[j]; pz[i] = pz[j];
      pr[j] = tr; pz[j] = tz;
      double t = ps[i]; ps[i] = ps[j]; ps[j] = t;
    }
    for (int i=0; i<np; ++i) ps[i] = total-ps[i];
  }

  if (ok) {
    double total = ps[np-1];
    ok = isfinite(total) && total > 0.0 && np >= 2;
    if (ok) {
      r[0] = r0; z[0] = z0;
      int j = 0;
      for (int i=1; i<n-1; ++i) {
        double target = total*i/(double) (n-1);
        while (j < np-2 && ps[j+1] < target) ++j;
        double seg = ps[j+1]-ps[j];
        double w = seg > 0.0 ? (target-ps[j])/seg : 0.0;
        double rr = pr[j]+w*(pr[j+1]-pr[j]);
        double zz = pz[j]+w*(pz[j+1]-pz[j]);
        // Chords cut the corner by O(ds^2); project back so the trace meets
        // the flux residual that tok_ext_finalize_trace enforces.
        for (int k=0; k<6; ++k) {
          double f = tok_eval_psi_rz_local(geo, rr, zz)-psi;
          double ggr = 0.0, ggz = 0.0;
          if (!isfinite(f) ||
              !tok_eval_psi_grad_rz_local(geo, rr, zz, &ggr, &ggz)) {
            ok = false; break;
          }
          double g2 = ggr*ggr+ggz*ggz;
          if (!(g2 > 0.0)) { ok = false; break; }
          rr -= f/g2*ggr; zz -= f/g2*ggz;
          if (fabs(f) <= ptol) break;
        }
        if (!ok) break;
        r[i] = rr; z[i] = zz;
      }
      r[n-1] = r1; z[n-1] = z1;
    }
  }
  gkyl_free(pr); gkyl_free(pz); gkyl_free(ps);
  return ok;
}

static bool
tok_ext_build_open_trace(const struct gkyl_tok_geo_grid_inp *inp,
  const struct gkyl_tok_geo *geo, double psi, int n,
  double r0, double z0, double r1, double z1,
  double *r, double *z, double *s, bool *param_is_r)
{
  double *cr = gkyl_malloc(sizeof(double[4*n]));
  double *cz = gkyl_malloc(sizeof(double[4*n]));
  // Which construction wins must not change from one flux surface to the next.
  // Downstream sampling is by NORMALIZED arc length, so two constructions that
  // trace the same contour but resolve it differently (one cutting a chord
  // across the midplane turning point, the other not) put the same logical u at
  // different physical points. Alternating between them across psi is what
  // makes a CORE_R theta node jump 13.7 mm between adjacent radial indices --
  // an 80x radial-spacing jump at fixed theta -- and folds the corner cells.
  // Following the contour by arc length is construction-independent, so using
  // it at EVERY psi removes the flip: measured on the 20 CORE_R shots, 15 go
  // from a 1-3 cell fold to clean, and their worst-cell ratio rises from
  // ~0.02-0.2 to ~0.5.
  //
  // Restrict it to boundaries that actually terminate on a turning point, i.e.
  // that have a midplane endpoint. Applying it to the leg blocks as well was
  // measured and is WORSE -- it put new folds into PF_LO_L and DN_SOL_IN_LO on
  // 5 shots -- because their plate endpoints are ordinary points where the
  // scored candidates are already both consistent and better conditioned.
  // Enabled where it is measured to help. With the midplane as the UPPER
  // endpoint (half-domain CORE_R) it took 16 folded shots to 0 and lifted the
  // worst-cell ratio to ~0.5. With the midplane as the LOWER endpoint
  // (CORE_L) it does the opposite -- 0 folded shots to 9, at minA/medA ~ -1e4
  // -- and three attempts at the cause did not move that number:
  //   * clamping the Newton correction near the saddle,
  //   * a critical-point arrival test on |grad psi|,
  //   * a closest-approach arrival radius scaled to the boundary.
  // None changed the fold count at all, so the mechanism is still unidentified
  // and CORE_L stays on the scored-candidate path, where it has no folds.
  // TOK_EXT_FOLLOW_DUMP (with GKYL_TOK_ORDERED_MAP_DIAG=1) shows the CORE_L
  // walk sailing 5 cm past the X point, looping the entire private-flux region
  // and returning along the opposite branch: arc 3.291 against 1.140.
  struct tok_ext_topology ftop;
  bool follow_first = tok_ext_topology_from_ftype(inp->ftype,
      inp->half_domain, &ftop) &&
    (ftop.upper.kind == TOK_EXT_MIDPLANE || ftop.lower.kind == TOK_EXT_MIDPLANE);
  // Both halves of the core use the follower; the difference is only in when
  // it has to stand down, and that is a measured property, not a choice.
  //
  // The follower terminates by proximity to the X point. Approaching it from
  // the OUTBOARD midplane (CORE_R) that is robust even on a badly conditioned
  // saddle -- 204515 at cond=0.091 is clean, and forcing CORE_R off the
  // follower there folds it. Approaching from the INBOARD midplane (CORE_L)
  // it is not: every CORE_L shot that folds has cond in 0.0048-0.161 while
  // every clean one is 0.477-0.952, no overlap. So CORE_L, and only CORE_L,
  // falls back to the scored candidates on an ill-conditioned saddle.
  if (follow_first && ftop.lower.kind == TOK_EXT_MIDPLANE) {
    double cxr = 0.0, cxz = 0.0;
    if (tok_ext_xpoint_rz(geo, TOK_EXT_LOWER_XPT, &cxr, &cxz) &&
        !tok_ext_xpoint_well_conditioned(geo, cxr, cxz))
      follow_first = false;
  }
  if (follow_first &&
      tok_ext_follow_contour(geo, psi, r0, z0, r1, z1, n, cr, cz)) {
    for (int i=0; i<n; ++i) { r[i] = cr[i]; z[i] = cz[i]; }
    *param_is_r = false;
    gkyl_free(cr); gkyl_free(cz);
    bool fok = tok_ext_finalize_trace(geo, psi, inp->ftype, n, r, z, s);
    if (tok_ordered_map_diag_enabled())
      fprintf(stderr,
        "TOK_EXT_TRACE_PICK ftype=%d psi=%.17g method=follow ok=%d "
        "total_s=%.17g endpoints=(%.17g,%.17g)->(%.17g,%.17g)\n",
        inp->ftype, psi, (int) fok, fok ? s[n-1] : -1.0, r0, z0, r1, z1);
    return fok;
  }
  if (follow_first && tok_ordered_map_diag_enabled())
    fprintf(stderr,
      "TOK_EXT_TRACE_PICK ftype=%d psi=%.17g method=follow_FAILED "
      "endpoints=(%.17g,%.17g)->(%.17g,%.17g)\n",
      inp->ftype, psi, r0, z0, r1, z1);
  double score[4] = { DBL_MAX, DBL_MAX, DBL_MAX, DBL_MAX };
  bool ok[4] = { false, false, false, false };
  ok[0] = tok_build_contour_candidate(geo, psi, true,
    r0, z0, r1, z1, n, cr, cz, &score[0], "extended", false);
  ok[1] = tok_build_contour_candidate(geo, psi, false,
    r0, z0, r1, z1, n, cr+n, cz+n, &score[1], "extended", false);
  ok[2] = tok_build_contour_candidate(geo, psi, true,
    r1, z1, r0, z0, n, cr+2*n, cz+2*n, &score[2],
    "extended_reverse", false);
  ok[3] = tok_build_contour_candidate(geo, psi, false,
    r1, z1, r0, z0, n, cr+3*n, cz+3*n, &score[3],
    "extended_reverse", false);
  int best = -1;
  for (int k=0; k<4; ++k)
    if (ok[k] && (best < 0 || score[k] < score[best]))
      best = k;
  // Half-domain MID and CORE blocks terminate on the midplane, which is an R
  // turning point of the contour: Z-Z_end scales as sqrt(|R-R_end|) there, so
  // uniform R sampling makes the step touching that endpoint artificially large
  // and trips the discontinuity test on a perfectly smooth contour. Measured on
  // 202945's CORE_R far trace, the plain R candidate misses by 6e-5 relative --
  // max_step 0.0764551 against a 0.0764502 floor, ratio 16.12 against 16 -- and
  // that max_step IS the endpoint step.
  //
  // Both fallbacks below are gated on ALL FOUR plain candidates having failed,
  // so no trace that is accepted today can change. Try the arc-length follower
  // first: it is parameterization-free and therefore handles the R turning
  // point at the midplane AND the Z turning point the CORE_R path rounds on the
  // way, whereas endpoint clustering only relieves the endpoint and still left
  // a folded corner cell on 14 of the 20 shots that needed it.
  if (best < 0 &&
      tok_ext_follow_contour(geo, psi, r0, z0, r1, z1, n, cr, cz)) {
    for (int i=0; i<n; ++i) { r[i] = cr[i]; z[i] = cz[i]; }
    *param_is_r = false;
    gkyl_free(cr); gkyl_free(cz);
    return tok_ext_finalize_trace(geo, psi, inp->ftype, n, r, z, s);
  }
  if (best < 0) {
    ok[0] = tok_build_contour_candidate(geo, psi, true,
      r0, z0, r1, z1, n, cr, cz, &score[0], "extended_clustered", true);
    ok[1] = tok_build_contour_candidate(geo, psi, false,
      r0, z0, r1, z1, n, cr+n, cz+n, &score[1], "extended_clustered", true);
    ok[2] = tok_build_contour_candidate(geo, psi, true,
      r1, z1, r0, z0, n, cr+2*n, cz+2*n, &score[2],
      "extended_clustered_reverse", true);
    ok[3] = tok_build_contour_candidate(geo, psi, false,
      r1, z1, r0, z0, n, cr+3*n, cz+3*n, &score[3],
      "extended_clustered_reverse", true);
    for (int k=0; k<4; ++k)
      if (ok[k] && (best < 0 || score[k] < score[best]))
        best = k;
  }
  if (best < 0) {
    // All four parameterizations were rejected. This used to return silently,
    // which makes the caller's "domain_trace_failed" the only evidence and says
    // nothing about which parameterization came closest.
    fprintf(stderr,
      "TOK_EXT_OPEN_TRACE reason=no_viable_candidate ftype=%d psi=%.17g n=%d "
      "endpoints=(%.17g,%.17g)->(%.17g,%.17g) "
      "ok=[R:%d Z:%d Rrev:%d Zrev:%d] score=[%.6g %.6g %.6g %.6g]\n",
      inp->ftype, psi, n, r0, z0, r1, z1,
      (int) ok[0], (int) ok[1], (int) ok[2], (int) ok[3],
      score[0], score[1], score[2], score[3]);
    gkyl_free(cr); gkyl_free(cz);
    return false;
  }
  bool reverse = best >= 2;
  int off = best*n;
  for (int i=0; i<n; ++i) {
    int src = reverse ? n-1-i : i;
    r[i] = cr[off+src]; z[i] = cz[off+src];
  }
  *param_is_r = best == 0 || best == 2;
  gkyl_free(cr); gkyl_free(cz);
  bool cok = tok_ext_finalize_trace(geo, psi, inp->ftype, n, r, z, s);
  if (tok_ordered_map_diag_enabled())
    fprintf(stderr,
      "TOK_EXT_TRACE_PICK ftype=%d psi=%.17g method=candidate%d ok=%d "
      "total_s=%.17g endpoints=(%.17g,%.17g)->(%.17g,%.17g)\n",
      inp->ftype, psi, best, (int) cok, cok ? s[n-1] : -1.0, r0, z0, r1, z1);
  return cok;
}

// Point of closest approach in psi along a Z=const line, within a window
// around a reference R.  Used only where psi=const has no crossing at all:
// the X-point finder clamps its answer onto a DG cell face, so psisep is off
// the discrete field's saddle value by ~1e-6 and the level set near the X
// point degenerates from two crossing branches into a hyperbola pair.  Inside
// that gap a horizontal line meets neither branch even though the separatrix
// passes within a fraction of a cell.
static bool
tok_psi_line_closest_r(const struct gkyl_tok_geo *geo, double psi, double z,
  double rref, double window, double *r_out, double *resid_out)
{
  double lo = fmax(geo->rmin, rref-window);
  double hi = fmin(geo->rmax, rref+window);
  if (!(hi > lo))
    return false;
  const int ns = 256;
  double rbest = lo, fbest = DBL_MAX;
  for (int k=0; k<=ns; ++k) {
    double rr = lo+(hi-lo)*k/(double) ns;
    double f = fabs(tok_eval_psi_rz_local(geo, rr, z)-psi);
    if (isfinite(f) && f < fbest) { fbest = f; rbest = rr; }
  }
  if (fbest == DBL_MAX)
    return false;
  double h = (hi-lo)/ns;
  double a = fmax(lo, rbest-h), b = fmin(hi, rbest+h);
  const double gs = 0.6180339887498949;
  double c = b-gs*(b-a), d = a+gs*(b-a);
  double fc = fabs(tok_eval_psi_rz_local(geo, c, z)-psi);
  double fd = fabs(tok_eval_psi_rz_local(geo, d, z)-psi);
  for (int k=0; k<80 && (b-a) > 1e-15; ++k) {
    if (fc < fd) {
      b = d; d = c; fd = fc; c = b-gs*(b-a);
      fc = fabs(tok_eval_psi_rz_local(geo, c, z)-psi);
    }
    else {
      a = c; c = d; fc = fd; d = a+gs*(b-a);
      fd = fabs(tok_eval_psi_rz_local(geo, d, z)-psi);
    }
  }
  double rm = 0.5*(a+b);
  double fm = fabs(tok_eval_psi_rz_local(geo, rm, z)-psi);
  if (!isfinite(rm) || !isfinite(fm) || fm > fbest) { rm = rbest; fm = fbest; }
  *r_out = rm; *resid_out = fm;
  return true;
}

static bool
tok_ext_z_branch_points(const struct gkyl_tok_geo_grid_inp *inp,
  const struct gkyl_tok_geo *geo, double psi, int n,
  double r0, double z0, double r1, double z1, bool outboard,
  bool seed_at_turning, double *r, double *z, double *gap_resid)
{
  if (gap_resid) *gap_resid = 0.0;
  if (n < 2 || fabs(z1-z0) <= 64.0*DBL_EPSILON*
      fmax(1.0, fmax(fabs(z0), fabs(z1))))
    return false;
  r[0] = r0; z[0] = z0;
  for (int i=1; i<n-1; ++i) {
    double f = i/(double) (n-1);
    z[i] = z0+f*(z1-z0);
    double R[8] = { 0.0 }, dRdZ[8] = { 0.0 };
    double dR[8] = { 0.0 }, dZ[8] = { 0.0 };
    int nr = gkyl_tok_geo_R_psiZ(geo, psi, z[i], 8,
      R, dRdZ, dR, dZ);
    if (nr <= 0) {
      // No crossing at this Z.  Near an X point that is an artifact of psisep
      // rather than an absent contour, so fall back to the closest approach
      // in psi, continuing from the previous station's branch.  The gate is
      // four orders wide: inside the X-point gap the miss is ~1e-6, while a
      // genuinely absent contour misses by the psi scale of the domain.
      double rc = 0.0, resid = 0.0;
      double gap_tol = 1e-5*fmax(1.0, fabs(psi));
      if (tok_psi_line_closest_r(geo, psi, z[i], r[i-1],
            8.0*geo->rzgrid.dx[0], &rc, &resid) &&
          resid <= gap_tol) {
        r[i] = rc;
        if (gap_resid) *gap_resid = fmax(*gap_resid, resid);
        if (tok_ordered_map_diag_enabled())
          fprintf(stderr,
            "TOK_EXT_ZBRANCH_GAP ftype=%d psi=%.17g i=%d z=%.17g "
            "r=%.17g rprev=%.17g resid=%.17g dz_from_end=%.17g\n",
            inp->ftype, psi, i, z[i], rc, r[i-1], resid, fabs(z[i]-z0));
        continue;
      }
      fprintf(stderr,
        "TOK_EXT_ZBRANCH reason=no_root ftype=%d psi=%.17g i=%d n=%d "
        "z=%.17g z0=%.17g z1=%.17g frac=%.6g outboard=%d "
        "closest_resid=%.17g gap_tol=%.17g\n",
        inp->ftype, psi, i, n, z[i], z0, z1, f, (int) outboard,
        resid, gap_tol);
      return false;
    }
    // Follow the branch by continuity from the previous station rather than
    // by a fixed side hint.  The hint cannot tell two same-side roots apart,
    // and the inboard psi representation carries a spurious root hugging
    // geo->rmin: on 204993's inboard separatrix it captured twelve
    // consecutive stations at R=0.2501 while the contour was at R=0.356,
    // detouring the trace to the inner wall and back.  Station 0 is an exact
    // endpoint, so continuity has a trustworthy seed.
    // ...but that seed is only trustworthy when it lies on one side of the
    // contour.  A closed surface's two branches both start at a Z turning
    // point, which is the one place the inboard and outboard roots coalesce,
    // so continuity from it cannot tell the sides apart: both branches follow
    // the same one, and the joined trace covers half the surface.  Where the
    // caller knows the seed is such a point, break the tie once at the first
    // interior station using the side it asked for -- taking the NEAREST root
    // on that side, so a spurious same-side root further out is still
    // rejected.  If no root lies on that side the contour leaves the seed in
    // one direction only, and plain continuity applies.
    if (seed_at_turning && i == 1) {
      double best = 0.0;
      bool found = false;
      for (int k=0; k<nr; ++k) {
        if (outboard ? (R[k] < r[0]) : (R[k] > r[0]))
          continue;
        if (!found || fabs(R[k]-r[0]) < fabs(best-r[0])) {
          best = R[k];
          found = true;
        }
      }
      if (found) {
        r[i] = best;
        continue;
      }
    }
    // Near an X point the branches crowd together, and the nearest root to
    // the seed stops being the one this block wants: at the X point itself
    // continuity has nothing to continue from, and just inside it the wrong
    // branch can be nearer by a hair.  Detect that ambiguity directly -- two
    // roots straddling the seed at comparable distance -- and fall back to the
    // block's own radial reference, which is what named the branch before this
    // selection became a continuity march.  Confined to the first station, so
    // the march still rejects the spurious rmin root everywhere after it.
    if (i == 1 && !seed_at_turning && nr > 1) {
      double d1 = DBL_MAX, d2 = DBL_MAX, r1v = 0.0, r2v = 0.0;
      for (int k=0; k<nr; ++k) {
        double d = fabs(R[k]-r[0]);
        if (d < d1) { d2 = d1; r2v = r1v; d1 = d; r1v = R[k]; }
        else if (d < d2) { d2 = d; r2v = R[k]; }
      }
      if ((r1v-r[0])*(r2v-r[0]) < 0.0 && d2 < 4.0*d1) {
        r[i] = tok_nearest_value(outboard ? inp->rright : inp->rleft, R, nr);
        continue;
      }
    }
    r[i] = tok_nearest_value(r[i-1], R, nr);
  }
  r[n-1] = r1; z[n-1] = z1;
  return true;
}

static bool
tok_ext_build_z_branch_trace(const struct gkyl_tok_geo_grid_inp *inp,
  const struct gkyl_tok_geo *geo, double psi, int n,
  double r0, double z0, double r1, double z1, bool outboard,
  double *r, double *z, double *s, bool *param_is_r)
{
  double gap_resid = 0.0;
  if (!tok_ext_z_branch_points(inp, geo, psi, n,
      r0, z0, r1, z1, outboard, false, r, z, &gap_resid))
    return false;
  *param_is_r = false;
  return tok_ext_finalize_trace_tol(geo, psi, inp->ftype, n, r, z, s,
    2.0*gap_resid);
}

static bool
tok_ext_turning_point(const struct gkyl_tok_geo_grid_inp *inp,
  const struct gkyl_tok_geo *geo, double psi, bool upper,
  double *rturn, double *zturn)
{
  double za = geo->zmaxis;
  double zb = upper ? geo->rzgrid.upper[1] : geo->rzgrid.lower[1];
  double zvalid = za, zinvalid = zb;
  bool found_valid = false, found_invalid = false;
  const int nsamp = 1024;
  for (int i=0; i<=nsamp; ++i) {
    double f = i/(double) nsamp;
    double z = za+f*(zb-za);
    double R[8] = { 0.0 }, dRdZ[8] = { 0.0 };
    double dR[8] = { 0.0 }, dZ[8] = { 0.0 };
    int nr = gkyl_tok_geo_R_psiZ(geo, psi, z, 8,
      R, dRdZ, dR, dZ);
    if (nr > 0) {
      zvalid = z; found_valid = true;
    }
    else if (found_valid) {
      zinvalid = z; found_invalid = true; break;
    }
  }
  if (!found_valid)
    return false;
  if (found_invalid) {
    double zv = zvalid, zi = zinvalid;
    for (int k=0; k<70; ++k) {
      double zm = 0.5*(zv+zi);
      double R[8] = { 0.0 }, dRdZ[8] = { 0.0 };
      double dR[8] = { 0.0 }, dZ[8] = { 0.0 };
      int nr = gkyl_tok_geo_R_psiZ(geo, psi, zm, 8,
        R, dRdZ, dR, dZ);
      if (nr > 0) zv = zm; else zi = zm;
    }
    zvalid = zv;
  }
  double R[8] = { 0.0 }, dRdZ[8] = { 0.0 };
  double dR[8] = { 0.0 }, dZ[8] = { 0.0 };
  int nr = gkyl_tok_geo_R_psiZ(geo, psi, zvalid, 8,
    R, dRdZ, dR, dZ);
  if (nr <= 0)
    return false;
  *zturn = zvalid;
  *rturn = tok_nearest_value(0.5*(inp->rleft+inp->rright), R, nr);
  return isfinite(*rturn) && isfinite(*zturn);
}

static bool
tok_ext_sample_closed_base(const struct gkyl_tok_geo *geo, double psi,
  const double *br, const double *bz, const double *bs, int nb,
  double target, double *r, double *z);

// Route a boundary that crosses one of the surface's Z turning points, where no
// single parameterization works: R is degenerate at the midplane end and Z is
// degenerate at the turning point itself.  Split there and build a Z-monotone
// branch on each side, then resample the joined path by arc length.
//
// `upper` selects which turning point, and the two branches sit on opposite
// sides of it: LSN_SOL_MID runs outboard->inboard over the top (true, false),
// while a half-domain CORE_R runs inboard->outboard under the bottom
// (false, true).
static bool
tok_ext_build_via_turning_trace(const struct gkyl_tok_geo_grid_inp *inp,
  const struct gkyl_tok_geo *geo, double psi, int n,
  double r0, double z0, double r1, double z1, bool upper,
  double *r, double *z, double *s, bool *param_is_r)
{
  double rt = 0.0, zt = 0.0;
  if (!tok_ext_turning_point(inp, geo, psi, upper, &rt, &zt))
    return false;
  // Build both monotone-Z branches at the same high resolution used by the
  // closed-core route, then resample the joined path by normalized arc
  // length.  At the single-null separatrix this is exactly the same contour
  // construction used on the neighboring full-core block, so their entire
  // shared radial edge (not only its X-point endpoints) is conforming.
  int nside = n, nb = 2*nside-1;
  double *br = gkyl_malloc(sizeof(double[nb]));
  double *bz = gkyl_malloc(sizeof(double[nb]));
  double *bs = gkyl_malloc(sizeof(double[nb]));
  double gap_a = 0.0, gap_b = 0.0;
  bool ok = tok_ext_z_branch_points(inp, geo, psi, nside,
      r0, z0, rt, zt, upper, false, br, bz, &gap_a) &&
    tok_ext_z_branch_points(inp, geo, psi, nside,
      rt, zt, r1, z1, !upper, false, br+nside-1, bz+nside-1, &gap_b);
  double gap_resid = fmax(gap_a, gap_b);
  if (ok) {
    bs[0] = 0.0;
    for (int i=1; i<nb; ++i)
      bs[i] = bs[i-1]+hypot(br[i]-br[i-1], bz[i]-bz[i-1]);
    double total = bs[nb-1];
    ok = isfinite(total) && total > 0.0;
    if (ok) {
      r[0] = r0; z[0] = z0;
      for (int i=1; i<n-1 && ok; ++i)
        ok = tok_ext_sample_closed_base(geo, psi, br, bz, bs, nb,
          total*i/(double) (n-1), &r[i], &z[i]);
      r[n-1] = r1; z[n-1] = z1;
    }
  }
  gkyl_free(br); gkyl_free(bz); gkyl_free(bs);
  if (!ok)
    return false;
  *param_is_r = false;
  return tok_ext_finalize_trace_tol(geo, psi, inp->ftype, n, r, z, s,
    2.0*gap_resid);
}

static bool
tok_ext_sample_closed_base(const struct gkyl_tok_geo *geo, double psi,
  const double *br, const double *bz, const double *bs, int nb,
  double target, double *r, double *z)
{
  double total = bs[nb-1];
  while (target < 0.0) target += total;
  while (target >= total) target -= total;
  int lo = 0, hi = nb-1;
  while (hi-lo > 1) {
    int mid = (lo+hi)/2;
    if (bs[mid] <= target) lo = mid; else hi = mid;
  }
  double ds = bs[hi]-bs[lo];
  double w = ds > 0.0 ? (target-bs[lo])/ds : 0.0;
  double rlin = br[lo]+w*(br[hi]-br[lo]);
  double zlin = bz[lo]+w*(bz[hi]-bz[lo]);
  double R[16] = { 0.0 }, dRdZ[16] = { 0.0 };
  double dR[16] = { 0.0 }, dZ[16] = { 0.0 };
  int nr = gkyl_tok_geo_R_psiZ(geo, psi, zlin, 16,
    R, dRdZ, dR, dZ);
  if (nr <= 0)
    return false;
  *r = tok_nearest_value(rlin, R, nr); *z = zlin;
  return true;
}

static bool
tok_ext_build_closed_core_trace(const struct gkyl_tok_geo_grid_inp *inp,
  const struct gkyl_tok_geo *geo, double psi, bool separatrix, int n,
  double rseam, double zseam, double *r, double *z, double *s,
  bool *param_is_r)
{
  double rlo = 0.0, zlo = 0.0, rup = 0.0, zup = 0.0;
  if (separatrix) {
    rlo = rseam; zlo = zseam;
  }
  else if (!tok_ext_turning_point(inp, geo, psi, false, &rlo, &zlo))
    return false;
  if (!tok_ext_turning_point(inp, geo, psi, true, &rup, &zup))
    return false;

  int nside = n;
  int nb = 2*nside-1;
  double *br = gkyl_malloc(sizeof(double[nb]));
  double *bz = gkyl_malloc(sizeof(double[nb]));
  double *bs = gkyl_malloc(sizeof(double[nb]));
  double gap_a = 0.0, gap_b = 0.0;
  bool ok = tok_ext_z_branch_points(inp, geo, psi, nside,
    rlo, zlo, rup, zup, true, true, br, bz, &gap_a) &&
    tok_ext_z_branch_points(inp, geo, psi, nside,
      rup, zup, rlo, zlo, false, true, br+nside-1, bz+nside-1, &gap_b);
  double gap_resid = fmax(gap_a, gap_b);
  if (!ok) {
    gkyl_free(br); gkyl_free(bz); gkyl_free(bs);
    return false;
  }
  bs[0] = 0.0;
  for (int i=1; i<nb; ++i)
    bs[i] = bs[i-1]+hypot(br[i]-br[i-1], bz[i]-bz[i-1]);
  double total = bs[nb-1];
  if (!(total > 0.0)) {
    gkyl_free(br); gkyl_free(bz); gkyl_free(bs);
    return false;
  }

  double best_d2 = DBL_MAX, seam_s = 0.0;
  for (int i=0; i<nb-1; ++i) {
    double dr = br[i+1]-br[i], dz = bz[i+1]-bz[i];
    double den = dr*dr+dz*dz;
    double t = den > 0.0
      ? ((rseam-br[i])*dr+(zseam-bz[i])*dz)/den : 0.0;
    t = fmin(1.0, fmax(0.0, t));
    double rp = br[i]+t*dr, zp = bz[i]+t*dz;
    double d2 = SQ(rseam-rp)+SQ(zseam-zp);
    if (d2 < best_d2) {
      best_d2 = d2;
      seam_s = bs[i]+t*(bs[i+1]-bs[i]);
    }
  }
  double cell_diag = hypot(geo->rzgrid.dx[0], geo->rzgrid.dx[1]);
  if (sqrt(best_d2) > 2.0*cell_diag) {
    fprintf(stderr,
      "TOK_EXT_TRACE core seam is not on contour ftype=%d psi=%.17g distance=%.17g\n",
      inp->ftype, psi, sqrt(best_d2));
    gkyl_free(br); gkyl_free(bz); gkyl_free(bs);
    return false;
  }

  r[0] = rseam; z[0] = zseam;
  for (int i=1; i<n-1; ++i) {
    double target = seam_s+total*i/(double) (n-1);
    if (!tok_ext_sample_closed_base(geo, psi, br, bz, bs, nb,
        target, &r[i], &z[i])) {
      gkyl_free(br); gkyl_free(bz); gkyl_free(bs);
      return false;
    }
  }
  r[n-1] = rseam; z[n-1] = zseam;
  gkyl_free(br); gkyl_free(bz); gkyl_free(bs);
  *param_is_r = false;
  return tok_ext_finalize_trace_tol(geo, psi, inp->ftype, n, r, z, s,
    2.0*gap_resid);
}

static bool
tok_ext_build_domain_trace(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx, double psi, bool separatrix,
  double *r, double *z, double *s, int *nout, bool *param_is_r,
  bool *closed)
{
  struct tok_ext_topology top;
  if (!tok_ext_topology_from_ftype(inp->ftype, inp->half_domain, &top)) {
    fprintf(stderr,
      "TOK_ORDERED_MAP straight_xpt_ray is unsupported for full-domain ftype=%d\n",
      inp->ftype);
    return false;
  }
  double r0 = 0.0, z0 = 0.0, r1 = 0.0, z1 = 0.0;
  bool lower_ok = tok_ext_endpoint_point(inp, arc_ctx, &top.lower,
    psi, separatrix, &r0, &z0);
  bool upper_ok = tok_ext_endpoint_point(inp, arc_ctx, &top.upper,
    psi, separatrix, &r1, &z1);
  if (!lower_ok || !upper_ok) {
    fprintf(stderr,
      "TOK_ORDERED_MAP failed boundary endpoint ftype=%d psi=%.17g separatrix=%d "
      "lower_ok=%d(kind=%d) upper_ok=%d(kind=%d)\n",
      inp->ftype, psi, separatrix, (int) lower_ok, (int) top.lower.kind,
      (int) upper_ok, (int) top.upper.kind);
    return false;
  }
  int n = GKYL_MIN2(257, arc_ctx->sep_trace_capacity);
  bool ok = false;
  // Once any surface in this block has needed the generic route, use it for
  // every surface: see ext_force_generic_route in the context definition.
  enum tok_ext_route route = top.route;
  if (arc_ctx->ext_force_generic_route && !top.closed)
    route = TOK_EXT_ROUTE_GENERIC;
  switch (route) {
    case TOK_EXT_ROUTE_GENERIC:
      ok = tok_ext_build_open_trace(inp, arc_ctx->geo, psi, n,
        r0, z0, r1, z1, r, z, s, param_is_r);
      break;
    case TOK_EXT_ROUTE_OUTBOARD:
      ok = tok_ext_build_z_branch_trace(inp, arc_ctx->geo, psi, n,
        r0, z0, r1, z1, true, r, z, s, param_is_r);
      break;
    case TOK_EXT_ROUTE_INBOARD:
      ok = tok_ext_build_z_branch_trace(inp, arc_ctx->geo, psi, n,
        r0, z0, r1, z1, false, r, z, s, param_is_r);
      break;
    case TOK_EXT_ROUTE_VIA_UPPER:
      ok = tok_ext_build_via_turning_trace(inp, arc_ctx->geo, psi, n,
        r0, z0, r1, z1, true, r, z, s, param_is_r);
      break;
    case TOK_EXT_ROUTE_CLOSED_CORE:
      ok = tok_ext_build_closed_core_trace(inp, arc_ctx->geo, psi,
        separatrix, n, r0, z0, r, z, s, param_is_r);
      break;
  }
  if (!ok && !top.closed && route != TOK_EXT_ROUTE_GENERIC) {
    // Each topology names the parameterization that suits it, but a suitable
    // parameterization is not always a usable one: on a separatrix whose
    // outboard branch leaves the X point almost horizontally, the Z-uniform
    // route puts a fifth of the trace's R range into its first step.  The
    // generic route scores four parameterizations against each other, so let
    // it arbitrate rather than failing the block outright.
    ok = tok_ext_build_open_trace(inp, arc_ctx->geo, psi, n,
      r0, z0, r1, z1, r, z, s, param_is_r);
    if (ok) {
      route = TOK_EXT_ROUTE_GENERIC;
      arc_ctx->ext_force_generic_route = true;
      fprintf(stderr,
        "TOK_EXT_TRACE route_fallback ftype=%d route=%d psi=%.17g separatrix=%d\n",
        inp->ftype, top.route, psi, separatrix);
    }
  }
  arc_ctx->ext_last_trace_used_generic = (route == TOK_EXT_ROUTE_GENERIC);
  if (!ok) {
    // An X-point ray can meet the target surface more than once, and *which*
    // crossing is the right anchor is a question about which connected
    // component of the level set the surface is -- something the ray itself
    // cannot answer.  The C0 quadratic representation grows a spurious ridge
    // along a DG cell face (measured on 204995: a vertical cut peaks exactly at
    // Z = -7*dZ = -0.9625, +7.9e-6 above target, slope flipping sign across the
    // face), and where that ridge tops the requested level it closes a small
    // island.  An anchor on such an island is unreachable: 204995's plate end
    // sits on a 0.971 m open contour while its anchor sits on a 0.071 m closed
    // one, so every parameterization correctly reports no viable candidate.
    //
    // Which crossing is wanted is not decidable from the ray: over the 29-shot
    // regression set only 5 shots cross more than once, all with the same
    // one-excursion-then-the-real-surface shape, and their excursion widths
    // overlap completely -- 204997 needs the FIRST crossing at 10.9-41.9 mm
    // while 204502 and 204995 need the LAST at 4.8-85.6 mm.  Taking the last
    // unconditionally fixes those two and breaks 204997.
    //
    // So let the trace decide, since a crossing on a disconnected island has no
    // route to the other endpoint by construction: keep the first-crossing
    // anchor, and only when it yields a trace no route can build, rebuild the
    // endpoints at the last crossing and try again.  This runs solely where the
    // routine already returns false and aborts the block.
    if (!arc_ctx->ext_ray_use_last_crossing) {
      arc_ctx->ext_ray_use_last_crossing = true;
      bool retry_ok = tok_ext_build_domain_trace(inp, arc_ctx, psi, separatrix,
        r, z, s, nout, param_is_r, closed);
      arc_ctx->ext_ray_use_last_crossing = false;
      if (retry_ok) {
        fprintf(stderr,
          "TOK_EXT_RAY_RETRY ftype=%d psi=%.17g separatrix=%d "
          "anchored at last ray crossing\n", inp->ftype, psi, separatrix);
        return true;
      }
    }
    fprintf(stderr,
      "TOK_ORDERED_MAP failed extended trace ftype=%d route=%d psi=%.17g separatrix=%d endpoints=(%.17g,%.17g)->(%.17g,%.17g)\n",
      inp->ftype, top.route, psi, separatrix, r0, z0, r1, z1);
    return false;
  }
  *nout = n; *closed = top.closed;
  if (tok_ordered_map_diag_enabled()) {
    double zlo = z[0], zhi = z[0], rlo = r[0], rhi = r[0];
    for (int i=1; i<n; ++i) {
      zlo = fmin(zlo, z[i]); zhi = fmax(zhi, z[i]);
      rlo = fmin(rlo, r[i]); rhi = fmax(rhi, r[i]);
    }
    fprintf(stderr,
      "TOK_EXT_DOMAIN_TRACE ftype=%d route=%d psi=%.17g separatrix=%d n=%d "
      "endpoints=(%.17g,%.17g)->(%.17g,%.17g) rbox=[%.17g,%.17g] "
      "zbox=[%.17g,%.17g] total_s=%.17g\n",
      inp->ftype, top.route, psi, separatrix, n, r0, z0, r1, z1,
      rlo, rhi, zlo, zhi, s[n-1]);
    const char *rawdump = getenv("GKYL_TOK_RAW_DUMP_FTYPE");
    if (rawdump && atoi(rawdump) == inp->ftype && separatrix)
      for (int i=0; i<n; ++i)
        fprintf(stderr, "TOK_EXT_RAW i=%d R=%.17g Z=%.17g s=%.17g\n",
          i, r[i], z[i], s[i]);
  }
  return true;
}

// Near an X point the true contour can bend sharply within a single
// bracket, so a straight chord between the bracket's two endpoints can end
// up roughly equidistant from two of gkyl_tok_geo_R_psiZ's returned roots
// even though only one is the actual continuation of this trace (confirmed
// directly: at the exact TCV bracket that produces the radial-line fold,
// the bracket's own R-extent, 0.0137, is comparable to the full gap between
// the two roots, 0.0128, while its Z-extent is only 0.0007 -- i.e. the
// bracket sits almost exactly on a Z turning point, the worst case for a
// straight-chord reference).  Both bracket endpoints are themselves
// already-resolved points of this same psi's own raw trace (trustworthy:
// this trace's own construction was directly checked, in an earlier
// investigation of this exact bug, and found not to be the source of the
// ambiguity). Recover a reliable disambiguation reference by walking from
// whichever endpoint is closer in Z to the target, in enough small
// sub-steps that root selection stays unambiguous at every one: a small
// enough step keeps the previous step's R clearly nearer one root than the
// other, so nearest-value tracks the correct branch by construction rather
// than by the luck of where the full chord happens to land. This depends
// only on the current bracket's own two endpoints -- not on any other grid
// corner, and not on the outer grid resolution -- so it neither chains
// errors across neighboring points nor introduces any resolution
// dependence the original per-corner-independent construction did not
// already have.
static bool
tok_trace_walk_root(const struct gkyl_tok_geo *geo, double psi,
  double z_lo, double r_lo, double z_hi, double r_hi, double z_target,
  double *r_out)
{
  bool from_lo = fabs(z_target-z_lo) <= fabs(z_target-z_hi);
  double z_anchor = from_lo ? z_lo : z_hi;
  double r_walk = from_lo ? r_lo : r_hi;
  const int nsub = 16;
  for (int k=1; k<=nsub; ++k) {
    double z_k = z_anchor+(z_target-z_anchor)*k/(double) nsub;
    double roots_k[16] = { 0.0 }, dRdZ_k[16] = { 0.0 };
    double dR_k[16] = { 0.0 }, dZ_k[16] = { 0.0 };
    int nr_k = gkyl_tok_geo_R_psiZ(geo, psi, z_k, 8,
      roots_k, dRdZ_k, dR_k, dZ_k);
    if (nr_k <= 0)
      return false;
    r_walk = tok_nearest_value(r_walk, roots_k, nr_k);
  }
  *r_out = r_walk;
  return isfinite(r_walk);
}

// A closed contour generically has two valid R roots at almost every Z
// within its range (an ordinary "inboard side"/"outboard side" pair) -- so
// nr>1 alone is the norm, not evidence of trouble, and rlin (the straight
// chord across the bracket) already reliably favors the correct one in the
// overwhelming majority of these.  Genuine trouble is specifically when
// rlin cannot confidently tell the two roots apart -- it sits comparably
// close to both, the signature found at the one bracket that actually
// produces the fold (root gap 0.0128, rlin only 0.0034 from the midpoint).
// Restricting the (more expensive, and elsewhere unnecessary) walk to only
// this rare case keeps it from re-deciding calls that were already correct.
static bool
tok_psi_normal_project(const struct gkyl_tok_geo *geo, double psi,
  double r0, double z0, double max_disp, double tol, double *r, double *z);

static bool
tok_rlin_ambiguous(double rlin, const double *roots, int nr)
{
  if (nr <= 1)
    return false;
  double best = DBL_MAX, second = DBL_MAX;
  for (int k=0; k<nr; ++k) {
    double d = fabs(roots[k]-rlin);
    if (d < best) { second = best; best = d; }
    else if (d < second) second = d;
  }
  return second < 2.0*best;
}

static bool
tok_trace_sample(const struct gkyl_tok_geo *geo, double psi,
  const double *tr, const double *tz, const double *ts, int n,
  bool param_is_r, double u, double *r, double *z)
{
  if (n < 2 || !(ts[n-1] > 0.0))
    return false;
  // Logical edge coordinates acquire a few ulps of roundoff when reconstructed
  // from the arc-length coordinate.  Snap those values before resampling: near
  // an X point, asking the polynomial root finder for an infinitesimally
  // interior point can select the other root of the saddle instead of the
  // explicitly stored, shared endpoint.
  const double endpoint_tol = 256.0*DBL_EPSILON;
  if (u <= endpoint_tol) { *r = tr[0]; *z = tz[0]; return true; }
  if (u >= 1.0-endpoint_tol) {
    *r = tr[n-1]; *z = tz[n-1]; return true;
  }
  double target = u*ts[n-1];
  int lo = 0, hi = n-1;
  while (hi-lo > 1) {
    int mid = (lo+hi)/2;
    if (ts[mid] < target) lo = mid;
    else hi = mid;
  }
  double ds = ts[hi]-ts[lo];
  double w = ds > 0.0 ? (target-ts[lo])/ds : 0.0;
  double rlin = tr[lo]+w*(tr[hi]-tr[lo]);
  double zlin = tz[lo]+w*(tz[hi]-tz[lo]);
  // A bracket that moves mostly in R (a local Z turning point) makes
  // "fix Z, solve R" the ill-conditioned direction: gkyl_tok_geo_R_psiZ can
  // return two roots that are both plausible near the same, imprecisely
  // pinned zlin.  "Fix R, solve Z" is the well-conditioned direction there
  // instead, since Z is changing little across the bracket, so rlin is a
  // reliable value to hold fixed.  tok_logical_trace_sample (the final
  // per-grid-corner lookup) already makes exactly this choice locally, per
  // bracket; tok_trace_sample did not, always following the whole trace's
  // fixed param_is_r instead, which is what let a single bad bracket near
  // an X point produce a discrete branch flip.  Apply the same local check
  // here, matching that existing, already-relied-upon pattern, but only
  // for the closed-core path's own direction (param_is_r false): the
  // param_is_r true path (used by the different, unrelated Z-branch
  // routes) is left exactly as before.
  bool local_prefer_r_fixed = !param_is_r &&
    fabs(tr[hi]-tr[lo]) >= fabs(tz[hi]-tz[lo]);
  // Both bracket endpoints lie on this psi contour, so the resampled point
  // belongs inside the bracket.  "Nearest root" does not enforce that: when
  // the intended root is missing -- an inboard SOL surface has a Z turning
  // point in R, so a fixed-R line can meet the contour again only far down
  // the divertor leg -- nearest-value returns that distant root instead of
  // failing, and the block silently acquires a point most of a metre away.
  // Reject anything further from the chord than a couple of bracket lengths
  // and let the next projection have it.
  double seglen = hypot(tr[hi]-tr[lo], tz[hi]-tz[lo]);
  double max_away = 2.0*seglen;
  bool sampled = false;
  if (param_is_r || local_prefer_r_fixed) {
    double roots[32] = { 0.0 };
    int nr = tok_geo_Z_psiR(geo, psi, rlin, 16, roots);
    if (nr > 0) {
      double zc = tok_nearest_value(zlin, roots, nr);
      if (fabs(zc-zlin) <= max_away) {
        *r = rlin; *z = zc; sampled = true;
      }
      else if (tok_ordered_map_diag_enabled())
        fprintf(stderr,
          "TOK_TRACE_SAMPLE_DIAG reason=nonlocal_root branch=fix_r u=%.17g psi=%.17g "
          "fixed_r=%.17g lin=(%.17g,%.17g) root_z=%.17g away=%.17g seglen=%.17g\n",
          u, psi, rlin, rlin, zlin, zc, fabs(zc-zlin), seglen);
    }
    else if (tok_ordered_map_diag_enabled())
      fprintf(stderr,
        "TOK_TRACE_SAMPLE_DIAG reason=no_roots branch=fix_r u=%.17g psi=%.17g "
        "fixed_r=%.17g lin=(%.17g,%.17g) bracket=[%d,%d] brk_r=(%.17g,%.17g) brk_z=(%.17g,%.17g)\n",
        u, psi, rlin, rlin, zlin, lo, hi, tr[lo], tr[hi], tz[lo], tz[hi]);
  }
  if (!sampled) {
    double roots[16] = { 0.0 }, dRdZ[16] = { 0.0 };
    double dR[16] = { 0.0 }, dZ[16] = { 0.0 };
    int nr = gkyl_tok_geo_R_psiZ(geo, psi, zlin, 8,
      roots, dRdZ, dR, dZ);
    if (nr > 0) {
      double r_choice = tok_nearest_value(rlin, roots, nr);
      if (tok_rlin_ambiguous(rlin, roots, nr)) {
        double r_walk = 0.0;
        if (tok_trace_walk_root(geo, psi, tz[lo], tr[lo], tz[hi], tr[hi],
            zlin, &r_walk))
          r_choice = r_walk;
      }
      if (fabs(r_choice-rlin) <= max_away) {
        *r = r_choice; *z = zlin; sampled = true;
      }
      else if (tok_ordered_map_diag_enabled())
        fprintf(stderr,
          "TOK_TRACE_SAMPLE_DIAG reason=nonlocal_root branch=r_of_z u=%.17g psi=%.17g "
          "fixed_z=%.17g lin=(%.17g,%.17g) root_r=%.17g away=%.17g seglen=%.17g\n",
          u, psi, zlin, rlin, zlin, r_choice, fabs(r_choice-rlin), seglen);
    }
    else if (tok_ordered_map_diag_enabled())
      fprintf(stderr,
        "TOK_TRACE_SAMPLE_DIAG reason=no_roots branch=r_of_z u=%.17g psi=%.17g "
        "fixed_z=%.17g lin=(%.17g,%.17g) bracket=[%d,%d] brk_r=(%.17g,%.17g) brk_z=(%.17g,%.17g)\n",
        u, psi, zlin, rlin, zlin, lo, hi, tr[lo], tr[hi], tz[lo], tz[hi]);
  }
  if (!sampled) {
    // Neither axis-aligned direction has a usable root.  Project the chord
    // point onto the contour along its own normal, bounded by the bracket.
    double pr = 0.0, pz = 0.0;
    if (seglen > 0.0 && tok_psi_normal_project(geo, psi, rlin, zlin,
        max_away, 1e-9*fmax(1.0, fabs(psi)), &pr, &pz)) {
      *r = pr; *z = pz; sampled = true;
    }
  }
  if (!sampled)
    return false;
  double residual = tok_eval_psi_rz_local(geo, *r, *z)-psi;
  bool sample_ok = isfinite(*r) && isfinite(*z) && isfinite(residual) &&
    fabs(residual) <= 1e-9*fmax(1.0, fabs(psi));
  if (!sample_ok && tok_ordered_map_diag_enabled())
    fprintf(stderr,
      "TOK_TRACE_SAMPLE_DIAG reason=residual_out_of_tolerance u=%.17g psi=%.17g "
      "residual=%.17g rz=(%.17g,%.17g) lin=(%.17g,%.17g) bracket=[%d,%d] "
      "param_is_r=%d local_prefer_r_fixed=%d\n",
      u, psi, residual, *r, *z, rlin, zlin, lo, hi,
      (int) param_is_r, (int) local_prefer_r_fixed);
  return sample_ok;
}

static bool
tok_build_far_trace(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx)
{
  if (arc_ctx->far_trace_initialized)
    return true;
  if (!arc_ctx->xpt_ray_initialized) {
    fprintf(stderr,
      "TOK_ORDERED_MAP far trace failed ftype=%d reason=ray_not_initialized\n",
      inp->ftype);
    return false;
  }
  double rf = 0.0, zf = 0.0;
  if (!tok_fixed_edge_point(inp, arc_ctx->geo, arc_ctx->xpt_ray_psi0,
      &rf, &zf)) {
    fprintf(stderr,
      "TOK_ORDERED_MAP far trace failed ftype=%d reason=no_fixed_edge_point psi=%.17g\n",
      inp->ftype, arc_ctx->xpt_ray_psi0);
    return false;
  }
  int n = GKYL_MIN2(257, arc_ctx->sep_trace_capacity);
  double *rr = gkyl_malloc(sizeof(double[n]));
  double *zr = gkyl_malloc(sizeof(double[n]));
  double *rz = gkyl_malloc(sizeof(double[n]));
  double *zz = gkyl_malloc(sizeof(double[n]));
  double score_r = DBL_MAX, score_z = DBL_MAX;
  bool ok_r = tok_build_contour_candidate(arc_ctx->geo,
    arc_ctx->xpt_ray_psi0, true, rf, zf, arc_ctx->xpt_ray_r0,
    arc_ctx->xpt_ray_z0, n, rr, zr, &score_r, "far_boundary", false);
  bool ok_z = tok_build_contour_candidate(arc_ctx->geo,
    arc_ctx->xpt_ray_psi0, false, rf, zf, arc_ctx->xpt_ray_r0,
    arc_ctx->xpt_ray_z0, n, rz, zz, &score_z, "far_boundary", false);
  // Preserve both original parameterizations whenever either works.  Only a
  // trace for which both were rejected and whose fixed edge is a known
  // midplane turning point receives the endpoint-clustered R retry.
  if (!ok_r && !ok_z && tok_fixed_edge_is_midplane(inp->ftype))
    ok_r = tok_build_contour_candidate(arc_ctx->geo,
      arc_ctx->xpt_ray_psi0, true, rf, zf, arc_ctx->xpt_ray_r0,
      arc_ctx->xpt_ray_z0, n, rr, zr, &score_r, "far_boundary", true);
  if (!ok_r && !ok_z) {
    fprintf(stderr,
      "TOK_ORDERED_MAP far trace failed ftype=%d psi=%.17g fixed=(%.17g,%.17g) ray=(%.17g,%.17g)\n",
      inp->ftype, arc_ctx->xpt_ray_psi0, rf, zf, arc_ctx->xpt_ray_r0,
      arc_ctx->xpt_ray_z0);
    gkyl_free(rr); gkyl_free(zr); gkyl_free(rz); gkyl_free(zz);
    return false;
  }
  bool use_r = ok_r && (!ok_z || score_r <= score_z);
  const double *src_r = use_r ? rr : rz;
  const double *src_z = use_r ? zr : zz;
  bool fixed_first = tok_sep_fixed_edge_is_first(inp->ftype);
  arc_ctx->far_trace_s[0] = 0.0;
  for (int i=0; i<n; ++i) {
    int src = fixed_first ? i : n-1-i;
    arc_ctx->far_trace_r[i] = src_r[src];
    arc_ctx->far_trace_z[i] = src_z[src];
    if (i > 0)
      arc_ctx->far_trace_s[i] = arc_ctx->far_trace_s[i-1]
        +hypot(arc_ctx->far_trace_r[i]-arc_ctx->far_trace_r[i-1],
          arc_ctx->far_trace_z[i]-arc_ctx->far_trace_z[i-1]);
  }
  arc_ctx->far_trace_n = n;
  arc_ctx->far_trace_param_is_r = use_r;
  arc_ctx->far_trace_initialized =
    isfinite(arc_ctx->far_trace_s[n-1]) && arc_ctx->far_trace_s[n-1] > 0.0;
  gkyl_free(rr); gkyl_free(zr); gkyl_free(rz); gkyl_free(zz);
  return arc_ctx->far_trace_initialized;
}

// Number of theta samples in a map trace.  The ladder must agree with
// tok_build_current_ordered_trace exactly, and the allocation sites size the
// ladder from it, so it lives in one place.
static int
tok_ext_map_trace_nodes(const struct gkyl_tok_geo_grid_inp *inp, int capacity)
{
  return GKYL_MIN2(capacity, GKYL_MAX2(49, 4*inp->cgrid.cells[2]+1));
}

// Rung count the march STARTS from.  Rungs must be close enough in psi that
// consecutive contours are near each other -- that is what makes the projection
// unambiguous -- so the block's own radial grid is the natural starting guess.
// It is only a starting guess: how fine the ladder must actually be is a
// property of the equilibrium, not of the grid, so the march refines from here
// until it resolves the correspondence (see tok_ext_build_theta_ladder).
static int
tok_ext_ladder_rows(const struct gkyl_tok_geo_grid_inp *inp)
{
  return GKYL_MAX2(1, inp->cgrid.cells[0]);
}

// Rungs must be an integer multiple of the block's radial cell count, and this
// is a correctness requirement rather than a convenience.
//
// A node row sits at radial fraction j/rows, and the lookup interpolates the
// correspondence linearly between the two rungs bracketing (j/rows)*m.  When m
// is a multiple of rows every row lands exactly ON a rung and is never
// interpolated.  Otherwise each row picks up an interpolation perturbation set
// by where it falls inside its rung interval; adjacent rows sitting at
// different phases within their intervals receive different perturbations, and
// the rows cross radially.  It is the phase DIFFERENCE between adjacent rows
// that matters, and holding m to a multiple of rows drives it to zero
// identically.
//
// Note the monotonicity argument made at the lookup covers the other axis only:
// a convex combination of two increasing rows is increasing, so an interpolated
// row is monotone in theta.  Nothing in that argument constrains RADIAL
// ordering between rows.
static int
tok_ext_ladder_rungs(const struct gkyl_tok_geo_grid_inp *inp, int k)
{
  return GKYL_MAX2(1, k)*tok_ext_ladder_rows(inp);
}

// Largest theta motion, in units of the map trace's own node spacing, that the
// ladder may leave between adjacent rungs.
//
// w is a normalized arc coordinate and the map trace carries n nodes, so 1/(n-1)
// is the finest theta feature the table can represent at all.  If a node's w
// moves further than that from one rung to the next, then interpolating between
// those rungs places intermediate surfaces to worse than the trace's own
// resolution, and adjacent node rows can cross.  Expressing the tolerance this
// way keeps it dimensionless and derived from the grid being built rather than
// from any equilibrium- or machine-specific flux scale.
static double
tok_ext_ladder_tol(int n)
{
  return n > 1 ? 1.0/(n-1) : 1.0;
}

// Weight of the identity map blended into each marched rung.  Projection can
// legitimately pile several nodes onto one place when a contour shortens (the
// plate strike of 204951 jumps 50 mm, leaving 48.8 mm of arc with nowhere to
// go), and PAVA then collapses them to equal values.  Blending beta of the
// uniform map back in guarantees every theta gap is at least beta/(n-1) of
// uniform, so the row stays strictly ordered and no cell degenerates.
//
// Smallest theta gap a marched rung may contain, as a fraction of the uniform
// gap 1/(n-1).  Projection can legitimately pile several nodes onto one place
// when a contour shortens (204951's plate strike jumps 50 mm, leaving 48.8 mm
// of arc with nowhere to go), and PAVA then collapses them to equal values, so
// some floor is required or a cell degenerates.
//
// This replaced a global blend w = beta*u + (1-beta)*w_proj, which enforced the
// same floor but paid for it everywhere: it pulled the ENTIRE row toward
// uniform whether or not that row was pinched.  That made beta a genuine
// trade-off between the block needing the most correspondence and the one
// wanting the least, measured (by node-file count, not by absence of output):
//
//   beta    204951   205004   205079   203963 b1 thinnest-cell ratio
//   0.10    clean    clean    clean    1.07e-02
//   0.40    clean    clean    clean    4.39e-02
//   0.55    clean    FAIL     clean    6.02e-02
//   0.90    clean    FAIL     FAIL     9.87e-02
//   none    FAIL     FAIL     clean    1.06e-01
//
// 205004 bound the top of that range and 203963's cell quality bound the
// bottom, with no setting good at both ends.  Applying the floor only where it
// binds removes the conflict outright, because it caps how extreme a row may
// get without destroying the SHAPE of the correspondence -- and shape is what
// ordering depends on.  Measured, same shots, adaptive floor:
//
//   min_gap  204102  205004  205079  204951   203963 b1 thinnest-cell ratio
//   0.40     clean   clean   clean   clean    4.95e-02
//   0.60     clean   clean   clean   clean    7.07e-02
//   0.70     clean    -       -       -
//   0.90     FAIL    clean   clean   clean    9.87e-02  (no ladder: 1.06e-01)
//   0.98      -      clean    -       -
//
// This WIDENS the safe window rather than removing it: the global blend already
// failed 205004 at 0.55, whereas the floor runs to ~0.8 before 204102's
// DN_SOL_OUT_MID crosses (ftype=6, the same class as 205004).  204102 binds
// here, not 205004.  0.6 sits ~25% below that edge and still recovers 67% of
// the un-laddered cell quality, against 41% for the blend at the same safety.
// The hard limit remains min_gap -> 1: the floor then claims the whole arc,
// every gap is forced uniform, and the march carries nothing.
static double
tok_ext_ladder_min_gap(void)
{
  double f = tok_geo_trace_double_env("GKYL_TOK_EXT_LADDER_MINGAP", 0.6);
  return isfinite(f) && f > 0.0 && f < 1.0 ? f : 0.6;
}

// Which radial boundary anchors the march.  The rung the march starts from
// keeps the identity map, so it imposes its own arc-length parameterization on
// everything downstream, and the separatrix is the wrong choice for that: it
// carries the X-point corner, where arc length concentrates hard while the
// neighbouring surface is still smooth.  The mismatch folds the FIRST radial
// cell -- 204997 folded exactly one cell, at (i_psi=0, j_theta=10), where the
// row-0 to row-1 spacing varies 260x across theta, and it was insensitive to
// min_gap (identical fold at 0.2, 0.4, 0.6) because the floor is not what
// produces it.  Anchoring at the far boundary, which carries no corner, lets
// the separatrix inherit a corrected placement instead of dictating a bad one:
// it fixes 204997 and leaves 204951, 204965, 205004 and 204102 clean.
// Set GKYL_TOK_EXT_LADDER_FROM_SEP=1 to go back to separatrix-anchored.
static bool
tok_ext_ladder_from_far(void)
{
  const char *e = getenv("GKYL_TOK_EXT_LADDER_FROM_SEP");
  return !(e && e[0] != '\0' && e[0] != '0');
}

static bool
tok_ext_theta_ladder_enabled(void)
{
  const char *e = getenv("GKYL_TOK_EXT_THETA_LADDER");
  return e && e[0] != '\0' && e[0] != '0';
}

// Core surfaces are closed flux surfaces cut at a seam: the two ends of the cut
// lie on top of each other, so a nearest-point projection can hop across the
// seam onto the far end of the same contour.  That is the ambiguity the
// `closed` branch below already guards against -- but a half-domain core arc is
// not flagged closed, so it needs saying explicitly.  Measured: giving CORE_L
// and CORE_R the ladder puts 21 reversals into CORE_R (min cos -0.985) on a
// shot whose cores were previously clean.  They keep the two-point blend, which
// is both validated and correct for them, since their theta boundaries are the
// X-point rays and vary smoothly with psi.  This also means the npsi=600 core
// blocks cost nothing.
static bool
tok_ext_ladder_applies(enum gkyl_tok_geo_type ftype)
{
  switch (ftype) {
    case GKYL_GEOMETRY_TOKAMAK_CORE:
    case GKYL_GEOMETRY_TOKAMAK_CORE_L:
    case GKYL_GEOMETRY_TOKAMAK_CORE_R:
    case GKYL_GEOMETRY_TOKAMAK_IWL:
      return false;
    default:
      return true;
  }
}

// Nearest point on a polyline, searched forward from *j_lo.  u increases
// monotonically along the ladder row, so its image must too; marching the
// window forward makes that true by construction and removes the backward jump
// between the two nearby legs of a single-null contour that a global argmin
// would allow.  Returns the normalized arc position.
static double
tok_ext_project_onto_trace(const double *tr, const double *tz, const double *ts,
  int n, double rp, double zp, int *j_lo)
{
  double best_d2 = DBL_MAX, best_v = 0.0;
  int best_j = *j_lo;
  double total = ts[n-1];
  for (int j=*j_lo; j<n-1; ++j) {
    double r0 = tr[j], z0 = tz[j];
    double dr = tr[j+1]-r0, dz = tz[j+1]-z0;
    double den = dr*dr+dz*dz;
    double t = den > 0.0 ? ((rp-r0)*dr+(zp-z0)*dz)/den : 0.0;
    t = fmin(1.0, fmax(0.0, t));
    double d2 = SQ(rp-(r0+t*dr))+SQ(zp-(z0+t*dz));
    if (d2 < best_d2) {
      best_d2 = d2;
      best_j = j;
      best_v = total > 0.0 ? (ts[j]+t*(ts[j+1]-ts[j]))/total : 0.0;
    }
  }
  *j_lo = best_j;
  return fmin(1.0, fmax(0.0, best_v));
}

// Force a row strictly increasing: isotonic regression on the interior (least
// squares, rather than greedily flattening after the first reversal), then the
// identity blend, then re-anchor the two radial-boundary nodes.
static void
tok_ext_ladder_condition_row(double *w, int n)
{
  double *bmean = gkyl_malloc(sizeof(double[n]));
  int *bcount = gkyl_malloc(sizeof(int[n]));
  int nblock = 0;
  for (int i=1; i<n-1; ++i) {
    bmean[nblock] = w[i];
    bcount[nblock] = 1;
    ++nblock;
    while (nblock > 1 && bmean[nblock-2] > bmean[nblock-1]) {
      int c = bcount[nblock-2]+bcount[nblock-1];
      bmean[nblock-2] = (bcount[nblock-2]*bmean[nblock-2]
        +bcount[nblock-1]*bmean[nblock-1])/c;
      bcount[nblock-2] = c;
      --nblock;
    }
  }
  int out = 1;
  for (int b=0; b<nblock; ++b)
    for (int j=0; j<bcount[b]; ++j)
      w[out++] = bmean[b];
  gkyl_free(bmean);
  gkyl_free(bcount);
  w[0] = 0.0;
  w[n-1] = 1.0;

  // Adaptive gap floor: the L1-cheapest way to satisfy gap >= g while keeping
  // the gaps summing to 1.  Raise every deficient gap to g, then take the
  // deficit back from the gaps that have room, in proportion to their EXCESS
  // over g -- so a gap already comfortably above the floor is barely touched,
  // and a row that is nowhere pinched comes through unchanged.  That is the
  // whole point: the old global blend distorted every row by beta regardless.
  //
  // One pass is exact.  Writing nd/nf for the deficient/free counts,
  //   before: 1 = (nd*g - deficit) + (nf*g + freesum) = (n-1)*g + freesum - deficit
  //   after:  nd*g + (nf*g + freesum*scale) with scale = (freesum-deficit)/freesum
  //         = (n-1)*g + freesum - deficit = 1
  // and every free gap stays >= g because it lands at g + (nonneg)*scale.  The
  // fallback can only trigger at g*(n-1) >= 1, i.e. min_gap >= 1, which the
  // accessor excludes.
  int ngap = n-1;
  double g = tok_ext_ladder_min_gap()/ngap;
  double *d = gkyl_malloc(sizeof(double[ngap]));
  double deficit = 0.0, freesum = 0.0;
  for (int i=0; i<ngap; ++i) {
    d[i] = w[i+1]-w[i];
    if (d[i] < g) deficit += g-d[i]; else freesum += d[i]-g;
  }
  if (deficit > 0.0) {
    if (freesum > deficit) {
      double scale = (freesum-deficit)/freesum;
      for (int i=0; i<ngap; ++i)
        d[i] = d[i] < g ? g : g+(d[i]-g)*scale;
    }
    else {
      for (int i=0; i<ngap; ++i) d[i] = 1.0/ngap;
    }
    double acc = 0.0;
    for (int i=0; i<ngap; ++i) { acc += d[i]; w[i+1] = acc; }
    w[0] = 0.0;
    w[n-1] = 1.0;
  }
  gkyl_free(d);
}

// March the ladder at a given rung count into `table`, and report the largest
// theta motion between adjacent rungs so the caller can judge whether that
// count resolved the correspondence.
static bool
tok_ext_march_theta_ladder(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx, int m, int n, double *table, double *dmax_out)
{
  const int cap = arc_ctx->sep_trace_capacity;
  const double psisep = arc_ctx->geo->psisep;
  const double span = arc_ctx->xpt_ray_psi0-psisep;
  *dmax_out = HUGE_VAL;
  if (n < 2 || m < 1 || !isfinite(span) || span == 0.0)
    return false;

  double *scratch = gkyl_malloc(sizeof(double[6*cap]));
  double *pr = scratch,      *pz = scratch+cap,   *ps = scratch+2*cap;
  double *cr = scratch+3*cap, *cz = scratch+4*cap, *cs = scratch+5*cap;
  double *wprev = gkyl_malloc(sizeof(double[n]));
  double *wcur = gkyl_malloc(sizeof(double[n]));

  const bool from_far = tok_ext_ladder_from_far();
  const int k_seed = from_far ? m : 0;
  int pn = from_far ? arc_ctx->far_trace_n : arc_ctx->sep_trace_n;
  bool pparam = from_far ? arc_ctx->far_trace_param_is_r
    : arc_ctx->sep_trace_param_is_r;
  double ppsi = from_far ? arc_ctx->xpt_ray_psi0 : psisep;
  const double *sr = from_far ? arc_ctx->far_trace_r : arc_ctx->sep_trace_r;
  const double *sz = from_far ? arc_ctx->far_trace_z : arc_ctx->sep_trace_z;
  const double *ss = from_far ? arc_ctx->far_trace_s : arc_ctx->sep_trace_s;
  for (int i=0; i<pn; ++i) { pr[i] = sr[i]; pz[i] = sz[i]; ps[i] = ss[i]; }
  for (int i=0; i<n; ++i) {
    wprev[i] = i/(double) (n-1);
    table[(size_t) k_seed*n+i] = wprev[i];
  }

  bool ok = true;
  double dmax = 0.0;
  for (int step=1; step<=m && ok; ++step) {
    int k = from_far ? m-step : step;
    double psi_k = psisep+span*(k/(double) m);
    int cn = 0;
    bool cparam = false, cclosed = false;
    if (!tok_ext_build_domain_trace(inp, arc_ctx, psi_k, k == 0,
        cr, cz, cs, &cn, &cparam, &cclosed) || cn < 2) {
      fprintf(stderr,
        "TOK_EXT_LADDER reason=rung_trace_failed ftype=%d k=%d psi=%.17g\n",
        inp->ftype, k, psi_k);
      ok = false;
      break;
    }
    int j_lo = 0;
    for (int i=0; i<n; ++i) {
      double rp = 0.0, zp = 0.0;
      if (!tok_trace_sample(arc_ctx->geo, ppsi, pr, pz, ps, pn, pparam,
          wprev[i], &rp, &zp)) {
        fprintf(stderr,
          "TOK_EXT_LADDER reason=prev_sample_failed ftype=%d k=%d i=%d w=%.17g\n",
          inp->ftype, k, i, wprev[i]);
        ok = false;
        break;
      }
      wcur[i] = tok_ext_project_onto_trace(cr, cz, cs, cn, rp, zp, &j_lo);
    }
    if (!ok) break;
    tok_ext_ladder_condition_row(wcur, n);
    for (int i=1; i<n; ++i)
      if (!(wcur[i] > wcur[i-1])) {
        fprintf(stderr,
          "TOK_EXT_LADDER reason=nonmonotonic ftype=%d k=%d i=%d prev=%.17g curr=%.17g\n",
          inp->ftype, k, i, wcur[i-1], wcur[i]);
        ok = false;
        break;
      }
    if (!ok) break;
    for (int i=0; i<n; ++i)
      table[(size_t) k*n+i] = wcur[i];
    for (int i=0; i<n; ++i)
      dmax = fmax(dmax, fabs(wcur[i]-wprev[i]));
    // This rung becomes the next one's reference.
    for (int i=0; i<cn; ++i) { pr[i] = cr[i]; pz[i] = cz[i]; ps[i] = cs[i]; }
    pn = cn; pparam = cparam; ppsi = psi_k;
    for (int i=0; i<n; ++i) wprev[i] = wcur[i];
  }

  gkyl_free(scratch);
  gkyl_free(wprev);
  gkyl_free(wcur);
  if (!ok) return false;
  *dmax_out = dmax;
  return true;
}

// Build the marched ladder, refining in psi until it actually resolves the
// correspondence.  Failure is not fatal: the caller falls back to the two-point
// blend, which is what shipped before.
//
// How many rungs a block needs is set by how fast its contours change with psi,
// which is a property of the equilibrium and not of the grid.  A block whose
// trace length steps sharply between adjacent surfaces needs several times the
// rungs its radial cell count alone would give, while most blocks are resolved
// by that count.  A single fixed count cannot serve both, and any fixed
// multiple of the grid size is only a new constant for a different equilibrium
// to defeat.  So march, measure, and refine until the ladder resolves itself,
// in steps that keep the rungs aligned with the rows.
//
// Two things stop the refinement, and neither is a tuned threshold: the
// tolerance is met, or doubling stops reducing the motion -- which is how a
// genuine geometric discontinuity announces itself, since no psi resolution can
// smooth a contour that really does jump.  Refining further would then only
// cost time, and the min_gap floor already converts such a jump into a stretched
// row rather than a crossed one.  Rungs are also capped by the trace capacity,
// beyond which the rungs would be finer than the traces they are built from.
static bool
tok_ext_build_theta_ladder(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx)
{
  const int cap = arc_ctx->sep_trace_capacity;
  const int n = tok_ext_map_trace_nodes(inp, cap);
  if (n < 2)
    return false;
  const double tol = tok_ext_ladder_tol(n);

  double *table = NULL;
  int m_done = 0;
  double dmax_prev = HUGE_VAL;
  int refinements = 0;
  // An explicit override pins the count: it exists to A/B a specific ladder
  // resolution, which refining away from would defeat.
  const int pinned = (int) tok_geo_trace_double_env("GKYL_TOK_EXT_LADDER_RUNGS", 0.0);

  for (int k = 1; ; ) {
    const int m = pinned > 0 ? pinned : tok_ext_ladder_rungs(inp, k);
    if (m > cap)
      break;
    double *cand = gkyl_malloc(sizeof(double)*(size_t) (m+1)*n);
    double dmax = HUGE_VAL;
    if (!tok_ext_march_theta_ladder(inp, arc_ctx, m, n, cand, &dmax)) {
      gkyl_free(cand);
      break;
    }
    gkyl_free(table);
    table = cand;
    m_done = m;
    if (pinned > 0 || dmax <= tol)
      break;
    // Not resolved -- but if refining did not reduce the motion, the contour
    // genuinely jumps and no rung spacing will smooth it.
    if (!(dmax < dmax_prev))
      break;
    dmax_prev = dmax;
    // Adjacent rungs are one rung spacing apart in psi, so to leading order the
    // motion between them falls off like 1/m, which makes the multiplier that
    // meets the tolerance predictable from the one just measured.  Each march
    // costs a full set of contour traces, so jumping straight there is worth
    // doing -- but only as far as the estimate can be trusted.  Far from that
    // asymptotic regime the ratio badly overestimates, so take the estimate
    // when it asks for less than a doubling and fall back to doubling when it
    // asks for more: the search then converges as fast as doubling in the worst
    // case, faster when the estimate is good, and can never skip a multiplier
    // that would have worked.
    //
    // Refining is safe here in a way it is not for an unconstrained count:
    // every candidate keeps the rungs aligned with the rows, so overshooting
    // costs time and nothing else.
    double predicted = k*(dmax/tol);
    int k_next = predicted < (double) (2*k) ? (int) ceil(predicted) : 2*k;
    k = GKYL_MAX2(k_next, k+1);
    ++refinements;
  }

  if (!table)
    return false;
  gkyl_free(arc_ctx->ext_ladder_w);
  arc_ctx->ext_ladder_w = table;
  arc_ctx->ext_ladder_m = m_done;
  arc_ctx->ext_ladder_n = n;
  arc_ctx->ext_ladder_initialized = true;
  if (tok_ordered_map_diag_enabled())
    fprintf(stderr,
      "TOK_EXT_LADDER built ftype=%d rungs=%d rows=%d k=%d nodes=%d "
      "min_gap=%.4g refinements=%d tol=%.4g\n",
      inp->ftype, m_done, tok_ext_ladder_rows(inp),
      m_done/GKYL_MAX2(1, tok_ext_ladder_rows(inp)), n,
      tok_ext_ladder_min_gap(), refinements, tol);
  return true;
}

// w(u_i, psi) by linear interpolation between the two bracketing rungs.  A
// convex combination of two increasing rows is increasing, so the interpolated
// row is monotone in THETA for every radial fraction.  That is the only axis
// that argument covers: it says nothing about RADIAL ordering between adjacent
// rows, which is why the rung count is held to a multiple of the row count --
// see tok_ext_ladder_rungs.
static double
tok_ext_ladder_w(const struct arc_length_ctx *arc_ctx, double rf, int i)
{
  int m = arc_ctx->ext_ladder_m, n = arc_ctx->ext_ladder_n;
  double x = fmin(1.0, fmax(0.0, rf))*m;
  int k = GKYL_MIN2(m-1, GKYL_MAX2(0, (int) floor(x)));
  double t = x-k;
  const double *w0 = arc_ctx->ext_ladder_w+(size_t) k*n;
  const double *w1 = w0+n;
  return (1.0-t)*w0[i]+t*w1[i];
}

static bool
tok_build_trace_correspondence(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx)
{
  const bool extended = tok_ext_construction(inp);
  if (arc_ctx->ordered_boundaries_initialized)
    return true;
  bool closed = false;
  if (extended) {
    bool sep_closed = false, far_closed = false;
    bool sep_built_generic = false, sep_built_here = false;
    if (!arc_ctx->sep_trace_initialized) {
      if (!tok_ext_build_domain_trace(inp, arc_ctx, arc_ctx->geo->psisep,
          true, arc_ctx->sep_trace_r, arc_ctx->sep_trace_z,
          arc_ctx->sep_trace_s, &arc_ctx->sep_trace_n,
          &arc_ctx->sep_trace_param_is_r, &sep_closed)) {
        fprintf(stderr,
          "TOK_TRACE_CORR reason=sep_domain_trace_failed ftype=%d psisep=%.17g\n",
          inp->ftype, arc_ctx->geo->psisep);
        return false;
      }
      arc_ctx->sep_trace_initialized = true;
      sep_built_generic = arc_ctx->ext_last_trace_used_generic;
      sep_built_here = true;
    }
    else {
      struct tok_ext_topology top;
      if (!tok_ext_topology_from_ftype(inp->ftype, inp->half_domain, &top)) {
        fprintf(stderr,
          "TOK_TRACE_CORR reason=sep_topology_unknown ftype=%d\n", inp->ftype);
        return false;
      }
      sep_closed = top.closed;
    }
    if (!arc_ctx->far_trace_initialized) {
      if (!tok_ext_build_domain_trace(inp, arc_ctx,
          arc_ctx->xpt_ray_psi0, false, arc_ctx->far_trace_r,
          arc_ctx->far_trace_z, arc_ctx->far_trace_s,
          &arc_ctx->far_trace_n, &arc_ctx->far_trace_param_is_r,
          &far_closed)) {
        fprintf(stderr,
          "TOK_TRACE_CORR reason=far_domain_trace_failed ftype=%d ray_psi0=%.17g\n",
          inp->ftype, arc_ctx->xpt_ray_psi0);
        return false;
      }
      arc_ctx->far_trace_initialized = true;
      // The far surface can be the one that trips the block onto the generic
      // route, in which case the separatrix trace above was already built with
      // the topology route.  Rebuild it so both radial boundaries -- and every
      // interior surface, which now follows the flag -- share one route.
      if (sep_built_here && !sep_built_generic &&
          arc_ctx->ext_force_generic_route && !sep_closed) {
        if (!tok_ext_build_domain_trace(inp, arc_ctx, arc_ctx->geo->psisep,
            true, arc_ctx->sep_trace_r, arc_ctx->sep_trace_z,
            arc_ctx->sep_trace_s, &arc_ctx->sep_trace_n,
            &arc_ctx->sep_trace_param_is_r, &sep_closed)) {
          fprintf(stderr,
            "TOK_TRACE_CORR reason=sep_regeneric_failed ftype=%d psisep=%.17g\n",
            inp->ftype, arc_ctx->geo->psisep);
          return false;
        }
        fprintf(stderr,
          "TOK_EXT_TRACE sep_rebuilt_generic ftype=%d psisep=%.17g\n",
          inp->ftype, arc_ctx->geo->psisep);
      }
    }
    else {
      struct tok_ext_topology top;
      if (!tok_ext_topology_from_ftype(inp->ftype, inp->half_domain, &top)) {
        fprintf(stderr,
          "TOK_TRACE_CORR reason=far_topology_unknown ftype=%d\n", inp->ftype);
        return false;
      }
      far_closed = top.closed;
    }
    if (sep_closed != far_closed) {
      fprintf(stderr,
        "TOK_TRACE_CORR reason=closedness_mismatch ftype=%d sep_closed=%d far_closed=%d\n",
        inp->ftype, (int) sep_closed, (int) far_closed);
      return false;
    }
    closed = sep_closed;
  }
  else {
    if (!arc_ctx->sep_trace_initialized) {
      double rf = 0.0, zf = 0.0;
      if (!tok_fixed_edge_point(inp, arc_ctx->geo, arc_ctx->geo->psisep,
          &rf, &zf)) {
        fprintf(stderr,
          "TOK_TRACE_CORR reason=fixed_edge_point_failed ftype=%d psisep=%.17g\n",
          inp->ftype, arc_ctx->geo->psisep);
        return false;
      }
      tok_build_sep_trace(inp, arc_ctx, zf, rf);
    }
    if (!tok_build_far_trace(inp, arc_ctx)) {
      fprintf(stderr,
        "TOK_TRACE_CORR reason=far_trace_failed ftype=%d ray_psi0=%.17g\n",
        inp->ftype, arc_ctx->xpt_ray_psi0);
      return false;
    }
  }

  int n = GKYL_MIN2(arc_ctx->sep_trace_n, arc_ctx->far_trace_n);
  if (n < 2) {
    fprintf(stderr,
      "TOK_TRACE_CORR reason=too_few_trace_points ftype=%d n=%d sep_n=%d far_n=%d\n",
      inp->ftype, n, arc_ctx->sep_trace_n, arc_ctx->far_trace_n);
    return false;
  }
  // The extended full-domain traces already share topology-matched endpoints
  // and direction.  Equal normalized contour length therefore supplies a
  // one-to-one coordinate on every radial surface.  In particular, it avoids
  // radial folds caused when nearest-point projection jumps between the two
  // nearby legs of a single-null contour.  A closed core also needs this
  // identity map because its duplicated seam makes nearest projection
  // ambiguous.  Keep the original projection/PAVA correspondence unchanged
  // for the established half-domain path.
  // A closed core keeps the identity map unconditionally: its duplicated seam
  // makes nearest projection ambiguous.  An open extended block keeps it only
  // while the projection correspondence is switched off.
  const bool ext_open = extended && !closed;
  const bool ext_ladder = ext_open && tok_ext_theta_ladder_enabled() &&
    tok_ext_ladder_applies(inp->ftype);
  if (extended || closed) {
    for (int i=0; i<n; ++i)
      arc_ctx->trace_corr_v[i] = i/(double) (n-1);
    arc_ctx->trace_corr_n = n;
    // Both radial boundary traces exist and the block's route has settled, so
    // this is the one point where the ladder can be marched without perturbing
    // route selection.  A failure here is not fatal -- the two-point blend
    // above remains valid, it just cannot resolve an interior discontinuity.
    if (ext_ladder && !arc_ctx->ext_ladder_initialized &&
        !arc_ctx->ext_ladder_failed) {
      if (!tok_ext_build_theta_ladder(inp, arc_ctx)) {
        arc_ctx->ext_ladder_failed = true;
        fprintf(stderr,
          "TOK_EXT_LADDER reason=build_failed ftype=%d falling back to blend\n",
          inp->ftype);
      }
    }
    arc_ctx->ordered_boundaries_initialized = true;
    return true;
  }
  arc_ctx->trace_corr_v[0] = 0.0;
  for (int i=1; i<n-1; ++i) {
    double u = i/(double) (n-1), rs = 0.0, zs = 0.0;
    if (!tok_trace_sample(arc_ctx->geo, arc_ctx->geo->psisep,
        arc_ctx->sep_trace_r, arc_ctx->sep_trace_z,
        arc_ctx->sep_trace_s, arc_ctx->sep_trace_n,
        arc_ctx->sep_trace_param_is_r, u, &rs, &zs)) {
      fprintf(stderr,
        "TOK_TRACE_CORR reason=sep_trace_sample_failed ftype=%d i=%d n=%d u=%.17g sep_n=%d\n",
        inp->ftype, i, n, u, arc_ctx->sep_trace_n);
      return false;
    }
    double best_d2 = DBL_MAX, best_v = 0.0;
    for (int j=0; j<arc_ctx->far_trace_n-1; ++j) {
      double r0 = arc_ctx->far_trace_r[j], z0 = arc_ctx->far_trace_z[j];
      double dr = arc_ctx->far_trace_r[j+1]-r0;
      double dz = arc_ctx->far_trace_z[j+1]-z0;
      double den = dr*dr+dz*dz;
      double t = den > 0.0 ? ((rs-r0)*dr+(zs-z0)*dz)/den : 0.0;
      t = fmin(1.0, fmax(0.0, t));
      double rp = r0+t*dr, zp = z0+t*dz;
      double d2 = SQ(rs-rp)+SQ(zs-zp);
      if (d2 < best_d2) {
        best_d2 = d2;
        best_v = (arc_ctx->far_trace_s[j]
          +t*(arc_ctx->far_trace_s[j+1]-arc_ctx->far_trace_s[j]))
          /arc_ctx->far_trace_s[arc_ctx->far_trace_n-1];
      }
    }
    arc_ctx->trace_corr_v[i] = best_v;
  }
  arc_ctx->trace_corr_v[n-1] = 1.0;

  // Isotonic regression (increasing PAVA) finds the least-squares monotone
  // correspondence, rather than greedily flattening every value after the
  // first local reversal.  Endpoints are already bounded by every projection,
  // so applying PAVA to the interior preserves the exact 0 and 1 anchors.
  int nblock = 0;
  double *block_mean = gkyl_malloc(sizeof(double[n]));
  int *block_count = gkyl_malloc(sizeof(int[n]));
  for (int i=1; i<n-1; ++i) {
    block_mean[nblock] = arc_ctx->trace_corr_v[i];
    block_count[nblock] = 1;
    ++nblock;
    while (nblock > 1 && block_mean[nblock-2] > block_mean[nblock-1]) {
      int count = block_count[nblock-2]+block_count[nblock-1];
      block_mean[nblock-2] =
        (block_count[nblock-2]*block_mean[nblock-2]
          +block_count[nblock-1]*block_mean[nblock-1])/count;
      block_count[nblock-2] = count;
      --nblock;
    }
  }
  int out = 1;
  for (int b=0; b<nblock; ++b)
    for (int j=0; j<block_count[b]; ++j)
      arc_ctx->trace_corr_v[out++] = block_mean[b];
  gkyl_free(block_mean);
  gkyl_free(block_count);
  // Nearest-point projection can legitimately collapse a finite interval of
  // the separatrix onto one far-boundary location.  Blend in a small identity
  // map so g remains well-conditioned and strictly order preserving without
  // discarding the nearest-normal correspondence.  Shot 203730 remains
  // fold-free over a broad range around this conservative one-percent blend.
  const double identity_fraction = tok_trace_corr_identity_fraction();
  for (int i=1; i<n-1; ++i) {
    double u = i/(double) (n-1);
    arc_ctx->trace_corr_v[i] = identity_fraction*u
      +(1.0-identity_fraction)*arc_ctx->trace_corr_v[i];
  }
  for (int i=1; i<n; ++i)
    if (!(arc_ctx->trace_corr_v[i] > arc_ctx->trace_corr_v[i-1])) {
      fprintf(stderr,
        "TOK_TRACE_CORR reason=nonmonotonic_correspondence ftype=%d i=%d n=%d v_prev=%.17g v_curr=%.17g\n",
        inp->ftype, i, n, arc_ctx->trace_corr_v[i-1], arc_ctx->trace_corr_v[i]);
      return false;
    }
  arc_ctx->trace_corr_n = n;
  arc_ctx->ordered_boundaries_initialized = true;
  return true;
}

static double
tok_trace_correspondence(const struct arc_length_ctx *arc_ctx, double u)
{
  const double endpoint_tol = 256.0*DBL_EPSILON;
  if (u <= endpoint_tol) return 0.0;
  if (u >= 1.0-endpoint_tol) return 1.0;
  double x = u*(arc_ctx->trace_corr_n-1);
  int i = GKYL_MIN2(arc_ctx->trace_corr_n-2, (int) floor(x));
  double w = x-i;
  return arc_ctx->trace_corr_v[i]
    +w*(arc_ctx->trace_corr_v[i+1]-arc_ctx->trace_corr_v[i]);
}

static bool
tok_ordered_chord_point(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx, double psi, double u, double *r, double *z)
{
  if (!tok_build_trace_correspondence(inp, arc_ctx))
    return false;
  u = fmin(1.0, fmax(0.0, u));
  double v = tok_trace_correspondence(arc_ctx, u);
  double ra = 0.0, za = 0.0, rb = 0.0, zb = 0.0;
  if (!tok_trace_sample(arc_ctx->geo, arc_ctx->geo->psisep,
      arc_ctx->sep_trace_r, arc_ctx->sep_trace_z,
      arc_ctx->sep_trace_s, arc_ctx->sep_trace_n,
      arc_ctx->sep_trace_param_is_r, u, &ra, &za) ||
      !tok_trace_sample(arc_ctx->geo, arc_ctx->xpt_ray_psi0,
        arc_ctx->far_trace_r, arc_ctx->far_trace_z,
        arc_ctx->far_trace_s, arc_ctx->far_trace_n,
        arc_ctx->far_trace_param_is_r, v, &rb, &zb))
    return false;
  double scale = fmax(1.0, fmax(fabs(psi), fabs(arc_ctx->geo->psisep)));
  if (tok_geo_same_flux(psi, arc_ctx->geo->psisep)) {
    *r = ra; *z = za; return true;
  }
  if (tok_geo_same_flux(psi, arc_ctx->xpt_ray_psi0)) {
    *r = rb; *z = zb; return true;
  }
  double fa = tok_eval_psi_rz_local(arc_ctx->geo, ra, za)-psi;
  double fb = tok_eval_psi_rz_local(arc_ctx->geo, rb, zb)-psi;
  if (!isfinite(fa) || !isfinite(fb) || fa*fb > 0.0) {
    fprintf(stderr,
      "TOK_ORDERED_MAP unbracketed chord ftype=%d psi=%.17g u=%.17g fa=%.17g fb=%.17g\n",
      inp->ftype, psi, u, fa, fb);
    return false;
  }
  // NOTE: psi is not always monotone along the chord -- where the chord leaves
  // the X point nearly tangent to a separatrix branch it can rise past psi, dip
  // back below and rise again, so psi=psi_curr has three roots and which one
  // this bisection converges to depends on the bracketing accident rather than
  // on a stated rule. Selecting the FIRST root instead was tried and is WRONG
  // here: on 204502's DN_SOL_IN_LO it puts the seam node of the first flux
  // surface at 0.031 m from the X point while the rest of that radial row sits
  // at ~0.15 m, because the near-X-point root lies on the piece of the contour
  // that wraps the saddle rather than on the arc that continues to the plate.
  // The outermost root, which this bisection happens to find, is the smooth
  // continuation (0.023 m from its theta neighbour, versus 0.119 m for the
  // first root). Leave the selection alone unless a case is measured where it
  // lands on a root that is not the neighbour-consistent one.
  double slo = 0.0, shi = 1.0, flo = fa;
  for (int k=0; k<44; ++k) {
    double smid = 0.5*(slo+shi);
    double rm = ra+smid*(rb-ra), zm = za+smid*(zb-za);
    double fm = tok_eval_psi_rz_local(arc_ctx->geo, rm, zm)-psi;
    if (!isfinite(fm)) return false;
    if (flo*fm <= 0.0)
      shi = smid;
    else { slo = smid; flo = fm; }
  }
  double s = 0.5*(slo+shi);
  *r = ra+s*(rb-ra); *z = za+s*(zb-za);
  double residual = tok_eval_psi_rz_local(arc_ctx->geo, *r, *z)-psi;
  return isfinite(*r) && isfinite(*z) && isfinite(residual) &&
    fabs(residual) <= 1e-9*scale;
}

static bool
tok_xpt_seam_endpoint(enum gkyl_tok_geo_type ftype, double u)
{
  const double endpoint_tol = 256.0*DBL_EPSILON;
  return tok_sep_fixed_edge_is_first(ftype)
    ? u >= 1.0-endpoint_tol : u <= endpoint_tol;
}

static bool
tok_eval_psi_grad_rz_local(const struct gkyl_tok_geo *geo,
  double R, double Z, double *dpsidR, double *dpsidZ)
{
  if (geo->use_cubics) {
    double xn[2] = { R, Z }, out[3] = { 0.0 };
    geo->efit->evf->eval_cubic_wgrad(0.0, xn, out, geo->efit->evf->ctx);
    *dpsidR = out[1]; *dpsidZ = out[2];
    return isfinite(*dpsidR) && isfinite(*dpsidZ);
  }
  int idx[2];
  idx[0] = GKYL_MIN2(geo->rzlocal.upper[0], GKYL_MAX2(geo->rzlocal.lower[0],
    geo->rzlocal.lower[0]
      +(int) floor((R-geo->rzgrid.lower[0])/geo->rzgrid.dx[0])));
  idx[1] = GKYL_MIN2(geo->rzlocal.upper[1], GKYL_MAX2(geo->rzlocal.lower[1],
    geo->rzlocal.lower[1]
      +(int) floor((Z-geo->rzgrid.lower[1])/geo->rzgrid.dx[1])));
  long loc = gkyl_range_idx(&geo->rzlocal, idx);
  const double *p = gkyl_array_cfetch(geo->psiRZ, loc);
  double xc[2];
  gkyl_rect_grid_cell_center(&geo->rzgrid, idx, xc);
  double x = (R-xc[0])/(0.5*geo->rzgrid.dx[0]);
  double y = (Z-xc[1])/(0.5*geo->rzgrid.dx[1]);
  if (geo->efit->rzbasis.poly_order == 1) {
    *dpsidR = (1.5*p[3]*y+0.8660254037844386*p[1])
      *2.0/geo->rzgrid.dx[0];
    *dpsidZ = (1.5*p[3]*x+0.8660254037844386*p[2])
      *2.0/geo->rzgrid.dx[1];
  }
  else {
    *dpsidR = (5.625*p[8]*(2.0*x*y*y-0.6666666666666666*x)
      +2.904737509655563*p[7]*(y*y-0.3333333333333333)
      +5.809475019311126*p[6]*x*y+1.5*p[3]*y
      +3.354101966249684*p[4]*x+0.8660254037844386*p[1])
      *2.0/geo->rzgrid.dx[0];
    *dpsidZ = (5.625*p[8]*(2.0*x*x*y-0.6666666666666666*y)
      +5.809475019311126*p[7]*x*y+3.354101966249684*p[5]*y
      +2.904737509655563*p[6]*(x*x-0.3333333333333333)
      +1.5*p[3]*x+0.8660254037844386*p[2])
      *2.0/geo->rzgrid.dx[1];
  }
  return isfinite(*dpsidR) && isfinite(*dpsidZ);
}

// Project a continuation predictor to the exact requested DG/cubic flux
// surface. Holding the coordinate that changes most along the contour keeps
// the root solve well conditioned at either an R or Z turning point.
static bool
tok_project_xpt_seam_candidate(const struct gkyl_tok_geo *geo, double psi,
  double rprev, double zprev, double rpred, double zpred,
  double tangent_r, double tangent_z, double step,
  double *r, double *z)
{
  const struct gkyl_rect_grid *grid = geo->use_cubics
    ? &geo->rzgrid_cubic : &geo->rzgrid;
  if (rpred < grid->lower[0] || rpred > grid->upper[0] ||
      zpred < grid->lower[1] || zpred > grid->upper[1])
    return false;

  double best_distance = DBL_MAX, second_distance = DBL_MAX;
  if (fabs(tangent_r) >= fabs(tangent_z)) {
    double roots[32] = { 0.0 };
    int nr = tok_geo_Z_psiR(geo, psi, rpred, 32, roots);
    if (nr <= 0)
      return false;
    int best = 0;
    for (int i=0; i<nr; ++i) {
      double distance = fabs(roots[i]-zpred);
      if (distance < best_distance) {
        second_distance = best_distance;
        best_distance = distance;
        best = i;
      }
      else if (distance < second_distance)
        second_distance = distance;
    }
    double ambiguity_scale = fmax(grid->dx[1], fabs(step));
    if (nr > 1 && second_distance-best_distance <=
        1e-8*fmax(1.0, ambiguity_scale))
      return false;
    *r = rpred; *z = roots[best];
  }
  else {
    double roots[32] = { 0.0 }, dRdZ[32] = { 0.0 };
    double dR[32] = { 0.0 }, dZ[32] = { 0.0 };
    int nr = gkyl_tok_geo_R_psiZ(geo, psi, zpred, 32,
      roots, dRdZ, dR, dZ);
    if (nr <= 0)
      return false;
    int best = 0;
    for (int i=0; i<nr; ++i) {
      double distance = fabs(roots[i]-rpred);
      if (distance < best_distance) {
        second_distance = best_distance;
        best_distance = distance;
        best = i;
      }
      else if (distance < second_distance)
        second_distance = distance;
    }
    double ambiguity_scale = fmax(grid->dx[0], fabs(step));
    if (nr > 1 && second_distance-best_distance <=
        1e-8*fmax(1.0, ambiguity_scale))
      return false;
    *r = roots[best]; *z = zpred;
  }

  double dr = *r-rprev, dz = *z-zprev;
  double forward = step*(dr*tangent_r+dz*tangent_z);
  double distance = hypot(dr, dz);
  double residual = tok_eval_psi_rz_local(geo, *r, *z)-psi;
  double flux_scale = fmax(1.0, fabs(psi));
  return isfinite(*r) && isfinite(*z) && isfinite(distance) &&
    isfinite(residual) && forward > 0.0 &&
    distance <= 2.5*fabs(step)+1e-12 &&
    fabs(residual) <= 1e-10*flux_scale;
}

static bool
tok_displace_xpt_seam_on_flux(const struct gkyl_tok_geo *geo, double psi,
  double delta_s, double *r, double *z, double *realized_out)
{
  *realized_out = 0.0;
  if (delta_s == 0.0)
    return true;
  const struct gkyl_rect_grid *grid = geo->use_cubics
    ? &geo->rzgrid_cubic : &geo->rzgrid;
  double max_step = 0.05*fmin(grid->dx[0], grid->dx[1]);
  if (!(max_step > 0.0) || !isfinite(max_step))
    return false;
  double required_steps = ceil(fabs(delta_s)/max_step);
  if (!isfinite(required_steps) || required_steps > 4096.0)
    return false;
  int nstep = GKYL_MAX2(16, (int) required_steps);
  double rstart = *r, zstart = *z;
  double arc_length_error = DBL_MAX, realized = 0.0;
  // Combined absolute+relative acceptance. Near the B1(q) taper's
  // endpoints the requested per-node delta_s shrinks toward zero, so a
  // purely relative criterion against fabs(delta_s) is unreachable no
  // matter how many times nstep is doubled; the absolute floor keeps that
  // from being amplified into spurious rejections. Measured empirically:
  // even at the nstep cap (4096) the achievable absolute arc-length error
  // plateaus around 1e-8 m regardless of delta_s magnitude (set by the
  // underlying per-step Newton projection's own precision, not by
  // discretization), so 1e-10 is unreachable in both framings -- 1e-7
  // comfortably covers the observed error with margin while still being
  // far tighter than anything physically meaningful at this length scale.
  const double arc_length_rel_tol = 1e-7;
  const double arc_length_abs_tol = 1e-7;
  for (;;) {
    *r = rstart; *z = zstart;
    double step = delta_s/nstep;
    realized = 0.0;

    for (int i=0; i<nstep; ++i) {
      double grad_r = 0.0, grad_z = 0.0;
      if (!tok_eval_psi_grad_rz_local(geo, *r, *z, &grad_r, &grad_z))
        return false;
      double grad = hypot(grad_r, grad_z);
      if (!(grad > 1e-14))
        return false;
      double tangent_r = -grad_z/grad, tangent_z = grad_r/grad;

      double rmid = *r+0.5*step*tangent_r;
      double zmid = *z+0.5*step*tangent_z;
      if (!tok_project_xpt_seam_candidate(geo, psi, *r, *z, rmid, zmid,
          tangent_r, tangent_z, 0.5*step, &rmid, &zmid)) {
        if (tok_ordered_map_diag_enabled())
          fprintf(stderr,
            "TOK_XPT_SEAM_DIAG reason=midpoint_projection_failed psi=%.17g step=%d/%d R=%.17g Z=%.17g delta_s=%.17g\n",
            psi, i, nstep, *r, *z, delta_s);
        return false;
      }
      if (!tok_eval_psi_grad_rz_local(geo, rmid, zmid, &grad_r, &grad_z))
        return false;
      grad = hypot(grad_r, grad_z);
      if (!(grad > 1e-14))
        return false;
      tangent_r = -grad_z/grad; tangent_z = grad_r/grad;

      double rnext = *r+step*tangent_r;
      double znext = *z+step*tangent_z;
      if (!tok_project_xpt_seam_candidate(geo, psi, *r, *z, rnext, znext,
          tangent_r, tangent_z, step, &rnext, &znext)) {
        if (tok_ordered_map_diag_enabled())
          fprintf(stderr,
            "TOK_XPT_SEAM_DIAG reason=full_projection_failed psi=%.17g step=%d/%d R=%.17g Z=%.17g delta_s=%.17g\n",
            psi, i, nstep, *r, *z, delta_s);
        return false;
      }
      realized += hypot(rnext-*r, znext-*z);
      *r = rnext; *z = znext;
    }
    double absolute_error = fabs(realized-fabs(delta_s));
    double error_tolerance = fmax(arc_length_rel_tol*fabs(delta_s),
      arc_length_abs_tol);
    arc_length_error = absolute_error;
    if (isfinite(absolute_error) && absolute_error <= error_tolerance) {
      *realized_out = realized;
      return true;
    }
    if (tok_ordered_map_diag_enabled())
      fprintf(stderr,
        "TOK_XPT_SEAM_DIAG reason=arc_length_refine psi=%.17g requested=%.17g realized=%.17g signed_absolute_error=%.17g error_tolerance=%.17g nstep=%d step_size=%.17g\n",
        psi, delta_s, realized, realized-fabs(delta_s), error_tolerance,
        nstep, fabs(step));
    if (nstep > 4096/2)
      break;
    nstep *= 2;
  }
  *r = rstart; *z = zstart;
  if (tok_ordered_map_diag_enabled())
    fprintf(stderr,
      "TOK_XPT_SEAM_DIAG reason=arc_length_not_converged psi=%.17g requested=%.17g realized=%.17g absolute_error=%.17g nstep=%d\n",
      psi, delta_s, realized, arc_length_error, nstep);
  return false;
}

// B1(q)=4q(1-q) is smooth, has unit peak, and uses explicit endpoint
// branches so both physical anchors remain exactly fixed.
static bool
tok_xpt_seam_delta_s(const struct gkyl_tok_geo_grid_inp *inp,
  const struct arc_length_ctx *arc_ctx, double psi, double *q,
  double *delta_s)
{
  double span = arc_ctx->xpt_ray_psi0-arc_ctx->geo->psisep;
  double scale = fmax(1.0, fmax(fabs(arc_ctx->xpt_ray_psi0),
    fabs(arc_ctx->geo->psisep)));
  if (!isfinite(span) || fabs(span) <= 256.0*DBL_EPSILON*scale)
    return false;
  if (tok_geo_same_flux(psi, arc_ctx->geo->psisep)) {
    *q = 0.0; *delta_s = 0.0;
    return true;
  }
  if (tok_geo_same_flux(psi, arc_ctx->xpt_ray_psi0)) {
    *q = 1.0; *delta_s = 0.0;
    return true;
  }
  *q = (psi-arc_ctx->geo->psisep)/span;
  if (*q < -1e-8 || *q > 1.0+1e-8)
    return false;
  *q = fmin(1.0, fmax(0.0, *q));
  if (*q == 0.0 || *q == 1.0) {
    *delta_s = 0.0;
    return true;
  }
  double coefficient = inp->relaxed_xpt_seam_delta_s_coeff;
  double bound = inp->relaxed_xpt_seam_delta_s_bound;
  if (coefficient == 0.0) {
    *delta_s = 0.0;
    return true;
  }
  if (!isfinite(coefficient) || !isfinite(bound) || !(bound > 0.0) ||
      fabs(coefficient) > bound*(1.0+64.0*DBL_EPSILON))
    return false;
  *delta_s = coefficient*4.0*(*q)*(1.0-*q);
  return isfinite(*delta_s) && fabs(*delta_s) <=
    bound*(1.0+64.0*DBL_EPSILON);
}

// Relax an already-validated seam-adjacent point (r,z), computed by
// whichever base construction was used for this ftype/domain mode (the
// half-domain chord point or the full-domain/extended ordered-map trace),
// by the bounded delta-s displacement. Contour-agnostic: the displacement
// primitive (tok_displace_xpt_seam_on_flux) only walks along the flux
// contour already passing through the given (r,z), so it composes safely
// regardless of which point-construction algorithm produced it -- this is
// what lets the same relaxation be reused for both half_domain=true and
// half_domain=false (extended) blocks instead of being wired only into the
// half-domain-only path it originated in.
static bool
tok_relax_xpt_seam_point(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx, double psi, double u, double *r, double *z)
{
  if (!inp->relaxed_xpt_seam)
    return true;

  double q = 0.0, delta_s = 0.0;
  if (!tok_xpt_seam_delta_s(inp, arc_ctx, psi, &q, &delta_s)) {
    fprintf(stderr,
      "TOK_XPT_SEAM_DIAG mode=rejected ftype=%d psi=%.17g q=%.17g coefficient=%.17g bound=%.17g reason=invalid_delta_s action=%s\n",
      inp->ftype, psi, q, inp->relaxed_xpt_seam_delta_s_coeff,
      inp->relaxed_xpt_seam_delta_s_bound,
      inp->relaxed_xpt_seam_sweep ? "reject_candidate" : "retain_straight");
    if (tok_xpt_seam_optimizer_trial(inp)) {
      tok_xpt_seam_trial_reject(inp,
        GKYL_XPT_SEAM_TRIAL_INVALID_PARAMETER);
      return true;
    }
    return !inp->relaxed_xpt_seam_sweep;
  }

  // Use signed local contour arc length, with s0=0 at the straight-ray
  // intersection. Return before doing coordinate arithmetic so zero mode is
  // bitwise identical to the existing construction.
  const double s0 = 0.0;
  if (delta_s == 0.0) {
    if (tok_ordered_map_diag_enabled() &&
        tok_xpt_seam_endpoint(inp->ftype, u))
      fprintf(stderr,
        "TOK_XPT_SEAM_DIAG mode=zero ftype=%d psi=%.17g q=%.17g s0=%.17g delta_s=%.17g seam_s=%.17g R=%.17g Z=%.17g construction=straight\n",
        inp->ftype, psi, q, s0, delta_s, s0, *r, *z);
    return true;
  }

  if (!inp->relaxed_xpt_seam_sweep || !inp->straight_xpt_ray) {
    fprintf(stderr,
      "TOK_XPT_SEAM_DIAG mode=rejected ftype=%d psi=%.17g q=%.17g s0=%.17g delta_s=%.17g reason=diagnostic_sweep_not_enabled fallback=straight\n",
      inp->ftype, psi, q, s0, delta_s);
    return true;
  }

  // The only free data are on the seam. A fixed linear arc weight carries
  // that boundary condition to the unchanged opposite edge, giving a smooth
  // diagnostic shadow map without introducing interior degrees of freedom.
  // A self-periodic single-null CORE block has no such opposite edge -- u=0
  // and u=1 are the same physical (X-point-adjacent) location, glued
  // together by the closed trace construction -- so a linear taper would
  // apply delta_s at one end and not the other, breaking the periodic
  // identification it is supposed to preserve. Apply it uniformly instead:
  // every point at this psi gets the same q-dependent radial bump
  // regardless of u, so u=0 and u=1 (which sample the same underlying raw
  // trace point) remain identically displaced and the wrap stays exact.
  double seam_weight = inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE
    ? 1.0 : tok_sep_fixed_edge_is_first(inp->ftype) ? u : 1.0-u;
  double point_delta_s = seam_weight*delta_s;
  if (point_delta_s == 0.0)
    return true;
  double candidate_r = *r, candidate_z = *z;
  double realized = 0.0;
  if (!tok_displace_xpt_seam_on_flux(arc_ctx->geo, psi, point_delta_s,
      &candidate_r, &candidate_z, &realized)) {
    fprintf(stderr,
      "TOK_XPT_SEAM_DIAG mode=rejected ftype=%d psi=%.17g q=%.17g u=%.17g delta_s=%.17g point_delta_s=%.17g reason=contour_tracking_failed action=reject_candidate fallback=none\n",
      inp->ftype, psi, q, u, delta_s, point_delta_s);
    if (tok_xpt_seam_optimizer_trial(inp)) {
      tok_xpt_seam_trial_reject(inp, GKYL_XPT_SEAM_TRIAL_CONTOUR);
      return true;
    }
    return false;
  }
  if (tok_xpt_seam_optimizer_trial(inp)) {
    struct gkyl_tok_geo_xpt_seam_trial_status *status =
      inp->relaxed_xpt_seam_trial_status;
    status->max_realized_displacement = fmax(
      status->max_realized_displacement, realized);
    double bound_tolerance = 64.0*DBL_EPSILON*fmax(1.0,
      inp->relaxed_xpt_seam_delta_s_bound);
    if (realized > inp->relaxed_xpt_seam_delta_s_bound+bound_tolerance) {
      tok_xpt_seam_trial_reject(inp, GKYL_XPT_SEAM_TRIAL_CONTOUR);
      return true;
    }
  }
  if (tok_xpt_seam_endpoint(inp->ftype, u) &&
      arc_ctx->xpt_ray_branch_valid) {
    bool resolved = false, on_right = false;
    if (!tok_xpt_classify_branch_at_point(inp, arc_ctx, candidate_r,
        candidate_z, psi, &resolved, &on_right) || !resolved ||
        on_right != arc_ctx->xpt_ray_on_right) {
      if (tok_xpt_seam_optimizer_trial(inp)) {
        tok_xpt_seam_trial_reject(inp, GKYL_XPT_SEAM_TRIAL_BRANCH);
        return true;
      }
      fprintf(stderr,
        "TOK_XPT_SEAM_DIAG mode=rejected ftype=%d psi=%.17g q=%.17g delta_s=%.17g reason=branch_identity_changed action=reject_candidate fallback=none\n",
        inp->ftype, psi, q, delta_s);
      return false;
    }
  }
  *r = candidate_r; *z = candidate_z;
  if (tok_ordered_map_diag_enabled() &&
      tok_xpt_seam_endpoint(inp->ftype, u))
    fprintf(stderr,
      "TOK_XPT_SEAM_DIAG mode=candidate ftype=%d psi=%.17g q=%.17g s0=%.17g delta_s=%.17g seam_s=%.17g R=%.17g Z=%.17g extension=fixed_edge_linear_arc_blend\n",
      inp->ftype, psi, q, s0, delta_s, s0+delta_s, *r, *z);
  return true;
}

// Half-domain (non-extended) seam point: the validated straight-ray chord
// point, relaxed by delta-s. Semantically unchanged from the pre-refactor
// implementation -- previously this function's own body started here, gated
// by "!inp->half_domain" in addition to "!inp->relaxed_xpt_seam"; since this
// call site is only ever reached when half_domain=true or straight_xpt_ray
// is false (see tok_build_current_ordered_trace/tok_ordered_map_lookup),
// that half_domain check was always true here and is preserved by construction
// rather than by an explicit condition.
static bool
tok_parameterized_xpt_seam_point(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx, double psi, double u, double *r, double *z)
{
  if (!tok_ordered_chord_point(inp, arc_ctx, psi, u, r, z))
    return false;
  return tok_relax_xpt_seam_point(inp, arc_ctx, psi, u, r, z);
}

static double
tok_fpol_at_psi(const struct gkyl_tok_geo *geo, double psi)
{
  double p = psi;
  if (p < geo->fgrid.lower[0] || p > geo->fgrid.upper[0])
    p = geo->sibry;
  int idx = GKYL_MIN2(geo->frange.upper[0], GKYL_MAX2(geo->frange.lower[0],
    geo->frange.lower[0]
      +(int) floor((p-geo->fgrid.lower[0])/geo->fgrid.dx[0])));
  long loc = gkyl_range_idx(&geo->frange, &idx);
  const double *coeffs = gkyl_array_cfetch(geo->fpoldg, loc);
  double xc;
  gkyl_rect_grid_cell_center(&geo->fgrid, &idx, &xc);
  double x = (p-xc)/(0.5*geo->fgrid.dx[0]);
  return geo->fbasis.eval_expand(&x, coeffs);
}

static bool
tok_ext_set_phi_reference(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx)
{
  struct tok_ext_topology top;
  if (!tok_ext_topology_from_ftype(inp->ftype, inp->half_domain, &top) ||
      arc_ctx->map_trace_n < 2)
    return false;
  if (top.phi_reference == TOK_EXT_PHI_LOWER) {
    arc_ctx->map_trace_phi_ref = arc_ctx->map_trace_phi[0];
    return true;
  }
  if (top.phi_reference == TOK_EXT_PHI_UPPER) {
    arc_ctx->map_trace_phi_ref =
      arc_ctx->map_trace_phi[arc_ctx->map_trace_n-1];
    return true;
  }

  double R[16] = { 0.0 }, dRdZ[16] = { 0.0 };
  double dR[16] = { 0.0 }, dZ[16] = { 0.0 };
  int nr = gkyl_tok_geo_R_psiZ(arc_ctx->geo, arc_ctx->psi,
    arc_ctx->geo->zmaxis, 16, R, dRdZ, dR, dZ);
  if (nr <= 0)
    return false;
  bool outboard = top.phi_reference == TOK_EXT_PHI_OUTBOARD_MIDPLANE;
  double target_r = tok_nearest_value(outboard ? inp->rright : inp->rleft,
    R, nr);
  double target_z = arc_ctx->geo->zmaxis;
  double best_d2 = DBL_MAX, best_phi = 0.0, max_step = 0.0;
  for (int i=0; i<arc_ctx->map_trace_n-1; ++i) {
    double r0 = arc_ctx->map_trace_r[i];
    double z0 = arc_ctx->map_trace_z[i];
    double dr = arc_ctx->map_trace_r[i+1]-r0;
    double dz = arc_ctx->map_trace_z[i+1]-z0;
    double den = dr*dr+dz*dz;
    if (!(den > 0.0))
      continue;
    max_step = fmax(max_step, sqrt(den));
    double w = ((target_r-r0)*dr+(target_z-z0)*dz)/den;
    w = fmin(1.0, fmax(0.0, w));
    double rp = r0+w*dr, zp = z0+w*dz;
    double d2 = SQ(target_r-rp)+SQ(target_z-zp);
    if (d2 < best_d2) {
      best_d2 = d2;
      best_phi = arc_ctx->map_trace_phi[i]
        +w*(arc_ctx->map_trace_phi[i+1]-arc_ctx->map_trace_phi[i]);
    }
  }
  if (!isfinite(best_phi) || !isfinite(best_d2) || !(max_step > 0.0) ||
      sqrt(best_d2) > 2.0*max_step) {
    fprintf(stderr,
      "TOK_ORDERED_MAP failed phi reference ftype=%d psi=%.17g side=%s distance=%.17g max_step=%.17g\n",
      inp->ftype, arc_ctx->psi, outboard ? "outboard" : "inboard",
      sqrt(best_d2), max_step);
    return false;
  }
  arc_ctx->map_trace_phi_ref = best_phi;
  return true;
}

// Project a point onto psi=const along grad psi -- the level set's own normal
// -- rather than along a coordinate axis.  Used where the axis-aligned solve
// has no solution at all: on the trace segment that leaves an X point the
// contour turns a corner, so the chord between two on-contour trace points
// bows onto the concave side and the line through it misses both branches.
// The displacement is bounded because |grad psi| collapses at the X point and
// an unclamped Newton step can land on a different separatrix branch.
static bool
tok_psi_normal_project(const struct gkyl_tok_geo *geo, double psi,
  double r0, double z0, double max_disp, double tol, double *r, double *z)
{
  double cr = r0, cz = z0;
  double step_cap = 0.25*max_disp;
  for (int k=0; k<32; ++k) {
    double f = tok_eval_psi_rz_local(geo, cr, cz)-psi;
    if (!isfinite(f))
      return false;
    if (fabs(f) <= tol) { *r = cr; *z = cz; return true; }
    double gr = 0.0, gz = 0.0;
    if (!tok_eval_psi_grad_rz_local(geo, cr, cz, &gr, &gz))
      return false;
    double g2 = gr*gr+gz*gz;
    if (!(g2 > 0.0))
      return false;
    double dr = -f/g2*gr, dz = -f/g2*gz;
    double dl = hypot(dr, dz);
    if (dl > step_cap) { dr *= step_cap/dl; dz *= step_cap/dl; }
    cr += dr; cz += dz;
    if (!isfinite(cr) || !isfinite(cz))
      return false;
    if (hypot(cr-r0, cz-z0) > max_disp)
      return false;
  }
  double f = tok_eval_psi_rz_local(geo, cr, cz)-psi;
  if (isfinite(f) && fabs(f) <= tol) { *r = cr; *z = cz; return true; }
  return false;
}

// Sample a trace whose array index, rather than physical arc length, is
// uniform in the logical block coordinate.  Project the local interpolation
// back to psi=constant using whichever physical coordinate varies most on the
// segment; this remains well conditioned at R and Z turning points.
static bool
tok_logical_trace_sample(const struct gkyl_tok_geo *geo, double psi,
  const double *tr, const double *tz, int n, double u,
  double *r, double *z)
{
  if (n < 2)
    return false;
  const double endpoint_tol = 256.0*DBL_EPSILON;
  if (u <= endpoint_tol) { *r = tr[0]; *z = tz[0]; return true; }
  if (u >= 1.0-endpoint_tol) {
    *r = tr[n-1]; *z = tz[n-1]; return true;
  }
  double x = u*(n-1);
  int i = GKYL_MIN2(n-2, GKYL_MAX2(0, (int) floor(x)));
  double w = x-i;
  double rlin = tr[i]+w*(tr[i+1]-tr[i]);
  double zlin = tz[i]+w*(tz[i+1]-tz[i]);
  double dr = tr[i+1]-tr[i], dz = tz[i+1]-tz[i];
  const double tol = 1e-9*fmax(1.0, fabs(psi));
  bool axis_ok = true;
  if (fabs(dr) >= fabs(dz)) {
    double roots[32] = { 0.0 };
    int nr = tok_geo_Z_psiR(geo, psi, rlin, 32, roots);
    if (nr <= 0)
      axis_ok = false;
    else { *r = rlin; *z = tok_nearest_value(zlin, roots, nr); }
  }
  else {
    double roots[16] = { 0.0 }, dRdZ[16] = { 0.0 };
    double dR[16] = { 0.0 }, dZ[16] = { 0.0 };
    int nr = gkyl_tok_geo_R_psiZ(geo, psi, zlin, 16,
      roots, dRdZ, dR, dZ);
    if (nr <= 0)
      axis_ok = false;
    else {
      double r_choice = tok_nearest_value(rlin, roots, nr);
      if (tok_rlin_ambiguous(rlin, roots, nr)) {
        double r_walk = 0.0;
        if (tok_trace_walk_root(geo, psi, tz[i], tr[i], tz[i+1], tr[i+1],
            zlin, &r_walk))
          r_choice = r_walk;
      }
      *r = r_choice; *z = zlin;
    }
  }
  double seglen = hypot(dr, dz);
  if (axis_ok) {
    // The axis-aligned solve takes the root nearest the chord point, but when
    // no root is near it returns a distant one instead of failing: on an
    // inboard SOL surface the same R column meets the contour again far down
    // the divertor leg, and that root silently replaces the intended one.
    // Both bracketing trace points are on the contour, so the answer cannot
    // be further from the chord than the segment itself.
    double away = hypot(*r-rlin, *z-zlin);
    if (away > 0.5*seglen) {
      if (tok_ordered_map_diag_enabled())
        fprintf(stderr,
          "TOK_LOGSAMP_NONLOCAL u=%.17g i=%d chord=(%.17g,%.17g) "
          "root=(%.17g,%.17g) away=%.17g seglen=%.17g psi=%.17g\n",
          u, i, rlin, zlin, *r, *z, away, seglen, psi);
      axis_ok = false;
    }
  }
  if (axis_ok) {
    double residual = tok_eval_psi_rz_local(geo, *r, *z)-psi;
    if (isfinite(*r) && isfinite(*z) && isfinite(residual) &&
        fabs(residual) <= tol)
      return true;
  }
  // The axis-aligned solve has no answer on the segment that leaves an X
  // point: the contour turns a corner there, so the chord between two
  // on-contour trace points bows onto the concave side and neither the
  // vertical nor the horizontal line through it meets psi=const.  Both trace
  // points are exact, so project the chord point back along the level set's
  // own normal instead, bounded by half the segment so the result stays on
  // this segment's branch.
  double pr = 0.0, pz = 0.0;
  if (seglen > 0.0 &&
      tok_psi_normal_project(geo, psi, rlin, zlin, 0.5*seglen, tol, &pr, &pz)) {
    if (tok_ordered_map_diag_enabled())
      fprintf(stderr,
        "TOK_LOGSAMP_NORMAL_PROJECT u=%.17g i=%d axis_ok=%d "
        "chord=(%.17g,%.17g) out=(%.17g,%.17g) disp=%.17g seglen=%.17g "
        "psi=%.17g\n",
        u, i, (int) axis_ok, rlin, zlin, pr, pz,
        hypot(pr-rlin, pz-zlin), seglen, psi);
    *r = pr; *z = pz;
    return true;
  }
  // Last resort, and the only one that survives a DG cell face.  The quadratic
  // psi representation is C0: psi is continuous across a face but grad psi is
  // not, so a Newton projection started on a face -- which is exactly where
  // 204980's segment sits, straddling Z = -7*dZ = -0.9625 to 13 digits -- has
  // no well-defined direction to move in and the step above fails.  Bisecting
  // for a sign change of psi-target along the segment normal never evaluates a
  // gradient, so the kink cannot defeat it.  Both bracketing trace points are
  // exact contour points, so a crossing exists within a segment length; bound
  // the search there to stay on this segment's branch.
  if (seglen > 0.0) {
    double nx = -dz/seglen, ny = dr/seglen;
    const int nstep = 32;
    double h = seglen/nstep;
    double f0 = tok_eval_psi_rz_local(geo, rlin, zlin)-psi;
    if (isfinite(f0)) {
      double fprev_p = f0, fprev_m = f0;
      for (int k=1; k<=nstep; ++k) {
        for (int side=0; side<2; ++side) {
          double sgn = side == 0 ? 1.0 : -1.0;
          double t = sgn*k*h;
          double fk = tok_eval_psi_rz_local(geo,
            rlin+t*nx, zlin+t*ny)-psi;
          double fprev = side == 0 ? fprev_p : fprev_m;
          if (isfinite(fk) && isfinite(fprev) && (fk < 0.0) != (fprev < 0.0)) {
            double lo = sgn*(k-1)*h, hi = t, flo = fprev;
            for (int m=0; m<80; ++m) {
              double mid = 0.5*(lo+hi);
              double fm = tok_eval_psi_rz_local(geo,
                rlin+mid*nx, zlin+mid*ny)-psi;
              if (!isfinite(fm)) break;
              if ((fm < 0.0) == (flo < 0.0)) { lo = mid; flo = fm; } else hi = mid;
            }
            double tb = 0.5*(lo+hi);
            *r = rlin+tb*nx; *z = zlin+tb*ny;
            if (tok_ordered_map_diag_enabled())
              fprintf(stderr,
                "TOK_LOGSAMP_PERP_BISECT u=%.17g i=%d chord=(%.17g,%.17g) "
                "out=(%.17g,%.17g) disp=%.17g seglen=%.17g psi=%.17g\n",
                u, i, rlin, zlin, *r, *z, fabs(tb), seglen, psi);
            return true;
          }
          if (side == 0) fprev_p = fk; else fprev_m = fk;
        }
      }
    }
  }
  // Fourth resort: the requested level set does not exist anywhere near this
  // segment, so there is no crossing for any of the three searches above to
  // find and widening their windows cannot help.  psisep is psi evaluated at
  // whatever point the X-point finder returned, and when that point is not a
  // critical point of *this* (quadratic) representation the psi=psisep level
  // set is severed near the X point rather than crossing itself there.
  // Measured on 203997, whose quadratic rep has no in-vessel critical point at
  // all so the cubic location is substituted (|grad psi| there is 5.8e-4, not
  // ~0): {psi >= psisep} splits into two components separated by a 10.6 mm gap
  // straddling the R = 0.67875 cell face, where |grad psi| jumps 9x across the
  // C0 seam.  The trace's last segment spans that gap, so the level set is
  // absent over its whole length and psi stays below target across +/-20 mm of
  // normal offset.
  //
  // The two bracketing trace points remain the best available representatives
  // of the contour, so fall back to the linear interpolation between them --
  // which is exactly what the trace polyline already is between its stations.
  // This only runs where the routine previously returned false and aborted the
  // block, so it cannot move a node that any other path could place; and the
  // projection it replaces is a sub-micron refinement of the chord point
  // wherever the level set does exist.  The resulting grid is still subject to
  // the fold and surface-crossing checks, which is what decides whether the
  // interpolated node is acceptable.
  if (!tok_logsamp_chord_fallback_disabled()) {
    double residual = tok_eval_psi_rz_local(geo, rlin, zlin)-psi;
    if (isfinite(rlin) && isfinite(zlin) && isfinite(residual)) {
      fprintf(stderr,
        "TOK_LOGSAMP_CHORD_FALLBACK u=%.17g i=%d chord=(%.17g,%.17g) "
        "residual=%.17g seglen=%.17g psi=%.17g\n",
        u, i, rlin, zlin, residual, seglen, psi);
      *r = rlin; *z = zlin;
      return true;
    }
  }
  fprintf(stderr,
    "TOK_LOGSAMP reason=no_projection axis_ok=%d u=%.17g i=%d w=%.17g "
    "chord=(%.17g,%.17g) p0=(%.17g,%.17g) p1=(%.17g,%.17g) seglen=%.17g "
    "psi=%.17g\n",
    (int) axis_ok, u, i, w, rlin, zlin, tr[i], tz[i], tr[i+1], tz[i+1],
    seglen, psi);
  return false;
}

static bool
tok_build_current_ordered_trace(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx)
{
  double psi = arc_ctx->psi;
  if (arc_ctx->map_trace_initialized && tok_geo_same_flux(
      psi, arc_ctx->map_trace_psi))
    return true;
  if (!tok_build_trace_correspondence(inp, arc_ctx)) {
    fprintf(stderr,
      "TOK_ORDERED_MAP_TRACE reason=correspondence_failed ftype=%d psi=%.17g\n",
      inp->ftype, psi);
    return false;
  }
  int n = tok_ext_map_trace_nodes(inp, arc_ctx->sep_trace_capacity);
  const bool extended = tok_ext_construction(inp);
  double *raw_r = 0, *raw_z = 0, *raw_s = 0;
  int raw_n = 0;
  bool raw_param_is_r = false, raw_closed = false;
  double radial_fraction = 0.0;
  if (extended) {
    int raw_capacity = arc_ctx->sep_trace_capacity;
    raw_r = gkyl_malloc(sizeof(double[raw_capacity]));
    raw_z = gkyl_malloc(sizeof(double[raw_capacity]));
    raw_s = gkyl_malloc(sizeof(double[raw_capacity]));
    bool at_sep = tok_geo_same_flux(psi, arc_ctx->geo->psisep);
    if (!tok_ext_build_domain_trace(inp, arc_ctx, psi, at_sep,
        raw_r, raw_z, raw_s, &raw_n, &raw_param_is_r, &raw_closed)) {
      fprintf(stderr,
        "TOK_ORDERED_MAP_TRACE reason=domain_trace_failed ftype=%d psi=%.17g at_sep=%d\n",
        inp->ftype, psi, (int) at_sep);
      gkyl_free(raw_r); gkyl_free(raw_z); gkyl_free(raw_s);
      return false;
    }
    struct tok_ext_topology top;
    bool have_top = tok_ext_topology_from_ftype(inp->ftype, inp->half_domain, &top);
    if (!have_top || raw_closed != top.closed) {
      fprintf(stderr,
        "TOK_ORDERED_MAP_TRACE reason=topology_mismatch ftype=%d psi=%.17g have_topology=%d raw_closed=%d expected_closed=%d\n",
        inp->ftype, psi, (int) have_top, (int) raw_closed,
        have_top ? (int) top.closed : -1);
      gkyl_free(raw_r); gkyl_free(raw_z); gkyl_free(raw_s);
      return false;
    }
    double span = arc_ctx->xpt_ray_psi0-arc_ctx->geo->psisep;
    if (!isfinite(span) || fabs(span) <= 256.0*DBL_EPSILON*
        fmax(1.0, fabs(arc_ctx->geo->psisep))) {
      fprintf(stderr,
        "TOK_ORDERED_MAP_TRACE reason=degenerate_radial_span ftype=%d psi=%.17g span=%.17g ray_psi0=%.17g psisep=%.17g\n",
        inp->ftype, psi, span, arc_ctx->xpt_ray_psi0, arc_ctx->geo->psisep);
      gkyl_free(raw_r); gkyl_free(raw_z); gkyl_free(raw_s);
      return false;
    }
    radial_fraction = (psi-arc_ctx->geo->psisep)/span;
    if (radial_fraction < -1e-8 || radial_fraction > 1.0+1e-8) {
      fprintf(stderr,
        "TOK_ORDERED_MAP_TRACE reason=radial_fraction_out_of_range ftype=%d psi=%.17g radial_fraction=%.17g span=%.17g psisep=%.17g\n",
        inp->ftype, psi, radial_fraction, span, arc_ctx->geo->psisep);
      gkyl_free(raw_r); gkyl_free(raw_z); gkyl_free(raw_s);
      return false;
    }
    radial_fraction = fmin(1.0, fmax(0.0, radial_fraction));
  }
  double fpol = tok_fpol_at_psi(arc_ctx->geo, psi);
  arc_ctx->map_trace_s[0] = 0.0;
  arc_ctx->map_trace_phi[0] = 0.0;
  for (int i=0; i<n; ++i) {
    double u = i/(double) (n-1);
    bool point_ok = false;
    if (extended) {
      double w;
      if (arc_ctx->ext_ladder_initialized && n == arc_ctx->ext_ladder_n) {
        w = tok_ext_ladder_w(arc_ctx, radial_fraction, i);
      }
      else {
        double v = tok_trace_correspondence(arc_ctx, u);
        w = (1.0-radial_fraction)*u+radial_fraction*v;
      }
      point_ok = tok_trace_sample(arc_ctx->geo, psi,
        raw_r, raw_z, raw_s, raw_n, raw_param_is_r, w,
        &arc_ctx->map_trace_r[i], &arc_ctx->map_trace_z[i]);
      // Relaxed delta-s was originally wired only into the half-domain chord
      // point below. Apply the same bounded relaxation here so full-domain
      // (production) geometry can use it too; a no-op unless relaxed_xpt_seam
      // is set (see tok_relax_xpt_seam_point).
      if (point_ok)
        point_ok = tok_relax_xpt_seam_point(inp, arc_ctx, psi, u,
          &arc_ctx->map_trace_r[i], &arc_ctx->map_trace_z[i]);
    }
    else {
      point_ok = tok_parameterized_xpt_seam_point(inp, arc_ctx, psi, u,
        &arc_ctx->map_trace_r[i], &arc_ctx->map_trace_z[i]);
    }
    if (!point_ok) {
      fprintf(stderr,
        "TOK_ORDERED_MAP_TRACE reason=trace_point_failed ftype=%d psi=%.17g i=%d n=%d u=%.17g extended=%d radial_fraction=%.17g\n",
        inp->ftype, psi, i, n, u, (int) extended, radial_fraction);
      if (extended) {
        gkyl_free(raw_r); gkyl_free(raw_z); gkyl_free(raw_s);
      }
      return false;
    }
    if (i > 0) {
      double ds = hypot(arc_ctx->map_trace_r[i]-arc_ctx->map_trace_r[i-1],
        arc_ctx->map_trace_z[i]-arc_ctx->map_trace_z[i-1]);
      if (!(ds > 0.0) || !isfinite(ds)) {
        fprintf(stderr,
          "TOK_ORDERED_MAP_TRACE reason=nonpositive_arc_step ftype=%d psi=%.17g i=%d n=%d ds=%.17g prev=(%.17g,%.17g) curr=(%.17g,%.17g)\n",
          inp->ftype, psi, i, n, ds,
          arc_ctx->map_trace_r[i-1], arc_ctx->map_trace_z[i-1],
          arc_ctx->map_trace_r[i], arc_ctx->map_trace_z[i]);
        if (extended) {
          gkyl_free(raw_r); gkyl_free(raw_z); gkyl_free(raw_s);
        }
        return false;
      }
      arc_ctx->map_trace_s[i] = arc_ctx->map_trace_s[i-1]+ds;
      double rm = 0.5*(arc_ctx->map_trace_r[i]
        +arc_ctx->map_trace_r[i-1]);
      double zm = 0.5*(arc_ctx->map_trace_z[i]
        +arc_ctx->map_trace_z[i-1]);
      double gr = 0.0, gz = 0.0;
      if (!tok_eval_psi_grad_rz_local(arc_ctx->geo, rm, zm, &gr, &gz)) {
        fprintf(stderr,
          "TOK_ORDERED_MAP_TRACE reason=psi_grad_eval_failed ftype=%d psi=%.17g i=%d mid=(%.17g,%.17g)\n",
          inp->ftype, psi, i, rm, zm);
        if (extended) {
          gkyl_free(raw_r); gkyl_free(raw_z); gkyl_free(raw_s);
        }
        return false;
      }
      double grad = hypot(gr, gz);
      if (!(grad > 1e-14) || !(rm > 0.0)) {
        fprintf(stderr,
          "TOK_ORDERED_MAP_TRACE reason=degenerate_grad_or_radius ftype=%d psi=%.17g i=%d grad=%.17g mid=(%.17g,%.17g)\n",
          inp->ftype, psi, i, grad, rm, zm);
        if (extended) {
          gkyl_free(raw_r); gkyl_free(raw_z); gkyl_free(raw_s);
        }
        return false;
      }
      arc_ctx->map_trace_phi[i] = arc_ctx->map_trace_phi[i-1]
        +fpol*ds/(rm*grad);
    }
  }
  if (extended) {
    gkyl_free(raw_r); gkyl_free(raw_z); gkyl_free(raw_s);
  }
  arc_ctx->map_trace_n = n;
  arc_ctx->map_trace_psi = psi;
  if (tok_ext_construction(inp) &&
      !tok_ext_set_phi_reference(inp, arc_ctx)) {
    fprintf(stderr,
      "TOK_ORDERED_MAP_TRACE reason=phi_reference_failed ftype=%d psi=%.17g\n",
      inp->ftype, psi);
    return false;
  }
  arc_ctx->map_trace_initialized = true;
  return true;
}

static bool
tok_ordered_map_lookup(const struct gkyl_tok_geo_grid_inp *inp,
  struct arc_length_ctx *arc_ctx, double theta, double alpha,
  struct tok_ordered_point *out)
{
  if (!tok_xpt_mapping_requested(inp))
    return false;
  const struct gkyl_tok_geo_grid_inp *effective_inp = inp;
  struct gkyl_tok_geo_grid_inp straight_inp;
  if (!tok_build_current_ordered_trace(effective_inp, arc_ctx)) {
    if (tok_xpt_seam_optimizer_trial(inp)) {
      tok_xpt_seam_trial_reject(inp, GKYL_XPT_SEAM_TRIAL_TRACE_ORDERING);
      straight_inp = *inp;
      straight_inp.relaxed_xpt_seam_delta_s_coeff = 0.0;
      straight_inp.relaxed_xpt_seam_sweep = false;
      straight_inp.relaxed_xpt_seam_optimizer_trial = false;
      straight_inp.relaxed_xpt_seam_trial_status = 0;
      arc_ctx->map_trace_initialized = false;
      if (!tok_build_current_ordered_trace(&straight_inp, arc_ctx)) {
        fprintf(stderr,
          "TOK_ORDERED_MAP straight trial fallback failed ftype=%d psi=%.17g\n",
          inp->ftype, arc_ctx->psi);
        abort();
      }
      effective_inp = &straight_inp;
    }
    else {
      fprintf(stderr,
        "TOK_ORDERED_MAP initialization failed ftype=%d psi=%.17g\n",
        inp->ftype, arc_ctx->psi);
      abort();
    }
  }
  double dtheta = inp->cgrid.upper[2]-inp->cgrid.lower[2];
  double u = (theta-inp->cgrid.lower[2])/dtheta;
  static int dump_ftype = -2;
  if (dump_ftype == -2) {
    const char *e = getenv("GKYL_TOK_LOOKUP_DUMP_FTYPE");
    dump_ftype = e ? atoi(e) : -1;
  }
  if (u < -1e-10 || u > 1.0+1e-10) {
    fprintf(stderr,
      "TOK_ORDERED_MAP_LOOKUP reason=u_out_of_range ftype=%d psi=%.17g theta=%.17g u=%.17g "
      "cgrid_ndim=%d cgrid_lo=(%.17g,%.17g,%.17g) cgrid_up=(%.17g,%.17g,%.17g) "
      "cgrid_cells=(%d,%d,%d) dtheta=%.17g\n",
      inp->ftype, arc_ctx->psi, theta, u,
      inp->cgrid.ndim,
      inp->cgrid.lower[0], inp->cgrid.lower[1], inp->cgrid.lower[2],
      inp->cgrid.upper[0], inp->cgrid.upper[1], inp->cgrid.upper[2],
      inp->cgrid.cells[0], inp->cgrid.cells[1], inp->cgrid.cells[2],
      dtheta);
    return false;
  }
  if (u <= 256.0*DBL_EPSILON) u = 0.0;
  else if (u >= 1.0-256.0*DBL_EPSILON) u = 1.0;
  else u = fmin(1.0, fmax(0.0, u));
  if (tok_ext_construction(effective_inp)) {
    if (!tok_logical_trace_sample(arc_ctx->geo, arc_ctx->psi,
        arc_ctx->map_trace_r, arc_ctx->map_trace_z,
        arc_ctx->map_trace_n, u, &out->r, &out->z)) {
      fprintf(stderr,
        "TOK_ORDERED_MAP_LOOKUP reason=logical_trace_sample_failed ftype=%d psi=%.17g u=%.17g n=%d\n",
        inp->ftype, arc_ctx->psi, u, arc_ctx->map_trace_n);
      return false;
    }
  }
  else if (!tok_parameterized_xpt_seam_point(effective_inp, arc_ctx,
      arc_ctx->psi, u,
      &out->r, &out->z)) {
    fprintf(stderr,
      "TOK_ORDERED_MAP_LOOKUP reason=seam_point_failed ftype=%d psi=%.17g u=%.17g\n",
      inp->ftype, arc_ctx->psi, u);
    return false;
  }
  if (dump_ftype == inp->ftype) {
    double mzlo = arc_ctx->map_trace_z[0], mzhi = mzlo;
    for (int k=1; k<arc_ctx->map_trace_n; ++k) {
      mzlo = fmin(mzlo, arc_ctx->map_trace_z[k]);
      mzhi = fmax(mzhi, arc_ctx->map_trace_z[k]);
    }
    double xx = u*(arc_ctx->map_trace_n-1);
    int ii = GKYL_MIN2(arc_ctx->map_trace_n-2, GKYL_MAX2(0, (int) floor(xx)));
    fprintf(stderr,
      "TOK_LOOKUP_DUMP ftype=%d psi=%.17g theta=%.17g u=%.17g out=(%.17g,%.17g) "
      "ext=%d mtrace_zbox=[%.17g,%.17g] n=%d i=%d p0=(%.17g,%.17g) "
      "p1=(%.17g,%.17g)\n",
      inp->ftype, arc_ctx->psi, theta, u, out->r, out->z,
      (int) tok_ext_construction(effective_inp), mzlo, mzhi,
      arc_ctx->map_trace_n, ii,
      arc_ctx->map_trace_r[ii], arc_ctx->map_trace_z[ii],
      arc_ctx->map_trace_r[ii+1], arc_ctx->map_trace_z[ii+1]);
  }
  double x = u*(arc_ctx->map_trace_n-1);
  int i = GKYL_MIN2(arc_ctx->map_trace_n-2, GKYL_MAX2(0, (int) floor(x)));
  double du = 1.0/(arc_ctx->map_trace_n-1);
  double dr = arc_ctx->map_trace_r[i+1]-arc_ctx->map_trace_r[i];
  double dz = arc_ctx->map_trace_z[i+1]-arc_ctx->map_trace_z[i];
  double speed_u = hypot(dr, dz)/du;
  double gr = 0.0, gz = 0.0;
  if (!tok_eval_psi_grad_rz_local(arc_ctx->geo, out->r, out->z,
      &gr, &gz)) {
    fprintf(stderr,
      "TOK_ORDERED_MAP_LOOKUP reason=grad_eval_failed ftype=%d psi=%.17g rz=(%.17g,%.17g)\n",
      inp->ftype, arc_ctx->psi, out->r, out->z);
    return false;
  }
  double tr = dr, tz = dz, grad = hypot(gr, gz);
  if (tok_ordered_map_diag_enabled() &&
      tok_ext_construction(inp) && u == 0.0) {
    fprintf(stderr,
      "TOK_ORDERED_MAP_DIAG endpoint_tangent ftype=%d psi=%.17g R=%.17g Z=%.17g dRdu=%.17g dZdu=%.17g gradR=%.17g gradZ=%.17g canonical_dot=%.17g\n",
      inp->ftype, arc_ctx->psi, out->r, out->z, dr/du, dz/du,
      gr, gz, (-gz)*dr+gr*dz);
  }
  if (grad > 1e-14) {
    tr = -gz/grad; tz = gr/grad;
    if (tr*dr+tz*dz < 0.0) { tr = -tr; tz = -tz; }
  }
  else {
    double tmag = hypot(tr, tz);
    if (!(tmag > 0.0)) {
      fprintf(stderr,
        "TOK_ORDERED_MAP_LOOKUP reason=degenerate_tangent ftype=%d psi=%.17g u=%.17g grad=%.17g\n",
        inp->ftype, arc_ctx->psi, u, grad);
      return false;
    }
    tr /= tmag; tz /= tmag;
  }
  out->dr_dtheta = tr*speed_u/dtheta;
  out->dz_dtheta = tz*speed_u/dtheta;
  double w = x-i;
  double path_phi = arc_ctx->map_trace_phi[i]
    +w*(arc_ctx->map_trace_phi[i+1]-arc_ctx->map_trace_phi[i]);
  double ref_phi = tok_ext_construction(effective_inp)
    ? arc_ctx->map_trace_phi_ref
    : (tok_sep_fixed_edge_is_first(effective_inp->ftype) ? 0.0
      : arc_ctx->map_trace_phi[arc_ctx->map_trace_n-1]);
  out->phi = alpha+path_phi-ref_phi;
  out->dphi_dtheta = (arc_ctx->map_trace_phi[i+1]
    -arc_ctx->map_trace_phi[i])/du/dtheta;
  return isfinite(out->r) && isfinite(out->z) && isfinite(out->phi) &&
    isfinite(out->dr_dtheta) && isfinite(out->dz_dtheta) &&
    isfinite(out->dphi_dtheta);
}

// Function to calculate phi given alpha
double
phi_func(double alpha_curr, double Z, void *ctx)
{
  struct arc_length_ctx *actx = ctx;
  double *arc_memo = actx->arc_memo;
  double psi = actx->psi, rclose = actx->rclose, zmin = actx->zmin, arcL = actx->arcL, zmax = actx->zmax;

  // Here we will abandon conventions about alpha and phi except for full core and full SN cases
  // The convention for phi only affects b_x - it does not affect any quantities used in axisymmetric simulations
  // I have not quite figured out full 3D yet. b_x presents a serious problem as of now. Akash Shukla 1/20/2024
  // The idea for axisymmetry is that I am avoiding starting integrals at the x-point to minimize issues
  double ival = 0;
  double phi_ref = 0.0;
  if (actx->ftype==GKYL_GEOMETRY_TOKAMAK_CORE){ // phi = alpha at outboard midplane
    if(actx->right==true){
      if(Z<actx->zmaxis)
        ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmaxis, rclose, false, false, arc_memo);
      else
        ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmaxis, Z, rclose, false, false, arc_memo);
    }
    else{
      ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmax, rclose, false, false, arc_memo);
      phi_ref = actx->phi_right;
    }
  }
  else if (actx->ftype==GKYL_GEOMETRY_TOKAMAK_CORE_L){ // alpha = phi at inboard midplane
    if(Z<actx->zmaxis)
      ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmaxis, rclose, false, false, arc_memo);
    else
      ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmaxis, Z, rclose, false, false, arc_memo);
  }

  else if (actx->ftype==GKYL_GEOMETRY_TOKAMAK_CORE_R){ // alpha = phi at outboard midplane
    if(Z<actx->zmaxis)
      ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmaxis, rclose, false, false, arc_memo);
    else
      ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmaxis, Z, rclose, false, false, arc_memo);
  }

  else if( (actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT) || (actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID)){ // alpha = phi at outboard midplane
    if (Z<actx->zmaxis)
      ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmaxis, rclose, false, false, arc_memo);
    else
      ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmaxis, Z, rclose, false, false, arc_memo);
  }
  else if(actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO){ // alpha = phi at lower plate and increases towards xpt
    ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, false, false, arc_memo);
  }
  else if(actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP){ //alpha = phi at upper plate and decreases towards xpt
    ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, false, false, arc_memo);
  }
  if( (actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN) || (actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID) ){ // alpha = phi at inboard midplane
    if (Z<actx->zmaxis)
      ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmaxis, rclose, false, false, arc_memo);
    else
      ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmaxis, Z, rclose, false, false, arc_memo);
  }
  else if(actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO){ // alpha = phi at lower plate and decreases towards xpt
    ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, false, false, arc_memo);
  }
  else if(actx->ftype==GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP){// alpha = phi at upper plate and increases towards xpt
    ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, false, false, arc_memo);
  }
  else if(actx->ftype==GKYL_GEOMETRY_TOKAMAK_LSN_SOL  || actx->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID){
    // alpha = phi at outboard midplane
    if (actx->right==true){
      if (Z<actx->zmaxis)
        ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmaxis, rclose, false, false, arc_memo);
      else
        ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmaxis, Z, rclose, false, false, arc_memo);
    }
    else{
      ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmax, rclose, false, false, arc_memo);
      phi_ref = actx->phi_right;
    }
  }
  else if(actx->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO){ //alpha = phi at outer plate and increases towards xpt
    ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, false, false, arc_memo);
  }
  else if(actx->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP){ //alpha = phi at inner plate and decreases towards xpt
    ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, false, false, arc_memo);
  }
  else if(actx->ftype==GKYL_GEOMETRY_TOKAMAK_PF_LO_R){ // alpha = phi at outer plate and increases towards xpt
      ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, false, false, arc_memo);
  }
  else if(actx->ftype==GKYL_GEOMETRY_TOKAMAK_PF_LO_L){ //alpha = phi at inner plate and decreases towards xpt
      ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, zmin, Z, rclose, false, false, arc_memo);
  }
  else if(actx->ftype==GKYL_GEOMETRY_TOKAMAK_PF_UP_R){ // alpha = phi at outer plate and decreases towards Xpt
      ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, false, false, arc_memo);
  }
  else if(actx->ftype==GKYL_GEOMETRY_TOKAMAK_PF_UP_L){ // alpha = phi at inner plate and increases towards xpt
      ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, zmax, rclose, false, false, arc_memo);
  }
  else if (actx->ftype==GKYL_GEOMETRY_TOKAMAK_IWL) {
    // phi = alpha at outboard midplane
    if(actx->right==true){
      if(Z<actx->zmaxis)
        ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmaxis, rclose, false, false, arc_memo);
      else
        ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmaxis, Z, rclose, false, false, arc_memo);
    }
    else{
      if (Z<actx->zmaxis) {
        ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmin, Z, rclose, false, false, arc_memo) ;
        phi_ref  = -actx->phi_right;
      }
      else {
        ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, Z, actx->zmax, rclose, false, false, arc_memo) ;
        phi_ref  = actx->phi_right;
      }
    }
  }
  // Now multiply by fpol
  double R[4] = {0};
  double dRdZ[4] = {0};
  double dR[4] = {0};
  double dZ[4] = {0};
  int nr = gkyl_tok_geo_R_psiZ(actx->geo, psi, Z, 4, R, dRdZ, dR, dZ);
  double r_curr = nr == 1 ? R[0] : choose_closest(rclose, R, R, nr);
  double psi_fpol = psi;
  if ( (psi_fpol < actx->geo->fgrid.lower[0]) || (psi_fpol > actx->geo->fgrid.upper[0]) ) // F = F(psi_sep) in the SOL.
    psi_fpol = actx->geo->sibry;
  int idx = fmin(actx->geo->frange.lower[0] + (int) floor((psi_fpol - actx->geo->fgrid.lower[0])/actx->geo->fgrid.dx[0]), actx->geo->frange.upper[0]);
  long loc = gkyl_range_idx(&actx->geo->frange, &idx);
  const double *coeffs = gkyl_array_cfetch(actx->geo->fpoldg,loc);
  double fxc;
  gkyl_rect_grid_cell_center(&actx->geo->fgrid, &idx, &fxc);
  double fx = (psi_fpol-fxc)/(actx->geo->fgrid.dx[0]*0.5);
  double fpol = actx->geo->fbasis.eval_expand(&fx, coeffs);
  ival = ival*fpol;

  //while(ival < -M_PI){
  //  ival +=2*M_PI;
  //}
  //while(ival > M_PI){
  //  ival -=2*M_PI;
  //}
  return alpha_curr + ival + phi_ref;
}

double
qprofile_func(void *ctx)
{
  // Function to calculate the flux surface averaged q profile.

  struct arc_length_ctx *actx = ctx;
  double *arc_memo = actx->arc_memo;
  double psi = actx->psi, rclose = actx->rclose, zmin = actx->zmin, arcL = actx->arcL, zmax = actx->zmax;
  double rleft =  actx->rleft, rright = actx->rright;

  // Calculate q(psi) = -F(psi)/2pi * integral_zmin^zmax 1/Rgrad(psi).

  double ival = 0;
  double phi_ref = 0.0;
  if (actx->ftype == GKYL_GEOMETRY_TOKAMAK_CORE) {
    double ival1 = integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmin, actx->zmax, rright, false, false, arc_memo);
    double ival2 = -integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmin, actx->zmax, rleft, false, false, arc_memo);
    ival = ival1 + ival2;
  }

  if (actx->ftype == GKYL_GEOMETRY_TOKAMAK_IWL) {
    double ival1 = integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmin, actx->zmax, rright, false, false, arc_memo);
    double ival2 = -integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmin, actx->zmax, rleft, false, false, arc_memo);
    ival = ival1 + ival2;
  }

  if (actx->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL) {
    double ival1 = integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmin_right, actx->zmax, rright, false, false, arc_memo);
    double ival2 = -integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmin_left, actx->zmax, rleft, false, false, arc_memo);
    ival = ival1 + ival2;
  }

  if (actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT || actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID || actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO || actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP){
    ival = integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmin, actx->zmax, rclose, false, false, arc_memo);
  }

  if (actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN || actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID || actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO || actx->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP) {
    ival = -integrate_phi_along_psi_contour_memo(actx->geo, psi, actx->zmin, actx->zmax, rclose, false, false, arc_memo);
  }

  // Now multiply by fpol/2pi.
  double R[4] = {0};
  double dR[4] = {0};
  double psi_fpol = psi;
  if ( (psi_fpol < actx->geo->fgrid.lower[0]) || (psi_fpol > actx->geo->fgrid.upper[0]) ) // F = F(psi_sep) in the SOL.
    psi_fpol = actx->geo->sibry;
  int idx = fmin(actx->geo->frange.lower[0] + (int) floor((psi_fpol - actx->geo->fgrid.lower[0])/actx->geo->fgrid.dx[0]), actx->geo->frange.upper[0]);
  long loc = gkyl_range_idx(&actx->geo->frange, &idx);
  const double *coeffs = gkyl_array_cfetch(actx->geo->fpoldg,loc);
  double fxc;
  gkyl_rect_grid_cell_center(&actx->geo->fgrid, &idx, &fxc);
  double fx = (psi_fpol-fxc)/(actx->geo->fgrid.dx[0]*0.5);
  double fpol = actx->geo->fbasis.eval_expand(&fx, coeffs);
  double qout = -ival*fpol/M_PI;

  // AS 1/15/25: The 3 lines below are a useful check to compare against q from efit.
  //coeffs = gkyl_array_cfetch(actx->geo->qdg,loc);
  //double q_efit = actx->geo->fbasis.eval_expand(&fx, coeffs);
  // printf("psi_curr = %g, my q = %g, efit q = %g\n", psi_fpol, qout, q_efit);

  return qout;
}

static double
dphidtheta_func(double Z, void *ctx)
{
  struct arc_length_ctx *actx = ctx;
  double *arc_memo = actx->arc_memo;
  double psi = actx->psi, rclose = actx->rclose, zmin = actx->zmin, arcL = actx->arcL, zmax = actx->zmax;

  // Get the integrand
  double integrand = 0.0;
  struct contour_ctx cctx = {
    .geo = actx->geo,
    .psi = psi,
    .ncall = 0,
    .last_R = rclose
  };
  integrand = dphidtheta_integrand(Z, &cctx);
  // Now multiply by fpol
  double R[4] = {0};
  double dRdZ[4] = {0};
  double dR[4] = {0};
  double dZ[4] = {0};
  int nr = gkyl_tok_geo_R_psiZ(actx->geo, psi, Z, 4, R, dRdZ, dR, dZ);
  double r_curr = nr == 1 ? R[0] : choose_closest(rclose, R, R, nr);
  double psi_fpol = psi;
  if ( (psi_fpol < actx->geo->fgrid.lower[0]) || (psi_fpol > actx->geo->fgrid.upper[0]) ) // F = F(psi_sep) in the SOL.
    psi_fpol = actx->geo->sibry;
  int idx = fmin(actx->geo->frange.lower[0] + (int) floor((psi_fpol - actx->geo->fgrid.lower[0])/actx->geo->fgrid.dx[0]), actx->geo->frange.upper[0]);
  long loc = gkyl_range_idx(&actx->geo->frange, &idx);
  const double *coeffs = gkyl_array_cfetch(actx->geo->fpoldg,loc);
  double fxc;
  gkyl_rect_grid_cell_center(&actx->geo->fgrid, &idx, &fxc);
  double fx = (psi_fpol-fxc)/(actx->geo->fgrid.dx[0]*0.5);
  double fpol = actx->geo->fbasis.eval_expand(&fx, coeffs);
  integrand = integrand*fpol;
  double darc_dtheta = actx->xpt_map_valid
    ? actx->xpt_map_darc_dtheta : actx->arcL_tot/(2.0*M_PI);
  integrand = integrand*darc_dtheta;
  return integrand;
}

static double
bmag_func(double r_curr, double Z, void *ctx)
{
  struct arc_length_ctx *actx = ctx;
  double *arc_memo = actx->arc_memo;
  double psi = actx->psi, rclose = actx->rclose, zmin = actx->zmin, arcL = actx->arcL, zmax = actx->zmax;
  // Calculate fpol
  double psi_fpol = psi;
  if ( (psi_fpol < actx->geo->fgrid.lower[0]) || (psi_fpol > actx->geo->fgrid.upper[0]) ) // F = F(psi_sep) in the SOL.
    psi_fpol = actx->geo->sibry;
  int idx = fmin(actx->geo->frange.lower[0] + (int) floor((psi_fpol - actx->geo->fgrid.lower[0])/actx->geo->fgrid.dx[0]), actx->geo->frange.upper[0]);
  long loc = gkyl_range_idx(&actx->geo->frange, &idx);
  const double *coeffs = gkyl_array_cfetch(actx->geo->fpoldg,loc);
  double fxc;
  gkyl_rect_grid_cell_center(&actx->geo->fgrid, &idx, &fxc);
  double fx = (psi_fpol-fxc)/(actx->geo->fgrid.dx[0]*0.5);
  double fpol = actx->geo->fbasis.eval_expand(&fx, coeffs);
  double Bphi = fpol/r_curr;
  double Br = 0.0, Bz = 0.0, bmag = 0.0;

  if (actx->geo->use_cubics) {
  double xn[2] = {r_curr, Z};
  double fout[3];
  actx->geo->efit->evf->eval_cubic_wgrad(0.0, xn, fout, actx->geo->efit->evf->ctx);
  double dpsidR = fout[1];
  double dpsidZ = fout[2];

  Br = 1.0/r_curr*dpsidZ;
  Bz = -1.0/r_curr*dpsidR;
  }
  else {
    int rzidx[2];
    int idxtemp = actx->geo->rzlocal.lower[0] + (int) floor((r_curr - actx->geo->rzgrid.lower[0])/actx->geo->rzgrid.dx[0]);
    idxtemp = GKYL_MIN2(idxtemp, actx->geo->rzlocal.upper[0]);
    idxtemp = GKYL_MAX2(idxtemp, actx->geo->rzlocal.lower[0]);
    rzidx[0] = idxtemp;
    idxtemp = actx->geo->rzlocal.lower[1] + (int) floor((Z - actx->geo->rzgrid.lower[1])/actx->geo->rzgrid.dx[1]);
    idxtemp = GKYL_MIN2(idxtemp, actx->geo->rzlocal.upper[1]);
    idxtemp = GKYL_MAX2(idxtemp, actx->geo->rzlocal.lower[1]);
    rzidx[1] = idxtemp;

    long loc = gkyl_range_idx((&actx->geo->rzlocal), rzidx);
    const double *psih = gkyl_array_cfetch(actx->geo->psiRZ, loc);

    double xc[2];
    gkyl_rect_grid_cell_center((&actx->geo->rzgrid), rzidx, xc);
    double x = (r_curr-xc[0])/(actx->geo->rzgrid.dx[0]*0.5);
    double y = (Z-xc[1])/(actx->geo->rzgrid.dx[1]*0.5);

    double dpsidx = 5.625*psih[8]*(2.0*x*SQ(y)-0.6666666666666666*x)+2.904737509655563*psih[7]*(SQ(y)-0.3333333333333333)+5.809475019311126*psih[6]*x*y+1.5*psih[3]*y+3.354101966249684*psih[4]*x+0.8660254037844386*psih[1];
    double dpsidy = 5.625*psih[8]*(2.0*SQ(x)*y-0.6666666666666666*y)+5.809475019311126*psih[7]*x*y+3.354101966249684*psih[5]*y+2.904737509655563*psih[6]*(SQ(x)-0.3333333333333333)+1.5*psih[3]*x+0.8660254037844386*psih[2];
    double dpsidR = dpsidx*2.0/actx->geo->rzgrid.dx[0];
    double dpsidZ = dpsidy*2.0/actx->geo->rzgrid.dx[1];
    Br = 1.0/r_curr*dpsidZ;
    Bz = -1.0/r_curr*dpsidR;
  }
  bmag = sqrt(Br*Br+Bz*Bz+Bphi*Bphi);
  return bmag;
}

static void
curlbhat_func(double psi, double r_curr, double Z, double phi, double *curlbhat, void *ctx)
{
  struct arc_length_ctx *actx = ctx;
  double *arc_memo = actx->arc_memo;
  // Calculate fpol and fpolprime
  // First get the location on the flux grid
  double psi_fpol = psi;
  if ( (psi_fpol < actx->geo->fgrid.lower[0]) || (psi_fpol > actx->geo->fgrid.upper[0]) ) // F = F(psi_sep) in the SOL.
    psi_fpol = actx->geo->sibry;
  int idx = fmin(actx->geo->frange.lower[0] + (int) floor((psi_fpol - actx->geo->fgrid.lower[0])/actx->geo->fgrid.dx[0]), actx->geo->frange.upper[0]);
  long loc = gkyl_range_idx(&actx->geo->frange, &idx);
  const double *coeffs = gkyl_array_cfetch(actx->geo->fpoldg,loc);
  double fxc;
  gkyl_rect_grid_cell_center(&actx->geo->fgrid, &idx, &fxc);
  double fx = (psi_fpol-fxc)/(actx->geo->fgrid.dx[0]*0.5);
  // Second calculate fpol and bphi 
  double fpol = actx->geo->fbasis.eval_expand(&fx, coeffs);
  double Bphi = fpol/r_curr;
  // Third calculate fpolprime
  coeffs = gkyl_array_cfetch(actx->geo->fpolprimedg,loc);
  double fpolprime = actx->geo->fbasis.eval_expand(&fx, coeffs);


  // Now calculate psi's various derivatives
  double Br = 0.0, Bz = 0.0, bmag = 0.0;
  double dpsidR = 0.0, dpsidZ = 0.0;
  double d2psidR2 = 0.0, d2psidZ2 = 0.0, d2psidRdZ = 0.0;
  double dBdR = 0.0, dBdZ = 0.0;
  double dBrdR = 0.0, dBrdZ = 0.0;
  double dBzdR = 0.0, dBzdZ = 0.0;

  if (actx->geo->use_cubics) {
    double xn[2] = {r_curr, Z};
    double fout[4];
    actx->geo->efit->evf->eval_cubic_wgrad(0.0, xn, fout, actx->geo->efit->evf->ctx);
    dpsidR = fout[1];
    dpsidZ = fout[2];
    actx->geo->efit->evf->eval_cubic_wgrad2(0.0, xn, fout, actx->geo->efit->evf->ctx);
    d2psidR2 = fout[1];
    d2psidZ2 = fout[2];
    d2psidRdZ = fout[3];
  }
  else {
    int rzidx[2];
    int idxtemp = actx->geo->rzlocal.lower[0] + (int) floor((r_curr - actx->geo->rzgrid.lower[0])/actx->geo->rzgrid.dx[0]);
    idxtemp = GKYL_MIN2(idxtemp, actx->geo->rzlocal.upper[0]);
    idxtemp = GKYL_MAX2(idxtemp, actx->geo->rzlocal.lower[0]);
    rzidx[0] = idxtemp;
    idxtemp = actx->geo->rzlocal.lower[1] + (int) floor((Z - actx->geo->rzgrid.lower[1])/actx->geo->rzgrid.dx[1]);
    idxtemp = GKYL_MIN2(idxtemp, actx->geo->rzlocal.upper[1]);
    idxtemp = GKYL_MAX2(idxtemp, actx->geo->rzlocal.lower[1]);
    rzidx[1] = idxtemp;

    long loc = gkyl_range_idx((&actx->geo->rzlocal), rzidx);
    const double *psih = gkyl_array_cfetch(actx->geo->psiRZ, loc);

    double xc[2];
    gkyl_rect_grid_cell_center((&actx->geo->rzgrid), rzidx, xc);
    double x = (r_curr-xc[0])/(actx->geo->rzgrid.dx[0]*0.5);
    double y = (Z-xc[1])/(actx->geo->rzgrid.dx[1]*0.5);

    double dpsidx = 5.625*psih[8]*(2.0*x*SQ(y)-0.6666666666666666*x)+2.904737509655563*psih[7]*(SQ(y)-0.3333333333333333)+5.809475019311126*psih[6]*x*y+1.5*psih[3]*y+3.354101966249684*psih[4]*x+0.8660254037844386*psih[1];
    double dpsidy = 5.625*psih[8]*(2.0*SQ(x)*y-0.6666666666666666*y)+5.809475019311126*psih[7]*x*y+3.354101966249684*psih[5]*y+2.904737509655563*psih[6]*(SQ(x)-0.3333333333333333)+1.5*psih[3]*x+0.8660254037844386*psih[2];
    dpsidR = dpsidx*2.0/actx->geo->rzgrid.dx[0];
    dpsidZ = dpsidy*2.0/actx->geo->rzgrid.dx[1];
    double d2psidx2 = 11.25*psih[8]*y*y+5.809475019311125*psih[6]*y-3.75*psih[8]+3.354101966249685*psih[4];
    double d2psidy2 = 11.25*psih[8]*x*x+5.809475019311125*psih[7]*x-3.75*psih[8]+3.354101966249685*psih[5];
    double d2psidxdy = 22.5*psih[8]*x*y+5.809475019311125*psih[7]*y+5.809475019311125*psih[6]*x+1.5*psih[3]; 
    d2psidR2 = d2psidx2*(2.0/actx->geo->rzgrid.dx[0])*(2.0/actx->geo->rzgrid.dx[0]);
    d2psidZ2 = d2psidy2*(2.0/actx->geo->rzgrid.dx[1])*(2.0/actx->geo->rzgrid.dx[1]);
    d2psidRdZ = d2psidxdy*(2.0/actx->geo->rzgrid.dx[0])*(2.0/actx->geo->rzgrid.dx[1]);
  }

  Br = 1.0/r_curr*dpsidZ;
  Bz = -1.0/r_curr*dpsidR;
  bmag = sqrt(Br*Br+Bz*Bz+Bphi*Bphi);

  dBrdR  = 1.0/r_curr*d2psidRdZ;
  dBrdZ  = 1.0/r_curr*d2psidZ2;
  dBzdR  = -1.0/r_curr*d2psidR2;
  dBzdZ  = -1.0/r_curr*d2psidRdZ;


  double dFdR = fpolprime*dpsidR;
  double dFdZ = fpolprime*dpsidZ;
  dBdR = 1/bmag*(Br*dBrdR + Bz*dBzdR + fpol/r_curr*(dFdR/r_curr - fpol/r_curr/r_curr));
  dBdZ = 1/bmag*(Br*dBrdZ + Bz*dBzdZ + fpol/r_curr*dFdZ);

  // Get the polar components (contravariant, upperscript components on tangent basis)
  double polar_comp[3] = {0.0};
  polar_comp[0] = 1.0/bmag*-1.0/r_curr*dFdZ - 1.0/bmag*1.0/bmag*dBdZ*fpol/r_curr; // R component ^1
  polar_comp[1] = 1.0/bmag*1.0/r_curr*(dBrdZ - dBzdR) + (-dBdR*Bz/r_curr + dBdZ*Br/r_curr); // Phi component ^2
  polar_comp[2] = 1.0/bmag*1.0/r_curr*dFdR + 1.0/bmag*1.0/bmag*dBdR*fpol/r_curr;

  // Convert to cartesian
  curlbhat[0] = polar_comp[0] * cos(phi) - polar_comp[1] * sin(phi) * r_curr;
  curlbhat[1] = polar_comp[0] * sin(phi) + polar_comp[1] * cos(phi) * r_curr;
  curlbhat[2] = polar_comp[2];
}






struct gkyl_tok_geo*
gkyl_tok_geo_new(const struct gkyl_efit_inp *inp, const struct gkyl_tok_geo_grid_inp *ginp)
{
  struct gkyl_tok_geo *geo = gkyl_malloc(sizeof(*geo));
  *geo = (struct gkyl_tok_geo) {};

  geo->efit = gkyl_efit_new(inp);

  geo->plate_spec = ginp->plate_spec;
  geo->plate_func_lower = ginp->plate_func_lower;
  geo->plate_func_upper = ginp->plate_func_upper;

  geo->rzbasis = geo->efit->rzbasis;
  geo->rzbasis_cubic = geo->efit->rzbasis_cubic;
  geo->rzgrid = geo->efit->rzgrid;
  geo->rzgrid_cubic = geo->efit->rzgrid_cubic;
  geo->psiRZ = gkyl_array_acquire(geo->efit->psizr);
  geo->psiRZ_cubic = gkyl_array_acquire(geo->efit->psizr_cubic);

  geo->num_rzbasis = geo->rzbasis.num_basis;
  geo->rzlocal = geo->efit->rzlocal;
  geo->rzlocal_ext = geo->efit->rzlocal_ext;
  geo->rzlocal_cubic = geo->efit->rzlocal_cubic;
  geo->rzlocal_cubic_ext = geo->efit->rzlocal_cubic_ext;
  geo->fgrid = geo->efit->fluxgrid;
  geo->fbasis = geo->efit->fluxbasis;
  geo->frange = geo->efit->fluxlocal;
  geo->frange_ext = geo->efit->fluxlocal_ext;
  geo->fpoldg= gkyl_array_acquire(geo->efit->fpolflux);
  geo->fpolprimedg= gkyl_array_acquire(geo->efit->fpolprimeflux);
  geo->qdg= gkyl_array_acquire(geo->efit->qflux);
  geo->sibry = geo->efit->sibry;
  geo->psisep = geo->efit->psisep;
  geo->zmaxis = geo->efit->zmaxis;

  geo->use_cubics = ginp->use_cubics;
  geo->use_hyperbolic_numbers = ginp->use_hyperbolic_numbers;
  geo->root_param.eps =
    ginp->root_param.eps > 0 ? ginp->root_param.eps : 1e-10;
  geo->root_param.max_iter =
    ginp->root_param.max_iter > 0 ? ginp->root_param.max_iter : 100;

  geo->quad_param.max_level =
    ginp->quad_param.max_levels > 0 ? ginp->quad_param.max_levels : 10;
  geo->quad_param.eps =
    ginp->quad_param.eps > 0 ? ginp->quad_param.eps : 1e-10;

  if (geo->use_cubics) {
    if(geo->use_hyperbolic_numbers)
      geo->calc_roots = calc_RdR_p3_hyperbolic;
    else
      geo->calc_roots = calc_RdR_p3;
    geo->calc_grad_psi = calc_grad_psi_p3;
  }
  else if (geo->efit->rzbasis.poly_order == 1) {
    geo->calc_roots = calc_RdR_p1;
    geo->calc_grad_psi = calc_grad_psi_p1;
  }
  else if (geo->efit->rzbasis.poly_order == 2){
    geo->calc_roots = calc_RdR_p2_tensor_nrc;
    geo->calc_grad_psi = calc_grad_psi_p2_tensor;
  }

  geo->stat = (struct gkyl_tok_geo_stat) { };

  
  return geo;
}

double
gkyl_tok_geo_integrate_psi_contour(const struct gkyl_tok_geo *geo, double psi,
  double zmin, double zmax, double rclose)
{
  return integrate_psi_contour_memo(geo, psi, zmin, zmax, rclose,
    false, false, 0);
}

int
gkyl_tok_geo_R_psiZ(const struct gkyl_tok_geo *geo, double psi, double Z, int nmaxroots,
  double *R, double *dRdZ, double *dR, double *dZ)
{
  if(geo->use_cubics)
    return R_psiZ_cubic(geo, psi, Z, nmaxroots, R, dRdZ, dR, dZ);
  else
    return R_psiZ(geo, psi, Z, nmaxroots, R, dRdZ, dR, dZ);
}

void gkyl_tok_geo_calc(struct gk_geometry* up, struct gkyl_range *nrange, struct gkyl_tok_geo *geo, 
  struct gkyl_tok_geo_grid_inp *inp, struct gkyl_position_map *position_map)
{

  geo->rleft = inp->rleft;
  geo->rright = inp->rright;

  geo->inexact_roots = inp->inexact_roots;

  geo->rmax = inp->rmax;
  geo->rmin = inp->rmin;

  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
  
  double dtheta = inp->cgrid.dx[TH_IDX],
    dpsi = inp->cgrid.dx[PSI_IDX],
    dalpha = inp->cgrid.dx[AL_IDX];
  
  double theta_lo = up->grid.lower[TH_IDX] + (up->local.lower[TH_IDX] - up->global.lower[TH_IDX])*up->grid.dx[TH_IDX],
    psi_lo = up->grid.lower[PSI_IDX] + (up->local.lower[PSI_IDX] - up->global.lower[PSI_IDX])*up->grid.dx[PSI_IDX],
    alpha_lo = up->grid.lower[AL_IDX] + (up->local.lower[AL_IDX] - up->global.lower[AL_IDX])*up->grid.dx[AL_IDX];
    
  double dx_fact = up->basis.poly_order == 1.0/up->basis.poly_order;
  dtheta *= dx_fact; dpsi *= dx_fact; dalpha *= dx_fact;

  double rclose = inp->rclose;
  double rright = inp->rright;
  double rleft = inp->rleft;

  int nzcells;
  if(geo->use_cubics)
    nzcells = geo->rzgrid_cubic.cells[1];
  else
    nzcells = geo->rzgrid.cells[1];
  double *arc_memo = gkyl_malloc(sizeof(double[nzcells]));
  double *arc_memo_left = gkyl_malloc(sizeof(double[nzcells]));
  double *arc_memo_right = gkyl_malloc(sizeof(double[nzcells]));
  int sep_trace_capacity = 16*nzcells+1;
  double *sep_trace_r = gkyl_malloc(sizeof(double[sep_trace_capacity]));
  double *sep_trace_z = gkyl_malloc(sizeof(double[sep_trace_capacity]));
  double *sep_trace_s = gkyl_malloc(sizeof(double[sep_trace_capacity]));
  // The marched theta ladder owns its own block: how many rungs it needs is
  // discovered while marching, not known here, so it cannot be sized with the
  // rest.  It is freed alongside this one.
  double *ordered_trace_storage =
    gkyl_malloc(sizeof(double)*8*(size_t) sep_trace_capacity);

  struct arc_length_ctx arc_ctx = {
    .geo = geo,
    .arc_memo = arc_memo,
    .arc_memo_right = arc_memo_right,
    .arc_memo_left = arc_memo_left,
    .sep_trace_r = sep_trace_r,
    .sep_trace_z = sep_trace_z,
    .sep_trace_s = sep_trace_s,
    .sep_trace_capacity = sep_trace_capacity,
    .far_trace_r = ordered_trace_storage,
    .far_trace_z = ordered_trace_storage+sep_trace_capacity,
    .far_trace_s = ordered_trace_storage+2*sep_trace_capacity,
    .trace_corr_v = ordered_trace_storage+3*sep_trace_capacity,
    .map_trace_r = ordered_trace_storage+4*sep_trace_capacity,
    .map_trace_z = ordered_trace_storage+5*sep_trace_capacity,
    .map_trace_s = ordered_trace_storage+6*sep_trace_capacity,
    .map_trace_phi = ordered_trace_storage+7*sep_trace_capacity,
    .ext_ladder_w = NULL,
    .ftype = inp->ftype,
    .zmaxis = geo->zmaxis
  };
  struct plate_ctx pctx = {
    .geo = geo
  };
  tok_init_xpt_ray_target(inp, position_map, &arc_ctx);

  int cidx[3] = { 0 };
  for (int ia=nrange->lower[AL_IDX]; ia<nrange->lower[AL_IDX]+1; ++ia){
    cidx[AL_IDX] = ia;
    double alpha_curr = alpha_lo + ia*dalpha;

    for (int ip=nrange->lower[PSI_IDX]; ip<=nrange->upper[PSI_IDX]; ++ip) {
      double psi_curr = psi_lo + ip*dpsi;

      // Non-uniform psi. Finite differences are calculated in calc_metric.c
      position_map->maps[0](0.0, &psi_curr,  &psi_curr,  position_map->ctxs[0]);

      double darcL, arcL_curr, arcL_lo;

      // For double null blocks this should set arc_ctx :
      // zmin, zmax, rclose, arcL_tot for all blocks. No left and right
      // For a full core case:
      // also set phi_right and arcL_right
      // For a single null case:
      // also set zmin_left and zmin_right 
      if (tok_xpt_mapping_requested(inp))
        tok_prepare_ordered_map(inp, &arc_ctx, psi_curr);
      else
        tok_find_endpoints(inp, geo, &arc_ctx, &pctx, psi_curr, alpha_curr,
          arc_memo, arc_memo_left, arc_memo_right);

      darcL = arc_ctx.arcL_tot/(up->basis.poly_order*inp->cgrid.cells[TH_IDX])
        * (inp->cgrid.upper[TH_IDX] - inp->cgrid.lower[TH_IDX])/2/M_PI;
      // At the beginning of each theta loop we need to reset things.
      cidx[PSI_IDX] = ip;
      arcL_curr = 0.0;
      arcL_lo = (theta_lo + M_PI)/2/M_PI*arc_ctx.arcL_tot;
      double ridders_min, ridders_max;
      // Set node coordinates.
      for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
        int it_delta = 0;
        arcL_curr = arcL_lo + it*darcL;
        double theta_curr = arcL_curr*(2*M_PI/arc_ctx.arcL_tot) - M_PI ; 

        // Calculate derivatives using finite difference for ddtheta,
        // as well as transform the computational coordiante to the non-uniform field-aligned value

        // We cannot do non-uniform alpha because we are modeling axisymmetric systems
        // Non-uniform theta
        double Theta_curr;
        position_map->maps[2](0.0, &theta_curr,  &Theta_curr,  position_map->ctxs[2]);
        theta_curr = Theta_curr;

        struct tok_ordered_point ordered = { 0.0 };
        bool ordered_mapping = tok_ordered_map_lookup(inp, &arc_ctx,
          theta_curr, alpha_curr, &ordered);
        if (tok_xpt_mapping_requested(inp) && !ordered_mapping) {
          fprintf(stderr,
            "TOK_ORDERED_MAP lookup failed ftype=%d psi=%.17g theta=%.17g\n",
            inp->ftype, psi_curr, theta_curr);
          abort();
        }
        double r_curr = 0.0, z_curr = 0.0, phi_curr = 0.0;
        if (ordered_mapping) {
          r_curr = ordered.r; z_curr = ordered.z; phi_curr = ordered.phi;
        }
        else {
          arcL_curr = tok_xpt_theta_to_arc(inp, &arc_ctx, theta_curr);

          tok_set_ridders(inp, &arc_ctx, psi_curr, arcL_curr, &rclose, &ridders_min, &ridders_max);

          struct gkyl_qr_res res = gkyl_ridders(arc_length_func, &arc_ctx,
            arc_ctx.zmin, arc_ctx.zmax, ridders_min, ridders_max,
            geo->root_param.max_iter, 1e-10);
          z_curr = res.res;
          ((struct gkyl_tok_geo *)geo)->stat.nroot_cont_calls += res.nevals;

        if (inp->half_domain) { //Alternative for half domain
          if (tok_geo_same_flux(psi_curr, geo->psisep)) {
            if (it == nrange->upper[TH_IDX] && (up->local.upper[TH_IDX]== up->global.upper[TH_IDX])) {
              if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_R || inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L || inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO|| inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID)
                z_curr = geo->efit->Zxpt[0];
              else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO || inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID || inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE)
                z_curr = geo->efit->Zxpt[0];
            }
            if (it == nrange->lower[TH_IDX] && (up->local.lower[TH_IDX]== up->global.lower[TH_IDX])) {
              if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_L || inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R || inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID|| inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO)
                z_curr = geo->efit->Zxpt[0];
              else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP || inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID || inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE)
                z_curr = geo->efit->Zxpt[0];
            }
          }
        }
        else { // For full domain
          if (tok_geo_same_flux(psi_curr, geo->psisep)) {
            if (it == nrange->upper[TH_IDX] && (up->local.upper[TH_IDX]== up->global.upper[TH_IDX])) {
              if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_UP_L || inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R || inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID || inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_UP)
                z_curr = (geo->efit->num_xpts > 1 ? geo->efit->Zxpt[1] : geo->efit->Zxpt[0]);
              else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_R || inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L || inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_LO|| inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID)
                z_curr = geo->efit->Zxpt[0];
              else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_LO || inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID || inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE)
                z_curr = geo->efit->Zxpt[0];
            }
            if (it == nrange->lower[TH_IDX] && (up->local.lower[TH_IDX]== up->global.lower[TH_IDX])) {
              if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_UP_R || inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L || inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_UP|| inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_MID)
                z_curr = (geo->efit->num_xpts > 1 ? geo->efit->Zxpt[1] : geo->efit->Zxpt[0]);
              else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_PF_LO_L || inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R || inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_OUT_MID|| inp->ftype == GKYL_GEOMETRY_TOKAMAK_DN_SOL_IN_LO)
                z_curr = geo->efit->Zxpt[0];
              else if (inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_UP || inp->ftype == GKYL_GEOMETRY_TOKAMAK_LSN_SOL_MID || inp->ftype == GKYL_GEOMETRY_TOKAMAK_CORE)
                z_curr = geo->efit->Zxpt[0];
            }
          }
        }

        double sep_r_curr = 0.0;
        bool at_sep_trace = tok_half_domain_sep_rz(inp, &arc_ctx,
          theta_curr, &sep_r_curr, &z_curr);

        double fixed_root_ref = rclose;
        bool at_fixed_edge = tok_xpt_at_fixed_edge(inp, &arc_ctx, it,
          nrange,
          up->local.lower[TH_IDX] == up->global.lower[TH_IDX],
          up->local.upper[TH_IDX] == up->global.upper[TH_IDX],
          &z_curr, &fixed_root_ref);
        bool at_xpt_anchor = arc_ctx.xpt_anchor_valid &&
          tok_xpt_at_seam(inp, it, nrange,
            up->local.lower[TH_IDX] == up->global.lower[TH_IDX],
            up->local.upper[TH_IDX] == up->global.upper[TH_IDX]);
        if (at_xpt_anchor)
          z_curr = arc_ctx.xpt_anchor_z;

        double R[4] = { 0 }, dRdZ[4] = { 0 };
        double dR[4] = { 0 }, dZ[4] = { 0 };
        int nr = gkyl_tok_geo_R_psiZ(geo, psi_curr, z_curr, 4, R, dRdZ, dR, dZ);
        double root_ref = at_sep_trace ? sep_r_curr
          : (at_xpt_anchor ? arc_ctx.xpt_anchor_r
            : (at_fixed_edge ? fixed_root_ref : rclose));
        r_curr = choose_closest(root_ref, R, R, nr);
        double drdz_curr = choose_closest(root_ref, R, dRdZ, nr);
        double dr_curr = choose_closest(root_ref, R, dR, nr);
        double dz_curr = choose_closest(root_ref, R, dZ, nr);
        if (at_xpt_anchor)
          r_curr = arc_ctx.xpt_anchor_r;
        if (at_sep_trace)
          r_curr = sep_r_curr;

        if (tok_geo_same_flux(psi_curr, geo->psisep)) {
          if (z_curr == geo->efit->Zxpt[0]) {
            nr = 1;
            r_curr = geo->efit->Rxpt[0];
          }
          if (z_curr == (geo->efit->num_xpts > 1 ? geo->efit->Zxpt[1] : geo->efit->Zxpt[0])) {
            nr = 1;
            r_curr = (geo->efit->num_xpts > 1 ? geo->efit->Rxpt[1] : geo->efit->Rxpt[0]);
          }
        }

        if (nr==0) {
          printf(" ip = %d, it = %d, ia = %d\n", ip, it, ia);
          printf("Block Type = %d | Failed to find a root at psi = %g, Z = %1.16f\n", inp->ftype, psi_curr, z_curr);
          assert(false);
        }

          phi_curr = phi_func(alpha_curr, z_curr, &arc_ctx);
        }
        cidx[TH_IDX] = it;
        double *mc2p_n = gkyl_array_fetch(up->geo_corn.mc2p_nodal, gkyl_range_idx(nrange, cidx));
        double *mc2nu_n = gkyl_array_fetch(up->geo_corn.mc2nu_pos_nodal, gkyl_range_idx(nrange, cidx));
        double *bmag_n = gkyl_array_fetch(up->geo_corn.bmag_nodal, gkyl_range_idx(nrange, cidx));

        mc2p_n[X_IDX] = r_curr;
        mc2p_n[Y_IDX] = z_curr;
        mc2p_n[Z_IDX] = phi_curr;
        mc2nu_n[X_IDX] = psi_curr;
        mc2nu_n[Y_IDX] = alpha_curr;
        mc2nu_n[Z_IDX] = theta_curr;
        bmag_n[0] = bmag_func(r_curr, z_curr, &arc_ctx);
      }
    }
  }

  // Populate other alpha indices by using axisymmetry
  for (int ia=nrange->lower[AL_IDX]+1; ia<=nrange->upper[AL_IDX]; ++ia){
    cidx[AL_IDX] = ia;
    double alpha_curr = alpha_lo + ia*dalpha;
    double alpha_donor = alpha_lo + nrange->lower[AL_IDX]*dalpha;
    double alpha_diff = alpha_curr -  alpha_donor;
    for (int ip=nrange->lower[PSI_IDX]; ip<=nrange->upper[PSI_IDX]; ++ip) {
      cidx[PSI_IDX] = ip;
      for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
        cidx[TH_IDX] = it;

        double *mc2p_n = gkyl_array_fetch(up->geo_corn.mc2p_nodal, gkyl_range_idx(nrange, cidx));
        double *mc2nu_n = gkyl_array_fetch(up->geo_corn.mc2nu_pos_nodal, gkyl_range_idx(nrange, cidx));
        double *bmag_n = gkyl_array_fetch(up->geo_corn.bmag_nodal, gkyl_range_idx(nrange, cidx));

        int donor_cidx[3] ;
        donor_cidx[AL_IDX] = nrange->lower[AL_IDX];
        donor_cidx[PSI_IDX] = ip;
        donor_cidx[TH_IDX] = it;

        double *donor_mc2p_n = gkyl_array_fetch(up->geo_corn.mc2p_nodal, gkyl_range_idx(nrange, donor_cidx));
        double *donor_mc2nu_n = gkyl_array_fetch(up->geo_corn.mc2nu_pos_nodal, gkyl_range_idx(nrange, donor_cidx));
        double *donor_bmag_n = gkyl_array_fetch(up->geo_corn.bmag_nodal, gkyl_range_idx(nrange, donor_cidx));

        mc2p_n[X_IDX] = donor_mc2p_n[X_IDX];
        mc2p_n[Y_IDX] = donor_mc2p_n[Y_IDX];
        mc2p_n[Z_IDX] = donor_mc2p_n[Z_IDX] + alpha_diff;
        mc2nu_n[X_IDX] = donor_mc2nu_n[X_IDX];
        mc2nu_n[Y_IDX] = donor_mc2nu_n[AL_IDX] + alpha_diff;
        mc2nu_n[Z_IDX] = donor_mc2nu_n[Z_IDX];
        bmag_n[0] = donor_bmag_n[0];
      }
    }
  }

  // A folded cell is a grid that is wrong, not merely poor: the map from
  // computational to physical coordinates has reversed orientation there, so
  // the Jacobian changes sign inside the block. Such a block used to be written
  // out and "pass", and was only caught post-hoc by check_grid_quality.py.
  // Check it here instead, on the corner nodes, and fail the build.
  //
  // Signed area of each (psi,theta) quad in the R-Z plane, at every alpha. The
  // sign convention is per block, so compare against the block's own dominant
  // sign rather than assuming positive.
  {
    int nfold = 0, nquad = 0;
    double worst = 0.0, worst_r = 0.0, worst_z = 0.0;
    int worst_ip = -1, worst_it = -1;
    long npos = 0, nneg = 0;
    int cidx0[3] = { 0, 0, 0 };
    for (int pass=0; pass<2; ++pass) {
      for (int ia=nrange->lower[AL_IDX]; ia<=nrange->upper[AL_IDX]; ++ia) {
        for (int ip=nrange->lower[PSI_IDX]; ip<nrange->upper[PSI_IDX]; ++ip) {
          for (int it=nrange->lower[TH_IDX]; it<nrange->upper[TH_IDX]; ++it) {
            double rq[4], zq[4];
            const int dp[4] = { 0, 1, 1, 0 }, dt[4] = { 0, 0, 1, 1 };
            for (int k=0; k<4; ++k) {
              cidx0[PSI_IDX] = ip+dp[k];
              cidx0[AL_IDX] = ia;
              cidx0[TH_IDX] = it+dt[k];
              const double *p = gkyl_array_cfetch(up->geo_corn.mc2p_nodal,
                gkyl_range_idx(nrange, cidx0));
              rq[k] = p[X_IDX]; zq[k] = p[Y_IDX];
            }
            double a = 0.0;
            for (int k=0; k<4; ++k) {
              int k1 = (k+1)%4;
              a += rq[k]*zq[k1]-rq[k1]*zq[k];
            }
            a *= 0.5;
            if (pass == 0) {
              if (a > 0.0) ++npos; else if (a < 0.0) ++nneg;
              continue;
            }
            ++nquad;
            double signed_a = npos >= nneg ? a : -a;
            if (signed_a <= 0.0) {
              ++nfold;
              if (signed_a < worst) {
                worst = signed_a; worst_ip = ip; worst_it = it;
                worst_r = rq[0]; worst_z = zq[0];
              }
            }
          }
        }
      }
    }
    if (nfold > 0) {
      fprintf(stderr,
        "TOK_GEO_FOLDED_CELLS ftype=%d nfold=%d of %d quads "
        "worst_area=%.17g at (ip=%d,it=%d) node=(%.17g,%.17g)\n",
        inp->ftype, nfold, nquad, worst, worst_ip, worst_it, worst_r, worst_z);
      const char *allow = getenv("GKYL_TOK_ALLOW_FOLDED_CELLS");
      if (!(allow && allow[0] != '\0' && allow[0] != '0')) {
        fprintf(stderr,
          "TOK_GEO_FOLDED_CELLS ftype=%d aborting: the coordinate map reverses "
          "orientation inside this block, so its Jacobian changes sign. Set "
          "GKYL_TOK_ALLOW_FOLDED_CELLS=1 to emit the grid anyway.\n",
          inp->ftype);
        abort();
      }
    }
  }

  // Consecutive flux surfaces crossing is the other way a block can be wrong
  // while still "passing": every quad keeps its orientation, so the fold check
  // above sees nothing, but the radial ordering reverses and the grid is
  // physically meaningless there. 205004's DN_SOL_OUT_MID emitted 44 reversed
  // radial steps (min cos -0.995) and was caught only post-hoc by
  // check_grid_quality.py.
  //
  // Same criterion as that script: walk the radial node line at fixed theta and
  // compare successive step directions. A pair counts as a reversal only well
  // past orthogonal, so ordinary curvature near an X point is not reported.
  // Like the fold count above this sums over both alpha nodes, so it is 2x what
  // check_grid_quality.py reports from a single alpha slice.
  {
    const double reversal_cos = -0.5;
    int nrev = 0, npair = 0;
    double worst_cos = 1.0, worst_r = 0.0, worst_z = 0.0;
    int worst_ip = -1, worst_it = -1;
    int cidx0[3] = { 0, 0, 0 };
    for (int ia=nrange->lower[AL_IDX]; ia<=nrange->upper[AL_IDX]; ++ia) {
      for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
        for (int ip=nrange->lower[PSI_IDX]; ip+2<=nrange->upper[PSI_IDX]; ++ip) {
          double rp[3], zp[3];
          for (int k=0; k<3; ++k) {
            cidx0[PSI_IDX] = ip+k;
            cidx0[AL_IDX] = ia;
            cidx0[TH_IDX] = it;
            const double *p = gkyl_array_cfetch(up->geo_corn.mc2p_nodal,
              gkyl_range_idx(nrange, cidx0));
            rp[k] = p[X_IDX]; zp[k] = p[Y_IDX];
          }
          double ar = rp[1]-rp[0], az = zp[1]-zp[0];
          double br = rp[2]-rp[1], bz = zp[2]-zp[1];
          double la = hypot(ar, az), lb = hypot(br, bz);
          if (!(la > 0.0) || !(lb > 0.0))
            continue;
          ++npair;
          double c = (ar*br+az*bz)/(la*lb);
          if (c < reversal_cos) {
            ++nrev;
            if (c < worst_cos) {
              worst_cos = c; worst_ip = ip; worst_it = it;
              worst_r = rp[1]; worst_z = zp[1];
            }
          }
        }
      }
    }
    if (nrev > 0) {
      fprintf(stderr,
        "TOK_GEO_SURFACE_CROSS ftype=%d nrev=%d of %d radial step pairs "
        "worst_cos=%.17g at (ip=%d,it=%d) node=(%.17g,%.17g)\n",
        inp->ftype, nrev, npair, worst_cos, worst_ip, worst_it,
        worst_r, worst_z);
      const char *allow = getenv("GKYL_TOK_ALLOW_SURFACE_CROSS");
      if (!(allow && allow[0] != '\0' && allow[0] != '0')) {
        fprintf(stderr,
          "TOK_GEO_SURFACE_CROSS ftype=%d aborting: consecutive flux surfaces "
          "cross inside this block, so its radial ordering reverses. Set "
          "GKYL_TOK_ALLOW_SURFACE_CROSS=1 to emit the grid anyway.\n",
          inp->ftype);
        abort();
      }
    }
  }

  struct gkyl_nodal_ops *n2m =  gkyl_nodal_ops_new(&inp->cbasis, &inp->cgrid, false);
  gkyl_nodal_ops_n2m(n2m, &inp->cbasis, &inp->cgrid, nrange, &up->local, 3, up->geo_corn.mc2p_nodal, up->geo_corn.mc2p, false);
  gkyl_nodal_ops_n2m(n2m, &inp->cbasis, &inp->cgrid, nrange, &up->local, 3, up->geo_corn.mc2nu_pos_nodal, up->geo_corn.mc2nu_pos, false);
  gkyl_nodal_ops_n2m(n2m, &inp->cbasis, &inp->cgrid, nrange, &up->local, 1, up->geo_corn.bmag_nodal, up->geo_corn.bmag, false);
  gkyl_nodal_ops_release(n2m);

  // Need 1/B for LBO collisions, computed weakly.
  gkyl_dg_inv_op_range(&inp->cbasis, 0, up->geo_corn.bmag_inv, 0, up->geo_corn.bmag, &up->local);

  gkyl_free(arc_memo);
  gkyl_free(arc_memo_left);
  gkyl_free(arc_memo_right);
  gkyl_free(sep_trace_r);
  gkyl_free(sep_trace_z);
  gkyl_free(sep_trace_s);
  gkyl_free(arc_ctx.ext_ladder_w);
  gkyl_free(ordered_trace_storage);
}

void gkyl_tok_geo_calc_interior(struct gk_geometry* up, struct gkyl_range *nrange, double dzc[3], 
    struct gkyl_tok_geo *geo, struct gkyl_tok_geo_grid_inp *inp, struct gkyl_position_map *position_map)
{

  geo->rleft = inp->rleft;
  geo->rright = inp->rright;

  geo->inexact_roots = inp->inexact_roots;

  geo->rmax = inp->rmax;
  geo->rmin = inp->rmin;

  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  enum { X_IDX, Y_IDX, Z_IDX }; // arrangement of cartesian coordinates
  
  double dtheta = inp->cgrid.dx[TH_IDX],
    dpsi = inp->cgrid.dx[PSI_IDX],
    dalpha = inp->cgrid.dx[AL_IDX];
  
  double theta_lo = up->grid.lower[TH_IDX] + (up->local.lower[TH_IDX] - up->global.lower[TH_IDX])*up->grid.dx[TH_IDX],
    psi_lo = up->grid.lower[PSI_IDX] + (up->local.lower[PSI_IDX] - up->global.lower[PSI_IDX])*up->grid.dx[PSI_IDX],
    alpha_lo = up->grid.lower[AL_IDX] + (up->local.lower[AL_IDX] - up->global.lower[AL_IDX])*up->grid.dx[AL_IDX];

  double dels[2] = {1.0/sqrt(3), 1.0-1.0/sqrt(3) };
  theta_lo = theta_lo + dels[1]*dtheta/2.0;
  psi_lo = psi_lo + dels[1]*dpsi/2.0;
  alpha_lo = alpha_lo + dels[1]*dalpha/2.0;
    
  double dx_fact = up->basis.poly_order == 1.0/up->basis.poly_order;
  dtheta *= dx_fact; dpsi *= dx_fact; dalpha *= dx_fact;

  // used for finite differences 
  double delta_alpha = dalpha*1e-2;
  double delta_psi = dpsi*1e-2;
  double delta_theta = dtheta*1e-2;
  dzc[0] = delta_psi;
  dzc[1] = delta_alpha;
  dzc[2] = delta_theta;
  int modifiers[5] = {0, -1, 1, -2, 2};

  double rclose = inp->rclose;
  double rright = inp->rright;
  double rleft = inp->rleft;

  int nzcells;
  if(geo->use_cubics)
    nzcells = geo->rzgrid_cubic.cells[1];
  else
    nzcells = geo->rzgrid.cells[1];
  double *arc_memo = gkyl_malloc(sizeof(double[nzcells]));
  double *arc_memo_left = gkyl_malloc(sizeof(double[nzcells]));
  double *arc_memo_right = gkyl_malloc(sizeof(double[nzcells]));
  int sep_trace_capacity = 16*nzcells+1;
  double *sep_trace_r = gkyl_malloc(sizeof(double[sep_trace_capacity]));
  double *sep_trace_z = gkyl_malloc(sizeof(double[sep_trace_capacity]));
  double *sep_trace_s = gkyl_malloc(sizeof(double[sep_trace_capacity]));
  // The marched theta ladder owns its own block: how many rungs it needs is
  // discovered while marching, not known here, so it cannot be sized with the
  // rest.  It is freed alongside this one.
  double *ordered_trace_storage =
    gkyl_malloc(sizeof(double)*8*(size_t) sep_trace_capacity);

  struct arc_length_ctx arc_ctx = {
    .geo = geo,
    .arc_memo = arc_memo,
    .arc_memo_right = arc_memo_right,
    .arc_memo_left = arc_memo_left,
    .sep_trace_r = sep_trace_r,
    .sep_trace_z = sep_trace_z,
    .sep_trace_s = sep_trace_s,
    .sep_trace_capacity = sep_trace_capacity,
    .far_trace_r = ordered_trace_storage,
    .far_trace_z = ordered_trace_storage+sep_trace_capacity,
    .far_trace_s = ordered_trace_storage+2*sep_trace_capacity,
    .trace_corr_v = ordered_trace_storage+3*sep_trace_capacity,
    .map_trace_r = ordered_trace_storage+4*sep_trace_capacity,
    .map_trace_z = ordered_trace_storage+5*sep_trace_capacity,
    .map_trace_s = ordered_trace_storage+6*sep_trace_capacity,
    .map_trace_phi = ordered_trace_storage+7*sep_trace_capacity,
    .ext_ladder_w = NULL,
    .ftype = inp->ftype,
    .zmaxis = geo->zmaxis
  };
  struct plate_ctx pctx = {
    .geo = geo
  };
  tok_init_xpt_ray_target(inp, position_map, &arc_ctx);

  // Temporary array to store nodal q profile.
  struct gkyl_array *qprofile_nodal = gkyl_array_new(GKYL_DOUBLE, up->geo_int.bmag_nodal->ncomp, up->geo_int.bmag_nodal->size);

  int cidx[3] = { 0 };
  for(int ia=nrange->lower[AL_IDX]; ia<nrange->lower[AL_IDX]+1; ++ia){
    cidx[AL_IDX] = ia;
    double alpha_curr = calc_running_coord(alpha_lo, ia-nrange->lower[AL_IDX], dalpha);

    for (int ip=nrange->lower[PSI_IDX]; ip<=nrange->upper[PSI_IDX]; ++ip) {
      int ip_delta_max = 3;
      for (int ip_delta = 0; ip_delta < ip_delta_max; ip_delta++) {

        double psi_curr = calc_running_coord(psi_lo, ip-nrange->lower[PSI_IDX], dpsi) + modifiers[ip_delta]*delta_psi;

        // Non-uniform psi. Finite differences are calculated in calc_metric.c
        double Psi_curr;
        position_map->maps[0](0.0, &psi_curr,  &Psi_curr,  position_map->ctxs[0]);
        double dPsi_dpsi = gkyl_position_map_slope(position_map, 0, psi_curr, delta_psi, ip, nrange);
        psi_curr = Psi_curr;
        
        double darcL, arcL_curr, arcL_lo;

        // For double null blocks this should set arc_ctx :
        // zmin, zmax, rclose, arcL_tot for all blocks. No left and right
        // For a full core case:
        // also set phi_right and arcL_right
        // For a single null case:
        // also set zmin_left and zmin_right 
        double qprofile = 0.0;
        if (tok_xpt_mapping_requested(inp)) {
          if (!inp->half_domain) {
            // The ordered map changes only the coordinate construction.  Use
            // the established contour setup, with the feature flags disabled
            // on a local copy, to preserve the pre-feature q profile exactly.
            struct gkyl_tok_geo_grid_inp legacy_inp = *inp;
            struct arc_length_ctx qctx = arc_ctx;
            legacy_inp.straight_xpt_ray = false;
            legacy_inp.straight_core_xpt_ray = false;
            tok_find_endpoints(&legacy_inp, geo, &qctx, &pctx,
              psi_curr, alpha_curr, arc_memo, arc_memo_left,
              arc_memo_right);
            qprofile = qprofile_func(&qctx);
          }
          tok_prepare_ordered_map(inp, &arc_ctx, psi_curr);
        }
        else {
          tok_find_endpoints(inp, geo, &arc_ctx, &pctx, psi_curr, alpha_curr,
            arc_memo, arc_memo_left, arc_memo_right);
          qprofile = qprofile_func(&arc_ctx);
        }

        // Calculate the q profile
        // qhat = - F(psi) * s(psi) / (R * grad(psi))
        // q = integral_0^2pi qhat dtheta
        //   = F(psi)*s(psi) * integral 1/(R*grad(psi)) dl
        //   = 1/s(psi) * integral (dphidtheta) ; dphidtheta = F(psi)/(R*grad(psi))
        // The legacy half-domain block types return zero from qprofile_func;
        // the full-domain branch above evaluates the unchanged legacy path.

        darcL = arc_ctx.arcL_tot/(up->basis.poly_order*inp->cgrid.cells[TH_IDX])
          * (inp->cgrid.upper[TH_IDX] - inp->cgrid.lower[TH_IDX])/2/M_PI;
        // at the beginning of each theta loop we need to reset things
        cidx[PSI_IDX] = ip;
        arcL_curr = 0.0;
        arcL_lo = (theta_lo + M_PI)/2/M_PI*arc_ctx.arcL_tot;
        double ridders_min, ridders_max;

        for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
          arcL_curr = calc_running_coord(arcL_lo, it-nrange->lower[TH_IDX], darcL);
          double theta_curr = arcL_curr*(2*M_PI/arc_ctx.arcL_tot) - M_PI ; 

          // Calculate derivatives using finite difference for ddtheta,
          // as well as transform the computational coordiante to the non-uniform field-aligned value

          // We cannot do non-uniform alpha because we are modeling axisymmetric systems
          // Non-uniform theta
          double Theta_curr;
          position_map->maps[2](0.0, &theta_curr,  &Theta_curr,  position_map->ctxs[2]);
          double dTheta_dtheta = gkyl_position_map_slope(position_map, 2, theta_curr,\
            delta_theta, it, nrange);
          theta_curr = Theta_curr;

          struct tok_ordered_point ordered = { 0.0 };
          bool ordered_mapping = tok_ordered_map_lookup(inp, &arc_ctx,
            theta_curr, alpha_curr, &ordered);
          if (tok_xpt_mapping_requested(inp) && !ordered_mapping) {
            fprintf(stderr,
              "TOK_ORDERED_MAP lookup failed ftype=%d psi=%.17g theta=%.17g\n",
              inp->ftype, psi_curr, theta_curr);
            abort();
          }
          double r_curr = 0.0, z_curr = 0.0, phi_curr = 0.0;
          double drdz_curr = 0.0, dr_curr = 0.0, dz_curr = 0.0;
          if (ordered_mapping) {
            r_curr = ordered.r; z_curr = ordered.z; phi_curr = ordered.phi;
          }
          else {
            arcL_curr = tok_xpt_theta_to_arc(inp, &arc_ctx, theta_curr);

            tok_set_ridders(inp, &arc_ctx, psi_curr, arcL_curr, &rclose, &ridders_min, &ridders_max);

          bool trace_this_block = tok_geo_trace_enabled() && tok_geo_trace_block_selected();
          if (trace_this_block && tok_geo_trace_int_env("GKYL_TOK_GEO_TRACE_BEFORE", 1))
            tok_geo_trace_surface_row("before_ridders", inp->ftype, -1, ip,
              it, ia, ip_delta, psi_curr, alpha_curr, theta_curr, arcL_curr,
              arc_ctx.zmin, arc_ctx.zmax, rclose, ridders_min, ridders_max,
              0.0, -1, 0.0);

          clock_t trace_t0 = trace_this_block ? clock() : 0;
          struct gkyl_qr_res res = gkyl_ridders(arc_length_func, &arc_ctx,
            arc_ctx.zmin, arc_ctx.zmax, ridders_min, ridders_max,
            geo->root_param.max_iter, 1e-10);
          double trace_elapsed = trace_this_block ?
            ((double) (clock() - trace_t0))/CLOCKS_PER_SEC : 0.0;
          z_curr = res.res;
          ((struct gkyl_tok_geo *)geo)->stat.nroot_cont_calls += res.nevals;

          if (trace_this_block) {
            int neval_threshold = tok_geo_trace_int_env("GKYL_TOK_GEO_TRACE_NEVAL_THRESHOLD", 25);
            double time_threshold = tok_geo_trace_double_env("GKYL_TOK_GEO_TRACE_TIME_THRESHOLD", 0.01);
            if (res.nevals >= neval_threshold || trace_elapsed >= time_threshold ||
              tok_geo_trace_int_env("GKYL_TOK_GEO_TRACE_ALL_AFTER", 0))
              tok_geo_trace_surface_row("after_ridders", inp->ftype, -1, ip,
                it, ia, ip_delta, psi_curr, alpha_curr, theta_curr, arcL_curr,
                arc_ctx.zmin, arc_ctx.zmax, rclose, ridders_min, ridders_max,
                z_curr, res.nevals, trace_elapsed);
          }

          double sep_r_curr = 0.0;
          bool at_sep_trace = tok_half_domain_sep_rz(inp, &arc_ctx,
            theta_curr, &sep_r_curr, &z_curr);

          // These are volume-interior Gauss nodes, even when their nodal
          // indices are the first or last in the local range.  The affine
          // arc map already places them correctly; only a true block-edge
          // node may be snapped exactly to the shared ray.
          bool at_xpt_anchor = false;

          double R[4] = { 0 }, dRdZ[4] = { 0 };
          double dR[4] = { 0 }, dZ[4] = { 0 };
          int nr = gkyl_tok_geo_R_psiZ(geo, psi_curr, z_curr, 4, R, dRdZ, dR, dZ);
          double root_ref = at_sep_trace ? sep_r_curr
            : (at_xpt_anchor ? arc_ctx.xpt_anchor_r : rclose);
          r_curr = choose_closest(root_ref, R, R, nr);
          drdz_curr = choose_closest(root_ref, R, dRdZ, nr);
          dr_curr = choose_closest(root_ref, R, dR, nr);
          dz_curr = choose_closest(root_ref, R, dZ, nr);
          if (at_xpt_anchor)
            r_curr = arc_ctx.xpt_anchor_r;
          if (at_sep_trace)
            r_curr = sep_r_curr;

          if (tok_geo_same_flux(psi_curr, geo->psisep) && ip_delta==0) {
            if (z_curr == geo->efit->Zxpt[0]) {
              nr = 1;
              r_curr = geo->efit->Rxpt[0];
            }
            if (z_curr == (geo->efit->num_xpts > 1 ? geo->efit->Zxpt[1] : geo->efit->Zxpt[0])) {
              nr = 1;
              r_curr = (geo->efit->num_xpts > 1 ? geo->efit->Rxpt[1] : geo->efit->Rxpt[0]);
            }
          }

          if (nr==0) {
            printf("ip = %d, it = %d, ia = %d, ip_delta = %d\n", ip, it, ia, ip_delta);
            printf("Block Type = %d | Failed to find a root at psi = %g, Z = %1.16f\n", inp->ftype, psi_curr, z_curr);
            assert(false);
          }

            phi_curr = phi_func(alpha_curr, z_curr, &arc_ctx);
          }
          cidx[TH_IDX] = it;
          int lidx = 0;
          if (ip_delta != 0)
            lidx = 3 + 3*(ip_delta-1);
          double *mc2p_fd_n = gkyl_array_fetch(up->geo_int.mc2p_nodal_fd, gkyl_range_idx(nrange, cidx));
          double *ddtheta_n = gkyl_array_fetch(up->geo_int.ddtheta_nodal, gkyl_range_idx(nrange, cidx));
          double *ddpsi_n = gkyl_array_fetch(up->geo_int.ddpsi_nodal, gkyl_range_idx(nrange, cidx));
          double *mc2p_n = gkyl_array_fetch(up->geo_int.mc2p_nodal, gkyl_range_idx(nrange, cidx));
          double *bmag_n = gkyl_array_fetch(up->geo_int.bmag_nodal, gkyl_range_idx(nrange, cidx));
          double *curlbhat_n = gkyl_array_fetch(up->geo_int.curlbhat_nodal, gkyl_range_idx(nrange, cidx));
          double *qprofile_n = gkyl_array_fetch(qprofile_nodal, gkyl_range_idx(nrange, cidx));

          mc2p_fd_n[lidx+X_IDX] = r_curr;
          mc2p_fd_n[lidx+Y_IDX] = z_curr;
          mc2p_fd_n[lidx+Z_IDX] = phi_curr;

          if(ip_delta==0){
            if (ordered_mapping) {
              ddtheta_n[0] = ordered.dr_dtheta*dTheta_dtheta;
              ddtheta_n[1] = ordered.dz_dtheta*dTheta_dtheta;
              ddtheta_n[2] = ordered.dphi_dtheta*dTheta_dtheta;
            }
            else {
              double darc_dtheta = arc_ctx.xpt_map_valid
                ? arc_ctx.xpt_map_darc_dtheta : arc_ctx.arcL_tot/(2.0*M_PI);
              ddtheta_n[0] = sin(atan2(dr_curr, dz_curr))*darc_dtheta*dTheta_dtheta;
              ddtheta_n[1] = cos(atan2(dr_curr, dz_curr))*darc_dtheta*dTheta_dtheta;
              ddtheta_n[2] = dphidtheta_func(z_curr, &arc_ctx)*dTheta_dtheta;
            }
            ddpsi_n[0] = dPsi_dpsi;
            mc2p_n[lidx+X_IDX] = r_curr;
            mc2p_n[lidx+Y_IDX] = z_curr;
            mc2p_n[lidx+Z_IDX] = phi_curr;
            bmag_n[0] = bmag_func(r_curr, z_curr, &arc_ctx);
            qprofile_n[0] = qprofile;
            curlbhat_func(psi_curr, r_curr, z_curr, phi_curr, curlbhat_n, &arc_ctx);
          }
        }
      }
    }
  }

  // Populate other alpha indices by using axisymmetry
  for (int ia=nrange->lower[AL_IDX]+1; ia<=nrange->upper[AL_IDX]; ++ia){
    cidx[AL_IDX] = ia;
    double alpha_curr = calc_running_coord(alpha_lo, ia-nrange->lower[AL_IDX], dalpha);
    double alpha_donor = calc_running_coord(alpha_lo, 0, dalpha);
    double alpha_diff = alpha_curr -  alpha_donor;
    for (int ip=nrange->lower[PSI_IDX]; ip<=nrange->upper[PSI_IDX]; ++ip) {
      cidx[PSI_IDX] = ip;
      int ip_delta_max = 3;
      for(int ip_delta = 0; ip_delta < ip_delta_max; ip_delta++){
        double psi_curr = calc_running_coord(psi_lo, ip-nrange->lower[PSI_IDX], dpsi) + modifiers[ip_delta]*delta_psi;
        for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
          cidx[TH_IDX] = it;
            int lidx = 0;
            if (ip_delta != 0)
              lidx = 3 + 3*(ip_delta-1);

          double *mc2p_fd_n = gkyl_array_fetch(up->geo_int.mc2p_nodal_fd, gkyl_range_idx(nrange, cidx));
          double *ddtheta_n = gkyl_array_fetch(up->geo_int.ddtheta_nodal, gkyl_range_idx(nrange, cidx));
          double *ddpsi_n = gkyl_array_fetch(up->geo_int.ddpsi_nodal, gkyl_range_idx(nrange, cidx));
          double *mc2p_n = gkyl_array_fetch(up->geo_int.mc2p_nodal, gkyl_range_idx(nrange, cidx));
          double *bmag_n = gkyl_array_fetch(up->geo_int.bmag_nodal, gkyl_range_idx(nrange, cidx));
          double *curlbhat_n = gkyl_array_fetch(up->geo_int.curlbhat_nodal, gkyl_range_idx(nrange, cidx));

          int donor_cidx[3] ;
          donor_cidx[AL_IDX] = nrange->lower[AL_IDX];
          donor_cidx[PSI_IDX] = ip;
          donor_cidx[TH_IDX] = it;

          double *donor_mc2p_fd_n = gkyl_array_fetch(up->geo_int.mc2p_nodal_fd, gkyl_range_idx(nrange, donor_cidx));
          double *donor_ddtheta_n = gkyl_array_fetch(up->geo_int.ddtheta_nodal, gkyl_range_idx(nrange, donor_cidx));
          double *donor_ddpsi_n = gkyl_array_fetch(up->geo_int.ddpsi_nodal, gkyl_range_idx(nrange, donor_cidx));
          double *donor_mc2p_n = gkyl_array_fetch(up->geo_int.mc2p_nodal, gkyl_range_idx(nrange, donor_cidx));
          double *donor_bmag_n = gkyl_array_fetch(up->geo_int.bmag_nodal, gkyl_range_idx(nrange, donor_cidx));
          double *donor_curlbhat_n = gkyl_array_fetch(up->geo_int.curlbhat_nodal, gkyl_range_idx(nrange, donor_cidx));

          mc2p_fd_n[lidx+X_IDX] =  donor_mc2p_fd_n[lidx+X_IDX];
          mc2p_fd_n[lidx+Y_IDX] =  donor_mc2p_fd_n[lidx+Y_IDX];
          mc2p_fd_n[lidx+Z_IDX] =  donor_mc2p_fd_n[lidx+Z_IDX] + alpha_diff;
          if(ip_delta==0){
            ddtheta_n[0] = donor_ddtheta_n[0];
            ddtheta_n[1] = donor_ddtheta_n[1];
            ddtheta_n[2] = donor_ddtheta_n[2];
            ddpsi_n[0] = donor_ddpsi_n[0];
            mc2p_n[lidx+X_IDX] = donor_mc2p_n[lidx+X_IDX];
            mc2p_n[lidx+Y_IDX] = donor_mc2p_n[lidx+Y_IDX];
            mc2p_n[lidx+Z_IDX] = donor_mc2p_n[lidx+Z_IDX] + alpha_diff;
            bmag_n[0] = donor_bmag_n[0];
            curlbhat_func(psi_curr, mc2p_n[lidx+X_IDX], mc2p_n[lidx+Y_IDX], mc2p_n[lidx+Z_IDX], curlbhat_n, &arc_ctx);
          }
        }
      }
    }
  }

  struct gkyl_nodal_ops *n2m =  gkyl_nodal_ops_new(&inp->cbasis, &inp->cgrid, false);
  gkyl_nodal_ops_n2m(n2m, &inp->cbasis, &inp->cgrid, nrange, &up->local, 3, up->geo_int.mc2p_nodal, up->geo_int.mc2p, true);
  gkyl_nodal_ops_n2m(n2m, &inp->cbasis, &inp->cgrid, nrange, &up->local, 1, up->geo_int.bmag_nodal, up->geo_int.bmag, true);
  gkyl_nodal_ops_n2m(n2m, &inp->cbasis, &inp->cgrid, nrange, &up->local, 1, qprofile_nodal, up->geo_int.qprofile, true);
  gkyl_nodal_ops_release(n2m);
  gkyl_array_release(qprofile_nodal);

  gkyl_free(arc_memo);
  gkyl_free(arc_memo_left);
  gkyl_free(arc_memo_right);
  gkyl_free(sep_trace_r);
  gkyl_free(sep_trace_z);
  gkyl_free(sep_trace_s);
  gkyl_free(arc_ctx.ext_ladder_w);
  gkyl_free(ordered_trace_storage);
}

void gkyl_tok_geo_calc_surface(struct gk_geometry* up, int dir, struct gkyl_range *nrange, double dzc[3], 
    struct gkyl_tok_geo *geo, struct gkyl_tok_geo_grid_inp *inp, struct gkyl_position_map *position_map)
{

  geo->rleft = inp->rleft;
  geo->rright = inp->rright;

  geo->inexact_roots = inp->inexact_roots;

  geo->rmax = inp->rmax;
  geo->rmin = inp->rmin;

  enum { PSI_IDX, AL_IDX, TH_IDX }; // Arrangement of computational coordinates.
  enum { X_IDX, Y_IDX, Z_IDX }; // Arrangement of cartesian coordinates.
  
  double dtheta = inp->cgrid.dx[TH_IDX],
    dpsi = inp->cgrid.dx[PSI_IDX],
    dalpha = inp->cgrid.dx[AL_IDX];
  
  double theta_lo = up->grid.lower[TH_IDX] + (up->local.lower[TH_IDX] - up->global.lower[TH_IDX])*up->grid.dx[TH_IDX],
    psi_lo = up->grid.lower[PSI_IDX] + (up->local.lower[PSI_IDX] - up->global.lower[PSI_IDX])*up->grid.dx[PSI_IDX],
    alpha_lo = up->grid.lower[AL_IDX] + (up->local.lower[AL_IDX] - up->global.lower[AL_IDX])*up->grid.dx[AL_IDX];

  double dels[2] = {1.0/sqrt(3), 1.0-1.0/sqrt(3) };
  theta_lo += dir == 2 ? 0.0 : dels[1]*dtheta/2.0;
  psi_lo += dir == 0 ? 0.0 : dels[1]*dpsi/2.0;
  alpha_lo += dir == 1 ? 0. : dels[1]*dalpha/2.0;

    
  double dx_fact = up->basis.poly_order == 1.0/up->basis.poly_order;
  dtheta *= dx_fact; dpsi *= dx_fact; dalpha *= dx_fact;

  // Used for finite differences.
  double delta_alpha = dalpha*1e-2;
  double delta_psi = dpsi*1e-2;
  double delta_theta = dtheta*1e-2;
  dzc[0] = delta_psi;
  dzc[1] = delta_alpha;
  dzc[2] = delta_theta;
  int modifiers[5] = {0, -1, 1, -2, 2};

  double rclose = inp->rclose;
  double rright = inp->rright;
  double rleft = inp->rleft;


  int nzcells;
  if (geo->use_cubics)
    nzcells = geo->rzgrid_cubic.cells[1];
  else
    nzcells = geo->rzgrid.cells[1];
  double *arc_memo = gkyl_malloc(sizeof(double[nzcells]));
  double *arc_memo_left = gkyl_malloc(sizeof(double[nzcells]));
  double *arc_memo_right = gkyl_malloc(sizeof(double[nzcells]));
  int sep_trace_capacity = 16*nzcells+1;
  double *sep_trace_r = gkyl_malloc(sizeof(double[sep_trace_capacity]));
  double *sep_trace_z = gkyl_malloc(sizeof(double[sep_trace_capacity]));
  double *sep_trace_s = gkyl_malloc(sizeof(double[sep_trace_capacity]));
  // The marched theta ladder owns its own block: how many rungs it needs is
  // discovered while marching, not known here, so it cannot be sized with the
  // rest.  It is freed alongside this one.
  double *ordered_trace_storage =
    gkyl_malloc(sizeof(double)*8*(size_t) sep_trace_capacity);

  struct arc_length_ctx arc_ctx = {
    .geo = geo,
    .arc_memo = arc_memo,
    .arc_memo_right = arc_memo_right,
    .arc_memo_left = arc_memo_left,
    .sep_trace_r = sep_trace_r,
    .sep_trace_z = sep_trace_z,
    .sep_trace_s = sep_trace_s,
    .sep_trace_capacity = sep_trace_capacity,
    .far_trace_r = ordered_trace_storage,
    .far_trace_z = ordered_trace_storage+sep_trace_capacity,
    .far_trace_s = ordered_trace_storage+2*sep_trace_capacity,
    .trace_corr_v = ordered_trace_storage+3*sep_trace_capacity,
    .map_trace_r = ordered_trace_storage+4*sep_trace_capacity,
    .map_trace_z = ordered_trace_storage+5*sep_trace_capacity,
    .map_trace_s = ordered_trace_storage+6*sep_trace_capacity,
    .map_trace_phi = ordered_trace_storage+7*sep_trace_capacity,
    .ext_ladder_w = NULL,
    .ftype = inp->ftype,
    .zmaxis = geo->zmaxis
  };
  struct plate_ctx pctx = {
    .geo = geo
  };
  tok_init_xpt_ray_target(inp, position_map, &arc_ctx);

  int cidx[3] = { 0 };
  for(int ia=nrange->lower[AL_IDX]; ia<nrange->lower[AL_IDX]+1; ++ia){
    cidx[AL_IDX] = ia;
    double alpha_curr = dir==1 ? alpha_lo + ia*dalpha : calc_running_coord(alpha_lo, ia-nrange->lower[AL_IDX], dalpha);

    for (int ip=nrange->lower[PSI_IDX]; ip<=nrange->upper[PSI_IDX]; ++ip) {
      int ip_delta_max = 5;
      for(int ip_delta = 0; ip_delta < ip_delta_max; ip_delta++){
        if((ip == nrange->lower[PSI_IDX]) && (up->local.lower[PSI_IDX]== up->global.lower[PSI_IDX]) && dir==0){
          if(ip_delta == 1 || ip_delta == 3)
            continue; // one sided stencils at edge
        }
        else if((ip == nrange->upper[PSI_IDX]) && (up->local.upper[PSI_IDX]== up->global.upper[PSI_IDX]) && dir==0){
          if(ip_delta == 2 || ip_delta == 4)
            continue; // one sided stencils at edge
        }
        else{ // interior 
          if( ip_delta == 3 || ip_delta == 4)
            continue;
        }

        double psi_curr = dir == 0 ? psi_lo + ip*dpsi : calc_running_coord(psi_lo, ip-nrange->lower[PSI_IDX], dpsi) ;
        psi_curr += modifiers[ip_delta]*delta_psi;
        // Non-uniform psi. Finite differences are calculated in calc_metric.c
        double Psi_curr;
        position_map->maps[0](0.0, &psi_curr,  &Psi_curr,  position_map->ctxs[0]);
        double dPsi_dpsi = gkyl_position_map_slope(position_map, 0, psi_curr, delta_psi, ip, nrange);
        psi_curr = Psi_curr;
        
        double darcL, arcL_curr, arcL_lo;

        // For double null blocks this should set arc_ctx :
        // zmin, zmax, rclose, arcL_tot for all blocks. No left and right
        // For a full core case:
        // also set phi_right and arcL_right
        // For a single null case:
        // also set zmin_left and zmin_right 
        if (tok_xpt_mapping_requested(inp))
          tok_prepare_ordered_map(inp, &arc_ctx, psi_curr);
        else
          tok_find_endpoints(inp, geo, &arc_ctx, &pctx, psi_curr, alpha_curr,
            arc_memo, arc_memo_left, arc_memo_right);

        darcL = arc_ctx.arcL_tot/(up->basis.poly_order*inp->cgrid.cells[TH_IDX]) * (inp->cgrid.upper[TH_IDX] - inp->cgrid.lower[TH_IDX])/2/M_PI;
        // at the beginning of each theta loop we need to reset things
        cidx[PSI_IDX] = ip;
        arcL_curr = 0.0;
        arcL_lo = (theta_lo + M_PI)/2/M_PI*arc_ctx.arcL_tot;
        double ridders_min, ridders_max;

        for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
          arcL_curr = dir==2 ? arcL_lo + it*darcL: calc_running_coord(arcL_lo, it-nrange->lower[TH_IDX], darcL);
          double theta_curr = arcL_curr*(2*M_PI/arc_ctx.arcL_tot) - M_PI ; 

          // Calculate derivatives using finite difference for ddtheta,
          // as well as transform the computational coordiante to the non-uniform field-aligned value

          // We cannot do non-uniform alpha because we are modeling axisymmetric systems
          // Non-uniform theta
          double Theta_curr;
          position_map->maps[2](0.0, &theta_curr,  &Theta_curr,  position_map->ctxs[2]);
          double dTheta_dtheta = gkyl_position_map_slope(position_map, 2, theta_curr,\
            delta_theta, it, nrange);
          theta_curr = Theta_curr;

          struct tok_ordered_point ordered = { 0.0 };
          bool ordered_mapping = tok_ordered_map_lookup(inp, &arc_ctx,
            theta_curr, alpha_curr, &ordered);
          if (tok_xpt_mapping_requested(inp) && !ordered_mapping) {
            fprintf(stderr,
              "TOK_ORDERED_MAP lookup failed ftype=%d psi=%.17g theta=%.17g\n",
              inp->ftype, psi_curr, theta_curr);
            abort();
          }
          double r_curr = 0.0, z_curr = 0.0, phi_curr = 0.0;
          double drdz_curr = 0.0, dr_curr = 0.0, dz_curr = 0.0;
          if (ordered_mapping) {
            r_curr = ordered.r; z_curr = ordered.z; phi_curr = ordered.phi;
          }
          else {
            arcL_curr = tok_xpt_theta_to_arc(inp, &arc_ctx, theta_curr);

            tok_set_ridders(inp, &arc_ctx, psi_curr, arcL_curr, &rclose, &ridders_min, &ridders_max);

          bool trace_this_block = tok_geo_trace_enabled() && tok_geo_trace_block_selected();
          if (trace_this_block && tok_geo_trace_int_env("GKYL_TOK_GEO_TRACE_BEFORE", 1))
            tok_geo_trace_surface_row("before_ridders", inp->ftype, dir, ip,
              it, ia, ip_delta, psi_curr, alpha_curr, theta_curr, arcL_curr,
              arc_ctx.zmin, arc_ctx.zmax, rclose, ridders_min, ridders_max,
              0.0, -1, 0.0);

          clock_t trace_t0 = trace_this_block ? clock() : 0;
          struct gkyl_qr_res res = gkyl_ridders(arc_length_func, &arc_ctx,
            arc_ctx.zmin, arc_ctx.zmax, ridders_min, ridders_max,
            geo->root_param.max_iter, 1e-10);
          double trace_elapsed = trace_this_block ?
            ((double) (clock() - trace_t0))/CLOCKS_PER_SEC : 0.0;
          z_curr = res.res;
          ((struct gkyl_tok_geo *)geo)->stat.nroot_cont_calls += res.nevals;

          if (trace_this_block) {
            int neval_threshold = tok_geo_trace_int_env("GKYL_TOK_GEO_TRACE_NEVAL_THRESHOLD", 25);
            double time_threshold = tok_geo_trace_double_env("GKYL_TOK_GEO_TRACE_TIME_THRESHOLD", 0.01);
            if (res.nevals >= neval_threshold || trace_elapsed >= time_threshold ||
              tok_geo_trace_int_env("GKYL_TOK_GEO_TRACE_ALL_AFTER", 0))
              tok_geo_trace_surface_row("after_ridders", inp->ftype, dir, ip,
                it, ia, ip_delta, psi_curr, alpha_curr, theta_curr, arcL_curr,
                arc_ctx.zmin, arc_ctx.zmax, rclose, ridders_min, ridders_max,
                z_curr, res.nevals, trace_elapsed);
          }

          double sep_r_curr = 0.0;
          bool at_sep_trace = tok_half_domain_sep_rz(inp, &arc_ctx,
            theta_curr, &sep_r_curr, &z_curr);

          // Only a surface normal to theta contains true theta-boundary
          // nodes.  On radial and alpha surfaces the theta coordinates are
          // Gauss points and must not be snapped to the seam.
          double fixed_root_ref = rclose;
          bool at_fixed_edge = dir == TH_IDX &&
            tok_xpt_at_fixed_edge(inp, &arc_ctx, it, nrange,
              up->local.lower[TH_IDX] == up->global.lower[TH_IDX],
              up->local.upper[TH_IDX] == up->global.upper[TH_IDX],
              &z_curr, &fixed_root_ref);
          bool at_xpt_anchor = dir == TH_IDX && arc_ctx.xpt_anchor_valid &&
            tok_xpt_at_seam(inp, it, nrange,
              up->local.lower[TH_IDX] == up->global.lower[TH_IDX],
              up->local.upper[TH_IDX] == up->global.upper[TH_IDX]);
          if (at_xpt_anchor)
            z_curr = arc_ctx.xpt_anchor_z;

          double R[4] = { 0 }, dRdZ[4] = { 0 };
          double dR[4] = { 0 }, dZ[4] = { 0 };
          int nr = gkyl_tok_geo_R_psiZ(geo, psi_curr, z_curr, 4, R, dRdZ, dR, dZ);
          double root_ref = at_sep_trace ? sep_r_curr
            : (at_xpt_anchor ? arc_ctx.xpt_anchor_r
              : (at_fixed_edge ? fixed_root_ref : rclose));
          r_curr = choose_closest(root_ref, R, R, nr);
          drdz_curr = choose_closest(root_ref, R, dRdZ, nr);
          dr_curr = choose_closest(root_ref, R, dR, nr);
          dz_curr = choose_closest(root_ref, R, dZ, nr);
          if (at_xpt_anchor)
            r_curr = arc_ctx.xpt_anchor_r;
          if (at_sep_trace)
            r_curr = sep_r_curr;

          if (tok_geo_same_flux(psi_curr, geo->psisep) && ip_delta==0) {
            if (z_curr == geo->efit->Zxpt[0]) {
              nr = 1;
              r_curr = geo->efit->Rxpt[0];
            }
            if (z_curr == (geo->efit->num_xpts > 1 ? geo->efit->Zxpt[1] : geo->efit->Zxpt[0])) {
              nr = 1;
              r_curr = (geo->efit->num_xpts > 1 ? geo->efit->Rxpt[1] : geo->efit->Rxpt[0]);
            }
          }

          if(nr==0){
            printf("ip = %d, it = %d, ia = %d, ip_delta = %d\n", ip, it, ia, ip_delta);
            printf("Block Type = %d | Failed to find a root at psi = %g, Z = %1.16f\n", inp->ftype, psi_curr, z_curr);
            assert(false);
          }

            phi_curr = phi_func(alpha_curr, z_curr, &arc_ctx);
          }
          cidx[TH_IDX] = it;
          int lidx = 0;
          if (ip_delta != 0)
            lidx = 3 + 3*(ip_delta-1);
          double *mc2p_fd_n = gkyl_array_fetch(up->geo_surf[dir].mc2p_nodal_fd, gkyl_range_idx(nrange, cidx));
          double *ddtheta_n = gkyl_array_fetch(up->geo_surf[dir].ddtheta_nodal, gkyl_range_idx(nrange, cidx));
          double *ddpsi_n = gkyl_array_fetch(up->geo_surf[dir].ddpsi_nodal, gkyl_range_idx(nrange, cidx));
          double *bmag_n = gkyl_array_fetch(up->geo_surf[dir].bmag_nodal, gkyl_range_idx(nrange, cidx));
          double *curlbhat_n = gkyl_array_fetch(up->geo_surf[dir].curlbhat_nodal, gkyl_range_idx(nrange, cidx));
          double *deltats_n = gkyl_array_fetch(up->geo_surf[dir].deltats_nodal, gkyl_range_idx(nrange, cidx));

          mc2p_fd_n[lidx+X_IDX] = r_curr;
          mc2p_fd_n[lidx+Y_IDX] = z_curr;
          mc2p_fd_n[lidx+Z_IDX] = phi_curr;

          if(ip_delta==0){
            if (ordered_mapping) {
              ddtheta_n[0] = ordered.dr_dtheta*dTheta_dtheta;
              ddtheta_n[1] = ordered.dz_dtheta*dTheta_dtheta;
              ddtheta_n[2] = ordered.dphi_dtheta*dTheta_dtheta;
            }
            else {
              double darc_dtheta = arc_ctx.xpt_map_valid
                ? arc_ctx.xpt_map_darc_dtheta : arc_ctx.arcL_tot/(2.0*M_PI);
              ddtheta_n[0] = sin(atan2(dr_curr,dz_curr))*darc_dtheta*dTheta_dtheta;
              ddtheta_n[1] = cos(atan2(dr_curr,dz_curr))*darc_dtheta*dTheta_dtheta;
              ddtheta_n[2] = dphidtheta_func(z_curr, &arc_ctx)*dTheta_dtheta;
            }
            ddpsi_n[0] = dPsi_dpsi;
            bmag_n[0] = bmag_func(r_curr, z_curr, &arc_ctx);
            curlbhat_func(psi_curr, r_curr, z_curr, phi_curr, curlbhat_n, &arc_ctx);
            deltats_n[0] = phi_curr - alpha_curr;
          }
        }
      }
    }
  }

  // Populate other alpha indices by using axisymmetry
  for (int ia=nrange->lower[AL_IDX]+1; ia<=nrange->upper[AL_IDX]; ++ia){
    cidx[AL_IDX] = ia;
    double alpha_curr = dir==1 ? alpha_lo + ia*dalpha : calc_running_coord(alpha_lo, ia-nrange->lower[AL_IDX], dalpha);
    double alpha_donor= dir==1 ? alpha_lo + nrange->lower[AL_IDX]*dalpha : calc_running_coord(alpha_lo, 0, dalpha);
    double alpha_diff = alpha_curr -  alpha_donor;

    for (int ip=nrange->lower[PSI_IDX]; ip<=nrange->upper[PSI_IDX]; ++ip) {
      cidx[PSI_IDX] = ip;
      int ip_delta_max = 5;
      for(int ip_delta = 0; ip_delta < ip_delta_max; ip_delta++){
        if((ip == nrange->lower[PSI_IDX]) && (up->local.lower[PSI_IDX]== up->global.lower[PSI_IDX]) && dir==0){
          if(ip_delta == 1 || ip_delta == 3)
            continue; // one sided stencils at edge
        }
        else if((ip == nrange->upper[PSI_IDX]) && (up->local.upper[PSI_IDX]== up->global.upper[PSI_IDX]) && dir==0){
          if(ip_delta == 2 || ip_delta == 4)
            continue; // one sided stencils at edge
        }
        else{ // interior 
          if( ip_delta == 3 || ip_delta == 4)
            continue;
        }
        double psi_curr = dir == 0 ? psi_lo + ip*dpsi : calc_running_coord(psi_lo, ip-nrange->lower[PSI_IDX], dpsi) ;
        psi_curr += modifiers[ip_delta]*delta_psi;

        for (int it=nrange->lower[TH_IDX]; it<=nrange->upper[TH_IDX]; ++it) {
          cidx[TH_IDX] = it;
            int lidx = 0;
            if (ip_delta != 0)
              lidx = 3 + 3*(ip_delta-1);

          double *mc2p_fd_n = gkyl_array_fetch(up->geo_surf[dir].mc2p_nodal_fd, gkyl_range_idx(nrange, cidx));
          double *ddtheta_n = gkyl_array_fetch(up->geo_surf[dir].ddtheta_nodal, gkyl_range_idx(nrange, cidx));
          double *ddpsi_n = gkyl_array_fetch(up->geo_surf[dir].ddpsi_nodal, gkyl_range_idx(nrange, cidx));
          double *bmag_n = gkyl_array_fetch(up->geo_surf[dir].bmag_nodal, gkyl_range_idx(nrange, cidx));
          double *curlbhat_n = gkyl_array_fetch(up->geo_surf[dir].curlbhat_nodal, gkyl_range_idx(nrange, cidx));
          double *deltats_n = gkyl_array_fetch(up->geo_surf[dir].deltats_nodal, gkyl_range_idx(nrange, cidx));

          int donor_cidx[3] ;
          donor_cidx[AL_IDX] = nrange->lower[AL_IDX];
          donor_cidx[PSI_IDX] = ip;
          donor_cidx[TH_IDX] = it;

          double *donor_mc2p_fd_n = gkyl_array_fetch(up->geo_surf[dir].mc2p_nodal_fd, gkyl_range_idx(nrange, donor_cidx));
          double *donor_ddtheta_n = gkyl_array_fetch(up->geo_surf[dir].ddtheta_nodal, gkyl_range_idx(nrange, donor_cidx));
          double *donor_ddpsi_n = gkyl_array_fetch(up->geo_surf[dir].ddpsi_nodal, gkyl_range_idx(nrange, donor_cidx));
          double *donor_bmag_n = gkyl_array_fetch(up->geo_surf[dir].bmag_nodal, gkyl_range_idx(nrange, donor_cidx));

          mc2p_fd_n[lidx+X_IDX] =  donor_mc2p_fd_n[lidx+X_IDX];
          mc2p_fd_n[lidx+Y_IDX] =  donor_mc2p_fd_n[lidx+Y_IDX];
          mc2p_fd_n[lidx+Z_IDX] =  donor_mc2p_fd_n[lidx+Z_IDX] + alpha_diff;
          if(ip_delta==0){
            ddtheta_n[0] = donor_ddtheta_n[0];
            ddtheta_n[1] = donor_ddtheta_n[1];
            ddtheta_n[2] = donor_ddtheta_n[2];
            ddpsi_n[0] = donor_ddpsi_n[0];
            bmag_n[0] = donor_bmag_n[0];
            curlbhat_func(psi_curr, mc2p_fd_n[X_IDX], mc2p_fd_n[Y_IDX], mc2p_fd_n[Z_IDX], curlbhat_n, &arc_ctx);
            deltats_n[0] = mc2p_fd_n[Z_IDX] - alpha_curr;
          }
        }
      }
    }
  }

  gkyl_free(arc_memo);
  gkyl_free(arc_memo_left);
  gkyl_free(arc_memo_right);
  gkyl_free(sep_trace_r);
  gkyl_free(sep_trace_z);
  gkyl_free(sep_trace_s);
  gkyl_free(arc_ctx.ext_ladder_w);
  gkyl_free(ordered_trace_storage);
}


void
gkyl_tok_geo_set_extent(struct gkyl_tok_geo_grid_inp* inp, struct gkyl_tok_geo *geo, double *theta_lo, double *theta_up)
{
  tok_geo_set_extent(inp, geo, theta_lo, theta_up);
}

struct gkyl_tok_geo_stat
gkyl_tok_geo_get_stat(const struct gkyl_tok_geo *geo)
{
  return geo->stat;
}

void
gkyl_tok_geo_release(struct gkyl_tok_geo *geo)
{
  gkyl_array_release(geo->psiRZ);
  gkyl_array_release(geo->psiRZ_cubic);
  gkyl_array_release(geo->fpoldg);
  gkyl_array_release(geo->fpolprimedg);
  gkyl_array_release(geo->qdg);
  gkyl_efit_release(geo->efit);
  gkyl_free(geo);
}
