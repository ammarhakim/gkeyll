// Operator-decomposed energy budget diagnostic for the DG GR Maxwell solver.
//
// Stage 1 deliverable for debugging the ergosphere checkerboard mode in
// rt_dg_gr_maxwell_wald_kerr.lua. For a programmatically-constructed
// perturbation dQ we compute the partitioned semi-discrete growth proxy
//
//     <dQ, L_partition dQ> / <dQ, dQ>
//
// where L is the assembled DG operator and the partitions are
//
//   volume, r-interior faces, theta-interior faces,
//   r-inner BC face, r-outer BC face,
//   theta-pole-lo face, theta-pole-up face.
//
// We exploit the fact that the GR Maxwell DG operator is linear in the field
// for fixed metric: L(Q_wald + dQ) - L(Q_wald) == L(dQ). So this diagnostic
// applies L to dQ directly and never needs the Wald background, only the
// metric (lapse, shift, h_ij, det_h) projected onto the same nodal points the
// production solver uses.
//
// Norm: Euclidean modal sum over (cell, component, mode). The curved-space
// quadratic form 1/2 J_c (E_i D^i + H_i B^i) is a follow-up; that integrand
// is indefinite inside the ergosphere and needs a separate utility.

#include <acutest.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_dg_gr_maxwell_conf_flux_surf.h>
#include <gkyl_dg_gr_maxwell_conf_flux_surf_priv.h>
#include <gkyl_dg_gr_maxwell_divide_Jc.h>
#include <gkyl_dg_gr_maxwell_geom.h>
#include <gkyl_dg_gr_maxwell_kernels.h>
#include <gkyl_dg_gr_maxwell_surf_and_vol_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

// ---- Configuration ----
//
// Defaults match a cut-down version of rt_dg_gr_maxwell_wald_kerr.lua:
// same r/theta extents, same Kerr spin, ser_p1 basis. Resolution is small
// (32x32) by default so the unit test is fast; the failure-eliciting run
// uses 192x192 and is too large to iterate on inside acutest.

#define MASS_BH         1.0
#define SPIN_BH         0.95
#define R_LOWER         1.5
#define R_UPPER         5.0
#define THETA_LOWER     0.0
#define THETA_UPPER     M_PI

// Polynomial order is fixed to 1 in this file; the surface-flux kernels and
// surfx/surfy kernels are picked up by name (ser_p1 variants) below.
#define POLY_ORDER      1
#define NUM_BASIS       4    // 2D serendipity p=1: {1, xi, eta, xi*eta}
#define NUM_FIELD_COMPS 8    // JD^{r,theta,phi}, JB^{r,theta,phi}, phi, psi
#define NUM_QUAD_FACE   2    // (p+1) gauss nodes per surface in 2D p=1
#define CONF_FLUX_PER_CELL (2*NUM_FIELD_COMPS*NUM_QUAD_FACE) // = 32

// Partition tags.
typedef enum {
  PART_VOLUME = 0,
  PART_R_INTERIOR,
  PART_THETA_INTERIOR,
  PART_R_INNER_BC,
  PART_R_OUTER_BC,
  PART_THETA_POLE_LO,
  PART_THETA_POLE_UP,
  PART_COUNT,
} partition_t;

static const char *partition_name[PART_COUNT] = {
  "volume",
  "r-interior faces",
  "theta-interior faces",
  "r-inner BC",
  "r-outer BC",
  "theta-pole lo",
  "theta-pole up",
};

// ---- Bench (encapsulates state for a single (Nr, Ntheta) configuration) ----

struct budget_bench {
  bool use_lax;                          // true: LLF; false: Roe.
  int nr, ntheta;
  struct gkyl_rect_grid grid;
  struct gkyl_basis basis;
  struct gkyl_range local;
  struct gkyl_range local_ext;
  int ghost[GKYL_MAX_DIM];
  struct gkyl_dg_gr_maxwell_geom_ctx geom_ctx;
  // Geometry projected at vol+surf nodes.
  struct gkyl_surf_and_vol_node_arrays *lapse, *shift, *h_ij, *det_h;
  // Per-direction BC configuration. Default (set by bench_init) matches the
  // production input: reflective theta-pole at θ ∈ {0, π}, no outflow.
  // bench_init_full can override to test BC dependence of the spectrum.
  int theta_pole_lo[GKYL_MAX_CDIM];
  int theta_pole_up[GKYL_MAX_CDIM];
  int outflow_lo[GKYL_MAX_CDIM];
  int outflow_up[GKYL_MAX_CDIM];
  // Diagnostic flag: when set, after conf_flux_surf_advance fills conf_flux_surf
  // for all faces (using the standard theta_pole=1 path that handles the
  // coordinate singularity safely), we explicitly zero the dir=1 slice for
  // theta-pole face cells. This produces an "insulating-wall" theta-pole
  // (no flux exchange) — a non-reflective alternative to the reflective
  // dissipation-zeroed LLF that the production setup uses. Useful for
  // diagnosing whether the theta-pole reflection is what's driving the
  // unstable spectrum (vs ergoregion-interior).
  int zero_theta_pole_flux;
  // Conf-flux-surf updater (we use up->conf_flux_surf directly per face).
  struct gkyl_dg_gr_maxwell_conf_flux_surf *conf_flux_up;
  // State arrays. dQ is the perturbation we apply L to.
  struct gkyl_array *dQ_with_J;          // [JD, JB, phi, psi] modes per cell.
  struct gkyl_array *dQ_no_J;            // [D, B, phi, psi] modes per cell.
  // Staging arrays.
  struct gkyl_array *conf_flux_surf;     // populated per partition.
  struct gkyl_array *RHS;                // partition's contribution.
  struct gkyl_array *cflrate;            // per-cell scalar; required by conf_flux_surf_advance.
};

// ---- Helpers ----

static struct gkyl_array *
mkarr(long nc, long size)
{
  return gkyl_array_new(GKYL_DOUBLE, nc, size);
}

// Apply gkyl_dg_gr_maxwell_divide_Jc to the dQ_with_J state to populate dQ_no_J.
static void
recompute_dQ_no_J(struct budget_bench *bench)
{
  gkyl_dg_gr_maxwell_divide_Jc(&bench->basis, &bench->local,
    bench->det_h->nodal_arr_vol, bench->dQ_with_J, bench->dQ_no_J, false);
}

// Zero an array on the local extended range (modes set to 0).
static void
zero_array(struct gkyl_array *arr)
{
  gkyl_array_clear(arr, 0.0);
}

// Euclidean modal inner product over confLocal: sum over cells, components,
// modes of a*b. Each modal coefficient counts equally; this is uniform across
// the domain regardless of geometry. Indefiniteness inside the ergosphere is
// not visible in this norm — that is by design (the curved-space norm is the
// follow-up).
static double
euclidean_inner(const struct budget_bench *bench,
  const struct gkyl_array *a, const struct gkyl_array *b)
{
  double s = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &bench->local);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(&bench->local, iter.idx);
    const double *ad = gkyl_array_cfetch(a, lidx);
    const double *bd = gkyl_array_cfetch(b, lidx);
    for (int k = 0; k < NUM_FIELD_COMPS*NUM_BASIS; ++k) {
      s += ad[k]*bd[k];
    }
  }
  return s;
}

// ---- Volume RHS application ----
// Iterate confLocal cells, calling gr_maxwell_vol_2x_ser_p1 to accumulate
// volume contribution into bench->RHS. Caller is responsible for zeroing RHS
// beforehand.
static void
apply_volume_to_RHS(struct budget_bench *bench)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &bench->local);
  while (gkyl_range_iter_next(&iter)) {
    double xc[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&bench->grid, iter.idx, xc);
    long cidx = gkyl_range_idx(&bench->local, iter.idx);

    const double *lapse_d = gkyl_array_cfetch(bench->lapse->nodal_arr_vol, cidx);
    const double *shift_d = gkyl_array_cfetch(bench->shift->nodal_arr_vol, cidx);
    const double *h_ij_d  = gkyl_array_cfetch(bench->h_ij->nodal_arr_vol, cidx);
    const double *det_h_d = gkyl_array_cfetch(bench->det_h->nodal_arr_vol, cidx);

    const double *fields = gkyl_array_cfetch(bench->dQ_with_J, cidx);
    double *out = gkyl_array_fetch(bench->RHS, cidx);

    // The volume kernel is the production kernel; using dQ_with_J as input
    // gives the linearization L_volume(dQ) by linearity.
    gr_maxwell_vol_2x_ser_p1(xc, bench->grid.dx, lapse_d, shift_d,
      h_ij_d, det_h_d, fields, out);
  }
}

// ---- Surface RHS application ----
// Given a populated conf_flux_surf array, accumulate surfx + surfy
// contributions into bench->RHS for every cell in confLocal.
//
// Each cell reads its own lower-face flux (cidxC) and its right-neighbor's
// lower-face flux (cidxR = upper-face flux of current cell). We rely on the
// extended range so the right ghost cells contain populated data when we are
// at confLocal.upper[dir].
static void
apply_surface_to_RHS_using_flux(struct budget_bench *bench)
{
  // dir=0 (r-direction): surfx_2x_ser_p1
  {
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &bench->local);
    while (gkyl_range_iter_next(&iter)) {
      double xc[GKYL_MAX_DIM];
      gkyl_rect_grid_cell_center(&bench->grid, iter.idx, xc);
      long cidxC = gkyl_range_idx(&bench->local, iter.idx);
      int idx_r[GKYL_MAX_DIM];
      gkyl_copy_int_arr(bench->grid.ndim, iter.idx, idx_r);
      idx_r[0] += 1;
      long cidxR = gkyl_range_idx(&bench->local_ext, idx_r);

      const double *flux_l = gkyl_array_cfetch(bench->conf_flux_surf, cidxC);
      const double *flux_r = gkyl_array_cfetch(bench->conf_flux_surf, cidxR);
      double *out = gkyl_array_fetch(bench->RHS, cidxC);

      gr_maxwell_surfx_2x_ser_p1(xc, bench->grid.dx, flux_l, flux_r, out);
    }
  }
  // dir=1 (theta-direction): surfy_2x_ser_p1
  {
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &bench->local);
    while (gkyl_range_iter_next(&iter)) {
      double xc[GKYL_MAX_DIM];
      gkyl_rect_grid_cell_center(&bench->grid, iter.idx, xc);
      long cidxC = gkyl_range_idx(&bench->local, iter.idx);
      int idx_r[GKYL_MAX_DIM];
      gkyl_copy_int_arr(bench->grid.ndim, iter.idx, idx_r);
      idx_r[1] += 1;
      long cidxR = gkyl_range_idx(&bench->local_ext, idx_r);

      const double *flux_l = gkyl_array_cfetch(bench->conf_flux_surf, cidxC);
      const double *flux_r = gkyl_array_cfetch(bench->conf_flux_surf, cidxR);
      double *out = gkyl_array_fetch(bench->RHS, cidxC);

      gr_maxwell_surfy_2x_ser_p1(xc, bench->grid.dx, flux_l, flux_r, out);
    }
  }
}

// ---- Partition predicate ----
// Returns true if the face *at the lower-direction boundary of cell idx*,
// in direction dir, belongs to the requested partition. Boundary faces
// (whose owning idx lies on the lower edge) are classified as boundary;
// the upper-edge faces, owned by the right ghost cell, are classified by
// the ghost-cell idx.
static bool
face_in_partition(const struct budget_bench *bench, const int *idx, int dir,
  partition_t which)
{
  // Determine the face class first.
  bool is_inner_bc = (dir == 0) && (idx[0] == bench->local.lower[0]);
  bool is_outer_bc = (dir == 0) && (idx[0] == bench->local.upper[0] + 1);
  bool is_pole_lo  = (dir == 1) && (idx[1] == bench->local.lower[1]) && bench->theta_pole_lo[1];
  bool is_pole_up  = (dir == 1) && (idx[1] == bench->local.upper[1] + 1) && bench->theta_pole_up[1];
  bool is_r_face   = (dir == 0);
  bool is_t_face   = (dir == 1);

  switch (which) {
    case PART_R_INNER_BC:    return is_inner_bc;
    case PART_R_OUTER_BC:    return is_outer_bc;
    case PART_THETA_POLE_LO: return is_pole_lo;
    case PART_THETA_POLE_UP: return is_pole_up;
    case PART_R_INTERIOR:
      return is_r_face && !is_inner_bc && !is_outer_bc;
    case PART_THETA_INTERIOR:
      return is_t_face && !is_pole_lo && !is_pole_up;
    case PART_VOLUME:
    case PART_COUNT:
      return false;
  }
  return false;
}

// ---- Populate conf_flux_surf for a single surface partition ----
// Mirrors the dispatch loop in dg_gr_maxwell_conf_flux_surf_advance, but only
// invokes up->conf_flux_surf for faces in the requested partition. All other
// face slots in conf_flux_surf are left zero (the array is cleared first).
static void
populate_conf_flux_for_partition(struct budget_bench *bench, partition_t which)
{
  zero_array(bench->conf_flux_surf);
  if (which == PART_VOLUME || which == PART_COUNT) return;

  struct gkyl_dg_gr_maxwell_conf_flux_surf *up = bench->conf_flux_up;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &bench->local);
  while (gkyl_range_iter_next(&iter)) {
    int idx[GKYL_MAX_DIM];
    gkyl_copy_int_arr(bench->grid.ndim, iter.idx, idx);
    long cidxC = gkyl_range_idx(&bench->local, idx);

    double xcC[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&bench->grid, idx, xcC);

    const double *fc_c = gkyl_array_cfetch(bench->dQ_with_J, cidxC);
    const double *fnc_c = gkyl_array_cfetch(bench->dQ_no_J, cidxC);
    double *fl = gkyl_array_fetch(bench->conf_flux_surf, cidxC);

    // For each direction, the cell owns its lower face in that direction.
    for (int dir = 0; dir < bench->grid.ndim; ++dir) {
      // Geometry pulled from the surface arrays at this cell, this direction.
      const double *lapse_d, *shift_d, *h_ij_d, *det_h_d;
      if (dir == 0) {
        lapse_d = gkyl_array_cfetch(bench->lapse->nodal_arr_surf_x, cidxC);
        shift_d = gkyl_array_cfetch(bench->shift->nodal_arr_surf_x, cidxC);
        h_ij_d  = gkyl_array_cfetch(bench->h_ij->nodal_arr_surf_x, cidxC);
        det_h_d = gkyl_array_cfetch(bench->det_h->nodal_arr_surf_x, cidxC);
      } else {
        lapse_d = gkyl_array_cfetch(bench->lapse->nodal_arr_surf_y, cidxC);
        shift_d = gkyl_array_cfetch(bench->shift->nodal_arr_surf_y, cidxC);
        h_ij_d  = gkyl_array_cfetch(bench->h_ij->nodal_arr_surf_y, cidxC);
        det_h_d = gkyl_array_cfetch(bench->det_h->nodal_arr_surf_y, cidxC);
      }

      int idx_l[GKYL_MAX_DIM];
      gkyl_copy_int_arr(bench->grid.ndim, idx, idx_l);
      idx_l[dir] -= 1;
      long cidxL = gkyl_range_idx(&bench->local_ext, idx_l);
      const double *fc_l  = gkyl_array_cfetch(bench->dQ_with_J, cidxL);
      const double *fnc_l = gkyl_array_cfetch(bench->dQ_no_J, cidxL);

      int theta_pole = 0;
      if (idx[dir] == bench->local.lower[dir] && bench->theta_pole_lo[dir]) {
        theta_pole = 1;
      }

      if (face_in_partition(bench, idx, dir, which)) {
        up->conf_flux_surf(up, dir, xcC, bench->grid.dx, theta_pole,
          lapse_d, shift_d, h_ij_d, det_h_d, fc_l, fc_c, fnc_l, fnc_c, fl);
      }

      // Upper-edge face: owned by the right ghost.
      if (idx[dir] == bench->local.upper[dir]) {
        int idx_r[GKYL_MAX_DIM];
        gkyl_copy_int_arr(bench->grid.ndim, idx, idx_r);
        idx_r[dir] += 1;
        long cidxR = gkyl_range_idx(&bench->local_ext, idx_r);

        const double *fc_r  = gkyl_array_cfetch(bench->dQ_with_J, cidxR);
        const double *fnc_r = gkyl_array_cfetch(bench->dQ_no_J, cidxR);
        double *flR = gkyl_array_fetch(bench->conf_flux_surf, cidxR);

        // Right-edge geometry comes from the right ghost's surface arrays.
        const double *lapse_r, *shift_r, *h_ij_r, *det_h_r;
        if (dir == 0) {
          lapse_r = gkyl_array_cfetch(bench->lapse->nodal_arr_surf_x, cidxR);
          shift_r = gkyl_array_cfetch(bench->shift->nodal_arr_surf_x, cidxR);
          h_ij_r  = gkyl_array_cfetch(bench->h_ij->nodal_arr_surf_x, cidxR);
          det_h_r = gkyl_array_cfetch(bench->det_h->nodal_arr_surf_x, cidxR);
        } else {
          lapse_r = gkyl_array_cfetch(bench->lapse->nodal_arr_surf_y, cidxR);
          shift_r = gkyl_array_cfetch(bench->shift->nodal_arr_surf_y, cidxR);
          h_ij_r  = gkyl_array_cfetch(bench->h_ij->nodal_arr_surf_y, cidxR);
          det_h_r = gkyl_array_cfetch(bench->det_h->nodal_arr_surf_y, cidxR);
        }

        int theta_pole_r = 0;
        if (idx[dir] == bench->local.upper[dir] && bench->theta_pole_up[dir]) {
          theta_pole_r = 1;
        }

        double xcR[GKYL_MAX_DIM];
        gkyl_rect_grid_cell_center(&bench->grid, idx_r, xcR);

        if (face_in_partition(bench, idx_r, dir, which)) {
          up->conf_flux_surf(up, dir, xcR, bench->grid.dx, theta_pole_r,
            lapse_r, shift_r, h_ij_r, det_h_r, fc_c, fc_r, fnc_c, fnc_r, flR);
        }
      }
    }
  }
}

// ---- One partition, one growth proxy ----
// Returns <dQ, L_partition dQ>. Caller normalizes by <dQ, dQ>.
static double
partition_inner(struct budget_bench *bench, partition_t which)
{
  zero_array(bench->RHS);

  if (which == PART_VOLUME) {
    apply_volume_to_RHS(bench);
  } else {
    populate_conf_flux_for_partition(bench, which);
    apply_surface_to_RHS_using_flux(bench);
  }
  return euclidean_inner(bench, bench->dQ_with_J, bench->RHS);
}

// ---- Bench setup / release ----

// Full bench setup with BC configuration.
//   theta_pole_active: 1 -> use pole regularity (sign-flip + zero LLF jump at pole).
//                      0 -> treat θ-pole as a plain face (no special handling).
//   outflow_r:         1 -> use bcOutflow at both r boundaries (instead of LLF/Roe).
//   outflow_theta:     1 -> use bcOutflow at both θ boundaries (replaces pole regularity).
//                      Setting outflow_theta=1 implies theta_pole_active=0.
static void
bench_init_full_v2(struct budget_bench *bench, int nr, int ntheta, bool use_lax,
  int theta_pole_active, int outflow_r, int outflow_theta, bool use_curved_norm);

static void
bench_init_full(struct budget_bench *bench, int nr, int ntheta, bool use_lax,
  int theta_pole_active, int outflow_r, int outflow_theta)
{
  bench_init_full_v2(bench, nr, ntheta, use_lax,
    theta_pole_active, outflow_r, outflow_theta, /*use_curved_norm=*/false);
}

static void
bench_init_full_v2(struct budget_bench *bench, int nr, int ntheta, bool use_lax,
  int theta_pole_active, int outflow_r, int outflow_theta, bool use_curved_norm)
{
  bench->use_lax = use_lax;
  bench->nr = nr;
  bench->ntheta = ntheta;
  bench->geom_ctx.mass_bh = MASS_BH;
  bench->geom_ctx.spin_bh = SPIN_BH;

  // Grid + basis.
  double lower[2] = { R_LOWER, THETA_LOWER };
  double upper[2] = { R_UPPER, THETA_UPPER };
  int cells[2] = { nr, ntheta };
  gkyl_rect_grid_init(&bench->grid, 2, lower, upper, cells);
  gkyl_cart_modal_serendip(&bench->basis, 2, POLY_ORDER);

  // Local ranges with one ghost in each direction.
  bench->ghost[0] = 1; bench->ghost[1] = 1;
  gkyl_create_grid_ranges(&bench->grid, bench->ghost,
    &bench->local_ext, &bench->local);

  // Resolve BC flags.
  //
  // The theta-pole is a coordinate singularity: h_zz = r^2 sin^2(theta) -> 0
  // there, so det_h -> 0 and the alpha_quad kernel's J-weighted formula
  // computes alpha*sqrt(h_xx h_zz - h_xz^2)/det_h = 0/0 = NaN. The kernel's
  // theta_pole=1 path bypasses this by using a non-J-weighted formula AND
  // zeroing the dissipative jump (which makes the pole reflective).
  //
  // For our diagnostic "non-reflective theta-pole" we leave theta_pole=1
  // active in the kernel (so it doesn't NaN) but, AFTER conf_flux_surf_advance
  // runs, explicitly zero the dir=1 face flux at the theta-pole face cells.
  // That gives an insulating-wall pole (no flux through) without the
  // reflective dissipation-zeroing; effectively removes the perfect
  // reflection while keeping the pole numerically stable.
  bench->zero_theta_pole_flux = outflow_theta ? 1 : 0;
  bench->theta_pole_lo[0] = 0;
  bench->theta_pole_lo[1] = theta_pole_active ? 1 : 0;
  bench->theta_pole_up[0] = 0;
  bench->theta_pole_up[1] = theta_pole_active ? 1 : 0;
  bench->outflow_lo[0] = outflow_r ? 1 : 0;
  bench->outflow_up[0] = outflow_r ? 1 : 0;
  bench->outflow_lo[1] = 0; // see note above; we don't use outflow_flux_y at the pole
  bench->outflow_up[1] = 0;

  // Geometry projection at vol+surf nodes.
  enum gkyl_triad_preset_geom_type t = GKYL_TRIAD_GR_KERR_SCHILD_RTHETA;
  bool use_gpu = false;

  struct gkyl_dg_gr_maxwell_surf_and_vol_nodes *p;
  p = gkyl_dg_gr_maxwell_surf_and_vol_nodes_new(&bench->grid, &bench->basis,
    1, POLY_ORDER, gkyl_dg_gr_maxwell_preset_lapse(t), &bench->geom_ctx);
  bench->lapse = gkyl_surf_and_vol_node_arrays_new(p, bench->local_ext.volume, use_gpu);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance(p, 0.0, &bench->local_ext, bench->lapse);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_release(p);

  p = gkyl_dg_gr_maxwell_surf_and_vol_nodes_new(&bench->grid, &bench->basis,
    3, POLY_ORDER, gkyl_dg_gr_maxwell_preset_shift(t), &bench->geom_ctx);
  bench->shift = gkyl_surf_and_vol_node_arrays_new(p, bench->local_ext.volume, use_gpu);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance(p, 0.0, &bench->local_ext, bench->shift);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_release(p);

  p = gkyl_dg_gr_maxwell_surf_and_vol_nodes_new(&bench->grid, &bench->basis,
    6, POLY_ORDER, gkyl_dg_gr_maxwell_preset_h_ij(t), &bench->geom_ctx);
  bench->h_ij = gkyl_surf_and_vol_node_arrays_new(p, bench->local_ext.volume, use_gpu);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance(p, 0.0, &bench->local_ext, bench->h_ij);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_release(p);

  p = gkyl_dg_gr_maxwell_surf_and_vol_nodes_new(&bench->grid, &bench->basis,
    1, POLY_ORDER, gkyl_dg_gr_maxwell_preset_det_h(t), &bench->geom_ctx);
  bench->det_h = gkyl_surf_and_vol_node_arrays_new(p, bench->local_ext.volume, use_gpu);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_advance(p, 0.0, &bench->local_ext, bench->det_h);
  gkyl_dg_gr_maxwell_surf_and_vol_nodes_release(p);

  // Conf-flux-surf updater.
  struct gkyl_dg_gr_maxwell_conf_flux_surf_inp inp = {
    .conf_basis = &bench->basis,
    .conf_grid = &bench->grid,
    .field_id = GKYL_FIELD_GR_D_B,
    .theta_pole_lo = bench->theta_pole_lo,
    .theta_pole_up = bench->theta_pole_up,
    .outflow_lo = bench->outflow_lo,
    .outflow_up = bench->outflow_up,
    .use_lax = use_lax,
    .use_curved_norm = use_curved_norm,
    .use_gpu = use_gpu,
  };
  bench->conf_flux_up = gkyl_dg_gr_maxwell_conf_flux_surf_inew(&inp);

  // State arrays.
  bench->dQ_with_J = mkarr(NUM_FIELD_COMPS*NUM_BASIS, bench->local_ext.volume);
  bench->dQ_no_J   = mkarr(NUM_FIELD_COMPS*NUM_BASIS, bench->local_ext.volume);
  bench->conf_flux_surf = mkarr(CONF_FLUX_PER_CELL, bench->local_ext.volume);
  bench->RHS       = mkarr(NUM_FIELD_COMPS*NUM_BASIS, bench->local_ext.volume);
  bench->cflrate   = mkarr(1, bench->local_ext.volume);

  zero_array(bench->dQ_with_J);
  zero_array(bench->dQ_no_J);
  zero_array(bench->conf_flux_surf);
  zero_array(bench->RHS);
  zero_array(bench->cflrate);
}

// Production-default BC config (theta-pole reflective, no outflow). Existing
// callers retain their behavior unchanged.
static void
bench_init(struct budget_bench *bench, int nr, int ntheta, bool use_lax)
{
  bench_init_full(bench, nr, ntheta, use_lax,
    /*theta_pole_active=*/1, /*outflow_r=*/0, /*outflow_theta=*/0);
}

static void
bench_release(struct budget_bench *bench)
{
  gkyl_array_release(bench->dQ_with_J);
  gkyl_array_release(bench->dQ_no_J);
  gkyl_array_release(bench->conf_flux_surf);
  gkyl_array_release(bench->RHS);
  gkyl_array_release(bench->cflrate);
  gkyl_dg_gr_maxwell_conf_flux_surf_release(bench->conf_flux_up);
  gkyl_surf_and_vol_node_arrays_release(bench->lapse);
  gkyl_surf_and_vol_node_arrays_release(bench->shift);
  gkyl_surf_and_vol_node_arrays_release(bench->h_ij);
  gkyl_surf_and_vol_node_arrays_release(bench->det_h);
}

// ---- Seeded perturbations (programmatic, cell-average mode only) ----
//
// Each seed writes to dQ_with_J, then divides by J_c to get dQ_no_J. The
// "cell-average mode" means we set basis coefficient 0 (the constant mode)
// for the chosen components and zero the slope/cross modes. This is the
// most aggressive checkerboard pattern: per-cell jump of ±amp at every
// face, zero gradient inside cells.

// Component indices (with-J): 0=JD^r, 1=JD^theta, 2=JD^phi, 3=JB^r, 4=JB^theta, 5=JB^phi.
// 6,7 = phi, psi cleaning (we leave these at zero).

static void
seed_zero_dQ(struct budget_bench *bench)
{
  zero_array(bench->dQ_with_J);
  zero_array(bench->dQ_no_J);
}

// Set dQ_with_J's component-`comp` cell-average mode to amp*sign(idx) per cell,
// for all cells in confLocal. Other modes and components untouched.
typedef double (*sign_fn_t)(const int *idx, void *ctx);

struct gate_ctx {
  int lower_r, lower_theta;
  double r_lower, dr;       // for r-position gating
  bool gate_r_lt; double r_max; // include only r < r_max
  bool gate_r_gt; double r_min; // include only r > r_min
};

static double sign_alt_r(const int *idx, void *ctx_) {
  struct gate_ctx *ctx = ctx_;
  // Map (idx[0] - lower_r) -> ± alternating sign; gate by r-window if requested.
  double r_center = ctx->r_lower + (idx[0] - ctx->lower_r + 0.5)*ctx->dr;
  if (ctx->gate_r_lt && r_center >= ctx->r_max) return 0.0;
  if (ctx->gate_r_gt && r_center <= ctx->r_min) return 0.0;
  return ((idx[0] - ctx->lower_r) & 1) ? -1.0 : 1.0;
}
static double sign_alt_theta(const int *idx, void *ctx_) {
  struct gate_ctx *ctx = ctx_;
  return ((idx[1] - ctx->lower_theta) & 1) ? -1.0 : 1.0;
}
static double sign_alt_rt(const int *idx, void *ctx_) {
  struct gate_ctx *ctx = ctx_;
  double r_center = ctx->r_lower + (idx[0] - ctx->lower_r + 0.5)*ctx->dr;
  if (ctx->gate_r_lt && r_center >= ctx->r_max) return 0.0;
  if (ctx->gate_r_gt && r_center <= ctx->r_min) return 0.0;
  int rt = (idx[0] - ctx->lower_r) + (idx[1] - ctx->lower_theta);
  return (rt & 1) ? -1.0 : 1.0;
}

// Seed: cell-average mode on a single component, with sign pattern from sign_fn.
// A pure cell-average seed has zero spatial gradient *inside* each cell, so the
// volume term vanishes by construction; all contributions come from face jumps.
static void
seed_dQ_pattern(struct budget_bench *bench, int comp, double amp,
  sign_fn_t sf, void *sf_ctx)
{
  zero_array(bench->dQ_with_J);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &bench->local);
  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(&bench->local, iter.idx);
    double *q = gkyl_array_fetch(bench->dQ_with_J, cidx);
    double s = sf(iter.idx, sf_ctx);
    q[comp*NUM_BASIS + 0] = amp*s;
  }
  recompute_dQ_no_J(bench);
}

// Seed: per-cell slope mode (xi, eta, or xi*eta) on a single component.
// `slope`: 0 = xi (r-slope, basis idx 1), 1 = eta (theta-slope, basis idx 2),
//          2 = xi*eta (cross slope, basis idx 3).
// This excites the volume term directly because intra-cell gradients are
// nonzero. Combined with a sign pattern (sf) it produces a slope-checkerboard.
static void
seed_dQ_slope_pattern(struct budget_bench *bench, int comp, int slope,
  double amp, sign_fn_t sf, void *sf_ctx)
{
  zero_array(bench->dQ_with_J);
  int basis_idx = slope + 1;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &bench->local);
  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(&bench->local, iter.idx);
    double *q = gkyl_array_fetch(bench->dQ_with_J, cidx);
    double s = sf(iter.idx, sf_ctx);
    q[comp*NUM_BASIS + basis_idx] = amp*s;
  }
  recompute_dQ_no_J(bench);
}

// ---- Driver: print one budget table for the seeded dQ ----
static void
print_partition_table(struct budget_bench *bench, const char *label)
{
  double dQ_norm2 = euclidean_inner(bench, bench->dQ_with_J, bench->dQ_with_J);
  if (dQ_norm2 == 0.0) {
    printf("[%s] dQ is identically zero; skipping.\n", label);
    return;
  }
  printf("\n=== %s   (||dQ||^2 = %1.6e, useLax=%d, %dx%d) ===\n",
    label, dQ_norm2, (int)bench->use_lax, bench->nr, bench->ntheta);
  printf("  %-24s %18s %18s\n", "partition", "<dQ,L_p dQ>", "growth_proxy");
  double total_inner = 0.0;
  for (int p = 0; p < PART_COUNT; ++p) {
    double inner = partition_inner(bench, (partition_t)p);
    total_inner += inner;
    printf("  %-24s %18.6e %18.6e\n",
      partition_name[p], inner, inner/dQ_norm2);
  }
  printf("  %-24s %18.6e %18.6e\n",
    "(sum of partitions)", total_inner, total_inner/dQ_norm2);
}

// ---- Tests ----

// Smoke test: bench builds and tears down without crash; volume kernel runs
// on a zero perturbation and returns zero RHS.
static void
test_bench_smoke(void)
{
  struct budget_bench bench;
  bench_init(&bench, 8, 8, true);

  zero_array(bench.dQ_with_J);
  zero_array(bench.dQ_no_J);
  zero_array(bench.RHS);

  apply_volume_to_RHS(&bench);
  double s = euclidean_inner(&bench, bench.dQ_with_J, bench.RHS);
  TEST_CHECK(s == 0.0);

  // All-partition sanity: still zero.
  for (int p = 0; p < PART_COUNT; ++p) {
    double v = partition_inner(&bench, (partition_t)p);
    TEST_CHECK(v == 0.0);
  }

  bench_release(&bench);
}

// Seeded sweep: alternating-r checkerboard, alternating-theta, 2D checkerboard,
// each restricted to JB^phi-only and to all components, gated to r<2 / r>2.
// Prints the partition tables. The most positive growth_proxy across all
// seeds is a candidate for the dominant unstable mode structure.
static void
test_seeded_sweep(void)
{
  // 192x192 matches the failing run in rt_dg_gr_maxwell_wald_kerr.lua. The user
  // reports lower resolutions inhibit the mode, so we run the diagnostic at
  // production resolution. (~50 MB of arrays; runs in seconds.)
  struct budget_bench bench;
  bench_init(&bench, 192, 192, /*use_lax=*/true);

  struct gate_ctx ctx = {
    .lower_r = bench.local.lower[0],
    .lower_theta = bench.local.lower[1],
    .r_lower = R_LOWER,
    .dr = (R_UPPER - R_LOWER) / bench.nr,
    .gate_r_lt = false, .r_max = 0.0,
    .gate_r_gt = false, .r_min = 0.0,
  };
  const double amp = 1.0;

  // Seed names follow: (sign-pattern) x (component-set) x (gate).
  // Component 5 = JB^phi, component 0..5 = all six EM components.
  struct gate_ctx ctx_lt2 = ctx; ctx_lt2.gate_r_lt = true; ctx_lt2.r_max = 2.0;
  struct gate_ctx ctx_gt2 = ctx; ctx_gt2.gate_r_gt = true; ctx_gt2.r_min = 2.0;

  // 1. Alternating-r checkerboard, JB^phi only.
  seed_dQ_pattern(&bench, /*comp=*/5, amp, sign_alt_r, &ctx);
  print_partition_table(&bench, "alt-r, JB^phi, all-r");
  // 2. Alternating-theta, JB^phi only.
  seed_dQ_pattern(&bench, /*comp=*/5, amp, sign_alt_theta, &ctx);
  print_partition_table(&bench, "alt-theta, JB^phi, all-r");
  // 3. Alternating-r+theta (2D checkerboard), JB^phi only.
  seed_dQ_pattern(&bench, /*comp=*/5, amp, sign_alt_rt, &ctx);
  print_partition_table(&bench, "alt-r+theta, JB^phi, all-r");
  // 4. 2D checkerboard, JB^phi, gated r<2 (inside ergosphere).
  seed_dQ_pattern(&bench, /*comp=*/5, amp, sign_alt_rt, &ctx_lt2);
  print_partition_table(&bench, "alt-r+theta, JB^phi, r<2");
  // 5. 2D checkerboard, JB^phi, gated r>2 (outside ergosphere).
  seed_dQ_pattern(&bench, /*comp=*/5, amp, sign_alt_rt, &ctx_gt2);
  print_partition_table(&bench, "alt-r+theta, JB^phi, r>2");
  // 6. 2D checkerboard, JD^r only, gated r<2 (control).
  seed_dQ_pattern(&bench, /*comp=*/0, amp, sign_alt_rt, &ctx_lt2);
  print_partition_table(&bench, "alt-r+theta, JD^r, r<2");

  // Slope-mode seeds: these excite the volume term. xi-slope (r-derivative
  // slope), eta-slope (theta-derivative slope), and the xi*eta cross slope.
  // Sign-alternating in r+theta puts a "slope-checkerboard" on top of the
  // cell-average pattern.
  //
  // The xi-slope sweep is run inside, outside, and across the whole domain
  // so that the localization of the (positive) volume contribution to the
  // ergoregion can be checked directly.
  // 7a. xi-slope, JB^phi, r<2 (inside ergosphere).
  seed_dQ_slope_pattern(&bench, /*comp=*/5, /*slope=xi*/0, amp, sign_alt_rt, &ctx_lt2);
  print_partition_table(&bench, "xi-slope, JB^phi, r<2 (in ergo)");
  // 7b. xi-slope, JB^phi, r>2 (outside ergosphere; equatorial radius of ergosphere).
  seed_dQ_slope_pattern(&bench, /*comp=*/5, /*slope=xi*/0, amp, sign_alt_rt, &ctx_gt2);
  print_partition_table(&bench, "xi-slope, JB^phi, r>2 (out ergo)");
  // 7c. xi-slope, JB^phi, all-r (control: tests whether the volume sign is
  // a property of the kernel itself or a property of the ergoregion).
  seed_dQ_slope_pattern(&bench, /*comp=*/5, /*slope=xi*/0, amp, sign_alt_rt, &ctx);
  print_partition_table(&bench, "xi-slope, JB^phi, all-r");
  // 8. eta-slope checkerboard, JB^phi, r<2.
  seed_dQ_slope_pattern(&bench, /*comp=*/5, /*slope=eta*/1, amp, sign_alt_rt, &ctx_lt2);
  print_partition_table(&bench, "eta-slope, JB^phi, r<2");
  // 9. xi*eta cross-slope checkerboard, JB^phi, r<2.
  seed_dQ_slope_pattern(&bench, /*comp=*/5, /*slope=xi*eta*/2, amp, sign_alt_rt, &ctx_lt2);
  print_partition_table(&bench, "xi*eta-slope, JB^phi, r<2");

  // Repeat the most diagnostic seeds with Roe flux (Stage 0.1 pairing).
  // If the inner-BC dissipation differs substantially under Roe, the LLF
  // dissipation structure at the boundary face is implicated.
  bench_release(&bench);
  bench_init(&bench, 192, 192, /*use_lax=*/false);
  ctx.lower_r = bench.local.lower[0];
  ctx.lower_theta = bench.local.lower[1];
  ctx.dr = (R_UPPER - R_LOWER) / bench.nr;
  ctx_lt2 = ctx; ctx_lt2.gate_r_lt = true; ctx_lt2.r_max = 2.0;
  seed_dQ_pattern(&bench, /*comp=*/5, amp, sign_alt_rt, &ctx_lt2);
  print_partition_table(&bench, "ROE alt-r+theta, JB^phi, r<2");
  seed_dQ_pattern(&bench, /*comp=*/5, amp, sign_alt_rt, &ctx);
  print_partition_table(&bench, "ROE alt-r+theta, JB^phi, all-r");
  seed_dQ_slope_pattern(&bench, /*comp=*/5, /*slope=xi*/0, amp, sign_alt_rt, &ctx_lt2);
  print_partition_table(&bench, "ROE xi-slope, JB^phi, r<2");

  bench_release(&bench);
}

// ---- Power iteration: find the dominant eigenmode of the assembled DG operator ----
//
// Uses SSP-RK3 stepping (the same scheme the production app uses) on a
// renormalized perturbation. For a stable scheme on stable physical modes
// and unstable on a few unstable modes, the renormalized iterate converges to
// the eigenmode with the largest growth rate. After convergence,
//   growth_per_step ~ exp(lambda * dt),
//   lambda_Rayleigh = <dQ, L dQ> / <dQ, dQ>
// agree with the largest unstable eigenvalue.

// Apply theta-pole regularity BC to the ghost rows of `arr`. For r-theta
// Kerr-Schild axisymmetric runs the spherical pole demands per-component
// sign-flips of the skin field as it is reflected into the ghost: components
// 0/2/3/5 (Dr, Dphi, Br, Bphi) flip on odd modes; components 1/4 (Dtheta,
// Btheta) flip on even modes. Components 6/7 (cleaning) are left unchanged.
//
// Without this, the lax flux at the pole — which zeros the dissipative jump
// when theta_pole=1 — leaves the perturbation's pole face un-dissipated and
// admits a fast spurious eigenmode that swamps the physical instability the
// production simulation experiences.
//
// Skipped when bench->theta_pole_lo[1]/up[1] are 0 (e.g., when outflow_theta
// is configured); ghost rows are then left at zero, which is consistent with
// the outflow flux kernel never reading them.
static void
apply_theta_pole_bc(struct budget_bench *bench, struct gkyl_array *arr)
{
  if (!bench->theta_pole_lo[1] && !bench->theta_pole_up[1]) return;
  int dir = 1;
  // Lower theta-pole: ghost row at idx[1] = local.lower[1] - 1.
  for (int i = bench->local.lower[0]; i <= bench->local.upper[0]; ++i) {
    int idx_skin[2]  = { i, bench->local.lower[1]     };
    int idx_ghost[2] = { i, bench->local.lower[1] - 1 };
    long lidx_s = gkyl_range_idx(&bench->local_ext, idx_skin);
    long lidx_g = gkyl_range_idx(&bench->local_ext, idx_ghost);
    const double *s = gkyl_array_cfetch(arr, lidx_s);
    double       *g = gkyl_array_fetch (arr, lidx_g);
    bench->basis.flip_odd_sign (dir, &s[NUM_BASIS*0], &g[NUM_BASIS*0]); // Dr
    bench->basis.flip_even_sign(dir, &s[NUM_BASIS*1], &g[NUM_BASIS*1]); // Dtheta
    bench->basis.flip_odd_sign (dir, &s[NUM_BASIS*2], &g[NUM_BASIS*2]); // Dphi
    bench->basis.flip_odd_sign (dir, &s[NUM_BASIS*3], &g[NUM_BASIS*3]); // Br
    bench->basis.flip_even_sign(dir, &s[NUM_BASIS*4], &g[NUM_BASIS*4]); // Btheta
    bench->basis.flip_odd_sign (dir, &s[NUM_BASIS*5], &g[NUM_BASIS*5]); // Bphi
  }
  // Upper theta-pole: ghost row at idx[1] = local.upper[1] + 1.
  for (int i = bench->local.lower[0]; i <= bench->local.upper[0]; ++i) {
    int idx_skin[2]  = { i, bench->local.upper[1]     };
    int idx_ghost[2] = { i, bench->local.upper[1] + 1 };
    long lidx_s = gkyl_range_idx(&bench->local_ext, idx_skin);
    long lidx_g = gkyl_range_idx(&bench->local_ext, idx_ghost);
    const double *s = gkyl_array_cfetch(arr, lidx_s);
    double       *g = gkyl_array_fetch (arr, lidx_g);
    bench->basis.flip_odd_sign (dir, &s[NUM_BASIS*0], &g[NUM_BASIS*0]);
    bench->basis.flip_even_sign(dir, &s[NUM_BASIS*1], &g[NUM_BASIS*1]);
    bench->basis.flip_odd_sign (dir, &s[NUM_BASIS*2], &g[NUM_BASIS*2]);
    bench->basis.flip_odd_sign (dir, &s[NUM_BASIS*3], &g[NUM_BASIS*3]);
    bench->basis.flip_even_sign(dir, &s[NUM_BASIS*4], &g[NUM_BASIS*4]);
    bench->basis.flip_odd_sign (dir, &s[NUM_BASIS*5], &g[NUM_BASIS*5]);
  }
}

// Apply the full DG operator L to `in` (with-J representation), writing the
// result into `out`. Side-effects: bench->dQ_with_J / dQ_no_J / RHS /
// conf_flux_surf / cflrate are overwritten in the process.
static void
apply_full_L(struct budget_bench *bench, const struct gkyl_array *in,
  struct gkyl_array *out)
{
  // Stage input into bench's with-J / no-J arrays so the kernels see it.
  // Radial BCs use fixed-function around fixed Wald, so radial ghosts stay
  // zero. Theta-pole BCs are regularity sign-flips applied below.
  gkyl_array_clear(bench->dQ_with_J, 0.0);
  gkyl_array_copy_range(bench->dQ_with_J, in, &bench->local);
  apply_theta_pole_bc(bench, bench->dQ_with_J);
  recompute_dQ_no_J(bench);
  apply_theta_pole_bc(bench, bench->dQ_no_J);

  // Volume.
  zero_array(bench->RHS);
  apply_volume_to_RHS(bench);

  // Surface (all faces in one pass via the production updater).
  zero_array(bench->conf_flux_surf);
  gkyl_dg_gr_maxwell_conf_flux_surf_advance(bench->conf_flux_up,
    &bench->local, &bench->local_ext,
    bench->lapse, bench->shift, bench->h_ij, bench->det_h,
    bench->dQ_with_J, bench->dQ_no_J,
    bench->cflrate, bench->conf_flux_surf);

  // Diagnostic: convert reflective theta-pole to insulating-wall theta-pole
  // by zeroing the dir=1 face flux at the two theta-pole face cells. The
  // owning cells are: lower-pole flux at idx[1]=local.lower[1], upper-pole
  // flux at the right ghost idx[1]=local.upper[1]+1. dir=1 slice in
  // conf_flux_surf is offsets 16..31 (8 components * 2 quad nodes).
  if (bench->zero_theta_pole_flux) {
    int idx[GKYL_MAX_DIM];
    int ntheta_targets[2] = { bench->local.lower[1], bench->local.upper[1] + 1 };
    for (int t = 0; t < 2; ++t) {
      idx[1] = ntheta_targets[t];
      for (int i = bench->local.lower[0]; i <= bench->local.upper[0]; ++i) {
        idx[0] = i;
        long lidx = gkyl_range_idx(&bench->local_ext, idx);
        double *cf = gkyl_array_fetch(bench->conf_flux_surf, lidx);
        for (int k = CONF_FLUX_PER_CELL/2; k < CONF_FLUX_PER_CELL; ++k) cf[k] = 0.0;
      }
    }
  }

  apply_surface_to_RHS_using_flux(bench);

  // Hand back through `out`.
  gkyl_array_copy_range(out, bench->RHS, &bench->local);
}

// One SSP-RK3 step on `dQ`. `dQ_t` and `L_dQ` are scratch arrays. dQ is
// overwritten with dQ_{n+1}.
static void
ssp_rk3_step(struct budget_bench *bench, double dt,
  struct gkyl_array *dQ, struct gkyl_array *dQ_t, struct gkyl_array *L_dQ)
{
  // Stage 1: dQ_t = dQ + dt * L(dQ).
  apply_full_L(bench, dQ, L_dQ);
  gkyl_array_copy_range(dQ_t, dQ, &bench->local);
  gkyl_array_accumulate_range(dQ_t, dt, L_dQ, &bench->local);

  // Stage 2: dQ_t = (3/4) dQ + (1/4) dQ_t + (1/4) dt * L(dQ_t).
  apply_full_L(bench, dQ_t, L_dQ);
  gkyl_array_scale_range(dQ_t, 1.0/4.0, &bench->local);
  gkyl_array_accumulate_range(dQ_t, 3.0/4.0, dQ, &bench->local);
  gkyl_array_accumulate_range(dQ_t, dt/4.0, L_dQ, &bench->local);

  // Stage 3: dQ <- (1/3) dQ + (2/3) dQ_t + (2/3) dt * L(dQ_t).
  apply_full_L(bench, dQ_t, L_dQ);
  gkyl_array_scale_range(dQ, 1.0/3.0, &bench->local);
  gkyl_array_accumulate_range(dQ, 2.0/3.0, dQ_t, &bench->local);
  gkyl_array_accumulate_range(dQ, 2.0*dt/3.0, L_dQ, &bench->local);
}

// Run `num_steps` of SSP-RK3 power iteration with renormalization at each
// step. The starting iterate is bench->dQ_with_J (caller-set). On exit,
// bench->dQ_with_J holds the converged eigenmode (renormalized to ||dQ||^2 = 1
// in the Euclidean modal norm). If `out_fname` is non-NULL, the eigenmode is
// written to that file as a .gkyl frame for downstream visualization. If
// `out_lambda_eff` is non-NULL, the final lambda_eff is written there.
static void
power_iterate(struct budget_bench *bench, int num_steps, double dt,
  int print_every, const char *out_fname, double *out_lambda_eff)
{
  struct gkyl_array *dQ_n = mkarr(NUM_FIELD_COMPS*NUM_BASIS, bench->local_ext.volume);
  struct gkyl_array *dQ_t = mkarr(NUM_FIELD_COMPS*NUM_BASIS, bench->local_ext.volume);
  struct gkyl_array *L_dQ = mkarr(NUM_FIELD_COMPS*NUM_BASIS, bench->local_ext.volume);

  // Take seed from bench, normalize.
  gkyl_array_copy_range(dQ_n, bench->dQ_with_J, &bench->local);
  gkyl_array_copy_range(bench->dQ_with_J, dQ_n, &bench->local);
  recompute_dQ_no_J(bench);
  double seed_norm2 = euclidean_inner(bench, bench->dQ_with_J, bench->dQ_with_J);
  TEST_CHECK(seed_norm2 > 0.0);
  gkyl_array_scale_range(dQ_n, 1.0/sqrt(seed_norm2), &bench->local);

  printf("\n=== Power iteration: %d steps, dt=%g, useLax=%d, %dx%d ===\n",
    num_steps, dt, (int)bench->use_lax, bench->nr, bench->ntheta);
  printf("  step    growth_per_step    lambda_eff    lambda_Rayleigh\n");

  double lambda_eff_last = 0.0;
  for (int step = 0; step < num_steps; ++step) {
    // Rayleigh quotient pre-step (dQ_n is normalized to ||dQ_n||^2 = 1).
    apply_full_L(bench, dQ_n, L_dQ);
    gkyl_array_copy_range(bench->dQ_with_J, dQ_n, &bench->local);
    double rayleigh = euclidean_inner(bench, bench->dQ_with_J, L_dQ);

    // SSP-RK3 step on dQ_n in place.
    ssp_rk3_step(bench, dt, dQ_n, dQ_t, L_dQ);

    // Norm after one full step.
    gkyl_array_copy_range(bench->dQ_with_J, dQ_n, &bench->local);
    double norm2 = euclidean_inner(bench, bench->dQ_with_J, bench->dQ_with_J);
    double growth = sqrt(norm2);
    double lambda_eff = (growth > 0.0) ? log(growth)/dt : 0.0;
    if (step % print_every == 0 || step == num_steps - 1) {
      printf("  %4d   %1.6e   %1.6e   %1.6e\n",
        step, growth, lambda_eff, rayleigh);
    }
    lambda_eff_last = lambda_eff;

    // Renormalize.
    if (norm2 > 0.0) {
      gkyl_array_scale_range(dQ_n, 1.0/growth, &bench->local);
    }
  }

  // Stash converged mode.
  gkyl_array_copy_range(bench->dQ_with_J, dQ_n, &bench->local);
  recompute_dQ_no_J(bench);

  if (out_lambda_eff) *out_lambda_eff = lambda_eff_last;

  printf("Power iteration done. Final lambda_eff = %1.6e (compare to user's reported %1.6e).\n",
    lambda_eff_last, 11.8);

  if (out_fname) {
    enum gkyl_array_rio_status st = gkyl_grid_sub_array_write(&bench->grid,
      &bench->local, 0, bench->dQ_with_J, out_fname);
    if (st == GKYL_ARRAY_RIO_SUCCESS) {
      printf("Eigenmode (with-J representation) saved to %s\n", out_fname);
    } else {
      printf("WARNING: failed to write eigenmode file %s (status %d)\n",
        out_fname, (int)st);
    }
  }

  gkyl_array_release(dQ_n);
  gkyl_array_release(dQ_t);
  gkyl_array_release(L_dQ);
}

// ---- Diagnostics on the converged eigenmode ----

// Print the L2 norm (in Euclidean modal sum) of bench->dQ_with_J as a
// histogram in r, summing over theta. Useful to see whether the mode is
// localized at the inner boundary, theta-pole, or elsewhere.
static void
print_radial_localization(const struct budget_bench *bench, int nbins,
  const char *label)
{
  // nbins partitions of [R_LOWER, R_UPPER]. Each bin gets the L2 norm of all
  // cells whose r-index falls in that bin.
  if (nbins <= 0) return;
  double bin_norm[64] = { 0.0 }; // up to 64 bins for the diagnostic.
  if (nbins > 64) nbins = 64;
  int bin_for_idx[2048] = { 0 }; // up to 2048 r-cells.
  int Nr = bench->nr;
  for (int i = 0; i < Nr && i < 2048; ++i) {
    int b = (i*nbins) / Nr;
    if (b >= nbins) b = nbins - 1;
    bin_for_idx[i] = b;
  }
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &bench->local);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(&bench->local, iter.idx);
    const double *q = gkyl_array_cfetch(bench->dQ_with_J, lidx);
    double s = 0.0;
    for (int k = 0; k < NUM_FIELD_COMPS*NUM_BASIS; ++k) s += q[k]*q[k];
    int rb = bin_for_idx[iter.idx[0] - bench->local.lower[0]];
    bin_norm[rb] += s;
  }
  printf("%s: radial L2-norm histogram (Nr=%d, nbins=%d):\n", label, Nr, nbins);
  for (int b = 0; b < nbins; ++b) {
    double r_lo = R_LOWER + (R_UPPER - R_LOWER) * b / nbins;
    double r_hi = R_LOWER + (R_UPPER - R_LOWER) * (b + 1) / nbins;
    printf("  r in [%.3f, %.3f]: ||dQ||² = %1.6e\n", r_lo, r_hi, bin_norm[b]);
  }
}

// Print the L2 norm by component: which component of the eigenmode dominates?
static void
print_component_breakdown(const struct budget_bench *bench, const char *label)
{
  double comp_norm[NUM_FIELD_COMPS] = { 0.0 };
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &bench->local);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(&bench->local, iter.idx);
    const double *q = gkyl_array_cfetch(bench->dQ_with_J, lidx);
    for (int c = 0; c < NUM_FIELD_COMPS; ++c) {
      for (int m = 0; m < NUM_BASIS; ++m) {
        double v = q[c*NUM_BASIS + m];
        comp_norm[c] += v*v;
      }
    }
  }
  static const char *cn[NUM_FIELD_COMPS] = {
    "JD^r", "JD^theta", "JD^phi", "JB^r", "JB^theta", "JB^phi", "phi_clean", "psi_clean"
  };
  printf("%s: component breakdown of ||dQ||²:\n", label);
  double total = 0.0;
  for (int c = 0; c < NUM_FIELD_COMPS; ++c) total += comp_norm[c];
  for (int c = 0; c < NUM_FIELD_COMPS; ++c) {
    printf("  %-10s : %1.6e (%5.2f%%)\n", cn[c], comp_norm[c],
      100.0*comp_norm[c]/(total > 0 ? total : 1.0));
  }
}

// Test: power iteration on the LLF operator at production resolution. Seeds
// from a 2D-checkerboard JB^phi pattern, runs SSP-RK3 with renormalization
// for `num_steps` steps, prints the convergence trace, then prints the
// partition table for the converged eigenmode plus radial-localization and
// component-breakdown diagnostics. The eigenmode is saved to a .gkyl frame
// so the user can visualize it with their own tooling.
//
// dt is set to 0.4*Δr for SSP-RK3 stability margin (production cfl_frac=1.0
// gives dt~Δr; we shrink to 0.4 to be safe inside ergosphere where lapse is
// smaller and characteristic speeds may differ).
//
// Also runs a low-resolution comparison (32x32) since the user reports low
// resolution inhibits the observed mode in the simulation; this checks
// whether the operator's dominant eigenvalue is itself resolution-dependent.
static void
test_power_iteration(void)
{
  // High-resolution run (matches failing simulation).
  {
    struct budget_bench bench;
    bench_init(&bench, 192, 192, /*use_lax=*/true);

    struct gate_ctx ctx = {
      .lower_r = bench.local.lower[0],
      .lower_theta = bench.local.lower[1],
      .r_lower = R_LOWER,
      .dr = (R_UPPER - R_LOWER) / bench.nr,
      .gate_r_lt = false, .r_max = 0.0,
      .gate_r_gt = false, .r_min = 0.0,
    };
    seed_dQ_pattern(&bench, /*comp=*/5, /*amp=*/1.0, sign_alt_rt, &ctx);

    double dr = (R_UPPER - R_LOWER) / bench.nr;
    double dt = 0.4 * dr;
    power_iterate(&bench, /*num_steps=*/300, dt, /*print_every=*/20,
      "ctest_dg_gr_maxwell_eigenmode_LLF_192.gkyl", NULL);

    print_partition_table(&bench, "POWER ITER converged eigenmode (LLF, 192x192)");
    print_component_breakdown(&bench, "[LLF, 192x192]");
    print_radial_localization(&bench, /*nbins=*/8, "[LLF, 192x192]");

    bench_release(&bench);
  }

  // Low-resolution comparison: user reports low res inhibits the mode in the
  // simulation. Does the operator's dominant eigenvalue change with resolution?
  {
    struct budget_bench bench;
    bench_init(&bench, 32, 32, /*use_lax=*/true);
    struct gate_ctx ctx = {
      .lower_r = bench.local.lower[0],
      .lower_theta = bench.local.lower[1],
      .r_lower = R_LOWER,
      .dr = (R_UPPER - R_LOWER) / bench.nr,
      .gate_r_lt = false, .r_max = 0.0,
      .gate_r_gt = false, .r_min = 0.0,
    };
    seed_dQ_pattern(&bench, /*comp=*/5, /*amp=*/1.0, sign_alt_rt, &ctx);
    double dr = (R_UPPER - R_LOWER) / bench.nr;
    double dt = 0.4 * dr;
    power_iterate(&bench, /*num_steps=*/300, dt, /*print_every=*/20,
      "ctest_dg_gr_maxwell_eigenmode_LLF_032.gkyl", NULL);
    print_partition_table(&bench, "POWER ITER converged eigenmode (LLF, 32x32)");
    print_component_breakdown(&bench, "[LLF, 32x32]");
    print_radial_localization(&bench, /*nbins=*/8, "[LLF, 32x32]");
    bench_release(&bench);
  }
}

// ---- Phase 2A: Resolution scan ----
//
// Run power iteration at several N to map lambda_max(N). The diagnostic
// signal we want is whether
//   lambda_max(N) -> lambda_continuum   (mode is converging, physical),
//   lambda_max(N) -> infinity as 1/Δr   (grid-dependent spurious branch),
// or some mixture. The scan prints lambda_max, Δr, and the product
// lambda_max * Δr. If that product is approximately constant, the mode is
// scaling like 1/Δr (discrete artifact).
static void
test_resolution_scan(void)
{
  // Extended down to N=8 to cross-check with the explicit-spectrum test, and
  // to see whether the unstable branch has a critical resolution below which
  // it does not form (the user reports low resolution inhibits the mode in
  // the simulation).
  int N_list[] = { 8, 12, 16, 24, 32, 48, 64, 96, 128, 192 };
  const int n_N = (int)(sizeof(N_list)/sizeof(N_list[0]));

  printf("\n=== Phase 2A: Resolution scan, lambda_max(N) (LLF, 200 RK3 steps) ===\n");
  printf("    N        Delta_r       lambda_max     lambda*Delta_r       lambda/lambda(16)\n");

  double lambda_at_16 = 0.0;
  for (int i = 0; i < n_N; ++i) {
    int N = N_list[i];
    struct budget_bench bench;
    bench_init(&bench, N, N, /*use_lax=*/true);

    struct gate_ctx ctx = {
      .lower_r = bench.local.lower[0],
      .lower_theta = bench.local.lower[1],
      .r_lower = R_LOWER,
      .dr = (R_UPPER - R_LOWER) / bench.nr,
      .gate_r_lt = false, .r_max = 0.0,
      .gate_r_gt = false, .r_min = 0.0,
    };
    seed_dQ_pattern(&bench, /*comp=*/5, /*amp=*/1.0, sign_alt_rt, &ctx);

    double dr = (R_UPPER - R_LOWER) / bench.nr;
    double dt = 0.4 * dr;
    double lambda_eff = 0.0;
    // Run power iteration silently (print_every larger than num_steps).
    power_iterate(&bench, /*num_steps=*/200, dt, /*print_every=*/10000,
      /*out_fname=*/NULL, /*out_lambda_eff=*/&lambda_eff);

    if (i == 0) lambda_at_16 = lambda_eff;
    double ratio = lambda_at_16 != 0.0 ? lambda_eff/lambda_at_16 : 0.0;
    printf("  %4d   %1.6e   %1.6e   %1.6e   %1.6e\n",
      N, dr, lambda_eff, lambda_eff*dr, ratio);

    bench_release(&bench);
  }
  printf("If lambda*Delta_r is approximately constant -> mode scales as 1/Delta_r (discrete).\n"
         "If lambda_max plateaus -> mode is converging to a continuum eigenvalue.\n");
}

// ---- Phase 2B: Explicit spectrum at low resolution ----

// Forward declaration of LAPACK's non-symmetric eigenvalue routine.
// Fortran ABI; Accelerate (Mac) and OpenBLAS both expose this symbol.
extern void dgeev_(const char *JOBVL, const char *JOBVR, const int *N,
  double *A, const int *LDA, double *WR, double *WI,
  double *VL, const int *LDVL, double *VR, const int *LDVR,
  double *WORK, const int *LWORK, int *INFO);

// Build the dense matrix representation of the assembled DG operator L on the
// interior cells (no ghosts). Layout: row-major linearization of
// (cell_iter_count, component, mode), where cell_iter_count is the order in
// which gkyl_range_iter walks confLocal. The matrix is stored in column-major
// (LAPACK convention).
//
// Each column j is computed by setting unit vector e_j and applying L. Slow
// but conceptually clean: O(N_total^2) cell visits. For the small grid sizes
// we use here (8x8 -> N_total = 2048), this finishes in well under a second.
static int
build_explicit_L_matrix(struct budget_bench *bench, double *L_mat)
{
  int Nr = bench->nr, Nt = bench->ntheta;
  int N_cells = Nr * Nt;
  int N_total = N_cells * NUM_FIELD_COMPS * NUM_BASIS;

  struct gkyl_array *unit  = mkarr(NUM_FIELD_COMPS*NUM_BASIS, bench->local_ext.volume);
  struct gkyl_array *L_unit = mkarr(NUM_FIELD_COMPS*NUM_BASIS, bench->local_ext.volume);

  // Cache cell index list (so we know the j -> (idx, comp, mode) mapping).
  long *cell_lidx = gkyl_malloc(sizeof(long) * N_cells);
  int idx_count = 0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &bench->local);
  while (gkyl_range_iter_next(&iter)) {
    cell_lidx[idx_count++] = gkyl_range_idx(&bench->local, iter.idx);
  }

  for (int j = 0; j < N_total; ++j) {
    int target_cell = j / (NUM_FIELD_COMPS * NUM_BASIS);
    int target_comp = (j / NUM_BASIS) % NUM_FIELD_COMPS;
    int target_mode = j % NUM_BASIS;

    gkyl_array_clear(unit, 0.0);
    double *q = gkyl_array_fetch(unit, cell_lidx[target_cell]);
    q[target_comp*NUM_BASIS + target_mode] = 1.0;

    apply_full_L(bench, unit, L_unit);

    // Read column j of L from L_unit.
    int row = 0;
    for (int c = 0; c < N_cells; ++c) {
      const double *out = gkyl_array_cfetch(L_unit, cell_lidx[c]);
      for (int comp = 0; comp < NUM_FIELD_COMPS; ++comp) {
        for (int m = 0; m < NUM_BASIS; ++m) {
          // Column-major: L_mat[row + j*N_total]
          L_mat[row + j*N_total] = out[comp*NUM_BASIS + m];
          ++row;
        }
      }
    }
  }

  gkyl_free(cell_lidx);
  gkyl_array_release(unit);
  gkyl_array_release(L_unit);

  return N_total;
}

// Diagonalize L via dgeev_ and print the top-K eigenvalues by Re(lambda) AND
// the top-K by |lambda|. The first ranking tells us whether the operator has
// real-positive eigenvalues (a true discrete-artifact unstable mode). The
// second tells us where the largest |lambda| lives — relevant because power
// iteration with SSP-RK3 at fixed CFL converges to that mode, not to the
// largest Re(lambda) mode if the largest |lambda| has |Im|·dt > sqrt(3).
//
// Returns the maximum real part via *out_max_re_lambda (NULL ok), the
// maximum magnitude via *out_max_abs_lambda (NULL ok), and the count of
// eigenvalues with Re > re_threshold via *out_n_unstable (NULL ok).
static void
diagonalize_and_print_spectrum(double *L_mat, int N_total, int top_K,
  double re_threshold,
  double *out_max_re_lambda, double *out_max_abs_lambda, int *out_n_unstable)
{
  // dgeev arguments. We only want eigenvalues; left/right eigenvectors are
  // skipped. Workspace is a query-then-allocate.
  double *WR = gkyl_malloc(sizeof(double) * N_total);
  double *WI = gkyl_malloc(sizeof(double) * N_total);
  double VL_dummy = 0.0, VR_dummy = 0.0;
  int LDVL = 1, LDVR = 1;
  int LDA = N_total;
  int INFO = 0;

  // Workspace query.
  double work_query = 0.0;
  int lwork_query = -1;
  dgeev_("N", "N", &N_total, L_mat, &LDA, WR, WI,
    &VL_dummy, &LDVL, &VR_dummy, &LDVR, &work_query, &lwork_query, &INFO);
  TEST_CHECK(INFO == 0);
  int LWORK = (int)work_query;
  double *WORK = gkyl_malloc(sizeof(double) * LWORK);

  // Actual diagonalization.
  dgeev_("N", "N", &N_total, L_mat, &LDA, WR, WI,
    &VL_dummy, &LDVL, &VR_dummy, &LDVR, WORK, &LWORK, &INFO);
  TEST_CHECK(INFO == 0);
  if (INFO != 0) {
    printf("dgeev_ failed with INFO=%d; aborting spectrum print.\n", INFO);
    gkyl_free(WORK); gkyl_free(WR); gkyl_free(WI);
    return;
  }

  // Selection-sort top K by Re(lambda) descending.
  int *order_re = gkyl_malloc(sizeof(int) * N_total);
  for (int i = 0; i < N_total; ++i) order_re[i] = i;
  for (int k = 0; k < top_K && k < N_total; ++k) {
    int best = k;
    for (int i = k+1; i < N_total; ++i) {
      if (WR[order_re[i]] > WR[order_re[best]]) best = i;
    }
    int tmp = order_re[k]; order_re[k] = order_re[best]; order_re[best] = tmp;
  }
  // Selection-sort top K by |lambda| descending.
  int *order_ab = gkyl_malloc(sizeof(int) * N_total);
  for (int i = 0; i < N_total; ++i) order_ab[i] = i;
  for (int k = 0; k < top_K && k < N_total; ++k) {
    int best = k;
    double bestmag = WR[order_ab[best]]*WR[order_ab[best]] + WI[order_ab[best]]*WI[order_ab[best]];
    for (int i = k+1; i < N_total; ++i) {
      double m = WR[order_ab[i]]*WR[order_ab[i]] + WI[order_ab[i]]*WI[order_ab[i]];
      if (m > bestmag) { best = i; bestmag = m; }
    }
    int tmp = order_ab[k]; order_ab[k] = order_ab[best]; order_ab[best] = tmp;
  }

  printf("\n  Top %d eigenvalues by Re(lambda):\n", top_K);
  printf("    rank        Re(lambda)        Im(lambda)        |lambda|\n");
  for (int k = 0; k < top_K && k < N_total; ++k) {
    int i = order_re[k];
    double mag = sqrt(WR[i]*WR[i] + WI[i]*WI[i]);
    printf("    %4d   %+1.6e   %+1.6e   %1.6e\n", k, WR[i], WI[i], mag);
  }
  printf("\n  Top %d eigenvalues by |lambda|:\n", top_K);
  printf("    rank        Re(lambda)        Im(lambda)        |lambda|\n");
  for (int k = 0; k < top_K && k < N_total; ++k) {
    int i = order_ab[k];
    double mag = sqrt(WR[i]*WR[i] + WI[i]*WI[i]);
    printf("    %4d   %+1.6e   %+1.6e   %1.6e\n", k, WR[i], WI[i], mag);
  }

  if (out_max_re_lambda)  *out_max_re_lambda  = WR[order_re[0]];
  if (out_max_abs_lambda) *out_max_abs_lambda =
    sqrt(WR[order_ab[0]]*WR[order_ab[0]] + WI[order_ab[0]]*WI[order_ab[0]]);

  // Count eigenvalues with Re > re_threshold. Useful to see whether the
  // unstable spectrum is filling out (continuum-like) or sparse (discrete).
  if (out_n_unstable) {
    int n_unstable = 0;
    for (int i = 0; i < N_total; ++i) if (WR[i] > re_threshold) n_unstable++;
    *out_n_unstable = n_unstable;
  }

  gkyl_free(order_re);
  gkyl_free(order_ab);
  gkyl_free(WORK);
  gkyl_free(WR);
  gkyl_free(WI);
}

// Build L explicitly at several small N, diagonalize, and report
//   max_Re(lambda)(N): the largest real part anywhere in the spectrum.
//   max_|lambda|(N):   the largest magnitude (typically a near-imaginary mode).
// Decision tree:
//   - max_Re(lambda) GROWS with N (e.g. ~0.004 at N=8, ~30 at N=12, ~200 at N=16):
//     genuine discrete-artifact unstable mode, scaling with 1/Δr.
//   - max_Re(lambda) STAYS small (~ 1e-3 at all N), but max_|lambda| ~ 1/Δr:
//     the operator is essentially continuum-stable. The "instability" the
//     simulation sees is a TIME-STEPPING problem (RK3 amplifying high-|Im|
//     modes when CFL is set near the imaginary-axis stability boundary).
// The two scenarios call for completely different fixes: scenario A is a DG
// scheme issue; scenario B is a CFL/timestepper issue (or both, if Re(lambda)
// is nonzero but small enough that timestepper amplifies the Im modes too).
//
// Resolutions: N=4 (matrix 512), N=6 (1152), N=8 (2048), N=12 (4608).
// At N=12 dgeev takes a couple of minutes but is still feasible.
static void
test_explicit_spectrum(void)
{
  // Scan N up to where dgeev is feasible:
  //   N=16 -> 8192^2  ~ 530 MB,  ~3 min
  //   N=20 -> 12800^2 ~ 1.3 GB,  ~12 min
  //   N=24 -> 18432^2 ~ 2.7 GB,  ~35 min
  // Total wall-clock for the scan is about 50 minutes. Larger N starts
  // running into laptop memory limits (N=32 -> 8.6 GB).
  int N_list[] = { 6, 8, 12, 16, 20, 24 };
  const int n_N = (int)(sizeof(N_list)/sizeof(N_list[0]));

  printf("\n=== Phase 2B: Explicit spectrum scan, max Re(lambda) and max |lambda| vs N ===\n");

  // Summary at the bottom.
  double max_re[16] = {0}, max_ab[16] = {0};
  int n_unstable[16] = {0};

  for (int idx = 0; idx < n_N; ++idx) {
    int N = N_list[idx];
    struct budget_bench bench;
    bench_init(&bench, N, N, /*use_lax=*/true);

    int N_cells = N*N;
    int N_total = N_cells * NUM_FIELD_COMPS * NUM_BASIS;

    printf("\n--- N = %d x %d  (matrix size = %d) ---\n", N, N, N_total);
    fflush(stdout);
    double *L_mat = gkyl_malloc(sizeof(double) * (long)N_total * (long)N_total);
    build_explicit_L_matrix(&bench, L_mat);

    double trace = 0.0;
    for (int i = 0; i < N_total; ++i) trace += L_mat[i + i*N_total];
    printf("  trace(L) = %+1.6e (sum of all eigenvalues' real parts)\n", trace);
    fflush(stdout);

    // Threshold for "non-trivially unstable": 1e-3. This filters out the
    // O(1e-15) noise eigenvalues and counts only modes that would matter
    // dynamically.
    diagonalize_and_print_spectrum(L_mat, N_total, /*top_K=*/15,
      /*re_threshold=*/1.0e-3,
      &max_re[idx], &max_ab[idx], &n_unstable[idx]);
    fflush(stdout);
    gkyl_free(L_mat);

    bench_release(&bench);
  }

  // Final summary, including a power-law fit to max Re(lambda) so we can see
  // whether the slope d(log Re)/d(log N) is decreasing (continuum
  // convergence) or stable/increasing (discrete artifact).
  printf("\n=== Phase 2B summary ===\n");
  printf("    N      Delta_r        max Re(lambda)        max |lambda|     Re*Delta_r        |L|*Delta_r     # Re>1e-3\n");
  for (int idx = 0; idx < n_N; ++idx) {
    int N = N_list[idx];
    double dr = (R_UPPER - R_LOWER) / N;
    printf("  %4d   %1.6e     %+1.6e      %1.6e     %+1.6e      %1.6e     %d\n",
      N, dr, max_re[idx], max_ab[idx], max_re[idx]*dr, max_ab[idx]*dr,
      n_unstable[idx]);
  }
  printf("\n  Local power-law slope p where Re_lambda(N) ~ N^p between consecutive N values:\n");
  for (int idx = 1; idx < n_N; ++idx) {
    if (max_re[idx] <= 0.0 || max_re[idx-1] <= 0.0) continue;
    double r = log(max_re[idx]/max_re[idx-1]) / log((double)N_list[idx]/(double)N_list[idx-1]);
    printf("    p(%d -> %d) = %1.3f\n", N_list[idx-1], N_list[idx], r);
  }
  printf("\n  Decreasing p sequence -> continuum mode being resolved (slope -> 0).\n");
  printf("  Constant p ~ 1 -> 1/Delta_r divergence (discrete artifact).\n");
  printf("  If # Re>1e-3 grows fast with N, the spectrum is filling like a continuum.\n");
}

// ---- Phase 2C: BC- and flux-dependence scan ----
//
// At a fixed N, build L explicitly and diagonalize for four BC configurations
// (A: ref r + ref theta-pole; B: outflow r + ref theta-pole; C: ref r +
// non-reflective theta-pole; D: outflow r + non-reflective theta-pole) under
// each of two flux choices (LLF / Roe). The flux comparison is informative
// because:
//   - LLF adds explicit dissipation |alpha_max|·(U_R-U_L). If this is mostly
//     masking an interior instability (rather than correctly damping it in
//     the right norm), Roe — which has no extra dissipation beyond the
//     characteristic-correct upwinding — should expose more of the underlying
//     real-positive spectrum, i.e., Roe would yield larger max Re(lambda)
//     than LLF for the same BC config.
//   - Conversely, if Roe gives smaller max Re(lambda) than LLF, the LLF
//     dissipation is itself the destabilizing piece (it has wrong sign in
//     the curved-space energy norm where the instability lives).
// Either outcome implicates a kernel-level fix (curved-norm dissipation),
// not a BC-level fix.
static void
test_explicit_spectrum_bc_scan(void)
{
  // N=12 keeps each dgeev under 30 sec; total runtime ~6-7 min for 12 configs
  // (4 BC variants x 3 flux variants: LLF, Roe, curved-norm LLF).
  const int N = 12;
  struct cfg { int theta_pole, outflow_r, outflow_theta; bool use_lax; bool use_curved_norm; const char *name; };
  static const struct cfg configs[] = {
    { 1, 0, 0, true,  false, "A_LLF      : ref r + ref theta-pole, LLF" },
    { 1, 0, 0, false, false, "A_Roe      : ref r + ref theta-pole, Roe" },
    { 1, 0, 0, true,  true,  "A_LLFcurved: ref r + ref theta-pole, LLF curved-norm" },
    { 1, 1, 0, true,  false, "B_LLF      : outflow r + ref theta-pole, LLF" },
    { 1, 1, 0, false, false, "B_Roe      : outflow r + ref theta-pole, Roe" },
    { 1, 1, 0, true,  true,  "B_LLFcurved: outflow r + ref theta-pole, LLF curved-norm" },
    { 1, 0, 1, true,  false, "C_LLF      : ref r + non-refl theta, LLF" },
    { 1, 0, 1, false, false, "C_Roe      : ref r + non-refl theta, Roe" },
    { 1, 0, 1, true,  true,  "C_LLFcurved: ref r + non-refl theta, LLF curved-norm" },
    { 1, 1, 1, true,  false, "D_LLF      : outflow r + non-refl theta, LLF" },
    { 1, 1, 1, false, false, "D_Roe      : outflow r + non-refl theta, Roe" },
    { 1, 1, 1, true,  true,  "D_LLFcurved: outflow r + non-refl theta, LLF curved-norm" },
  };
  const int n_cfg = (int)(sizeof(configs)/sizeof(configs[0]));

  printf("\n=== Phase 2C: BC- and flux-dependence scan at N=%d ===\n", N);
  double max_re[16] = {0}, max_ab[16] = {0};
  int n_unstable[16] = {0};

  for (int idx = 0; idx < n_cfg; ++idx) {
    printf("\n--- Config %s ---\n", configs[idx].name);
    fflush(stdout);

    struct budget_bench bench;
    bench_init_full_v2(&bench, N, N, configs[idx].use_lax,
      configs[idx].theta_pole, configs[idx].outflow_r, configs[idx].outflow_theta,
      configs[idx].use_curved_norm);

    int N_total = N*N * NUM_FIELD_COMPS * NUM_BASIS;
    double *L_mat = gkyl_malloc(sizeof(double) * (long)N_total * (long)N_total);
    build_explicit_L_matrix(&bench, L_mat);

    double trace = 0.0;
    for (int i = 0; i < N_total; ++i) trace += L_mat[i + i*N_total];
    printf("  trace(L) = %+1.6e\n", trace);
    fflush(stdout);

    diagonalize_and_print_spectrum(L_mat, N_total, /*top_K=*/10, /*re_threshold=*/1.0e-3,
      &max_re[idx], &max_ab[idx], &n_unstable[idx]);
    fflush(stdout);
    gkyl_free(L_mat);

    bench_release(&bench);
  }

  printf("\n=== Phase 2C summary (N=%d) ===\n", N);
  printf("    config                                                       max Re(lambda)        max |lambda|     # Re>1e-3\n");
  for (int idx = 0; idx < n_cfg; ++idx) {
    printf("    %-60s %+1.6e      %1.6e     %d\n",
      configs[idx].name, max_re[idx], max_ab[idx], n_unstable[idx]);
  }
  printf("\n  Reading the cross-product (rows = BC config, columns = flux choice):\n");
  printf("    BC      LLF max Re       Roe max Re       LLF-curved max Re      Roe/LLF      Curved/LLF\n");
  for (int b = 0; b < n_cfg; b += 3) {
    double r_llf = max_re[b], r_roe = max_re[b+1], r_curved = max_re[b+2];
    double rr_ratio = (r_llf > 0.0) ? r_roe/r_llf : 0.0;
    double rc_ratio = (r_llf > 0.0) ? r_curved/r_llf : 0.0;
    char tag[8] = "?"; tag[0] = configs[b].name[0]; tag[1] = '\0';
    printf("    %s    %+1.6e    %+1.6e    %+1.6e        %1.3f       %1.4f\n",
      tag, r_llf, r_roe, r_curved, rr_ratio, rc_ratio);
  }
  printf("\n  Interpretation:\n");
  printf("    - Roe/LLF > 1 -> standard LLF dissipation was masking interior instability\n");
  printf("                     (Roe with no extra dissipation exposes more).\n");
  printf("    - Curved/LLF << 1 -> curved-norm dissipation kills more of the unstable\n");
  printf("                          spectrum than standard LLF, confirming the\n");
  printf("                          interior instability lives in the curved norm.\n");
  printf("    - Curved max Re ~ 0 across BC configs -> the curved-norm patch fully\n");
  printf("                          stabilizes the interior; remaining max Re comes\n");
  printf("                          only from the BCs.\n");
}

TEST_LIST = {
  { "bench_smoke", test_bench_smoke },
  { "seeded_sweep", test_seeded_sweep },
  { "power_iteration", test_power_iteration },
  { "resolution_scan", test_resolution_scan },
  { "explicit_spectrum", test_explicit_spectrum },
  { "spectrum_bc_scan", test_explicit_spectrum_bc_scan },
  { NULL, NULL },
};
