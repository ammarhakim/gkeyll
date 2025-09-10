#include <gkyl_alloc.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_dg_basis_ops.h>
#include <gkyl_math.h>
#include <gkyl_mirror_grid_gen.h>
#include <gkyl_rect_decomp.h>

struct gkyl_mirror_grid_gen_x {
  enum gkyl_mirror_grid_gen_field_line_coord fl_coord; // field-line coordinate to use
  bool include_axis; // add nodes on r=0 axis (the axis is assumed be psi=0)
};

// context for use in root finder
struct psirz_ctx {
  struct gkyl_basis_ops_evalf *evcub; // cubic eval functions
  double Z; // local Z value
  double psi; // psi to match
};

static inline double
floor_sqrt(double x)
{
  return sqrt( fmax(x, 1e-14) );
}

static
double psirz(double R, void *ctx)
{
  struct psirz_ctx *rctx = ctx;
  double Z = rctx->Z;
  double xn[2] = { R, Z };
  double fout[1];
  rctx->evcub->eval_cubic(0, xn, fout, rctx->evcub->ctx);
  return fout[0] - rctx->psi;
}

struct gkyl_mirror_grid_gen *
gkyl_mirror_grid_gen_inew(const struct gkyl_mirror_grid_gen_inp *inp)
{
  struct gkyl_mirror_grid_gen *geo = gkyl_malloc(sizeof *geo);
  geo->gg_x = gkyl_malloc(sizeof *geo->gg_x);

  geo->gg_x->fl_coord = inp->fl_coord;
  geo->gg_x->include_axis = inp->include_axis;

  int nr = inp->nrcells, nz = inp->nzcells;
  int cells[] = { nr, nz };
  double lower[2] = { inp->R[0], inp->Z[0] };
  double upper[2] = { inp->R[1], inp->Z[1] };

  struct gkyl_rect_grid gridRZ;
  gkyl_rect_grid_init(&gridRZ, 2, lower, upper, cells);

  struct gkyl_basis_ops_evalf *evcub =
    gkyl_dg_basis_ops_evalf_new(&gridRZ, inp->psiRZ);

  do {
    const char *fname = inp->psi_cubic_fname ? inp->psi_cubic_fname : "psi_cubic.gkyl";
    if (inp->write_psi_cubic)
      gkyl_dg_basis_ops_evalf_write_cubic(evcub, fname);
  } while (0);

  // Construct grid in RZ plane
  
  enum { NPSI, NZ };
  int nc[2];
  nc[NPSI] = inp->comp_grid->cells[0]+1;
  nc[NZ] = inp->comp_grid->cells[2]+1;
  
  long nctot = nc[NPSI]*nc[NZ];
  geo->nodes_rz = gkyl_array_new(GKYL_DOUBLE, 2, nctot);
  geo->nodes_psi = gkyl_array_new(GKYL_DOUBLE, 1, nctot);
  geo->nodes_geom = gkyl_array_new(GKYL_USER, sizeof(struct gkyl_mirror_grid_gen_geom), nctot);

  struct gkyl_range node_rng;
  gkyl_range_init_from_shape(&node_rng, 2, nc);

  double z_lo = inp->comp_grid->lower[2];
  double z_up = inp->comp_grid->upper[2];
  double dz = (z_up-z_lo)/(nc[NZ]-1);

  double psi_lo = inp->comp_grid->lower[0];
  double psi_up = inp->comp_grid->upper[0];

  bool inc_axis = inp->include_axis;

  // Adjust if we are using sqrt(psi) as radial coordinate
  double psic_lo = psi_lo, psic_up = psi_up;
  if (inp->fl_coord == GKYL_MIRROR_GRID_GEN_SQRT_PSI_CART_Z) {
    // When we include axis, psi_lo is ignored
    psic_lo = inc_axis ? 0.0 : sqrt(psi_lo);
    psic_up = sqrt(psi_up);
  }
  
  double dpsi = (psic_up-psic_lo)/(nc[NPSI]-1);

  double rlow = lower[0], rup = upper[0];
  double rmin = rlow + 1e-8*(rup-rlow);

  struct psirz_ctx pctx = { .evcub = evcub };

  // Compute node locations
  bool status = true;
  for (int iz=0; iz<nc[NZ]; ++iz) {
    double zcurr = z_lo + iz*dz;

    double psi_min[1], psi_max[1];
    evcub->eval_cubic(0.0, (double[2]) { rmin, zcurr }, psi_min, evcub->ctx);
    evcub->eval_cubic(0.0, (double[2]) { rup, zcurr }, psi_max, evcub->ctx);

    for (int ipsi=0; ipsi<nc[NPSI]; ++ipsi) {

      if (inc_axis && (ipsi == 0)) {
        int idx[2] = { ipsi, iz };
        double *rz = gkyl_array_fetch(geo->nodes_rz, gkyl_range_idx(&node_rng, idx));
        rz[0] = 0.0; rz[1] = zcurr;
      }
      else {
        double psic_curr = psic_lo + ipsi*dpsi;
        pctx.Z = zcurr;

        // We continue to do root-finding for psi and not sqrt(psi)
        double psi_curr = psic_curr;
        if (inp->fl_coord == GKYL_MIRROR_GRID_GEN_SQRT_PSI_CART_Z)
          psi_curr = psic_curr*psic_curr;
        
        pctx.psi =  psi_curr;
        
        struct gkyl_qr_res root = gkyl_ridders(psirz, &pctx, rmin, rup,
          psi_min[0]-psi_curr, psi_max[0]-psi_curr,
          100, 1e-10);

        if (root.status) {
          status = false;
          goto cleanup;
        }
        
        int idx[2] = { ipsi, iz };
        double *rz = gkyl_array_fetch(geo->nodes_rz, gkyl_range_idx(&node_rng, idx));
        rz[0] = root.res; rz[1] = zcurr;
        double *psi_coord = gkyl_array_fetch(geo->nodes_psi, gkyl_range_idx(&node_rng, idx));
        psi_coord[0] = psi_curr; // store psi coordinate
      }
    }
  }

  enum { PSI_I, DPSI_R_I, DPSI_Z_I };
  
  // Compute geometry at nodes
  for (int iz=0; iz<nc[NZ]; ++iz) {
    
    for (int ipsi=0; ipsi<nc[NPSI]; ++ipsi) {
      int idx[2] = { ipsi, iz };
      long loc = gkyl_range_idx(&node_rng, idx);
      
      const double *rzp = gkyl_array_cfetch(geo->nodes_rz, loc);
      double rz[2] = { rzp[0], rzp[1] };
      
      struct gkyl_mirror_grid_gen_geom *g = gkyl_array_fetch(geo->nodes_geom, loc);
      
      if (inc_axis && (ipsi == 0)) {
        double fout2[4]; // second derivative of psi is needed
        evcub->eval_cubic_wgrad2(0.0, rz, fout2, evcub->ctx);

        // On-axis the coordinate system breaks down. Below we choose
        // some reasonable defaults for the tnagent and
        // duals. However, the Jacobians and magnetic field are
        // correct and computed using the estimated asymptotic
        // behavior of psi as r -> 0.
        
        g->dual[0].x[0] = 1.0;
        g->dual[0].x[1] = 0.0;
        g->dual[0].x[2] = 0.0;

        g->dual[1].x[0] = 0;
        g->dual[1].x[1] = 1.0;
        g->dual[1].x[2] = 0.0;

        g->dual[2].x[0] = 0;
        g->dual[2].x[1] = 0.0;
        g->dual[2].x[2] = 1.0;        
        
        g->tang[0].x[0] = 1.0;
        g->tang[0].x[1] = 0.0;
        g->tang[0].x[2] = 0.0;

        g->tang[1].x[0] = 0.0;
        g->tang[1].x[1] = 1.0;
        g->tang[1].x[2] = 0.0;        
        
        g->tang[2].x[0] = 0;
        g->tang[2].x[1] = 0.0;
        g->tang[2].x[2] = 1.0;
        
        if (inp->fl_coord == GKYL_MIRROR_GRID_GEN_SQRT_PSI_CART_Z)
          g->Jc = 0; // assumes asymptotics of psi ~ r^2 as r -> 0
        else
          g->Jc = 1/fout2[DPSI_R_I];

        g->B.x[0] = 0.0; // no radial component
        g->B.x[1] = 0.0;
        g->B.x[2] = fout2[DPSI_R_I]; // diff(psi,r,2)
      }
      else {
        double fout[3]; // first derivative of psi is needed
        evcub->eval_cubic_wgrad(0.0, rz, fout, evcub->ctx);
      
        // e^1
        g->dual[0].x[0] = fout[DPSI_R_I]; // dpsi/dr
        g->dual[0].x[1] = 0.0; // no toroidal component
        g->dual[0].x[2] = fout[DPSI_Z_I]; // dspi/dz
      
        if (inp->fl_coord == GKYL_MIRROR_GRID_GEN_SQRT_PSI_CART_Z) {
          // for sqrt(psi) as radial coordinate e^1 = grad(psi)/2*sqrt(psi)
          g->dual[0].x[0] = g->dual[0].x[0]/(2*floor_sqrt(fout[0]));
          g->dual[0].x[2] = g->dual[0].x[2]/(2*floor_sqrt(fout[0]));
        }

        // e^2 is just e^phi
        g->dual[1].x[0] = 0;
        g->dual[1].x[1] = 1.0/(rz[0]*rz[0]);
        g->dual[1].x[2] = 0.0;

        // e^3 is just sigma_3
        g->dual[2].x[0] = 0;
        g->dual[2].x[1] = 0.0;
        g->dual[2].x[2] = 1.0;

        // e_1 points along the radial direction
        g->tang[0].x[0] = 1/g->dual[0].x[0];
        g->tang[0].x[1] = 0.0;
        g->tang[0].x[2] = 0.0;

        // e_2
        g->tang[1].x[0] = 0;
        g->tang[1].x[1] = 1.0;
        g->tang[1].x[2] = 0.0;

        // e_3
        g->tang[2].x[0] = -fout[DPSI_Z_I]/fout[DPSI_R_I];
        g->tang[2].x[1] = 0.0;
        g->tang[2].x[2] = 1.0;
        
        if (inp->fl_coord == GKYL_MIRROR_GRID_GEN_SQRT_PSI_CART_Z)
          g->Jc = 2*floor_sqrt(fout[PSI_I])*rz[0]/fout[DPSI_R_I];
        else
          g->Jc = rz[0]/fout[DPSI_R_I];
        
        g->B.x[0] = -fout[DPSI_Z_I]/rz[0];
        g->B.x[1] = 0.0;
        g->B.x[2] = fout[DPSI_R_I]/rz[0];
      }
    }
  }
  
  cleanup:

  if (true != status) {
    gkyl_mirror_grid_gen_release(geo);
    geo = 0;
    fprintf(stderr, "gkyl_mirror_grid_gen_inew failed to generate a grid\n");
  }
  
  gkyl_dg_basis_ops_evalf_release(evcub);
  
  return geo;
}

bool
gkyl_mirror_grid_gen_is_include_axis(const struct gkyl_mirror_grid_gen *geom)
{
  return geom->gg_x->include_axis;
}

enum gkyl_mirror_grid_gen_field_line_coord
  gkyl_mirror_grid_gen_fl_coord(const struct gkyl_mirror_grid_gen *geom)
{
  return geom->gg_x->fl_coord;
}
  
void
gkyl_mirror_grid_gen_release(struct gkyl_mirror_grid_gen *geom)
{
  gkyl_array_release(geom->nodes_rz);
  gkyl_array_release(geom->nodes_psi);
  gkyl_array_release(geom->nodes_geom);
  gkyl_free(geom->gg_x);
  gkyl_free(geom);
}
