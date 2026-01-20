#include <gkyl_vlasov_priv.h>
#include <gkyl_vlasov_triad_geom.h>

void
eval_annulus_vierbein_2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Downstairs components of the Vierbeins
  double q_r = xn[0];

  // e_ij = e_i(x) . \sigma_j(x)
  double e_rr = 1.0; // Vierbein Coefficients (r-r coefficient).
  double e_rt = 0.0; // Vierbein Coefficients (r-theta coefficient).
  double e_tr = 0.0; // Vierbein Coefficients (theta-r coefficient).
  double e_tt = q_r; // Vierbein Coefficients (theta-theta coefficient).
  
  fout[0] = e_rr;
  fout[1] = e_rt;
  fout[2] = e_tr;
  fout[3] = e_tt;
}

void
eval_annulus_vierbein_gradient_2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Downstairs components of the Vierbein Gradients 
  double q_r = xn[0];

  // d (e_ij) / dr 
  double d_e_rr_dr = 0.0; // Vierbein Gradient Coefficients (r-r coefficient).
  double d_e_rt_dr = 0.0; // Vierbein Gradient Coefficients (r-theta coefficient).
  double d_e_tr_dr = 0.0; // Vierbein Gradient Coefficients (theta-r coefficient).
  double d_e_tt_dr = 1.0; // Vierbein Gradient Coefficients (theta-theta coefficient).

  // d (e_ij) / dtheta 
  double d_e_rr_dt = 0.0; // Vierbein Gradient Coefficients (r-r coefficient).
  double d_e_rt_dt = 0.0; // Vierbein Gradient Coefficients (r-theta coefficient).
  double d_e_tr_dt = 0.0; // Vierbein Gradient Coefficients (theta-r coefficient).
  double d_e_tt_dt = 0.0; // Vierbein Gradient Coefficients (theta-theta coefficient).
  
  fout[0] = d_e_rr_dr;
  fout[1] = d_e_rt_dr;
  fout[2] = d_e_tr_dr;
  fout[3] = d_e_tt_dr;

  fout[4] = d_e_rr_dt;
  fout[5] = d_e_rt_dt;
  fout[6] = d_e_tr_dt;
  fout[7] = d_e_tt_dt;

}

void
eval_ks_rphi_vierbein_2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Parameters
  const struct vm_geom *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;

  // Coordinates
  double r = xn[0];
  double phi = xn[1];

  // Intermediate Variables
  double rho_sq = r * r ;

  // Metric spatial covariant components
  double h_rr = ( 1.0 + 2.0 * M * r / rho_sq);
  double h_pp = ( rho_sq + a * a * ( 1.0 + 2.0 * M * r / rho_sq )   );
  double h_thth = rho_sq;
  double h_rp = - a * ( 1.0 + 2.0 * M * r / rho_sq );

  // Vierbein: e_i^a = g_i . sigma^a
  fout[0] = sqrt( h_rr );
  fout[1] = 0.0;
  fout[2] = h_rp / sqrt( h_rr );
  fout[3] = sqrt( h_pp - h_rp * h_rp / h_rr );
}

void
eval_ks_rphi_vierbein_gradient_2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  const struct vm_geom *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;

  double r = xn[0];
  double phi = xn[1];

  // Intermediate Variables
  double rho_sq = r * r ;

  // Gradient w.r.t. r: d(e_i^a)/dr
  fout[0] = (M * ( - r * r)) / (pow(rho_sq, 1.5) * sqrt(rho_sq + 2.0 * M * r));
  fout[1] = 0.0;
  fout[2] = -(M * a * ( - r * r)) / (pow(rho_sq, 1.5) * sqrt(rho_sq + 2.0 * M * r));
  fout[3] = r / sqrt(rho_sq);

  // Gradient w.r.t. theta: d(e_i^a)/dphi
  fout[4] = 0.0;
  fout[5] = 0.0;
  fout[6] = 0.0;
  fout[7] = 0.0;

}

void 
eval_ks_vierbein_3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  // Parameters
  const struct vm_geom *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;

  // Coordinates
  double r = xn[0];
  double theta = xn[1];
  double phi = xn[2];

  // Intermediate Variables
  double rho_sq = r * r + a * a * cos(theta) * cos(theta);

  // Metric spatial covariant components
  double h_rr = ( 1.0 + 2.0 * M * r / rho_sq);
  double h_pp = sin(theta) * sin(theta) * ( rho_sq + a * a * ( 1.0 + 2.0 * M * r / rho_sq ) * sin(theta) * sin(theta)  );
  double h_thth = rho_sq;
  double h_rp = - a * ( 1.0 + 2.0 * M * r / rho_sq ) * sin(theta) * sin(theta);

  // Vierbein: e_i^a = g_i . sigma^a
  fout[0] = sqrt( h_rr );
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 0.0;
  fout[4] = sqrt( rho_sq );
  fout[5] = 0.0;
  fout[6] = h_rp / sqrt( h_rr );
  fout[7] = 0.0;
  fout[8] = sqrt( h_pp - h_rp * h_rp / h_rr );
}

void 
eval_ks_vierbein_gradient_3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  const struct vm_geom *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;

  double r = xn[0];
  double theta = xn[1];
  double phi = xn[2];

  // Intermediate Variables
  double rho_sq = r * r + a * a * cos(theta) * cos(theta);

  // Gradient w.r.t. r: d(e_i^a)/dr
  fout[0] = (M * (a * a * cos(theta) * cos(theta) - r * r)) / (pow(rho_sq, 1.5) * sqrt(rho_sq + 2.0 * M * r));
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 0.0;
  fout[4] = r / sqrt(rho_sq);
  fout[5] = 0.0;
  fout[6] = -(M * a * sin(theta) * sin(theta) * (a * a * cos(theta) * cos(theta) - r * r)) / (pow(rho_sq, 1.5) * sqrt(rho_sq + 2.0 * M * r));
  fout[7] = 0.0;
  fout[8] = (r * sin(theta)) / sqrt(rho_sq);

  // Gradient w.r.t. theta: d(e_i^a)/dtheta
  fout[9]  = (M * a * a * r * sin(2.0 * theta)) / (pow(rho_sq, 1.5) * sqrt(rho_sq + 2.0 * M * r));
  fout[10] = 0.0;
  fout[11] = 0.0;
  fout[12] = 0.0;
  fout[13] = -(a * a * sin(2.0 * theta)) / (2.0 * sqrt(rho_sq));
  fout[14] = 0.0;
  fout[15] = -(2.0 * a * cos(theta) * sin(theta) * (a * a * a * a * cos(theta) * cos(theta) * cos(theta) * cos(theta) + 2.0 * M * r * r * r + r * r * r * r + 2.0 * a * a * r * r * cos(theta) * cos(theta) + 2.0 * M * a * a * r * cos(theta) * cos(theta) + M * a * a * r * sin(theta) * sin(theta))) / (pow(rho_sq, 1.5) * sqrt(rho_sq + 2.0 * M * r));
  fout[16] = 0.0;
  fout[17] = (sin(4.0 * theta) * a * a + 2.0 * sin(2.0 * theta) * r * r) / (2.0 * sqrt(1.0 - cos(2.0 * theta)) * sqrt(a * a * cos(2.0 * theta) + a * a + 2.0 * r * r));

  // Gradient w.r.t. theta: d(e_i^a)/dphi
  fout[18] = 0.0;
  fout[19] = 0.0;
  fout[20] = 0.0;
  fout[21] = 0.0;
  fout[22] = 0.0;
  fout[23] = 0.0;
  fout[24] = 0.0;
  fout[25] = 0.0;
  fout[26] = 0.0;
}


void
eval_rz_cylindrical_vierbein_3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Downstairs components of the Vierbeins
  double q_r = xn[0];

  // e_ij = e_i(x) . \sigma_j(x)
  double e_rr = 1.0; // Vierbein Coefficients (r-r coefficient).
  double e_rz = 0.0;
  double e_rt = 0.0; // Vierbein Coefficients (r-theta coefficient).
  
  double e_zr = 0.0;
  double e_zz = 1.0;
  double e_zt = 0.0; 

  double e_tr = 0.0; // Vierbein Coefficients (theta-r coefficient).
  double e_tz = 0.0; 
  double e_tt = q_r; // Vierbein Coefficients (theta-theta coefficient).
  
  fout[0] = e_rr;
  fout[1] = e_rz;
  fout[2] = e_rt;

  fout[3] = e_zr;
  fout[4] = e_zz;
  fout[5] = e_zt;

  fout[6] = e_tr;
  fout[7] = e_tz;
  fout[8] = e_tt;
}

void
eval_rz_cylindrical_vierbein_gradient_3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Downstairs components of the Vierbein Gradients 
  double q_r = xn[0];

  // d (e_ij) / dr 
  double d_e_rr_dr = 0.0; // Vierbein Gradient Coefficients (r-r coefficient).
  double d_e_rz_dr = 0.0;
  double d_e_rt_dr = 0.0; // Vierbein Gradient Coefficients (r-theta coefficient).
  
  double d_e_zr_dr = 0.0;
  double d_e_zz_dr = 0.0;
  double d_e_zt_dr = 0.0;

  double d_e_tr_dr = 0.0; // Vierbein Gradient Coefficients (theta-r coefficient).
  double d_e_tz_dr = 0.0;
  double d_e_tt_dr = 1.0; // Vierbein Gradient Coefficients (theta-theta coefficient).

  // d(e_ij) / dz
  double d_e_rr_dz = 0.0; // Vierbein Gradient Coefficients (r-r coefficient).
  double d_e_rz_dz = 0.0;
  double d_e_rt_dz = 0.0; // Vierbein Gradient Coefficients (r-theta coefficient).
  
  double d_e_zr_dz = 0.0;
  double d_e_zz_dz = 0.0;
  double d_e_zt_dz = 0.0;

  double d_e_tr_dz = 0.0; // Vierbein Gradient Coefficients (theta-r coefficient).
  double d_e_tz_dz = 0.0;
  double d_e_tt_dz = 0.0; // Vierbein Gradient Coefficients (theta-theta coefficient).


  // d (e_ij) / dtheta 
  double d_e_rr_dt = 0.0; // Vierbein Gradient Coefficients (r-r coefficient).
  double d_e_rz_dt = 0.0;
  double d_e_rt_dt = 0.0; // Vierbein Gradient Coefficients (r-theta coefficient).

  double d_e_zr_dt = 0.0;
  double d_e_zz_dt = 0.0;
  double d_e_zt_dt = 0.0;

  double d_e_tr_dt = 0.0; // Vierbein Gradient Coefficients (theta-r coefficient).
  double d_e_tz_dt = 0.0;
  double d_e_tt_dt = 0.0; // Vierbein Gradient Coefficients (theta-theta coefficient).
  
  fout[0] = d_e_rr_dr;
  fout[1] = d_e_rz_dr;
  fout[2] = d_e_rt_dr;
  fout[3] = d_e_zr_dr;
  fout[4] = d_e_zz_dr;
  fout[5] = d_e_zt_dr;
  fout[6] = d_e_tr_dr;
  fout[7] = d_e_tz_dr;
  fout[8] = d_e_tt_dr;

  fout[9+0] = d_e_rr_dz;
  fout[9+1] = d_e_rz_dz;
  fout[9+2] = d_e_rt_dz;
  fout[9+3] = d_e_zr_dz;
  fout[9+4] = d_e_zz_dz;
  fout[9+5] = d_e_zt_dz;
  fout[9+6] = d_e_tr_dz;
  fout[9+7] = d_e_tz_dz;
  fout[9+8] = d_e_tt_dz;

  fout[18+0] = d_e_rr_dt;
  fout[18+1] = d_e_rz_dt;
  fout[18+2] = d_e_rt_dt;
  fout[18+3] = d_e_zr_dt;
  fout[18+4] = d_e_zz_dt;
  fout[18+5] = d_e_zt_dt;
  fout[18+6] = d_e_tr_dt;
  fout[18+7] = d_e_tz_dt;
  fout[18+8] = d_e_tt_dt;

}

void
eval_ks_rphi_hamil_2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Grab the geometry
  const struct vm_geom *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;
  
  double r = xn[0];
  double phi = xn[1];

  double pr_hat = xn[2];
  double pphi_hat = xn[3];

  // Build the Hamiltonian for Kerr-Schild in 3D
  double rho_sq = r * r ;
  double H = (1.0/ sqrt(1 + 2*M*r/rho_sq)) * (
        sqrt(1 + pr_hat * pr_hat + pphi_hat * pphi_hat) 
        - (2*M*r/rho_sq) * pr_hat );

  fout[0] = H;
}


void
eval_ks_hamil_3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Grab the geometry
  const struct vm_geom *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;
  
  double r = xn[0];
  double theta = xn[1];
  double phi = xn[2];

  double pr_hat = xn[3];
  double ptheta_hat = xn[4];
  double pphi_hat = xn[5];

  // Build the Hamiltonian for Kerr-Schild in 3D
  double rho_sq = r * r + a * a * cos(theta) * cos(theta);
  double H = (1.0/ sqrt(1 + 2*M*r/rho_sq)) * (
        sqrt(1 + pr_hat * pr_hat + ptheta_hat * ptheta_hat + pphi_hat * pphi_hat) 
        - (2*M*r/rho_sq) * pr_hat );

  fout[0] = H;
}

// for use in kernel tables
typedef struct { evalf_t kernels[3]; } hamil_kern_list;
typedef struct { evalf_t kernels[3]; } vierbein_kern_list;
typedef struct { evalf_t kernels[3]; } vierbein_gradient_kern_list;

static const hamil_kern_list ks_rphi_hamil_list[] = {
  { NULL },
  { eval_ks_rphi_hamil_2v },
  { NULL }
};

static const hamil_kern_list ks_hamil_3v_list[] = {
  { NULL },
  { NULL },
  { eval_ks_hamil_3v }
};

static const vierbein_kern_list annulus_vierbein_list[] = {
  { NULL },
  { eval_annulus_vierbein_2v },
  { NULL }
};

static const vierbein_gradient_kern_list annulus_vierbein_gradient_list[] = {
  { NULL },
  { eval_annulus_vierbein_gradient_2v },
  { NULL }
};

static const vierbein_kern_list ks_rphi_vierbein_list[] = {
  { NULL },
  { eval_ks_rphi_vierbein_2v },
  { NULL }
};

static const vierbein_gradient_kern_list ks_rphi_vierbein_gradient_list[] = {
  { NULL },
  { eval_ks_rphi_vierbein_gradient_2v },
  { NULL }
};

static const vierbein_kern_list ks_vierbein_list[] = {
  { NULL },
  { NULL},
  { eval_ks_vierbein_3v  }
};

static const vierbein_gradient_kern_list ks_vierbein_gradient_list[] = {
  { NULL },
  { eval_ks_vierbein_gradient_3v },
  { NULL }
};

static const vierbein_kern_list rz_cylindrical_vierbein_list[] = {
  { NULL },
  { NULL },
  { eval_rz_cylindrical_vierbein_3v }
};

static const vierbein_gradient_kern_list rz_cylindrical_vierbein_gradient_list[] = {
  { NULL },
  { NULL },
  { eval_rz_cylindrical_vierbein_gradient_3v }
};

static evalf_t
choose_vierbein_kern(enum gkyl_triad_preset_geom_type type, int vdim)
{
  switch(type) {
    case GKYL_TRIAD_ANNULUS:
      return annulus_vierbein_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_CYLINDRICAL_RZ:
      return rz_cylindrical_vierbein_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_GR_KERR_SCHILD_RPHI:
      return ks_rphi_vierbein_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_GR_KERR_SCHILD_3V:
      return ks_vierbein_list[vdim-1].kernels[0];
      break;
    default:
      assert(false);
  }
}

static evalf_t
choose_vierbein_gradient_kern(enum gkyl_triad_preset_geom_type type, int vdim)
{
  switch(type) {
    case GKYL_TRIAD_ANNULUS:
      return annulus_vierbein_gradient_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_CYLINDRICAL_RZ:
      return rz_cylindrical_vierbein_gradient_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_GR_KERR_SCHILD_RPHI:
      return ks_rphi_vierbein_gradient_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_GR_KERR_SCHILD_3V:
      return ks_vierbein_gradient_list[vdim-1].kernels[0];
      break;
    default:
      assert(false);
  }
}

static evalf_t
choose_hamil_kern(enum gkyl_triad_preset_geom_type type, int vdim)
{
  switch(type) {
    case GKYL_TRIAD_GR_KERR_SCHILD_RPHI:
      return ks_rphi_hamil_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_GR_KERR_SCHILD_3V:
      return ks_hamil_3v_list[vdim-1].kernels[0];
      break;
    default:
      assert(false);
  }
}
