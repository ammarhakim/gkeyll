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
eval_ks_rphi_hamil_1x2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Grab the geometry
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;
  
  double r = xn[0];

  double pr_hat = xn[1];
  double pphi_hat = xn[2];

  // Build the Hamiltonian for Kerr-Schild in 2D (Assumes theta = pi/2)
  double rho_sq = r * r ;
  double H = (1.0/ sqrt(1 + 2*M*r/rho_sq)) * (
        sqrt(1 + pr_hat * pr_hat + pphi_hat * pphi_hat) 
        - (2*M*r/rho_sq) * pr_hat );

  fout[0] = H;
}

void
eval_ks_rphi_hamil_2x2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Grab the geometry
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;
  
  double r = xn[0];
  double phi = xn[1];

  double pr_hat = xn[2];
  double pphi_hat = xn[3];

  // Build the Hamiltonian for Kerr-Schild in 2D (Assumes theta = pi/2)
  double rho_sq = r * r ;
  double H = (1.0/ sqrt(1 + 2*M*r/rho_sq)) * (
        sqrt(1 + pr_hat * pr_hat + pphi_hat * pphi_hat) 
        - (2*M*r/rho_sq) * pr_hat );

  fout[0] = H;
}

void
eval_ks_r_hamil_1x3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Grab the geometry
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;
  
  double r = xn[0];

  double pr_hat = xn[1];
  double ptheta_hat = xn[2];
  double pphi_hat = xn[3];

  // Build the Hamiltonian for Kerr-Schild in 3D
  double rho_sq = r * r;
  double H = (1.0/ sqrt(1 + 2*M*r/rho_sq)) * (
        sqrt(1 + pr_hat * pr_hat + ptheta_hat * ptheta_hat + pphi_hat * pphi_hat) 
        - (2*M*r/rho_sq) * pr_hat );

  fout[0] = H;
}


void
eval_ks_rtheta_hamil_2x3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Grab the geometry
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;
  
  double r = xn[0];
  double theta = xn[1];

  double pr_hat = xn[2];
  double ptheta_hat = xn[3];
  double pphi_hat = xn[4];

  // Build the Hamiltonian for Kerr-Schild in 3D
  double rho_sq = r * r + a * a * cos(theta) * cos(theta);
  double H = (1.0/ sqrt(1 + 2*M*r/rho_sq)) * (
        sqrt(1 + pr_hat * pr_hat + ptheta_hat * ptheta_hat + pphi_hat * pphi_hat) 
        - (2*M*r/rho_sq) * pr_hat );

  fout[0] = H;
}

void
eval_ks_rphi_vierbein_2v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Parameters
  const struct gkyl_triad_geom_ctx *geom = ctx;
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

  const struct gkyl_triad_geom_ctx *geom = ctx;
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

  // Gradient w.r.t. phi: d(e_i^a)/dphi
  fout[4] = 0.0;
  fout[5] = 0.0;
  fout[6] = 0.0;
  fout[7] = 0.0;

}

void 
eval_ks_r_vierbein_3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  // Parameters
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;

  // Coordinates
  double r = xn[0];

  // Intermediate Variables
  double rho_sq = r * r;

  // Metric spatial covariant components
  double h_rr = ( 1.0 + 2.0 * M * r / rho_sq);
  double h_pp = ( rho_sq + a * a * ( 1.0 + 2.0 * M * r / rho_sq )  );
  double h_thth = rho_sq;
  double h_rp = - a * ( 1.0 + 2.0 * M * r / rho_sq );

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
eval_ks_r_vierbein_gradient_3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;

  double r = xn[0];

  // Intermediate Variables
  double rho_sq = r * r ;

  // Gradient w.r.t. r: d(e_i^a)/dr
  fout[0] = (M * ( - r * r)) / (pow(rho_sq, 1.5) * sqrt(rho_sq + 2.0 * M * r));
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = 0.0;
  fout[4] = r / sqrt(rho_sq);
  fout[5] = 0.0;
  fout[6] = -(M * a * ( - r * r)) / (pow(rho_sq, 1.5) * sqrt(rho_sq + 2.0 * M * r));
  fout[7] = 0.0;
  fout[8] = (r) / sqrt(rho_sq);

  // Gradient w.r.t. theta: d(e_i^a)/dtheta
  fout[9]  = 0.0;
  fout[10] = 0.0;
  fout[11] = 0.0;
  fout[12] = 0.0;
  fout[13] = 0.0;
  fout[14] = 0.0;
  fout[15] = 0.0;
  fout[16] = 0.0;
  fout[17] = 0.0;

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
eval_ks_rtheta_vierbein_3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  // Parameters
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;

  // Coordinates
  double r = xn[0];
  double theta = xn[1];

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
eval_ks_rtheta_vierbein_gradient_3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;

  double r = xn[0];
  double theta = xn[1];

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
eval_ks_vierbein_3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  // Parameters
  const struct gkyl_triad_geom_ctx *geom = ctx;
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
  const struct gkyl_triad_geom_ctx *geom = ctx;
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
eval_ks_hamil_3x3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{

  // Grab the geometry
  const struct gkyl_triad_geom_ctx *geom = ctx;
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

void 
eval_ks_vierbein_cart_1x3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  // Parameters
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;

  // Coordinates
  double x = xn[0];
  double y = 0.0;
  double z = 0.0;

  // Cartesian Kerr-Schild radius.
  double x2 = x*x;
  double y2 = y*y;
  double z2 = z*z;
  double a2 = a*a;
  double R2_cart = x2 + y2 + z2;
  double disc = (R2_cart - a2)*(R2_cart - a2) + 4.0*a2*z2;
  double r = sqrt(0.5*((R2_cart - a2) + sqrt(disc)));

  // Compact Kerr-Schild scalar quantities.
  double r2 = r*r;
  double r3 = r2*r;
  double r4 = r2*r2;
  double RHO2 = r4 + a2*z2;
  double SIG = r2 + a2;
  double HK = M*r3/RHO2;

  // Spatial part of the Kerr-Schild null vector.
  double L1 = (r*x + a*y)/SIG;
  double L2 = (r*y - a*x)/SIG;
  double L3 = z/r;

  // Spatial covariant metric h_ij = delta_ij + 2 HK L_i L_j.
  double h11 = 1.0 + 2.0*HK*L1*L1;
  double h12 =       2.0*HK*L1*L2;
  double h13 =       2.0*HK*L1*L3;
  double h22 = 1.0 + 2.0*HK*L2*L2;
  double h23 =       2.0*HK*L2*L3;
  double h33 = 1.0 + 2.0*HK*L3*L3;

  // Spatial contravariant metric h^ij = delta^ij - kappa L^i L^j.
  double kappa = 2.0*HK/(1.0 + 2.0*HK);
  double hc11 = 1.0 - kappa*L1*L1;
  double hc12 =     - kappa*L1*L2;
  double hc13 =     - kappa*L1*L3;
  double hc22 = 1.0 - kappa*L2*L2;
  double hc23 =     - kappa*L2*L3;
  double hc33 = 1.0 - kappa*L3*L3;

  // Cholesky-like factors for the contravariant tetrad spatial block.
  double B = sqrt(hc11);
  double D = hc12/B;
  double E = sqrt(hc22 - hc12*hc12/hc11);
  double G = hc13/B;
  double Ch = (hc23 - hc13*hc12/hc11)/E;
  double I = sqrt(hc33 - hc13*hc13/hc11 - (hc23 - hc13*hc12/hc11)*(hc23 - hc13*hc12/hc11)/(hc22 - hc12*hc12/hc11));

  // Spatial block of e^i_a from the MATLAB construction.
  // The minus signs are inherited from vierb_con_mualpha.
  double vc11 = -B;
  double vc12 =  0.0;
  double vc13 =  0.0;
  double vc21 = -D;
  double vc22 = -E;
  double vc23 =  0.0;
  double vc31 = -G;
  double vc32 = -Ch;
  double vc33 = -I;

  // Raw lowered spatial block from the 4D covariant tetrad is h_ij e^j_a.
  // Gkeyll convention in the existing routines uses the opposite spatial sign,
  // so define e_i^a = - h_ij e^j_a.
  double e11 = -(h11*vc11 + h12*vc21 + h13*vc31);
  double e12 = -(h11*vc12 + h12*vc22 + h13*vc32);
  double e13 = -(h11*vc13 + h12*vc23 + h13*vc33);

  double e21 = -(h12*vc11 + h22*vc21 + h23*vc31);
  double e22 = -(h12*vc12 + h22*vc22 + h23*vc32);
  double e23 = -(h12*vc13 + h22*vc23 + h23*vc33);

  double e31 = -(h13*vc11 + h23*vc21 + h33*vc31);
  double e32 = -(h13*vc12 + h23*vc22 + h33*vc32);
  double e33 = -(h13*vc13 + h23*vc23 + h33*vc33);

  // Vierbein: e_i^a = g_i . sigma^a
  fout[0] = e11;
  fout[1] = e12;
  fout[2] = e13;
  fout[3] = e21;
  fout[4] = e22;
  fout[5] = e23;
  fout[6] = e31;
  fout[7] = e32;
  fout[8] = e33;
}

void 
eval_ks_vierbein_cart_2x3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  // Parameters
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;

  // Coordinates
  double x = xn[0];
  double y = xn[1];
  double z = 0.0;

  // Cartesian Kerr-Schild radius.
  double x2 = x*x;
  double y2 = y*y;
  double z2 = z*z;
  double a2 = a*a;
  double R2_cart = x2 + y2 + z2;
  double disc = (R2_cart - a2)*(R2_cart - a2) + 4.0*a2*z2;
  double r = sqrt(0.5*((R2_cart - a2) + sqrt(disc)));

  // Compact Kerr-Schild scalar quantities.
  double r2 = r*r;
  double r3 = r2*r;
  double r4 = r2*r2;
  double RHO2 = r4 + a2*z2;
  double SIG = r2 + a2;
  double HK = M*r3/RHO2;

  // Spatial part of the Kerr-Schild null vector.
  double L1 = (r*x + a*y)/SIG;
  double L2 = (r*y - a*x)/SIG;
  double L3 = z/r;

  // Spatial covariant metric h_ij = delta_ij + 2 HK L_i L_j.
  double h11 = 1.0 + 2.0*HK*L1*L1;
  double h12 =       2.0*HK*L1*L2;
  double h13 =       2.0*HK*L1*L3;
  double h22 = 1.0 + 2.0*HK*L2*L2;
  double h23 =       2.0*HK*L2*L3;
  double h33 = 1.0 + 2.0*HK*L3*L3;

  // Spatial contravariant metric h^ij = delta^ij - kappa L^i L^j.
  double kappa = 2.0*HK/(1.0 + 2.0*HK);
  double hc11 = 1.0 - kappa*L1*L1;
  double hc12 =     - kappa*L1*L2;
  double hc13 =     - kappa*L1*L3;
  double hc22 = 1.0 - kappa*L2*L2;
  double hc23 =     - kappa*L2*L3;
  double hc33 = 1.0 - kappa*L3*L3;

  // Cholesky-like factors for the contravariant tetrad spatial block.
  double B = sqrt(hc11);
  double D = hc12/B;
  double E = sqrt(hc22 - hc12*hc12/hc11);
  double G = hc13/B;
  double Ch = (hc23 - hc13*hc12/hc11)/E;
  double I = sqrt(hc33 - hc13*hc13/hc11 - (hc23 - hc13*hc12/hc11)*(hc23 - hc13*hc12/hc11)/(hc22 - hc12*hc12/hc11));

  // Spatial block of e^i_a from the MATLAB construction.
  // The minus signs are inherited from vierb_con_mualpha.
  double vc11 = -B;
  double vc12 =  0.0;
  double vc13 =  0.0;
  double vc21 = -D;
  double vc22 = -E;
  double vc23 =  0.0;
  double vc31 = -G;
  double vc32 = -Ch;
  double vc33 = -I;

  // Raw lowered spatial block from the 4D covariant tetrad is h_ij e^j_a.
  // Gkeyll convention in the existing routines uses the opposite spatial sign,
  // so define e_i^a = - h_ij e^j_a.
  double e11 = -(h11*vc11 + h12*vc21 + h13*vc31);
  double e12 = -(h11*vc12 + h12*vc22 + h13*vc32);
  double e13 = -(h11*vc13 + h12*vc23 + h13*vc33);

  double e21 = -(h12*vc11 + h22*vc21 + h23*vc31);
  double e22 = -(h12*vc12 + h22*vc22 + h23*vc32);
  double e23 = -(h12*vc13 + h22*vc23 + h23*vc33);

  double e31 = -(h13*vc11 + h23*vc21 + h33*vc31);
  double e32 = -(h13*vc12 + h23*vc22 + h33*vc32);
  double e33 = -(h13*vc13 + h23*vc23 + h33*vc33);

  // Vierbein: e_i^a = g_i . sigma^a
  fout[0] = e11;
  fout[1] = e12;
  fout[2] = e13;
  fout[3] = e21;
  fout[4] = e22;
  fout[5] = e23;
  fout[6] = e31;
  fout[7] = e32;
  fout[8] = e33;
}

void 
eval_ks_vierbein_cart_3x3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  // Parameters
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;

  // Coordinates
  double x = xn[0];
  double y = xn[1];
  double z = xn[2];

  // Cartesian Kerr-Schild radius.
  double x2 = x*x;
  double y2 = y*y;
  double z2 = z*z;
  double a2 = a*a;
  double R2_cart = x2 + y2 + z2;
  double disc = (R2_cart - a2)*(R2_cart - a2) + 4.0*a2*z2;
  double r = sqrt(0.5*((R2_cart - a2) + sqrt(disc)));

  // Compact Kerr-Schild scalar quantities.
  double r2 = r*r;
  double r3 = r2*r;
  double r4 = r2*r2;
  double RHO2 = r4 + a2*z2;
  double SIG = r2 + a2;
  double HK = M*r3/RHO2;

  // Spatial part of the Kerr-Schild null vector.
  double L1 = (r*x + a*y)/SIG;
  double L2 = (r*y - a*x)/SIG;
  double L3 = z/r;

  // Spatial covariant metric h_ij = delta_ij + 2 HK L_i L_j.
  double h11 = 1.0 + 2.0*HK*L1*L1;
  double h12 =       2.0*HK*L1*L2;
  double h13 =       2.0*HK*L1*L3;
  double h22 = 1.0 + 2.0*HK*L2*L2;
  double h23 =       2.0*HK*L2*L3;
  double h33 = 1.0 + 2.0*HK*L3*L3;

  // Spatial contravariant metric h^ij = delta^ij - kappa L^i L^j.
  double kappa = 2.0*HK/(1.0 + 2.0*HK);
  double hc11 = 1.0 - kappa*L1*L1;
  double hc12 =     - kappa*L1*L2;
  double hc13 =     - kappa*L1*L3;
  double hc22 = 1.0 - kappa*L2*L2;
  double hc23 =     - kappa*L2*L3;
  double hc33 = 1.0 - kappa*L3*L3;

  // Cholesky-like factors for the contravariant tetrad spatial block.
  double B = sqrt(hc11);
  double D = hc12/B;
  double E = sqrt(hc22 - hc12*hc12/hc11);
  double G = hc13/B;
  double Ch = (hc23 - hc13*hc12/hc11)/E;
  double I = sqrt(hc33 - hc13*hc13/hc11 - (hc23 - hc13*hc12/hc11)*(hc23 - hc13*hc12/hc11)/(hc22 - hc12*hc12/hc11));

  // Spatial block of e^i_a from the MATLAB construction.
  // The minus signs are inherited from vierb_con_mualpha.
  double vc11 = -B;
  double vc12 =  0.0;
  double vc13 =  0.0;
  double vc21 = -D;
  double vc22 = -E;
  double vc23 =  0.0;
  double vc31 = -G;
  double vc32 = -Ch;
  double vc33 = -I;

  // Raw lowered spatial block from the 4D covariant tetrad is h_ij e^j_a.
  // Gkeyll convention in the existing routines uses the opposite spatial sign,
  // so define e_i^a = - h_ij e^j_a.
  double e11 = -(h11*vc11 + h12*vc21 + h13*vc31);
  double e12 = -(h11*vc12 + h12*vc22 + h13*vc32);
  double e13 = -(h11*vc13 + h12*vc23 + h13*vc33);

  double e21 = -(h12*vc11 + h22*vc21 + h23*vc31);
  double e22 = -(h12*vc12 + h22*vc22 + h23*vc32);
  double e23 = -(h12*vc13 + h22*vc23 + h23*vc33);

  double e31 = -(h13*vc11 + h23*vc21 + h33*vc31);
  double e32 = -(h13*vc12 + h23*vc22 + h33*vc32);
  double e33 = -(h13*vc13 + h23*vc23 + h33*vc33);

  // Vierbein: e_i^a = g_i . sigma^a
  fout[0] = e11;
  fout[1] = e12;
  fout[2] = e13;
  fout[3] = e21;
  fout[4] = e22;
  fout[5] = e23;
  fout[6] = e31;
  fout[7] = e32;
  fout[8] = e33;
}

void
eval_ks_hamil_cart_1x3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Grab the geometry
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;


  // y and z assumed to be zero for this case, so only x and the momenta are used as inputs.
  double x = xn[0];
  double y = 0.0; 
  double z = 0.0;

  double px_hat = xn[1];
  double py_hat = xn[2];
  double pz_hat = xn[3];

  // Cartesian Kerr-Schild radius.
  double x2 = x*x;
  double y2 = y*y;
  double z2 = z*z;
  double a2 = a*a;
  double R2_cart = x2 + y2 + z2;
  double disc = (R2_cart - a2)*(R2_cart - a2) + 4.0*a2*z2;
  double r = sqrt(0.5*((R2_cart - a2) + sqrt(disc)));

  // Compact Kerr-Schild scalar quantities.
  double r2 = r*r;
  double r3 = r2*r;
  double r4 = r2*r2;
  double RHO2 = r4 + a2*z2;
  double SIG = r2 + a2;
  double HK = M*r3/RHO2;

  // Spatial part of the Kerr-Schild null vector.
  double L1 = (r*x + a*y)/SIG;
  double L2 = (r*y - a*x)/SIG;
  double L3 = z/r;

  // ADM lapse.
  double alpha = 1.0/sqrt(1.0 + 2.0*HK);

  // Spatial covariant metric h_ij = delta_ij + 2 HK L_i L_j.
  double h11 = 1.0 + 2.0*HK*L1*L1;
  double h12 =       2.0*HK*L1*L2;
  double h13 =       2.0*HK*L1*L3;
  double h22 = 1.0 + 2.0*HK*L2*L2;
  double h23 =       2.0*HK*L2*L3;
  double h33 = 1.0 + 2.0*HK*L3*L3;

  // Spatial contravariant metric h^ij = delta^ij - kappa L^i L^j.
  double kappa = 2.0*HK/(1.0 + 2.0*HK);
  double hc11 = 1.0 - kappa*L1*L1;
  double hc12 =     - kappa*L1*L2;
  double hc13 =     - kappa*L1*L3;
  double hc22 = 1.0 - kappa*L2*L2;
  double hc23 =     - kappa*L2*L3;
  double hc33 = 1.0 - kappa*L3*L3;

  // Cholesky-like factors for the contravariant tetrad spatial block.
  double B = sqrt(hc11);
  double D = hc12/B;
  double E = sqrt(hc22 - hc12*hc12/hc11);
  double G = hc13/B;
  double Ch = (hc23 - hc13*hc12/hc11)/E;
  double I = sqrt(hc33 - hc13*hc13/hc11 - (hc23 - hc13*hc12/hc11)*(hc23 - hc13*hc12/hc11)/(hc22 - hc12*hc12/hc11));

  // Spatial block of e^i_a from the MATLAB construction.
  double vc11 = -B;
  double vc12 =  0.0;
  double vc13 =  0.0;
  double vc21 = -D;
  double vc22 = -E;
  double vc23 =  0.0;
  double vc31 = -G;
  double vc32 = -Ch;
  double vc33 = -I;

  // Gkeyll-sign lowered spatial triad e_i^a = -h_ij e^j_a.
  double e11 = -(h11*vc11 + h12*vc21 + h13*vc31);
  double e12 = -(h11*vc12 + h12*vc22 + h13*vc32);
  double e13 = -(h11*vc13 + h12*vc23 + h13*vc33);

  double e21 = -(h12*vc11 + h22*vc21 + h23*vc31);
  double e22 = -(h12*vc12 + h22*vc22 + h23*vc32);
  double e23 = -(h12*vc13 + h22*vc23 + h23*vc33);

  double e31 = -(h13*vc11 + h23*vc21 + h33*vc31);
  double e32 = -(h13*vc12 + h23*vc22 + h33*vc32);
  double e33 = -(h13*vc13 + h23*vc23 + h33*vc33);

  // Shift vector beta^i.
  double beta1 = kappa*L1;
  double beta2 = kappa*L2;
  double beta3 = kappa*L3;

  // Local-frame shift projection b_a = beta^i e_i^a / alpha.
  // The Hamiltonian has the same structure as the existing spherical KS routines:
  // H = alpha * (gamma - b_a p_hat_a).
  double b1 = (beta1*e11 + beta2*e21 + beta3*e31)/alpha;
  double b2 = (beta1*e12 + beta2*e22 + beta3*e32)/alpha;
  double b3 = (beta1*e13 + beta2*e23 + beta3*e33)/alpha;

  double gamma = sqrt(1.0 + px_hat*px_hat + py_hat*py_hat + pz_hat*pz_hat);
  double H = alpha*(gamma - b1*px_hat - b2*py_hat - b3*pz_hat);

  fout[0] = H;
}

void
eval_ks_hamil_cart_2x3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Grab the geometry
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;


  // y and z assumed to be zero for this case, so only x, y and the momenta are used as inputs.
  double x = xn[0];
  double y = xn[1];
  double z = 0.0;

  double px_hat = xn[2];
  double py_hat = xn[3];
  double pz_hat = xn[4];

  // Cartesian Kerr-Schild radius.
  double x2 = x*x;
  double y2 = y*y;
  double z2 = z*z;
  double a2 = a*a;
  double R2_cart = x2 + y2 + z2;
  double disc = (R2_cart - a2)*(R2_cart - a2) + 4.0*a2*z2;
  double r = sqrt(0.5*((R2_cart - a2) + sqrt(disc)));

  // Compact Kerr-Schild scalar quantities.
  double r2 = r*r;
  double r3 = r2*r;
  double r4 = r2*r2;
  double RHO2 = r4 + a2*z2;
  double SIG = r2 + a2;
  double HK = M*r3/RHO2;

  // Spatial part of the Kerr-Schild null vector.
  double L1 = (r*x + a*y)/SIG;
  double L2 = (r*y - a*x)/SIG;
  double L3 = z/r;

  // ADM lapse.
  double alpha = 1.0/sqrt(1.0 + 2.0*HK);

  // Spatial covariant metric h_ij = delta_ij + 2 HK L_i L_j.
  double h11 = 1.0 + 2.0*HK*L1*L1;
  double h12 =       2.0*HK*L1*L2;
  double h13 =       2.0*HK*L1*L3;
  double h22 = 1.0 + 2.0*HK*L2*L2;
  double h23 =       2.0*HK*L2*L3;
  double h33 = 1.0 + 2.0*HK*L3*L3;

  // Spatial contravariant metric h^ij = delta^ij - kappa L^i L^j.
  double kappa = 2.0*HK/(1.0 + 2.0*HK);
  double hc11 = 1.0 - kappa*L1*L1;
  double hc12 =     - kappa*L1*L2;
  double hc13 =     - kappa*L1*L3;
  double hc22 = 1.0 - kappa*L2*L2;
  double hc23 =     - kappa*L2*L3;
  double hc33 = 1.0 - kappa*L3*L3;

  // Cholesky-like factors for the contravariant tetrad spatial block.
  double B = sqrt(hc11);
  double D = hc12/B;
  double E = sqrt(hc22 - hc12*hc12/hc11);
  double G = hc13/B;
  double Ch = (hc23 - hc13*hc12/hc11)/E;
  double I = sqrt(hc33 - hc13*hc13/hc11 - (hc23 - hc13*hc12/hc11)*(hc23 - hc13*hc12/hc11)/(hc22 - hc12*hc12/hc11));

  // Spatial block of e^i_a from the MATLAB construction.
  double vc11 = -B;
  double vc12 =  0.0;
  double vc13 =  0.0;
  double vc21 = -D;
  double vc22 = -E;
  double vc23 =  0.0;
  double vc31 = -G;
  double vc32 = -Ch;
  double vc33 = -I;

  // Gkeyll-sign lowered spatial triad e_i^a = -h_ij e^j_a.
  double e11 = -(h11*vc11 + h12*vc21 + h13*vc31);
  double e12 = -(h11*vc12 + h12*vc22 + h13*vc32);
  double e13 = -(h11*vc13 + h12*vc23 + h13*vc33);

  double e21 = -(h12*vc11 + h22*vc21 + h23*vc31);
  double e22 = -(h12*vc12 + h22*vc22 + h23*vc32);
  double e23 = -(h12*vc13 + h22*vc23 + h23*vc33);

  double e31 = -(h13*vc11 + h23*vc21 + h33*vc31);
  double e32 = -(h13*vc12 + h23*vc22 + h33*vc32);
  double e33 = -(h13*vc13 + h23*vc23 + h33*vc33);

  // Shift vector beta^i.
  double beta1 = kappa*L1;
  double beta2 = kappa*L2;
  double beta3 = kappa*L3;

  // Local-frame shift projection b_a = beta^i e_i^a / alpha.
  // The Hamiltonian has the same structure as the existing spherical KS routines:
  // H = alpha * (gamma - b_a p_hat_a).
  double b1 = (beta1*e11 + beta2*e21 + beta3*e31)/alpha;
  double b2 = (beta1*e12 + beta2*e22 + beta3*e32)/alpha;
  double b3 = (beta1*e13 + beta2*e23 + beta3*e33)/alpha;

  double gamma = sqrt(1.0 + px_hat*px_hat + py_hat*py_hat + pz_hat*pz_hat);
  double H = alpha*(gamma - b1*px_hat - b2*py_hat - b3*pz_hat);

  fout[0] = H;
}

void
eval_ks_hamil_cart_3x3v(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void* ctx)
{
  // Grab the geometry
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;


  double x = xn[0];
  double y = xn[1];
  double z = xn[2];

  double px_hat = xn[3];
  double py_hat = xn[4];
  double pz_hat = xn[5];

  // Cartesian Kerr-Schild radius.
  double x2 = x*x;
  double y2 = y*y;
  double z2 = z*z;
  double a2 = a*a;
  double R2_cart = x2 + y2 + z2;
  double disc = (R2_cart - a2)*(R2_cart - a2) + 4.0*a2*z2;
  double r = sqrt(0.5*((R2_cart - a2) + sqrt(disc)));

  // Compact Kerr-Schild scalar quantities.
  double r2 = r*r;
  double r3 = r2*r;
  double r4 = r2*r2;
  double RHO2 = r4 + a2*z2;
  double SIG = r2 + a2;
  double HK = M*r3/RHO2;

  // Spatial part of the Kerr-Schild null vector.
  double L1 = (r*x + a*y)/SIG;
  double L2 = (r*y - a*x)/SIG;
  double L3 = z/r;

  // ADM lapse.
  double alpha = 1.0/sqrt(1.0 + 2.0*HK);

  // Spatial covariant metric h_ij = delta_ij + 2 HK L_i L_j.
  double h11 = 1.0 + 2.0*HK*L1*L1;
  double h12 =       2.0*HK*L1*L2;
  double h13 =       2.0*HK*L1*L3;
  double h22 = 1.0 + 2.0*HK*L2*L2;
  double h23 =       2.0*HK*L2*L3;
  double h33 = 1.0 + 2.0*HK*L3*L3;

  // Spatial contravariant metric h^ij = delta^ij - kappa L^i L^j.
  double kappa = 2.0*HK/(1.0 + 2.0*HK);
  double hc11 = 1.0 - kappa*L1*L1;
  double hc12 =     - kappa*L1*L2;
  double hc13 =     - kappa*L1*L3;
  double hc22 = 1.0 - kappa*L2*L2;
  double hc23 =     - kappa*L2*L3;
  double hc33 = 1.0 - kappa*L3*L3;

  // Cholesky-like factors for the contravariant tetrad spatial block.
  double B = sqrt(hc11);
  double D = hc12/B;
  double E = sqrt(hc22 - hc12*hc12/hc11);
  double G = hc13/B;
  double Ch = (hc23 - hc13*hc12/hc11)/E;
  double I = sqrt(hc33 - hc13*hc13/hc11 - (hc23 - hc13*hc12/hc11)*(hc23 - hc13*hc12/hc11)/(hc22 - hc12*hc12/hc11));

  // Spatial block of e^i_a from the MATLAB construction.
  double vc11 = -B;
  double vc12 =  0.0;
  double vc13 =  0.0;
  double vc21 = -D;
  double vc22 = -E;
  double vc23 =  0.0;
  double vc31 = -G;
  double vc32 = -Ch;
  double vc33 = -I;

  // Gkeyll-sign lowered spatial triad e_i^a = -h_ij e^j_a.
  double e11 = -(h11*vc11 + h12*vc21 + h13*vc31);
  double e12 = -(h11*vc12 + h12*vc22 + h13*vc32);
  double e13 = -(h11*vc13 + h12*vc23 + h13*vc33);

  double e21 = -(h12*vc11 + h22*vc21 + h23*vc31);
  double e22 = -(h12*vc12 + h22*vc22 + h23*vc32);
  double e23 = -(h12*vc13 + h22*vc23 + h23*vc33);

  double e31 = -(h13*vc11 + h23*vc21 + h33*vc31);
  double e32 = -(h13*vc12 + h23*vc22 + h33*vc32);
  double e33 = -(h13*vc13 + h23*vc23 + h33*vc33);

  // Shift vector beta^i.
  double beta1 = kappa*L1;
  double beta2 = kappa*L2;
  double beta3 = kappa*L3;

  // Local-frame shift projection b_a = beta^i e_i^a / alpha.
  // The Hamiltonian has the same structure as the existing spherical KS routines:
  // H = alpha * (gamma - b_a p_hat_a).
  double b1 = (beta1*e11 + beta2*e21 + beta3*e31)/alpha;
  double b2 = (beta1*e12 + beta2*e22 + beta3*e32)/alpha;
  double b3 = (beta1*e13 + beta2*e23 + beta3*e33)/alpha;

  double gamma = sqrt(1.0 + px_hat*px_hat + py_hat*py_hat + pz_hat*pz_hat);
  double H = alpha*(gamma - b1*px_hat - b2*py_hat - b3*pz_hat);

  fout[0] = H;
}


void
eval_ks_vierbein_gradient_cart_1x3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;

  double x = xn[0];
  double y = 0.0;
  double z = 0.0;


  // Cartesian Kerr-Schild radius.
  double x2 = x*x;
  double y2 = y*y;
  double z2 = z*z;
  double a2 = a*a;
  double R2_cart = x2 + y2 + z2;
  double disc = (R2_cart - a2)*(R2_cart - a2) + 4.0*a2*z2;
  double r = sqrt(0.5*((R2_cart - a2) + sqrt(disc)));

  // Common powers and compact Kerr-Schild scalars.
  double r2 = r*r;
  double r3 = r2*r;
  double r4 = r2*r2;
  double RHO2 = r4 + a2*z2;
  double SIG  = r2 + a2;

  double HK = M*r3/RHO2;
  double L1 = (r*x + a*y)/SIG;
  double L2 = (r*y - a*x)/SIG;
  double L3 = z/r;

  // Implicit Cartesian derivatives of r.
  double rx = x*r3/RHO2;
  double ry = y*r3/RHO2;
  double rz = z*r*SIG/RHO2;

  // Quotient-rule derivatives of RHO2 and SIG.
  double RHO2_x = 4.0*r3*rx;
  double RHO2_y = 4.0*r3*ry;
  double RHO2_z = 4.0*r3*rz + 2.0*a2*z;

  double SIG_x = 2.0*r*rx;
  double SIG_y = 2.0*r*ry;
  double SIG_z = 2.0*r*rz;

  // Numerators for L1 and L2.
  double N1 = r*x + a*y;
  double N2 = r*y - a*x;

  double N1_x = rx*x + r;
  double N1_y = ry*x + a;
  double N1_z = rz*x;

  double N2_x = rx*y - a;
  double N2_y = ry*y + r;
  double N2_z = rz*y;

  // Compact scalar derivatives.
  double HK_x = M*(3.0*r2*rx*RHO2 - r3*RHO2_x)/(RHO2*RHO2);
  double HK_y = M*(3.0*r2*ry*RHO2 - r3*RHO2_y)/(RHO2*RHO2);
  double HK_z = M*(3.0*r2*rz*RHO2 - r3*RHO2_z)/(RHO2*RHO2);

  double L1_x = (N1_x*SIG - N1*SIG_x)/(SIG*SIG);
  double L1_y = (N1_y*SIG - N1*SIG_y)/(SIG*SIG);
  double L1_z = (N1_z*SIG - N1*SIG_z)/(SIG*SIG);

  double L2_x = (N2_x*SIG - N2*SIG_x)/(SIG*SIG);
  double L2_y = (N2_y*SIG - N2*SIG_y)/(SIG*SIG);
  double L2_z = (N2_z*SIG - N2*SIG_z)/(SIG*SIG);

  double L3_x = -z*rx/r2;
  double L3_y = -z*ry/r2;
  double L3_z = (r - z*rz)/r2;

  // Pack compact variables and their Cartesian derivatives.
  double L[3] = { L1, L2, L3 };
  double dL[3][3] = {
    { L1_x, L2_x, L3_x },
    { L1_y, L2_y, L3_y },
    { L1_z, L2_z, L3_z }
  };
  double dHK[3] = { HK_x, HK_y, HK_z };

  // Compact inverse spatial metric h^{ij} = delta^{ij} - kappa L^i L^j.
  double s = 1.0 + 2.0*HK;
  double kappa = 2.0*HK/s;

  double hcon[3][3];
  double hcov[3][3];
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) {
      double delta = (i == j) ? 1.0 : 0.0;
      hcon[i][j] = delta - kappa*L[i]*L[j];
      hcov[i][j] = delta + 2.0*HK*L[i]*L[j];
    }
  }

  // Cholesky-like factors for the contravariant spatial vierbein.
  double h11 = hcon[0][0];
  double h12 = hcon[0][1];
  double h13 = hcon[0][2];
  double h22 = hcon[1][1];
  double h23 = hcon[1][2];
  double h33 = hcon[2][2];

  double B = sqrt(h11);
  double D = h12/B;
  double Q = h22 - h12*h12/h11;
  double E = sqrt(Q);
  double G = h13/B;
  double Nch = h23 - h13*h12/h11;
  double Ch = Nch/E;
  double Iarg = h33 - h13*h13/h11 - Nch*Nch/Q;
  double I = sqrt(Iarg);

  // Spatial block of e^i_a. This matches the MATLAB sign convention.
  double V[3][3] = {
    { -B,  0.0, 0.0 },
    { -D, -E,  0.0 },
    { -G, -Ch, -I }
  };

  // Forward-mode derivative propagation for q = x,y,z.
  for (int q=0; q<3; ++q) {
    double dkappa = 2.0*dHK[q]/(s*s);

    double dhcon[3][3];
    double dhcov[3][3];
    for (int i=0; i<3; ++i) {
      for (int j=0; j<3; ++j) {
        dhcon[i][j] = -dkappa*L[i]*L[j]
                     - kappa*(dL[q][i]*L[j] + L[i]*dL[q][j]);
        dhcov[i][j] = 2.0*dHK[q]*L[i]*L[j]
                    + 2.0*HK*(dL[q][i]*L[j] + L[i]*dL[q][j]);
      }
    }

    double dh11 = dhcon[0][0];
    double dh12 = dhcon[0][1];
    double dh13 = dhcon[0][2];
    double dh22 = dhcon[1][1];
    double dh23 = dhcon[1][2];
    double dh33 = dhcon[2][2];

    double dB = dh11/(2.0*B);
    double dD = dh12/B - h12*dB/(B*B);
    double dQ = dh22 - (2.0*h12*dh12*h11 - h12*h12*dh11)/(h11*h11);
    double dE = dQ/(2.0*E);
    double dG = dh13/B - h13*dB/(B*B);
    double dNch = dh23 - (dh13*h12 + h13*dh12)/h11 + h13*h12*dh11/(h11*h11);
    double dCh = dNch/E - Nch*dE/(E*E);
    double dIarg = dh33
      - (2.0*h13*dh13*h11 - h13*h13*dh11)/(h11*h11)
      - (2.0*Nch*dNch*Q - Nch*Nch*dQ)/(Q*Q);
    double dI = dIarg/(2.0*I);

    double dV[3][3] = {
      { -dB,  0.0, 0.0 },
      { -dD, -dE,  0.0 },
      { -dG, -dCh, -dI }
    };

    // e_cov_spatial = -hcov * V, so d e = -(dhcov * V + hcov * dV).
    for (int i=0; i<3; ++i) {
      for (int alpha=0; alpha<3; ++alpha) {
        double de = 0.0;
        for (int j=0; j<3; ++j) {
          de += dhcov[i][j]*V[j][alpha] + hcov[i][j]*dV[j][alpha];
        }
        fout[9*q + 3*i + alpha] = -de;
      }
    }
  }
}


void
eval_ks_vierbein_gradient_cart_2x3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;

  double x = xn[0];
  double y = xn[1];
  double z = 0.0;


  // Cartesian Kerr-Schild radius.
  double x2 = x*x;
  double y2 = y*y;
  double z2 = z*z;
  double a2 = a*a;
  double R2_cart = x2 + y2 + z2;
  double disc = (R2_cart - a2)*(R2_cart - a2) + 4.0*a2*z2;
  double r = sqrt(0.5*((R2_cart - a2) + sqrt(disc)));

  // Common powers and compact Kerr-Schild scalars.
  double r2 = r*r;
  double r3 = r2*r;
  double r4 = r2*r2;
  double RHO2 = r4 + a2*z2;
  double SIG  = r2 + a2;

  double HK = M*r3/RHO2;
  double L1 = (r*x + a*y)/SIG;
  double L2 = (r*y - a*x)/SIG;
  double L3 = z/r;

  // Implicit Cartesian derivatives of r.
  double rx = x*r3/RHO2;
  double ry = y*r3/RHO2;
  double rz = z*r*SIG/RHO2;

  // Quotient-rule derivatives of RHO2 and SIG.
  double RHO2_x = 4.0*r3*rx;
  double RHO2_y = 4.0*r3*ry;
  double RHO2_z = 4.0*r3*rz + 2.0*a2*z;

  double SIG_x = 2.0*r*rx;
  double SIG_y = 2.0*r*ry;
  double SIG_z = 2.0*r*rz;

  // Numerators for L1 and L2.
  double N1 = r*x + a*y;
  double N2 = r*y - a*x;

  double N1_x = rx*x + r;
  double N1_y = ry*x + a;
  double N1_z = rz*x;

  double N2_x = rx*y - a;
  double N2_y = ry*y + r;
  double N2_z = rz*y;

  // Compact scalar derivatives.
  double HK_x = M*(3.0*r2*rx*RHO2 - r3*RHO2_x)/(RHO2*RHO2);
  double HK_y = M*(3.0*r2*ry*RHO2 - r3*RHO2_y)/(RHO2*RHO2);
  double HK_z = M*(3.0*r2*rz*RHO2 - r3*RHO2_z)/(RHO2*RHO2);

  double L1_x = (N1_x*SIG - N1*SIG_x)/(SIG*SIG);
  double L1_y = (N1_y*SIG - N1*SIG_y)/(SIG*SIG);
  double L1_z = (N1_z*SIG - N1*SIG_z)/(SIG*SIG);

  double L2_x = (N2_x*SIG - N2*SIG_x)/(SIG*SIG);
  double L2_y = (N2_y*SIG - N2*SIG_y)/(SIG*SIG);
  double L2_z = (N2_z*SIG - N2*SIG_z)/(SIG*SIG);

  double L3_x = -z*rx/r2;
  double L3_y = -z*ry/r2;
  double L3_z = (r - z*rz)/r2;

  // Pack compact variables and their Cartesian derivatives.
  double L[3] = { L1, L2, L3 };
  double dL[3][3] = {
    { L1_x, L2_x, L3_x },
    { L1_y, L2_y, L3_y },
    { L1_z, L2_z, L3_z }
  };
  double dHK[3] = { HK_x, HK_y, HK_z };

  // Compact inverse spatial metric h^{ij} = delta^{ij} - kappa L^i L^j.
  double s = 1.0 + 2.0*HK;
  double kappa = 2.0*HK/s;

  double hcon[3][3];
  double hcov[3][3];
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) {
      double delta = (i == j) ? 1.0 : 0.0;
      hcon[i][j] = delta - kappa*L[i]*L[j];
      hcov[i][j] = delta + 2.0*HK*L[i]*L[j];
    }
  }

  // Cholesky-like factors for the contravariant spatial vierbein.
  double h11 = hcon[0][0];
  double h12 = hcon[0][1];
  double h13 = hcon[0][2];
  double h22 = hcon[1][1];
  double h23 = hcon[1][2];
  double h33 = hcon[2][2];

  double B = sqrt(h11);
  double D = h12/B;
  double Q = h22 - h12*h12/h11;
  double E = sqrt(Q);
  double G = h13/B;
  double Nch = h23 - h13*h12/h11;
  double Ch = Nch/E;
  double Iarg = h33 - h13*h13/h11 - Nch*Nch/Q;
  double I = sqrt(Iarg);

  // Spatial block of e^i_a. This matches the MATLAB sign convention.
  double V[3][3] = {
    { -B,  0.0, 0.0 },
    { -D, -E,  0.0 },
    { -G, -Ch, -I }
  };

  // Forward-mode derivative propagation for q = x,y,z.
  for (int q=0; q<3; ++q) {
    double dkappa = 2.0*dHK[q]/(s*s);

    double dhcon[3][3];
    double dhcov[3][3];
    for (int i=0; i<3; ++i) {
      for (int j=0; j<3; ++j) {
        dhcon[i][j] = -dkappa*L[i]*L[j]
                     - kappa*(dL[q][i]*L[j] + L[i]*dL[q][j]);
        dhcov[i][j] = 2.0*dHK[q]*L[i]*L[j]
                    + 2.0*HK*(dL[q][i]*L[j] + L[i]*dL[q][j]);
      }
    }

    double dh11 = dhcon[0][0];
    double dh12 = dhcon[0][1];
    double dh13 = dhcon[0][2];
    double dh22 = dhcon[1][1];
    double dh23 = dhcon[1][2];
    double dh33 = dhcon[2][2];

    double dB = dh11/(2.0*B);
    double dD = dh12/B - h12*dB/(B*B);
    double dQ = dh22 - (2.0*h12*dh12*h11 - h12*h12*dh11)/(h11*h11);
    double dE = dQ/(2.0*E);
    double dG = dh13/B - h13*dB/(B*B);
    double dNch = dh23 - (dh13*h12 + h13*dh12)/h11 + h13*h12*dh11/(h11*h11);
    double dCh = dNch/E - Nch*dE/(E*E);
    double dIarg = dh33
      - (2.0*h13*dh13*h11 - h13*h13*dh11)/(h11*h11)
      - (2.0*Nch*dNch*Q - Nch*Nch*dQ)/(Q*Q);
    double dI = dIarg/(2.0*I);

    double dV[3][3] = {
      { -dB,  0.0, 0.0 },
      { -dD, -dE,  0.0 },
      { -dG, -dCh, -dI }
    };

    // e_cov_spatial = -hcov * V, so d e = -(dhcov * V + hcov * dV).
    for (int i=0; i<3; ++i) {
      for (int alpha=0; alpha<3; ++alpha) {
        double de = 0.0;
        for (int j=0; j<3; ++j) {
          de += dhcov[i][j]*V[j][alpha] + hcov[i][j]*dV[j][alpha];
        }
        fout[9*q + 3*i + alpha] = -de;
      }
    }
  }
}

void
eval_ks_vierbein_gradient_cart_3x3v(double t, const double *xn, double* restrict fout, void *ctx)
{
  const struct gkyl_triad_geom_ctx *geom = ctx;
  double a = geom->spin_bh;
  double M = geom->mass_bh;

  double x = xn[0];
  double y = xn[1];
  double z = xn[2];


  // Cartesian Kerr-Schild radius.
  double x2 = x*x;
  double y2 = y*y;
  double z2 = z*z;
  double a2 = a*a;
  double R2_cart = x2 + y2 + z2;
  double disc = (R2_cart - a2)*(R2_cart - a2) + 4.0*a2*z2;
  double r = sqrt(0.5*((R2_cart - a2) + sqrt(disc)));

  // Common powers and compact Kerr-Schild scalars.
  double r2 = r*r;
  double r3 = r2*r;
  double r4 = r2*r2;
  double RHO2 = r4 + a2*z2;
  double SIG  = r2 + a2;

  double HK = M*r3/RHO2;
  double L1 = (r*x + a*y)/SIG;
  double L2 = (r*y - a*x)/SIG;
  double L3 = z/r;

  // Implicit Cartesian derivatives of r.
  double rx = x*r3/RHO2;
  double ry = y*r3/RHO2;
  double rz = z*r*SIG/RHO2;

  // Quotient-rule derivatives of RHO2 and SIG.
  double RHO2_x = 4.0*r3*rx;
  double RHO2_y = 4.0*r3*ry;
  double RHO2_z = 4.0*r3*rz + 2.0*a2*z;

  double SIG_x = 2.0*r*rx;
  double SIG_y = 2.0*r*ry;
  double SIG_z = 2.0*r*rz;

  // Numerators for L1 and L2.
  double N1 = r*x + a*y;
  double N2 = r*y - a*x;

  double N1_x = rx*x + r;
  double N1_y = ry*x + a;
  double N1_z = rz*x;

  double N2_x = rx*y - a;
  double N2_y = ry*y + r;
  double N2_z = rz*y;

  // Compact scalar derivatives.
  double HK_x = M*(3.0*r2*rx*RHO2 - r3*RHO2_x)/(RHO2*RHO2);
  double HK_y = M*(3.0*r2*ry*RHO2 - r3*RHO2_y)/(RHO2*RHO2);
  double HK_z = M*(3.0*r2*rz*RHO2 - r3*RHO2_z)/(RHO2*RHO2);

  double L1_x = (N1_x*SIG - N1*SIG_x)/(SIG*SIG);
  double L1_y = (N1_y*SIG - N1*SIG_y)/(SIG*SIG);
  double L1_z = (N1_z*SIG - N1*SIG_z)/(SIG*SIG);

  double L2_x = (N2_x*SIG - N2*SIG_x)/(SIG*SIG);
  double L2_y = (N2_y*SIG - N2*SIG_y)/(SIG*SIG);
  double L2_z = (N2_z*SIG - N2*SIG_z)/(SIG*SIG);

  double L3_x = -z*rx/r2;
  double L3_y = -z*ry/r2;
  double L3_z = (r - z*rz)/r2;

  // Pack compact variables and their Cartesian derivatives.
  double L[3] = { L1, L2, L3 };
  double dL[3][3] = {
    { L1_x, L2_x, L3_x },
    { L1_y, L2_y, L3_y },
    { L1_z, L2_z, L3_z }
  };
  double dHK[3] = { HK_x, HK_y, HK_z };

  // Compact inverse spatial metric h^{ij} = delta^{ij} - kappa L^i L^j.
  double s = 1.0 + 2.0*HK;
  double kappa = 2.0*HK/s;

  double hcon[3][3];
  double hcov[3][3];
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) {
      double delta = (i == j) ? 1.0 : 0.0;
      hcon[i][j] = delta - kappa*L[i]*L[j];
      hcov[i][j] = delta + 2.0*HK*L[i]*L[j];
    }
  }

  // Cholesky-like factors for the contravariant spatial vierbein.
  double h11 = hcon[0][0];
  double h12 = hcon[0][1];
  double h13 = hcon[0][2];
  double h22 = hcon[1][1];
  double h23 = hcon[1][2];
  double h33 = hcon[2][2];

  double B = sqrt(h11);
  double D = h12/B;
  double Q = h22 - h12*h12/h11;
  double E = sqrt(Q);
  double G = h13/B;
  double Nch = h23 - h13*h12/h11;
  double Ch = Nch/E;
  double Iarg = h33 - h13*h13/h11 - Nch*Nch/Q;
  double I = sqrt(Iarg);

  // Spatial block of e^i_a. This matches the MATLAB sign convention.
  double V[3][3] = {
    { -B,  0.0, 0.0 },
    { -D, -E,  0.0 },
    { -G, -Ch, -I }
  };

  // Forward-mode derivative propagation for q = x,y,z.
  for (int q=0; q<3; ++q) {
    double dkappa = 2.0*dHK[q]/(s*s);

    double dhcon[3][3];
    double dhcov[3][3];
    for (int i=0; i<3; ++i) {
      for (int j=0; j<3; ++j) {
        dhcon[i][j] = -dkappa*L[i]*L[j]
                     - kappa*(dL[q][i]*L[j] + L[i]*dL[q][j]);
        dhcov[i][j] = 2.0*dHK[q]*L[i]*L[j]
                    + 2.0*HK*(dL[q][i]*L[j] + L[i]*dL[q][j]);
      }
    }

    double dh11 = dhcon[0][0];
    double dh12 = dhcon[0][1];
    double dh13 = dhcon[0][2];
    double dh22 = dhcon[1][1];
    double dh23 = dhcon[1][2];
    double dh33 = dhcon[2][2];

    double dB = dh11/(2.0*B);
    double dD = dh12/B - h12*dB/(B*B);
    double dQ = dh22 - (2.0*h12*dh12*h11 - h12*h12*dh11)/(h11*h11);
    double dE = dQ/(2.0*E);
    double dG = dh13/B - h13*dB/(B*B);
    double dNch = dh23 - (dh13*h12 + h13*dh12)/h11 + h13*h12*dh11/(h11*h11);
    double dCh = dNch/E - Nch*dE/(E*E);
    double dIarg = dh33
      - (2.0*h13*dh13*h11 - h13*h13*dh11)/(h11*h11)
      - (2.0*Nch*dNch*Q - Nch*Nch*dQ)/(Q*Q);
    double dI = dIarg/(2.0*I);

    double dV[3][3] = {
      { -dB,  0.0, 0.0 },
      { -dD, -dE,  0.0 },
      { -dG, -dCh, -dI }
    };

    // e_cov_spatial = -hcov * V, so d e = -(dhcov * V + hcov * dV).
    for (int i=0; i<3; ++i) {
      for (int alpha=0; alpha<3; ++alpha) {
        double de = 0.0;
        for (int j=0; j<3; ++j) {
          de += dhcov[i][j]*V[j][alpha] + hcov[i][j]*dV[j][alpha];
        }
        fout[9*q + 3*i + alpha] = -de;
      }
    }
  }
}


// for use in kernel tables
typedef struct { evalf_t kernels[3]; } hamil_kern_list;
typedef struct { evalf_t kernels[3]; } vierbein_kern_list;
typedef struct { evalf_t kernels[3]; } vierbein_gradient_kern_list;

static const hamil_kern_list ks_rphi_hamil_list[] = {
  { NULL, NULL, NULL },
  { eval_ks_rphi_hamil_1x2v, eval_ks_rphi_hamil_2x2v, NULL },
  { NULL, NULL, NULL }
};

static const hamil_kern_list ks_r_hamil_list[] = {
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { eval_ks_r_hamil_1x3v, NULL, NULL }
};

static const hamil_kern_list ks_rtheta_hamil_list[] = {
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, eval_ks_rtheta_hamil_2x3v, eval_ks_hamil_3x3v }
};

static const hamil_kern_list ks_hamil_3v_list[] = {
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, eval_ks_hamil_3x3v }
};

static const hamil_kern_list ks_hamil_cart_3v_list[] = {
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { eval_ks_hamil_cart_1x3v, eval_ks_hamil_cart_2x3v, eval_ks_hamil_cart_3x3v }
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

static const vierbein_kern_list ks_r_vierbein_list[] = {
  { NULL },
  { NULL },
  { eval_ks_r_vierbein_3v }
};

static const vierbein_gradient_kern_list ks_r_vierbein_gradient_list[] = {
  { NULL },
  { NULL },
  { eval_ks_r_vierbein_gradient_3v }
};

static const vierbein_kern_list ks_rtheta_vierbein_list[] = {
  { NULL },
  { NULL },
  { eval_ks_rtheta_vierbein_3v }
};

static const vierbein_gradient_kern_list ks_rtheta_vierbein_gradient_list[] = {
  { NULL },
  { NULL },
  { eval_ks_rtheta_vierbein_gradient_3v }
};

static const vierbein_kern_list ks_vierbein_list[] = {
  { NULL },
  { NULL},
  { eval_ks_vierbein_3v  }
};

static const vierbein_kern_list ks_vierbein_cart_list[] = {
  { NULL , NULL , NULL },
  { NULL , NULL , NULL},
  { eval_ks_vierbein_cart_1x3v , eval_ks_vierbein_cart_2x3v , eval_ks_vierbein_cart_3x3v  }
};

static const vierbein_gradient_kern_list ks_vierbein_gradient_list[] = {
  { NULL },
  { NULL },
  { eval_ks_vierbein_gradient_3v }
};

static const vierbein_gradient_kern_list ks_vierbein_gradient_cart_list[] = {
  { NULL , NULL , NULL },
  { NULL , NULL , NULL},
  { eval_ks_vierbein_gradient_cart_1x3v , eval_ks_vierbein_gradient_cart_2x3v , eval_ks_vierbein_gradient_cart_3x3v  }
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
choose_vierbein_kern(enum gkyl_triad_preset_geom_type type, int cdim, int vdim)
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
    case GKYL_TRIAD_GR_KERR_SCHILD_R:
      return ks_r_vierbein_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_GR_KERR_SCHILD_RTHETA:
      return ks_rtheta_vierbein_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_GR_KERR_SCHILD_3V:
      return ks_vierbein_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_CART_GR_KERR_SCHILD_3V:
      return ks_vierbein_cart_list[vdim-1].kernels[cdim - 1];
      break;
    default:
      assert(false);
  }
}

static evalf_t
choose_vierbein_gradient_kern(enum gkyl_triad_preset_geom_type type, int cdim, int vdim)
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
    case GKYL_TRIAD_GR_KERR_SCHILD_R:
      return ks_r_vierbein_gradient_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_GR_KERR_SCHILD_RTHETA:
      return ks_rtheta_vierbein_gradient_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_GR_KERR_SCHILD_3V:
      return ks_vierbein_gradient_list[vdim-1].kernels[0];
      break;
    case GKYL_TRIAD_CART_GR_KERR_SCHILD_3V:
      return ks_vierbein_gradient_cart_list[vdim-1].kernels[cdim - 1];
      break;
    default:
      assert(false);
  }
}

static evalf_t
choose_hamil_kern(enum gkyl_triad_preset_geom_type type, int cdim, int vdim)
{
  switch(type) {
    case GKYL_TRIAD_GR_KERR_SCHILD_RPHI:
      return ks_rphi_hamil_list[vdim-1].kernels[cdim-1];
      break;
    case GKYL_TRIAD_GR_KERR_SCHILD_R:
      return ks_r_hamil_list[vdim-1].kernels[cdim-1];
      break;
    case GKYL_TRIAD_GR_KERR_SCHILD_RTHETA:
      return ks_rtheta_hamil_list[vdim-1].kernels[cdim-1];
      break;
    case GKYL_TRIAD_GR_KERR_SCHILD_3V:
      return ks_hamil_3v_list[vdim-1].kernels[cdim-1];
      break;
    case GKYL_TRIAD_CART_GR_KERR_SCHILD_3V:
      return ks_hamil_cart_3v_list[vdim-1].kernels[cdim-1];
      break;
    default:
      assert(false);
  }
}
