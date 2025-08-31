// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

static void 
kernel_metric_1x(const double *cov_tangent_basis, double* GKYL_RESTRICT h_ij) 
{ 
  // Finds the metric 
  // cov_tangent_basis: components of the covariant tangent basis 
  // h_ij: metric 

  h_ij[0] = cov_tangent_basis[0] * cov_tangent_basis[0];
} 

static void 
kernel_metric_2x(const double *cov_tangent_basis, double* GKYL_RESTRICT h_ij) 
{ 
  // Finds the metric 
  // cov_tangent_basis: components of the covariant tangent basis
  // h_ij: metric 

  int vdim = 2;
  for (int i=0; i<vdim; ++i) {
    for (int j=i; j<vdim; ++j) {
      int sym_tensor_index = (i*(2*vdim - i + 1))/2 + (j-i);
       h_ij[sym_tensor_index] = 0.0;
      for (int k=0; k<vdim; ++k) {
        h_ij[sym_tensor_index] += cov_tangent_basis[i*vdim+k] * cov_tangent_basis[j*vdim+k];
      }
    }
  }
} 

static void 
kernel_metric_3x(const double *cov_tangent_basis, double* GKYL_RESTRICT h_ij) 
{ 
  // Finds the metric 
  // cov_tangent_basis: components of the covariant tangent basis
  // h_ij: metric

  int vdim = 3;
  for (int i=0; i<vdim; ++i) {
    for (int j=i; j<vdim; ++j) {
      int sym_tensor_index = (i*(2*vdim - i + 1))/2 + (j-i);
       h_ij[sym_tensor_index] = 0.0;
      for (int k=0; k<vdim; ++k) {
        h_ij[sym_tensor_index] += cov_tangent_basis[i*vdim+k] * cov_tangent_basis[j*vdim+k];
      }
    }
  }
} 



static void 
kernel_metric_inv_1x(const double *h_ij, double* GKYL_RESTRICT h_ij_inv) 
{ 
  // Finds the inverse of the metric 
  // h_ij: metric 
  // h_ij_inv: metric inverse

  h_ij_inv[0] = 1.0/h_ij[0];
} 

static void 
kernel_metric_inv_2x(const double *h_ij, double* GKYL_RESTRICT h_ij_inv) 
{ 
  // Finds the inverse of the metric 
  // h_ij: metric 
  // h_ij_inv: metric inverse

  double det_h = h_ij[0]*h_ij[2] - h_ij[1]*h_ij[1];
  h_ij_inv[0] = h_ij[2]/det_h;
  h_ij_inv[1] = -h_ij[1]/det_h;
  h_ij_inv[2] = h_ij[0]/det_h;

} 

static void 
kernel_metric_inv_3x(const double *h_ij, double* GKYL_RESTRICT h_ij_inv) 
{ 
  // Finds the inverse of the metric 
  // h_ij: metric 
  // h_ij_inv: metric inverse

  double A =   h_ij[3] * h_ij[5] - h_ij[4] * h_ij[4];
  double D = -(h_ij[1] * h_ij[5] - h_ij[2] * h_ij[4]);
  double E =   h_ij[0] * h_ij[5] - h_ij[2] * h_ij[2];
  double G =   h_ij[1] * h_ij[4] - h_ij[2] * h_ij[3];
  double H = -(h_ij[0] * h_ij[4] - h_ij[2] * h_ij[1]);
  double I =   h_ij[0] * h_ij[3] - h_ij[1] * h_ij[1];

  double det_h = (h_ij[0] * A) + (h_ij[1] * D) + (h_ij[2] * G);

  h_ij_inv[0] = A/det_h;
  h_ij_inv[1] = D/det_h;
  h_ij_inv[2] = G/det_h;
  h_ij_inv[3] = E/det_h;
  h_ij_inv[4] = H/det_h;
  h_ij_inv[5] = I/det_h;

} 


static void 
kernel_metric_det_1x(const double *h_ij, double* GKYL_RESTRICT det_h) 
{ 
  // Finds the determinant of the metric 
  // h_ij: metric 
  // det_h: (squareroot) metric determinant

  det_h[0] = sqrt( h_ij[0] );

} 

static void 
kernel_metric_det_2x(const double *h_ij, double* GKYL_RESTRICT det_h) 
{ 
  // Finds the determinant of the metric 
  // h_ij: metric 
  // det_h: (squareroot) metric determinant

  det_h[0] = sqrt( h_ij[0]*h_ij[2] - h_ij[1]*h_ij[1] );

} 

static void 
kernel_metric_det_3x(const double *h_ij, double* GKYL_RESTRICT det_h) 
{ 
  // Finds the determinant of the metric 
  // h_ij: metric 
  // det_h: (squareroot) metric determinant

  double A =   h_ij[3] * h_ij[5] - h_ij[4] * h_ij[4];
  double D = -(h_ij[1] * h_ij[5] - h_ij[2] * h_ij[4]);
  double G =   h_ij[1] * h_ij[4] - h_ij[2] * h_ij[3];
  det_h[0] = sqrt( (h_ij[0] * A) + (h_ij[1] * D) + (h_ij[2] * G) );

} 


static void 
kernel_conf_poisson_tensor_1x1v(const double *h_ij_inv, const double *triad_basis,
  const double *cov_tangent_basis, const double *triad_basis_gradient, double* GKYL_RESTRICT conf_poisson_tensor) 
{ 
  // Finds the configuration space components of the Poisson tensor
  // h_ij_inv: inverse metric 
  // triad_basis: triad basis coefficients
  // cov_tangent_basis: covariant tangent basis coefficients
  // triad_basis_gradient: triad basis gradient coefficients
  // (out) conf_poisson_tensor:  configuration space components of the Poisson tensor

} 

static void 
kernel_conf_poisson_tensor_1x2v(const double *h_ij_inv, const double *triad_basis,
  const double *cov_tangent_basis, const double *triad_basis_gradient, double* GKYL_RESTRICT conf_poisson_tensor) 
{ 
  // Finds the configuration space components of the Poisson tensor
  // h_ij_inv: inverse metric 
  // triad_basis: triad basis coefficients
  // cov_tangent_basis: covariant tangent basis coefficients
  // triad_basis_gradient: triad basis gradient coefficients
  // (out) conf_poisson_tensor:  configuration space components of the Poisson tensor

} 

static void 
kernel_conf_poisson_tensor_1x3v(const double *h_ij_inv, const double *triad_basis,
  const double *cov_tangent_basis, const double *triad_basis_gradient, double* GKYL_RESTRICT conf_poisson_tensor) 
{ 
  // Finds the configuration space components of the Poisson tensor
  // h_ij_inv: inverse metric 
  // triad_basis: triad basis coefficients
  // cov_tangent_basis: covariant tangent basis coefficients
  // triad_basis_gradient: triad basis gradient coefficients
  // (out) conf_poisson_tensor:  configuration space components of the Poisson tensor

} 

static void 
kernel_conf_poisson_tensor_2x2v(const double *h_ij_inv, const double *triad_basis,
  const double *cov_tangent_basis, const double *triad_basis_gradient, double* GKYL_RESTRICT conf_poisson_tensor) 
{ 
  // Finds the configuration space components of the Poisson tensor
  // h_ij_inv: inverse metric 
  // triad_basis: triad basis coefficients
  // cov_tangent_basis: covariant tangent basis coefficients
  // triad_basis_gradient: triad basis gradient coefficients
  // (out) conf_poisson_tensor:  configuration space components of the Poisson tensor

} 

static void 
kernel_conf_poisson_tensor_2x3v(const double *h_ij_inv, const double *triad_basis,
  const double *cov_tangent_basis, const double *triad_basis_gradient, double* GKYL_RESTRICT conf_poisson_tensor) 
{ 
  // Finds the configuration space components of the Poisson tensor
  // h_ij_inv: inverse metric 
  // triad_basis: triad basis coefficients
  // cov_tangent_basis: covariant tangent basis coefficients
  // triad_basis_gradient: triad basis gradient coefficients
  // (out) conf_poisson_tensor:  configuration space components of the Poisson tensor

} 

static void 
compute_nu_inv_3x(const double *h_ij_inv, const double *triad_basis, const double *cov_tangent_basis,
  double* GKYL_RESTRICT nu_inv)
{
  // Compute nu_inv by raising the indicies:
  // nu^{ij} = \delta^{im} h^{jn} (\sigma_m \cdot g_n)
  //         = h^{jn} (\sigma^i \cdot g_n)
  //        = h^{jn} (\sigma^{i,m} g_{n,m})

  /* Shape of nu_inv.
  +----+----+----+
  | 00 | 01 | 02 |
  +----+----+----+
  | 03 | 04 | 05 |
  +----+----+----+
  | 06 | 07 | 08 |
  +----+----+----+
  */

  /* Shape of h_ij_inv
  +----+----+----+
  | 00 | 01 | 02 |
  +----+----+----+
  | -- | 03 | 04 |
  +----+----+----+
  | -- | -- | 05 |
  +----+----+----+
  */

  int vdim = 3;

  double h_ij_inv_local[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  h_ij_inv_local[0] = h_ij_inv[0];
  h_ij_inv_local[1] = h_ij_inv[1];
  h_ij_inv_local[2] = h_ij_inv[2];
  h_ij_inv_local[3] = h_ij_inv[1];
  h_ij_inv_local[4] = h_ij_inv[3];
  h_ij_inv_local[5] = h_ij_inv[4];
  h_ij_inv_local[6] = h_ij_inv[2];
  h_ij_inv_local[7] = h_ij_inv[4];
  h_ij_inv_local[8] = h_ij_inv[5];


  for (int i=0; i<vdim; ++i) {
    for (int j=0; j<vdim; ++j) {
      nu_inv[i*vdim + j] = 0;
      for (int n=0; n<vdim; ++n) {
        for (int m=0; m<vdim; ++m) {
          nu_inv[i*vdim + j] += h_ij_inv_local[j*vdim + n] * ( cov_tangent_basis[i*vdim+m] * cov_tangent_basis[n*vdim+m] );
        }
      }
    }
  }
}


static void 
kernel_conf_poisson_tensor_3x3v(const double *h_ij_inv, const double *triad_basis,
  const double *cov_tangent_basis, const double *triad_basis_gradient, double* GKYL_RESTRICT conf_poisson_tensor) 
{ 
  // Finds the configuration space components of the Poisson tensor
  // h_ij_inv: inverse metric 
  // triad_basis: triad basis coefficients
  // cov_tangent_basis: covariant tangent basis coefficients
  // triad_basis_gradient: triad basis gradient coefficients
  // (out) conf_poisson_tensor:  configuration space components of the Poisson tensor

  /* Shape of conf_poisson_tensor non-zero elements.
  +----+----+----+----+----+----+
  | -- | -- | -- | 00 | 01 | 02 |
  +----+----+----+----+----+----+
  | -- | -- | -- | 03 | 04 | 05 |
  +----+----+----+----+----+----+
  | -- | -- | -- | 06 | 07 | 08 |
  +----+----+----+----+----+----+
  | -- | -- | -- | -- |09+j|10+j|
  +----+----+----+----+----+----+
  | -- | -- | -- | -- | -- |11+j|
  +----+----+----+----+----+----+
  | -- | -- | -- | -- | -- | -- |
  +----+----+----+----+----+----+
  */

  /* Shape of nu_inv.
  +----+----+----+
  | 00 | 01 | 02 |
  +----+----+----+
  | 03 | 04 | 05 |
  +----+----+----+
  | 06 | 07 | 08 |
  +----+----+----+
  */


  int cdim = 3;
  int vdim = 3;
  int pdim = cdim + vdim;
  int n_sym_tensor_indices = 3;
  int n_nonzero_indices = vdim*vdim + n_sym_tensor_indices;

  // Zero out the Poisson Tensor, top left block is zero
  for (int i=0; i<n_nonzero_indices; ++i) {
    conf_poisson_tensor[i] = 0.0;
  }

  double nu_inv[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  // Compute nu_inv by raising the indicies:
  // nu^{ij} = \delta^{im} h^{jn} (\sigma_m \cdot g_n)
  //         = h^{jn} (\sigma^i \cdot g_n)
  compute_nu_inv_3x(h_ij_inv, triad_basis, cov_tangent_basis, nu_inv);

  // Top right block:
  conf_poisson_tensor[0] = nu_inv[0];
  conf_poisson_tensor[1] = nu_inv[1];
  conf_poisson_tensor[2] = nu_inv[2];
  conf_poisson_tensor[3] = nu_inv[3];
  conf_poisson_tensor[4] = nu_inv[4];
  conf_poisson_tensor[5] = nu_inv[5];  
  conf_poisson_tensor[6] = nu_inv[6];
  conf_poisson_tensor[7] = nu_inv[7];
  conf_poisson_tensor[8] = nu_inv[8];


  // Bottom right block 
  for (int i=0; i<vdim; ++i) {
    for (int j=0; j<vdim; ++j) {
      for (int k=0; k<vdim; ++k) {
        for (int p=0; p<3; ++p ){

         // index the free component of the antisymmetric tensor
         int asym_ten_indx = p;
         int m = 0;
         int n = 0;
         if (p == 0) {
          m = 0;
          n = 1;
         }
         else if (p == 1) {
          m = 0;
          n = 2;
         }
         else if (p == 2) {
          m = 1;
          n = 2;
         }

          // Diagonal elements are zero
          double gj_dot_d_sigma_i_dx_k = 0.0;
          double gi_dot_d_sigma_j_dx_k = 0.0;
          for (int p=0; p<vdim; ++p) {
            gj_dot_d_sigma_i_dx_k += cov_tangent_basis[vdim*j + p] * triad_basis_gradient[vdim*vdim*k + i*vdim + p];
            gi_dot_d_sigma_j_dx_k += cov_tangent_basis[vdim*i + p] * triad_basis_gradient[vdim*vdim*k + j*vdim + p];
          }
          // Offset by vdim^2
          conf_poisson_tensor[vdim*vdim + n_sym_tensor_indices*j + p] +=  nu_inv[m*vdim + k] * nu_inv[n*vdim + i] * 
            ( gj_dot_d_sigma_i_dx_k -  gi_dot_d_sigma_j_dx_k ); 
        }
      }
    }
  }
} 



typedef void (*metric_t)(const double *pn_cov_tangent_basis, double* GKYL_RESTRICT pn_h_ij);
typedef void (*metric_inv_t)(const double *pn_h_ij, double* GKYL_RESTRICT pn_h_ij_inv);
typedef void (*metric_det_t)(const double *pn_h_ij, double* GKYL_RESTRICT pn_det_h);
typedef void (*conf_poisson_tensor_t)(const double *pn_h_ij_inv, const double *pn_triad_basis,
  const double *pn_cov_tangent_basis, const double *pn_triad_basis_gradient,
  double* GKYL_RESTRICT pn_conf_poisson_tensor);


  
// for use in kernel tables
typedef struct { metric_t kernels[3]; } metric_kern_list;
typedef struct { metric_inv_t kernels[3]; } metric_inv_kern_list;
typedef struct { metric_det_t kernels[3]; } metric_det_kern_list;
typedef struct { conf_poisson_tensor_t kernels[6]; } conf_poisson_tensor_kern_list;

static const metric_kern_list metric_list[] = {
  { kernel_metric_1x },
  { kernel_metric_2x },
  { kernel_metric_3x }
};

static const metric_inv_kern_list metric_inv_list[] = {
  { kernel_metric_inv_1x },
  { kernel_metric_inv_2x },
  { kernel_metric_inv_3x }
};

static const metric_det_kern_list metric_det_list[] = {
  { kernel_metric_det_1x },
  { kernel_metric_det_2x },
  { kernel_metric_det_3x }
};

static const conf_poisson_tensor_kern_list conf_poisson_tensor_list[] = {
  { kernel_conf_poisson_tensor_1x1v, kernel_conf_poisson_tensor_1x2v, kernel_conf_poisson_tensor_1x3v },
  { NULL, kernel_conf_poisson_tensor_2x2v, kernel_conf_poisson_tensor_2x3v },
  { NULL, NULL, kernel_conf_poisson_tensor_3x3v }
};


static metric_t
choose_metric_kern(int vdim)
{
  return metric_list[vdim-1].kernels[0];
}

static metric_inv_t
choose_metric_inv_kern(int vdim)
{
  return metric_inv_list[vdim-1].kernels[0];
}

static metric_det_t
choose_metric_det_kern(int vdim)
{
  return metric_det_list[vdim-1].kernels[0];
}

static conf_poisson_tensor_t
choose_conf_poisson_tensor_kern(int cdim, int vdim)
{
  return conf_poisson_tensor_list[cdim-1].kernels[vdim-1];
}