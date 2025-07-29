#pragma once

// Private header, not for direct use in user code

#include <gkyl_array.h>
#include <gkyl_mom_type.h>
#include <gkyl_ref_count.h>
#include <gkyl_mom_vlasov_kernels.h>

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
GKYL_CU_D
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1,  3,  4,  5}, // 2x kernel indices
  {-1, -1, -1,  6}, // 3x kernel indices  
};

typedef void (*vlasov_momf_t)(const double *xc, const double *dx,
  const int *idx, const double* hamil, const double *fIn, double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { vlasov_momf_t kernels[4]; } gkyl_mom_kern_list;

// M0 kernel list (Serendipity basis)
GKYL_CU_D
static const gkyl_mom_kern_list ser_m0_kernels[] = {
  // 1x kernels
  { NULL, mom_vlasov_M0_1x1v_ser_p1, mom_vlasov_M0_1x1v_ser_p2, NULL }, // 0
  { NULL, mom_vlasov_M0_1x2v_ser_p1, mom_vlasov_M0_1x2v_ser_p2, NULL }, // 1
  { NULL, mom_vlasov_M0_1x3v_ser_p1, mom_vlasov_M0_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, mom_vlasov_M0_2x1v_ser_p1, mom_vlasov_M0_2x1v_ser_p2, NULL }, // 3
  { NULL, mom_vlasov_M0_2x2v_ser_p1, mom_vlasov_M0_2x2v_ser_p2, NULL }, // 4
  { NULL, mom_vlasov_M0_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, mom_vlasov_M0_3x3v_ser_p1, NULL, NULL }, // 6
};

// M0 kernel list (Tensor basis)
GKYL_CU_D
static const gkyl_mom_kern_list tensor_m0_kernels[] = {
  // 1x kernels
  { NULL, mom_vlasov_M0_1x1v_tensor_p1, mom_vlasov_M0_1x1v_tensor_p2, NULL }, // 0
  { NULL, mom_vlasov_M0_1x2v_tensor_p1, mom_vlasov_M0_1x2v_tensor_p2, NULL }, // 1
  { NULL, mom_vlasov_M0_1x3v_tensor_p1, mom_vlasov_M0_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, mom_vlasov_M0_2x1v_tensor_p1, mom_vlasov_M0_2x1v_tensor_p2, NULL }, // 3
  { NULL, mom_vlasov_M0_2x2v_tensor_p1, mom_vlasov_M0_2x2v_tensor_p2, NULL }, // 4
  { NULL, mom_vlasov_M0_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, mom_vlasov_M0_3x3v_tensor_p1, NULL, NULL }, // 6
};

// M1i (dH/dv moment, velocity-space Hamiltonian) kernel list (Serendipity basis)
GKYL_CU_D
static const gkyl_mom_kern_list ser_m1i_hamil_vel_kernels[] = {
  // 1x kernels
  { NULL, mom_vlasov_hamil_vel_M1i_1x1v_ser_p1, mom_vlasov_hamil_vel_M1i_1x1v_ser_p2, NULL }, // 0
  { NULL, mom_vlasov_hamil_vel_M1i_1x2v_ser_p1, mom_vlasov_hamil_vel_M1i_1x2v_ser_p2, NULL }, // 1
  { NULL, mom_vlasov_hamil_vel_M1i_1x3v_ser_p1, mom_vlasov_hamil_vel_M1i_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, mom_vlasov_hamil_vel_M1i_2x1v_ser_p1, mom_vlasov_hamil_vel_M1i_2x1v_ser_p2, NULL }, // 3
  { NULL, mom_vlasov_hamil_vel_M1i_2x2v_ser_p1, mom_vlasov_hamil_vel_M1i_2x2v_ser_p2, NULL }, // 4
  { NULL, mom_vlasov_hamil_vel_M1i_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, mom_vlasov_hamil_vel_M1i_3x3v_ser_p1, NULL, NULL }, // 6
};

// M1i (dH/dv moment, velocity-space Hamiltonian) kernel list (Tensor basis)
GKYL_CU_D
static const gkyl_mom_kern_list tensor_m1i_hamil_vel_kernels[] = {
  // 1x kernels
  { NULL, mom_vlasov_hamil_vel_M1i_1x1v_tensor_p1, mom_vlasov_hamil_vel_M1i_1x1v_tensor_p2, NULL }, // 0
  { NULL, mom_vlasov_hamil_vel_M1i_1x2v_tensor_p1, mom_vlasov_hamil_vel_M1i_1x2v_tensor_p2, NULL }, // 1
  { NULL, mom_vlasov_hamil_vel_M1i_1x3v_tensor_p1, mom_vlasov_hamil_vel_M1i_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, mom_vlasov_hamil_vel_M1i_2x1v_tensor_p1, mom_vlasov_hamil_vel_M1i_2x1v_tensor_p2, NULL }, // 3
  { NULL, mom_vlasov_hamil_vel_M1i_2x2v_tensor_p1, mom_vlasov_hamil_vel_M1i_2x2v_tensor_p2, NULL }, // 4
  { NULL, mom_vlasov_hamil_vel_M1i_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, mom_vlasov_hamil_vel_M1i_3x3v_tensor_p1, NULL, NULL }, // 6
};

// M1i (dH/dv moment, phase-space Hamiltonian) kernel list (Serendipity basis)
GKYL_CU_D
static const gkyl_mom_kern_list ser_m1i_hamil_gen_kernels[] = {
  // 1x kernels
  { NULL, mom_vlasov_hamil_gen_M1i_1x1v_ser_p1, mom_vlasov_hamil_gen_M1i_1x1v_ser_p2, NULL }, // 0
  { NULL, mom_vlasov_hamil_gen_M1i_1x2v_ser_p1, mom_vlasov_hamil_gen_M1i_1x2v_ser_p2, NULL }, // 1
  { NULL, mom_vlasov_hamil_gen_M1i_1x3v_ser_p1, mom_vlasov_hamil_gen_M1i_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, mom_vlasov_hamil_gen_M1i_2x1v_ser_p1, mom_vlasov_hamil_gen_M1i_2x1v_ser_p2, NULL }, // 3
  { NULL, mom_vlasov_hamil_gen_M1i_2x2v_ser_p1, mom_vlasov_hamil_gen_M1i_2x2v_ser_p2, NULL }, // 4
  { NULL, mom_vlasov_hamil_gen_M1i_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, mom_vlasov_hamil_gen_M1i_3x3v_ser_p1, NULL, NULL }, // 6
};

// M1i (dH/dv moment, phase-space Hamiltonian) kernel list (Tensor basis)
GKYL_CU_D
static const gkyl_mom_kern_list tensor_m1i_hamil_gen_kernels[] = {
  // 1x kernels
  { NULL, mom_vlasov_hamil_gen_M1i_1x1v_tensor_p1, mom_vlasov_hamil_gen_M1i_1x1v_tensor_p2, NULL }, // 0
  { NULL, mom_vlasov_hamil_gen_M1i_1x2v_tensor_p1, mom_vlasov_hamil_gen_M1i_1x2v_tensor_p2, NULL }, // 1
  { NULL, mom_vlasov_hamil_gen_M1i_1x3v_tensor_p1, mom_vlasov_hamil_gen_M1i_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, mom_vlasov_hamil_gen_M1i_2x1v_tensor_p1, mom_vlasov_hamil_gen_M1i_2x1v_tensor_p2, NULL }, // 3
  { NULL, mom_vlasov_hamil_gen_M1i_2x2v_tensor_p1, mom_vlasov_hamil_gen_M1i_2x2v_tensor_p2, NULL }, // 4
  { NULL, mom_vlasov_hamil_gen_M1i_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, mom_vlasov_hamil_gen_M1i_3x3v_tensor_p1, NULL, NULL }, // 6
};

// M2 (H moment, velocity-space Hamiltonian) kernel list (Serendipity basis)
GKYL_CU_D
static const gkyl_mom_kern_list ser_m2_hamil_vel_kernels[] = {
  // 1x kernels
  { NULL, mom_vlasov_hamil_vel_M2_1x1v_ser_p1, mom_vlasov_hamil_vel_M2_1x1v_ser_p2, NULL }, // 0
  { NULL, mom_vlasov_hamil_vel_M2_1x2v_ser_p1, mom_vlasov_hamil_vel_M2_1x2v_ser_p2, NULL }, // 1
  { NULL, mom_vlasov_hamil_vel_M2_1x3v_ser_p1, mom_vlasov_hamil_vel_M2_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, mom_vlasov_hamil_vel_M2_2x1v_ser_p1, mom_vlasov_hamil_vel_M2_2x1v_ser_p2, NULL }, // 3
  { NULL, mom_vlasov_hamil_vel_M2_2x2v_ser_p1, mom_vlasov_hamil_vel_M2_2x2v_ser_p2, NULL }, // 4
  { NULL, mom_vlasov_hamil_vel_M2_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, mom_vlasov_hamil_vel_M2_3x3v_ser_p1, NULL, NULL }, // 6
};

// M2 (H moment, velocity-space Hamiltonian) kernel list (Tensor basis)
GKYL_CU_D
static const gkyl_mom_kern_list tensor_m2_hamil_vel_kernels[] = {
  // 1x kernels
  { NULL, mom_vlasov_hamil_vel_M2_1x1v_tensor_p1, mom_vlasov_hamil_vel_M2_1x1v_tensor_p2, NULL }, // 0
  { NULL, mom_vlasov_hamil_vel_M2_1x2v_tensor_p1, mom_vlasov_hamil_vel_M2_1x2v_tensor_p2, NULL }, // 1
  { NULL, mom_vlasov_hamil_vel_M2_1x3v_tensor_p1, mom_vlasov_hamil_vel_M2_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, mom_vlasov_hamil_vel_M2_2x1v_tensor_p1, mom_vlasov_hamil_vel_M2_2x1v_tensor_p2, NULL }, // 3
  { NULL, mom_vlasov_hamil_vel_M2_2x2v_tensor_p1, mom_vlasov_hamil_vel_M2_2x2v_tensor_p2, NULL }, // 4
  { NULL, mom_vlasov_hamil_vel_M2_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, mom_vlasov_hamil_vel_M2_3x3v_tensor_p1, NULL, NULL }, // 6
};

// M2 (H moment, phase-space Hamiltonian) kernel list (Serendipity basis)
GKYL_CU_D
static const gkyl_mom_kern_list ser_m2_hamil_gen_kernels[] = {
  // 1x kernels
  { NULL, mom_vlasov_hamil_gen_M2_1x1v_ser_p1, mom_vlasov_hamil_gen_M2_1x1v_ser_p2, NULL }, // 0
  { NULL, mom_vlasov_hamil_gen_M2_1x2v_ser_p1, mom_vlasov_hamil_gen_M2_1x2v_ser_p2, NULL }, // 1
  { NULL, mom_vlasov_hamil_gen_M2_1x3v_ser_p1, mom_vlasov_hamil_gen_M2_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, mom_vlasov_hamil_gen_M2_2x1v_ser_p1, mom_vlasov_hamil_gen_M2_2x1v_ser_p2, NULL }, // 3
  { NULL, mom_vlasov_hamil_gen_M2_2x2v_ser_p1, mom_vlasov_hamil_gen_M2_2x2v_ser_p2, NULL }, // 4
  { NULL, mom_vlasov_hamil_gen_M2_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, mom_vlasov_hamil_gen_M2_3x3v_ser_p1, NULL, NULL }, // 6
};

// M2 (H moment, phase-space Hamiltonian) kernel list (Tensor basis)
GKYL_CU_D
static const gkyl_mom_kern_list tensor_m2_hamil_gen_kernels[] = {
  // 1x kernels
  { NULL, mom_vlasov_hamil_gen_M2_1x1v_tensor_p1, mom_vlasov_hamil_gen_M2_1x1v_tensor_p2, NULL }, // 0
  { NULL, mom_vlasov_hamil_gen_M2_1x2v_tensor_p1, mom_vlasov_hamil_gen_M2_1x2v_tensor_p2, NULL }, // 1
  { NULL, mom_vlasov_hamil_gen_M2_1x3v_tensor_p1, mom_vlasov_hamil_gen_M2_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, mom_vlasov_hamil_gen_M2_2x1v_tensor_p1, mom_vlasov_hamil_gen_M2_2x1v_tensor_p2, NULL }, // 3
  { NULL, mom_vlasov_hamil_gen_M2_2x2v_tensor_p1, mom_vlasov_hamil_gen_M2_2x2v_tensor_p2, NULL }, // 4
  { NULL, mom_vlasov_hamil_gen_M2_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, mom_vlasov_hamil_gen_M2_3x3v_tensor_p1, NULL, NULL }, // 6
};

// New energy flux kernels which are formally H*dH/dv are very dense, so not generating
// them yet. JJ 07/26/25
// // M3i (H*dH/dv moment, velocity-space Hamiltonian) kernel list (Serendipity basis)
// GKYL_CU_D
// static const gkyl_mom_kern_list ser_m3i_hamil_vel_kernels[] = {
//   // 1x kernels
//   { NULL, mom_vlasov_hamil_vel_M3i_1x1v_ser_p1, mom_vlasov_hamil_vel_M3i_1x1v_ser_p2, NULL }, // 0
//   { NULL, mom_vlasov_hamil_vel_M3i_1x2v_ser_p1, mom_vlasov_hamil_vel_M3i_1x2v_ser_p2, NULL }, // 1
//   { NULL, mom_vlasov_hamil_vel_M3i_1x3v_ser_p1, mom_vlasov_hamil_vel_M3i_1x3v_ser_p2, NULL }, // 2
//   // 2x kernels
//   { NULL, mom_vlasov_hamil_vel_M3i_2x1v_ser_p1, mom_vlasov_hamil_vel_M3i_2x1v_ser_p2, NULL }, // 3
//   { NULL, mom_vlasov_hamil_vel_M3i_2x2v_ser_p1, mom_vlasov_hamil_vel_M3i_2x2v_ser_p2, NULL }, // 4
//   { NULL, mom_vlasov_hamil_vel_M3i_2x3v_ser_p1, NULL, NULL }, // 5
//   // 3x kernels
//   { NULL, mom_vlasov_hamil_vel_M3i_3x3v_ser_p1, NULL, NULL }, // 6
// };

// // M3i (H*dH/dv moment, velocity-space Hamiltonian) kernel list (Tensor basis)
// GKYL_CU_D
// static const gkyl_mom_kern_list tensor_m3i_hamil_vel_kernels[] = {
//   // 1x kernels
//   { NULL, mom_vlasov_hamil_vel_M3i_1x1v_tensor_p1, mom_vlasov_hamil_vel_M3i_1x1v_tensor_p2, NULL }, // 0
//   { NULL, mom_vlasov_hamil_vel_M3i_1x2v_tensor_p1, mom_vlasov_hamil_vel_M3i_1x2v_tensor_p2, NULL }, // 1
//   { NULL, mom_vlasov_hamil_vel_M3i_1x3v_tensor_p1, mom_vlasov_hamil_vel_M3i_1x3v_tensor_p2, NULL }, // 2
//   // 2x kernels
//   { NULL, mom_vlasov_hamil_vel_M3i_2x1v_tensor_p1, mom_vlasov_hamil_vel_M3i_2x1v_tensor_p2, NULL }, // 3
//   { NULL, mom_vlasov_hamil_vel_M3i_2x2v_tensor_p1, mom_vlasov_hamil_vel_M3i_2x2v_tensor_p2, NULL }, // 4
//   { NULL, mom_vlasov_hamil_vel_M3i_2x3v_tensor_p1, NULL, NULL }, // 5
//   // 3x kernels
//   { NULL, mom_vlasov_hamil_vel_M3i_3x3v_tensor_p1, NULL, NULL }, // 6
// };

// // M3i (H*dH/dv moment, phase-space Hamiltonian) kernel list (Serendipity basis)
// GKYL_CU_D
// static const gkyl_mom_kern_list ser_m3i_hamil_gen_kernels[] = {
//   // 1x kernels
//   { NULL, mom_vlasov_hamil_gen_M3i_1x1v_ser_p1, mom_vlasov_hamil_gen_M3i_1x1v_ser_p2, NULL }, // 0
//   { NULL, mom_vlasov_hamil_gen_M3i_1x2v_ser_p1, mom_vlasov_hamil_gen_M3i_1x2v_ser_p2, NULL }, // 1
//   { NULL, mom_vlasov_hamil_gen_M3i_1x3v_ser_p1, mom_vlasov_hamil_gen_M3i_1x3v_ser_p2, NULL }, // 2
//   // 2x kernels
//   { NULL, mom_vlasov_hamil_gen_M3i_2x1v_ser_p1, mom_vlasov_hamil_gen_M3i_2x1v_ser_p2, NULL }, // 3
//   { NULL, mom_vlasov_hamil_gen_M3i_2x2v_ser_p1, mom_vlasov_hamil_gen_M3i_2x2v_ser_p2, NULL }, // 4
//   { NULL, mom_vlasov_hamil_gen_M3i_2x3v_ser_p1, NULL, NULL }, // 5
//   // 3x kernels
//   { NULL, mom_vlasov_hamil_gen_M3i_3x3v_ser_p1, NULL, NULL }, // 6
// };

// // M3i (H*dH/dv moment, phase-space Hamiltonian) kernel list (Tensor basis)
// GKYL_CU_D
// static const gkyl_mom_kern_list tensor_m3i_hamil_gen_kernels[] = {
//   // 1x kernels
//   { NULL, mom_vlasov_hamil_gen_M3i_1x1v_tensor_p1, mom_vlasov_hamil_gen_M3i_1x1v_tensor_p2, NULL }, // 0
//   { NULL, mom_vlasov_hamil_gen_M3i_1x2v_tensor_p1, mom_vlasov_hamil_gen_M3i_1x2v_tensor_p2, NULL }, // 1
//   { NULL, mom_vlasov_hamil_gen_M3i_1x3v_tensor_p1, mom_vlasov_hamil_gen_M3i_1x3v_tensor_p2, NULL }, // 2
//   // 2x kernels
//   { NULL, mom_vlasov_hamil_gen_M3i_2x1v_tensor_p1, mom_vlasov_hamil_gen_M3i_2x1v_tensor_p2, NULL }, // 3
//   { NULL, mom_vlasov_hamil_gen_M3i_2x2v_tensor_p1, mom_vlasov_hamil_gen_M3i_2x2v_tensor_p2, NULL }, // 4
//   { NULL, mom_vlasov_hamil_gen_M3i_2x3v_tensor_p1, NULL, NULL }, // 5
//   // 3x kernels
//   { NULL, mom_vlasov_hamil_gen_M3i_3x3v_tensor_p1, NULL, NULL }, // 6
// };

// M2ij kernel list (Serendipity basis)
GKYL_CU_D
static const gkyl_mom_kern_list ser_m2ij_kernels[] = {
  // 1x kernels
  { NULL, mom_vlasov_M2ij_1x1v_ser_p1, mom_vlasov_M2ij_1x1v_ser_p2, NULL }, // 0
  { NULL, mom_vlasov_M2ij_1x2v_ser_p1, mom_vlasov_M2ij_1x2v_ser_p2, NULL }, // 1
  { NULL, mom_vlasov_M2ij_1x3v_ser_p1, mom_vlasov_M2ij_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, mom_vlasov_M2ij_2x1v_ser_p1, mom_vlasov_M2ij_2x1v_ser_p2, NULL }, // 3
  { NULL, mom_vlasov_M2ij_2x2v_ser_p1, mom_vlasov_M2ij_2x2v_ser_p2, NULL }, // 4
  { NULL, mom_vlasov_M2ij_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, mom_vlasov_M2ij_3x3v_ser_p1, NULL, NULL }, // 6
};

// M2ij kernel list (Tensor basis)
GKYL_CU_D
static const gkyl_mom_kern_list tensor_m2ij_kernels[] = {
  // 1x kernels
  { NULL, mom_vlasov_M2ij_1x1v_tensor_p1, mom_vlasov_M2ij_1x1v_tensor_p2, NULL }, // 0
  { NULL, mom_vlasov_M2ij_1x2v_tensor_p1, mom_vlasov_M2ij_1x2v_tensor_p2, NULL }, // 1
  { NULL, mom_vlasov_M2ij_1x3v_tensor_p1, mom_vlasov_M2ij_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, mom_vlasov_M2ij_2x1v_tensor_p1, mom_vlasov_M2ij_2x1v_tensor_p2, NULL }, // 3
  { NULL, mom_vlasov_M2ij_2x2v_tensor_p1, mom_vlasov_M2ij_2x2v_tensor_p2, NULL }, // 4
  { NULL, mom_vlasov_M2ij_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, mom_vlasov_M2ij_3x3v_tensor_p1, NULL, NULL }, // 6
};

// M3ijk kernel list (Serendipity basis)
GKYL_CU_D
static const gkyl_mom_kern_list ser_m3ijk_kernels[] = {
  // 1x kernels
  { NULL, mom_vlasov_M3ijk_1x1v_ser_p1, mom_vlasov_M3ijk_1x1v_ser_p2, NULL }, // 0
  { NULL, mom_vlasov_M3ijk_1x2v_ser_p1, mom_vlasov_M3ijk_1x2v_ser_p2, NULL }, // 1
  { NULL, mom_vlasov_M3ijk_1x3v_ser_p1, mom_vlasov_M3ijk_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, mom_vlasov_M3ijk_2x1v_ser_p1, mom_vlasov_M3ijk_2x1v_ser_p2, NULL }, // 3
  { NULL, mom_vlasov_M3ijk_2x2v_ser_p1, mom_vlasov_M3ijk_2x2v_ser_p2, NULL }, // 4
  { NULL, mom_vlasov_M3ijk_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, mom_vlasov_M3ijk_3x3v_ser_p1, NULL, NULL }, // 6
};

// M3ijk kernel list (Tensor basis)
GKYL_CU_D
static const gkyl_mom_kern_list tensor_m3ijk_kernels[] = {
  // 1x kernels
  { NULL, mom_vlasov_M3ijk_1x1v_tensor_p1, mom_vlasov_M3ijk_1x1v_tensor_p2, NULL }, // 0
  { NULL, mom_vlasov_M3ijk_1x2v_tensor_p1, mom_vlasov_M3ijk_1x2v_tensor_p2, NULL }, // 1
  { NULL, mom_vlasov_M3ijk_1x3v_tensor_p1, mom_vlasov_M3ijk_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, mom_vlasov_M3ijk_2x1v_tensor_p1, mom_vlasov_M3ijk_2x1v_tensor_p2, NULL }, // 3
  { NULL, mom_vlasov_M3ijk_2x2v_tensor_p1, mom_vlasov_M3ijk_2x2v_tensor_p2, NULL }, // 4
  { NULL, mom_vlasov_M3ijk_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, mom_vlasov_M3ijk_3x3v_tensor_p1, NULL, NULL }, // 6
};

// Five moments (Zeroth, First, and Second moment together) kernel list (Serendipity basis)
// {1, dH/dv, H} moments for velocity-space Hamiltonian. 
GKYL_CU_D
static const gkyl_mom_kern_list ser_five_moments_hamil_vel_kernels[] = {
  // 1x kernels
  { NULL, mom_vlasov_hamil_vel_five_moments_1x1v_ser_p1, mom_vlasov_hamil_vel_five_moments_1x1v_ser_p2, NULL }, // 0
  { NULL, mom_vlasov_hamil_vel_five_moments_1x2v_ser_p1, mom_vlasov_hamil_vel_five_moments_1x2v_ser_p2, NULL }, // 1
  { NULL, mom_vlasov_hamil_vel_five_moments_1x3v_ser_p1, mom_vlasov_hamil_vel_five_moments_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, mom_vlasov_hamil_vel_five_moments_2x1v_ser_p1, mom_vlasov_hamil_vel_five_moments_2x1v_ser_p2, NULL }, // 3
  { NULL, mom_vlasov_hamil_vel_five_moments_2x2v_ser_p1, mom_vlasov_hamil_vel_five_moments_2x2v_ser_p2, NULL }, // 4
  { NULL, mom_vlasov_hamil_vel_five_moments_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, mom_vlasov_hamil_vel_five_moments_3x3v_ser_p1, NULL, NULL }, // 6
};

// Five moments (Zeroth, First, and Second moment together) kernel list (Tensor basis)
// {1, dH/dv, H} moments for velocity-space Hamiltonian. 
GKYL_CU_D
static const gkyl_mom_kern_list tensor_five_moments_hamil_vel_kernels[] = {
  // 1x kernels
  { NULL, mom_vlasov_hamil_vel_five_moments_1x1v_tensor_p1, mom_vlasov_hamil_vel_five_moments_1x1v_tensor_p2, NULL }, // 0
  { NULL, mom_vlasov_hamil_vel_five_moments_1x2v_tensor_p1, mom_vlasov_hamil_vel_five_moments_1x2v_tensor_p2, NULL }, // 1
  { NULL, mom_vlasov_hamil_vel_five_moments_1x3v_tensor_p1, mom_vlasov_hamil_vel_five_moments_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, mom_vlasov_hamil_vel_five_moments_2x1v_tensor_p1, mom_vlasov_hamil_vel_five_moments_2x1v_tensor_p2, NULL }, // 3
  { NULL, mom_vlasov_hamil_vel_five_moments_2x2v_tensor_p1, mom_vlasov_hamil_vel_five_moments_2x2v_tensor_p2, NULL }, // 4
  { NULL, mom_vlasov_hamil_vel_five_moments_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, mom_vlasov_hamil_vel_five_moments_3x3v_tensor_p1, NULL, NULL }, // 6
};

// Five moments (Zeroth, First, and Second moment together) kernel list (Serendipity basis)
// {1, dH/dv, H} moments for phase-space Hamiltonian. 
GKYL_CU_D
static const gkyl_mom_kern_list ser_five_moments_hamil_gen_kernels[] = {
  // 1x kernels
  { NULL, mom_vlasov_hamil_gen_five_moments_1x1v_ser_p1, mom_vlasov_hamil_gen_five_moments_1x1v_ser_p2, NULL }, // 0
  { NULL, mom_vlasov_hamil_gen_five_moments_1x2v_ser_p1, mom_vlasov_hamil_gen_five_moments_1x2v_ser_p2, NULL }, // 1
  { NULL, mom_vlasov_hamil_gen_five_moments_1x3v_ser_p1, mom_vlasov_hamil_gen_five_moments_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, mom_vlasov_hamil_gen_five_moments_2x1v_ser_p1, mom_vlasov_hamil_gen_five_moments_2x1v_ser_p2, NULL }, // 3
  { NULL, mom_vlasov_hamil_gen_five_moments_2x2v_ser_p1, mom_vlasov_hamil_gen_five_moments_2x2v_ser_p2, NULL }, // 4
  { NULL, mom_vlasov_hamil_gen_five_moments_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, mom_vlasov_hamil_gen_five_moments_3x3v_ser_p1, NULL, NULL }, // 6
};

// Five moments (Zeroth, First, and Second moment together) kernel list (Tensor basis)
// {1, dH/dv, H} moments for phase-space Hamiltonian. 
GKYL_CU_D
static const gkyl_mom_kern_list tensor_five_moments_hamil_gen_kernels[] = {
  // 1x kernels
  { NULL, mom_vlasov_hamil_gen_five_moments_1x1v_tensor_p1, mom_vlasov_hamil_gen_five_moments_1x1v_tensor_p2, NULL }, // 0
  { NULL, mom_vlasov_hamil_gen_five_moments_1x2v_tensor_p1, mom_vlasov_hamil_gen_five_moments_1x2v_tensor_p2, NULL }, // 1
  { NULL, mom_vlasov_hamil_gen_five_moments_1x3v_tensor_p1, mom_vlasov_hamil_gen_five_moments_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, mom_vlasov_hamil_gen_five_moments_2x1v_tensor_p1, mom_vlasov_hamil_gen_five_moments_2x1v_tensor_p2, NULL }, // 3
  { NULL, mom_vlasov_hamil_gen_five_moments_2x2v_tensor_p1, mom_vlasov_hamil_gen_five_moments_2x2v_tensor_p2, NULL }, // 4
  { NULL, mom_vlasov_hamil_gen_five_moments_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, mom_vlasov_hamil_gen_five_moments_3x3v_tensor_p1, NULL, NULL }, // 6
};

// Integrated moments {1, dH/dv, H}, velocity-space Hamiltonian, kernel list (Serendipity basis)
GKYL_CU_D
static const gkyl_mom_kern_list ser_int_five_moments_hamil_vel_kernels[] = {
  // 1x kernels
  { NULL, mom_vlasov_hamil_vel_int_five_moments_1x1v_ser_p1, mom_vlasov_hamil_vel_int_five_moments_1x1v_ser_p2, NULL }, // 0
  { NULL, mom_vlasov_hamil_vel_int_five_moments_1x2v_ser_p1, mom_vlasov_hamil_vel_int_five_moments_1x2v_ser_p2, NULL }, // 1
  { NULL, mom_vlasov_hamil_vel_int_five_moments_1x3v_ser_p1, mom_vlasov_hamil_vel_int_five_moments_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, mom_vlasov_hamil_vel_int_five_moments_2x1v_ser_p1, mom_vlasov_hamil_vel_int_five_moments_2x1v_ser_p2, NULL }, // 3
  { NULL, mom_vlasov_hamil_vel_int_five_moments_2x2v_ser_p1, mom_vlasov_hamil_vel_int_five_moments_2x2v_ser_p2, NULL }, // 4
  { NULL, mom_vlasov_hamil_vel_int_five_moments_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, mom_vlasov_hamil_vel_int_five_moments_3x3v_ser_p1, NULL, NULL }, // 6
};

// Integrated moments {1, dH/dv, H}, velocity-space Hamiltonian, kernel list (Tensor basis)
GKYL_CU_D
static const gkyl_mom_kern_list tensor_int_five_moments_hamil_vel_kernels[] = {
  // 1x kernels
  { NULL, mom_vlasov_hamil_vel_int_five_moments_1x1v_tensor_p1, mom_vlasov_hamil_vel_int_five_moments_1x1v_tensor_p2, NULL }, // 0
  { NULL, mom_vlasov_hamil_vel_int_five_moments_1x2v_tensor_p1, mom_vlasov_hamil_vel_int_five_moments_1x2v_tensor_p2, NULL }, // 1
  { NULL, mom_vlasov_hamil_vel_int_five_moments_1x3v_tensor_p1, mom_vlasov_hamil_vel_int_five_moments_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, mom_vlasov_hamil_vel_int_five_moments_2x1v_tensor_p1, mom_vlasov_hamil_vel_int_five_moments_2x1v_tensor_p2, NULL }, // 3
  { NULL, mom_vlasov_hamil_vel_int_five_moments_2x2v_tensor_p1, mom_vlasov_hamil_vel_int_five_moments_2x2v_tensor_p2, NULL }, // 4
  { NULL, mom_vlasov_hamil_vel_int_five_moments_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, mom_vlasov_hamil_vel_int_five_moments_3x3v_tensor_p1, NULL, NULL }, // 6
};

// Integrated moments {1, dH/dv, H}, velocity-space Hamiltonian, kernel list (Serendipity basis)
GKYL_CU_D
static const gkyl_mom_kern_list ser_int_five_moments_hamil_gen_kernels[] = {
  // 1x kernels
  { NULL, mom_vlasov_hamil_gen_int_five_moments_1x1v_ser_p1, mom_vlasov_hamil_gen_int_five_moments_1x1v_ser_p2, NULL }, // 0
  { NULL, mom_vlasov_hamil_gen_int_five_moments_1x2v_ser_p1, mom_vlasov_hamil_gen_int_five_moments_1x2v_ser_p2, NULL }, // 1
  { NULL, mom_vlasov_hamil_gen_int_five_moments_1x3v_ser_p1, mom_vlasov_hamil_gen_int_five_moments_1x3v_ser_p2, NULL }, // 2
  // 2x kernels
  { NULL, mom_vlasov_hamil_gen_int_five_moments_2x1v_ser_p1, mom_vlasov_hamil_gen_int_five_moments_2x1v_ser_p2, NULL }, // 3
  { NULL, mom_vlasov_hamil_gen_int_five_moments_2x2v_ser_p1, mom_vlasov_hamil_gen_int_five_moments_2x2v_ser_p2, NULL }, // 4
  { NULL, mom_vlasov_hamil_gen_int_five_moments_2x3v_ser_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, mom_vlasov_hamil_gen_int_five_moments_3x3v_ser_p1, NULL, NULL }, // 6
};

// Integrated moments {1, dH/dv, H}, velocity-space Hamiltonian, kernel list (Tensor basis)
GKYL_CU_D
static const gkyl_mom_kern_list tensor_int_five_moments_hamil_gen_kernels[] = {
  // 1x kernels
  { NULL, mom_vlasov_hamil_gen_int_five_moments_1x1v_tensor_p1, mom_vlasov_hamil_gen_int_five_moments_1x1v_tensor_p2, NULL }, // 0
  { NULL, mom_vlasov_hamil_gen_int_five_moments_1x2v_tensor_p1, mom_vlasov_hamil_gen_int_five_moments_1x2v_tensor_p2, NULL }, // 1
  { NULL, mom_vlasov_hamil_gen_int_five_moments_1x3v_tensor_p1, mom_vlasov_hamil_gen_int_five_moments_1x3v_tensor_p2, NULL }, // 2
  // 2x kernels
  { NULL, mom_vlasov_hamil_gen_int_five_moments_2x1v_tensor_p1, mom_vlasov_hamil_gen_int_five_moments_2x1v_tensor_p2, NULL }, // 3
  { NULL, mom_vlasov_hamil_gen_int_five_moments_2x2v_tensor_p1, mom_vlasov_hamil_gen_int_five_moments_2x2v_tensor_p2, NULL }, // 4
  { NULL, mom_vlasov_hamil_gen_int_five_moments_2x3v_tensor_p1, NULL, NULL }, // 5
  // 3x kernels
  { NULL, mom_vlasov_hamil_gen_int_five_moments_3x3v_tensor_p1, NULL, NULL }, // 6
};

struct mom_type_vlasov {
  struct gkyl_mom_type momt;
  int hamil_dim; // Dimensionality of Hamiltonian. 
  int hamil_offset; // Offset for indexing Hamiltonian from phase-space index. 
  struct gkyl_range hamil_range; // Range for indexing Hamiltonian (either velocity-space range or full phase-space range).
  const struct gkyl_array *hamil; // Hamiltonian utilized to compute certain moments. 
  vlasov_momf_t kernel; // moment calculation kernel
};

/**
 * Free moment object.
 *
 * @param ref Reference counter for moment to free
 */
void gkyl_mom_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static void
kernel(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov *mom_vlasov = container_of(momt, struct mom_type_vlasov, momt);
  
  int idx_hamil[GKYL_MAX_DIM];
  for (int i=0; i<mom_vlasov->hamil_dim; ++i) {
    idx_hamil[i] = idx[mom_vlasov->hamil_offset+i];
  }
  long hidx = gkyl_range_idx(&mom_vlasov->hamil_range, idx_hamil);  

  return mom_vlasov->kernel(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vlasov->hamil, hidx), f, out);
}
