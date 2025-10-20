GKYL_CU_DH static inline void 
ser_1x_p1_inv(const double *A, double *A_inv) 
{ 
  // A:     Input DG field. 
  // A_inv: Output DG field (expansion of 1/A). 
 
  const double A0R2 = pow(A[0],2);
  const double A1R2 = pow(A[1],2);

  double det = -0.5*(A1R2-1.0*A0R2); 
 
  A_inv[0] = A[0]/det; 
  A_inv[1] = -(1.0*A[1])/det; 
} 
 
