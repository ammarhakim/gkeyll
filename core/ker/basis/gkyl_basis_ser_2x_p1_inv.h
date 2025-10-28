GKYL_CU_DH static inline void 
ser_2x_p1_inv(const double *A, double *A_inv) 
{ 
  // A:     Input DG field. 
  // A_inv: Output DG field (expansion of 1/A). 
 
  const double A0R2 = pow(A[0],2);
  const double A0R3 = pow(A[0],3);
  const double A0R4 = pow(A[0],4);
  const double A1R2 = pow(A[1],2);
  const double A1R3 = pow(A[1],3);
  const double A1R4 = pow(A[1],4);
  const double A2R2 = pow(A[2],2);
  const double A2R3 = pow(A[2],3);
  const double A2R4 = pow(A[2],4);
  const double A3R2 = pow(A[3],2);
  const double A3R3 = pow(A[3],3);
  const double A3R4 = pow(A[3],4);

  double det = 0.0625*(A3R4+((-2.0*A2R2)-2.0*A1R2-2.0*A0R2)*A3R2+8.0*A[0]*A[1]*A[2]*A[3]+A2R4+((-2.0*A1R2)-2.0*A0R2)*A2R2+A1R4-2.0*A0R2*A1R2+A0R4); 
 
  A_inv[0] = -(0.25*(A[0]*A3R2-2.0*A[1]*A[2]*A[3]+A[0]*A2R2+A[0]*A1R2-1.0*A0R3))/det; 
  A_inv[1] = -(0.25*(A[1]*A3R2-2.0*A[0]*A[2]*A[3]+A[1]*A2R2-1.0*A1R3+A0R2*A[1]))/det; 
  A_inv[2] = -(0.25*(A[2]*A3R2-2.0*A[0]*A[1]*A[3]-1.0*A2R3+(A1R2+A0R2)*A[2]))/det; 
  A_inv[3] = (0.25*(A3R3+((-1.0*A2R2)-1.0*A1R2-1.0*A0R2)*A[3]+2.0*A[0]*A[1]*A[2]))/det; 
} 
 
