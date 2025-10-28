GKYL_CU_DH static inline void 
ser_1x_p2_inv(const double *A, double *A_inv) 
{ 
  // A:     Input DG field. 
  // A_inv: Output DG field (expansion of 1/A). 
 
  const double A0R2 = pow(A[0],2);
  const double A0R3 = pow(A[0],3);
  const double A1R2 = pow(A[1],2);
  const double A1R3 = pow(A[1],3);
  const double A2R2 = pow(A[2],2);
  const double A2R3 = pow(A[2],3);

  double det = (-0.3162277660168379*A2R3)-0.1515228816828315*A[0]*A2R2+0.4065785563073631*A1R2*A[2]+0.5421047417431506*A0R2*A[2]-0.6363961030678927*A[0]*A1R2+0.3535533905932737*A0R3; 
 
  A_inv[0] = (0.4040610178208844*A2R2+1.084209483486302*A[0]*A[2]-0.5656854249492381*A1R2+0.7071067811865475*A0R2)/det; 
  A_inv[1] = (0.1807015805810502*A[1]*A[2]-0.7071067811865475*A[0]*A[1])/det; 
  A_inv[2] = ((-0.6324555320336759*A2R2)-0.7071067811865475*A[0]*A[2]+0.6324555320336759*A1R2)/det; 
} 
 
