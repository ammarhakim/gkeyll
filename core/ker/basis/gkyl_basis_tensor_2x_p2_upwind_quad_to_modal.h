GKYL_CU_DH static inline void 
tensor_2x_p2_upwind_quad_to_modal(const double* fUpwindQuad, double* GKYL_RESTRICT fUpwind) { 
  fUpwind[0] = 0.392837100659193*fUpwindQuad[2]+0.6285393610547091*fUpwindQuad[1]+0.392837100659193*fUpwindQuad[0]; 
  fUpwind[1] = 0.5270462766947298*fUpwindQuad[2]-0.5270462766947298*fUpwindQuad[0]; 
  fUpwind[2] = 0.3513641844631533*fUpwindQuad[2]-0.7027283689263066*fUpwindQuad[1]+0.3513641844631533*fUpwindQuad[0]; 

} 
