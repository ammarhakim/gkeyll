GKYL_CU_DH static inline void 
hyb_2x1v_p1_xdir_upwind_quad_to_modal(const double* fUpwindQuad, double* GKYL_RESTRICT fUpwind) { 
  fUpwind[0] = 0.2777777777777778*fUpwindQuad[5]+0.4444444444444444*fUpwindQuad[4]+0.2777777777777778*fUpwindQuad[3]+0.2777777777777778*fUpwindQuad[2]+0.4444444444444444*fUpwindQuad[1]+0.2777777777777778*fUpwindQuad[0]; 
  fUpwind[1] = 0.2777777777777778*fUpwindQuad[5]+0.4444444444444444*fUpwindQuad[4]+0.2777777777777778*fUpwindQuad[3]-0.2777777777777778*fUpwindQuad[2]-0.4444444444444444*fUpwindQuad[1]-0.2777777777777778*fUpwindQuad[0]; 
  fUpwind[2] = 0.3726779962499649*fUpwindQuad[5]-0.3726779962499649*fUpwindQuad[3]+0.3726779962499649*fUpwindQuad[2]-0.3726779962499649*fUpwindQuad[0]; 
  fUpwind[3] = 0.3726779962499649*fUpwindQuad[5]-0.3726779962499649*fUpwindQuad[3]-0.3726779962499649*fUpwindQuad[2]+0.3726779962499649*fUpwindQuad[0]; 
  fUpwind[4] = 0.2484519974999766*fUpwindQuad[5]-0.4969039949999532*fUpwindQuad[4]+0.2484519974999766*fUpwindQuad[3]+0.2484519974999766*fUpwindQuad[2]-0.4969039949999532*fUpwindQuad[1]+0.2484519974999766*fUpwindQuad[0]; 
  fUpwind[5] = 0.2484519974999767*fUpwindQuad[5]-0.4969039949999535*fUpwindQuad[4]+0.2484519974999767*fUpwindQuad[3]-0.2484519974999767*fUpwindQuad[2]+0.4969039949999535*fUpwindQuad[1]-0.2484519974999767*fUpwindQuad[0]; 

} 

GKYL_CU_DH static inline void 
hyb_2x1v_p1_vdir_upwind_quad_to_modal(const double* fUpwindQuad, double* GKYL_RESTRICT fUpwind) { 
  fUpwind[0] = 0.5*fUpwindQuad[3]+0.5*fUpwindQuad[2]+0.5*fUpwindQuad[1]+0.5*fUpwindQuad[0]; 
  fUpwind[1] = 0.5*fUpwindQuad[3]+0.5*fUpwindQuad[2]-0.5*fUpwindQuad[1]-0.5*fUpwindQuad[0]; 
  fUpwind[2] = 0.5*fUpwindQuad[3]-0.5*fUpwindQuad[2]+0.5*fUpwindQuad[1]-0.5*fUpwindQuad[0]; 
  fUpwind[3] = 0.5*fUpwindQuad[3]-0.5*fUpwindQuad[2]-0.5*fUpwindQuad[1]+0.5*fUpwindQuad[0]; 

} 

