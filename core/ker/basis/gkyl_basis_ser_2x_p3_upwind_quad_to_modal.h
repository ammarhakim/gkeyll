GKYL_CU_DH static inline void 
ser_2x_p3_upwind_quad_to_modal(const double* fUpwindQuad, double* GKYL_RESTRICT fUpwind) { 
  fUpwind[0] = 0.2459705198652899*fUpwindQuad[3]+0.4611362613212574*fUpwindQuad[2]+0.4611362613212574*fUpwindQuad[1]+0.2459705198652899*fUpwindQuad[0]; 
  fUpwind[1] = 0.366872863045464*fUpwindQuad[3]+0.2715467467935446*fUpwindQuad[2]-0.2715467467935446*fUpwindQuad[1]-0.366872863045464*fUpwindQuad[0]; 
  fUpwind[2] = 0.3367876570272815*fUpwindQuad[3]-0.3367876570272815*fUpwindQuad[2]-0.3367876570272815*fUpwindQuad[1]+0.3367876570272815*fUpwindQuad[0]; 
  fUpwind[3] = 0.1983222754244991*fUpwindQuad[3]-0.5023295150965306*fUpwindQuad[2]+0.5023295150965306*fUpwindQuad[1]-0.1983222754244991*fUpwindQuad[0]; 

} 
