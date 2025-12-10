#pragma once

#include <gkyl_wv_eqn.h>

// input packaged as a struct
struct gkyl_wv_coldfluid_inp {
  struct gkyl_wv_embed_geo *embed_geo; // embedded geometry
};
/**
 * Create a new cold-fluid equation object.
 * 
 * @return Pointer to cold-fluid equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_coldfluid_new(void);
/**
 * Create a new cold-fluid equation object.
 * 
 * @param inp Input parameters
 * @return Pointer to cold-fluid equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_coldfluid_inew(const struct gkyl_wv_coldfluid_inp *inp);

