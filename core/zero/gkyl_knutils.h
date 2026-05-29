#pragma once

#include <gkyl_ref_count.h>

// Struct to hold data passed to KANN functions: nvec float vectors
// with N elements each
struct gkyl_kn_vec {
  int nvec; // number of vectors
  int N; // length of each vector
  float **vals; // for passing into kann functions
  float *data; // actual data (do not use directly)

  struct gkyl_ref_count ref_count;
};

/**
 * Create a new kn_vec to hold @a nvec vectors, each of size @a N.
 *
 * @param nvec Number of vectors to create
 * @param N Length of each vector.
 * @return New kn_vec
 */
struct gkyl_kn_vec* gkyl_kn_vec_new(int nvec, int N);

/**
 * Acquire pointer to kn_vec.
 *
 * @param vec Vector to which reference is needed
 * @return Pointer to acquired array
 */
struct gkyl_kn_vec* gkyl_kn_vec_acquire(const struct gkyl_kn_vec* vec);

/**
 * Release kn_vec.
 *
 * @param vec Vector to release.
 */
void gkyl_kn_vec_release(struct gkyl_kn_vec *vec);
