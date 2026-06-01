#pragma once

// Struct to hold data passed to KANN functions: nvec float vectors
// with N elements each
struct kn_vec {
  int nvec; // number of vectors
  int N; // length of each vector
  float **vals; // for passing into kann functions
  float *data; // actual data (do not use directly)
};

/**
 * Create a new kn_vec to hold @a nvec vectors, each of size @a N.
 *
 * @param nvec Number of vectors to create
 * @param N Length of each vector.
 * @return New kn_vec
 */
struct kn_vec* kn_vec_new(int nvec, int N);

/**
 * Release kn_vec.
 *
 * @param vec Vector to release.
 */
void kn_vec_release(struct kn_vec *vec);
