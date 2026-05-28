#include <knutils.h>

#include <stdio.h>
#include <stdlib.h>

struct
kn_vec* kn_vec_new(int nvec, int N)
{
  struct kn_vec *vec = malloc(sizeof(*vec));
  vec->nvec = nvec;
  vec->N = N;

  vec->vals = malloc(nvec*sizeof(float*));
  vec->data = calloc(nvec*N, sizeof(float));

  vec->vals[0] = &vec->data[0];
  for (int i=1; i<nvec; ++i)
    vec->vals[i] = vec->vals[i-1] + N;

  return vec;
}

void
kn_vec_release(struct kn_vec *vec)
{
  free(vec->vals);
  free(vec->data);
  free(vec);
}
