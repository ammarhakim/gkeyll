#include <gkyl_alloc.h>
#include <gkyl_knutils.h>

#include <stdio.h>
#include <stdlib.h>

static void
kn_vec_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_kn_vec *vec = container_of(ref, struct gkyl_kn_vec, ref_count);
  
  gkyl_free(vec->vals);
  gkyl_free(vec->data);
  gkyl_free(vec);
}

struct
gkyl_kn_vec* gkyl_kn_vec_new(int nvec, int N)
{
  struct gkyl_kn_vec *vec = gkyl_malloc(sizeof(*vec));
  vec->nvec = nvec;
  vec->N = N;

  vec->vals = gkyl_malloc(nvec*sizeof(float*));
  vec->data = gkyl_calloc(nvec*N, sizeof(float));

  vec->vals[0] = &vec->data[0];
  for (int i=1; i<nvec; ++i)
    vec->vals[i] = vec->vals[i-1] + N;

  vec->ref_count = gkyl_ref_count_init(kn_vec_free);

  return vec;
}

struct gkyl_kn_vec *
gkyl_kn_vec_acquire(const struct gkyl_kn_vec* vec)
{
  gkyl_ref_count_inc(&vec->ref_count);
  return (struct gkyl_kn_vec*) vec;  
}
    
void
gkyl_kn_vec_release(struct gkyl_kn_vec *vec)
{
  gkyl_ref_count_dec(&vec->ref_count);
}
