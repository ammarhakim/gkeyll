#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include <gkyl_util.h>
#include <gkyl_alloc.h>

#include <mpack.h>

int
gkyl_search_str_int_pair_by_str(const struct gkyl_str_int_pair pairs[], const char *str, int def)
{
  for (int i=0; pairs[i].str != 0; ++i) {
    if (strcmp(pairs[i].str, str) == 0)
      return pairs[i].val;
  }
  return def;  
}

const char *
gkyl_search_str_int_pair_by_int(const struct gkyl_str_int_pair pairs[], int val, const char *def)
{
  for (int i=0; pairs[i].str != 0; ++i) {
    if (pairs[i].val == val)
      return pairs[i].str;
  }
  return def;  
}

int
gkyl_tm_trigger_check_and_bump(struct gkyl_tm_trigger *tmt, double tcurr)
{
  int status = 0;
  if (tcurr >= tmt->tcurr) {
    status = 1;
    tmt->curr += 1;
    tmt->tcurr += tmt->dt;
  }
  return status;
}

void
gkyl_exit(const char* msg)
{
  fprintf(stderr, "Error: %s\n", msg);
  exit(EXIT_FAILURE);
}

int
gkyl_compare_float(float a, float b, float eps)
{
  //if (isnanf(a) || isnanf(b)) return 0;
  
  float absa = fabs(a), absb = fabs(b), diff = fabs(a-b);

  if (a == b) return 1;
  if (a == 0 || b == 0 || (absa+absb < FLT_MIN)) return diff < eps;
  if (absa < eps) return diff < eps;
  if (absb < eps) return diff < eps;
  return diff/fminf(absa+absb, FLT_MAX) < eps;
}

int
gkyl_compare_double(double a, double b, double eps)
{
  if (isnan(a) || isnan(b)) return 0;
  
  double absa = fabs(a), absb = fabs(b), diff = fabs(a-b);
  if (a == b) return 1;
  if (a == 0 || b == 0 || (absa+absb < DBL_MIN)) return diff < eps;
  if (absa < eps) return diff < eps;
  if (absb < eps) return diff < eps;
  return diff/fmin(absa+absb, DBL_MAX) < eps;
}

struct timespec
gkyl_wall_clock(void)
{
  struct timespec tm = { 0 };
#ifdef GKYL_HAVE_CUDA
  cudaDeviceSynchronize();
#endif
  // we were using CLOCK_REALTIME here
  clock_gettime(CLOCK_MONOTONIC, &tm);
  return tm;
}

struct timespec
gkyl_time_diff(struct timespec start, struct timespec end)
{
  struct timespec tm;
  if ((end.tv_nsec-start.tv_nsec)<0) {
    tm.tv_sec = end.tv_sec-start.tv_sec-1;
    tm.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
  }
  else {
    tm.tv_sec = end.tv_sec-start.tv_sec;
    tm.tv_nsec = end.tv_nsec-start.tv_nsec;
  }
  return tm;  
}

double
gkyl_time_diff_now_sec(struct timespec tm)
{
  return gkyl_time_sec(gkyl_time_diff(tm, gkyl_wall_clock()));
}
   
double
gkyl_time_sec(struct timespec tm)
{
  return tm.tv_sec + 1e-9*tm.tv_nsec;
}

double
gkyl_time_now(void)
{
  return gkyl_time_sec( gkyl_wall_clock() );
}

pcg32_random_t
gkyl_pcg32_init(bool nd_seed)
{
  pcg32_random_t rng;
  int rounds = 5;

  if (nd_seed)
    // seed with external entropy -- the time and some program addresses
    // (which will actually be somewhat random on most modern systems).
    pcg32_srandom_r(&rng, time(NULL) ^ (intptr_t)&printf, 
      (intptr_t)&rounds);
  else
    // seed with a fixed constant
    pcg32_srandom_r(&rng, 42u, 54u);

  return rng;
}

uint32_t
gkyl_pcg32_rand_uint32(pcg32_random_t* rng)
{
  return pcg32_random_r(rng);
}

double
gkyl_pcg32_rand_double(pcg32_random_t* rng)
{
  return ldexp(pcg32_random_r(rng), -32);
}

static void
pcg64_srandom_r(pcg64_random_t* rng, uint64_t seed1, uint64_t seed2,
  uint64_t seq1,  uint64_t seq2)
{
  uint64_t mask = ~0ull >> 1;
  // stream for each generators *must* be distinct
  if ((seq1 & mask) == (seq2 & mask)) 
    seq2 = ~seq2;
  pcg32_srandom_r(rng->gen,   seed1, seq1);
  pcg32_srandom_r(rng->gen+1, seed2, seq2);
}

static int _dummy_global = 0; // just to provide address for use in seed

pcg64_random_t
gkyl_pcg64_init(bool nd_seed)
{
  pcg64_random_t rng;
  int rounds = 5;

  if (nd_seed)
    pcg64_srandom_r(&rng,
      time(NULL) ^ (intptr_t)&printf, ~time(NULL) ^ (intptr_t)&pcg32_random_r,
      (intptr_t)&rounds, (intptr_t)&_dummy_global);
  else
    pcg64_srandom_r(&rng, 42u, 42u, 54u, 54u);

  return rng;
}

uint64_t
gkyl_pcg64_rand_uint64(pcg64_random_t* rng)
{
  return ((uint64_t)(pcg32_random_r(rng->gen)) << 32) | pcg32_random_r(rng->gen+1);
}

double
gkyl_pcg64_rand_double(pcg64_random_t* rng)
{
  return ldexp(gkyl_pcg64_rand_uint64(rng), -64);
}

bool
gkyl_check_file_exists(const char *fname)
{
  return access(fname, F_OK) == 0;
}

int64_t
gkyl_file_size(const char *fname)
{
  struct stat st;
  stat(fname, &st);
  return st.st_size;
}

char*
gkyl_load_file(const char *fname, int64_t *sz)
{
  int64_t msz = gkyl_file_size(fname);
  char *buff = gkyl_malloc(msz);
  FILE *fp = fopen(fname, "r");
  int n = fread(buff, msz, 1, fp);
  *sz = msz;
  fclose(fp);
  return buff;
}

struct gkyl_msgpack_data *
gkyl_msgpack_create(int nvals, const struct gkyl_msgpack_map_elem *elist)
{
  struct gkyl_msgpack_data *mdata = gkyl_malloc(sizeof *mdata);
  mdata->meta_sz = 0;
  mdata->meta = 0;

  mpack_writer_t writer;
  mpack_writer_init_growable(&writer, &mdata->meta, &mdata->meta_sz);

  mpack_build_map(&writer);  
  for (int i=0; i<nvals; ++i) {
    mpack_write_cstr(&writer, elist[i].key);

    switch (elist[i].elem_type) {
      case GKYL_MP_INT:
        mpack_write_i64(&writer, elist[i].ival);
        break;

      case GKYL_MP_DOUBLE:
        mpack_write_double(&writer, elist[i].dval);
        break;

      case GKYL_MP_STRING:
        mpack_write_cstr(&writer, elist[i].cval);
        break;
    }
  }

  mpack_complete_map(&writer);

  int status = mpack_writer_destroy(&writer);

  if (status != mpack_ok) {
    MPACK_FREE(mdata->meta); // we need to use free here as mpack does its own malloc
    gkyl_free(mdata);
    mdata = 0;
  }

  return mdata;
}

void
gkyl_msgpack_data_release(struct gkyl_msgpack_data *mdata)
{
  if (!mdata) return;
  if (mdata->meta_sz > 0)
    MPACK_FREE(mdata->meta);
  gkyl_free(mdata);
}

