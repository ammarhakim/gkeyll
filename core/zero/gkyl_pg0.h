/* gkyl_pg0.h — the postgkyl support API (the "pg0 shim").
 *
 * A small, stable C surface over the core library for the postgkyl
 * post-processing tool. It deliberately exposes NO gkyl types: only opaque
 * handles, scalars, and plain buffers, so a consumer compiled against this
 * header holds no knowledge of struct layouts or calling conventions.
 * Everything layout- and ABI-sensitive lives in zero/pg0.c, compiled into
 * libg0core.so — so shim and library can never drift apart, and any core
 * API change that affects the shim fails THIS build, at the producer.
 *
 * Consumers (postgkyl's _g0py CPython extension) include only this header
 * and link -lg0core. Bump PG0_API_VERSION when this header changes shape;
 * consumers check it at load time to reject a stale pairing.
 */
#pragma once

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Bump when this header changes shape; checked by consumers at load time. */
#define PG0_API_VERSION 1
int pg0_api_version(void);

/* Enough for GKYL_MAX_DIM; fixed here so callers can size buffers without
 * seeing gkyl headers. pg0.c statically asserts it matches GKYL_MAX_DIM. */
#define PG0_MAX_DIM 7

typedef struct pg0_array pg0_array; /* wraps a gkyl_array* (double data)   */
typedef struct pg0_basis pg0_basis; /* wraps a struct gkyl_basis           */

/* ---- arrays: create / copy / introspect / release ---------------------- */
pg0_array *pg0_array_new(size_t ncomp, size_t size); /* zero-initialized   */
/* Zero-copy view of caller-owned memory laid out (size x ncomp) doubles.
 * The buffer must outlive the returned handle. */
pg0_array *pg0_array_from_buff(size_t ncomp, size_t size, double *data);
pg0_array *pg0_array_clone(const pg0_array *a);
void pg0_array_release(pg0_array *a);
size_t pg0_array_ncomp(const pg0_array *a);
size_t pg0_array_size(const pg0_array *a);
double *pg0_array_data(pg0_array *a); /* contiguous (size x ncomp) buffer  */

/* ---- file I/O (gkyl_array_rio) ------------------------------------------
 * Grid out-parameters are caller-allocated buffers of length PG0_MAX_DIM.
 * Status codes are gkyl_array_rio_status values; 0 == success.            */
int pg0_file_type(const char *fname); /* gkyl file type, or -1 */
int pg0_read_header(const char *fname, int *ndim, double *lower,
    double *upper, int *cells, int *file_type, size_t *esznc,
    size_t *tot_cells, char **meta, size_t *meta_sz); /* *meta: malloc'd   */
void pg0_meta_release(char *meta);
pg0_array *pg0_read_field(const char *fname, int *ndim, double *lower,
    double *upper, int *cells); /* NULL on failure */
const char *pg0_status_msg(int status);

/* ---- basis objects ------------------------------------------------------
 * type: "serendipity" | "tensor" (the long names stored in .gkyl metadata).
 * Returns NULL for unknown types.                                          */
pg0_basis *pg0_basis_new(const char *type, int ndim, int poly_order);
void pg0_basis_release(pg0_basis *b);
int pg0_basis_ndim(const pg0_basis *b);
int pg0_basis_poly_order(const pg0_basis *b);
int pg0_basis_num_basis(const pg0_basis *b);
const char *pg0_basis_id(const pg0_basis *b);
/* Evaluate all num_basis functions at reference point z[ndim].            */
void pg0_basis_eval(const pg0_basis *b, const double *z, double *bvals);
/* Node coordinates: coords[num_basis * ndim].                              */
void pg0_basis_node_list(const pg0_basis *b, double *coords);
/* Exact nodal -> modal change of basis on one cell's num_basis values.     */
void pg0_basis_nodal_to_modal(const pg0_basis *b, const double *fnodal,
    double *fmodal);

/* ---- weak (DG) algebra --------------------------------------------------
 * Operands must have ncomp == nfields * num_basis; the per-field loop runs
 * here. Returns 0 on success, nonzero on shape mismatch.                   */
int pg0_dg_mul(const pg0_basis *b, pg0_array *out, const pg0_array *a1,
    const pg0_array *a2);
int pg0_dg_div(const pg0_basis *b, pg0_array *out, const pg0_array *a1,
    const pg0_array *a2);
int pg0_dg_inv(const pg0_basis *b, pg0_array *out, const pg0_array *a1);

/* ---- linear coefficient ops / reductions ------------------------------- */
void pg0_array_set(pg0_array *out, double c, const pg0_array *a);
void pg0_array_accumulate(pg0_array *out, double c, const pg0_array *a);
void pg0_array_scale(pg0_array *a, double c);
void pg0_array_shiftc(pg0_array *a, double val, unsigned comp);
/* op: 0 min, 1 max, 2 sum; out[ncomp].                                     */
void pg0_array_reduce(double *out, const pg0_array *a, int op);

/* ---- integration (gkyl_array_integrate) ---------------------------------
 * op: 0 none, 1 abs, 2 sq. out[nfields]. Returns 0 on success, nonzero if
 * the grid does not cover the array.                                        */
int pg0_array_integrate(int ndim, const double *lower, const double *upper,
    const int *cells, const pg0_basis *b, int nfields, int op, double factor,
    const pg0_array *a, double *out);

#ifdef __cplusplus
}
#endif
