/* gkyl_gpython.h — the postgkyl support API (the "gpython shim").
 *
 * A small, stable C surface over the core library for the postgkyl
 * post-processing tool. It deliberately exposes NO gkyl types: only opaque
 * handles, scalars, and plain buffers, so a consumer compiled against this
 * header holds no knowledge of struct layouts or calling conventions.
 * Everything layout- and ABI-sensitive lives in zero/gpython.c, compiled into
 * libg0core.so — so shim and library can never drift apart, and any core
 * API change that affects the shim fails THIS build, at the producer.
 *
 * Consumers (postgkyl's _gpython CPython extension) include only this header
 * and link -lg0core. Bump GPYTHON_API_VERSION when this header changes shape;
 * consumers check it at load time to reject a stale pairing.
 */
#pragma once

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Bump when this header changes shape; checked by consumers at load time.
 * v2: added pg0_array_dg_reduce, pg0_write_field, pg0_dynvec_{read,write}.
 * v3: added pg0_basis_new_hybrid (hybrid/gkhybrid bases).
 * v4: added pg0_dg_mul_conf_phase (conf-space x phase-space weak product).
 * v5: added pg0_array_average (weighted/plain average over a dim subset). */
#define GPYTHON_API_VERSION 5
int pg0_api_version(void);

/* Enough for GKYL_MAX_DIM; fixed here so callers can size buffers without
 * seeing gkyl headers. gpython.c statically asserts it matches GKYL_MAX_DIM. */
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
 * Returns NULL for unknown types. Use pg0_basis_new_hybrid for the two
 * mixed-order bases below.                                                 */
pg0_basis *pg0_basis_new(const char *type, int ndim, int poly_order);
/* type: "hybrid" | "gkhybrid" — mixed-order bases with a fixed poly_order
 * of 1, parameterized by configuration/velocity dimension (ndim = cdim +
 * vdim) rather than by an explicit poly_order, mirroring Gkeyll's own
 * gkyl_cart_modal_hybrid/gkhybrid constructors exactly. Returns NULL for
 * unknown types; out-of-range cdim/vdim is a process abort in the
 * underlying constructor (see gkyl_basis.h), not a clean failure — callers
 * must respect its bounds.                                                 */
pg0_basis *pg0_basis_new_hybrid(const char *type, int cdim, int vdim);
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
/* Conf-space x phase-space weak product: pout = cop * pop, where cop lives
 * on a conf-space grid of cbasis->ndim dimensions and pop/pout live on a
 * phase-space grid of pbasis->ndim (> cbasis->ndim) dimensions sharing the
 * conf-space grid's leading dimensions (gkyl_dg_mul_conf_phase_op_range).
 * Single-field only: cop.ncomp == cbasis num_basis, pop/pout.ncomp == pbasis
 * num_basis -- the underlying kernel takes no field-index arguments.
 * conf_cells/phase_cells give each grid's per-dimension cell count (length
 * cbasis->ndim / pbasis->ndim), used to build the index ranges Gkeyll maps
 * each phase cell's conf cell through. Returns 0 on success, nonzero if the
 * arrays don't match the bases or cells don't cover the arrays.             */
int pg0_dg_mul_conf_phase(const pg0_basis *cbasis, const pg0_basis *pbasis,
    pg0_array *pout, const pg0_array *cop, const pg0_array *pop,
    const int *conf_cells, const int *phase_cells);

/* ---- linear coefficient ops / reductions ------------------------------- */
void pg0_array_set(pg0_array *out, double c, const pg0_array *a);
void pg0_array_accumulate(pg0_array *out, double c, const pg0_array *a);
void pg0_array_scale(pg0_array *a, double c);
void pg0_array_shiftc(pg0_array *a, double val, unsigned comp);
/* op: 0 min, 1 max, 2 sum; out[ncomp]. Reduces the raw DG COEFFICIENTS —
 * exact for "sum" but not the field's true min/max; see
 * pg0_array_dg_reduce for the field-aware version.                        */
void pg0_array_reduce(double *out, const pg0_array *a, int op);

/* Field-aware reduction (gkyl_array_dg_reducec): evaluates the DG expansion
 * of field `comp` (0-based, ncomp must be a multiple of num_basis) at each
 * cell's Gauss-Legendre quadrature nodes and reduces THOSE — the true
 * min/max/sum of the represented field, not the coefficients. op: 0/1/2
 * min/max/sum. Returns 0 on success, nonzero if comp is out of range.      */
int pg0_array_dg_reduce(double *out, const pg0_basis *b, const pg0_array *a,
    int comp, int op);

/* ---- integration (gkyl_array_integrate) ---------------------------------
 * op: 0 none, 1 abs, 2 sq. out[nfields]. Returns 0 on success, nonzero if
 * the grid does not cover the array.                                        */
int pg0_array_integrate(int ndim, const double *lower, const double *upper,
    const int *cells, const pg0_basis *b, int nfields, int op, double factor,
    const pg0_array *a, double *out);

/* ---- averaging (gkyl_array_average) -------------------------------------
 * Single-field weighted (or plain) average of `a` over the donor dims
 * flagged in `avg_dim` (length ndim, 1 = averaged, 0 = kept):
 *   weight != NULL:  int f(x) w(x) dx^avg / int w(x) dx^avg
 *   weight == NULL:  int f(x) dx^avg / int dx^avg
 * `b` is the full donor basis (ndim/poly_order); `b_avg` the reduced-
 * dimension (same poly_order) target basis for the surviving `ndim_avg`
 * dims -- `ndim_avg` is 1 with `cells_avg = {1}` for a full reduction
 * (Gkeyll keeps a single-cell target rather than a true 0-dim basis, same
 * convention `pg0_array_integrate`'s caller uses at the Python layer).
 * `a`/`weight`/`out` are single-field: their ncomp must equal the relevant
 * basis's num_basis exactly (gkyl_array_average takes no field-index
 * argument, like pg0_dg_mul_conf_phase; a multi-field caller loops in
 * Python). Returns 0 on success, nonzero if `a`/`weight`/`out` don't match
 * the bases or the donor/target cell counts don't cover them.               */
int pg0_array_average(int ndim, const double *lower, const double *upper,
    const int *cells, const pg0_basis *b, const pg0_basis *b_avg,
    int ndim_avg, const int *cells_avg, const int *avg_dim,
    const pg0_array *weight, const pg0_array *a, pg0_array *out);

/* ---- writing (gkyl_array_rio) --------------------------------------------
 * Mirrors pg0_read_field: writes the FULL array over a uniform grid built
 * from lower/upper/cells (the same convention pg0_array_integrate uses).
 * meta/meta_sz: a raw msgpack byte blob, or meta_sz == 0 for none — encoding
 * policy lives in postgkyl's io/ layer, exactly as decoding does for reads.
 * Returns -1 if the grid does not cover the array (checked here to avoid an
 * out-of-bounds C loop); otherwise a gkyl_array_rio_status (0 == success,
 * rendered by pg0_status_msg).                                             */
int pg0_write_field(const char *fname, int ndim, const double *lower,
    const double *upper, const int *cells, const char *meta, size_t meta_sz,
    const pg0_array *a);

/* ---- dynvector (time-series) I/O (gkyl_dynvec) ---------------------------
 * Double-precision dynvectors only (the only kind postgkyl ever produces).
 * pg0_dynvec_read allocates *tm (ncomp=1) and *data (ncomp=*ncomp), both
 * sized to the number of time samples found in the file — the caller
 * releases them like any other pg0_array. Returns 0 on success; nonzero if
 * the file is missing/unreadable or not a double-typed dynvector.          */
int pg0_dynvec_read(const char *fname, size_t *ncomp, pg0_array **tm,
    pg0_array **data);
/* Writes a dynvector built by zipping tm[n] with data[n*ncomp:(n+1)*ncomp].
 * Returns 0 on success.                                                    */
int pg0_dynvec_write(const char *fname, size_t ncomp, size_t n,
    const double *tm, const double *data);

#ifdef __cplusplus
}
#endif
