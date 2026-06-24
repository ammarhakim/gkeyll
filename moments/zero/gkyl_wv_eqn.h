#pragma once

#include <gkyl_eqn_type.h>
#include <gkyl_ref_count.h>
#include <gkyl_util.h>
#include <gkyl_evalf_def.h>

// Flux type for use in wave/qfluct methods
enum gkyl_wv_flux_type { GKYL_WV_HIGH_ORDER_FLUX, GKYL_WV_LOW_ORDER_FLUX };

// Forward declare for use in function pointers
struct gkyl_wv_eqn;

// Function pointer to compute waves from RP solver
typedef double (*wv_waves_t)(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, double phil, double phir, 
  double *waves, double *speeds);

// Function pointer to compute q-fluctuations from waves
typedef void (*wv_qfluct_t)(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, double phil, double phir, const double *waves, const double *speeds,
  double *amdq, double *apdq);

// Function pointer to compute jump in flux. Returns absolute maximum
// wave-speed
typedef double (*wv_flux_jump_t)(const struct gkyl_wv_eqn *eqn,
  const double *ql, const double *qr, double *flux_jump);

// Function pointer to check if invariant domain is preserved
typedef bool (*wv_check_inv)(const struct gkyl_wv_eqn *eqn, const double *q);

// Function pointer to compute maximum speed given local state
typedef double (*wv_max_speed_t)(const struct gkyl_wv_eqn *eqn, const double *q);
// Optional per-direction variant: max speed in coordinate direction @a dir.
// Equations that provide it get a directionally-aware initial-dt seed in
// gkyl_wave_prop_max_dt (dt ≤ cfl·dx_dir/λ_dir); others fall back to the
// direction-agnostic max_speed_func.
typedef double (*wv_max_speed_dir_t)(const struct gkyl_wv_eqn *eqn,
  const double *q, int dir);

// Function pointer to rotate conserved variables to local
// tangent-normal frame: tau1 X tau2 = norm
typedef void (*wv_rotate_to_local)(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm,
  const double* qglobal, double* qlocal);

// Function pointer to rotate conserved variables to local
// tangent-normal frame: tau1 X tau2 = norm
typedef void (*wv_rotate_to_global)(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm,
  const double* qlocal, double* qglobal);

// Function pointer to convert conserved variables to Riemann
// variables, given an input state 'qstate'
typedef void (*wv_cons_to_riem)(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *qin, double *wout);

// Function pointer to convert Riemann variables back to conserved
// variables, given an input state 'qstate'
typedef void (*wv_riem_to_cons)(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *win, double *qout);

// Function pointer to compute diagnostic variables from conserved
// variables
typedef void (*wv_cons_to_diag)(const struct gkyl_wv_eqn *eqn,
    const double *qin, double *diag);
  
// Function pointer to compute the forcing/source term vector.
typedef void (*wv_source_func_t)(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout);

// Optional function pointer returning a scalar limiter alpha in [0,1] for an
// explicit source step q_new = q_old + alpha*dt*S(q_old), so an equation can
// keep the limited update inside its admissible set (e.g. GR-Euler positivity).
// NULL means no limiting (alpha = 1).
typedef double (*wv_source_limiter_func_t)(const struct gkyl_wv_eqn* eqn,
  const double* q_old, const double* sout, double dt);

// Optional function pointers used by wave_prop (and similar drivers) to
// inform the equation object of the relevant cell index before each callback
// fires. Equations that need to fetch per-cell auxfield data implement these;
// others leave them NULL.
//
// set_interface_idx_func: called once per interface before the rotate/waves/
// qfluct/flux_jump callbacks, supplying both bracketing cells so the equation
// can keep left- and right-side auxfields available.
typedef void (*wv_set_interface_idx_t)(const struct gkyl_wv_eqn* eqn,
  const int* idxl, const int* idxr);

// set_cell_idx_func: called once per cell before single-state callbacks
// (check_inv, max_speed, source, cons_to_diag).
typedef void (*wv_set_cell_idx_t)(const struct gkyl_wv_eqn* eqn,
  const int* idx);

// repair_state_func: optional last-resort hook called by wave_prop (after
// the positivity-redo) and by moment_spacetime_coupling (after each
// forward-Euler source stage) when a cell remains outside the convex
// admissibility set. The equation projects q in-place onto the boundary
// of the admissibility set, fixing one constraint per call. q_prev is the
// cell's pre-update state — guaranteed to be in the admissibility set
// (it's the state the failing update started from). Equations that want
// history-informed repair (anchor target on the cell's last-valid
// primitives) read it; equations that don't can ignore it. NULL for
// equations that do not implement state repair.
typedef void (*wv_repair_state_t)(const struct gkyl_wv_eqn *eqn,
  const double *q_prev, double *q);

struct gkyl_wv_eqn {
  enum gkyl_eqn_type type; // Equation type
  int num_equations; // number of equations in system
  int num_waves; // number of waves in system
  int num_diag; // number of diagnostic variables

  wv_waves_t waves_func; // function to compute waves and speeds
  wv_qfluct_t qfluct_func; // function to compute q-fluctuations
  wv_qfluct_t ffluct_func; // function to compute f-fluctuations

  wv_flux_jump_t flux_jump; // function to compute jump in flux

  wv_check_inv check_inv_func; // function to check invariant domains
  wv_max_speed_t max_speed_func; // function to compute max-speed
  wv_max_speed_dir_t max_speed_dir_func; // optional per-direction max-speed
                                         // (NULL → max_speed_func fallback)
  wv_rotate_to_local rotate_to_local_func; // function to rotate to local frame
  wv_rotate_to_global rotate_to_global_func; // function to rotate to global frame

  wv_cons_to_riem cons_to_riem; // function to convert cons to Riemann vars
  wv_cons_to_riem riem_to_cons; // function to convert Riemann vars to cons

  wv_bc_func_t wall_bc_func; // function to apply wall BC
  wv_bc_func_t no_slip_bc_func; // function to apply no-slip BC

  wv_cons_to_diag cons_to_diag; // function for diagnostic variables

  wv_source_func_t source_func; // function for computing the forcing/source term vector.
  wv_source_limiter_func_t source_limiter_func; // optional positivity limiter for explicit source steps (NULL => no limiting).

  // Optional setters invoked by wave_prop (and similar drivers) before each
  // callback to communicate the relevant cell index/indices. NULL for
  // equations that do not need per-cell auxfield lookup.
  wv_set_interface_idx_t set_interface_idx_func;
  wv_set_cell_idx_t set_cell_idx_func;

  // Optional last-resort state-repair hook. NULL for equations that do
  // not implement it; callers must check before invoking.
  wv_repair_state_t repair_state_func;

  // Diagnostic context hint set by the caller before invoking
  // repair_state_func. Lets equations that maintain per-call-site repair
  // counters dispatch to the right tally. Convention:
  //   0 = source-step (or unspecified)
  //   1 = wave_prop post-positivity-sweep
  // Equations that don't track per-context stats can ignore this field.
  // Mutable because callers reach in through const struct gkyl_wv_eqn*.
  int cur_repair_ctx;

  struct gkyl_wv_embed_geo *embed_geo;

  uint32_t flags;  
  struct gkyl_ref_count ref_count; // reference count
  struct gkyl_wv_eqn *on_dev; // pointer to itself or device data
};

/**
 * Check if equation is on device.
 *
 * @param eqn Equation to check
 * @return true if eqn on device, false otherwise
 */
bool gkyl_wv_eqn_is_cu_dev(const struct gkyl_wv_eqn *eqn);

/**
 * Acquire pointer to equation object. Delete using the release()
 * method
 *
 * @param eqn Equation object.
 * @return Acquired eqn obj pointer
 */
struct gkyl_wv_eqn *gkyl_wv_eqn_acquire(const struct gkyl_wv_eqn *eqn);

/**
 * Default function to convert conserved vars to diagostics: for many
 * eqn systems the conserved vara are the diagnostics one wishes to
 * compute.
 */
GKYL_CU_DH
static inline void
gkyl_default_cons_to_diag(const struct gkyl_wv_eqn *eqn,
  const double *qin, double *diag)
{
  for (int i=0; i<eqn->num_equations; ++i) diag[i] = qin[i];
}

/**
* Default function to compute forcing/source term vector: assumes that the system of equations being solved is strictly homogeneous
* (i.e. source-free).
*
* @param eqn Base equation object.
* @param qin Conserved variable vector (input).
* @param sout Forcing/source term vector (output).
*/
GKYL_CU_DH
static inline void
gkyl_default_source_func(const struct gkyl_wv_eqn *eqn, const double *qin, double *sout)
{
  for (int i = 0; i < eqn->num_equations; i++) {
    sout[i] = 0.0;
  }
}

/**
 * Compute waves and speeds from left/right conserved variables. The
 * 'waves' array has size num_equations X num_waves in length. The 'm'
 * wave (m = 0 ... num_waves-1) is stored starting at location
 * waves[m*num_equations].
 *
 * @param eqn Equation object
 * @param delta Jump across interface to split
 * @param ql Conserved variables on left of interface
 * @param qr Conserved variables on right of interface
 * @param waves On output, waves
 * @param speeds On output wave speeds[num_wave]
 * @return Maximum wave speed.
 */
GKYL_CU_DH
static inline double
gkyl_wv_eqn_waves(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, const double phil, const double phir, 
  double *waves, double *speeds)
{
  return eqn->waves_func(eqn, type, delta, ql, qr, phil, phir, waves, speeds);
}

/**
 * Compute waves and speeds from left/right conserved variables. The
 * 'waves' array has size num_equations X num_waves in length. The 'm'
 * wave (m = 0 ... num_waves-1) is stored starting at location
 * waves[m*num_equations].
 *
 * @param eqn Equation object
 * @param ql Conserved variables on left of interface
 * @param qr Conserved variables on right of interface
 * @param waves Waves computed from waves() method
 * @param speeds Wave speeds[num_wave]
 * @param amdq On output, the left-going fluctuations.
 * @param apdq On output, the right-going fluctuations.
 */
GKYL_CU_DH
static inline void
gkyl_wv_eqn_qfluct(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double phil, const double phir, const double *waves, const double *speeds,
  double *amdq, double *apdq)
{
  eqn->qfluct_func(eqn, type, ql, qr, phil, phir, waves, speeds, amdq, apdq);
}

/**
 * See signature for gkyl_wv_eqn_qfluct. This function computes the
 * fluctuations using f-waves rather than q-waves.
 */
GKYL_CU_DH
static inline void
gkyl_wv_eqn_ffluct(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double phil, const double phir, const double *waves, const double *speeds,
  double *amdq, double *apdq)
{
  eqn->ffluct_func(eqn, type, ql, qr, phil, phir, waves, speeds, amdq, apdq);
}

/**
 * Compute jump in flux given two conserved variable states.
 *
 * @param eqn Equation object
 * @param ql Conserved variables on left
 * @param qr Conserved variables on right
 * @param flux_jump Jump in flux (F(qr)-F(ql))
 * @return Maximum wave speed for states qr and ql.
 */
GKYL_CU_DH
static inline double
gkyl_wv_eqn_flux_jump(const struct gkyl_wv_eqn *eqn,
  const double *ql, const double *qr, double *flux_jump)
{
  return eqn->flux_jump(eqn, ql, qr, flux_jump);
}

/**
 * Check invariant domain of equation system (e.g., pressure > 0.0)
 *
 * @param eqn Equation object
 * @param q Conserved variables
 * @return boolean (true if invariant domain is satisfied, false if not)
 */
GKYL_CU_DH
static inline bool
gkyl_wv_eqn_check_inv(const struct gkyl_wv_eqn *eqn, const double *q)
{
  return eqn->check_inv_func(eqn, q);
}

/**
 * Compute waves and speeds from left/right conserved variables. The
 * 'waves' array has size num_equations X num_waves in length. The 'm'
 * wave (m = 0 ... num_waves-1) is stored starting at location
 * waves[m*num_equations].
 *
 * @param eqn Equation object
 * @param q Conserved variables
 * @return maximum wave-speed in direction 'dir'
 */
GKYL_CU_DH
static inline double
gkyl_wv_eqn_max_speed(const struct gkyl_wv_eqn *eqn, const double *q)
{
  return eqn->max_speed_func(eqn, q);
}

/**
 * Rotate state (conserved/primitive) vector to local tangent-normal coordinate frame.
 *
 * @param eqn Equation object
 * @param tau1 Tangent vector
 * @param tau2 Tangent vector
 * @param norm Normal vector such that norm = tau1 x tau2
 * @param qglobal State vector in global coordinates
 * @param qlocal State vector in local coordinates
 */
GKYL_CU_DH
static inline void
gkyl_wv_eqn_rotate_to_local(const struct gkyl_wv_eqn* eqn,
  const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qglobal, double *GKYL_RESTRICT qlocal)
{
  eqn->rotate_to_local_func(eqn, tau1, tau2, norm, qglobal, qlocal);
}

/**
 * Rotate state (conserved/primitive) vector to global coordinate frame.
 *
 * @param eqn Equation object
 * @param tau1 Tangent vector
 * @param tau2 Tangent vector
 * @param norm Normal vector such that norm = tau1 x tau2
 * @param qlocal State vector in local coordinates
 * @param qglobal State vector in local coordinates
 */
GKYL_CU_DH
static inline void
gkyl_wv_eqn_rotate_to_global(const struct gkyl_wv_eqn* eqn,
  const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qlocal, double *GKYL_RESTRICT qglobal)
{
  eqn->rotate_to_global_func(eqn, tau1, tau2, norm, qlocal, qglobal);
}

/**
* Compute forcing/source term vector from conserved variables.
*
* @param eqn Base equation object.
* @param qin Conserved variable vector (input).
* @param sout Forcing/source term vector (output).
*/
GKYL_CU_DH
static inline void
gkyl_wv_eqn_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  eqn->source_func(eqn, qin, sout);
}

/**
 * Apply the equation's state-repair hook in place. Caller is responsible
 * for first checking that eqn->repair_state_func is non-NULL — equations
 * without a meaningful repair operation leave the slot NULL.
 *
 * @param eqn    Equation object
 * @param q_prev Cell's pre-update state (admissible). NULL is allowed —
 *               equations that don't need history-informed repair ignore it.
 * @param q      Conserved variables (modified in place)
 */
GKYL_CU_DH
static inline void
gkyl_wv_eqn_repair_state(const struct gkyl_wv_eqn *eqn,
  const double *q_prev, double *q)
{
  eqn->repair_state_func(eqn, q_prev, q);
}

/**
 * Delete equation object
 *
 * @param eqn Equation object to delete.
 */
void gkyl_wv_eqn_release(const struct gkyl_wv_eqn* eqn);
