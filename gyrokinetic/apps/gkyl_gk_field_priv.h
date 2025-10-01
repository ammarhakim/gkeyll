
// Function pointer for the rhs of the field solve
typedef void (*gk_field_fem_init_func_t)(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Project a DG field onto the parallel FEM basis to make it
 * continuous along z (or to solve a Poisson equation in 1x).
 *
 * @param app Application object
 * @param field Field object
 * @param arr_dg Input DG array (to be projected)
 * @param arr_fem Output FEM array (projected result)
 */
void
gk_field_fem_projection_par(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *arr_dg, struct gkyl_array *arr_fem);

/**
 * Initialize the FEM object for the Boltzmann field.
 *
 * @param app Application object
 * @param f Field object
 */
void
gk_field_fem_init_boltzmann(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Initialize the FEM object for the 1x field.
 *
 * @param app Application object
 * @param f Field object
 */
void
gk_field_fem_init_1x(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Initialize the FEM object for the 2x or 3x field.
 *
 * @param app Application object
 * @param f Field object
 */
void
gk_field_fem_init_2x3x(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Accumulate charge density for a Boltzmann response.
 *
 * @param app Application object
 * @param field Field object
 * @param s Species object
 */
void
gk_field_accumulate_rho_c_boltzmann(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gk_species *s);

/**
 * Accumulate charge density for an adiabatic response.
 *
 * @param app Application object
 * @param field Field object
 * @param s Species object
 */
void
gk_field_accumulate_rho_c_adiabatic(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gk_species *s);

/**
 * Accumulate charge density for a Poisson solve.
 *
 * @param app Application object
 * @param field Field object
 * @param s Species object
 */
void
gk_field_accumulate_rho_c_poisson(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gk_species *s);

/**
 * Solve for the field using a Boltzmann response.
 *
 * @param app Application object
 * @param field Field object
 */
void
gk_field_boltzmann_solve(struct gkyl_gyrokinetic_app *app, struct gk_field *field);

/**
 * Solve the Poisson equation in 1x with the parallel FEM projection.
 *
 * @param app Application object
 * @param field Field object
 */
void
gk_field_poisson_solve_1x(struct gkyl_gyrokinetic_app *app, struct gk_field *field);

/**
 * Solve the Poisson equation for IWL sims.
 *
 * @param app Application object
 * @param field Field object
 */
void
gk_field_poisson_deflate_solve_es_iwl(struct gkyl_gyrokinetic_app *app, struct gk_field *field);

/**
 * Solve the Poisson equation in 2x and 3x with the parallel FEM projection.
 *
 * @param app Application object
 * @param field Field object
 */
void
gk_field_poisson_perp_solve_2x3x(struct gkyl_gyrokinetic_app *app, struct gk_field *field);