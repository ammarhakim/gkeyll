
// Function pointer types for FEM object initialization
typedef void (*gk_field_fem_init_func_t)(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

void
gk_field_fem_projection_par(gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *arr_dg, struct gkyl_array *arr_fem);

void
gk_field_fem_init_boltzmann(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

void
gk_field_fem_init_1x(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

void
gk_field_fem_init_2x3x(struct gkyl_gyrokinetic_app *app, struct gk_field *f);


void
gk_field_boltzmann_solve(struct gkyl_gyrokinetic_app *app, struct gk_field *field);

void
gk_field_poisson_solve_1x(struct gkyl_gyrokinetic_app *app, struct gk_field *field);

void
gk_field_poisson_deflate_solve_es_iwl(struct gkyl_gyrokinetic_app *app, struct gk_field *field);

void
gk_field_poisson_perp_solve_2x3x(struct gkyl_gyrokinetic_app *app, struct gk_field *field);