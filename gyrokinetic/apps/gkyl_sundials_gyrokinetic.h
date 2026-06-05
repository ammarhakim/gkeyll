#pragma once

/**
 *
 * Gyrokinetic-specific methods for SUNDIALS time integrators.
 *
 */

#include <gkyl_sundials.h>

/**
 * Assign function pointers for methods needed in SUNDIALS' time integrators
 * that are specific to the gyrokinetic solver (e.g. the function that
 * calculates df/dt).
 *
 * @param gksun Gkeyll's Sundials object.
 */
void gkyl_sundials_gyrokinetic_assign_methods(struct gkyl_sundials *gksun);
