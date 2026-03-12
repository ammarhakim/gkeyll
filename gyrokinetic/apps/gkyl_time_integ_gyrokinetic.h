#pragma once
// Some definitions for objects used in time integration. 

#include <gkyl_array.h>

// Arguments to df/dt function.
struct gkyl_gyrokinetic_fdot_args {
  int num_species; // Number of species.
  int num_neut_species; // Number of neutral species.
  const struct gkyl_array **fin; // Input charged species state (e.g. distribution function).
  const struct gkyl_array **fin_neut; // Input neutral species state (e.g. distribution function).
  const struct gkyl_array **fieldin; // Input electromagnetic fields.
  struct gkyl_array **fout; // Time rate of change of fin (df/dt), or new charged f.
  struct gkyl_array **fout_neut; // Time rate of change of fin_neut (df/dt), or new neutral f.
  struct gkyl_array **fieldout; // Output electromagnetic fields.
  struct gkyl_array ***bflux_in; // Input moments of boundary fluxes of charged species.
  struct gkyl_array ***bflux_in_neut; // Input moments of boundary fluxes of neutral species.
  struct gkyl_array ***bflux_out; // Output moments of boundary fluxes of charged species.
  struct gkyl_array ***bflux_out_neut; // Output moments of boundary fluxes of neutral species.
  // Number of state vectors and offset in a vector of all state vectors (for SUNDIALS).
  int *num_arr_distf_charged, *offset_distf_charged;
  int *num_arr_bflux_charged, *offset_bflux_charged;
  int *num_arr_distf_neut, *offset_distf_neut;
  int *num_arr_bflux_neut, *offset_bflux_neut;
  int *num_arr_field, *offset_field;
};

// Arguments to df/dt function.
struct gkyl_gyrokinetic_multib_fdot_args {
  int nblocks_local; // Number of blocks in this MPI process.
  int num_species; // Number of species.
  int num_neut_species; // Number of neutral species.
  const struct gkyl_array **fin; // Input charged species state (e.g. distribution function).
  const struct gkyl_array **fin_neut; // Input neutral species state (e.g. distribution function).
  const struct gkyl_array **fieldin; // Input electromagnetic fields.
  struct gkyl_array **fout; // Time rate of change of fin (df/dt), or new charged f.
  struct gkyl_array **fout_neut; // Time rate of change of fin_neut (df/dt), or new neutral f.
  struct gkyl_array **fieldout; // Output electromagnetic fields.
  struct gkyl_array ***bflux_in; // Input moments of boundary fluxes of charged species.
  struct gkyl_array ***bflux_in_neut; // Input moments of boundary fluxes of neutral species.
  struct gkyl_array ***bflux_out; // Output moments of boundary fluxes of charged species.
  struct gkyl_array ***bflux_out_neut; // Output moments of boundary fluxes of neutral species.
};

