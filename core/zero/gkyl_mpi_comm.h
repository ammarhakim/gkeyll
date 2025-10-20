#pragma once

#ifdef GKYL_HAVE_MPI

#include <mpi.h>
#include <gkyl_comm.h>
#include <gkyl_rect_decomp.h>

// input to create new MPI communicator
struct gkyl_mpi_comm_inp {
  MPI_Comm mpi_comm; // MPI communicator to use
  const struct gkyl_rect_decomp *decomp; // pre-computed decomposition
  bool sync_corners; // should we sync corners?
};

/**
 * Return a new MPI communicator. The methods in the communicator work
 * on arrays defined on the decomp specified in the input. This
 * contract is not checked (and can't be checked) and so the user of
 * this object must ensure consistency of the arrays and the
 * decomposition.
 *
 * @param inp Input struct to use for initialization
 * @return New MPI communicator
 */
struct gkyl_comm *gkyl_mpi_comm_new(const struct gkyl_mpi_comm_inp *inp);

#endif
