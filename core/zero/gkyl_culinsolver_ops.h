#ifdef GKYL_HAVE_CUDSS

// Code was compiled with cuDSS.
#include <gkyl_cudss_ops.h>

#else

// Code was compiled without cuDSS, use cuSolver/cuSparse instead.
#include <gkyl_cusolver_ops.h>

#endif
