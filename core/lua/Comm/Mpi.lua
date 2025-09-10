-- Gkyl ------------------------------------------------------------------------
--
-- Lua wrapper for MPI
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- don't bother if MPI if not built in
assert(GKYL_HAVE_MPI, "Gkyl was not built with MPI!")

local ffi  = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"
local new, typeof = xsys.from(ffi,
     "new, typeof")

local Lin = require "Lib.Linalg"
local Time = require "Lib.Time"

-- global count of barriers
local numMpiBarrier = 0
-- time spent in barrier function
local timeMpiBarrier = 0.0

local _M = {}

ffi.cdef [[
  // Opaque types
  typedef struct MPI_Comm_type *MPI_Comm;
  typedef struct MPI_Datatype_type *MPI_Datatype;
  typedef struct MPI_Op_type *MPI_Op;
  typedef struct MPI_Status_type *MPI_Status;
  typedef struct MPI_Group_type *MPI_Group;
  typedef struct MPI_Request_type *MPI_Request;
  typedef struct MPI_Info_type *MPI_Info;
  typedef struct MPI_Win_type *MPI_Win;
  typedef ptrdiff_t MPI_Aint;

  // size of various objects
  int sizeof_MPI_Status();
  int sizeof_MPI_Group();
  int sizeof_MPI_Comm();
  int sizeof_MPI_Request();
  int sizeof_MPI_Aint();
  // size of various object pointers
  int sizeof_ptr_MPI_Status();
  int sizeof_ptr_MPI_Request();

  // Pre-defined objects and constants
  MPI_Comm get_MPI_COMM_WORLD();
  MPI_Comm get_MPI_COMM_NULL();
  MPI_Comm get_MPI_COMM_SELF();
  MPI_Request get_MPI_REQUEST_NULL();
  MPI_Status *getPtr_MPI_STATUS_IGNORE();
  MPI_Info get_MPI_INFO_NULL();
  int get_MPI_PROC_NULL();
  int get_MPI_COMM_TYPE_SHARED();
  int get_MPI_UNDEFINED();
  int get_MPI_ORDER_C();
  int get_MPI_ORDER_FORTRAN();
  int * get_MPI_UNWEIGHTED();
  void * get_MPI_BOTTOM();

  // Datatypes
  MPI_Datatype get_MPI_C_BOOL();
  MPI_Datatype get_MPI_CHAR();
  MPI_Datatype get_MPI_BYTE();
  MPI_Datatype get_MPI_SHORT();
  MPI_Datatype get_MPI_INT();
  MPI_Datatype get_MPI_LONG();
  MPI_Datatype get_MPI_FLOAT();
  MPI_Datatype get_MPI_DOUBLE();
  MPI_Datatype get_MPI_UNSIGNED_CHAR();
  MPI_Datatype get_MPI_UNSIGNED_SHORT();
  MPI_Datatype get_MPI_UNSIGNED();
  MPI_Datatype get_MPI_UNSIGNED_LONG();
  MPI_Datatype get_MPI_LONG_DOUBLE();
  MPI_Datatype get_MPI_LONG_LONG_INT();
  MPI_Datatype get_MPI_FLOAT_INT();
  MPI_Datatype get_MPI_LONG_INT();
  MPI_Datatype get_MPI_DOUBLE_INT();
  MPI_Datatype get_MPI_SHORT_INT();
  MPI_Datatype get_MPI_2INT();
  MPI_Datatype get_MPI_LONG_DOUBLE_INT();
  MPI_Datatype get_MPI_PACKED();

  // Operators
  MPI_Op get_MPI_MAX();
  MPI_Op get_MPI_MIN();
  MPI_Op get_MPI_SUM();
  MPI_Op get_MPI_PROD();
  MPI_Op get_MPI_LAND();
  MPI_Op get_MPI_BAND();
  MPI_Op get_MPI_LOR();
  MPI_Op get_MPI_BOR();
  MPI_Op get_MPI_LXOR();
  MPI_Op get_MPI_BXOR();
  MPI_Op get_MPI_MINLOC();
  MPI_Op get_MPI_MAXLOC();

  // Errors
  int get_MPI_SUCCESS();

  // Communicators
  int MPI_Comm_rank(MPI_Comm comm, int *rank);
  int MPI_Comm_size(MPI_Comm comm, int *size);
  int MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm);
  int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm);
  int MPI_Comm_free(MPI_Comm *comm);

  // MPI_Datatype handling
  int MPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype *newtype);
  int MPI_Type_vector(int count,
    int blocklength, int stride, MPI_Datatype oldtype, MPI_Datatype * newtype);
  int MPI_Type_indexed(int count, const int array_of_blocklengths[],
    const int array_of_displacements[], MPI_Datatype oldtype, MPI_Datatype *newtype);
  int MPI_Type_create_subarray(int ndims, const int array_of_sizes[], const
    int array_of_subsizes[], const int array_of_starts[], int order, MPI_Datatype
    oldtype, MPI_Datatype *newtype);
  int MPI_Type_commit(MPI_Datatype * datatype);
  int MPI_Type_free(MPI_Datatype *datatype);
  int MPI_Get_address(const void *location, MPI_Aint *address);

  // Win & SHM calls
  int MPI_Comm_split_type(MPI_Comm comm, int split_type, int key, MPI_Info info, MPI_Comm *newcomm);
  int MPI_Win_allocate_shared (MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm, void *baseptr, MPI_Win *win);
  int MPI_Win_shared_query(MPI_Win win, int rank, MPI_Aint *size, int *disp_unit, void *baseptr);
  int MPI_Win_lock_all(int assert, MPI_Win win);
  int MPI_Win_unlock_all(MPI_Win win);
  int MPI_Win_free(MPI_Win *win);

  // Point-to-point communication
  int MPI_Get_count(const MPI_Status *status, MPI_Datatype datatype, int *count);
  int MPI_Allgather(const void *sendbuf, int sendcount, MPI_Datatype senddatatype,
    const void *recvbuf, int recvcount, MPI_Datatype recvdatatype, MPI_Comm comm);
  int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
    MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
  int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root,
    MPI_Comm comm);
  int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
    MPI_Comm comm);
  int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
    MPI_Comm comm, MPI_Status *status);
  int MPI_Wait(MPI_Request *request, MPI_Status *status);
  int MPI_Waitall(int count, MPI_Request array_of_requests[],
    MPI_Status *array_of_statuses);
  int MPI_Waitany(int count, MPI_Request array_of_requests[],
    int *index, MPI_Status *status);

  // Nonblocking
  int MPI_Iallgather(const void *sendbuf, int sendcount, MPI_Datatype senddatatype,
    const void *recvbuf, int recvcount, MPI_Datatype recvdatatype, MPI_Comm comm);
  int MPI_Iallreduce(const void *sendbuf, void *recvbuf, int count,
    MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request);
  int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest,
    int tag, MPI_Comm comm, MPI_Request *request);
  int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source,
    int tag, MPI_Comm comm, MPI_Request *request);

  int MPI_Request_free(MPI_Request *request);

  // Groups
  int MPI_Comm_group(MPI_Comm comm, MPI_Group *group);
  int MPI_Group_rank(MPI_Group group, int *rank);
  int MPI_Group_size(MPI_Group group, int *size);
  int MPI_Group_incl(MPI_Group group, int n, const int ranks[], MPI_Group *newgroup);
  int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm);
  int MPI_Group_translate_ranks(MPI_Group group1, int n, const int ranks1[], MPI_Group group2, int ranks2[]);
  int MPI_Comm_create_group(MPI_Comm comm, MPI_Group group, int tag, MPI_Comm *newcomm);
  int MPI_Group_free(MPI_Group *group);

  // Cartesian communicator functions.
  int MPI_Cart_coords(MPI_Comm comm, int rank, int maxdims, int coords[]);
  int MPI_Cart_create(MPI_Comm comm_old, int ndims, const int dims[], const int periods[], int reorder, MPI_Comm *comm_cart);
  int MPI_Cart_get(MPI_Comm comm, int maxdims, int dims[], int periods[], int coords[]);
  int MPI_Cart_rank(MPI_Comm comm, int coords[], int *rank);
  int MPI_Cart_shift(MPI_Comm comm, int direction, int disp, int *rank_source, int *rank_dest);
  int MPI_Cart_sub(MPI_Comm comm, const int remain_dims[], MPI_Comm *comm_new);
  int MPI_Cartdim_get(MPI_Comm comm, int *ndims);

  // Graph constructors and query functions.
  int MPI_Dist_graph_create(MPI_Comm comm_old, int n, const int sources[],
                            const int degrees[], const int destinations[], const int weights[],
                            MPI_Info info, int reorder, MPI_Comm *comm_dist_graph);
  int MPI_Dist_graph_create_adjacent(MPI_Comm comm_old, int indegree, const int sources[],
                                     const int sourceweights[], int outdegree, const int destinations[],
                                     const int destweights[], MPI_Info info, int reorder,
                                     MPI_Comm *comm_dist_graph);
  int MPI_Dist_graph_neighbors(MPI_Comm comm, int maxindegree, int sources[],
                               int sourceweights[], int maxoutdegree, int destinations[],
                               int destweights[]);
  int MPI_Dist_graph_neighbors_count(MPI_Comm comm, int *indegree, int *outdegree, int *weighted);

  // Neighbor collectives.
  int MPI_Neighbor_allgather(const void *sendbuf, int  sendcount, MPI_Datatype sendtype, void *recvbuf,
                             int recvcount, MPI_Datatype recvtype, MPI_Comm comm);
  int MPI_Ineighbor_allgather(const void *sendbuf, int  sendcount, MPI_Datatype sendtype, void *recvbuf,
                              int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request req);
  int MPI_Neighbor_allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                              const int recvcounts[], const int displs[], MPI_Datatype recvtype, MPI_Comm comm);
  int MPI_Ineighbor_allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                               const int recvcounts[], const int displs[], MPI_Datatype recvtype, MPI_Comm comm,
                               MPI_Request *request);
  int MPI_Neighbor_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
                            MPI_Datatype recvtype, MPI_Comm comm);
  int MPI_Ineighbor_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
                             MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request);
  int MPI_Neighbor_alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[], MPI_Datatype sendtype,
                             void *recvbuf, const int recvcounts[], const int rdispls[], MPI_Datatype recvtype,
                             MPI_Comm comm);
  int MPI_Ineighbor_alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[], MPI_Datatype sendtype,
                              void *recvbuf, const int recvcounts[], const int rdispls[], MPI_Datatype recvtype,
                              MPI_Comm comm, MPI_Request *request);
  int MPI_Neighbor_alltoallw(const void *sendbuf, const int sendcounts[], const MPI_Aint sdispls[],
                             const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                             const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm);
  int MPI_Ineighbor_alltoallw(const void *sendbuf, const int sendcounts[], const MPI_Aint sdispls[],
                              const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                              const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Request *request);



  // Global operators
  int MPI_Barrier(MPI_Comm comm);
  int MPI_Abort(MPI_Comm comm, int errorcode);

  // CUDA test 
  int MPIX_Query_cuda_support();

  // Gkyl utility functions
  void GkMPI_fillStatus(const MPI_Status* inStatus, int *outStatus);

  // Gkyl structs holding status and requests.
  typedef struct {
     MPI_Request *req;
  } gkyl_MPI_Request;

  typedef struct {
     MPI_Status *stat;
  } gkyl_MPI_Status;

  typedef struct {
     MPI_Request *req;
     MPI_Status *stat;
  } gkyl_MPI_Request_Status;

  // Functions to allocate and free structs holding requests and statuses.
  void gkyl_MPI_Request_alloc(gkyl_MPI_Request *rs, int num);
  void gkyl_MPI_Request_release(gkyl_MPI_Request *rs);
  void gkyl_MPI_Status_alloc(gkyl_MPI_Status *ss, int num);
  void gkyl_MPI_Status_release(gkyl_MPI_Status *ss);
  void gkyl_MPI_Request_Status_alloc(gkyl_MPI_Request_Status *rss, int num);
  void gkyl_MPI_Request_Status_release(gkyl_MPI_Request_Status *rss);

  // Functions to fetch members of status.
  int gkyl_mpi_get_status_SOURCE(const MPI_Status* instat, int off);
  int gkyl_mpi_get_status_TAG(const MPI_Status* instat, int off);
  int gkyl_mpi_get_status_ERROR(const MPI_Status* instat, int off);

  // Get count from a status (which may be one of several in an array of
  // statuses).
  int gkyl_mpi_get_status_count(const MPI_Status *instat, MPI_Datatype datatype, int *count, int off);
]]

-- Predefined objects and constants
_M.COMM_WORLD    = ffiC.get_MPI_COMM_WORLD()
_M.COMM_NULL     = ffiC.get_MPI_COMM_NULL()
_M.COMM_SELF     = ffiC.get_MPI_COMM_SELF()
_M.REQUEST_NULL  = ffiC.get_MPI_REQUEST_NULL()
_M.STATUS_IGNORE = ffiC.getPtr_MPI_STATUS_IGNORE()
_M.PROC_NULL     = ffiC.get_MPI_PROC_NULL()

_M.INFO_NULL        = ffiC.get_MPI_INFO_NULL()
_M.COMM_TYPE_SHARED = ffiC.get_MPI_COMM_TYPE_SHARED()
_M.UNDEFINED        = ffiC.get_MPI_UNDEFINED()
_M.ORDER_C          = ffiC.get_MPI_ORDER_C()
_M.ORDER_FORTRAN    = ffiC.get_MPI_ORDER_FORTRAN()
_M.UNWEIGHTED       = ffiC.get_MPI_UNWEIGHTED()
_M.BOTTOM           = ffiC.get_MPI_BOTTOM()

-- Object sizes
_M.SIZEOF_AINT        = ffiC.sizeof_MPI_Aint()
_M.SIZEOF_STATUS      = ffiC.sizeof_MPI_Status()
_M.SIZEOF_REQUEST     = ffiC.sizeof_MPI_Request()
_M.SIZEOF_STATUS_PTR  = ffiC.sizeof_ptr_MPI_Status()
_M.SIZEOF_REQUEST_PTR = ffiC.sizeof_ptr_MPI_Request()

-- Datatypes
_M.C_BOOL          = ffiC.get_MPI_C_BOOL()
_M.CHAR            = ffiC.get_MPI_CHAR()
_M.BYTE            = ffiC.get_MPI_BYTE()
_M.SHORT           = ffiC.get_MPI_SHORT()
_M.INT             = ffiC.get_MPI_INT()
_M.LONG            = ffiC.get_MPI_LONG()
_M.FLOAT           = ffiC.get_MPI_FLOAT()
_M.DOUBLE          = ffiC.get_MPI_DOUBLE()
_M.UNSIGNED_CHAR   = ffiC.get_MPI_UNSIGNED_CHAR()
_M.UNSIGNED_SHORT  = ffiC.get_MPI_UNSIGNED_SHORT()
_M.UNSIGNED        = ffiC.get_MPI_UNSIGNED()
_M.UNSIGNED_LONG   = ffiC.get_MPI_UNSIGNED_LONG()
_M.LONG_DOUBLE     = ffiC.get_MPI_LONG_DOUBLE()
_M.LONG_LONG_INT   = ffiC.get_MPI_LONG_LONG_INT()
_M.FLOAT_INT       = ffiC.get_MPI_FLOAT_INT()
_M.LONG_INT        = ffiC.get_MPI_LONG_INT()
_M.DOUBLE_INT      = ffiC.get_MPI_DOUBLE_INT()
_M.SHORT_INT       = ffiC.get_MPI_SHORT_INT()
_M.TWOINT          = ffiC.get_MPI_2INT()
_M.LONG_DOUBLE_INT = ffiC.get_MPI_LONG_DOUBLE_INT()
_M.PACKED          = ffiC.get_MPI_PACKED()

-- Operators
_M.MAX    = ffiC.get_MPI_MAX()
_M.MIN    = ffiC.get_MPI_MIN()
_M.SUM    = ffiC.get_MPI_SUM()
_M.PROD   = ffiC.get_MPI_PROD()
_M.LAND   = ffiC.get_MPI_LAND()
_M.BAND   = ffiC.get_MPI_BAND()
_M.LOR    = ffiC.get_MPI_LOR()
_M.BOR    = ffiC.get_MPI_BOR()
_M.LXOR   = ffiC.get_MPI_LXOR()
_M.BXOR   = ffiC.get_MPI_BXOR()
_M.MINLOC = ffiC.get_MPI_MINLOC()
_M.MAXLOC = ffiC.get_MPI_MAXLOC()

-- Error codes
_M.SUCCESS = ffiC.get_MPI_SUCCESS()

-- some types for use in MPI functions
local int_1  = typeof("int[1]")
local uint_1 = typeof("unsigned[1]")
local voidp  = typeof("void *[1]")

-- ctors for various MPI objects.
-- MF 2022/10/04: we keep requests and statuses in our own structs
-- because we encountered some seg faults when using non-blocking
-- send and receives and couldn't find another way to do it safely.
-- Probably something to do with Lua's GC not handling these objects
-- robustly, in general.
local function new_MPI_Request(sz)
   local sz = sz or 1
   local rs = ffi.gc(new("gkyl_MPI_Request"), ffiC.gkyl_MPI_Request_release)
   ffiC.gkyl_MPI_Request_alloc(rs, sz)
   return rs
end
local function new_MPI_Status(sz)
   local sz = sz or 1
   local ss = ffi.gc(new("gkyl_MPI_Status"), ffiC.gkyl_MPI_Status_release)
   ffiC.gkyl_MPI_Status_alloc(ss, sz)
   return ss
end
local function new_MPI_RequestStatus(sz)
   local sz = sz or 1
   local rss = ffi.gc(new("gkyl_MPI_Request_Status"), ffiC.gkyl_MPI_Request_Status_release) 
   ffiC.gkyl_MPI_Request_Status_alloc(rss, sz)
   return rss
end
local function new_MPI_Comm()
   return new("MPI_Comm[1]")
end
local function new_MPI_Group()
   return new("MPI_Group[1]")
end
local function new_MPI_Win()
   return new("MPI_Win[1]")
end
local function new_MPI_Datatype()
   return new("MPI_Datatype[1]")
end
local function new_MPI_Datatype_vec(sz)
   local sz = sz or 1
   return new("MPI_Datatype[?]", sz)
end
local function new_MPI_Aint()
   return new("MPI_Aint[1]")
end
local function new_MPI_Aint_vec(sz)
   local sz = sz or 1
   return new("MPI_Aint[?]", sz)
end

-- de-reference if object is a pointer
local function getObj(obj, ptyp)
   return ffi.istype(typeof(ptyp), obj) and obj[0] or obj
end
local function getRequest(obj)
   return (ffi.istype(typeof("gkyl_MPI_Request"), obj)
           or ffi.istype(typeof("gkyl_MPI_Request_Status"), obj))
          and obj.req or obj
end
local function getStatus(obj)
   return (ffi.istype(typeof("gkyl_MPI_Status"), obj)
           or ffi.istype(typeof("gkyl_MPI_Request_Status"), obj))
          and obj.stat or obj
end

-- Extract object from (potentially) a pointer
function _M.getComm(obj)
   return getObj(obj, "MPI_Comm[1]")
end

-- MPI_Status object
function _M.Status(sz) return new_MPI_Status(sz) end 
-- MPI_Request types.
function _M.Request(sz) return new_MPI_Request(sz) end 
-- MPI_Request/MPI_Status pair.
function _M.RequestStatus(sz) return new_MPI_RequestStatus(sz) end 

-- Functions to obtain members of MPI_Status.
function _M.getStatusSOURCE(stat, off)
   local off = off or 0
   return ffiC.gkyl_mpi_get_status_SOURCE(getStatus(stat), off);
end
function _M.getStatusTAG(stat, off)
   local off = off or 0
   return ffiC.gkyl_mpi_get_status_TAG(getStatus(stat), off);
end
function _M.getStatusERROR(stat, off)
   local off = off or 0
   return ffiC.gkyl_mpi_get_status_ERROR(getStatus(stat), off);
end

-- MPI_Comm object
function _M.Comm()  return new_MPI_Comm() end
-- MPI_Group object
function _M.Group()  return new_MPI_Group() end
-- MPI_Datatype object
function _M.MPI_Datatype() return new_MPI_Datatype() end
function _M.MPI_Datatype_vec(sz) return new_MPI_Datatype_vec(sz) end
-- MPI_Aint object
function _M.MPI_Aint() return new_MPI_Aint() end
function _M.MPI_Aint_vec(sz) return new_MPI_Aint_vec(sz) end

-- MPI_Comm_rank
function _M.Comm_rank(comm)
   local r = int_1()
   local _ = ffiC.MPI_Comm_rank(getObj(comm, "MPI_Comm[1]"), r)
   return r[0]
end
-- MPI_Comm_size
function _M.Comm_size(comm)
   local r = int_1()
   local _ = ffiC.MPI_Comm_size(getObj(comm, "MPI_Comm[1]"), r)
   return r[0]
end
-- MPI_Comm_dup
function _M.Comm_dup(comm)
   local c = new_MPI_Comm()
   local _ = ffiC.MPI_Comm_dup(getObj(comm, "MPI_Comm[1]"), c)
   return c
end
-- MPI_Comm_split
function _M.Comm_split(comm, color, key)
   local c = new_MPI_Comm()
   local _ = ffiC.MPI_Comm_split(getObj(comm, "MPI_Comm[1]"), color, key, c)
   return c
end
-- MPI_Comm_free
function _M.Comm_free(comm)
   local _ = ffiC.MPI_Comm_free(comm)
end

-- MPI_Comm_split_type
function _M.Comm_split_type(comm, split_type, key, info)
   local c = new_MPI_Comm()
   local _ = ffiC.MPI_Comm_split_type(getObj(comm, "MPI_Comm[1]"), split_type, key, info, c)
   return c
end
-- MPI_Win_allocate_shared
function _M.Win_allocate_shared(sz, disp_unit, info, comm)
   local w = new_MPI_Win()
   local baseptr = voidp()
   local _ = ffiC.MPI_Win_allocate_shared(sz, disp_unit, info, getObj(comm, "MPI_Comm[1]"), baseptr, w)
   return baseptr[0], w
end
-- MPI_Win_shared_query
function _M.Win_shared_query(win, rank)
   local sz, du = uint_1(), int_1()
   local baseptr = voidp()
   local _ = ffiC.MPI_Win_shared_query(getObj(win, "MPI_Win[1]"), rank, sz, du, baseptr)
   return sz[0], du[0], baseptr[0]
end
-- MPI_Win_lock_all
function _M.Win_lock_all(assert, win)
   local _ = ffiC.MPI_Win_lock_all(assert, getObj(win, "MPI_Win[1]"))
end
-- MPI_Win_unlock_all
function _M.Win_unlock_all(win)
   local _ = ffiC.MPI_Win_unlock_all(getObj(win, "MPI_Win[1]"))
end
-- MPI_Win_free
function _M.Win_free(win)
   local _ = ffiC.MPI_Win_free(win)
end

-- MPI_Type_contiguous
function _M.Type_contiguous(count, oldtype)
   local t = new_MPI_Datatype()
   local _ = ffiC.MPI_Type_contiguous(count, getObj(oldtype, "MPI_Datatype[1]"), t)
   return t
end
-- MPI_Type_vector
function _M.Type_vector(count, blocklength, stride, oldtype)
   local t = new_MPI_Datatype()
   local _ = ffiC.MPI_Type_vector(
      count, blocklength, stride,
      getObj(oldtype, "MPI_Datatype[1]"), t)
   return t
end
-- MPI_Type_indexed
function _M.Type_indexed(array_of_blocklengths, array_of_displacements, oldtype)
   local t = new_MPI_Datatype()
   local count = #array_of_blocklengths
   local _ = ffiC.MPI_Type_indexed(
      count, array_of_blocklengths:data(), array_of_displacements:data(),
      getObj(oldtype, "MPI_Datatype[1]"), t)
   return t
end
-- MPI_Type_create_subarray
function _M.Type_create_subarray(array_of_sizes, array_of_subsizes, array_of_starts, order, oldtype)
   local t = new_MPI_Datatype()
   local ndims = #array_of_sizes
   local _ = ffiC.MPI_Type_create_subarray(
      ndims, array_of_sizes,
      array_of_subsizes, array_of_starts, order,
      getObj(oldtype, "MPI_Datatype[1]"), t)
   return t
end
-- MPI_Type_commit
function _M.Type_commit(datatype)
   local _ = ffiC.MPI_Type_commit(datatype)
   return datatype
end
-- MPI_Type_free
function _M.Type_free(datatype)
   local _ = ffiC.MPI_Type_free(datatype)
end
-- int MPI_Get_address(const void *location, MPI_Aint *address)
function _M.Get_address(location)
   local newAddress = new_MPI_Aint()
   local _ = ffiC.MPI_Get_address(location, newAddress)
   return newAddress
end

-- MPI_Get_count
function _M.Get_count(status, datatype, off)
   local off = off or 0
   local r = int_1()
   local _ = ffiC.gkyl_mpi_get_status_count(getStatus(status), datatype, r, off)
   return r[0]
end
-- MPI_Allgather
function _M.Allgather(sendbuff, sendcount, senddatatype, recvbuf, recvcount, recvdatatype, comm)
   local _ = ffiC.MPI_Allgather(sendbuff, sendcount, senddatatype,
                                recvbuf, recvcount, recvdatatype, getObj(comm, "MPI_Comm[1]"))
end
-- MPI_Allreduce
function _M.Allreduce(sendbuf, recvbuf, count, datatype, op, comm)
   local _ = ffiC.MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, getObj(comm, "MPI_Comm[1]"))
end
-- MPI_Bcast
function _M.Bcast(buffer, count, datatype, root, comm)
   local _ = ffiC.MPI_Bcast(buffer, count, datatype, root, getObj(comm, "MPI_Comm[1]"))
end
-- MPI_Send
function _M.Send(buf, count, datatype, dest, tag, comm)
   local _ = ffiC.MPI_Send(
      buf, count, getObj(datatype, "MPI_Datatype[1]"), dest, tag, getObj(comm, "MPI_Comm[1]"))
end
-- MPI_Recv
function _M.Recv(buf, count, datatype, source, tag, comm, status)
   local st = status or _M.STATUS_IGNORE
   local _ = ffiC.MPI_Recv(
      buf, count, getObj(datatype, "MPI_Datatype[1]"), source, tag, getObj(comm, "MPI_Comm[1]"), getStatus(st))
end
-- MPI_Wait
function _M.Wait(request, status)
   local st = status or _M.STATUS_IGNORE
   local err = ffiC.MPI_Wait(getRequest(request), getStatus(st))
   return err
end
-- MPI_Waitall
function _M.Waitall(count, request, status)
   local st = status or _M.STATUS_IGNORE
   local err = ffiC.MPI_Waitall(count, getRequest(request), getStatus(st))
   return err
end
---- MPI_Waitany
----  int MPI_Waitany(int count, MPI_Request array_of_requests[],
----                  int *index, MPI_Status *status);
--function _M.Waitany(count, requestArray, ind, status)
--   local st = status and status.mpiStatus or _M.STATUS_IGNORE
--   local _ = ffiC.MPI_Waitany(count, request, st)
--   -- store MPI_Status
--   if statusArray ~= nil then
--      local gks = new("int[3]")
--      ffiC.GkMPI_fillStatus(st, gks)
--      status.SOURCE, status.TAG, status.ERROR = gks[0], gks[1], gks[2]
--   end
--end

-- MPI_Iallgather
function _M.Iallgather(sendbuff, sendcount, senddatatype, recvbuf, recvcount, recvdatatype, comm)
   return ffiC.MPI_Iallgather(sendbuff, sendcount, senddatatype, recvbuf, recvcount, recvdatatype, getObj(comm, "MPI_Comm[1]")) --not sure which comm to getObj of
end
-- MPI_Iallreduce
function _M.Iallreduce(sendbuf, recvbuf, count, datatype, op, tag, comm, req)
   return ffiC.MPI_Iallreduce(sendbuf, recvbuf, count, datatype, op, getObj(comm, "MPI_Comm[1]"), getRequest(req))
end
-- MPI_Isend
function _M.Isend(buf, count, datatype, dest, tag, comm, req)
   return ffiC.MPI_Isend(
      buf, count, getObj(datatype, "MPI_Datatype[1]"), dest, tag, getObj(comm, "MPI_Comm[1]"), getRequest(req))
end
-- MPI_Irecv
function _M.Irecv(buf, count, datatype, source, tag, comm, req)
   return ffiC.MPI_Irecv(
      buf, count, getObj(datatype, "MPI_Datatype[1]"), source, tag, getObj(comm, "MPI_Comm[1]"), getRequest(req))
end

-- MPI_Request_free
function _M.Request_free(req)
  return ffiC.MPI_Request_free(getRequest(req))
end

-- MPI_Barrier
function _M.Barrier(comm)
   numMpiBarrier = numMpiBarrier+1
   local tm = Time.clock()
   local _ = ffiC.MPI_Barrier(getObj(comm, "MPI_Comm[1]"))
   timeMpiBarrier = timeMpiBarrier + (Time.clock() - tm)
end
-- MPI_Abort
function _M.Abort(comm, errCode)
   local _ = ffiC.MPI_Abort(getObj(comm, "MPI_Comm[1]"), errCode)
end

-- MPI_Comm_group
function _M.Comm_group(comm)
   local grp = new_MPI_Group()
   local _ = ffiC.MPI_Comm_group(getObj(comm, "MPI_Comm[1]"), grp)
   return grp
end
-- MPI_Group_rank
function _M.Group_rank(group)
   local r = int_1()
   local _ = ffiC.MPI_Group_rank(getObj(group, "MPI_Group[1]"), r)
   return r[0]
end
-- MPI_Group_size
function _M.Group_size(group)
   local r = int_1()
   local _ = ffiC.MPI_Group_size(getObj(group, "MPI_Group[1]"), r)
   return r[0]
end
-- MPI_Group_incl
function _M.Group_incl(group, n, ranks)
   local newgroup = new_MPI_Group()
   local _ = ffiC.MPI_Group_incl(getObj(group, "MPI_Group[1]"), n, ranks, newgroup)
   return newgroup
end
-- MPI_Comm_create
function _M.Comm_create(comm, group)
   local c = new_MPI_Comm()
   local _ = ffiC.MPI_Comm_create(
      getObj(comm, "MPI_Comm[1]"), getObj(group, "MPI_Group[1]"), c)
   return c
end
-- MPI_Group_translate_ranks
function _M.Group_translate_ranks(group1, ranks1, group2)
   local n = #ranks1
   local ranks2 = Lin.IntVec(n)
   local _ = ffiC.MPI_Group_translate_ranks(
      getObj(group1, "MPI_Group[1]"), n, ranks1:data(), getObj(group2, "MPI_Group[1]"), ranks2:data())
   return ranks2
end
-- MPI_Comm_create_group.
function _M.Comm_create_group(comm, group, tag)
   local c = new_MPI_Comm()
   local _ = ffiC.MPI_Comm_create_group(getObj(comm, "MPI_Comm[1]"), getObj(group, "MPI_Group[1]"), tag, c)
   return c
end
-- MPI_Group_free
function _M.Group_free(group)
   local _ = ffiC.MPI_Group_free(group)
end

-- MPI_Cart_coords
function _M.Cart_coords(comm, rank, maxdims)
   local coords = Lin.IntVec(maxdims)
   local _ = ffiC.MPI_Cart_coords(getObj(comm, "MPI_Comm[1]"), rank, maxdims, coords:data()) 
   return coords
end
-- MPI_Cart_create
function _M.Cart_create(commOld, ndims, dims, isDirPeriodic, reorder)
   local commCart = new_MPI_Comm()
   local _ = ffiC.MPI_Cart_create(getObj(commOld, "MPI_Comm[1]"), ndims, dims:data(), isDirPeriodic:data(), reorder, commCart)
   return commCart
end
-- MPI_Cart_get
function _M.Cart_get(comm, maxdims)
   local dims = Lin.IntVec(maxdims)
   local isDirPeriodic = Lin.IntVec(maxdims)
   local coords = Lin.IntVec(maxdims)
   local _ = ffiC.MPI_Cart_get(getObj(comm, "MPI_Comm[1]"), maxdims, dims:data(), isDirPeriodic:data(), coords:data())
   return dims, isDirPeriodic, coords
end
-- MPI_Cart_rank
function _M.Cart_rank(comm, coords)
   local r = int_1()
   local _ = ffiC.MPI_Cart_rank(getObj(comm, "MPI_Comm[1]"), coords:data(), r)
   return r[0]
end
-- MPI_Cart_shift
function _M.Cart_shift(comm, direction, disp)
   local rSrc, rDest = int_1(), int_1()
   local _ = ffiC.MPI_Cart_shift(getObj(comm, "MPI_Comm[1]"), direction, disp, rSrc, rDest)
   return rSrc[0], rDest[0]
end
-- MPI_Cart_sub
function _M.Cart_sub(comm, remainDims)
   local subComm = new_MPI_Comm()
   local err = ffiC.MPI_Cart_sub(getObj(comm, "MPI_Comm[1]"), remainDims:data(), subComm)
   return subComm
end
-- MPI_Cartdim_get
function _M.Cartdim_get(comm)
   local ndims = int_1()
   local _ = ffiC.MPI_Cartdim_get(getObj(comm, "MPI_Comm[1]"), ndims)
   return ndims[0]
end

-- MPI_Dist_graph_create_adjacent.
function _M.Dist_graph_create_adjacent(comm_old, indegree, sources, sourceweights, outdegree, destinations,
                                       destweights, info, reorder)
   local comm_dist_graph = new_MPI_Comm()
   local _ = ffiC.MPI_Dist_graph_create_adjacent(getObj(comm_old, "MPI_Comm[1]"), indegree, sources:data(),
      sourceweights, outdegree, destinations:data(), destweights, info, reorder,
      comm_dist_graph)
   return comm_dist_graph
end
-- MPI_Dist_graph_neighbors.
function _M.Dist_graph_neighbors(comm, maxindegree, maxoutdegree)
   local sources, sourceweights = Lin.IntVec(maxindegree), Lin.IntVec(maxindegree)
   local destinations, destinationweights = Lin.IntVec(maxoutdegree), Lin.IntVec(maxoutdegree)
   local _ = ffiC.MPI_Dist_graph_neighbors(getObj(comm, "MPI_Comm[1]"),
     maxindegree, sources:data(), sourceweights:data(),
     maxoutdegree, destinations:data(), destinationweights:data());
   return sources, sourceweights, destinations, destinationweights
end
-- MPI_Dist_graph_neighbors_count.
function _M.Dist_graph_neighbors_count(comm)
  local indegree, outdegree, weighted = int_1(), int_1(), int_1()
  local _ = ffiC.MPI_Dist_graph_neighbors_count(getObj(comm, "MPI_Comm[1]"), indegree, outdegree, weighted);
  return indegree[0], outdegree[0], weighted[0]==1 and true or false
end
-- MPI_Neighbor_allgather.
function _M.Neighbor_allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm)
  local _ = ffiC.MPI_Neighbor_allgather(sendbuf, sendcount, getObj(sendtype, "MPI_Datatype[1]"), recvbuf,
                                        recvcount, getObj(recvtype, "MPI_Datatype[1]"), getObj(comm, "MPI_Comm[1]"));
end
-- MPI_Neighbor_allgatherv.
function _M.Neighbor_allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm)
  local _ = ffiC.MPI_Neighbor_allgatherv(sendbuf, sendcount, getObj(sendtype, "MPI_Datatype[1]"), recvbuf,
                                         recvcounts, displs, getObj(recvtype, "MPI_Datatype[1]"), getObj(comm, "MPI_Comm[1]"
));
end
-- MPI_Neighbor_alltoall.
function _M.Neighbor_alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm)
  local _ = ffiC.MPI_Neighbor_alltoall(sendbuf, sendcount, getObj(sendtype, "MPI_Datatype[1]"),
                                       recvbuf, recvcount, getObj(recvtype, "MPI_Datatype[1]"),
                                       getObj(comm, "MPI_Comm[1]"));
end
-- MPI_Neighbor_alltoallv.
function _M.Neighbor_alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, recvtype, comm)
  local _ = ffiC.MPI_Neighbor_alltoallv(sendbuf, sendcounts, sdispls, getObj(sendtype, "MPI_Datatype[1]"),
                                        recvbuf, recvcounts, rdispls, getObj(recvtype, "MPI_Datatype[1]"),
                                        getObj(comm, "MPI_Comm[1]"));
end
-- MPI_Neighbor_alltoallw.
function _M.Neighbor_alltoallw(sendbuf, sendcounts, sdispls, sendtypes, recvbuf, recvcounts, rdispls, recvtypes, comm)
  -- Here we assume sendtypes and recvtypes is an array of MPI_Datatype, so we don't need to getObj
  -- as in previous functions, otherwise it would pick the first type only and not the whole array.
  local _ = ffiC.MPI_Neighbor_alltoallw(sendbuf, sendcounts, sdispls, sendtypes,
                                        recvbuf, recvcounts, rdispls, recvtypes,
                                        getObj(comm, "MPI_Comm[1]"));
end


-- Convenience functions (these are not wrappers over MPI but make
-- some things a little cleaner)

-- Collect 'ranks' from 'comm' and create a new communicator with just
-- those ranks. The 'ranks' parameter must be of a Linalg.IntVec object
function _M.Split_comm(comm, ranks)
   local grp = _M.Comm_group(comm)
   local newGrp = _M.Group_incl(grp, #ranks, ranks:data())
   return _M.Comm_create(comm, newGrp)
end

-- Check if the communicator is NULL
function _M.Is_comm_valid(comm)
   return getObj(comm, "MPI_Comm[1]") ~=  _M.COMM_NULL
end

-- Makes a MPI_Datatype object from block sizes and offsets
function _M.createDataTypeFromBlockSizeAndOffset(blockSize, blockOffset, oldtype)
   if #blockSize == 1 then
      return _M.Type_contiguous(blockSize[1], oldtype)
   end

   -- check if sizes of all blocks and relative offsets are the same
   local isVectorType = true
   local sz = blockSize[1]
   for i = 2, #blockSize do
      if blockSize[i] ~= sz then
	 isVectorType = false
      end
   end
   if #blockSize > 1 then
      local stride = blockOffset[2]-blockOffset[1]
      for i = 2, #blockSize do
	 if stride ~= blockOffset[i]-blockOffset[i-1] then
	    isVectorType = false
	 end
      end
   end

   if isVectorType then
      return _M.Type_vector(
	 #blockOffset,
	 blockSize[1],
	 #blockOffset>1 and blockOffset[2]-blockOffset[1] or 0,
	 oldtype)
   else
      -- must be contructed as an indexed type
      local array_of_blocklengths, array_of_displacements
	 = Lin.IntVec(#blockSize), Lin.IntVec(#blockSize)
      array_of_blocklengths:setValues(blockSize)
      array_of_displacements:setValues(blockOffset)
      return _M.Type_indexed(array_of_blocklengths, array_of_displacements, oldtype)
   end
end

-- Constructs block offsets and sizes from range and a specified
-- sub-range. Sub-range must be completely inside the range
function _M.createBlockInfoFromRangeAndSubRange(subRange, range, numComponent, ordering)
   local indexer = range:genIndexer(ordering)
   local lo = subRange:lowerAsVec()
   local up = subRange:upperAsVec()

   -- check if subRange is contiguous subset of range
   local dist = indexer(up)-indexer(lo)+1
   if dist == subRange:volume() then
      return { subRange:volume()*numComponent }, { 0 }
   end

   local firstOffset = indexer(lo)
   local blockSize, blockOffset = {}, {}   
   local currBlockSize, currOffset = 1, indexer(subRange:lowerAsVec())
   local count = 0
   local lastLinIdx = indexer(lo)
   
   for idx in subRange:iter(ordering) do
      if count == 0 then
	 count = 1
      else
	 local linIdx = indexer(idx)
	 if linIdx-lastLinIdx == 1 then -- indices are contiguous
	    currBlockSize = currBlockSize+1
	 else
	    -- store size and index
	    blockSize[count] = currBlockSize*numComponent
	    blockOffset[count] = (currOffset-firstOffset)*numComponent
	    
	    -- prep for next round
	    currBlockSize = 1
	    currOffset = linIdx
	    count = count+1
	 end
	 lastLinIdx = linIdx -- for next round
	 blockSize[count] = currBlockSize*numComponent
	 blockOffset[count] = (currOffset-firstOffset)*numComponent
      end
   end
   return blockSize, blockOffset
end

-- Creates an MPI_Datatype object from a range object in a specified
-- direction and ordering. 'ordering' must be one of Range.rowMajor or
-- Range.colMajor.
--
-- 'nlayer' is number of layers in range to include in datatype and
-- 'numComponent' is the number of components in field (usually 1)
function _M.createDataTypeFromRange(dir, range, nlayer, numComponent, ordering, oldtype)
   local blockSize, blockOffset = _M.createBlockInfoFromRangeAndSubRange(
	 range:shorten(dir, nlayer), range, numComponent, ordering)
   return _M.Type_commit(
      _M.createDataTypeFromBlockSizeAndOffset(blockSize, blockOffset, oldtype)
   )
end

-- Creates an MPI_Datatype object from a range and its subRange object.
-- 'ordering' must be one of Range.rowMajor or Range.colMajor.
--
-- 'numComponent' is the number of components in field (usually 1)
function _M.createDataTypeFromRangeAndSubRange(subRange, range, numComponent, ordering, oldtype)
   local blockSize, blockOffset = _M.createBlockInfoFromRangeAndSubRange(
      subRange, range, numComponent, ordering)
   return _M.Type_commit(
      _M.createDataTypeFromBlockSizeAndOffset(blockSize, blockOffset, oldtype)
   )
end

-- Gets total number of barriers called in system
function _M.getNumBarriers()  return numMpiBarrier end

-- Gets total time spent in barriers
function _M.getTimeBarriers()  return timeMpiBarrier end

-- MPI_Abort
function _M.Query_cuda_support()
   if ffiC.MPIX_Query_cuda_support then
      return ffiC.MPIX_Query_cuda_support() == 1 and true or false
   end
   return false
end

return _M
