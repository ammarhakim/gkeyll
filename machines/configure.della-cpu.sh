module load intel/2021.1.2
module load openmpi/intel-2021.1/4.1.0 

: "${PREFIX:=$HOME/gkylsoft}"

./configure CC=cc --prefix=$PREFIX --lapack-inc=$PREFIX/OpenBLAS/include --lapack-lib=$PREFIX/OpenBLAS/lib/libopenblas.a --superlu-inc=$PREFIX/superlu/include --superlu-lib=$PREFIX/superlu/lib/libsuperlu.a --use-mpi=yes --mpi-inc=$MPI_HOME/include --mpi-lib=$MPI_HOME/lib64 --use-lua=yes;
