module load intel/2021.1.2 
module load openmpi/intel-2021.1/4.1.2 

: "${PREFIX:=$HOME/gkylsoft}"
./configure CC=cc --prefix=$PREFIX --use-mpi=yes --mpi-inc=$MPI_HOME/include --mpi-lib=$MPI_HOME/lib64 --use-lua=yes

