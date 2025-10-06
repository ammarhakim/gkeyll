: "${PREFIX:=$HOME/gkylsoft}"
: "${MPI_HOME:=$HOME/gkylsoft/openmpi}"
./configure CC=gcc --prefix=$PREFIX --use-lua=yes --use-mpi=yes --mpi-inc=$MPI_HOME/include --app=vlasov