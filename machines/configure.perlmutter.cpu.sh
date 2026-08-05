module unload darshan
module load cray-mpich/9.0.1
module load cray-libsci/25.09.0

: "${PREFIX:=$HOME/gkylsoft}"

./configure CC=cc --prefix=$PREFIX --lapack-lib-name=sci_gnu --lapack-inc=$CRAY_LIBSCI_PREFIX/include --lapack-lib=$CRAY_LIBSCI_PREFIX/lib --use-mpi=yes --mpi-inc=$CRAY_MPICH_DIR/include --mpi-lib=$CRAY_MPICH_DIR/lib --use-lua=yes;


