module unload darshan
module unload cudatoolkit
module unload gpu
module load cray-mpich/9.0.1
module load cray-libsci/25.09.0
module load craype-accel-host

: "${PREFIX:=/global/homes/m/mana/perlmutter/gkeyll/code/gkeyll_cpu_v2/gkylsoft}"

cd install-deps
./mkdeps.sh CC=cc MPICC=$CRAY_MPICH_DIR/bin/mpicc MPICXX=$CRAY_MPICH_DIR/bin/mpicxx MPIFC=$CRAY_MPICH_DIR/bin/mpifort MPIRUN=srun --build-superlu=yes --build-luajit=yes --prefix=$PREFIX --build-sundials=yes --lapack-lib-name=sci_gnu_mpi --lapack-inc=$CRAY_LIBSCI_PREFIX/include --lapack-lib=$CRAY_LIBSCI_PREFIX/lib
