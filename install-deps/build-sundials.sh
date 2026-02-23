#!/bin/bash

source ./build-opts.sh

# Install prefix
PREFIX=$GKYLSOFT/sundials_install
# Location where dependency sources will be downloaded
DEP_SOURCES=$GKYLSOFT/dep_src/

mkdir -p $DEP_SOURCES
cd $DEP_SOURCES

if [ "$DOWNLOAD_PKGS" = "yes" ]
then
    echo "Downloading SUNDIALS .."
    # delete old checkout and builds
    rm -rf sundials
    git clone --branch=feature/fes-get-stage-index https://github.com/LLNL/sundials
fi

if [ "$BUILD_PKGS" = "yes" ]
then
    echo "Building SUNDIALS .."
    mkdir -p sundials/build
    cd sundials/build
    echo 'C compiler: $MPICC'
    echo 'C++ compiler: $MPICXX'
    echo 'Fortran compiler: $MPIFC'
    echo 'MPI executor: $MPIEXEC'

    cmake -DCMAKE_INSTALL_PREFIX=$PREFIX -DENABLE_MPI=ON -DSUNDIALS_INDEX_SIZE=32 -DCMAKE_C_COMPILER=$MPICC -DCMAKE_CXX_COMPILER=$MPICXX -DCMAKE_CXX_COMPILER=$MPICXX -DMPI_C_COMPILER=$MPICC -DMPIEXEC_EXECUTABLE=$MPIEXEC -DENABLE_LAPACK=ON -DLAPACK_LIBRARIES=$GKYLSOFT/OpenBLAS/lib/libopenblas.a -DMPI_Fortran_WORKS=ON -DMPI_Fortran_COMPILER=$MPIFC ..

    make -j 4 install

    # soft-link 
    ln -sfn $PREFIX $GKYLSOFT/sundials
fi
