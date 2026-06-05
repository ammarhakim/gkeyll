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
    git clone --branch=develop https://github.com/LLNL/sundials
fi

if [ "$BUILD_PKGS" = "yes" ]
then
    echo "Building SUNDIALS .."
    mkdir -p sundials/build
    cd sundials/build
    echo $MPICC
    echo $MPICXX

    # Choose LAPACK/BLAS option depending on platform
    if [[ "$(uname -s)" == "Darwin" ]]; then
        LAPACK_OPTION="-DBLA_VENDOR=Apple"
    else
        LAPACK_OPTION="-DLAPACK_LIBRARIES=${GKYLSOFT}/OpenBLAS/lib/libopenblas.a"
    fi

    cmake -DCMAKE_INSTALL_PREFIX=$PREFIX -DSUNDIALS_ENABLE_MPI=ON -DSUNDIALS_INDEX_SIZE=32 -DCMAKE_C_COMPILER=$MPICC -DCMAKE_CXX_COMPILER=$MPICXX -DMPI_C_COMPILER=$MPICC -DMPIEXEC_EXECUTABLE=$MPIEXEC -DSUNDIALS_ENABLE_LAPACK=ON ${LAPACK_OPTION} -DMPI_Fortran_WORKS=ON -DMPI_Fortran_COMPILER=$GKYLSOFT/openmpi/bin/mpifort ..

    make -j 4 install

    # soft-link 
    ln -sfn $PREFIX $GKYLSOFT/sundials
fi
