#!/bin/bash

source ./build-opts.sh

# Edit to suite your system
# Location where dependency sources will be downloaded
DEP_SOURCES=$GKYLSOFT/dep_src/

CUDA_VERSION_THRESH="13.0"
CUDA_VERSION=$(nvcc --version 2>/dev/null | grep release | cut -d ' ' -f 5 | cut -c1- | cut -d ',' -f 1)

mkdir -p $DEP_SOURCES
cd $DEP_SOURCES

if [ "$DOWNLOAD_PKGS" = "yes" ]
then
    echo "Downloading cuDSS .."
    # delete old checkout and builds
    rm -rf cuDSS-*

    if awk -v v1="$CUDA_VERSION" -v v2="$CUDA_VERSION_THRESH" 'BEGIN {
    split(v1, a, "."); split(v2, b, ".");
    for (i=1; i<=3; i++) {
        if (a[i] > b[i]) {exit 0} else if (a[i] < b[i]) {exit 1}
    }
    exit 0
}'; then
        # Cuda >= 13.0
        PREFIX=$GKYLSOFT/cuDSS-0.7.1.4
        wget https://developer.download.nvidia.com/compute/cudss/redist/libcudss/linux-x86_64/libcudss-linux-x86_64-0.7.1.4_cuda13-archive.tar.xz
    else
        # Cuda < 13.0
        PREFIX=$GKYLSOFT/cuDSS-0.7.0.20
        wget https://developer.download.nvidia.com/compute/cudss/redist/libcudss/linux-x86_64/libcudss-linux-x86_64-0.7.0.20_cuda12-archive.tar.xz
    fi
fi

if [ "$BUILD_PKGS" = "yes" ]
then
    echo "Building cuDSS .."
    if awk -v v1="$CUDA_VERSION" -v v2="$CUDA_VERSION_THRESH" 'BEGIN {
    split(v1, a, "."); split(v2, b, ".");
    for (i=1; i<=3; i++) {
        if (a[i] > b[i]) {exit 0} else if (a[i] < b[i]) {exit 1}
    }
    exit 0
}'; then
        # Cuda >= 13.0
        tar xf libcudss-linux-x86_64-0.7.1.4_cuda13-archive.tar.xz
        mv libcudss-linux-x86_64-0.7.1.4_cuda13-archive $PREFIX
    else
        # Cuda < 13.0
        tar xf libcudss-linux-x86_64-0.7.0.20_cuda12-archive.tar.xz
        mv libcudss-linux-x86_64-0.7.0.20_cuda12-archive $PREFIX
    fi

    # soft-link 
    ln -sfn $PREFIX $GKYLSOFT/cuDSS
fi
