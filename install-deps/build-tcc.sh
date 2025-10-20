#!/bin/bash

source ./build-opts.sh

# Install prefix
PREFIX=$GKYLSOFT/tcc-git
# Location where dependency sources will be downloaded
DEP_SOURCES=$GKYLSOFT/dep_src/

mkdir -p $DEP_SOURCES
cd $DEP_SOURCES

if [ "$DOWNLOAD_PKGS" = "yes" ]
then
    echo "Downloading tcc .."
    # delete old checkout and builds
    rm -rf tinycc
    git clone https://github.com/TinyCC/tinycc.git
fi

if [ "$BUILD_PKGS" = "yes" ]
then
    echo "Building tcc .."
    cd tinycc
    ./configure --prefix=$PREFIX
    make -j6 
    make install

    # soft-link 
    ln -sfn $PREFIX $GKYLSOFT/tcc
fi
