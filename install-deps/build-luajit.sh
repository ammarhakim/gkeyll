#!/bin/bash

source ./build-opts.sh

# Install prefix
PREFIX=$GKYLSOFT/luajit-git
# Location where dependency sources will be downloaded
DEP_SOURCES=$GKYLSOFT/dep_src/

mkdir -p $DEP_SOURCES
cd $DEP_SOURCES

if [ "$DOWNLOAD_PKGS" = "yes" ]
then
    echo "Downloading LuaJIT .."
    # delete old checkout and builds
    rm -rf luajit
    git clone https://github.com/LuaJIT/LuaJIT.git luajit
fi

if [ "$BUILD_PKGS" = "yes" ]
then
    echo "Building LuaJIT .."
    cd luajit

    make -j6 PREFIX=$PREFIX install 

    # soft-link 
    ln -sfn $PREFIX $GKYLSOFT/luajit
fi
