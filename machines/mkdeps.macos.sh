cd install-deps
: "${PREFIX:=/Users/mfrancis/Documents/gkeyll/code/gkeyll_v0/gkylsoft}"
export MACOSX_DEPLOYMENT_TARGET=14.7
./mkdeps.sh --build-openblas=no --build-superlu=no --build-luajit=no --prefix=$PREFIX --build-openmpi=yes
