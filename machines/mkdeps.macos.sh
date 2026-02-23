cd install-deps
: "${PREFIX:=/Users/mfrancis/Documents/gkeyll/code/gkeyll_v3/gkylsoft}"
export MACOSX_DEPLOYMENT_TARGET=14.7
./mkdeps.sh --build-openblas=no --build-superlu=yes --build-luajit=yes --build-openmpi=yes --prefix=$PREFIX
