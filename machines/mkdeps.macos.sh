cd install-deps
: "${PREFIX:=/Users/mfrancis/Documents/gkeyll/code/gkeyll_v2/gkylsoft}"
export MACOSX_DEPLOYMENT_TARGET=14.7
./mkdeps.sh --prefix=$PREFIX --build-openblas=no --build-superlu=yes --build-luajit=yes --build-openmpi=yes --build-sundials=yes
