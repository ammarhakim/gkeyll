cd install-deps
: "${PREFIX:=/opt/gkylsoft}"
./mkdeps.sh --build-openblas=yes --build-superlu=yes --build-openmpi=yes --build-luajit=yes --build-cudss=yes --prefix=$PREFIX
