module unload darshan
module load cray-mpich/9.0.1
module load cray-libsci/25.09.0

: "${PREFIX:=$HOME/gkylsoft}"

cd install-deps
./mkdeps.sh --build-superlu=yes --build-luajit=yes --prefix=$PREFIX
