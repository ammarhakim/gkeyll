module load nvhpc/25.3
: "${PREFIX:=$HOME/gkylsoft}"

cd install-deps
./mkdeps.sh --build-openblas=yes --build-superlu=yes --build-cudss=yes --prefix=$PREFIX --build-luajit=yes MPICC=mpicc  MPICXX=mpicxx
